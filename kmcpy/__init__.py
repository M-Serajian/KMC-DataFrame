"""
kmcpy v1.0.0
============
Python bindings for the KMC k-mer counter.
Stage 2 output is always in RAM. Stage 1 bins spill to disk only when
the data exceeds the available RAM budget (hybrid mode).

Built on KMC 3.2.4 (https://github.com/refresh-bio/KMC).

Internally, all k-mers are stored in 2-bit packed format (A=0, C=1, G=2, T=3).
The ``decoded`` parameter controls what you get back:

    decoded=True  (default) — ACGT letter strings  →  pandas DataFrame
    decoded=False           — kmer_0..kmer_{n-1} (n=ceil(2k/64) uint64 cols) → space-optimal, all k

Quick start
-----------
    import kmcpy

    # Decoded ACGT strings (default)
    df = kmcpy.count_kmers("assembly.fasta", k=15,
                           input_fmt="mfasta",
                           tmp_dir="/scratch/kmc_tmp")

    # Packed uint64 key — zero overhead, for cuDF merge
    df = kmcpy.count_kmers("assembly.fasta", k=15,
                           input_fmt="mfasta",
                           tmp_dir="/scratch/kmc_tmp",
                           decoded=False)

    # Show help
    kmcpy.help()

For full documentation:  help(kmcpy.count_kmers)
"""

from ._version import __version__, __kmc_engine_version__
from . import _core as _ext

import os
import sys
import uuid
import shutil
import signal
import atexit
import time
import numpy as np
import pandas as pd
from typing import Union, List, Optional


# ---------------------------------------------------------------------------
# Input format map  (mirrors KMC -f flags)
# ---------------------------------------------------------------------------
_INPUT_FMT = {
    "fastq"   : 0,   # -fq   (KMC default)
    "fq"      : 0,
    "fastq.gz": 0,
    "fasta"   : 1,   # -fa   single-line FASTA
    "fa"      : 1,
    "fna"     : 1,
    "fasta.gz": 1,
    "mfasta"  : 2,   # -fm   multi-line FASTA (most assembler outputs)
    "fm"      : 2,
    "bam"     : 3,   # -fbam
    "kmc"     : 4,   # -fkmc (use existing KMC database as input)
}

# ---------------------------------------------------------------------------
# Module-level stats from the most recent run
# ---------------------------------------------------------------------------
_last_stats: dict = {}

# ---------------------------------------------------------------------------
# Global registry of active subdirectories owned by this process.
# Used by atexit and signal handlers to clean up on abnormal exit.
# Protected by the GIL — no explicit lock needed for CPython.
# ---------------------------------------------------------------------------
_active_subdirs: set = set()


# ---------------------------------------------------------------------------
# Temp directory helpers
# ---------------------------------------------------------------------------

def _register_subdir(path: str) -> None:
    """Register a subdir for cleanup on exit."""
    _active_subdirs.add(path)


def _unregister_subdir(path: str) -> None:
    """Unregister a subdir (already cleaned up normally)."""
    _active_subdirs.discard(path)


def _cleanup_subdir(path: str) -> None:
    """
    Remove a run subdir safely on both local and NFS filesystems.

    Problem: KMC uses internal C++ threads that may still be flushing
    bin files when the Python-visible call returns.  Calling unlink()
    on a file that a thread still has open raises EBUSY on some
    filesystems (NFS, GPFS, Lustre).

    Solution — atomic rename then deferred delete:
      1. Atomically rename the subdir to a hidden "tombstone" name.
         This is instant and always succeeds even on busy files.
         It frees the original path immediately so the next call can
         reuse the slot.  KMC threads keep writing to the tombstone
         path harmlessly — the kernel keeps the inode alive.
      2. Register an atexit hook to delete the tombstone.
         By the time atexit runs, all KMC threads have exited and
         all file descriptors are closed.  rmtree succeeds cleanly.
         On NFS, any .nfs* silly-rename files will also be gone by then.

    This approach never calls unlink() on an open file, so EBUSY
    cannot occur.
    """
    if not os.path.isdir(path):
        return

    # Step 1: atomically rename to a tombstone path in the same directory.
    # Same-directory rename is guaranteed atomic on POSIX.
    parent    = os.path.dirname(path)
    tombstone = os.path.join(parent,
                             ".kmcpy_tomb_{}_{}".format(
                                 os.getpid(), uuid.uuid4().hex))
    try:
        os.rename(path, tombstone)
    except OSError:
        # rename failed (cross-device, permissions, etc.) —
        # fall back to best-effort rmtree
        tombstone = path

    # Step 2: deferred delete via atexit.
    # All KMC threads and file descriptors are closed by then.
    def _deferred_rmtree(p):
        try:
            shutil.rmtree(p, ignore_errors=True)
        except Exception:
            pass

    atexit.register(_deferred_rmtree, tombstone)


def _cleanup_all_active() -> None:
    """
    atexit / signal handler — clean up all subdirs owned by this process.
    Called on: normal exit, SIGINT, SIGTERM.
    Not called on: SIGKILL (unavoidable).
    """
    for path in list(_active_subdirs):
        _cleanup_subdir(path)
    _active_subdirs.clear()


def _signal_handler(signum, frame) -> None:
    """Handle SIGINT and SIGTERM — clean up then re-raise."""
    _cleanup_all_active()
    # Re-raise as default to allow normal termination behaviour
    signal.signal(signum, signal.SIG_DFL)
    os.kill(os.getpid(), signum)


# Register cleanup hooks at import time — runs once per process
atexit.register(_cleanup_all_active)
signal.signal(signal.SIGINT,  _signal_handler)
signal.signal(signal.SIGTERM, _signal_handler)


def _purge_stale_subdirs(tmp_dir: str, max_age_hours: float = 24.0) -> None:
    """
    Scan tmp_dir for orphaned kmcpy_<PID>_<UUID>/ subdirectories left
    by processes that were SIGKILL'd (the only case we cannot clean up
    in the signal handler).

    A subdir is considered stale when ALL of:
      1. Its name matches the kmcpy_<PID>_<UUID> pattern.
      2. The PID embedded in its name is no longer alive.
      3. Its mtime is older than max_age_hours.

    Condition 3 is a safety guard — we never delete a dir whose PID
    happens to have been recycled by a new process.
    """
    now = time.time()
    cutoff = max_age_hours * 3600.0
    try:
        for entry in os.scandir(tmp_dir):
            if not entry.is_dir():
                continue
            name = entry.name
            if not name.startswith("kmcpy_"):
                continue
            # Expected pattern: kmcpy_<PID>_<UUID>
            parts = name.split("_", 2)   # ['kmcpy', '<PID>', '<UUID>']
            if len(parts) != 3:
                continue
            try:
                pid = int(parts[1])
            except ValueError:
                continue
            # Check age
            try:
                mtime = entry.stat().st_mtime
            except OSError:
                continue
            if (now - mtime) < cutoff:
                continue   # Too young — might still be active
            # Check if PID is alive
            pid_alive = True
            try:
                os.kill(pid, 0)   # Signal 0: existence check, no actual signal
            except ProcessLookupError:
                pid_alive = False
            except PermissionError:
                pid_alive = True  # Exists but owned by another user — leave it
            if not pid_alive:
                # Safe to remove — PID is gone, NFS locks are released
                shutil.rmtree(entry.path, ignore_errors=True)
                # If .nfs* files remain after rmtree, wait briefly and retry
                if os.path.isdir(entry.path):
                    time.sleep(1.0)
                    shutil.rmtree(entry.path, ignore_errors=True)
    except OSError:
        pass   # tmp_dir not accessible — silently skip


def _make_run_subdir(tmp_dir: str) -> str:
    """
    Create and return a unique subdirectory for this specific KMC run:

        tmp_dir/kmcpy_<PID>_<UUID>/

    PID makes it easy to identify ownership.
    UUID guarantees uniqueness even if PID is recycled or multiple
    calls run in parallel within the same process.
    """
    name    = "kmcpy_{}_{}".format(os.getpid(), uuid.uuid4().hex)
    subdir  = os.path.join(tmp_dir, name)
    os.makedirs(subdir, exist_ok=False)
    return subdir


# ---------------------------------------------------------------------------
# help()
# ---------------------------------------------------------------------------

def help() -> None:
    """Print a formatted help reference for the kmcpy package."""
    _BOLD  = "\033[1m"  if sys.stdout.isatty() else ""
    _CYAN  = "\033[36m" if sys.stdout.isatty() else ""
    _GREEN = "\033[32m" if sys.stdout.isatty() else ""
    _RESET = "\033[0m"  if sys.stdout.isatty() else ""
    _DIM   = "\033[2m"  if sys.stdout.isatty() else ""

    def header(text):
        print("\n{}{}{}{}".format(_BOLD, _CYAN, text, _RESET))

    def row(flag, pyarg, default, desc):
        print("  {:<10}  {:<20}  {:<12}  {}".format(
            _GREEN + flag + _RESET,
            _DIM   + pyarg + _RESET,
            str(default),
            desc))

    print()
    print("{}{}kmcpy {}{}  —  KMC k-mer counting, hybrid RAM/disk{}".format(
        _BOLD, _CYAN, __version__, _RESET, _RESET))
    print("{}KMC engine {}{}".format(_DIM, __kmc_engine_version__, _RESET))
    print()
    print("  {}import kmcpy{}".format(_BOLD, _RESET))
    print("  {}df = kmcpy.count_kmers(input_files, tmp_dir=..., **options){}".format(
        _BOLD, _RESET))

    header("Input")
    print("  {:<10}  {:<20}  {:<12}  {}".format(
        "argument", "Python parameter", "default", "description"))
    print("  " + "-"*72)
    row("(pos)",  "input_files",   "",        "str or list[str] — FASTQ/FASTA/BAM paths")
    row("(req)",  "tmp_dir",       "",        "base directory for temporary files (REQUIRED)")
    row("-fq",    "input_fmt",     "'fastq'", "FASTQ input (gzipped or not)")
    row("-fa",    "input_fmt",     "",        "single-line FASTA input")
    row("-fm",    "input_fmt",     "",        "multi-line FASTA (assembler output)")
    row("-fbam",  "input_fmt",     "",        "BAM input")
    row("-fkmc",  "input_fmt",     "",        "existing KMC database as input")

    header("k-mer options")
    print("  {:<10}  {:<20}  {:<12}  {}".format(
        "KMC flag", "Python parameter", "default", "description"))
    print("  " + "-"*72)
    row("-k<len>", "k",            25,    "k-mer length (1–256; k<=13 uses small-k path)")
    row("-b",      "canonical",    True,  "False = count non-canonical k-mers (-b)")
    row("-hc",     "homopolymer",  False, "homopolymer-compressed k-mers (experimental)")
    row("-ci<n>",  "min_count",    1,     "exclude k-mers with count < n  (default 1 = keep all)")
    row("-cx<n>",  "max_count",    "1e9", "exclude k-mers with count > n")
    row("-cs<n>",  "counter_max",  65535, "saturate counter at n (max 65535)")

    header("Performance options")
    print("  {:<10}  {:<20}  {:<12}  {}".format(
        "KMC flag", "Python parameter", "default", "description"))
    print("  " + "-"*72)
    row("-m<GB>",  "max_ram_gb",    12,    "RAM budget in gigabytes")
    row("-t<n>",   "threads",       0,     "total threads (0 = all CPU cores)")
    row("-sm",     "strict_mem",    True,  "strict memory mode — never exceed -m limit")
    row("-p<n>",   "signature_len", 9,     "signature length 5–11 (advanced tuning)")
    row("-n<n>",   "n_bins",        0,     "number of bins (0 = auto)")
    row("-sf<n>",  "n_readers",     0,     "FASTQ reader threads (0 = auto)")
    row("-sp<n>",  "n_splitters",   0,     "splitter threads (0 = auto)")

    header("Disk / memory options")
    print()
    print("  {}tmp_dir{}  (REQUIRED)".format(_BOLD, _RESET))
    print("    Base directory for temporary files.")
    print("    On HPC: use /scratch or /localscratch, not /tmp.")
    print("    A unique subdirectory kmcpy_<PID>_<UUID>/ is created")
    print("    inside tmp_dir for each call and removed automatically")
    print("    on completion, exception, SIGINT, or SIGTERM.")
    print("    On SIGKILL: stale dirs are cleaned on the next call.")
    print()
    print("  {}Hybrid spill strategy:{}".format(_BOLD, _RESET))
    print("    Stage 1 bins stay in RAM while memory allows (ramOnlyMode).")
    print("    When RAM budget (max_ram_gb) is exceeded, KMC automatically")
    print("    spills bins to tmp_dir/kmcpy_<PID>_<UUID>/.")
    print("    Stage 2 output (k-mer counts) always comes back to RAM.")
    print("    strict_mem=True (default) ensures max_ram_gb is never exceeded.")

    header("Output representation  ( decoded parameter )")
    print()
    print("  k-mers are always stored internally as 2-bit integer codes:")
    print()
    print("    {}A = 0    C = 1    G = 2    T = 3{}".format(_BOLD, _RESET))
    print()
    print("  The {}decoded{} parameter controls what you receive:".format(
        _BOLD, _RESET))
    print()
    row("",  "drop_count",  False, "True = skip counts at C++ level — no allocation, no transfer, no column")
    print()
    print("  {}decoded=True{}  (default)".format(_BOLD, _RESET))
    print("    Each 2-bit code is converted to its ACGT letter.")
    print("    Returns a pandas DataFrame with columns:")
    print("      kmer   (StringDtype)  — ACGT string of length k")
    print("      count  (uint32)       — occurrence count")
    print("    Sorted lexicographically.  Ready for export or inspection.")
    print()
    print("  {}decoded=False{}".format(_BOLD, _RESET))
    print("    Space-optimal packed key — n_words = ceil(k*2/64) uint64 columns.")
    print("    Returns a pandas DataFrame with columns:")
    print("      kmer_0          (uint64) — most significant word (base[0] side)")
    print("      kmer_1          (uint64) — next word")
    print("      ...             ...")
    print("      kmer_{n_words-1} (uint64) — least significant word (base[k-1] side)")
    print("      count           (uint32) — occurrence count")
    print("    Bit layout: base[k-1] at bits 1:0 of kmer_{{n-1}}; base[0] at MSB.")
    print("    Unused high bits of kmer_0 are always zero.")
    print("    A=00, C=01, G=10, T=11 (exact KMC 2-bit encoding, zero transformation).")
    print("    Space efficiency (minimum columns for k bits):")
    print("      k<=32:  1 col   k<=64:  2 cols  k<=128: 4 cols  k<=256: 8 cols")
    print("    Valid for ALL k (1-256). No arbitrary cutoff.")
    print("    cuDF merge (same syntax for all k):")
    print("      key_cols = [f\"kmer_{{i}}\" for i in range(n_words)]")
    print("      merged = gdf1.merge(gdf2, on=key_cols)")
    print()
    print("  {}drop_count=False{}  (default)".format(_BOLD, _RESET))
    print("    Count column is included in the returned DataFrame.")
    print("    C++ allocates, fills, and transfers counts_arr (n x uint32).")
    print()
    print("  {}drop_count=True{}".format(_BOLD, _RESET))
    print("    True memory optimisation at the C++->Python boundary.")
    print("    counts_arr is allocated with size 0 — zero heap allocation.")
    print("    Count bytes in bin records are never read or decoded.")
    print("    No count data crosses the C++->Python boundary.")
    print("    count column is absent from the returned DataFrame.")
    print("    Memory saved at C++->Python boundary:")
    print("      10M  k-mers ->   40 MB saved")
    print("      100M k-mers ->  400 MB saved")
    print("      1B   k-mers -> 4000 MB saved")
    print("      8 parallel jobs x 100M k-mers -> 3.2 GB saved")
    print("    Use when counts are not needed downstream:")
    print("      set intersection, Jaccard similarity, SBWT construction,")
    print("      presence/absence matrix, any count-agnostic analysis.")
    print()
    print("  {}2-bit encoding reference (KMC internal format):{}".format(_DIM, _RESET))
    print("    Bases packed MSB-first, 4 per byte.")
    print("    Suffix: ceil((k - lut_prefix_len) / 4) bytes.")
    print("    Prefix: stored as uint64 LUT row index.")

    header("Returns")
    print("  decoded=True,  drop_count=False  -> kmer (StringDtype), count (uint32)")
    print("  decoded=True,  drop_count=True   -> kmer (StringDtype) only")
    print("  decoded=False, drop_count=False  -> kmer_0..kmer_{{n-1}} (uint64), count (uint32)")
    print("  decoded=False, drop_count=True   -> kmer_0..kmer_{{n-1}} (uint64) only")
    print("                   n = ceil(k*2/64) — space-optimal, valid for all k 1-256")

    header("Temp directory lifecycle")
    print("  Each call to count_kmers() creates:")
    print("    tmp_dir/kmcpy_<PID>_<UUID>/")
    print("  Removed automatically on:")
    print("    normal return    -> finally: block")
    print("    exception        -> finally: block")
    print("    SIGINT (Ctrl+C)  -> signal handler + finally: block")
    print("    SIGTERM          -> signal handler + finally: block")
    print("    SIGKILL          -> cannot catch; cleaned on next call (stale PID detection)")

    header("Parallel usage")
    print("  Each call is fully independent — safe to run concurrently.")
    print("  Each call gets its own unique subdir — no cross-contamination.")
    print("  Budget: n_parallel_calls x max_ram_gb must fit in total RAM.")
    print()
    print("  from concurrent.futures import ProcessPoolExecutor")
    print("  import functools, kmcpy")
    print("  fn = functools.partial(kmcpy.count_kmers, k=15,")
    print("                         input_fmt='mfasta', tmp_dir='/scratch')")
    print("  with ProcessPoolExecutor(max_workers=4) as pool:")
    print("      results = list(pool.map(fn, genome_files))")

    header("Other functions")
    print("  kmcpy.last_run_stats()  — dict of timing and count statistics")
    print("  kmcpy.kmc_version()     — KMC engine version string")
    print("  kmcpy.help()            — this message")

    header("Examples")
    examples = [
        ("# ACGT strings (default)",
         "df = kmcpy.count_kmers('assembly.fasta', k=15,\n"
         "                        input_fmt='mfasta',\n"
         "                        tmp_dir='/scratch/kmc_tmp')"),
        ("# Packed uint64 key — zero overhead, for cuDF merge",
         "df = kmcpy.count_kmers('assembly.fasta', k=15,\n"
         "                        input_fmt='mfasta',\n"
         "                        tmp_dir='/scratch/kmc_tmp',\n"
         "                        decoded=False)\n"
         "# n_words = ceil(k*2/64); same syntax for all k:\n"
         "# key_cols = [f'kmer_{{i}}' for i in range(n_words)]\n"
         "# merged = gdf1.merge(gdf2, on=key_cols)"),
        ("# Paired FASTQ, k=31, 8 threads, 32 GB RAM",
         "df = kmcpy.count_kmers(['r1.fastq.gz', 'r2.fastq.gz'], k=31,\n"
         "                        tmp_dir='/scratch/kmc_tmp',\n"
         "                        threads=8, max_ram_gb=32)"),
        ("# Equivalent to: kmc -k27 -m24 -t8 -ci1 -cx500 -fm sample.fasta",
         "df = kmcpy.count_kmers('sample.fasta', k=27, input_fmt='mfasta',\n"
         "                        tmp_dir='/scratch/kmc_tmp',\n"
         "                        max_ram_gb=24, threads=8,\n"
         "                        min_count=1, max_count=500)"),
    ]
    for comment, code in examples:
        print()
        print("  {}{}{}".format(_DIM, comment, _RESET))
        for line in code.split("\n"):
            print("  {}{}{}".format(_BOLD, line, _RESET))

    print()
    print("{}Note:{} tmp_dir is REQUIRED — use a fast, large filesystem "
          "(e.g. /scratch on HPC).".format(_BOLD, _RESET))
    print("{}Note:{} min_count defaults to 1 (keep all k-mers). "
          "KMC CLI default is 2.".format(_BOLD, _RESET))
    print("{}Note:{} Most assembler FASTA files are multi-line. "
          "Use input_fmt='mfasta', not 'fasta'.".format(_BOLD, _RESET))
    print("{}Note:{} k <= 13 uses KMC small-k optimisation path — "
          "fully supported.".format(_BOLD, _RESET))
    print()


# ---------------------------------------------------------------------------
# count_kmers()
# ---------------------------------------------------------------------------

def count_kmers(
    input_files  : Union[str, List[str]],
    tmp_dir      : str,                        # REQUIRED — no default
    # ---- Stage 1 (KMC -f / -k / -m / -t / -p / -n / -sf / -sp / -b / -hc) ----
    k            : int   = 25,
    input_fmt    : str   = "fastq",
    max_ram_gb   : int   = 12,
    threads      : int   = 0,
    canonical    : bool  = True,
    homopolymer  : bool  = False,
    signature_len: int   = 9,
    n_bins       : int   = 0,
    n_readers    : int   = 0,
    n_splitters  : int   = 0,
    # ---- Stage 2 (KMC -ci / -cx / -cs / -sm) --------------------------------
    min_count    : int   = 1,
    max_count    : int   = 1_000_000_000,
    counter_max  : int   = 65535,
    strict_mem   : bool  = True,
    # ---- Output ---------------------------------------------------------------
    decoded      : bool  = True,
    drop_count   : bool  = False,
) -> pd.DataFrame:
    """
    Count k-mers from sequencing data.

    Stage 1 bins are kept in RAM while the memory budget allows.  When
    ``max_ram_gb`` is exceeded, KMC automatically spills bins to a unique
    subdirectory inside ``tmp_dir``.  Stage 2 output (k-mer counts) always
    comes back to RAM as a pandas DataFrame.

    The temporary subdirectory is created on entry and removed on exit —
    whether the call succeeds, raises an exception, receives SIGINT, or
    receives SIGTERM.  Subdirectories orphaned by SIGKILL are detected and
    removed on the next call to this function.

    Parameters
    ----------
    input_files : str or list[str]
        Path to one or more input files.  All files must share the same
        format (see ``input_fmt``).  Gzipped files are supported.

    tmp_dir : str  (REQUIRED)
        Base directory for temporary files.  Must exist and be writable.
        On HPC systems use a dedicated scratch filesystem (e.g. /scratch),
        not /tmp which is often RAM-backed or very small.

        A unique subdirectory ``kmcpy_<PID>_<UUID>/`` is created inside
        ``tmp_dir`` for each call and removed automatically on completion.

    k : int, default 25
        K-mer length (1–256).  k <= 13 uses KMC small-k optimisation.
        (KMC flag: ``-k<len>``)

    input_fmt : str, default "fastq"
        Input file format:

        ============  ==============================================
        ``"fastq"``   FASTQ, gzipped or plain       (``-fq``)
        ``"fasta"``   Single-line FASTA              (``-fa``)
        ``"mfasta"``  Multi-line FASTA (assemblies)  (``-fm``)
        ``"bam"``     BAM                           (``-fbam``)
        ``"kmc"``     Existing KMC database         (``-fkmc``)
        ============  ==============================================

        .. warning::
            Most assemblers produce multi-line FASTA.  Use
            ``input_fmt='mfasta'`` or k-mers will be silently truncated.

    max_ram_gb : int, default 12
        RAM budget in gigabytes.  KMC respects this limit — when bins
        exceed it, they spill to ``tmp_dir``.  (KMC flag: ``-m<size>``)

    threads : int, default 0
        Total parallel threads.  0 = all CPU cores.
        (KMC flag: ``-t<value>``)

    canonical : bool, default True
        Count canonical k-mers (k-mer + reverse complement together).
        (KMC flag: ``-b`` disables)

    homopolymer : bool, default False
        Homopolymer-compressed k-mers.  Experimental.
        (KMC flag: ``-hc``)

    signature_len : int, default 9
        Internal binning signature length (5–11).  Rarely needs changing.
        (KMC flag: ``-p<par>``)

    n_bins : int, default 0
        Number of bins.  0 = auto.  (KMC flag: ``-n<value>``)

    n_readers : int, default 0
        Reader threads.  0 = auto.  (KMC flag: ``-sf<value>``)

    n_splitters : int, default 0
        Splitter threads.  0 = auto.  (KMC flag: ``-sp<value>``)

    min_count : int, default 1
        Exclude k-mers appearing fewer than this many times.
        Default 1 keeps all k-mers including singletons.
        KMC CLI default is 2.  (KMC flag: ``-ci<value>``)

    max_count : int, default 1_000_000_000
        Exclude k-mers appearing more than this many times.
        (KMC flag: ``-cx<value>``)

    counter_max : int, default 65535
        Saturate counter at this value.  Max 65535.
        Default 65535 = no effective saturation for typical datasets.
        (KMC flag: ``-cs<value>``)

    strict_mem : bool, default True
        Enforce strict memory limit during Stage 2 sorting.
        When True, KMC splits oversized bins into sub-bins rather than
        exceeding ``max_ram_gb``.  Recommended for large datasets.
        (KMC flag: ``-sm``)

    drop_count : bool, default False
        If True, the occurrence count is neither computed nor transferred
        from C++ to Python.

        **This is a true memory optimisation at the C++→Python boundary,
        not merely a cosmetic column drop.**

        When ``drop_count=False`` (default), the C++ layer allocates a
        ``uint32`` array of length *n* (one entry per unique k-mer),
        fills it with count values decoded from the KMC bin format, and
        transfers it to Python as a numpy array.

        When ``drop_count=True``:

        * The ``counts_arr`` numpy array is allocated with size **0** —
          no heap memory is reserved for counts at all.
        * Count bytes in the bin records are never read or decoded.
        * No count data crosses the C++→Python boundary.
        * The ``count`` column is absent from the returned DataFrame.

        Memory saved at the C++→Python boundary:

        ============  ===========  ==============================
        Unique k-mers Saved (MB)   Example scenario
        ============  ===========  ==============================
        10 M          40 MB        small genome, k=15
        100 M         400 MB       large genome, k=25
        1 B           4,000 MB     100k genomes combined
        ============  ===========  ==============================

        For 8 parallel jobs each with 100M k-mers: **3.2 GB saved**.

        Use ``drop_count=True`` when only k-mer presence or identity
        matters, not frequency — for example:

        * Set intersection / union across samples
        * Jaccard similarity between genomes
        * Presence/absence matrix construction
        * SBWT index construction
        * Any downstream tool that ignores counts

    decoded : bool, default True
        Controls k-mer representation in the returned DataFrame.

        k-mers are always stored **internally as 2-bit codes**:

        .. code-block:: text

            A = 0    C = 1    G = 2    T = 3
            Packed MSB-first, 4 bases per byte.

        ``decoded=True`` *(default)*
            Returns DataFrame with columns:
            - **kmer**  (StringDtype) — ACGT string of length k
            - **count** (uint32)     — occurrence count
            Sorted lexicographically.

        ``decoded=False``
            Space-optimal packed key representation.  Valid for all k (1–256).

            Returns DataFrame with columns:
            - **kmer_0** … **kmer_{n-1}** (uint64) — packed key words,
              where n = ``ceil(k * 2 / 64)``.
              kmer_0 holds the most significant bits (base[0] side).
              kmer_{n-1} holds the least significant bits (base[k-1] side).
              base[k-1] occupies bits 1:0 of kmer_{n-1}.
              Unused high bits of kmer_0 are always zero.
              A=00, C=01, G=10, T=11 (exact KMC 2-bit encoding).
            - **count** (uint32) — occurrence count

            Space efficiency — minimum columns for k:

            =========  =======  =====================
            k range    n_words  bytes per k-mer key
            =========  =======  =====================
            k ≤  32    1        8  (≤ 2 bits wasted)
            k ≤  64    2        16 (≤ 2 bits wasted)
            k ≤  96    3        24
            k ≤ 128    4        32
            k ≤ 192    6        48
            k ≤ 256    8        64
            =========  =======  =====================

            cuDF merge (identical syntax for all k):

            .. code-block:: python

                n_words  = len([c for c in df.columns if c.startswith("kmer_")])
                key_cols = [f"kmer_{i}" for i in range(n_words)]
                merged   = gdf1.merge(gdf2, on=key_cols)

    Returns
    -------
    pandas.DataFrame
        drop_count=False (default):
          decoded=True:  kmer (StringDtype), count (uint32)
          decoded=False: kmer_0..kmer_{n-1} (uint64), count (uint32)
        drop_count=True:
          decoded=True:  kmer (StringDtype) only
          decoded=False: kmer_0..kmer_{n-1} (uint64) only
        where n = ceil(k * 2 / 64).  Valid for all k 1–256.

    Raises
    ------
    TypeError
        If ``tmp_dir`` is not provided.
    ValueError
        If ``tmp_dir`` does not exist or is not writable.
    FileNotFoundError
        If any input file does not exist.
    ValueError
        If ``k`` < 1 or ``input_fmt`` is not recognised.

    Examples
    --------
    Count 15-mers from a multi-line FASTA assembly:

    >>> df = kmcpy.count_kmers("assembly.fasta", k=15,
    ...                         input_fmt="mfasta",
    ...                         tmp_dir="/scratch/kmc_tmp")

    Packed uint64 key for cuDF merge (zero overhead):

    >>> df = kmcpy.count_kmers("assembly.fasta", k=15,
    ...                         input_fmt="mfasta",
    ...                         tmp_dir="/scratch/kmc_tmp",
    ...                         decoded=False)
    >>> gdf = cudf.DataFrame.from_pandas(df)
    >>> # Identical syntax for all k:
    >>> n_words  = len([c for c in df.columns if c.startswith("kmer_")])
    >>> key_cols = [f"kmer_{i}" for i in range(n_words)]
    >>> merged   = gdf1.merge(gdf2, on=key_cols)

    Parallel execution over multiple genomes:

    >>> from concurrent.futures import ProcessPoolExecutor
    >>> import functools
    >>> fn = functools.partial(kmcpy.count_kmers, k=15,
    ...                        input_fmt="mfasta",
    ...                        tmp_dir="/scratch/kmc_tmp",
    ...                        threads=4, max_ram_gb=8)
    >>> with ProcessPoolExecutor(max_workers=8) as pool:
    ...     results = list(pool.map(fn, genome_files))

    Notes
    -----
    * ``tmp_dir`` must be a fast, large filesystem — not RAM-backed /tmp.
    * Stage 1 bins spill to ``tmp_dir/kmcpy_<PID>_<UUID>/`` only when
      ``max_ram_gb`` is exceeded — otherwise everything stays in RAM.
    * Stage 2 output always comes to RAM regardless of data size.
    * ``strict_mem=True`` (default) guarantees ``max_ram_gb`` is never
      exceeded during Stage 2 sorting.
    * Each call creates and removes its own unique subdirectory — fully
      safe for concurrent parallel calls.
    * ``drop_count=True`` saves significant memory when processing many
      samples — e.g. 4M k-mers × 4 bytes = 16 MB saved per sample.
    """
    global _last_stats

    # ------------------------------------------------------------------
    # 1. Validate tmp_dir — fail early and clearly
    # ------------------------------------------------------------------
    if not isinstance(tmp_dir, str) or not tmp_dir.strip():
        raise TypeError(
            "tmp_dir is required and must be a non-empty string path.\n"
            "Provide a writable directory on a fast, large filesystem.\n"
            "On HPC: use /scratch or /localscratch, not /tmp.\n"
            "Example: count_kmers('file.fasta', tmp_dir='/scratch/kmc_tmp', ...)"
        )
    if not os.path.isdir(tmp_dir):
        raise ValueError(
            "tmp_dir '{}' does not exist. "
            "Create it first: os.makedirs('{}', exist_ok=True)".format(
                tmp_dir, tmp_dir)
        )
    if not os.access(tmp_dir, os.W_OK):
        raise ValueError(
            "tmp_dir '{}' is not writable by this process.".format(tmp_dir)
        )

    # ------------------------------------------------------------------
    # 2. Purge stale subdirs left by SIGKILL'd processes (best-effort)
    # ------------------------------------------------------------------
    _purge_stale_subdirs(tmp_dir)

    # ------------------------------------------------------------------
    # 3. Validate inputs
    # ------------------------------------------------------------------
    if isinstance(input_files, str):
        input_files = [input_files]

    for f in input_files:
        if not os.path.isfile(f):
            raise FileNotFoundError("Input file not found: {}".format(f))

    if k < 1:
        raise ValueError("k must be >= 1, got {}.".format(k))

    if threads <= 0:
        threads = os.cpu_count() or 4

    # Normalise format string
    fmt_key = (input_fmt.lower()
               .lstrip("-")
               .replace("fq",   "fastq")
               .replace(".gz",  "")
               .replace(".fna", "fasta")
               .replace(".faa", "fasta"))
    if fmt_key not in _INPUT_FMT:
        raise ValueError(
            "Unknown input_fmt '{}'. "
            "Choose from: fastq, fasta, mfasta, bam, kmc.\n"
            "Call kmcpy.help() for details.".format(input_fmt)
        )

    # Warn if "fasta" is used on a multi-line FASTA file
    if _INPUT_FMT.get(fmt_key, -1) == 1:
        try:
            with open(input_files[0], "r") as _fh:
                _in_seq = False
                _lines  = 0
                for _line in _fh:
                    if _line.startswith(">"):
                        if _in_seq and _lines > 1:
                            import warnings
                            warnings.warn(
                                "input_fmt='fasta' was specified but '{}' appears "
                                "to be multi-line FASTA. Automatically switching to "
                                "'mfasta'. Pass input_fmt='mfasta' to suppress "
                                "this warning.".format(input_files[0]),
                                UserWarning, stacklevel=2
                            )
                            fmt_key = "mfasta"
                            break
                        _in_seq = True
                        _lines  = 0
                    else:
                        _lines += 1
                        if _lines > 50:
                            break
        except Exception:
            pass

    # ------------------------------------------------------------------
    # 4. Create unique run subdirectory
    #    Registered globally so atexit/signal handlers can clean it up
    #    if the process exits before the finally: block runs.
    # ------------------------------------------------------------------
    subdir = _make_run_subdir(tmp_dir)
    _register_subdir(subdir)

    try:
        # --------------------------------------------------------------
        # 5. Run KMC
        #    - Stage 1 bins live in subdir if they spill from RAM
        #    - Stage 2 output comes back to RAM (packedKmers/prefixArray)
        # --------------------------------------------------------------
        raw = _ext._count_kmers_internal(
            input_files   = input_files,
            tmp_path      = subdir,          # KMC writes spill bins here
            kmer_len      = k,
            max_ram_gb    = max_ram_gb,
            n_threads     = threads,
            n_readers     = n_readers,
            n_splitters   = n_splitters,
            signature_len = signature_len,
            n_bins        = n_bins,
            canonical     = canonical,
            homopolymer   = homopolymer,
            strict_mem    = strict_mem,
            input_type    = _INPUT_FMT[fmt_key],
            cutoff_min    = min_count,
            cutoff_max    = max_count,
            counter_max   = counter_max,
            drop_count    = drop_count,      # Fork: skip counts_arr allocation in C++
        )

    finally:
        # --------------------------------------------------------------
        # 6. Always clean up the subdir — success, exception, or signal.
        #    ignore_errors=True so a partial cleanup never masks the
        #    original exception.
        # --------------------------------------------------------------
        _cleanup_subdir(subdir)
        _unregister_subdir(subdir)

    # ------------------------------------------------------------------
    # 7. Store run statistics
    # ------------------------------------------------------------------
    _last_stats = {
        "n_unique"    : int(raw["n_unique"]),
        "n_total"     : int(raw["n_total"]),
        "n_sequences" : int(raw["n_sequences"]),
        "below_cutoff": int(raw["below_cutoff"]),
        "above_cutoff": int(raw["above_cutoff"]),
        "stage1_time" : float(raw["stage1_time"]),
        "stage2_time" : float(raw["stage2_time"]),
        "kmer_len"    : int(raw["kmer_len"]),
        "drop_count"  : bool(drop_count),
    }

    n = raw["kmers"].shape[0]


    # ------------------------------------------------------------------
    # 8. decoded=False — space-optimal packed key columns.
    #
    # n_words = ceil(k * 2 / 64) — minimum uint64 words needed.
    # Column names: kmer_0 (MSB) ... kmer_{n_words-1} (LSB).
    #
    # Space efficiency:
    #   k ≤  32: n_words=1  →  1 column  (at most  2 bits wasted)
    #   k ≤  64: n_words=2  →  2 columns (at most  2 bits wasted)
    #   k ≤  96: n_words=3  →  3 columns
    #   k ≤ 128: n_words=4  →  4 columns
    #   k ≤ 192: n_words=6  →  6 columns
    #   k ≤ 256: n_words=8  →  8 columns (KMC maximum)
    #
    # Key layout per row:
    #   kmer_0          holds the most significant bits (base[0] side)
    #   kmer_{n_words-1} holds the least significant bits (base[k-1] side)
    #   base[k-1] in bits 1:0 of kmer_{n_words-1}
    #   Unused high bits of kmer_0 are always zero.
    #   A=00, C=01, G=10, T=11 (exact KMC 2-bit encoding, zero transformation)
    #
    # cuDF merge (works for all k):
    #   key_cols = [f"kmer_{i}" for i in range(n_words)]
    #   merged = gdf1.merge(gdf2, on=key_cols)
    # ------------------------------------------------------------------
    if not decoded:
        n_words    = int(raw["n_words"])
        key_words  = raw["key_words"]          # (n, n_words) uint64
        key_cols   = ["kmer_{}".format(i) for i in range(n_words)]

        if n == 0:
            col_dict = {col: pd.Series(dtype="uint64") for col in key_cols}
            if not drop_count:
                col_dict["count"] = pd.Series(dtype="uint32")
            return pd.DataFrame(col_dict)

        col_dict = {
            "kmer_{}".format(i): key_words[:, i].astype("uint64")
            for i in range(n_words)
        }
        if not drop_count:
            col_dict["count"] = raw["counts"].astype("uint32")
        return pd.DataFrame(col_dict)

    # ------------------------------------------------------------------
    # 9. decoded=True — map 2-bit codes → ACGT strings.
    # ------------------------------------------------------------------
    if n == 0:
        col_dict = {"kmer": pd.Series(dtype="string")}
        if not drop_count:
            col_dict["count"] = pd.Series(dtype="uint32")
        return pd.DataFrame(col_dict)

    _DECODE   = np.array([65, 67, 71, 84], dtype=np.uint8)   # A C G T
    ascii_arr = _DECODE[raw["kmers"]]                          # (n, k) uint8
    kmer_bytes = ascii_arr.view("S{}".format(k)).reshape(n)
    col_dict = {
        "kmer": pd.array(kmer_bytes.astype("U{}".format(k)), dtype="string"),
    }
    if not drop_count:
        col_dict["count"] = raw["counts"].astype("uint32")
    df = pd.DataFrame(col_dict)
    df.sort_values("kmer", ignore_index=True, inplace=True)
    return df


# ---------------------------------------------------------------------------
# last_run_stats()
# ---------------------------------------------------------------------------

def last_run_stats() -> dict:
    """
    Return statistics from the most recent count_kmers() call.

    Returns
    -------
    dict with keys:
        n_unique     : int   — unique k-mers returned (after cutoff filtering)
        n_total      : int   — total k-mer occurrences in the input
        n_sequences  : int   — number of sequences in the input
        below_cutoff : int   — k-mers discarded because count < min_count
        above_cutoff : int   — k-mers discarded because count > max_count
        stage1_time  : float — seconds in Stage 1 (splitting + binning)
        stage2_time  : float — seconds in Stage 2 (sorting + counting)
        kmer_len     : int   — k-mer length used in the run

    Examples
    --------
    >>> df = kmcpy.count_kmers("sample.fasta", k=15,
    ...                         input_fmt="mfasta",
    ...                         tmp_dir="/scratch/kmc_tmp")
    >>> stats = kmcpy.last_run_stats()
    >>> print(f"Counted {stats['n_unique']:,} unique k-mers "
    ...       f"in {stats['stage1_time'] + stats['stage2_time']:.2f}s")
    """
    return dict(_last_stats)


# ---------------------------------------------------------------------------
# kmc_version()
# ---------------------------------------------------------------------------

def kmc_version() -> str:
    """Return the KMC engine version string (e.g. '3.2.4')."""
    return __kmc_engine_version__


# ---------------------------------------------------------------------------
# Module exports
# ---------------------------------------------------------------------------

__all__ = [
    "count_kmers",
    "last_run_stats",
    "kmc_version",
    "help",
    "__version__",
    "__kmc_engine_version__",
]