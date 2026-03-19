#!/usr/bin/env python3
"""
diagnostics/check_install.py
=============================
kmcpy installation and runtime diagnostic suite.

Checks in order:
  1. Python environment and sys.path shadowing
  2. Package metadata (version, location)
  3. C extension import and symbol integrity
  4. KMC engine version
  5. Temp directory access
  6. Functional smoke test (tiny synthetic FASTA)
  7. decoded=True / decoded=False output shape
  8. drop_count behaviour
  9. min_count / max_count filters
 10. Parallel safety (two concurrent calls)

Usage:
    python diagnostics/check_install.py
    python diagnostics/check_install.py --verbose
    python diagnostics/check_install.py --tmp-dir /scratch/kmc_tmp
"""

import argparse
import os
import sys
import shutil
import tempfile
import textwrap
import traceback
import math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

# ---------------------------------------------------------------------------
# ANSI colours
# ---------------------------------------------------------------------------
_TTY = sys.stdout.isatty()
G = "\033[32m" if _TTY else ""   # green
R = "\033[31m" if _TTY else ""   # red
Y = "\033[33m" if _TTY else ""   # yellow
B = "\033[34m" if _TTY else ""   # blue (dim)
W = "\033[1m"  if _TTY else ""   # bold
D = "\033[2m"  if _TTY else ""   # dim
X = "\033[0m"  if _TTY else ""   # reset

PASS  = f"{G}PASS{X}"
FAIL  = f"{R}FAIL{X}"
WARN  = f"{Y}WARN{X}"
SKIP  = f"{D}SKIP{X}"

_results = []

def _header(title):
    width = 64
    print(f"\n{W}{'─' * width}{X}")
    print(f"{W}  {title}{X}")
    print(f"{W}{'─' * width}{X}")

def _check(name, status, detail="", verbose=False):
    tag = {"pass": PASS, "fail": FAIL, "warn": WARN, "skip": SKIP}[status]
    print(f"  {tag}  {name}")
    if detail and (verbose or status in ("fail", "warn")):
        for line in detail.strip().splitlines():
            print(f"        {D}{line}{X}")
    _results.append((name, status))

def _run_check(name, fn, verbose=False):
    try:
        result = fn()
        if result is None:
            _check(name, "pass", verbose=verbose)
        elif isinstance(result, tuple):
            status, detail = result
            _check(name, status, detail, verbose=verbose)
        else:
            _check(name, "pass", str(result), verbose=verbose)
    except Exception as exc:
        _check(name, "fail", traceback.format_exc(), verbose=verbose)


# ---------------------------------------------------------------------------
# Synthetic FASTA helper
# ---------------------------------------------------------------------------
def _write_fasta(path: Path, sequences: list):
    """Write a multi-line FASTA file."""
    with open(path, "w") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f">seq{i}\n")
            # wrap at 60 chars (multi-line FASTA)
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j+60] + "\n")


# Deterministic synthetic genome: 1000 bp, reproducible k-mer content
_SYNTHETIC = (
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "TTTAAACCCGGGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "GCGCGCGCGCGCGCGCGCGCGCGCGCATCATCATCATCATCATCATCATCAT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAAC"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTA"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAAC"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTA"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAAC"
)


# ---------------------------------------------------------------------------
# Check functions
# ---------------------------------------------------------------------------

def check_sys_path(verbose):
    """Detect local kmcpy/ folder shadowing the installed package."""
    cwd = os.path.abspath(os.getcwd())
    script_dir = os.path.dirname(os.path.abspath(__file__))
    shadowed = []
    for p in sys.path:
        ap = os.path.abspath(p) if p else cwd
        candidate = os.path.join(ap, "kmcpy")
        if os.path.isdir(candidate):
            init = os.path.join(candidate, "__init__.py")
            core_so = [f for f in os.listdir(candidate) if "_core" in f and f.endswith(".so")]
            shadowed.append((ap, bool(core_so)))
    if not shadowed:
        return "warn", "No kmcpy directory found anywhere on sys.path — is it installed?"
    installed = [(p, has_so) for p, has_so in shadowed if has_so]
    shadow    = [(p, has_so) for p, has_so in shadowed if not has_so]
    detail = ""
    if shadow:
        detail += "Source folders (no _core.so) on sys.path:\n"
        for p, _ in shadow:
            detail += f"  {p}/kmcpy/\n"
        detail += "These will shadow the installed package if they appear first in sys.path.\n"
        detail += "Run scripts from outside the project root to avoid this.\n"
        return "warn", detail
    return "pass", f"Installed package at: {installed[0][0]}/kmcpy/"


def check_import(verbose):
    import kmcpy  # noqa
    loc = Path(kmcpy.__file__).parent
    return "pass", f"kmcpy {kmcpy.__version__} at {loc}"


def check_core_extension(verbose):
    import kmcpy._core as _core  # noqa
    return "pass", f"_core extension loaded from {_core.__file__}"


def check_kmc_version(verbose):
    import kmcpy
    v = kmcpy.kmc_version()
    if v != "3.2.4":
        return "warn", f"Expected KMC 3.2.4, got {v}"
    return "pass", f"KMC engine {v}"


def check_tmp_dir(tmp_dir, verbose):
    if not os.path.isdir(tmp_dir):
        return "fail", f"tmp_dir does not exist: {tmp_dir}"
    if not os.access(tmp_dir, os.W_OK):
        return "fail", f"tmp_dir is not writable: {tmp_dir}"
    # write test
    probe = os.path.join(tmp_dir, ".kmcpy_probe")
    try:
        with open(probe, "w") as fh:
            fh.write("probe")
        os.remove(probe)
    except Exception as e:
        return "fail", f"Cannot write to tmp_dir: {e}"
    # check free space
    stat = shutil.disk_usage(tmp_dir)
    free_gb = stat.free / 1e9
    detail = f"{tmp_dir}  ({free_gb:.1f} GB free)"
    if free_gb < 1.0:
        return "warn", detail + " — less than 1 GB free"
    return "pass", detail


def check_smoke(tmp_dir, fasta_path, verbose):
    import kmcpy
    df = kmcpy.count_kmers(
        str(fasta_path),
        k         = 15,
        input_fmt = "mfasta",
        tmp_dir   = tmp_dir,
        min_count = 1,
    )
    if df is None or len(df) == 0:
        return "fail", "count_kmers() returned empty DataFrame"
    if list(df.columns) != ["kmer", "count"]:
        return "fail", f"Unexpected columns: {list(df.columns)}"
    if df["kmer"].str.len().max() != 15:
        return "fail", "k-mer length mismatch"
    n = len(df)
    detail = f"{n:,} unique 15-mers  |  columns: {list(df.columns)}  |  dtypes: kmer={df['kmer'].dtype} count={df['count'].dtype}"
    return "pass", detail


def check_decoded_false(tmp_dir, fasta_path, k, verbose):
    import kmcpy
    import math
    df = kmcpy.count_kmers(
        str(fasta_path),
        k         = k,
        input_fmt = "mfasta",
        tmp_dir   = tmp_dir,
        min_count = 1,
        decoded   = False,
    )
    n_words_expected = math.ceil(k * 2 / 64)
    key_cols = [f"kmer_{i}" for i in range(n_words_expected)]
    expected_cols = key_cols + ["count"]
    if list(df.columns) != expected_cols:
        return "fail", f"Expected {expected_cols}, got {list(df.columns)}"
    if df[key_cols[0]].dtype != "uint64":
        return "fail", f"Expected uint64, got {df[key_cols[0]].dtype}"
    detail = (f"k={k}  n_words={n_words_expected}  "
              f"{len(df):,} k-mers  columns={list(df.columns)}")
    return "pass", detail


def check_drop_count(tmp_dir, fasta_path, verbose):
    import kmcpy
    df_with    = kmcpy.count_kmers(str(fasta_path), k=15, input_fmt="mfasta",
                                    tmp_dir=tmp_dir, min_count=1, drop_count=False)
    df_without = kmcpy.count_kmers(str(fasta_path), k=15, input_fmt="mfasta",
                                    tmp_dir=tmp_dir, min_count=1, drop_count=True)
    if "count" in df_without.columns:
        return "fail", "count column present despite drop_count=True"
    if "count" not in df_with.columns:
        return "fail", "count column missing with drop_count=False"
    if len(df_with) != len(df_without):
        return "fail", f"Row count mismatch: {len(df_with)} vs {len(df_without)}"
    detail = (f"drop_count=False: {list(df_with.columns)}  "
              f"drop_count=True: {list(df_without.columns)}  "
              f"rows match: {len(df_with):,}")
    return "pass", detail


def check_min_count_filter(tmp_dir, fasta_path, verbose):
    import kmcpy
    df_all  = kmcpy.count_kmers(str(fasta_path), k=15, input_fmt="mfasta",
                                  tmp_dir=tmp_dir, min_count=1)
    df_filt = kmcpy.count_kmers(str(fasta_path), k=15, input_fmt="mfasta",
                                  tmp_dir=tmp_dir, min_count=2)
    if len(df_filt) > len(df_all):
        return "fail", "min_count=2 returned MORE k-mers than min_count=1"
    if len(df_filt) > 0 and df_filt["count"].min() < 2:
        return "fail", f"min_count=2 filter violated: min count = {df_filt['count'].min()}"
    detail = (f"min_count=1: {len(df_all):,} k-mers  "
              f"min_count=2: {len(df_filt):,} k-mers  "
              f"singletons removed: {len(df_all)-len(df_filt):,}")
    return "pass", detail


def check_max_count_filter(tmp_dir, fasta_path, verbose):
    import kmcpy
    df_all  = kmcpy.count_kmers(str(fasta_path), k=15, input_fmt="mfasta",
                                  tmp_dir=tmp_dir, min_count=1)
    max_c   = int(df_all["count"].median())
    df_filt = kmcpy.count_kmers(str(fasta_path), k=15, input_fmt="mfasta",
                                  tmp_dir=tmp_dir, min_count=1, max_count=max_c)
    if len(df_filt) > 0 and df_filt["count"].max() > max_c:
        return "fail", f"max_count={max_c} filter violated: max count = {df_filt['count'].max()}"
    detail = (f"max_count={max_c}: {len(df_filt):,} k-mers  "
              f"(from {len(df_all):,} total)")
    return "pass", detail


def check_stats(tmp_dir, fasta_path, verbose):
    import kmcpy
    kmcpy.count_kmers(str(fasta_path), k=15, input_fmt="mfasta",
                       tmp_dir=tmp_dir, min_count=1)
    s = kmcpy.last_run_stats()
    required = {"n_unique", "n_total", "n_sequences", "stage1_time", "stage2_time", "kmer_len"}
    missing = required - set(s.keys())
    if missing:
        return "fail", f"Missing stats keys: {missing}"
    detail = (f"n_unique={s['n_unique']:,}  n_total={s['n_total']:,}  "
              f"stage1={s['stage1_time']:.3f}s  stage2={s['stage2_time']:.3f}s")
    return "pass", detail


def _parallel_worker(args):
    """Worker for parallel safety test — must be top-level for pickling."""
    fasta_path, tmp_dir, worker_id = args
    import kmcpy
    df = kmcpy.count_kmers(
        fasta_path,
        k         = 15,
        input_fmt = "mfasta",
        tmp_dir   = tmp_dir,
        min_count = 1,
    )
    return len(df)


def check_parallel(tmp_dir, fasta_path, verbose):
    args = [(str(fasta_path), tmp_dir, i) for i in range(2)]
    with ProcessPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(_parallel_worker, args))
    if len(set(results)) != 1:
        return "fail", f"Parallel calls returned different row counts: {results}"
    detail = f"2 concurrent calls both returned {results[0]:,} k-mers"
    return "pass", detail


def check_sorted_output(tmp_dir, fasta_path, verbose):
    import kmcpy
    df = kmcpy.count_kmers(str(fasta_path), k=15, input_fmt="mfasta",
                             tmp_dir=tmp_dir, min_count=1)
    kmers = df["kmer"].tolist()
    if kmers != sorted(kmers):
        return "fail", "decoded=True output is not sorted lexicographically"
    return "pass", "Output is lexicographically sorted"


def check_large_k(tmp_dir, fasta_path, verbose):
    """Test k=63 (2 uint64 words) and k=128 (4 uint64 words)."""
    import kmcpy
    import math
    issues = []
    for k in [31, 63, 100]:
        if len(_SYNTHETIC) < k:
            continue
        try:
            df = kmcpy.count_kmers(str(fasta_path), k=k, input_fmt="mfasta",
                                    tmp_dir=tmp_dir, min_count=1, decoded=False)
            n_words = math.ceil(k * 2 / 64)
            expected = [f"kmer_{i}" for i in range(n_words)] + ["count"]
            if list(df.columns) != expected:
                issues.append(f"k={k}: wrong columns {list(df.columns)}")
        except Exception as e:
            issues.append(f"k={k}: {e}")
    if issues:
        return "fail", "\n".join(issues)
    return "pass", "k=31, k=63, k=100 all produce correct column layouts"


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
def _summary():
    total  = len(_results)
    passed = sum(1 for _, s in _results if s == "pass")
    warned = sum(1 for _, s in _results if s == "warn")
    failed = sum(1 for _, s in _results if s == "fail")
    skipped= sum(1 for _, s in _results if s == "skip")

    print(f"\n{'─'*64}")
    print(f"{W}  Summary{X}")
    print(f"{'─'*64}")
    print(f"  {G}{passed} passed{X}  "
          f"{Y}{warned} warnings{X}  "
          f"{R}{failed} failed{X}  "
          f"{D}{skipped} skipped{X}  "
          f"/ {total} total")
    print(f"{'─'*64}\n")

    if failed:
        print(f"{R}  Some checks failed. See details above.{X}\n")
        return 1
    if warned:
        print(f"{Y}  All checks passed with warnings.{X}\n")
        return 0
    print(f"{G}  All checks passed.{X}\n")
    return 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--tmp-dir", default=tempfile.gettempdir(),
        help="Scratch directory for KMC temp files (default: system temp)",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Show detail for passing checks too",
    )
    args = parser.parse_args()
    tmp_dir = args.tmp_dir
    verbose = args.verbose

    print(f"\n{W}kmcpy diagnostic suite{X}")
    print(f"{D}Run from outside the project source tree for accurate results.{X}")

    # ── 1. Environment ──────────────────────────────────────────────────────
    _header("1. Python environment")
    _run_check("Python version",
               lambda: ("pass", f"{sys.version}"), verbose)
    _run_check("sys.path shadowing",
               lambda: check_sys_path(verbose), verbose)

    # ── 2. Package ──────────────────────────────────────────────────────────
    _header("2. Package installation")
    _run_check("import kmcpy",         lambda: check_import(verbose), verbose)
    _run_check("_core C extension",    lambda: check_core_extension(verbose), verbose)
    _run_check("KMC engine version",   lambda: check_kmc_version(verbose), verbose)

    # ── 3. Temp directory ───────────────────────────────────────────────────
    _header("3. Temp directory")
    _run_check(f"tmp_dir access  [{tmp_dir}]",
               lambda: check_tmp_dir(tmp_dir, verbose), verbose)

    # ── 4. Functional tests ─────────────────────────────────────────────────
    _header("4. Functional tests")

    # Write synthetic FASTA
    with tempfile.TemporaryDirectory() as td:
        fasta = Path(td) / "synthetic.fasta"
        _write_fasta(fasta, [_SYNTHETIC, _SYNTHETIC[:500], _SYNTHETIC[200:700]])

        _run_check("Smoke test  (decoded=True, k=15)",
                   lambda: check_smoke(tmp_dir, fasta, verbose), verbose)
        _run_check("decoded=False  k=15  (1 uint64 word)",
                   lambda: check_decoded_false(tmp_dir, fasta, 15, verbose), verbose)
        _run_check("decoded=False  k=33  (2 uint64 words)",
                   lambda: check_decoded_false(tmp_dir, fasta, 33, verbose), verbose)
        _run_check("Large k  (k=31, k=63, k=100)",
                   lambda: check_large_k(tmp_dir, fasta, verbose), verbose)
        _run_check("drop_count=True  (C++ boundary saving)",
                   lambda: check_drop_count(tmp_dir, fasta, verbose), verbose)
        _run_check("min_count filter",
                   lambda: check_min_count_filter(tmp_dir, fasta, verbose), verbose)
        _run_check("max_count filter",
                   lambda: check_max_count_filter(tmp_dir, fasta, verbose), verbose)
        _run_check("last_run_stats()",
                   lambda: check_stats(tmp_dir, fasta, verbose), verbose)
        _run_check("Output lexicographic sort",
                   lambda: check_sorted_output(tmp_dir, fasta, verbose), verbose)

    # ── 5. Parallel safety ──────────────────────────────────────────────────
    _header("5. Parallel safety")
    with tempfile.TemporaryDirectory() as td:
        fasta = Path(td) / "synthetic.fasta"
        _write_fasta(fasta, [_SYNTHETIC])
        _run_check("2 concurrent ProcessPoolExecutor workers",
                   lambda: check_parallel(tmp_dir, fasta, verbose), verbose)

    return _summary()


if __name__ == "__main__":
    sys.exit(main())