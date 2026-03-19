#!/usr/bin/env python3
"""
tests/kmcpy/sample_usage.py
============================
End-to-end usage examples for kmcpy.

This script is self-contained -- it uses the sample FASTA file committed
alongside it (tests/kmcpy/sample.fasta) and writes all KMC temp files to
tests/kmcpy/tmp/ which is pre-created in the repository.
No external data or configuration is required.

Demonstrates:
  1. Basic k-mer counting          (decoded=True, ACGT strings)
  2. Filtering                     (min_count, max_count)
  3. Packed integer keys           (decoded=False, for GPU merge)
  4. Count-free mode               (drop_count=True, maximum memory efficiency)
  5. Multiple input files          (list of FASTA paths)

Run from the project root:
    python tests/kmcpy/sample_usage.py

Or from inside tests/kmcpy/:
    python sample_usage.py

IMPORTANT -- pip install users:
    Do NOT run this script from the project root with:
        cd KMC-DataFrame && python -c "import kmcpy"
    Python will find the local kmcpy/ source folder before the installed
    package and raise ImportError. Always run from outside the project root
    or from inside tests/kmcpy/ as shown above.
"""

import math
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Resolve paths relative to this script's location.
# This guarantees the correct paths regardless of where the script is invoked.
# ---------------------------------------------------------------------------
SCRIPT_DIR  = Path(__file__).resolve().parent          # tests/kmcpy/
SAMPLE_FASTA = SCRIPT_DIR / "sample.fasta"             # committed sample input
TMP_DIR      = SCRIPT_DIR / "tmp"                      # pre-created scratch dir

# ---------------------------------------------------------------------------
# Guard: detect sys.path shadowing before importing kmcpy.
# If the project root is on sys.path AND the local kmcpy/ folder exists there,
# Python will import from source (no _core.so) instead of the installed package.
# ---------------------------------------------------------------------------
def _check_shadowing():
    project_root = SCRIPT_DIR.parent.parent             # KMC-DataFrame/
    local_kmcpy  = project_root / "kmcpy" / "__init__.py"
    if local_kmcpy.exists():
        for p in sys.path:
            if Path(p).resolve() == project_root:
                print(
                    "\n"
                    "WARNING: The project root is on sys.path and the local\n"
                    "  kmcpy/ source folder may shadow the installed package.\n"
                    "  If you get ImportError, run this script from outside\n"
                    "  the project root:\n"
                    "\n"
                    "    cd tests/kmcpy && python sample_usage.py\n"
                    "  or\n"
                    "    cd ~ && python /path/to/KMC-DataFrame/tests/kmcpy/sample_usage.py\n"
                )
                break

_check_shadowing()

# ---------------------------------------------------------------------------
# Import kmcpy -- must be installed via pip install .
# ---------------------------------------------------------------------------
try:
    import kmcpy
except ImportError as e:
    print(
        f"\nERROR: {e}\n"
        "\nkmcpy is not installed. Install it first:\n"
        "  pip install scikit-build-core pybind11 cmake ninja\n"
        "  cd KMC-DataFrame && pip install .\n"
    )
    sys.exit(1)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def section(title: str) -> None:
    width = 64
    print(f"\n{'─' * width}")
    print(f"  {title}")
    print(f"{'─' * width}")


def check_prerequisites():
    """Verify sample.fasta exists and tmp/ is writable."""
    if not SAMPLE_FASTA.exists():
        print(f"\nERROR: Sample FASTA not found: {SAMPLE_FASTA}")
        print("Make sure you cloned the full repository including tests/kmcpy/sample.fasta")
        sys.exit(1)

    TMP_DIR.mkdir(parents=True, exist_ok=True)
    probe = TMP_DIR / ".write_test"
    try:
        probe.write_text("ok")
        probe.unlink()
    except OSError as e:
        print(f"\nERROR: tmp directory is not writable: {TMP_DIR}\n{e}")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    check_prerequisites()

    print(f"\n{'═' * 64}")
    print(f"  kmcpy usage examples")
    print(f"{'═' * 64}")
    print(f"  kmcpy version  : {kmcpy.__version__}")
    print(f"  KMC engine     : {kmcpy.kmc_version()}")
    print(f"  Input FASTA    : {SAMPLE_FASTA}")
    print(f"  Scratch dir    : {TMP_DIR}")
    print(f"  Input size     : {SAMPLE_FASTA.stat().st_size:,} bytes")

    # ── 1. Basic counting (decoded=True) ─────────────────────────────────────
    section("1. Basic k-mer counting  (decoded=True, k=15)")

    df = kmcpy.count_kmers(
        str(SAMPLE_FASTA),
        k         = 15,
        input_fmt = "mfasta",   # multi-line FASTA
        tmp_dir   = str(TMP_DIR),
        min_count = 1,          # keep all k-mers including singletons
    )

    print(f"  Unique 15-mers  : {len(df):,}")
    print(f"  Columns         : {list(df.columns)}")
    print(f"  dtypes          : kmer={df['kmer'].dtype}  count={df['count'].dtype}")
    print(f"\n  Top 5 most frequent 15-mers:")
    print(df.nlargest(5, "count").to_string(index=False))

    stats = kmcpy.last_run_stats()
    print(f"\n  Run statistics:")
    print(f"    n_unique    = {stats['n_unique']:,}")
    print(f"    n_total     = {stats['n_total']:,}")
    print(f"    n_sequences = {stats['n_sequences']}")
    print(f"    stage1      = {stats['stage1_time']:.4f}s")
    print(f"    stage2      = {stats['stage2_time']:.4f}s")

    # ── 2. Filtering ──────────────────────────────────────────────────────────
    section("2. Filtering  (min_count / max_count)")

    df_filt = kmcpy.count_kmers(
        str(SAMPLE_FASTA),
        k         = 15,
        input_fmt = "mfasta",
        tmp_dir   = str(TMP_DIR),
        min_count = 2,    # discard singletons
        max_count = 50,   # discard highly repetitive k-mers
    )

    singletons = len(df) - len(df_filt)
    print(f"  min_count=1            : {len(df):,} k-mers")
    print(f"  min_count=2, max=50    : {len(df_filt):,} k-mers")
    print(f"  singletons removed     : {singletons:,}")
    if len(df_filt) > 0:
        print(f"  count range            : [{df_filt['count'].min()}, "
              f"{df_filt['count'].max()}]")

    # ── 3. Packed integer keys (decoded=False) ────────────────────────────────
    section("3. Packed integer keys  (decoded=False)  — GPU-ready")

    df_packed = kmcpy.count_kmers(
        str(SAMPLE_FASTA),
        k         = 15,
        input_fmt = "mfasta",
        tmp_dir   = str(TMP_DIR),
        min_count = 1,
        decoded   = False,    # uint64 columns instead of ACGT strings
    )

    n_words  = math.ceil(15 * 2 / 64)
    key_cols = [f"kmer_{i}" for i in range(n_words)]
    print(f"  k=15 → n_words={n_words}  ({n_words * 8} bytes per k-mer key)")
    print(f"  Columns  : {list(df_packed.columns)}")
    print(f"  dtypes   : {dict(df_packed.dtypes)}")
    print(f"\n  First 3 rows:")
    print(df_packed.head(3).to_string(index=False))
    print(f"\n  cuDF merge pattern (same syntax for all k):")
    print(f"    n_words  = {n_words}")
    print(f"    key_cols = {key_cols}")
    print(f"    merged   = gdf1.merge(gdf2, on=key_cols)")

    # ── 4. drop_count — maximum memory efficiency ─────────────────────────────
    section("4. drop_count=True  — presence only, no count column")

    df_keys = kmcpy.count_kmers(
        str(SAMPLE_FASTA),
        k          = 15,
        input_fmt  = "mfasta",
        tmp_dir    = str(TMP_DIR),
        min_count  = 1,
        decoded    = False,
        drop_count = True,    # counts never allocated at C++ level
    )

    mem_keys  = df_keys.memory_usage(deep=True).sum()
    mem_saved = len(df_keys) * 4   # uint32 per k-mer never allocated
    print(f"  Columns          : {list(df_keys.columns)}")
    print(f"  Memory (keys)    : {mem_keys / 1e6:.3f} MB")
    print(f"  Memory saved     : {mem_saved / 1e6:.3f} MB  "
          f"(counts never allocated in C++)")
    print(f"\n  Use case: set intersection, Jaccard similarity, SBWT construction")

    # ── 5. Multiple input files ───────────────────────────────────────────────
    section("5. Multiple input files  — one combined counting pass")

    # Use the same file twice to simulate two samples
    df_multi = kmcpy.count_kmers(
        [str(SAMPLE_FASTA), str(SAMPLE_FASTA)],
        k         = 15,
        input_fmt = "mfasta",
        tmp_dir   = str(TMP_DIR),
        min_count = 1,
    )

    print(f"  Single file  : {len(df):,} unique 15-mers")
    print(f"  Two files    : {len(df_multi):,} unique 15-mers  "
          f"(same k-mers, higher counts)")
    print(f"  All input files are counted together in a single KMC pass.")

    print(f"\n{'═' * 64}")
    print(f"  All examples completed successfully.")
    print(f"{'═' * 64}\n")


if __name__ == "__main__":
    main()