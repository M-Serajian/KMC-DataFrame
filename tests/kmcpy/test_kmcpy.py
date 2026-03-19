#!/usr/bin/env python3
"""
tests/test_kmcpy.py
===================
Regression test suite for kmcpy.

Run from outside the project source tree:
    python tests/test_kmcpy.py
    python tests/test_kmcpy.py --tmp-dir /scratch/kmc_tmp

Or via pytest (from outside the project root):
    pytest tests/test_kmcpy.py
"""

import argparse
import math
import os
import sys
import tempfile
from pathlib import Path

try:
    import pytest
except ImportError:
    print("pytest is required to run the test suite.")
    print("Install it with:  pip install pytest")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------
_GENOME = (
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
) * 3


def _write_fasta(path, sequences):
    with open(path, "w") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f">seq{i}\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j+60] + "\n")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def tmp_dir():
    d = str(Path(__file__).parent / "tmp")
    os.makedirs(d, exist_ok=True)
    return d


@pytest.fixture(scope="session")
def fasta(tmp_dir):
    p = Path(tmp_dir) / "test_input.fasta"
    _write_fasta(p, [_GENOME, _GENOME[:600], _GENOME[200:800]])
    return str(p)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestImport:
    def test_import(self):
        import kmcpy
        assert kmcpy.__version__

    def test_kmc_version(self):
        import kmcpy
        assert kmcpy.kmc_version() == "3.2.4"

    def test_core_extension(self):
        import kmcpy._core
        assert hasattr(kmcpy._core, "_count_kmers_internal")


class TestBasicCounting:
    def test_returns_dataframe(self, fasta, tmp_dir):
        import kmcpy, pandas as pd
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1)
        assert isinstance(df, pd.DataFrame)

    def test_columns_decoded_true(self, fasta, tmp_dir):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1)
        assert list(df.columns) == ["kmer", "count"]

    def test_kmer_length(self, fasta, tmp_dir):
        import kmcpy
        for k in [15, 21, 31]:
            df = kmcpy.count_kmers(fasta, k=k, input_fmt="mfasta",
                                    tmp_dir=tmp_dir, min_count=1)
            assert df["kmer"].str.len().max() == k, f"k={k} length mismatch"

    def test_sorted_output(self, fasta, tmp_dir):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1)
        kmers = df["kmer"].tolist()
        assert kmers == sorted(kmers), "Output is not lexicographically sorted"

    def test_valid_bases(self, fasta, tmp_dir):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1)
        invalid = df["kmer"].str.contains("[^ACGT]").any()
        assert not invalid, "Non-ACGT characters found in k-mers"

    def test_counts_positive(self, fasta, tmp_dir):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1)
        assert (df["count"] > 0).all()

    def test_dtypes(self, fasta, tmp_dir):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1)
        assert str(df["count"].dtype) == "uint32"


class TestDecodedFalse:
    @pytest.mark.parametrize("k,expected_words", [
        (15, 1), (31, 1), (32, 1), (33, 2), (63, 2), (64, 2), (65, 3),
    ])
    def test_n_words(self, fasta, tmp_dir, k, expected_words):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=k, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1, decoded=False)
        key_cols = [c for c in df.columns if c.startswith("kmer_")]
        assert len(key_cols) == expected_words, \
            f"k={k}: expected {expected_words} words, got {len(key_cols)}"

    def test_uint64_dtype(self, fasta, tmp_dir):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1, decoded=False)
        assert str(df["kmer_0"].dtype) == "uint64"

    def test_same_row_count_as_decoded(self, fasta, tmp_dir):
        import kmcpy
        df_dec = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                    tmp_dir=tmp_dir, min_count=1, decoded=True)
        df_raw = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                    tmp_dir=tmp_dir, min_count=1, decoded=False)
        assert len(df_dec) == len(df_raw)


class TestFilters:
    def test_min_count(self, fasta, tmp_dir):
        import kmcpy
        df1 = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                  tmp_dir=tmp_dir, min_count=1)
        df2 = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                  tmp_dir=tmp_dir, min_count=2)
        assert len(df2) <= len(df1)
        if len(df2) > 0:
            assert df2["count"].min() >= 2

    def test_max_count(self, fasta, tmp_dir):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                 tmp_dir=tmp_dir, min_count=1)
        max_c = int(df["count"].median())
        df_filt = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                     tmp_dir=tmp_dir, min_count=1,
                                     max_count=max_c)
        if len(df_filt) > 0:
            assert df_filt["count"].max() <= max_c

    def test_min_count_removes_singletons(self, fasta, tmp_dir):
        import kmcpy
        df1 = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                  tmp_dir=tmp_dir, min_count=1)
        df2 = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                  tmp_dir=tmp_dir, min_count=2)
        singletons = (df1["count"] == 1).sum()
        assert len(df1) - len(df2) == singletons


class TestDropCount:
    def test_no_count_column(self, fasta, tmp_dir):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1, drop_count=True)
        assert "count" not in df.columns

    def test_same_rows_with_and_without(self, fasta, tmp_dir):
        import kmcpy
        df_with    = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                        tmp_dir=tmp_dir, min_count=1,
                                        drop_count=False)
        df_without = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                        tmp_dir=tmp_dir, min_count=1,
                                        drop_count=True)
        assert len(df_with) == len(df_without)

    def test_decoded_false_drop_count(self, fasta, tmp_dir):
        import kmcpy
        df = kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                                tmp_dir=tmp_dir, min_count=1,
                                decoded=False, drop_count=True)
        assert "count" not in df.columns
        assert "kmer_0" in df.columns


class TestMultipleInputFiles:
    def test_list_input(self, tmp_dir):
        import kmcpy
        f1 = Path(tmp_dir) / "input1.fasta"
        f2 = Path(tmp_dir) / "input2.fasta"
        _write_fasta(f1, [_GENOME])
        _write_fasta(f2, [_GENOME[:600]])
        df = kmcpy.count_kmers([str(f1), str(f2)], k=15,
                                 input_fmt="mfasta", tmp_dir=tmp_dir,
                                 min_count=1)
        assert len(df) > 0


class TestStats:
    def test_last_run_stats_keys(self, fasta, tmp_dir):
        import kmcpy
        kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                           tmp_dir=tmp_dir, min_count=1)
        s = kmcpy.last_run_stats()
        required = {"n_unique", "n_total", "n_sequences",
                    "stage1_time", "stage2_time", "kmer_len"}
        assert required.issubset(set(s.keys()))

    def test_stats_kmer_len(self, fasta, tmp_dir):
        import kmcpy
        kmcpy.count_kmers(fasta, k=21, input_fmt="mfasta",
                           tmp_dir=tmp_dir, min_count=1)
        assert kmcpy.last_run_stats()["kmer_len"] == 21

    def test_n_unique_positive(self, fasta, tmp_dir):
        import kmcpy
        kmcpy.count_kmers(fasta, k=15, input_fmt="mfasta",
                           tmp_dir=tmp_dir, min_count=1)
        assert kmcpy.last_run_stats()["n_unique"] > 0


# ---------------------------------------------------------------------------
# Standalone runner (without pytest)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tmp-dir",
        default=str(Path(__file__).parent / "tmp"),
        help="Scratch directory for KMC temp files")
    args = parser.parse_args()

    sys.exit(pytest.main([__file__, "-v",
        f"--tmp-dir={args.tmp_dir}", "--tb=short"]))