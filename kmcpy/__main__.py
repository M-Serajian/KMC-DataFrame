"""
kmcpy.__main__
--------------
Command-line interface for kmcpy.

Usage
-----
    conda activate KMC-env
    python -m kmcpy --help
    python -m kmcpy reads.fastq --k 27
    python -m kmcpy assembly.fasta --fmt fasta --k 15 --min-count 1
    python -m kmcpy r1.fq r2.fq --k 31 --threads 8 --out kmers.tsv
"""

import argparse
import sys
import os
import textwrap

# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog        = "python -m kmcpy",
        description = textwrap.dedent("""\
            kmcpy — KMC k-mer counting, fully in-memory.
            No intermediate files written to disk.
            Equivalent to the KMC CLI but returns results to stdout or a file.
        """),
        formatter_class = argparse.RawDescriptionHelpFormatter,
        epilog = textwrap.dedent("""\
            KMC flag equivalents
            --------------------
              -k<len>   →  --k
              -fq/-fa   →  --fmt fastq/fasta
              -m<GB>    →  --max-ram
              -t<n>     →  --threads
              -ci<n>    →  --min-count
              -cx<n>    →  --max-count
              -cs<n>    →  --counter-max
              -b        →  --non-canonical
              -hc       →  --homopolymer
              -sm       →  --strict-mem
              -p<n>     →  --signature-len

            Examples
            --------
              # FASTQ, k=27, KMC defaults
              python -m kmcpy reads.fastq --k 27

              # FASTA assembly, all k-mers (singletons included)
              python -m kmcpy assembly.fasta --fmt fasta --k 15 --min-count 1

              # Paired-end FASTQ, save to TSV
              python -m kmcpy r1.fastq r2.fastq --k 31 --threads 8 --out kmers.tsv

              # Full options
              python -m kmcpy sample.fasta --fmt fasta --k 27 --max-ram 24 \\
                     --threads 8 --min-count 2 --max-count 500 --out out.tsv
        """),
    )

    # Input
    parser.add_argument(
        "input_files", nargs="+", metavar="FILE",
        help="Input file(s). All must be the same format."
    )

    # k-mer options
    g_kmer = parser.add_argument_group("k-mer options")
    g_kmer.add_argument("--k", type=int, default=25, metavar="LEN",
        help="K-mer length (default: 25)  [-k]")
    g_kmer.add_argument("--fmt", default="fastq", metavar="FORMAT",
        choices=["fastq","fq","fasta","fa","mfasta","fm","bam","kmc"],
        help="Input format: fastq/fasta/mfasta/bam/kmc (default: fastq)  [-f]")
    g_kmer.add_argument("--min-count", type=int, default=2, metavar="N",
        help="Exclude k-mers with count < N (default: 2)  [-ci]")
    g_kmer.add_argument("--max-count", type=int, default=1_000_000_000, metavar="N",
        help="Exclude k-mers with count > N (default: 1e9)  [-cx]")
    g_kmer.add_argument("--counter-max", type=int, default=255, metavar="N",
        help="Saturate counter at N (default: 255)  [-cs]")
    g_kmer.add_argument("--non-canonical", action="store_true",
        help="Count non-canonical k-mers  [-b]")
    g_kmer.add_argument("--homopolymer", action="store_true",
        help="Homopolymer-compressed k-mers (experimental)  [-hc]")

    # Performance options
    g_perf = parser.add_argument_group("performance options")
    g_perf.add_argument("--max-ram", type=int, default=12, metavar="GB",
        help="RAM budget in GB (default: 12)  [-m]")
    g_perf.add_argument("--threads", type=int, default=0, metavar="N",
        help="Total threads (default: 0 = all cores)  [-t]")
    g_perf.add_argument("--strict-mem", action="store_true",
        help="Strict memory mode — never exceed --max-ram  [-sm]")
    g_perf.add_argument("--signature-len", type=int, default=9, metavar="N",
        help="Signature length 5-11 (default: 9)  [-p]")
    g_perf.add_argument("--n-bins", type=int, default=0, metavar="N",
        help="Number of bins (default: 0 = auto)  [-n]")
    g_perf.add_argument("--n-readers", type=int, default=0, metavar="N",
        help="Reader threads (default: 0 = auto)  [-sf]")
    g_perf.add_argument("--n-splitters", type=int, default=0, metavar="N",
        help="Splitter threads (default: 0 = auto)  [-sp]")

    # Output options
    g_out = parser.add_argument_group("output options")
    g_out.add_argument("--out", default=None, metavar="FILE",
        help="Write TSV output to FILE instead of stdout")
    g_out.add_argument("--no-header", action="store_true",
        help="Omit the header line from TSV output")
    g_out.add_argument("--stats", action="store_true",
        help="Print run statistics to stderr after counting")
    g_out.add_argument("--version", action="store_true",
        help="Print kmcpy and KMC engine versions and exit")

    return parser


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import kmcpy

    parser = build_parser()
    args   = parser.parse_args()

    if args.version:
        print("kmcpy {}".format(kmcpy.__version__))
        print("KMC engine {}".format(kmcpy.kmc_version()))
        sys.exit(0)

    # Run
    try:
        df = kmcpy.count_kmers(
            input_files   = args.input_files,
            k             = args.k,
            input_fmt     = args.fmt,
            max_ram_gb    = args.max_ram,
            threads       = args.threads,
            canonical     = not args.non_canonical,
            homopolymer   = args.homopolymer,
            signature_len = args.signature_len,
            n_bins        = args.n_bins,
            n_readers     = args.n_readers,
            n_splitters   = args.n_splitters,
            min_count     = args.min_count,
            max_count     = args.max_count,
            counter_max   = args.counter_max,
            strict_mem    = args.strict_mem,
        )
    except (FileNotFoundError, ValueError) as e:
        print("Error: {}".format(e), file=sys.stderr)
        sys.exit(1)

    # Output
    sep    = "\t"
    header = not args.no_header

    if args.out:
        df.to_csv(args.out, sep=sep, index=False, header=header)
        print("Wrote {:,} k-mers to {}".format(len(df), args.out),
              file=sys.stderr)
    else:
        df.to_csv(sys.stdout, sep=sep, index=False, header=header)

    # Stats
    if args.stats:
        s = kmcpy.last_run_stats()
        print("\n--- run statistics ---", file=sys.stderr)
        for k, v in s.items():
            print("  {:15s} {}".format(k, v), file=sys.stderr)


if __name__ == "__main__":
    main()