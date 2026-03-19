<div align="center">

[![License: GPL v3](https://img.shields.io/badge/License-GPL--3.0-bd0000?style=flat-square)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-3776ab?style=flat-square&logo=python&logoColor=white)](https://www.python.org/downloads/)
[![KMC Engine](https://img.shields.io/badge/KMC%20Engine-3.2.4-1d4ed8?style=flat-square&logo=cplusplus&logoColor=white)](https://github.com/refresh-bio/KMC/releases/tag/v3.2.4)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-64748b?style=flat-square&logo=linux&logoColor=white)](https://github.com/M-Serajian/KMC-DataFrame#installation)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtx304-16a34a?style=flat-square)](https://doi.org/10.1093/bioinformatics/btx304)
[![Issues](https://img.shields.io/badge/Issues-Welcome-f59e0b?style=flat-square&logo=github&logoColor=white)](https://github.com/M-Serajian/KMC-DataFrame/issues)

<br/><br/>

# kmcpy

### High-performance k-mer counting for Python — fully in-memory, GPU-ready

**kmcpy** is a Python package built on a fork of [KMC 3.2.4](https://github.com/refresh-bio/KMC/releases/tag/v3.2.4),
the fastest k-mer counter available.  
It eliminates all disk I/O for output, returning k-mer counts directly as **pandas DataFrames** — ready for downstream analysis, machine learning, and GPU-accelerated merges via **cuDF**.

[Installation](#installation) •
[Quick Start](#quick-start) •
[API Reference](#api-reference) •
[GPU Usage](#gpu-accelerated-merges-with-cudf) •
[Examples](tests/kmcpy/sample_usage.py) •
[Citing](#citing)

</div>

---

## What is this fork?

**kmcpy** is installable as a standard Python package:

```bash
pip install .
```

```python
import kmcpy

df = kmcpy.count_kmers("genomes.fasta", k=25, tmp_dir="./tmp")
# Returns a pandas DataFrame directly in RAM
#          kmer  count
# 0  ACGTACGT...      3
# 1  ACGTACGT...      1
```

[KMC](https://github.com/refresh-bio/KMC) is one of the fastest and most memory-efficient k-mer counters available. kmcpy is a Python-friendly extension of KMC 3.2.4 that modifies the Stage 2 output pipeline to return counted k-mers directly as a **pandas DataFrame**, with no disk writes for the final output. Stage 1 (splitting + binning) still spills to a configurable scratch directory when the RAM budget is exceeded, but the final k-mer table always stays in RAM — making kmcpy well suited for large-scale Python workflows where downstream analysis happens entirely in memory.



### Pipeline overview

```
Input files (FASTQ / FASTA / BAM)
        │
        ▼
┌─────────────────────────────────────────────────────┐
│  Stage 1: Splitting & Binning          (KMC, unmodified) │
│  Reads input, distributes k-mers into signature bins  │
│  Bins stay in RAM; spill to tmp_dir when max_ram_gb   │
│  is exceeded (hybrid RAM/disk, configurable)          │
└─────────────────────────────────────────────────────┘
        │  bin files (RAM or tmp_dir)
        ▼
┌─────────────────────────────────────────────────────┐
│  Stage 2: Sorting & Counting           (Fork: modified)  │
│  Each bin is radix-sorted and counted                │
│  withoutOutput=true → results go to packedKmers[]    │
│  and prefixArray[] in RAM instead of .kmc_pre/.suf   │
│                                                       │
│  Filters applied here (KMC, unmodified):             │
│    min_count  — discard k-mers below threshold        │
│    max_count  — discard k-mers above threshold        │
│    counter_max — saturate counter at this value       │
└─────────────────────────────────────────────────────┘
        │  packedKmers (flat byte array)
        │  prefixArray (uint64 per k-mer)
        ▼
┌─────────────────────────────────────────────────────┐
│  _core.cpp pybind11 binding            (Fork: new)       │
│  Decodes packed binary format into numpy arrays:     │
│    kmers_arr   (n, k)       uint8  — 2-bit base codes│
│    key_words   (n, n_words) uint64 — packed int keys │
│    counts_arr  (n,)         uint32 — occurrence counts│
│  drop_count=True → counts_arr never allocated        │
└─────────────────────────────────────────────────────┘
        │  numpy arrays (zero-copy to Python)
        ▼
┌─────────────────────────────────────────────────────┐
│  __init__.py Python API                (Fork: new)       │
│  decoded=True  → ACGT strings + count               │
│  decoded=False → packed uint64 keys + count         │
│  drop_count=True → keys only, no count column       │
└─────────────────────────────────────────────────────┘
        │
        ▼
  pandas DataFrame  →  cuDF (GPU)  →  your analysis
```

### Filters

All filters are applied by KMC during Stage 2 — before any data reaches Python.
Only k-mers that pass all filters are transferred across the C++/Python boundary.

| Parameter | Default | Effect |
|---|---|---|
| `min_count` | 1 | Discard k-mers with count < n. Default 1 keeps all k-mers including singletons. KMC CLI default is 2. |
| `max_count` | 1,000,000,000 | Discard k-mers with count > n. Useful for removing highly repetitive elements. |
| `counter_max` | 65,535 | Saturate the counter at this value. Max 65,535 (2-byte bin format). |

### Key changes from KMC 3.2.4

| File | Change |
|---|---|
| `kmc_core/kmc_runner.h` | `Stage2Results` carries `packedKmers` / `prefixArray` buffers instead of `kmerTable`. `withoutOutput=true` by default. |
| `kmc_core/kb_completer.h/cpp` | Stage 2 accumulates raw packed bytes into `packedKmers` and `prefixArray` instead of writing `.kmc_pre` / `.kmc_suf`. Both normal-k and small-k (k≤13) paths modified. |
| `kmc_core/kmc.h` | `GetKmerTable()` writes into RAM buffers. Both normal-k and small-k paths. |
| `kmc_core/kxmer_set.h` | Multithreaded k-mer merger modified to write records to `out_buffer` for `OutputType::KMC` regardless of `without_output`, ensuring all k-mers are captured in the in-memory path. |
| `kmc_core/kb_completer.cpp` | Stage 2 output accumulates into `packedKmers` / `prefixArray` in RAM. Added `lut_prefix_len=0` safety path. |
| `kmc_core/raduls_impl.h` | Explicit template instantiations replace `InstantiateTempl` trick (GCC 11+ dead-code elimination fix). |
| `Makefile` | SIMD isolation flags (`-mno-avx`, `-mno-sse4.1` etc.) prevent CPU feature bleed between raduls translation units. |
| `kmcpy/_core.cpp` | New pybind11 binding. Decodes packed binary format into numpy arrays. Supports `drop_count` to skip counts allocation at C++ level. |
| `kmcpy/__init__.py` | New Python API. Temp dir management, signal handlers, stale dir cleanup, decoded/packed output. |

All C++ modifications are marked `// Fork:` in the source.

---

## Features

- **Zero output disk I/O** — k-mer counts come back as a pandas DataFrame, no files written
- **Hybrid RAM/disk Stage 1** — bins spill to scratch only when the RAM budget is exceeded
- **Space-optimal packed keys** — `decoded=False` returns `ceil(k*2/64)` uint64 columns per k-mer, valid for all k (1–256)
- **GPU-ready** — packed integer keys enable O(n) merge across thousands of genomes via cuDF
- **`drop_count=True`** — skip count allocation entirely at the C++ level; saves 400 MB per 100M k-mers
- **Parallel safe** — each call uses an isolated temp subdirectory; fully concurrent
- **Automatic cleanup** — temp dirs removed on normal exit, exception, SIGINT, SIGTERM; stale dirs from SIGKILL cleaned on next call

---

## Installation

### Requirements

| Requirement | Minimum |
|---|---|
| C++ compiler | g++ ≥ 7 or clang ≥ 5 (C++14) |
| CMake | ≥ 3.17 |
| Python | ≥ 3.8 |
| pybind11 | ≥ 2.10 |
| zlib | ≥ 1.2 |

---

### pip install (recommended)

Works with any Python 3.8+ environment. No conda required.

```bash
# 1. Install build tools (one-time)
pip install scikit-build-core pybind11 cmake ninja

# 2. Clone the repository
git clone https://github.com/M-Serajian/KMC-DataFrame.git
cd KMC-DataFrame

# 3. Build and install
pip install .

# 4. Optional: install pytest to run the test suite
pip install pytest
```

### Uninstalling

```bash
pip uninstall kmcpy

# Also clean build artifacts from the source tree
make clean
rm -rf bin/ include/ _skbuild/
```

---

### Verify installation

> **Important:** run this from **outside** the project root directory.
> Running `import kmcpy` from inside `KMC-DataFrame/` causes Python to find
> the local `kmcpy/` source folder before the installed package and raises
> `ImportError: cannot import name '_core'`.

```bash
# Move out of the project root first
cd ~

# Then verify
python -c "import kmcpy; print(kmcpy.__version__)"
# 1.0.0
```

If you see `ImportError: cannot import name '_core'`, see [Troubleshooting](#troubleshooting).

---

### conda_setup.py (advanced — fully isolated build)

For HPC environments where you need a fully isolated build with a
conda-managed compiler and no dependency on the system g++.
Requires `conda` on `PATH`.

```bash
# Install — creates KMC-env conda environment with its own
# g++, zlib, pybind11, numpy, pandas and builds everything.
python conda_setup.py install

# Activate before use
conda activate KMC-env

# Verify (from outside the project root)
cd ~
python -c "import kmcpy; print(kmcpy.__version__)"
# 1.0.0

# Uninstall — deactivate first, then run from base environment
conda deactivate
python conda_setup.py uninstall
```

---

## Quick Start

> **Full working examples:** [`tests/kmcpy/sample_usage.py`](tests/kmcpy/sample_usage.py) — see also [`tests/kmcpy/README.md`](tests/kmcpy/README.md)
>
> Uses the committed [`tests/kmcpy/sample.fasta`](tests/kmcpy/sample.fasta)
> and writes temp files to `tests/kmcpy/tmp/` — **no setup required**.
>
> ```bash
> # From the project root
> python tests/kmcpy/sample_usage.py
>
> # From inside tests/kmcpy/
> cd tests/kmcpy && python sample_usage.py
> ```
>
> **Important:** do not run `python -c "import kmcpy"` from the project root —
> Python will find the local `kmcpy/` source folder before the installed package.
> Always run scripts from outside the project root or from inside `tests/kmcpy/`.

```python
import kmcpy

# Count 25-mers from a multi-line FASTA assembly
df = kmcpy.count_kmers(
    "assembly.fasta",
    k          = 25,
    input_fmt  = "mfasta",   # multi-line FASTA (most assembler outputs)
    tmp_dir    = "/scratch/$USER/kmc_tmp",
)
# Returns pandas DataFrame:
#          kmer  count
# 0  ACGTACGT...      3
# 1  ACGTACGT...      1
# ...

# Count from multiple files (treated as one combined dataset)
df = kmcpy.count_kmers(
    ["genome1.fasta", "genome2.fasta", "genome3.fasta"],
    k         = 25,
    input_fmt = "mfasta",
    tmp_dir   = "/scratch/$USER/kmc_tmp",
)

# Paired-end FASTQ, 8 threads, 32 GB RAM
df = kmcpy.count_kmers(
    ["r1.fastq.gz", "r2.fastq.gz"],
    k          = 31,
    tmp_dir    = "/scratch/$USER/kmc_tmp",
    threads    = 8,
    max_ram_gb = 32,
)

# Show full help
kmcpy.help()
```

---

## API Reference

```python
kmcpy.count_kmers(
    input_files,          # str or list[str] — FASTQ/FASTA/BAM paths
    tmp_dir,              # str — REQUIRED. Scratch directory for spill bins.
    k             = 25,   # k-mer length (1–256)
    input_fmt     = "fastq",  # fastq | fasta | mfasta | bam | kmc
    max_ram_gb    = 12,   # RAM budget in GB; bins spill to tmp_dir above this
    threads       = 0,    # 0 = all CPU cores
    canonical     = True, # count canonical k-mers
    min_count     = 1,    # exclude k-mers with count < n (default: keep all)
    max_count     = 1_000_000_000,
    counter_max   = 65535,
    strict_mem    = True, # never exceed max_ram_gb during Stage 2 sorting
    decoded       = True, # True = ACGT strings; False = packed uint64 keys
    drop_count    = False,# True = omit count column (saves n*4 bytes at C++ boundary)
) -> pandas.DataFrame
```

### Output formats

| `decoded` | `drop_count` | Columns | Use case |
|---|---|---|---|
| `True` | `False` | `kmer` (str), `count` (uint32) | human-readable, export |
| `True` | `True` | `kmer` (str) only | presence/absence |
| `False` | `False` | `kmer_0`…`kmer_{n-1}` (uint64), `count` (uint32) | cuDF merge with counts |
| `False` | `True` | `kmer_0`…`kmer_{n-1}` (uint64) only | cuDF merge, max memory efficiency |

`n = ceil(k * 2 / 64)` — the minimum number of uint64 columns to represent k bases.

| k range | n columns | bytes/k-mer key |
|---|---|---|
| k ≤ 32 | 1 | 8 |
| k ≤ 64 | 2 | 16 |
| k ≤ 128 | 4 | 32 |
| k ≤ 256 | 8 | 64 |

### `drop_count` memory savings

`drop_count=True` eliminates the counts array **at the C++ level** — the buffer is never allocated, filled, or transferred across the Python/C++ boundary.

| Unique k-mers | Saved |
|---|---|
| 10 M | 40 MB |
| 100 M | 400 MB |
| 1 B | 4 GB |
| 8 parallel jobs × 100 M | 3.2 GB |

---

## GPU-accelerated merges with cuDF

`decoded=False` returns space-optimal packed integer keys that enable O(n) hash joins on GPU via [cuDF](https://docs.rapids.ai/api/cudf/stable/).

```python
import kmcpy
import cudf

# Count k-mers for two samples
df1 = kmcpy.count_kmers("sample1.fasta", k=25, input_fmt="mfasta",
                          tmp_dir="/scratch", decoded=False)
df2 = kmcpy.count_kmers("sample2.fasta", k=25, input_fmt="mfasta",
                          tmp_dir="/scratch", decoded=False)

# Transfer to GPU
gdf1 = cudf.DataFrame.from_pandas(df1)
gdf2 = cudf.DataFrame.from_pandas(df2)

# Merge on composite integer key — identical syntax for ALL k (1–256)
n_words  = len([c for c in df1.columns if c.startswith("kmer_")])
key_cols = [f"kmer_{i}" for i in range(n_words)]
merged   = gdf1.merge(gdf2, on=key_cols, suffixes=("_s1", "_s2"))
```

The composite key `(kmer_0, kmer_1, ...)` preserves lexicographic order because
A=00 < C=01 < G=10 < T=11 matches A < C < G < T alphabetically.

---

## Parallel execution

Each `count_kmers()` call creates an isolated temp subdirectory and is safe to run concurrently.
Budget `max_workers × max_ram_gb` must fit in total available RAM.

```python
from concurrent.futures import ProcessPoolExecutor
import functools
import kmcpy

fn = functools.partial(
    kmcpy.count_kmers,
    k          = 25,
    input_fmt  = "mfasta",
    tmp_dir    = "/scratch/$USER/kmc_tmp",
    threads    = 4,
    max_ram_gb = 8,
    decoded    = False,
    drop_count = True,
)

genome_files = ["g1.fasta", "g2.fasta", ...]   # thousands of genomes

with ProcessPoolExecutor(max_workers=8) as pool:
    results = list(pool.map(fn, genome_files))
```

---

## Temp directory lifecycle

`tmp_dir` is **required** and must be:
- A directory you have **write permission** to
- On a filesystem with **sufficient free space** (at minimum: input file size × 5)
- Ideally a **fast local or scratch filesystem** — avoid RAM-backed `/tmp` for large datasets

The default in all examples is `tests/kmcpy/tmp/` which is pre-created in the
repository and requires no setup. For production use with large genomes:

```bash
# Use any writable directory you own
df -h                              # check available space on filesystems
mkdir -p ~/kmc_tmp                 # home directory (always works, may be slow)
mkdir -p /tmp/$USER/kmc_tmp        # local /tmp (fast, usually wiped on reboot)

# On HPC with a dedicated scratch filesystem
mkdir -p /scratch/$USER/kmc_tmp   # check your HPC docs for the correct path
```

```
count_kmers() creates:    tmp_dir/kmcpy_<PID>_<UUID>/
Removed on:               normal return, exception, SIGINT, SIGTERM
Not removed on:           SIGKILL — cleaned automatically on the next call
```

---

## Requirements

| Requirement | Minimum |
|---|---|
| C++ compiler | g++ ≥ 7 or clang ≥ 5 (C++14) |
| CMake | ≥ 3.17 |
| Python | ≥ 3.8 |
| pybind11 | ≥ 2.10 |
| zlib | ≥ 1.2 |
| numpy | ≥ 1.24 |
| pandas | ≥ 2.0 |
| cuDF | ≥ 23.0 (optional, GPU only) |

---

## Project structure

```
KMC-DataFrame/
├── kmc_core/           KMC engine (forked from KMC 3.2.4, modifications marked // Fork:)
├── kmc_api/            KMC C++ API (unmodified)
├── kmc_CLI/            KMC command-line interface (unmodified)
├── 3rd_party/          cloudflare zlib (unmodified)
├── kmcpy/
│   ├── __init__.py     Public Python API
│   ├── _core.cpp       pybind11 binding
│   └── _version.py     Version constants
├── testing/            Test scripts
├── CMakeLists.txt      pip build system (scikit-build-core)
├── pyproject.toml      Python packaging metadata
├── environment.yml     conda environment specification
└── conda_setup.py      Fully isolated conda build (legacy method)
```

---

## Citing

If you use kmcpy in your research, please cite the original KMC papers:

**KMC 3** (recommended):
> Marek Kokot, Maciej Długosz, Sebastian Deorowicz,
> *KMC 3: counting and manipulating k-mer statistics*,
> Bioinformatics, Volume 33, Issue 17, 2017, Pages 2759–2761.
> https://doi.org/10.1093/bioinformatics/btx304

**KMC 2**:
> Sebastian Deorowicz, Marek Kokot, Szymon Grabowski, Agnieszka Debudaj-Grabysz,
> *KMC 2: fast and resource-frugal k-mer counting*,
> Bioinformatics, Volume 31, Issue 10, 2015, Pages 1569–1576.
> https://doi.org/10.1093/bioinformatics/btv022

**KMC 1**:
> Deorowicz, S., Debudaj-Grabysz, A. & Grabowski, S.,
> *Disk-based k-mer counting on a PC*,
> BMC Bioinformatics 14, 160 (2013).
> https://doi.org/10.1186/1471-2105-14-160

---

## License

kmcpy is a fork of KMC and is distributed under the
**GNU General Public License v3.0** — the same license as the original KMC software.

| Component | License | Authors |
|---|---|---|
| KMC engine (`kmc_core/`, `kmc_api/`, `kmc_CLI/`) | GPL-3.0 | Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot |
| kmcpy modifications (`kmcpy/`, `Makefile`, `CMakeLists.txt`) | GPL-3.0 | kmcpy contributors |
| pybind11 ≥ 2.10 (`py_kmc_api/libs/pybind11/`) | BSD-3-Clause | Wenzel Jakob et al. |

See [LICENSE](LICENSE) for the full GPL v3 text.

> **Note:** Because KMC is GPL v3, any fork or derivative work — including kmcpy —
> must also be distributed under GPL v3. This is a legal requirement of the
> GPL v3 copyleft terms, not a choice.

### Disclaimer

This software is provided **as-is**, without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose, and non-infringement. In no event shall the
authors or copyright holders be liable for any claim, damages, or other
liability, whether in an action of contract, tort, or otherwise, arising from,
out of, or in connection with the software or the use or other dealings in the
software.

For the full warranty disclaimer see the
[GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.html),
sections 15 and 16.

---

## Acknowledgements

kmcpy is a fork of [KMC](https://github.com/refresh-bio/KMC) developed by the
[REFRESH Bioinformatics Group](http://sun.aei.polsl.pl/REFRESH/) at Silesian University of Technology.
The core counting engine, SIMD-accelerated sorting kernels, and bin management are
entirely the work of the original KMC authors.
This fork adds only the RAM-output pipeline and Python bindings.

---

---

## Troubleshooting

### `ImportError: undefined symbol: RadixSortMSD_SSE2`

The SIMD raduls kernels were not compiled with the correct isolation flags.
This happens when the build machine's CPU supports AVX2 and the compiler
defines `__AVX2__` even when compiling the SSE2 translation unit.

```bash
# Clean all build artifacts and rebuild from scratch
make clean
rm -rf bin/ include/
pip uninstall kmcpy -y
pip install .
```

If the problem persists, verify the Makefile contains `-mno-avx -mno-avx2`
on the `raduls_sse2` compile rule:

```bash
grep "mno-avx" Makefile
# Expected output:
# $(CC) $(CFLAGS) -msse2 -mno-avx -mno-avx2 -c $< -o $@
```

### `ImportError: cannot import name '_core' from 'kmcpy'`

Python is finding the local `kmcpy/` source folder before the installed package.
This happens when you run a script from inside the project root directory —
Python adds the current directory to `sys.path` and finds `kmcpy/__init__.py`
there (which has no compiled `_core.so`) before finding the installed package.

```bash
# WRONG -- runs from project root, finds local kmcpy/ source folder
cd KMC-DataFrame
python -c "import kmcpy"           # fails
python my_script.py                # fails

# CORRECT -- run from outside the project root
cd ~
python -c "import kmcpy"           # works
python tests/kmcpy/sample_usage.py # works (resolves paths internally)

# CORRECT -- run from inside tests/kmcpy/
cd KMC-DataFrame/tests/kmcpy
python sample_usage.py             # works
python -c "import kmcpy"           # works

# Or run the diagnostic to confirm what is happening
cd ~
python /path/to/KMC-DataFrame/diagnostics/check_install.py
```

### `PermissionError: Cannot write to tmp_dir`

You do not have write permission to the specified directory.
Use a directory you own:

```bash
# Always works — your home directory
python tests/kmcpy/sample_usage.py --tmp-dir ~/kmc_tmp

# Local temp (fast, wiped on reboot)
mkdir -p /tmp/$USER/kmc_tmp
python tests/kmcpy/sample_usage.py --tmp-dir /tmp/$USER/kmc_tmp

# On HPC — check your system docs for the correct scratch path
# Common patterns: /scratch/$USER, /localscratch/$USER, $SCRATCH
```

The default (`tests/kmcpy/tmp/`) is pre-created in the repo and always writable —
no `--tmp-dir` argument needed for the built-in examples.

### `ValueError: tmp_dir does not exist`

Create the directory before calling `count_kmers()`:

```python
import os, kmcpy
os.makedirs("/scratch/$USER/kmc_tmp", exist_ok=True)
df = kmcpy.count_kmers(..., tmp_dir="/scratch/$USER/kmc_tmp")
```

### Build fails with `relocation R_X86_64_32 ... recompile with -fPIC`

The object files were compiled without `-fPIC`. Clean and rebuild:

```bash
make clean
rm -rf bin/ include/
pip uninstall kmcpy -y
pip install .
```

### Running the diagnostic suite

If you encounter any issue, run the full diagnostic first and share the output
when reporting a bug:

```bash
python diagnostics/check_install.py --verbose --tmp-dir /scratch/$USER/kmc_tmp
```

---

## Contributing and reporting issues

Bug reports, feature requests, and pull requests are welcome.

- **Report a bug:** [open an issue](https://github.com/M-Serajian/KMC-DataFrame/issues) and include the output of `python diagnostics/check_install.py --verbose`
- **Feature request:** open an issue describing the use case
- **Pull request:** fork the repo, make your changes (mark all C++ modifications with `// Fork:`), and open a PR against `master`

When reporting a bug please include:
- Output of `python diagnostics/check_install.py --verbose`
- Your OS, Python version, and compiler version (`g++ --version`)
- Whether you used pip or conda_setup.py to install

---

<div align="center">

Built on [KMC 3.2.4](https://github.com/refresh-bio/KMC/releases/tag/v3.2.4) •
[Report an issue](https://github.com/M-Serajian/KMC-DataFrame/issues) •
[KMC original repository](https://github.com/refresh-bio/KMC)

</div>