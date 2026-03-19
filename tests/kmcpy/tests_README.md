# tests/kmcpy/

Test suite and usage examples for kmcpy.

## Structure

```
tests/kmcpy/
├── sample.fasta        Sample input — 3 synthetic sequences, committed to repo
├── sample_usage.py     End-to-end usage examples (start here)
├── test_kmcpy.py       Regression test suite (requires pytest)
├── README.md           This file
└── tmp/                Scratch directory for KMC temp files (pre-created)
```

## Running the examples

No setup required. Uses `sample.fasta` (committed) and `tmp/` (pre-created):

```bash
# From the project root
python tests/kmcpy/sample_usage.py

# From inside tests/kmcpy/
cd tests/kmcpy
python sample_usage.py
```

For large real datasets, point to a faster scratch filesystem:

```bash
python tests/kmcpy/sample_usage.py --tmp-dir /path/to/scratch
```

## Running the tests

Install pytest if not already installed:

```bash
pip install pytest
```

Then run:

```bash
# From the project root
pytest tests/kmcpy/test_kmcpy.py -v

# From inside tests/kmcpy/
cd tests/kmcpy
pytest test_kmcpy.py -v
```

## Important: where to run from

**Do not run scripts from the project root with `python -c "import kmcpy"`.**
Python adds the current directory to `sys.path` and finds the local `kmcpy/`
source folder before the installed package, causing `ImportError`.

```bash
# WRONG
cd KMC-DataFrame
python -c "import kmcpy"        # ImportError

# CORRECT
cd KMC-DataFrame/tests/kmcpy
python -c "import kmcpy"        # works

# CORRECT
cd ~
python -c "import kmcpy"        # works
```
