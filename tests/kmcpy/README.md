# tests/kmcpy/

Test suite and usage examples for kmcpy.

## Structure

```
tests/kmcpy/
├── sample_usage.py     End-to-end usage examples (start here)
├── test_kmcpy.py       Regression test suite (requires pytest)
├── README.md           This file
└── tmp/                Scratch directory for KMC temp files (pre-created)
```

## Running the examples

No setup required — uses `tests/kmcpy/tmp/` automatically:

```bash
# From the project root
python tests/kmcpy/sample_usage.py

# From inside tests/kmcpy/
python sample_usage.py

# For large datasets, use a faster scratch filesystem
python tests/kmcpy/sample_usage.py --tmp-dir /path/to/scratch
```

## Running the tests

First install pytest if not already installed:

```bash
pip install pytest
```

Then run:

```bash
# From the project root
pytest tests/kmcpy/test_kmcpy.py -v

# From inside tests/kmcpy/
pytest test_kmcpy.py -v
```

## Note

Always run from **outside** the project source root, or from inside
`tests/kmcpy/`. Running `python -c "import kmcpy"` from the project root
causes Python to find the local `kmcpy/` source folder before the installed
package — run `cd ..` first.