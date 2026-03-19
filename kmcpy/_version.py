# kmcpy/_version.py
# Single source of truth for the kmcpy version.
# Referenced by:
#   - kmcpy/__init__.py  (runtime: kmcpy.__version__)
#   - pyproject.toml     (build: dynamic version)

__version__ = "1.0.0"

# KMC engine version this package was built against
__kmc_engine_version__ = "3.2.4"