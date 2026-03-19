#!/usr/bin/env python3
"""
conda_setup.py
--------------
Manages the KMC-env conda environment and builds the KMC project.
This is the legacy fully-isolated build method using conda.
For most users, the recommended installation is:  pip install .

Design principle: fully isolated
---------------------------------
The conda environment provides EVERYTHING needed to build KMC:

  - g++ compiler         (gxx_linux-64)
  - ar, ld, as           (binutils)
  - make                 (make)
  - git                  (git)
  - zlib headers + libz  (zlib)
  - libstdc++.a (static) (libstdcxx-static)
  - libgcc.a    (static) (libgcc-ng)

No system library or tool is used.  The build is fully reproducible
regardless of what the HPC system provides.

The Makefile is patched to:
  1. Use the conda zlib instead of the cloudflare submodule.
  2. Use the conda g++/ld/ar instead of whatever 'g++' is on PATH.
  3. Point the static linker at the conda static libraries.

Usage
-----
    python setup.py install     -- create environment, patch Makefile, build
    python setup.py uninstall   -- remove environment, restore Makefile, clean

Requirements
------------
    - conda (Anaconda or Miniconda) available on PATH
    - Script must be run from the root of the KMC repository
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ENV_NAME        = "KMC-env"
MAKEFILE        = Path("Makefile")
MAKEFILE_BACKUP = Path("Makefile.original")
CONDA_CHANNEL   = "conda-forge"

# Everything the build needs — fully isolated from the system.
#
#   gxx_linux-64      g++ C++14 compiler and its driver wrappers
#   binutils          ar, ld, as  (ar rcs libkmc_core.a; final link)
#   make              GNU make
#   git               git submodule update --init (ntHash)
#   zlib              zlib.h + libz.so/.a  (replaces cloudflare submodule)
#   libgcc-ng         libgcc shared + import stubs
#   libstdcxx-ng      libstdc++ shared
#   libstdcxx-static  libstdc++.a  — required by -static in STATIC_LFLAGS
#   glibc             libc.a, libpthread.a, libm.a — required by -static
#   glibc-static      static versions of glibc (libc.a, libpthread.a, etc.)


CONDA_PACKAGES = [
    "gxx_linux-64",   # g++ compiler + wrapper scripts
                      # auto-pulls: libstdcxx-devel_linux-64 (libstdc++.a)
                      #             libgcc-devel_linux-64 (libgcc.a)
                      #             sysroot_linux-64 (libc.a, libpthread.a)
    "binutils",       # ar, ld, as
    "make",           # GNU make
    "git",            # git submodule update --init
    "zlib",           # zlib.h + libz  (replaces cloudflare submodule)
    "libgcc-ng",      # libgcc shared runtime
    "libstdcxx-ng",   # libstdc++ shared runtime
    "python=3.12",    
    "pybind11=2.13.6", 
    "numpy>=1.24",     # ← add this
    "pandas>=2.0"
]


# ---------------------------------------------------------------------------
# Terminal colours
# ---------------------------------------------------------------------------
class Colour:
    if sys.stdout.isatty() and platform.system() != "Windows":
        RED    = "\033[0;31m"
        GREEN  = "\033[0;32m"
        YELLOW = "\033[1;33m"
        CYAN   = "\033[0;36m"
        RESET  = "\033[0m"
    else:
        RED = GREEN = YELLOW = CYAN = RESET = ""


def log_info(msg):  print("{c}[INFO]{r}  {m}".format(c=Colour.CYAN,   r=Colour.RESET, m=msg))
def log_ok(msg):    print("{c}[OK]{r}    {m}".format(c=Colour.GREEN,   r=Colour.RESET, m=msg))
def log_warn(msg):  print("{c}[WARN]{r}  {m}".format(c=Colour.YELLOW,  r=Colour.RESET, m=msg))
def log_error(msg): print("{c}[ERROR]{r} {m}".format(c=Colour.RED,     r=Colour.RESET, m=msg),
                          file=sys.stderr)


def abort(msg):
    log_error(msg)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Subprocess helper
# ---------------------------------------------------------------------------
def run(cmd, check=True, capture=False, env=None):
    return subprocess.run(
        cmd,
        check=check,
        capture_output=capture,
        text=True,
        env=env,
    )


# ---------------------------------------------------------------------------
# Conda helpers
# ---------------------------------------------------------------------------
def conda_exe():
    exe = shutil.which("conda")
    if not exe:
        abort("conda not found on PATH.")
    return exe


def conda_env_exists(name):
    result = run([conda_exe(), "env", "list"], capture=True)
    for line in result.stdout.splitlines():
        parts = line.split()
        if parts and parts[0] == name:
            return True
    return False


def conda_prefix(name):
    """
    Parse the environment path from 'conda env list'.
    The path is always the last whitespace-separated token on the line.
    """
    result = run([conda_exe(), "env", "list"], capture=True)
    for line in result.stdout.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if parts[0] == name:
            prefix = Path(parts[-1])
            if prefix.is_dir():
                return prefix
            abort("Environment '{}' listed but prefix missing: {}".format(name, prefix))
    abort("Environment '{}' not found in 'conda env list'.".format(name))


def check_repo_root():
    if not MAKEFILE.exists() or not Path("kmc_core").is_dir():
        abort(
            "setup.py must be run from the root of the KMC repository "
            "(Makefile and kmc_core/ not found)."
        )


def cpu_count():
    return os.cpu_count() or 1


# ---------------------------------------------------------------------------
# Locate tools inside the conda environment
# ---------------------------------------------------------------------------
def find_conda_tool(prefix, *candidates):
    """
    Search for a binary inside the conda environment's bin/ directory.
    Tries each candidate name in order.  Aborts if none is found.
    """
    bin_dir = prefix / "bin"
    for name in candidates:
        path = bin_dir / name
        if path.exists():
            return str(path)
    abort(
        "Could not find any of {} in {}.\n"
        "The conda environment may not have installed correctly.".format(
            candidates, bin_dir
        )
    )


# ---------------------------------------------------------------------------
# Makefile patching
# ---------------------------------------------------------------------------
def _check(text, old, label):
    if old not in text:
        abort(
            "{} — expected string not found in Makefile:\n"
            "  '{}'\n"
            "The Makefile may already be patched or is a different "
            "version.".format(label, old)
        )


def patch_makefile(prefix):
    """
    Patch the Makefile to use the conda environment exclusively.

    Patch 1  LIB_ZLIB cleared; ZLIB_INC and ZLIB_LIB point at conda.
    Patch 2  CLINK gets -lz -L$(ZLIB_LIB) for linking.
    Patch 3  Compile rule: -I 3rd_party/cloudflare -> -I$(ZLIB_INC).
    Patch 4  kmc_tools link rule: same include replacement.
    Patch 5  Cloudflare zlib build target commented out.
    Patch 6  CC overridden to use the conda g++ wrapper.
             The conda gxx_linux-64 package installs a wrapper script
             (x86_64-conda-linux-gnu-g++) that sets correct sysroot,
             include paths, and library paths inside the conda env.
             Using it instead of bare 'g++' gives the linker full access
             to libstdc++.a, libc.a, libpthread.a etc. from the conda env.
    """
    zlib_inc = prefix / "include"
    zlib_lib = prefix / "lib"

    # The conda gxx_linux-64 package installs two wrappers:
    #   x86_64-conda-linux-gnu-g++   (cross-compile style name)
    #   $prefix/bin/g++              (may or may not exist depending on version)
    # We prefer the explicit triple-prefixed one as it is always present.
    conda_gxx = find_conda_tool(
        prefix,
        "x86_64-conda-linux-gnu-g++",
        "g++",
    )
    conda_ar  = find_conda_tool(
        prefix,
        "x86_64-conda-linux-gnu-ar",
        "ar",
    )

    log_info("Conda prefix : {}".format(prefix))
    log_info("Conda g++    : {}".format(conda_gxx))
    log_info("Conda ar     : {}".format(conda_ar))
    log_info("zlib headers : {}".format(zlib_inc / "zlib.h"))
    log_info("zlib library : {}".format(zlib_lib))

    if not (zlib_inc / "zlib.h").exists():
        abort("zlib.h not found at {}.".format(zlib_inc))

    # Locate the GCC internal library directory inside the conda env.
    # gxx_linux-64 auto-installs libstdcxx-devel_linux-64 which places
    # libstdc++.a at:  prefix/lib/gcc/x86_64-conda-linux-gnu/VERSION/
    # We pass this directory to the linker via -L so the -static flag
    # can find libstdc++.a without relying on any system path.
    import glob as _glob
    gcc_lib_matches = _glob.glob(
        str(prefix / "lib/gcc/x86_64-conda-linux-gnu/*/libstdc++.a")
    )
    if not gcc_lib_matches:
        abort(
            "libstdc++.a not found under {}/lib/gcc/.\n"
            "Try reinstalling the environment: python setup.py uninstall && "
            "python setup.py install".format(prefix)
        )
    gcc_lib_dir = str(Path(gcc_lib_matches[0]).parent)
    log_ok("libstdc++.a  : {}".format(gcc_lib_matches[0]))
    log_ok("gcc lib dir  : {}".format(gcc_lib_dir))

    # Back up original Makefile
    if not MAKEFILE_BACKUP.exists():
        shutil.copy(MAKEFILE, MAKEFILE_BACKUP)
        log_info("Original Makefile backed up as {}".format(MAKEFILE_BACKUP))
    else:
        log_warn("Backup {} already exists — skipping.".format(MAKEFILE_BACKUP))

    text = MAKEFILE_BACKUP.read_text()

    # ------------------------------------------------------------------
    # Patch 1 — conda zlib paths
    # ------------------------------------------------------------------
    old = "LIB_ZLIB=3rd_party/cloudflare/libz.a"
    new = (
        "# Fork: cloudflare zlib replaced by conda zlib (see setup.py)\n"
        "LIB_ZLIB=\n"
        "ZLIB_INC={inc}\n"
        "ZLIB_LIB={lib}"
    ).format(inc=zlib_inc, lib=zlib_lib)
    _check(text, old, "Patch 1")
    text = text.replace(old, new, 1)

    # ------------------------------------------------------------------
    # Patch 2 — add -lz and -L$(ZLIB_LIB) to linker.
    # ------------------------------------------------------------------
    old = "CLINK\t= -lm $(STATIC_LFLAGS) -O3 -std=c++14"
    new = "CLINK\t= -lm -lz -L$(ZLIB_LIB) $(STATIC_LFLAGS) -O3 -std=c++14"
    _check(text, old, "Patch 2")
    text = text.replace(old, new, 1)

    # ------------------------------------------------------------------
    # Patch 3 — cloudflare include -> conda zlib include (compile rule)
    # ------------------------------------------------------------------
    old = "\t$(CC) $(CFLAGS) -I 3rd_party/cloudflare -c $< -o $@"
    new = "\t$(CC) $(CFLAGS) -I$(ZLIB_INC) -c $< -o $@"
    _check(text, old, "Patch 3")
    text = text.replace(old, new, 1)

    # ------------------------------------------------------------------
    # Patch 4 — cloudflare include -> conda zlib include (kmc_tools link)
    # ------------------------------------------------------------------
    old = "\t$(CC) $(CLINK) -I 3rd_party/cloudflare -o $(OUT_BIN_DIR)/$@ $^"
    new = "\t$(CC) $(CLINK) -I$(ZLIB_INC) -o $(OUT_BIN_DIR)/$@ $^"
    _check(text, old, "Patch 4")
    text = text.replace(old, new, 1)

    # ------------------------------------------------------------------
    # Patch 5 — disable cloudflare zlib build target
    # ------------------------------------------------------------------
    old = "\tcd 3rd_party/cloudflare; ./configure; make libz.a"
    new = "#\tcd 3rd_party/cloudflare; ./configure; make libz.a  # Fork: disabled"
    _check(text, old, "Patch 5")
    text = text.replace(old, new, 1)

    # ------------------------------------------------------------------
    # Patch 6 — override CC to the conda g++ wrapper
    #
    # The Makefile sets:  CC = g++  (on Linux)
    # We override it to the conda triple-prefixed wrapper which carries
    # the correct --sysroot, -L paths, and finds libstdc++.a inside the
    # conda environment instead of on the system.
    # ------------------------------------------------------------------
    old = "\tCC \t= g++"
    new = (
        "\t# Fork: use conda g++ for fully isolated build (see setup.py)\n"
        "\tCC \t= {gxx}"
    ).format(gxx=conda_gxx)
    _check(text, old, "Patch 6")
    text = text.replace(old, new, 1)

    # ------------------------------------------------------------------
    # Patch 7 — replace -static link flags with -pthread only.
    #
    # The Makefile uses -static on Linux x86-64 to produce a fully static
    # binary.  This requires libstdc++.a, libc.a, libpthread.a, etc. as
    # static archives AND requires that no dynamic library (.so) is linked.
    #
    # The conda toolchain's sysroot does not support a fully static link:
    # its libstdc++.a references __tls_get_addr from ld-linux.so (a DSO),
    # and conda zlib ships only libz.so (no libz.a).  Attempting -static
    # therefore always fails with:
    #   "attempted static link of dynamic object libz.so"
    #   "undefined reference to __tls_get_addr@@GLIBC_2.3"
    #
    # The fix: replace -static with just -pthread.  The conda g++ wrapper
    # already carries the correct --sysroot and rpath settings for full
    # isolation — no -static is needed.  The resulting binary links
    # dynamically against libz.so and libstdc++.so from the conda env
    # and is fully portable across nodes with the same conda env.
    # ------------------------------------------------------------------
    old_cflags = "\t\tSTATIC_CFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive"
    new_cflags  = "\t\tSTATIC_CFLAGS = -pthread  # Fork: -static removed, conda g++ handles isolation"
    _check(text, old_cflags, "Patch 7a (STATIC_CFLAGS)")
    text = text.replace(old_cflags, new_cflags, 1)

    old_lflags = "\t\tSTATIC_LFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive"
    new_lflags  = "\t\tSTATIC_LFLAGS = -pthread  # Fork: -static removed, conda g++ handles isolation"
    _check(text, old_lflags, "Patch 7b (STATIC_LFLAGS)")
    text = text.replace(old_lflags, new_lflags, 1)

    # ------------------------------------------------------------------
    # Patch 8 — use -fPIC (not -fPIE) in CFLAGS.
    #
    # -fPIE (Position Independent Executable) is only for final executables.
    # -fPIC (Position Independent Code) works for BOTH executables and .so.
    #
    # libkmc_core.a objects must be linked into:
    #   bin/kmc                       — an executable  (needs PIE or PIC)
    #   kmcpy/_core.cpython-*.so      — a shared library (REQUIRES fPIC)
    #
    # -fPIC is a strict superset of -fPIE: it satisfies the conda ld PIE
    # default AND allows the objects to be linked into a shared library.
    # Using -fPIE here causes the .so build to fail with:
    #   relocation R_X86_64_PC32 can not be used when making a shared object
    # ------------------------------------------------------------------
    old = "CFLAGS\t= -Wall -O3 -fsigned-char $(CPU_FLAGS) $(STATIC_CFLAGS) -std=c++14"
    new = "CFLAGS\t= -Wall -O3 -fsigned-char -fPIC $(CPU_FLAGS) $(STATIC_CFLAGS) -std=c++14  # Fork: -fPIC for both executable PIE and .so linking"
    _check(text, old, "Patch 8 (CFLAGS -fPIC)")
    text = text.replace(old, new, 1)

    MAKEFILE.write_text(text)
    log_ok("Makefile patched successfully (8 patches applied).")
    return conda_gxx


# ---------------------------------------------------------------------------
# install
# ---------------------------------------------------------------------------
def install():
    print()
    log_info("=" * 60)
    log_info("  Installing {}  (fully isolated)".format(ENV_NAME))
    log_info("=" * 60)

    check_repo_root()
    conda = conda_exe()
    log_info("conda    : {}".format(conda))
    log_info("channel  : {}".format(CONDA_CHANNEL))
    log_info("packages : {}".format(", ".join(CONDA_PACKAGES)))

    # ------------------------------------------------------------------
    # Step 1 — Create conda environment
    # ------------------------------------------------------------------
    if conda_env_exists(ENV_NAME):
        log_warn("Environment '{}' already exists — skipping creation.".format(ENV_NAME))
    else:
        log_info("Creating conda environment '{}' ...".format(ENV_NAME))
        run([
            conda, "create", "-y",
            "-n", ENV_NAME,
            "-c", CONDA_CHANNEL,
        ] + CONDA_PACKAGES)
        log_ok("Conda environment '{}' created.".format(ENV_NAME))

    # ------------------------------------------------------------------
    # Step 2 — Patch the Makefile
    # ------------------------------------------------------------------
    log_info("Resolving conda environment prefix ...")
    prefix = conda_prefix(ENV_NAME)
    log_ok("Environment prefix: {}".format(prefix))

    log_info("Patching Makefile ...")
    conda_gxx = patch_makefile(prefix)

    # ------------------------------------------------------------------
    # Step 3 — Initialise git submodules (ntHash)
    # ------------------------------------------------------------------
    log_info("Initialising git submodules ...")
    result = run(
        ["git", "submodule", "update", "--init", "--recursive"],
        check=False,
    )
    if result.returncode != 0:
        log_warn("git submodule update failed — continuing.")
    else:
        log_ok("Git submodules initialised.")

    # ------------------------------------------------------------------
    # Step 4 — Build using the conda make (fully isolated)
    # ------------------------------------------------------------------
    log_info("Cleaning stale artefacts ...")
    conda_make = find_conda_tool(prefix, "make")
    run([conda_make, "clean"], check=False)
    # The Makefile's clean target omits kmc_CLI/*.o — remove it explicitly
    # so stale objects compiled without -fPIE are never reused.
    import glob as _glob2
    for stale in _glob2.glob("kmc_CLI/*.o"):
        os.remove(stale)
        log_info("Removed stale: {}".format(stale))
    log_ok("Clean complete.")

    log_info("Building KMC using {} parallel jobs ...".format(cpu_count()))
    run([
        conda_make,
        "-j{}".format(cpu_count()),
        "kmc", "kmc_tools", "kmc_dump",
    ])

    # ------------------------------------------------------------------
    # Step 5 — Build the kmcpy Python extension (_core.so)
    #
    # Compiles kmcpy/_core.cpp against libkmc_core.a and the conda
    # pybind11/Python headers into a .so that Python can import.
    # The output filename includes the Python ABI tag automatically
    # via `python3-config --extension-suffix`, e.g.:
    #   kmcpy/_core.cpython-312-x86_64-linux-gnu.so
    # ------------------------------------------------------------------
    log_info("Building kmcpy Python extension ...")

    # Resolve python3 and python3-config from the conda env
    conda_python      = find_conda_tool(prefix, "python3", "python")
    conda_python_cfg  = find_conda_tool(prefix, "python3-config", "python-config")

    import subprocess as _sp

    # Get Python include flags:  -I/path/to/python3.12/include
    py_includes = _sp.check_output(
        [conda_python_cfg, "--includes"], text=True
    ).strip()

    # Get the .so extension suffix:  .cpython-312-x86_64-linux-gnu.so
    ext_suffix = _sp.check_output(
        [conda_python_cfg, "--extension-suffix"], text=True
    ).strip()

    core_src = Path("kmcpy") / "_core.cpp"
    core_out = Path("kmcpy") / ("_core" + ext_suffix)

    if not core_src.exists():
        abort(
            "kmcpy/_core.cpp not found. "
            "Make sure the kmcpy/ directory exists at the project root."
        )

    cmd = (
        "{gxx} -O3 -shared -fPIC -std=c++14 "
        "-fvisibility=hidden "
        "-I include "
        "-I {conda_inc} "
        "{py_inc} "
        "-o {out} "
        "{src} "
        "bin/libkmc_core.a "
        "-lz -L{conda_lib} "
        "-pthread "
        "-Wl,-rpath,{conda_lib}"
    ).format(
        gxx       = conda_gxx,
        conda_inc = str(prefix / "include"),
        py_inc    = py_includes,
        out       = str(core_out),
        src       = str(core_src),
        conda_lib = str(prefix / "lib"),
    )

    log_info("Compile: {}".format(cmd))
    result = run(cmd.split(), check=False)
    if result.returncode != 0:
        abort("kmcpy/_core build failed (see above).")

    log_ok("kmcpy extension: {}".format(core_out.resolve()))

    print()
    log_ok("=" * 60)
    log_ok("  Build complete.")
    log_ok("  Binaries : {}".format(Path("bin/kmc").resolve()))
    log_ok("  Library  : {}".format(Path("bin/libkmc_core.a").resolve()))
    log_ok("  Headers  : {}".format(Path("include/kmc_runner.h").resolve()))
    log_ok("  Python   : {}".format(core_out.resolve()))
    log_ok("=" * 60)
    log_info("Usage:")
    log_info("  conda activate KMC-env")
    log_info("  python -c \"import kmcpy; print(kmcpy.__version__)\"")
    print()


# ---------------------------------------------------------------------------
# uninstall
# ---------------------------------------------------------------------------
def uninstall():
    print()
    log_info("=" * 60)
    log_info("  Uninstalling {}".format(ENV_NAME))
    log_info("=" * 60)

    conda = conda_exe()

    # Step 1 — clean build artefacts
    if MAKEFILE.exists() and Path("kmc_core").is_dir():
        log_info("Cleaning build artefacts ...")
        run(["make", "clean"], check=False)
        log_ok("Build artefacts removed.")
    else:
        log_warn("Not in KMC repo root — skipping make clean.")

    # Remove compiled kmcpy extension
    import glob as _g
    for so in _g.glob("kmcpy/_core*.so"):
        os.remove(so)
        log_ok("Removed: {}".format(so))

    # Step 2 — restore original Makefile
    if MAKEFILE_BACKUP.exists():
        shutil.move(str(MAKEFILE_BACKUP), str(MAKEFILE))
        log_ok("Makefile restored from {}.".format(MAKEFILE_BACKUP))
    else:
        log_warn("No Makefile backup found — Makefile left as-is.")

    # Step 3 — remove conda environment
    if conda_env_exists(ENV_NAME):
        log_info("Removing conda environment '{}' ...".format(ENV_NAME))
        run([conda, "env", "remove", "-y", "-n", ENV_NAME])
        log_ok("Conda environment '{}' removed.".format(ENV_NAME))
    else:
        log_warn("Conda environment '{}' not found.".format(ENV_NAME))

    print()
    log_ok("=" * 60)
    log_ok("  Uninstall complete.")
    log_ok("=" * 60)
    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        prog="conda_setup.py",
        description=(
            "Manage the KMC-env conda environment and build KMC.\n"
            "The environment is fully isolated — no system libraries are used."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python setup.py install     # create environment and build\n"
            "  python setup.py uninstall   # remove environment and clean\n"
        ),
    )
    parser.add_argument(
        "command",
        choices=["install", "uninstall"],
        help="'install' to build, 'uninstall' to remove everything.",
    )
    args = parser.parse_args()

    if args.command == "install":
        install()
    elif args.command == "uninstall":
        uninstall()


if __name__ == "__main__":
    main()