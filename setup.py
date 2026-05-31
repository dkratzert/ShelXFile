"""
setup.py — responsible only for the optional sdm_cpp C++ extension.
All project metadata lives in pyproject.toml.

The sdm_cpp extension is built only when pybind11 is available.  If it is
absent the package still installs and works correctly using the pure-Python
SDM implementation.  To enable the C++ acceleration:

    pip install pybind11
    pip install -e . --no-build-isolation

OpenMP (optional, for multi-threaded acceleration)
----------------------------------------------------
macOS  : brew install libomp
Linux  : installed with GCC by default
Windows: included with MSVC
"""
from __future__ import annotations

import subprocess
import sys
from setuptools import setup

try:
    from setuptools import Extension
    import pybind11

    pybind11_include = pybind11.get_include()

    # ── OpenMP detection ──────────────────────────────────────────────────────

    def _find_openmp() -> tuple[list[str], list[str]]:
        """Return (extra_compile_args, extra_link_args) for OpenMP, or ([], [])."""
        if sys.platform == "darwin":
            import os
            archflags = os.environ.get("ARCHFLAGS", "")
            if "x86_64" in archflags:
                print("[sdm_cpp] Cross-compilation detected — skipping OpenMP.")
                return [], []
            try:
                result = subprocess.run(
                    ["brew", "--prefix", "libomp"],
                    capture_output=True, text=True, timeout=15,
                )
                if result.returncode == 0:
                    prefix = result.stdout.strip()
                    omp_h = f"{prefix}/include/omp.h"
                    omp_lib = f"{prefix}/lib/libomp.dylib"
                    import os.path
                    if os.path.isfile(omp_h) and os.path.isfile(omp_lib):
                        print(f"[sdm_cpp] OpenMP found via Homebrew libomp: {prefix}")
                        return (
                            ["-Xpreprocessor", "-fopenmp", f"-I{prefix}/include"],
                            [f"-L{prefix}/lib", "-lomp"],
                        )
                    print(f"[sdm_cpp] Homebrew libomp found but headers/lib missing — skipping.")
            except (FileNotFoundError, subprocess.TimeoutExpired):
                pass
            print("[sdm_cpp] libomp not found — building without OpenMP (single-threaded).")
            return [], []

        if sys.platform.startswith("linux"):
            return ["-fopenmp"], ["-fopenmp"]

        if sys.platform == "win32":
            return ["/openmp"], []

        return [], []

    omp_compile, omp_link = _find_openmp()

    if sys.platform == "win32":
        base_compile = ["/O2", "/std:c++17"] + omp_compile
    else:
        base_compile = ["-O3", "-std=c++17"] + omp_compile

    sdm_cpp_ext = Extension(
        name="sdm_cpp",
        sources=["src/sdm_cpp/sdm_cpp.cpp"],
        include_dirs=[pybind11_include],
        extra_compile_args=base_compile,
        extra_link_args=omp_link,
        language="c++",
    )

    setup(ext_modules=[sdm_cpp_ext])
    print("[sdm_cpp] C++ extension will be compiled.")

except ImportError:
    print(
        "[sdm_cpp] pybind11 not found — skipping C++ extension build.\n"
        "          Install pybind11 and re-run to enable the C++ acceleration:\n"
        "              pip install pybind11 && pip install -e . --no-build-isolation"
    )
    setup()

