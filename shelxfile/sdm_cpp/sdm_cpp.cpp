// ----------------------------------------------------------------------------
// "THE BEER-WARE LICENSE" (Revision 42):
// <dkratzert@gmx.de> wrote this file. As long as you retain
// this notice you can do whatever you want with this stuff. If we meet some day,
// and you think this stuff is worth it, you can buy me a beer in return.
// Daniel Kratzert
// ----------------------------------------------------------------------------
//
// C++ acceleration for the Shortest-Distance-Matrix algorithm.
// This module is optional; sdm.py falls back to pure Python when it cannot be
// imported (HAS_CPP = False).
//
// Build:
//   macOS  : uv pip install pybind11 && uv pip install -e . --no-build-isolation
//            (optionally: brew install libomp  for multi-threaded acceleration)
//   Linux  : uv pip install pybind11 && uv pip install -e . --no-build-isolation
//   Windows: uv pip install pybind11 && uv pip install -e . --no-build-isolation
//
// Interface (mirrors the Python call in sdm.py):
//
//   calc_sdm_cpp(
//       coords   : list[list[float]],           # N × 3,  fractional coordinates
//       symm_m   : list[list[list[float]]],      # S × 3×3 symmetry rotation matrices
//       symm_t   : list[list[float]],            # S × 3  symmetry translation vectors
//       aga, bbe, cal, asq, bsq, csq : float,   # metric-tensor cross-terms / squares
//       radii    : list[float],                  # N covalent radii (Å), from Python
//       is_h     : list[bool],                   # N hydrogen flags
//       parts    : list[float],                  # N disorder-part numbers
//   ) -> list[tuple[int, int, int, float, float, bool]]
//          i   j  best_n  mind   dddd  covalent

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <cmath>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

// Result type: (atom_i, atom_j, symmetry_index, min_distance, covalent_cutoff, is_covalent)
using ResultTuple = std::tuple<int, int, int, double, double, bool>;

// ---------------------------------------------------------------------------
// calc_sdm_cpp
//
// Direct C++17 translation of the pure-Python fallback in sdm.py::SDM.calc_sdm().
// The outer i-loop is parallelised with OpenMP; each thread accumulates results
// in a thread-local vector that is merged serially after the parallel region.
// When OpenMP is not available at compile time the loop runs single-threaded.
// ---------------------------------------------------------------------------
std::vector<ResultTuple> calc_sdm_cpp(
    const std::vector<std::vector<double>>& coords,       // N × 3
    const std::vector<std::vector<std::vector<double>>>& symm_m,  // S × 3×3
    const std::vector<std::vector<double>>& symm_t,       // S × 3
    double aga, double bbe, double cal,
    double asq, double bsq, double csq,
    const std::vector<double>& radii,
    const std::vector<bool>& is_h,
    const std::vector<double>& parts
) {
    const int N = static_cast<int>(coords.size());
    const int S = static_cast<int>(symm_m.size());

    // Pre-compute at2 + 0.5 (at2_plushalf in Python)
    std::vector<std::array<double, 3>> at2_ph(N);
    for (int j = 0; j < N; ++j) {
        at2_ph[j] = {coords[j][0] + 0.5, coords[j][1] + 0.5, coords[j][2] + 0.5};
    }

    // Flatten symmetry data for cache-friendly access
    // symm_m_flat[s][row][col] and symm_t_flat[s][0..2]
    std::vector<std::array<std::array<double, 3>, 3>> sm(S);
    std::vector<std::array<double, 3>>                st(S);
    for (int s = 0; s < S; ++s) {
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                sm[s][r][c] = symm_m[s][r][c];
        for (int c = 0; c < 3; ++c)
            st[s][c] = symm_t[s][c];
    }

#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
    std::vector<std::vector<ResultTuple>> thread_results(nthreads);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        std::vector<ResultTuple>& local = thread_results[omp_get_thread_num()];
#else
    std::vector<ResultTuple> serial_results;
    serial_results.reserve(static_cast<size_t>(N) * N / 4);

    for (int i = 0; i < N; ++i) {
        std::vector<ResultTuple>& local = serial_results;
#endif
        const double x1 = coords[i][0];
        const double y1 = coords[i][1];
        const double z1 = coords[i][2];

        // Apply all symmetry operations to atom i → prime_array
        std::vector<std::array<double, 3>> prime(S);
        for (int s = 0; s < S; ++s) {
            prime[s][0] = x1 * sm[s][0][0] + y1 * sm[s][1][0] + z1 * sm[s][2][0] + st[s][0];
            prime[s][1] = x1 * sm[s][0][1] + y1 * sm[s][1][1] + z1 * sm[s][2][1] + st[s][1];
            prime[s][2] = x1 * sm[s][0][2] + y1 * sm[s][1][2] + z1 * sm[s][2][2] + st[s][2];
        }

        for (int j = 0; j < N; ++j) {
            double mind  = 1.0e6;
            int    best_n = -1;

            const double atp_x = at2_ph[j][0];
            const double atp_y = at2_ph[j][1];
            const double atp_z = at2_ph[j][2];

            for (int n = 0; n < S; ++n) {
                const double dx = prime[n][0] - atp_x;
                const double dy = prime[n][1] - atp_y;
                const double dz = prime[n][2] - atp_z;

                // Reduce to [-0.5, 0.5) range  (equivalent to Python's  dx - floor(dx) - 0.5)
                const double dpx = dx - std::floor(dx) - 0.5;
                const double dpy = dy - std::floor(dy) - 0.5;
                const double dpz = dz - std::floor(dz) - 0.5;

                const double A   = 2.0 * (dpx * dpy * aga + dpx * dpz * bbe + dpy * dpz * cal);
                const double dk2 = dpx * dpx * asq + dpy * dpy * bsq + dpz * dpz * csq + A;

                if (dk2 > 16.0) continue;  // distance > 4 Å — skip

                double dk = std::sqrt(dk2);
                if (n > 0) dk += 0.0001;   // symmetry-generated copies get a tiny penalty

                if (dk > 0.01 && mind >= dk) {
                    mind   = dk;
                    best_n = n;
                }
            }

            if (best_n < 0) continue;  // no valid distance found → skip pair

            // Compute covalent cutoff (dddd) — mirrors Python logic exactly:
            //   if (neither is H  AND  part product is 0) OR  same part  → use radii sum * 1.2
            //   else → 0.0 (non-bonded)
            const bool   ih_i   = is_h[i];
            const bool   ih_j   = is_h[j];
            const double part_i = parts[i];
            const double part_j = parts[j];

            double dddd;
            if ((!ih_i && !ih_j && part_i * part_j == 0.0) || part_i == part_j) {
                dddd = (radii[i] + radii[j]) * 1.2;
            } else {
                dddd = 0.0;
            }

            const bool covalent = mind < dddd;
            local.emplace_back(i, j, best_n, mind, dddd, covalent);
        }
    }  // end parallel for

#ifdef _OPENMP
    // Merge thread-local vectors into a single result
    std::vector<ResultTuple> results;
    for (auto& tr : thread_results)
        results.insert(results.end(), tr.begin(), tr.end());
    return results;
#else
    return serial_results;
#endif
}

// ---------------------------------------------------------------------------
// Module registration
// ---------------------------------------------------------------------------
PYBIND11_MODULE(sdm_cpp, m) {
    m.doc() = R"doc(
Fast C++ implementation of the Shortest-Distance-Matrix (SDM) algorithm
used by fastmolwidget to grow crystal structures.

Compiled with OpenMP support when libomp is available (macOS: optionally brew install libomp).
Falls back gracefully to single-threaded execution when OpenMP is absent.
)doc";

    m.def(
        "calc_sdm_cpp",
        &calc_sdm_cpp,
        py::arg("coords"),
        py::arg("symm_m"),
        py::arg("symm_t"),
        py::arg("aga"), py::arg("bbe"), py::arg("cal"),
        py::arg("asq"), py::arg("bsq"), py::arg("csq"),
        py::arg("radii"),
        py::arg("is_h"),
        py::arg("parts"),
        R"doc(
Calculate the Shortest Distance Matrix for all atom pairs.

Parameters
----------
coords  : list of [x, y, z] (fractional coordinates), length N
symm_m  : list of 3×3 symmetry rotation matrices, length S
symm_t  : list of 3-element translation vectors, length S
aga     : cell[0]*cell[1]*cos(gamma)
bbe     : cell[0]*cell[2]*cos(beta)
cal     : cell[1]*cell[2]*cos(alpha)
asq     : cell[0]**2
bsq     : cell[1]**2
csq     : cell[2]**2
radii   : covalent radii for each atom (Å), length N
is_h    : True if the atom is H or D, length N
parts   : SHELX disorder-part number for each atom, length N

Returns
-------
list of (i, j, best_n, mind, dddd, covalent)
  i, j     : atom indices
  best_n   : symmetry-operation index that gave the shortest distance
  mind     : shortest distance found (Å)
  dddd     : covalent-bond cutoff (0.0 if the pair cannot be bonded)
  covalent : True when mind < dddd
)doc"
    );

#ifdef _OPENMP
    m.attr("has_openmp") = true;
#else
    m.attr("has_openmp") = false;
#endif
}

