# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
import time
from math import sqrt, radians, sin, floor
from pathlib import Path

import numpy as np

try:
    import sdm_cpp

    HAS_CPP = True
except ImportError:
    HAS_CPP = False

from shelxfile.atoms.atom import Atom
from shelxfile.misc.dsrmath import vol_unitcell
from shelxfile.misc.misc import wrap_line
from shelxfile.shelx.cards import AFIX, RESI


def _unique_legal_atom_name(element: str, used_names: set) -> str:
    """Return a unique, legal SHELXL atom name (≤4 alphanumeric chars).

    The name follows the pattern ``El1``, ``El2``, … up to the digit limit
    imposed by the 4-character maximum:
    - 1-char element symbol  → up to ``E999``
    - 2-char element symbol  → up to ``EE99``

    The generated name is added to *used_names* (upper-case) before returning
    so subsequent calls automatically avoid it.

    :raises ValueError: when the numbering range for the element is exhausted.
    """
    prefix = element.capitalize()
    max_digits = 4 - len(prefix)
    if max_digits < 1:
        raise ValueError(
            f"Element symbol '{element}' is too long to form a legal SHELXL atom name."
        )
    max_num = 10 ** max_digits - 1
    for n in range(1, max_num + 1):
        candidate = f'{prefix}{n}'
        if candidate.upper() not in used_names:
            used_names.add(candidate.upper())
            return candidate
    raise ValueError(
        f"Cannot generate a unique legal SHELXL atom name for element '{element}': "
        f"all {max_num} slots are taken."
    )


class SDMItem(object):
    __slots__ = ['dist', 'atom1', 'atom2', 'a1', 'a2', 'symmetry_number', 'covalent', 'dddd']

    def __init__(self):
        self.dist = 0.0
        self.atom1 = None
        self.a1 = 0
        self.atom2 = None
        self.a2 = 0
        self.symmetry_number = 0
        self.covalent = True
        self.dddd = 0

    def __lt__(self, a2):
        return True if self.dist < a2.dist else False

    def __eq__(self, other: 'SDMItem'):
        if other.a1 == self.a2 and other.a2 == self.a1:
            return True
        return False

    def __repr__(self):
        return '{} {} dist: {:.6} coval: {} sn: {}  {}\n'.format(
            self.a1, self.a2, self.dist, self.covalent, self.symmetry_number, self.dddd)


class SDM():
    """
    Calculates the shortest distance matrix and creates a completed (grown) structure
    by crystal symmetry.

    Algorithm improvements (matching Fastmolwidget):
    - Pre-computed symmetry matrices (column-major tuples) avoid per-iteration numpy calls.
    - Squared distance cutoff (dk2 > 16.0 Å²) eliminates the sqrt for far pairs.
    - Union-Find (path-halving + union-by-rank) replaces the O(K·|sdm_list|) molindex loop.
    - numpy is used for U-value transformations.
    """

    def __init__(self, shx: 'Shelxfile'):
        self.shx = shx
        self.cell = (
            shx.cell.a, shx.cell.b, shx.cell.c,
            shx.cell.alpha, shx.cell.beta, shx.cell.gamma)
        self.aga = shx.cell.a * shx.cell.b * shx.cell.cosga
        self.bbe = shx.cell.a * shx.cell.c * shx.cell.cosbe
        self.cal = shx.cell.b * shx.cell.c * shx.cell.cosal
        self.asq = shx.cell[0] ** 2
        self.bsq = shx.cell[1] ** 2
        self.csq = shx.cell[2] ** 2
        # reciprocal lattice vectors (used for U-value transforms and vector_length)
        self.astar = (shx.cell.b * shx.cell.c * sin(radians(shx.cell.alpha))) / shx.cell.V
        self.bstar = (shx.cell.c * shx.cell.a * sin(radians(shx.cell.beta))) / shx.cell.V
        self.cstar = (shx.cell.a * shx.cell.b * sin(radians(shx.cell.gamma))) / shx.cell.V
        self.sdm_list: list[SDMItem] = []
        self.maxmol = 1
        self.sdmtime = 0
        self.bondlist = []
        # Pre-computed symmetry matrices (set in _build_symm_arrays, reused everywhere)
        self._symm_m: list[tuple] = []
        self._symm_t: list[tuple] = []

    def _build_symm_arrays(self) -> None:
        """Pre-compute column-major symmetry matrices as nested tuples for fast inner loops.

        Each symmcard.matrix is a (3,3) numpy array stored row-major (matrix[i,j] is the
        coefficient of input axis j for output axis i).  Taking .T converts to column-major
        so the inner-loop formula
            px = x*m[0][0] + y*m[1][0] + z*m[2][0]
        correctly evaluates  matrix[0,:] · [x,y,z].
        """
        self._symm_m = []
        self._symm_t = []
        for s in self.shx.symmcards:
            self._symm_m.append(tuple(map(tuple, s.matrix.T)))
            self._symm_t.append(tuple(s.trans))

    def calc_sdm(self) -> list:
        t1 = time.perf_counter()
        self._build_symm_arrays()
        symm_m = self._symm_m
        symm_t = self._symm_t
        nlen = len(symm_m)

        all_atoms = self.shx.atoms.all_atoms
        self.sdm_list.clear()
        self.bondlist.clear()

        # Pre-extract primitive data to avoid per-iteration attribute lookups
        coords = [(at.x, at.y, at.z) for at in all_atoms]
        radii = [at.radius for at in all_atoms]
        is_h = [at.ishydrogen for at in all_atoms]
        parts = [at.part.n for at in all_atoms]

        aga, bbe, cal = self.aga, self.bbe, self.cal
        asq, bsq, csq = self.asq, self.bsq, self.csq

        # ── C++ fast path ──────────────────────────────────────────────────
        if HAS_CPP:
            cpp_results = sdm_cpp.calc_sdm_cpp(
                coords, symm_m, symm_t,
                aga, bbe, cal,
                asq, bsq, csq,
                radii, is_h, parts,
            )
            for (i, j, best_n, mind, dddd, covalent) in cpp_results:
                sdm_item = SDMItem()
                sdm_item.dist = mind
                sdm_item.atom1 = all_atoms[i]
                sdm_item.atom2 = all_atoms[j]
                sdm_item.a1 = i
                sdm_item.a2 = j
                sdm_item.symmetry_number = best_n
                sdm_item.dddd = dddd
                sdm_item.covalent = covalent
                self.sdm_list.append(sdm_item)

        # ── Pure-Python fallback ────────────────────────────────────────────
        else:
            at2_plushalf = [(x + 0.5, y + 0.5, z + 0.5) for (x, y, z) in coords]

            for i, at1 in enumerate(all_atoms):
                x1, y1, z1 = coords[i]

                # Apply all symmetry operations to at1 once (outer loop)
                prime_array = []
                for m, t in zip(symm_m, symm_t):
                    px = x1 * m[0][0] + y1 * m[1][0] + z1 * m[2][0] + t[0]
                    py = x1 * m[0][1] + y1 * m[1][1] + z1 * m[2][1] + t[1]
                    pz = x1 * m[0][2] + y1 * m[1][2] + z1 * m[2][2] + t[2]
                    prime_array.append((px, py, pz))

                for j, at2 in enumerate(all_atoms):
                    mind = 1_000_000.0
                    hma = False
                    atp_x, atp_y, atp_z = at2_plushalf[j]
                    sdm_item = SDMItem()

                    for n in range(nlen):
                        px, py, pz = prime_array[n]

                        dx = px - atp_x
                        dy = py - atp_y
                        dz = pz - atp_z

                        dpx = dx - floor(dx) - 0.5
                        dpy = dy - floor(dy) - 0.5
                        dpz = dz - floor(dz) - 0.5

                        A = 2.0 * (dpx * dpy * aga + dpx * dpz * bbe + dpy * dpz * cal)
                        dk2 = dpx * dpx * asq + dpy * dpy * bsq + dpz * dpz * csq + A

                        if dk2 > 16.0:  # 4 Å hard cutoff (squared) – skip sqrt
                            continue

                        dk = sqrt(dk2)
                        if n:
                            dk += 0.0001  # slight penalty for symmetry-generated images

                        if (dk > 0.01) and (mind >= dk):
                            mind = dk
                            sdm_item.dist = mind
                            sdm_item.atom1 = at1
                            sdm_item.atom2 = at2
                            sdm_item.a1 = i
                            sdm_item.a2 = j
                            sdm_item.symmetry_number = n
                            hma = True

                    if not sdm_item.atom1:
                        continue

                    if ((not is_h[i] and not is_h[j]) and parts[i] * parts[j] == 0) \
                            or parts[i] == parts[j]:
                        dddd = (radii[i] + radii[j]) * 1.2
                        sdm_item.dddd = dddd
                    else:
                        dddd = 0.0

                    if sdm_item.dist < dddd:
                        if hma:
                            sdm_item.covalent = True
                    else:
                        sdm_item.covalent = False

                    if hma:
                        self.sdm_list.append(sdm_item)

        t2 = time.perf_counter()
        self.sdmtime = t2 - t1
        if self.shx.debug:
            print(f'SDM {"(C++)" if HAS_CPP else "(Python fallback)"}: {self.sdmtime:.4f} s')

        self.sdm_list.sort()
        self.calc_molindex(all_atoms)
        need_symm = self.collect_needed_symmetry()
        if self.shx.debug:
            print(f"The asymmetric unit contains {self.maxmol} fragments.")
        return need_symm

    def collect_needed_symmetry(self) -> list:
        need_symm = []
        symm_m = self._symm_m
        symm_t = self._symm_t

        aga, bbe, cal = self.aga, self.bbe, self.cal
        asq, bsq, csq = self.asq, self.bsq, self.csq

        for sdm_item in self.sdm_list:
            if not sdm_item.covalent:
                continue
            if sdm_item.atom1.molindex < 1 or sdm_item.atom1.molindex > 6:
                continue

            x1, y1, z1 = sdm_item.atom1.x, sdm_item.atom1.y, sdm_item.atom1.z
            x2, y2, z2 = sdm_item.atom2.x, sdm_item.atom2.y, sdm_item.atom2.z
            part1 = sdm_item.atom1.part.n
            part2 = sdm_item.atom2.part.n
            is_h1 = sdm_item.atom1.ishydrogen
            is_h2 = sdm_item.atom2.ishydrogen

            for n, (m, t) in enumerate(zip(symm_m, symm_t)):
                if part1 != 0 and part2 != 0 and part1 != part2:
                    continue
                # Skip H–H pairs of the same atomic number
                if is_h1 and is_h2 and sdm_item.atom1.an == sdm_item.atom2.an:
                    continue

                px = x1 * m[0][0] + y1 * m[1][0] + z1 * m[2][0] + t[0]
                py = x1 * m[0][1] + y1 * m[1][1] + z1 * m[2][1] + t[1]
                pz = x1 * m[0][2] + y1 * m[1][2] + z1 * m[2][2] + t[2]

                Dx = px - x2 + 0.5
                Dy = py - y2 + 0.5
                Dz = pz - z2 + 0.5

                fDx = floor(Dx)
                fDy = floor(Dy)
                fDz = floor(Dz)

                if n == 0 and fDx == 0 and fDy == 0 and fDz == 0:
                    continue

                dpx = Dx - fDx - 0.5
                dpy = Dy - fDy - 0.5
                dpz = Dz - fDz - 0.5

                A = 2.0 * (dpx * dpy * aga + dpx * dpz * bbe + dpy * dpz * cal)
                dk2 = dpx * dpx * asq + dpy * dpy * bsq + dpz * dpz * csq + A

                if dk2 <= 1e-6:
                    continue

                dk = sqrt(dk2)
                dddd = sdm_item.dist + 0.2
                if is_h1 and is_h2:
                    dddd = 1.8

                if dk <= dddd:
                    bs = [n + 1, int(5 - fDx), int(5 - fDy), int(5 - fDz), sdm_item.atom1.molindex]
                    if bs not in need_symm:
                        need_symm.append(bs)

        return need_symm

    def calc_molindex(self, all_atoms: list) -> None:
        """Assign a molecule index to every atom using Union-Find (path-halving +
        union-by-rank).  Replaces the original O(K|sdm_list|) repeated-scan loop
        with an essentially linear O(N + M·α(N)) algorithm.
        """
        n = len(all_atoms)
        for at in all_atoms:
            at.molindex = -1

        parent = list(range(n))
        rank = [0] * n

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]  # path halving
                x = parent[x]
            return x

        def union(x: int, y: int) -> None:
            rx, ry = find(x), find(y)
            if rx == ry:
                return
            if rank[rx] < rank[ry]:
                rx, ry = ry, rx
            parent[ry] = rx
            if rank[rx] == rank[ry]:
                rank[rx] += 1

        for sdm_item in self.sdm_list:
            if sdm_item.covalent:
                union(sdm_item.a1, sdm_item.a2)

        root_to_mol: dict[int, int] = {}
        mol_counter = 0
        for i in range(n):
            root = find(i)
            if root not in root_to_mol:
                mol_counter += 1
                root_to_mol[root] = mol_counter
            all_atoms[i].molindex = root_to_mol[root]

        self.maxmol = mol_counter

    def vector_length(self, x: float, y: float, z: float) -> float:
        """
        Calculates the vector length given in fractional coordinates.
        """
        A = 2.0 * (x * y * self.aga + x * z * self.bbe + y * z * self.cal)
        return sqrt(x ** 2 * self.asq + y ** 2 * self.bsq + z ** 2 * self.csq + A)

    def packer(self, sdm: 'SDM', need_symm: list, with_qpeaks=False) -> list:
        """
        Packs atoms of the asymmetric unit to real molecules.
        """
        if not self._symm_m:
            self._build_symm_arrays()
        symm_m = self._symm_m
        symm_t = self._symm_t

        asymm = self.shx.atoms.all_atoms
        if with_qpeaks:
            showatoms = asymm[:]
        else:
            showatoms = [at for at in asymm if not at.qpeak]

        # Seed used_names with all names already present so generated names
        # never clash with the asymmetric-unit atoms or with each other.
        used_names: set = {at.name.upper() for at in showatoms}

        for symm in need_symm:
            symm_num, h, k, l, symmgroup = symm  # noqa: E741
            h -= 5
            k -= 5
            l -= 5  # noqa: E741
            symm_num -= 1

            m = symm_m[symm_num]
            t = symm_t[symm_num]

            for atom in asymm:
                if not with_qpeaks and atom.qpeak:
                    continue
                if atom.molindex != symmgroup:
                    continue

                new_atom = Atom(self.shx)
                if atom.qpeak:
                    continue

                uvals = list(atom.uvals)
                # Transform anisotropic U-values so ADPs are oriented correctly
                # after the symmetry rotation.
                if sum(abs(u) for u in uvals[2:]) > 1e-5:
                    uvals = list(self.transform_uvalues(uvals, symm_num))

                # Apply symmetry operation and lattice shift
                x1, y1, z1 = atom.x, atom.y, atom.z
                px = x1 * m[0][0] + y1 * m[1][0] + z1 * m[2][0] + t[0] + h
                py = x1 * m[0][1] + y1 * m[1][1] + z1 * m[2][1] + t[1] + k
                pz = x1 * m[0][2] + y1 * m[1][2] + z1 * m[2][2] + t[2] + l

                # Assign a unique, legal SHELXL atom name (≤4 alphanumeric chars).
                element = self.shx.sfac2elem(atom.sfac_num)
                new_name = _unique_legal_atom_name(element, used_names)

                new_atom.set_atom_parameters(
                    name=new_name,
                    sfac_num=atom.sfac_num,
                    coords=[px, py, pz],
                    part=atom.part,
                    afix=AFIX(self.shx, (atom.afix).split()) if atom.afix else None,
                    resi=RESI(self.shx, (f'RESI {atom.resinum} {atom.resiclass}').split()) if atom.resi else None,
                    site_occupation=atom.sof,
                    uvals=uvals,
                    symmgen=True,
                )

                isthere = False
                if new_atom.part.n >= 0:
                    for existing in showatoms:
                        if existing.part.n != new_atom.part.n:
                            continue
                        length = sdm.vector_length(
                            new_atom.x - existing.x,
                            new_atom.y - existing.y,
                            new_atom.z - existing.z)
                        if length < 0.2:
                            isthere = True
                            break
                if not isthere:
                    showatoms.append(new_atom)
                else:
                    # Atom was a duplicate — release the name we just reserved
                    used_names.discard(new_name.upper())

        return showatoms

    def pack_unit_cell(
            self,
            symmop_indices: list[int] | None = None,
            *,
            cart_tolerance: float = 0.2,
            with_qpeaks: bool = False,
    ) -> list:
        """Pack all symmetry-equivalent positions into one unit cell.

        For every atom in the asymmetric unit every selected symmetry operation
        is applied and the result is folded back into [0, 1) fractional
        coordinates.  Positions already occupied within *cart_tolerance* Å
        (with periodic boundary conditions) are discarded as duplicates.

        This can be called on a fresh :class:`SDM` object before
        :meth:`calc_sdm` — it does **not** require the SDM to have been run.

        :param symmop_indices: 0-based indices into the symmetry-card list
            (identity is index 0).  ``None`` applies all operations.
        :param cart_tolerance: Cartesian duplicate threshold in Å (default 0.2).
        :param with_qpeaks: include Q-peaks in the output (default False).
        :returns: List of :class:`~shelxfile.atoms.atom.Atom` objects with
            updated fractional coordinates.
        """
        if not self._symm_m:
            self._build_symm_arrays()
        symm_m = self._symm_m
        symm_t = self._symm_t

        selected: list[int] = (
            list(range(len(symm_m)))
            if symmop_indices is None
            else list(symmop_indices)
        )

        a, b, c, alpha, beta, gamma = self.cell
        aga, bbe, cal = self.aga, self.bbe, self.cal
        asq, bsq, csq = self.asq, self.bsq, self.csq
        tol_sq = cart_tolerance * cart_tolerance

        # Spatial grid for O(1) average-case duplicate detection.
        # Each bin spans ~cart_tolerance; at least 3 bins per axis so the
        # 3×3×3 neighbour check covers exactly one cell-image distance.
        grid_nx = max(3, int(a / cart_tolerance))
        grid_ny = max(3, int(b / cart_tolerance))
        grid_nz = max(3, int(c / cart_tolerance))

        # grid maps (ix, iy, iz) → [(fx, fy, fz, part, packed_idx), ...]
        grid: dict[tuple[int, int, int], list[tuple]] = {}
        # packed entries: (orig_atom, px, py, pz, symm_idx)
        packed: list[tuple] = []

        asymm = self.shx.atoms.all_atoms
        for at in asymm:
            if not with_qpeaks and at.qpeak:
                continue
            x1, y1, z1 = at.x, at.y, at.z
            part = at.part.n

            for idx in selected:
                m = symm_m[idx]
                t = symm_t[idx]

                # Apply symmetry operation (column-major) and fold to [0, 1)
                px = (x1 * m[0][0] + y1 * m[1][0] + z1 * m[2][0] + t[0]) % 1.0
                py = (x1 * m[0][1] + y1 * m[1][1] + z1 * m[2][1] + t[1]) % 1.0
                pz = (x1 * m[0][2] + y1 * m[1][2] + z1 * m[2][2] + t[2]) % 1.0

                # Determine grid bin for this candidate position
                ix = int(px * grid_nx) % grid_nx
                iy = int(py * grid_ny) % grid_ny
                iz = int(pz * grid_nz) % grid_nz

                # Check 3×3×3 neighbouring bins for a duplicate
                is_dup = False
                for dix in (-1, 0, 1):
                    if is_dup:
                        break
                    nix = (ix + dix) % grid_nx
                    for diy in (-1, 0, 1):
                        if is_dup:
                            break
                        niy = (iy + diy) % grid_ny
                        for diz in (-1, 0, 1):
                            if is_dup:
                                break
                            niz = (iz + diz) % grid_nz
                            bucket = grid.get((nix, niy, niz))
                            if bucket is None:
                                continue
                            for (efx, efy, efz, epart, _) in bucket:
                                # Different non-zero disorder parts are never duplicates
                                if epart != 0 and part != 0 and epart != part:
                                    continue
                                # Fractional difference folded to [-0.5, 0.5]
                                ddx = px - efx
                                ddx -= int(ddx + (0.5 if ddx >= 0.0 else -0.5))
                                ddy = py - efy
                                ddy -= int(ddy + (0.5 if ddy >= 0.0 else -0.5))
                                ddz = pz - efz
                                ddz -= int(ddz + (0.5 if ddz >= 0.0 else -0.5))
                                d2 = (ddx * ddx * asq + ddy * ddy * bsq
                                      + ddz * ddz * csq
                                      + 2.0 * (ddx * ddy * aga
                                               + ddx * ddz * bbe
                                               + ddy * ddz * cal))
                                if d2 < tol_sq:
                                    is_dup = True
                                    break

                if not is_dup:
                    idx_packed = len(packed)
                    packed.append((at, px, py, pz, idx))
                    key = (ix, iy, iz)
                    bucket = grid.get(key)
                    if bucket is None:
                        grid[key] = [(px, py, pz, part, idx_packed)]
                    else:
                        bucket.append((px, py, pz, part, idx_packed))

        # Build Atom objects with the packed fractional coordinates.
        # Assign unique, legal SHELXL names: identity-operation atoms keep
        # their original names; symmetry-generated copies get new sequential
        # names so no two atoms in the result share the same label.
        result: list[Atom] = []
        used_names: set = set()
        for (orig_at, px, py, pz, symm_num) in packed:
            new_atom = Atom(self.shx)
            uvals = list(orig_at.uvals)
            if sum(abs(u) for u in uvals[2:]) > 1e-5:
                uvals = list(self.transform_uvalues(uvals, symm_num))
            # For the identity operation the original name is available unless
            # already taken (shouldn't happen, but guard anyway).
            orig_name_up = orig_at.name.upper()
            if orig_name_up not in used_names:
                atom_name = orig_at.name
                used_names.add(orig_name_up)
            else:
                element = self.shx.sfac2elem(orig_at.sfac_num)
                atom_name = _unique_legal_atom_name(element, used_names)
            new_atom.set_atom_parameters(
                name=atom_name,
                sfac_num=orig_at.sfac_num,
                coords=[px, py, pz],
                part=orig_at.part,
                afix=AFIX(self.shx, orig_at.afix.split()) if orig_at.afix else None,
                resi=RESI(self.shx, f'RESI {orig_at.resinum} {orig_at.resiclass}'.split())
                if orig_at.resi else None,
                site_occupation=orig_at.sof,
                uvals=uvals,
                symmgen=(symm_num != 0),
            )
            result.append(new_atom)

        return result

    def transform_uvalues(self, uvals: list | tuple, symm_num: int) -> tuple:
        """
        Transforms the Uij values according to the given symmetry rotation.

        Reference: R. W. Grosse-Kunstleve, P. D. Adams (2002). J. Appl. Cryst. 35, 477–480.
        http://dx.doi.org/10.1107/S0021889802008580

        Steps:
          U(star) = N  @ U(cif) @ N.T      (CIF → star parameterisation)
          U(star) = R  @ U(star) @ R.T     (symmetry rotation)
          U(cif)  = N⁻¹ @ U(star) @ N⁻¹.T  (star → CIF)

        where N = diag(a*, b*, c*) and R = symmcards[symm_num].matrix (numpy, row-major).
        """
        U11, U22, U33, U23, U13, U12 = uvals
        Ucif = np.array([[U11, U12, U13],
                         [U12, U22, U23],
                         [U13, U23, U33]], dtype=float)
        N = np.diag([self.astar, self.bstar, self.cstar])
        N_inv = np.linalg.inv(N)
        R = self.shx.symmcards[symm_num].matrix  # numpy (3,3), not transposed

        Ustar = N @ Ucif @ N.T
        Ustar = R @ Ustar @ R.T
        Ucif_new = N_inv @ Ustar @ N_inv.T

        return (float(Ucif_new[0, 0]), float(Ucif_new[1, 1]), float(Ucif_new[2, 2]),
                float(Ucif_new[1, 2]), float(Ucif_new[0, 2]), float(Ucif_new[0, 1]))


def ufrac_to_ucart(A, cell: tuple, uvals: list) -> np.ndarray:
    """
    Converts anisotropic displacement parameters from fractional (CIF) to
    Cartesian coordinates.

    :param A:     OrthogonalMatrix (or any object with a .values attribute)
    :param cell:  (a, b, c, alpha, beta, gamma)
    :param uvals: [U11, U22, U33, U23, U13, U12] in fractional parameterisation
    :returns:     3×3 numpy array of U(cart)
    """
    U11, U22, U33, U23, U13, U12 = uvals
    Uij = np.array([[U11, U12, U13],
                    [U12, U22, U23],
                    [U13, U23, U33]], dtype=float)
    a, b, c, alpha, beta, gamma = cell
    V = vol_unitcell(*cell)
    astar = (b * c * sin(radians(alpha))) / V
    bstar = (c * a * sin(radians(beta))) / V
    cstar = (a * b * sin(radians(gamma))) / V
    N = np.diag([astar, bstar, cstar])
    A_np = np.array(A.values, dtype=float)
    return A_np @ N @ Uij @ N.T @ A_np.T


if __name__ == "__main__":
    from shelxfile.shelx.shelx import Shelxfile

    shx = Shelxfile()
    shx.read_file('tests/resources/p-31c.res')
    t1 = time.perf_counter()
    sdm = SDM(shx)
    needsymm = sdm.calc_sdm()
    print('sdmlist:', len(sdm.sdm_list))
    print(needsymm)
    packed_atoms = sdm.packer(sdm, needsymm)
    print('Zeit für sdm:', round(time.perf_counter() - t1, 3), 's')

    head = """
REM Solution 1  R1  0.081,  Alpha = 0.0146  in P31c
REM Flack x = -0.072 ( 0.041 ) from Parsons' quotients
REM C19.667 N2.667 P4
TITL p-31c-neu in P-31c
CELL  0.71073  12.5067  12.5067  24.5615   90.000   90.000  120.000
ZERR   2   0.0043   0.0043   0.0085   0.000   0.000   0.000
LATT -1
SFAC C H N P Cl
UNIT 120 186 14 12 12
TEMP -173.000
OMIT 0   0   2
L.S. 10
BOND $H
ACTA
CONF

LIST 4
FMAP 2
PLAN 40
WGHT    0.034600    0.643600
FVAR       0.22604   0.76052   0.85152\n"""

    tail = """\n
HKLF 4
 
REM  p-31c-neu in P-31c
REM R1 =  0.0308 for    4999 Fo > 4sig(Fo)  and  0.0343 for all    5352 data
REM    287 parameters refined using    365 restraints
 
END 
 
WGHT      0.0348      0.6278 
"""
    for at in packed_atoms:
        if at.qpeak:
            continue
        head += wrap_line(str(at)) + '\n'
    head += tail
    p = Path('./test.res')
    p.write_text(head)
    print(len(shx.atoms))
    print(len(packed_atoms))
    print('Zeit für sdm:', round(sdm.sdmtime, 3), 's')
