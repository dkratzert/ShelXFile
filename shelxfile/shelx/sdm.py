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
from string import ascii_letters

import numpy as np

from shelxfile.atoms.atom import Atom
from shelxfile.misc.dsrmath import vol_unitcell
from shelxfile.misc.misc import wrap_line
from shelxfile.shelx.cards import AFIX, RESI


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
        at2_plushalf = [(x + 0.5, y + 0.5, z + 0.5) for (x, y, z) in coords]
        radii = [at.radius for at in all_atoms]
        is_h = [at.ishydrogen for at in all_atoms]
        parts = [at.part.n for at in all_atoms]

        aga, bbe, cal = self.aga, self.bbe, self.cal
        asq, bsq, csq = self.asq, self.bsq, self.csq

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

                    if dk2 > 16.0:      # 4 Å hard cutoff (squared) – skip sqrt
                        continue

                    dk = sqrt(dk2)
                    if n:
                        dk += 0.0001   # slight penalty for symmetry-generated images

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
            print('Zeit sdm_calc:', self.sdmtime)

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

                new_atom.set_atom_parameters(
                    name=atom.name[:3] + ">>" + str(symm_num) + '_' + ascii_letters[atom.part.n],
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

        return showatoms

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
        R = self.shx.symmcards[symm_num].matrix   # numpy (3,3), not transposed

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
