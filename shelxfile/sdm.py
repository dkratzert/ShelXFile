# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
import time
from math import sqrt, radians, sin
from pathlib import Path
from string import ascii_letters

from atoms import Atom
from cards import AFIX, RESI
from dsrmath import Array, Matrix, vol_unitcell
from misc import DEBUG, wrap_line


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
        return '{} {} {} {} dist: {} coval: {} sn: {} {}'.format(self.atom1.name, self.atom2.name, self.a1, self.a2,
                                                                 self.dist, self.covalent,
                                                                 self.symmetry_number, self.dddd)


class SDM():
    """
    This class calculates the shortest distance matrix and creates a completed (grown) structure by crystal symmetry.
    """

    def __init__(self, shx: 'Shelxfile'):
        self.shx = shx
        self.aga = self.shx.cell.a * self.shx.cell.b * self.shx.cell.cosga
        self.bbe = self.shx.cell.a * self.shx.cell.c * self.shx.cell.cosbe
        self.cal = self.shx.cell.b * self.shx.cell.c * self.shx.cell.cosal
        self.sdm_list = []  # list of sdmitems
        self.maxmol = 1
        self.sdmtime = 0
        self.bondlist = []
        self.asq = self.shx.cell[0] ** 2
        self.bsq = self.shx.cell[1] ** 2
        self.csq = self.shx.cell[2] ** 2
        # calculate reciprocal lattice vectors:
        self.astar = (self.shx.cell.b * self.shx.cell.c * sin(radians(self.shx.cell.al))) / self.shx.cell.V
        self.bstar = (self.shx.cell.c * self.shx.cell.a * sin(radians(self.shx.cell.be))) / self.shx.cell.V
        self.cstar = (self.shx.cell.a * self.shx.cell.b * sin(radians(self.shx.cell.ga))) / self.shx.cell.V

    def calc_sdm(self) -> list:
        t1 = time.perf_counter()
        all_atoms = self.shx.atoms.all_atoms
        self.bondlist.clear()
        for i, at1 in enumerate(all_atoms):
            prime_array = [Array(at1.frac_coords) * symop.matrix + symop.trans for symop in self.shx.symmcards]
            for j, at2 in enumerate(all_atoms):
                mind = 1000000
                hma = False
                at2_plushalf = Array(at2.frac_coords) + 0.5
                sdmItem = SDMItem()
                for n, symop in enumerate(self.shx.symmcards):
                    D = prime_array[n] - at2_plushalf
                    dp = [v - 0.5 for v in D - D.floor]
                    dk = self.vector_length(*dp)
                    if dk > 5:
                        continue
                    if n:
                        dk += 0.0001
                    if (dk > 0.01) and (mind >= dk):
                        mind = min(dk, mind)
                        sdmItem.dist = mind
                        sdmItem.atom1 = at1
                        sdmItem.atom2 = at2
                        sdmItem.a1 = i
                        sdmItem.a2 = j
                        sdmItem.symmetry_number = n
                        hma = True
                if not sdmItem.atom1:
                    # Do not grow grown atoms:
                    continue
                if (not sdmItem.atom1.ishydrogen and not sdmItem.atom2.ishydrogen) and \
                        sdmItem.atom1.part.n * sdmItem.atom2.part.n == 0 \
                        or sdmItem.atom1.part.n == sdmItem.atom2.part.n:
                    dddd = (at1.radius + at2.radius) * 1.2
                    sdmItem.dddd = dddd
                else:
                    dddd = 0.0
                if sdmItem.dist < dddd:
                    if hma:
                        #self.bondlist.append((i, j, sdmItem.atom1.name, sdmItem.atom2.name, sdmItem.dist))
                        sdmItem.covalent = True
                else:
                    sdmItem.covalent = False
                if hma:
                    self.sdm_list.append(sdmItem)
        t2 = time.perf_counter()
        self.sdmtime = t2 - t1
        if DEBUG:
            print('Zeit sdm_calc:', self.sdmtime)
        self.sdm_list.sort()
        self.calc_molindex(all_atoms)
        need_symm = self.collect_needed_symmetry()
        if DEBUG:
            print("The asymmetric unit contains {} fragments.".format(self.maxmol))
        return need_symm

    def collect_needed_symmetry(self) -> list:
        need_symm = []
        # Collect needsymm list:
        for sdmItem in self.sdm_list:
            if sdmItem.covalent:
                if sdmItem.atom1.molindex < 1 or sdmItem.atom1.molindex > 6:
                    continue
                for n, symop in enumerate(self.shx.symmcards):
                    if sdmItem.atom1.part.n != 0 and sdmItem.atom2.part.n != 0 \
                            and sdmItem.atom1.part.n != sdmItem.atom2.part.n:
                        # both not part 0 and different part numbers
                        continue
                    # Both the same atomic number and number hydrogen:
                    if sdmItem.atom1.an == sdmItem.atom2.an and sdmItem.atom1.ishydrogen:
                        continue
                    prime = Array(sdmItem.atom1.frac_coords) * symop.matrix + symop.trans
                    D = prime - Array(sdmItem.atom2.frac_coords) + Array([0.5, 0.5, 0.5])
                    floorD = D.floor
                    dp = D - floorD - Array([0.5, 0.5, 0.5])
                    if n == 0 and Array([0, 0, 0]) == floorD:
                        continue
                    dk = self.vector_length(*dp)
                    dddd = sdmItem.dist + 0.2
                    if sdmItem.atom1.ishydrogen and sdmItem.atom2.ishydrogen:
                        dddd = 1.8
                    if (dk > 0.001) and (dddd >= dk):
                        bs = [n + 1, (5 - floorD[0]), (5 - floorD[1]), (5 - floorD[2]), sdmItem.atom1.molindex]
                        # Does not work:
                        # self.bondlist.append((sdmItem.a1, sdmItem.a2, sdmItem.atom1.name+'<',
                        #                      sdmItem.atom2.name+'<', sdmItem.dist))
                        if bs not in need_symm:
                            need_symm.append(bs)
        return need_symm

    def calc_molindex(self, all_atoms):
        # Start for George's "bring atoms together algorithm":
        someleft = 1
        nextmol = 1
        for at in all_atoms:
            at.molindex = -1
        all_atoms[0].molindex = 1
        while nextmol:
            someleft = 1
            nextmol = 0
            while someleft:
                someleft = 0
                for sdmItem in self.sdm_list:
                    if sdmItem.covalent and sdmItem.atom1.molindex * sdmItem.atom2.molindex < 0:
                        sdmItem.atom1.molindex = self.maxmol
                        sdmItem.atom2.molindex = self.maxmol
                        someleft += 1
            for ni, at in enumerate(all_atoms):
                if not at.ishydrogen and at.molindex < 0:
                    nextmol = ni
                    break
            if nextmol:
                self.maxmol += 1
                all_atoms[nextmol].molindex = self.maxmol

    def vector_length(self, x: float, y: float, z: float) -> float:
        """
        Calculates the vector length given in fractional coordinates.
        """
        A = 2.0 * (x * y * self.aga + x * z * self.bbe + y * z * self.cal)
        return sqrt(x ** 2 * self.asq + y ** 2 * self.bsq + z ** 2 * self.csq + A)

    def packer(self, sdm: 'SDM', need_symm: list, with_qpeaks=False):
        """
        Packs atoms of the asymmetric unit to real molecules.
        """
        # t1 = time.perf_counter()
        asymm = self.shx.atoms.all_atoms
        if with_qpeaks:
            showatoms = asymm[:]
        else:
            showatoms = [at for at in asymm if not at.qpeak]
        for symm in need_symm:
            symm_num, h, k, l, symmgroup = symm
            h -= 5
            k -= 5
            l -= 5
            symm_num -= 1
            for atom in asymm:
                if not with_qpeaks:
                    if atom.qpeak:
                        continue
                # Essential to have hydrogen atoms grown:
                # if not atom.ishydrogen and atom.molindex == symmgroup:
                if atom.molindex == symmgroup:
                    new_atom = Atom(self.shx)
                    if atom.qpeak:
                        continue
                    else:
                        pass
                        uvals = atom.uvals
                        # TODO: Transform u values according to symmetry:
                        # currently, the adps are directed in wrong directions after after applying symmetry to atoms.
                        # uvals = self.transform_uvalues(uvals, symm_num)
                    new_atom.set_atom_parameters(
                        name=atom.name[:3] + ">>" + str(symm_num) + '_' + ascii_letters[atom.part.n],
                        sfac_num=atom.sfac_num,
                        coords=Array(atom.frac_coords) * self.shx.symmcards[symm_num].matrix
                               + Array(self.shx.symmcards[symm_num].trans) + Array([h, k, l]),
                        part=atom.part,
                        afix=AFIX(self.shx, (atom.afix).split()) if atom.afix else None,
                        resi=RESI(self.shx, ('RESI ' + atom.resinum + atom.resiclass).split()) if atom.resi else None,
                        site_occupation=atom.sof,
                        uvals=uvals,
                        symmgen=True
                    )
                    isthere = False
                    if new_atom.part.n >= 0:
                        for atom in showatoms:
                            if atom.part.n != new_atom.part.n:
                                continue
                            length = sdm.vector_length(new_atom.frac_coords[0] - atom.frac_coords[0],
                                                       new_atom.frac_coords[1] - atom.frac_coords[1],
                                                       new_atom.frac_coords[2] - atom.frac_coords[2])
                            if length < 0.2:
                                isthere = True
                    if not isthere:
                        showatoms.append(new_atom)
                # elif grow_qpeaks:
                #    add q-peaks here
        # t2 = time.perf_counter()
        # print('packzeit:', t2-t1) # 0.04s
        return showatoms

    def transform_uvalues(self, uvals: (list, tuple), symm_num: int):
        """
        Transforms the Uij values according to local symmetry.
        R. W. Grosse-Kunstleve, P. D. Adams (2002). J. Appl. Cryst. 35, 477–480.
        http://dx.doi.org/10.1107/S0021889802008580

        U(star) = N * U(cif) * N.T
        U(star) = R * U(star) * R^t
        U(cif) = N^-1 * U(star) * (N^-1).T

        U(cart) = A * U(star) * A.T
        U(frac) = A^1 * U(cart) * (A^1).t
        U(star) = A^-1 * U(cart) * A^-1.t

        R is the rotation part of a given symmetry operation

        [ [U11, U12, U13]
          [U12, U22, U23]
          [U13, U23, U33]
        ]
        # Shelxl uses U* with a*,b*c*-parameterization
        atomname sfac x y z sof[11] U[0.05] or U11 U22 U33 U23 U13 U12

        Read:
        X-Ray Analysis and the Structure of Organic Molecules Second Edition, Dunitz, P240

        """
        U11, U22, U33, U23, U13, U12 = uvals
        U21 = U12
        U32 = U23
        U31 = U13
        Ucif = Matrix([[U11, U12, U13], [U21, U22, U23], [U31, U32, U33]])
        # matrix with the reciprocal lattice vectors:
        N = Matrix([[self.astar, 0, 0],
                    [0, self.bstar, 0],
                    [0, 0, self.cstar]])
        R = self.shx.symmcards[symm_num].matrix
        R_t = self.shx.symmcards[symm_num].matrix.transposed
        A = self.shx.orthogonal_matrix
        # U(star) = N * U(cif) * N.T
        # U(cart) = A * U(star) * A.T
        # U(star) = R * U(star) * R^t
        # U(cif) = N^-1 * U(star) * (N^-1).T
        # U(star) = A^-1 * U(cart) * A^-1.T
        Ustar = N * Ucif * N.T
        Ustar = R * Ustar * R_t
        Ucif = N.inversed * Ustar * N.inversed.T
        uvals = Ucif
        upper_diagonal = uvals.values[0][0], uvals.values[1][1], uvals.values[2][2], \
                         uvals.values[1][2], uvals.values[0][2], \
                         uvals.values[0][1]
        return upper_diagonal


def ufrac_to_ucart(A, cell, uvals):
    """
    #>>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    #>>> Ucart = Matrix([[ 0.0754483395556807,  0.030981701122469,  -0.0194466522033868], [  0.030981701122469,  0.03058, -0.01057], [-0.0194466522033868, -0.01057, 0.03216] ])
    #>>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    #>>> from shelxfile.dsrmath import OrthogonalMatrix
    #>>> A = OrthogonalMatrix(*cell)
    #>>> ufrac_to_ucart(A, cell, uvals)
    #TODO: test if this is right:
    | 0.0740  0.0310 -0.0299|
    | 0.0302  0.0306 -0.0148|
    |-0.0171 -0.0106  0.0346|
    <BLANKLINE>
    """
    U11, U22, U33, U23, U13, U12 = uvals
    U21 = U12
    U32 = U23
    U31 = U13
    Uij = Matrix([[U11, U12, U13], [U21, U22, U23], [U31, U32, U33]])
    a, b, c, alpha, beta, gamma = cell
    V = vol_unitcell(*cell)
    # calculate reciprocal lattice vectors:
    astar = (b * c * sin(radians(alpha))) / V
    bstar = (c * a * sin(radians(beta))) / V
    cstar = (a * b * sin(radians(gamma))) / V
    # matrix with the reciprocal lattice vectors:
    N = Matrix([[astar, 0, 0],
                [0, bstar, 0],
                [0, 0, cstar]])
    # Finally transform Uij values from fractional to cartesian axis system:
    Ucart = A * N * Uij * N.T * A.T
    return Ucart


if __name__ == "__main__":
    from .shelx import Shelxfile

    shx = Shelxfile('tests/p-31c.res')
    sdm = SDM(shx)
    needsymm = sdm.calc_sdm()
    packed_atoms = sdm.packer(sdm, needsymm)
    # print(needsymm)
    # [[8, 5, 5, 5, 1], [16, 5, 5, 5, 1], [7, 4, 5, 5, 3]]
    # print(len(shx.atoms))
    # print(len(packed_atoms))

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
        print(wrap_line(str(at)))
        head += wrap_line(str(at)) + '\n'
    head += tail
    p = Path('./test.res')
    p.write_text(head)
    print('Zeit für sdm:', round(sdm.sdmtime, 3), 's')
    # print(sdm.bondlist)
    print(len(sdm.bondlist), '(170) Atome in p-31c.res, (208) in I-43d.res')
