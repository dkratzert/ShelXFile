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
from dataclasses import dataclass
from math import sqrt, radians, sin
from string import ascii_letters
from typing import Union, Tuple

# noinspection PyUnresolvedReferences
from shelxfile.shelx.sdm_complete import calc_sdm

from shelxfile.atoms.atom import Atom
from shelxfile.misc.dsrmath import Array, Matrix, vol_unitcell
from shelxfile.misc.misc import DEBUG
from shelxfile.shelx.cards import AFIX, RESI


@dataclass
class RAtom:
    name: str
    x: float
    y: float
    z: float
    element: str
    part: int
    symmgen: bool
    molindex: int
    qpeak: bool
    radius: float


class SDMR():
    """
    This class calculates the shortest distance matrix and creates a completed (grown) structure by crystal symmetry.
    """

    def __init__(self, shx: 'Shelxfile'):
        self.shx = shx
        # self.sdm_list: List[SDMItem]
        self.all_atoms = self.get_atoms()
        self.sdm_list = []  # list of sdmitems
        self.maxmol: int = 1
        self.cell = (
            self.shx.cell.a, self.shx.cell.b, self.shx.cell.c, self.shx.cell.alpha, self.shx.cell.beta, self.shx.cell.gamma)
        self.aga = self.shx.cell.a * self.shx.cell.b * self.shx.cell.cosga
        self.bbe = self.shx.cell.a * self.shx.cell.c * self.shx.cell.cosbe
        self.cal = self.shx.cell.b * self.shx.cell.c * self.shx.cell.cosal
        self.maxmol = 1
        self.sdmtime = 0
        self.bondlist = []
        self.asq = self.shx.cell[0] ** 2
        self.bsq = self.shx.cell[1] ** 2
        self.csq = self.shx.cell[2] ** 2
        # calculate reciprocal lattice vectors:
        self.astar = (self.shx.cell.b * self.shx.cell.c * sin(radians(self.shx.cell.alpha))) / self.shx.cell.V
        self.bstar = (self.shx.cell.c * self.shx.cell.a * sin(radians(self.shx.cell.beta))) / self.shx.cell.V
        self.cstar = (self.shx.cell.a * self.shx.cell.b * sin(radians(self.shx.cell.gamma))) / self.shx.cell.V

    def get_atoms(self) -> Tuple[Union[RAtom, RAtom], ...]:
        atoms = []
        for atom in self.shx.atoms:
            atoms.append(
                RAtom(name=atom.name, x=atom.x, y=atom.y, z=atom.z, element=atom.element, part=atom.part.n,
                      symmgen=False, molindex=0, qpeak=atom.qpeak, radius=atom.radius))
        return tuple(atoms)

    def calc_sdm(self) -> list:
        t1 = time.perf_counter()
        symms = self.shx.symmcards._symmcards
        self.sdm_list = calc_sdm(self.all_atoms, symms, self.shx.cell)
        self.sdm_list.sort()
        t2 = time.perf_counter()
        print('Zeit für sdm:', round(t2 - t1, 3), 's')
        self.sdmtime = t2 - t1
        self.calc_molindex(self.shx.atoms.all_atoms)
        need_symm = self.collect_needed_symmetry(self.shx.atoms.all_atoms)
        if DEBUG:
            print("The asymmetric unit contains {} fragments.".format(self.maxmol))
        return need_symm

    def collect_needed_symmetry(self, all_atoms) -> list:
        need_symm = []
        # Collect needsymm list:
        for sdm_item in self.sdm_list:
            if sdm_item.covalent:
                at1: Atom = all_atoms[sdm_item.a1]
                at2: Atom = all_atoms[sdm_item.a2]
                if at1.molindex < 1 or at1.molindex > 6:
                    continue
                prime_array = [Array(at1.frac_coords) * symop.matrix + symop.trans for symop in self.shx.symmcards]
                for n, symop in enumerate(self.shx.symmcards):
                    if at1.part != 0 and at2.part != 0 and at1.part != at2.part:
                        # both not part 0 and different part numbers
                        continue
                    # Both the same atomic number and number hydrogen:
                    if at1.element == at2.element and at1.is_hydrogen:
                        continue
                    D = prime_array[n] - Array([at2.x, at2.y, at2.z]) + Array((0.5, 0.5, 0.5))
                    floor_d = D.floor
                    if n == 0 and D.floor.values == (0, 0, 0):
                        continue
                    dp = D - floor_d - Array((0.5, 0.5, 0.5))
                    dk = self.vector_length(*dp)
                    dddd = sdm_item.dist + 0.2
                    if at1.is_hydrogen and at2.is_hydrogen:
                        dddd = 1.8
                    if (dk > 0.001) and (dddd >= dk):
                        bs = [n + 1, (5 - floor_d[0]), (5 - floor_d[1]), (5 - floor_d[2]), at1.molindex]
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
                for sdm_item in self.sdm_list:
                    at1: Atom = all_atoms[sdm_item.a1]
                    at2: Atom = all_atoms[sdm_item.a2]
                    if sdm_item.covalent and at1.molindex * at2.molindex < 0:
                        at1.molindex = self.maxmol
                        at2.molindex = self.maxmol
                        someleft += 1
            for ni, at in enumerate(all_atoms):
                if not at.is_hydrogen and at.molindex < 0:
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

    def packer(self, sdm: 'SDMR', need_symm: list, with_qpeaks=False):
        """
        Packs atoms of the asymmetric unit to real molecules.
        """
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
                if not with_qpeaks and atom.qpeak:
                    continue
                # Essential to have hydrogen atoms grown:
                # if not atom.ishydrogen and atom.molindex == symmgroup:
                if atom.molindex == symmgroup:
                    new_atom = Atom(self.shx)
                    if atom.qpeak:
                        continue
                    else:
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
                        resi=RESI(self.shx, ('RESI ' + str(atom.resinum) + atom.resiclass).split()) if atom.resi else None,
                        site_occupation=atom.sof,
                        uvals=uvals,
                        symmgen=True
                    )
                    isthere = False
                    if new_atom.part.n >= 0:
                        for atom in showatoms:
                            if atom.part.n != new_atom.part.n:
                                continue
                            length = sdm.vector_length(new_atom.x - atom.x, new_atom.y - atom.y, new_atom.z - atom.z)
                            if length < 0.2:
                                isthere = True
                    if not isthere:
                        showatoms.append(new_atom)
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
    from shelxfile.shelx.shelx import Shelxfile

    shx = Shelxfile()
    shx.read_file('shelxfile/tests/resources/p-31c.res')
    t1 = time.perf_counter()
    sdm = SDMR(shx)
    needsymm = sdm.calc_sdm()
    print('Zeit für sdm:', round(time.perf_counter() - t1, 3), 's')
    print(needsymm)
    packed_atoms = sdm.packer(sdm, needsymm)
    print(len(shx.atoms))
    print(len(packed_atoms))
    assert str(packed_atoms[90]) == 'H2>>1_b 2   0.557744    0.080938   0.300634   21.00000   -1.30000'
    assert str(packed_atoms[
                   129]) == 'C15>>2_a  1    1.145639    0.497216    0.299794    11.00000    0.01803    0.01661      0.03458    0.00038   -0.00408    0.01187'
    # [[1, 5, 5, 5, 4], [1, 5, 5, 5, 5], [2, 6, 5, 5, 5], [3, 6, 6, 5, 5], [1, 5, 5, 5, 3], [2, 6, 6, 5, 3],
    # 88
    # 208

