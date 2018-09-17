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
from math import sqrt

from shelxfile.atoms import Atom
from shelxfile.cards import AFIX, RESI
from shelxfile.dsrmath import Array
from shelxfile.misc import DEBUG, wrap_line


class SDMItem(object):
    def __init__(self):
        self.dist = 0.0
        self.atom1 = None
        self.a1 = 0
        self.atom2 = None
        self.a2 = 0
        self.symmetry_number = 0
        self.floorD = 0
        self.covalent = True
        self.dddd = 0

    def __lt__(self, a2):
        return True if self.dist < a2.dist else False

    def __eq__(self, other: 'SDMItem'):
        if other.a1 == self.a2 and other.a2 == self.a1:
            return True
        return False

    def __repr__(self):
        return '{} {} {} {} dist: {} coval: {} sn: {} {} {}'.format(self.atom1.name, self.atom2.name, self.a1, self.a2,
                                                                    self.dist, self.covalent,
                                                                    self.symmetry_number, self.floorD, self.dddd)


class SDM():
    def __init__(self, shx: 'ShelXFile'):
        self.shx = shx
        self.aga = self.shx.cell[0] * self.shx.cell[1] * self.shx.cell.cosga
        self.bbe = self.shx.cell[0] * self.shx.cell[2] * self.shx.cell.cosbe
        self.cal = self.shx.cell[1] * self.shx.cell[2] * self.shx.cell.cosal
        self.dk, dddd = 0.0, 0.0
        self.prime = None
        self.dp = None
        self.D = None
        self.floorD = None
        # contact.clear()  # sdmitem list for hydrogen contacts (not needed)?
        self.sdm_list = []  # list of sdmitems
        self.knots = []
        self.maxmol = 1
        self.sdmtime = 0

    def __iter__(self):
        yield self.knots

    def calc_sdm(self) -> list:
        t1 = time.perf_counter()
        all_atoms = self.shx.atoms.all_atoms
        for i, at1 in enumerate(all_atoms):
            atneighb = []  # list of atom neigbors
            for j, at2 in enumerate(all_atoms):
                mind = 1000000
                hma = False
                sdmItem = SDMItem()
                for n, symop in enumerate(self.shx.symmcards):
                    prime = Array(at1.frac_coords) * symop.matrix + symop.trans
                    D = prime - Array(at2.frac_coords) + Array([0.5, 0.5, 0.5])
                    dp = D - D.floor - Array([0.5, 0.5, 0.5])
                    dk = self.vector_length(*dp)
                    if n:
                        dk += 0.0001
                    if (dk > 0.01) and (mind >= dk):
                        mind = min(dk, mind)
                        sdmItem.dist = mind
                        sdmItem.floorD = D.floor
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
                        atneighb.append(j)
                        sdmItem.covalent = True
                else:
                    sdmItem.covalent = False
                if hma:
                    self.sdm_list.append(sdmItem)
                # print(sdmItem)
            self.knots.append(atneighb)
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
                    if sdmItem.atom1.part != 0 and sdmItem.atom2.part != 0 \
                            and sdmItem.atom1.part != sdmItem.atom2.part:
                        # both not part 0 and different part numbers
                        continue
                    # Both the same atomic number and number 0 (hydrogen)
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
                    # TODO: Do I need this?
                    if sdmItem.atom1.ishydrogen and sdmItem.atom2.ishydrogen:
                        dddd = 1.8
                    if (dk > 0.001) and (dddd >= dk):
                        bs = [n + 1, (5 - floorD[0]), (5 - int(floorD[1])), (5 - int(floorD[2])),
                              sdmItem.atom1.molindex]
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
        a = 2.0 * x * y * self.aga
        b = 2.0 * x * z * self.bbe
        c = 2.0 * y * z * self.cal
        return sqrt(x ** 2 * self.shx.cell[0] ** 2 + y ** 2 * self.shx.cell[1] ** 2
                    + z ** 2 * self.shx.cell[2] ** 2 + a + b + c)

    def packer(self, sdm: 'SDM', need_symm: list, with_qpeaks=False):
        """
        Packs atoms of the asymmetric unit to real molecules.
        TODO: Support hydrogen atoms!
        """
        asymm = self.shx.atoms.all_atoms
        if with_qpeaks:
            showatoms = asymm[:]
        else:
            showatoms = [at for at in asymm if not at.qpeak]
        for symm in need_symm:
            s, h, k, l, symmgroup = symm
            h -= 5
            k -= 5
            l -= 5
            s -= 1
            for atom in asymm:
                if not with_qpeaks:
                    if atom.qpeak:
                        continue
                if not atom.ishydrogen and atom.molindex == symmgroup:
                #if atom.molindex == symmgroup:
                    new_atom = Atom(self.shx)
                    new_atom.set_atom_parameters(
                        name=atom.name + ">" + 'a',
                        sfac_num=atom.sfac_num,
                        coords=Array(atom.frac_coords) * self.shx.symmcards[s].matrix
                                + Array(self.shx.symmcards[s].trans) + Array([h, k, l]),
                        part=atom.part,
                        afix=AFIX(self.shx, ('AFIX ' + atom.afix).split()) if atom.afix else None,
                        resi=RESI(self.shx, ('RESI ' + atom.resinum + atom.resiclass).split()) if atom.resi else None,
                        site_occupation=atom.sof,
                        uvals=atom.uvals,
                        symmgen=True
                    )
                    # TODO: I have to transform the Uijs by symmetry here later.
                    isthere = False
                    if new_atom.part.n >= 0:
                        for atom in showatoms:
                            if atom.ishydrogen:
                                continue
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
        return showatoms


if __name__ == "__main__":
    from shelxfile.shelx import ShelXFile
    shx = ShelXFile('tests/I-43d.res')
    sdm = SDM(shx)
    needsymm = sdm.calc_sdm()
    packed_atoms = sdm.packer(sdm, needsymm)
    #print(needsymm)
    # [[8, 5, 5, 5, 1], [16, 5, 5, 5, 1], [7, 4, 5, 5, 3]]
    #print(len(shx.atoms))
    #print(len(packed_atoms))

    for at in packed_atoms:
        if at.qpeak:
            continue
        print(wrap_line(str(at)))

    print('Zeit für sdm:', round(sdm.sdmtime, 3), 's')