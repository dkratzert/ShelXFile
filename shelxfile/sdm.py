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
from shelxfile.atoms import Atom
from math import sqrt

from shelxfile.dsrmath import fmin, Array
from shelxfile.shelx import ShelXFile
from shelxfile.cards import AFIX, RESI


class SDMItem(object):
    def __init__(self):
        self.dist = 0.0
        self.atom1 = None
        self.a1 = 0
        self.atom2 = None
        self.a2 = 0
        self.symmetry_number = 0
        self.floor_dist = None
        self.covalent = True

    def __lt__(self, a2):
        return True if self.dist < a2.dist else False


class SDM():
    def __init__(self, shx):
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
        self.sdm = []  # list of sdmitems
        self.knots = []

    def calc_sdm(self):
        need_symm = []
        t1 = time.perf_counter()
        for n, at1 in enumerate(self.shx.atoms.all_atoms):
            atneighb = []  # list of atom neigbors
            for j, at2 in enumerate(self.shx.atoms):
                min = 1000000
                hma = False
                sdmItem = SDMItem()
                for n, symop in enumerate(self.shx.symmcards):
                    prime = Array(at1.frac_coords) * self.shx.symmcards[n].matrix + self.shx.symmcards[n].trans
                    D = prime - Array(at2.frac_coords) + Array([0.5, 0.5, 0.5])
                    dp = D - D.floor - Array([0.5, 0.5, 0.5])
                    dk = self.vector_length(*dp)
                    if n:
                        dk += 0.0001
                    if (dk > 0.01) and (min >= dk):
                        min = fmin(dk, min)
                        sdmItem.dist = min
                        sdmItem.floorD = D.floor
                        sdmItem.atom1 = at1
                        sdmItem.atom2 = at2
                        sdmItem.a1 = n
                        sdmItem.a2 = j
                        sdmItem.symmetry_number = n
                        hma = True
                if (not sdmItem.atom1.ishydrogen and not sdmItem.atom2.ishydrogen) and \
                        sdmItem.atom1.part.n * sdmItem.atom2.part.n == 0 \
                        or sdmItem.atom1.part.n == sdmItem.atom2.part.n:
                    dddd = (at1.radius + at2.radius) * 1.2
                else:
                    dddd = 0.0
                if sdmItem.dist < dddd:
                    if hma:
                        atneighb.append(j)
                        sdmItem.covalent = True
                else:
                    sdmItem.covalent = False
                if hma:
                   self.sdm.append(sdmItem)
            self.knots.append(atneighb)
        t2 = time.perf_counter()
        print('Zeit:', t2-t1)
        print('länge:', len(self.knots), self.knots)
        # return
        t3 = time.perf_counter()
        for sdmItem in self.sdm:
            if sdmItem.covalent:
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
                        bs = [n+1, (5-floorD[0]), (5-int(floorD[1])), (5-int(floorD[2]))]
                        if not bs in need_symm:
                            need_symm.append(bs)
                        #print(symop)
                        #print(bs)
        for x in need_symm:
            print(x)
        t4 = time.perf_counter()
        print('Zeit2 brauchsymm:', t4 - t3)
        t5 = time.perf_counter()
        self.sdm.sort()
        t6 = time.perf_counter()
        print('Zeit3 sort sdm:', t6 - t5)
        t7 = time.perf_counter()
        flags = []
        for i, atom in enumerate(self.shx.atoms.all_atoms):
            flags.append(-1 if atom.an == 0 else 1)
        for i, sdmitem in enumerate(self.sdm):
            if sdmitem.atom1.ishydrogen or sdmitem.atom2.ishydrogen:
                continue
            if flags[sdmitem.a1] * flags[sdmitem.a2] == -1:
                a = self.shx.atoms.all_atoms[sdmitem.a1]
                if a.ishydrogen or a.qpeak:
                    continue
            if sdmitem.symmetry_number == 0 and sdmitem.floorD == Array([0, 0, 0]):
                flags[sdmitem.a1] = 1
                continue
            if sdmitem.dist > 2.4:
                continue
            self.shx.atoms.all_atoms[i].frac_coords = self.shx.symmcards[sdmitem.symmetry_number].matrix * \
                                               Array(sdmitem.atom1.frac_coords) \
                                               + Array(self.shx.symmcards[sdmitem.symmetry_number].trans) \
                                               - Array(sdmitem.floorD)
            flags[sdmitem.a1] = 1
        t8 = time.perf_counter()
        print('Zeit4 Bring atoms together:', t8 - t7)
        return need_symm

    def vector_length(self, x: float, y: float, z: float) -> float:
        """
        Calculates the vector length given in fractional coordinates.
        """
        a = 0.0 if (self.shx.cell[5] == 90.0) else 2.0 * x * y * self.aga
        b = 0.0 if (self.shx.cell[4] == 90.0) else 2.0 * x * z * self.bbe
        c = 0.0 if (self.shx.cell[3] == 90.0) else 2.0 * y * z * self.cal
        return sqrt(x ** 2 * self.shx.cell[0] ** 2 + y ** 2 * self.shx.cell[1] ** 2 
                    + z ** 2 * self.shx.cell[2] ** 2 + a + b + c)

    def packer(self, sdm: 'SDM', need_symm: list):
        """
        Packs atoms of the asymmetric unit to real molecules.
        :param need_symm:
        :return:
        """
        t1 = time.perf_counter()
        asymm = self.shx.atoms.all_atoms
        showatoms = asymm[:]
        for j in need_symm:
            s, h, k, l = j
            h -= 5
            k -= 5
            l -= 5
            s = s - 1
            for atom in asymm:
                if not atom.ishydrogen:
                    new_atom = Atom(self.shx)
                    new_atom.set_atom_parameters(
                        name=atom.name + ">>" + 'new',
                        sfac_num=atom.sfac_num,
                        coords=Array(atom.frac_coords) * self.shx.symmcards[s].matrix \
                                           + Array(self.shx.symmcards[s].trans) + Array([h, k, l]),
                        part=atom.part,
                        afix=AFIX(self.shx, ('AFIX ' + atom.afix).split()) if atom.afix else None,
                        resi=RESI(self.shx, ('RESI ' + atom.resinum + atom.resiclass).split()) if atom.resi else None,
                        site_occupation=atom.sof,
                        uvals=atom.uvals
                        )
                    # TODO: I have to transform the Uijs by symmetry here later.
                    if new_atom.part.n >= 0:
                        for atom in showatoms:
                            if atom.ishydrogen:
                                continue
                    if atom.part.n != new_atom.part.n:
                        continue
                    if sdm.vector_length(new_atom.frac_coords[0] - atom.frac_coords[0],
                                         new_atom.frac_coords[1] - atom.frac_coords[1],
                                         new_atom.frac_coords[2] - atom.frac_coords[2]) > 0.2:
                        showatoms.append(new_atom)
                #elif grow_qpeaks:
                #    add q-peaks
        t2 = time.perf_counter()
        print('Packzeit:', t2-t1)
        return showatoms
    
if __name__ == "__main__":
    shx = ShelXFile('tests/p-31c.res')
    sdm = SDM(shx)
    needsymm = sdm.calc_sdm()
    packed_atoms = sdm.packer(sdm, needsymm)
    print(len(shx.atoms))
    print(len(packed_atoms))
    from shelxfile.misc import wrap_line
    for y in packed_atoms:
        pass
        print(wrap_line(str(y)))