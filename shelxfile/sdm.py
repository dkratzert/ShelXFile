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
from math import floor, sqrt
from shelxfile.dsrmath import vector_length, fmin, Matrix, Array
from shelxfile.shelx import ShelXFile


class SDMItem(object):
    def __init__(self):
        self.dist = 0.0
        self.atom1 = None
        self.atom2 = None
        self.symmetry_number = 0
        self.floor_dist = None
        self.covalent = True

    def __lt__(self, a2):
        #d1 = self.dist  # a1.a2 * 99999 + a1.d;
        #d2 = a2.dist  # a2.a2 * 99999 + a2.d;
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
        brauchSymm = []
        t1 = time.perf_counter()
        for at1 in self.shx.atoms.all_atoms:
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
                        sdmItem.symmetry_number = n
                        hma = True
                if (not sdmItem.atom1.ishydrogen and not sdmItem.atom2.ishydrogen) and \
                        sdmItem.atom1.part.n * sdmItem.atom2.part.n == 0 or sdmItem.atom1.part.n == sdmItem.atom2.part.n:
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
                    if sdmItem.atom1.an == sdmItem.atom2.an and sdmItem.atom1.an == 0:
                        continue
                    prime = Array(sdmItem.atom1.frac_coords) * symop.matrix + symop.trans
                    D = prime - Array(sdmItem.atom2.frac_coords) + Array([0.5, 0.5, 0.5])
                    floorD = D.floor
                    dp = D - floorD - Array([0.5, 0.5, 0.5])
                    if n == 0 and Array([0, 0, 0]) == floorD:
                        #print(floorD)
                        continue
                    dk = self.vector_length(*dp)
                    dddd = (sdmItem.dist + 0.2)
                    bs = ''
                    if (dk > 0.001) and (dddd >= dk):
                        #bs = "{}_{}{}{}".format(n+1, (5-int(floorD[0])), (5-int(floorD[1])), (5-int(floorD[2])) )
                        bs = [n+1, (5-floorD[0]), (5-int(floorD[1])), (5-int(floorD[2]))]
                        if not bs in brauchSymm:
                            brauchSymm.append(bs)
                        #print(symop)
                        #print(bs)
        for x in brauchSymm:
            print(x)
        t4 = time.perf_counter()
        print('Zeit2 brauchsymm:', t4 - t3)
        t5 = time.perf_counter()
        self.sdm.sort()
        t6 = time.perf_counter()
        print('Zeit3 sort sdm:', t6 - t5)
        return brauchSymm
        

    def vector_length(self, x: float, y: float, z: float) -> float:
        """
        Calculates the vector length given in fractional coordinates.

        >>> vector_length(1, 0, 0, [1, 1, 1, 90, 90, 90])
        1.0
        >>> round(vector_length(1.0, 1.0, 1.0, [5.773501, 5.773501, 5.773501, 90, 90, 90]), 5)
        10.0
        >>> round(vector_length(-0.269224, 0.349464, 0.0, [12.5067, 12.5067, 24.5615, 90.0, 90.0, 120.0]), 5)
        6.71984
        """
        a = 0.0 if (self.shx.cell[5] == 90.0) else 2.0 * x * y * self.aga
        b = 0.0 if (self.shx.cell[4] == 90.0) else 2.0 * x * z * self.bbe
        c = 0.0 if (self.shx.cell[3] == 90.0) else 2.0 * y * z * self.cal
        return sqrt(x ** 2 * self.shx.cell[0] ** 2 + y ** 2 * self.shx.cell[1] ** 2 
                    + z ** 2 * self.shx.cell[2] ** 2 + a + b + c)


    
if __name__ == "__main__":
    shx = ShelXFile('tests/p-31c.res')
    sdm = SDM(shx)
    sdm.calc_sdm()
    #print(shx.atoms)