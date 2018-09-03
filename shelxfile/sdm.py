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
from math import floor, sqrt
from shelxfile.dsrmath import vector_length, fmin
from shelxfile.shelx import ShelXFile

import numpy as np

class SDMItem(object):
    def __init__(self):
        self.dist = 0.0
        self.atom1 = 0
        self.atom2 = 0
        self.symmetry_number = 0
        self.floor_dist = None
        self.covalent = True

    def __lt__(self, a1, a2):
        d1 = a1.d  # a1.a2 * 99999 + a1.d;
        d2 = a2.d  # a2.a2 * 99999 + a2.d;
        return True if d1 < d2 else False


class SDM():
    def __init__(self, shx):
        self.shx = shx
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
        for i, at1 in enumerate(self.shx.atoms.all_atoms):
            atneighb = []  # list of atom neigbors
            for j, at2 in enumerate(self.shx.atoms):
                min = 1000000
                hma = False
                sdmItem = SDMItem()
                for n, symop in enumerate(self.shx.symmcards):
                    prime = np.array([at1.frac_coords]) * np.matrix(self.shx.symmcards[n].matrix.values) + np.array([self.shx.symmcards[n].trans])
                    D = prime - np.array(at2.frac_coords) + np.array([0.5, 0.5, 0.5])
                    dp = D - np.floor(D) - np.array((0.5, 0.5, 0.5))
                    dp = dp.reshape(3,1)
                    dk = self.vector_length(dp[0], dp[1], dp[2])
                    if n:
                        dk += 0.0001
                    if (dk > 0.01) and (min >= dk):
                        min = fmin(dk, min)
                        sdmItem.dist = min
                        sdmItem.floorD = np.floor(D)
                        sdmItem.a1 = at1
                        sdmItem.a2 = at2
                        sdmItem.sn = n
                        hma = True
                if (sdmItem.a1.part.n * sdmItem.a2.part.n == 0) or (sdmItem.a1.part.n == sdmItem.a2.part.n):
                    dddd = at1.radius + at2.radius * 0.012
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

        for k, sdmItem in enumerate(self.sdm):
            if sdmItem.covalent:
                for n, at in enumerate(self.shx.atoms):
                    if (self.shx.atoms[self.sdm[k].a1].part != 0) and (self.shx.atoms[self.sdm[k].a2].part != 0) \
                            and (self.shx.atoms[self.sdm[k].a1].part != self.shx.atoms[self.sdm[k].a2].part):
                        continue
                    if (self.shx.atoms[self.sdm[k].a1].an == self.shx.atoms[self.sdm[k].a2].an) \
                            and (self.shx.atoms[self.sdm[k].a1].an == 0):
                        continue
                    prime = self.shx.symmcards[n].matrix * self.shx.atoms[self.sdm[k].a1].frac_coords + self.shx.symmcards[n].trans
                    D = prime - self.shx.atoms[self.sdm[k].a2].frac_coords + np.array([0.5, 0.5, 0.5])
                    floorD = np.array([floor(D.x), floor(D.y), floor(D.z)])
                    dp = D - floorD - np.array([0.5, 0.5, 0.5])
                    if (n == 0) and (np.array([0., 0., 0.]) == floorD):
                        continue
                    dk = self.vector_length(dp.x, dp.y, dp.z)
                    dddd = (self.sdm[k].d + 0.2)
                    bs = ''
                    if (dk > 0.001) and (dddd >= dk):
                        bs = "%1_%2%3%4".format(n+1, (5-int(floorD[0])), (5-int(floorD[1])), (5-int(floorD[2])) )
                    if not bs in brauchSymm:
                        brauchSymm.append(bs)
                    print(bs)

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
        cell = self.shx.cell
        a = 0.0 if (cell[5] == 90.0) else 2.0 * x * y * cell[0] * cell[1] * cell.cosga
        b = 0.0 if (cell[4] == 90.0) else 2.0 * x * z * cell[0] * cell[2] * cell.cosbe
        c = 0.0 if (cell[3] == 90.0) else 2.0 * y * z * cell[1] * cell[2] * cell.cosal
        return sqrt(x ** 2 * cell[0] ** 2 + y ** 2 * cell[1] ** 2 + z ** 2 * cell[2] ** 2 + a + b + c)


    
if __name__ == "__main__":
    shx = ShelXFile('tests/p-31c.res')
    sdm = SDM(shx)
    sdm.calc_sdm()
    #print(shx.atoms)