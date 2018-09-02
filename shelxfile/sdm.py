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
from math import floor
from shelxfile.dsrmath import vector_length, fmin

import numpy as np

class SDMItem(object):
    def __init__(self):
        self.dist = 0.0
        self.atom1 = 0
        self.atom2 = 0
        self.symmetry_number = 0
        self.floor_dist = None
        self.covalent = True

    @property
    def dist(self):
        return self.dist

    @dist.setter
    def dist(self, d):
        self.dist = d

    def __lt__(self, a1, a2):
        d1 = a1.d  # a1.a2 * 99999 + a1.d;
        d2 = a2.d  # a2.a2 * 99999 + a2.d;
        return True if d1 < d2 else False


class SDM():
    def __init__(self, shx):
        self.shx = shx
        self.dk, dddd = 0.0, 0.0
        self.prime = np.array(3, 0)
        self.dp = np.array(3, 0)
        self.D = np.array(3, 0)
        self.floorD = np.array(3, 0)
        # contact.clear()  # sdmitem list for hydrogen contacts (not needed)?
        self.sdm = []  # list of sdmitems
        self.knots = []

    def calc_sdm(self):
        for i, at1 in enumerate(self.shx.atoms.all_atoms):
            kn = []  # list of atom neigbors
            for j, at2 in enumerate(self.shx.atoms):
                min = 1000000
                hma = False
                sdmItem = SDMItem()
                for n, symop in enumerate(self.shx.symmcards):
                    prime = np.matrix(self.shx.symmcards[n].matrix) * np.array(at1.frac_coords) + np.array(self.shx.symmcards[n].trans)
                    D = prime - at2.frac_coords + np.array([0.5, 0.5, 0.5])
                    floorD = np.array([np.floor(D.x), floor(D.y), floor(D.z)])
                    dp = D - floorD - np.array((0.5, 0.5, 0.5))
                    dk = vector_length(dp.x, dp.y, dp.z, self.shx.cell)
                    if n:
                        dk += 0.0001
                        if ((dk > 0.01) and (min >= dk)):
                            min = fmin(dk, min)
                            sdmItem.d = min
                            sdmItem.floorD = floorD
                            sdmItem.a1 = at1
                            sdmItem.a2 = at2
                            sdmItem.sn = n
                            hma = True
                if (sdmItem.a1.part * sdmItem.a2.part == 0) or (sdmItem.a1.part == sdmItem.a2.part):
                    dddd = at1.radius + at2.radius * 0.012
                else:
                    dddd = 0.0
                if sdmItem.d < dddd:
                    if hma:
                        self.knots.neighbors.append(j)
                        sdmItem.covalent = True
                else:
                    sdmItem.covalent = False
                if hma:
                   self.sdm.append(sdmItem)
        self.knots.append(kn)