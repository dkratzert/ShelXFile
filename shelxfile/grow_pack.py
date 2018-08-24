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

    def __lt__(self, a2):
        d1 = a1.d  # a1.a2 * 99999 + a1.d;
        d2 = a2.d  # a2.a2 * 99999 + a2.d;
        return True if d1 < d2 else False


class SDM():
    def __init__(self):
        self.dk, dddd = 0.0, 0.0
        self.prime = np.array(3, 0)
        self.dp = np.array(3, 0)
        self.D = np.array(3, 0)
        self.floorD = np.array(3, 0)
        # contact.clear()  # sdmitem list for hydrogen contacts (not needed)?
        self.sdm = []  # list of sdmitems

    def calc_sdm(self):
        pass