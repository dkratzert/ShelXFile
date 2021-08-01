import unittest
from collections import namedtuple
from pathlib import Path
from unittest import TestCase

from shelxfile.misc.misc import wrap_line
from shelxfile.shelx.sdm_rust import SDMR
from shelxfile.shelx.sdm_complete import get_vector_length
from shelxfile.shelx.shelx import Shelxfile
from shelxfile.shelx.sdm import SDM


class MySDMtest(TestCase):
    def setUp(self) -> None:
        self.maxDiff=None
        self.head = """
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

        self.tail = """\n
HKLF 4
 
REM  p-31c-neu in P-31c
REM R1 =  0.0308 for    4999 Fo > 4sig(Fo)  and  0.0343 for all    5352 data
REM    287 parameters refined using    365 restraints
 
END 
 
WGHT      0.0348      0.6278 
"""

    def test_SDM(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p-31c.res')
        sdm = SDM(shx)
        needsymm = sdm.calc_sdm()
        packed_atoms = sdm.packer(sdm, needsymm)
        # print(needsymm)
        # [[8, 5, 5, 5, 1], [16, 5, 5, 5, 1], [7, 4, 5, 5, 3]]
        # print(len(shx.atoms))
        # print(len(packed_atoms))

        for at in packed_atoms:
            if at.qpeak:
                continue
            # print(wrap_line(str(at)))
            self.head += wrap_line(str(at)) + '\n'
        self.head += self.tail
        p = Path('tests/resources/test-sdm1.res')
        print('Zeit für sdm:', round(sdm.sdmtime, 3), 's')
        self.assertEqual(p.read_text(), self.head)
        # print(sdm.bondlist)
        # print(len(sdm.bondlist), '(170) Atome in p-31c.res, (208) in I-43d.res')

    #@unittest.skip("Rust version does not work atm.")
    def test_SDM_rustversion(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p-31c.res')
        sdm = SDMR(shx)
        needsymm = sdm.calc_sdm()
        packed_atoms = sdm.packer(sdm, needsymm)

        for at in packed_atoms:
            if at.qpeak:
                print(at)
                continue
            self.head += wrap_line(str(at)) + '\n'
        self.head += self.tail
        p = Path('tests/resources/test-sdm1.res')
        print('Zeit für sdm:', round(sdm.sdmtime, 3), 's')
        self.assertEqual(p.read_text(), self.head)
        # print(sdm.bondlist)
        # print(len(sdm.bondlist), '(170) Atome in p-31c.res, (208) in I-43d.res')

    def test_vector_length(self):
        shx = Shelxfile()
        shx.read_string(self.head + self.tail)
        sdm = SDM(shx)
        self.assertEqual(7.343581289102655, sdm.vector_length(-0.3665069999999999, 0.293439, -0.06597900000000001))

    def test_vector_length2(self):
        Cell = namedtuple('UnitCell', 'a, b, c, alpha, beta, gamma')
        cell = Cell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
        print(get_vector_length(cell, 1.0, 1.0, 1.0))
        self.assertEqual(1.7320508075688774, get_vector_length(cell, 1.0, 1.0, 1.0))

    def test_vector_length3(self):
        Cell = namedtuple('UnitCell', 'a, b, c, alpha, beta, gamma')
        cell = Cell(12.5067, 12.5067, 24.5615, 90.0, 90.0, 120.0)
        print(get_vector_length(cell, 0.32436, 0.31780200000000003, 0.32296800000000003))
        self.assertEqual(8.891366756044848, get_vector_length(cell, 0.32436, 0.31780200000000003, 0.32296800000000003))