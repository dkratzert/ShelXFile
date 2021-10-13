from shutil import which
from unittest import TestCase

from shelxfile import Shelxfile
from shelxfile.refine.refine import get_xl_version_string


class TestRefine(TestCase):
    def test_get_xl_version_string_with_no_path(self):
        self.assertEqual('', get_xl_version_string(''))

    def test_get_xl_version_string_with_real_path(self):
        shx = which('shelxl') if which('shelxl') else which('xl')
        self.assertEqual('2018/3', get_xl_version_string(shx))

    def test_check_cycle_numbers(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/complete_run/p21c.res')
        self.assertEqual('L.S. 10', str(shx.cycles))

    def test_set_cycle_numbers_directly(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/complete_run/p21c.res')
        shx.cycles.number = 3
        self.assertEqual('L.S. 3', str(shx.cycles))

    def test_set_cycle_numbers_by_set_function(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/complete_run/p21c.res')
        shx.cycles.set_refine_cycles(4)
        self.assertEqual('L.S. 4', str(shx.cycles))
        self.assertEqual('L.S. 4', shx._reslist[shx.cycles.index].__repr__())

    def test_set_cycle_numbers_negative(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/complete_run/p21c.res')
        shx.cycles.set_refine_cycles(-1)
        self.assertEqual('L.S. -1', str(shx.cycles))

    def test_set_cycle_numbers_directly2(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/complete_run/p21c.res')
        shx.cycles.cycles = 3
        self.assertEqual('L.S. 10', str(shx.cycles))