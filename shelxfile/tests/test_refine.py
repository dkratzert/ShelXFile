import os
import shutil
from pathlib import Path
from shutil import which
from unittest import TestCase

from shelxfile import Shelxfile
from shelxfile.refine.refine import get_xl_version_string


def clean_refine_files(basename='p21c'):
    Path('{}.res'.format(basename)).unlink(missing_ok=True)
    Path('{}.ins'.format(basename)).unlink(missing_ok=True)
    Path('{}.lst'.format(basename)).unlink(missing_ok=True)
    Path('{}.fcf'.format(basename)).unlink(missing_ok=True)
    Path('{}.fcf6'.format(basename)).unlink(missing_ok=True)
    Path('{}.cif'.format(basename)).unlink(missing_ok=True)
    Path('{}.hkl'.format(basename)).unlink(missing_ok=True)


class TestRefine(TestCase):
    def setUp(self) -> None:
        os.chdir(Path(__file__).parent.parent)
        res = Path('tests/resources/complete_run/p21c.res')
        # hkl = Path('tests/resources/complete_run/p21c.hkl')
        shutil.copy(res, '.')
        # shutil.copy(hkl, '.')
        self.shx = Shelxfile()
        self.shx.read_file('./p21c.res')

    def tearDown(self) -> None:
        os.chdir(Path(__file__).parent.parent)
        clean_refine_files('p21c')

    def test_get_xl_version_string_with_no_path(self):
        self.assertEqual('', get_xl_version_string(''))

    def test_get_xl_version_string_with_real_path(self):
        shx = which('shelxl') if which('shelxl') else which('xl')
        self.assertEqual('2018/3', get_xl_version_string(shx))

    def test_check_cycle_numbers(self):
        self.assertEqual('L.S. 10', str(self.shx.cycles))

    def test_set_cycle_numbers_directly(self):
        self.shx.cycles.number = 3
        self.assertEqual('L.S. 3', str(self.shx.cycles))

    def test_set_cycle_numbers_by_set_function(self):
        self.shx.cycles.set_refine_cycles(4)
        self.assertEqual('L.S. 4', str(self.shx.cycles))
        self.assertEqual('L.S. 4', self.shx._reslist[self.shx.cycles.index].__repr__())

    def test_set_cycle_numbers_negative(self):
        self.shx.cycles.set_refine_cycles(-1)
        self.assertEqual('L.S. -1', str(self.shx.cycles))

    def test_set_cycle_numbers_directly2(self):
        self.shx.cycles.cycles = 3
        self.assertEqual('L.S. 10', str(self.shx.cycles))


class TestRefineFinishedmodel(TestCase):

    def setUp(self) -> None:
        os.chdir(Path(__file__).parent.parent)
        res = Path('tests/resources/model_finished/p21c.res')
        hkl = Path('tests/resources/model_finished/p21c.hkl')
        shutil.copy(res, '.')
        shutil.copy(hkl, '.')
        self.shx = Shelxfile()
        self.shx.read_file('./p21c.res')

    def tearDown(self) -> None:
        os.chdir(Path(__file__).parent.parent)
        clean_refine_files('p21c')

    def test_refine_with_cycle_number_set_to_4_in_LScycles(self):
        self.assertTrue('L.S. 10' in Path('p21c.res').read_text())
        self.shx.cycles.set_refine_cycles(4)
        self.shx.refine()
        txt = Path('p21c.res').read_text()
        self.assertTrue('L.S. 4' in txt)

    def test_refine_with_cycle_number_set_to_4_in_refine(self):
        self.assertTrue('L.S. 10' in Path('p21c.res').read_text())
        self.shx.refine(4)
        txt = Path('p21c.res').read_text()
        self.assertTrue('L.S. 4' in txt)

    def test_refine_with_ANIS_inserted_and_cycle_number_set_to_4_in_refine(self):
        self.assertFalse('ANIS' in Path('p21c.res').read_text())
        self.shx.insert_anis()
        self.shx.refine(3)
        self.assertTrue('ANIS' in Path('p21c.ins').read_text())
        self.assertFalse('ANIS' in Path('p21c.res').read_text())

    def test_refine_with_cycle_number_set_to_0_in_refine(self):
        self.assertTrue('L.S. 10' in Path('p21c.res').read_text())
        self.refine = self.shx.refine(0)
        txt = Path('p21c.res').read_text()
        self.assertTrue('L.S. 0' in txt)


if __name__ == '__main__':
    print('foo')
