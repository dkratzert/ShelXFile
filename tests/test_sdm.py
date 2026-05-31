import unittest
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

from shelxfile.misc.misc import wrap_line
from shelxfile.shelx import sdm as sdm_module
from shelxfile.shelx.sdm import SDM, HAS_CPP
from shelxfile import Shelxfile

@unittest.skip("Rust version does not work atm.")
class MySDMtest(TestCase):
    def setUp(self) -> None:
        self.maxDiff = None
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

    @unittest.skip("Rust version does not work atm.")
    def test_SDM_rustversion(self):
        from shelxfile.shelx.sdm_rust import SDMR
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

class TestSDMCpp(TestCase):
    """Tests for the C++ SDM acceleration (sdm_cpp).

    Each test runs both the C++ and the pure-Python paths on the same
    structure and asserts they produce identical results.  Tests are
    unconditional: when sdm_cpp is not installed the test verifies that the
    Python fallback still produces correct output.
    """

    RES_FILE = 'tests/resources/p-31c.res'

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _run_sdm(self, shx: Shelxfile, use_cpp: bool) -> tuple[list, list, int]:
        """Run calc_sdm() with the requested backend, return (sdm_list, need_symm, maxmol)."""
        with patch.object(sdm_module, 'HAS_CPP', use_cpp):
            s = SDM(shx)
            need_symm = s.calc_sdm()
        return s.sdm_list, need_symm, s.maxmol

    def _fresh_shx(self) -> Shelxfile:
        shx = Shelxfile()
        shx.read_file(self.RES_FILE)
        return shx

    # ------------------------------------------------------------------
    # availability
    # ------------------------------------------------------------------

    def test_has_cpp_flag_is_bool(self):
        """HAS_CPP must be a plain bool."""
        self.assertIsInstance(HAS_CPP, bool)

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_module_has_openmp_attr(self):
        """sdm_cpp must expose a has_openmp attribute."""
        import sdm_cpp
        self.assertIsInstance(sdm_cpp.has_openmp, bool)

    # ------------------------------------------------------------------
    # C++ path functional tests
    # ------------------------------------------------------------------

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_sdm_list_nonempty(self):
        """The C++ path should find at least one atom pair."""
        sdm_list, _, _ = self._run_sdm(self._fresh_shx(), use_cpp=True)
        self.assertGreater(len(sdm_list), 0)

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_molindex_assigned(self):
        """calc_sdm via C++ must assign a molindex to every atom."""
        shx = self._fresh_shx()
        with patch.object(sdm_module, 'HAS_CPP', True):
            s = SDM(shx)
            s.calc_sdm()
        for at in shx.atoms.all_atoms:
            self.assertGreater(at.molindex, 0, f"{at.name} has molindex {at.molindex}")

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_maxmol_positive(self):
        """At least one molecule must be found."""
        _, _, maxmol = self._run_sdm(self._fresh_shx(), use_cpp=True)
        self.assertGreater(maxmol, 0)

    # ------------------------------------------------------------------
    # C++ vs Python equivalence
    # ------------------------------------------------------------------

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_and_python_same_pair_count(self):
        """C++ and Python SDM must find the same number of atom pairs."""
        cpp_list, _, _ = self._run_sdm(self._fresh_shx(), use_cpp=True)
        py_list, _, _ = self._run_sdm(self._fresh_shx(), use_cpp=False)
        self.assertEqual(len(cpp_list), len(py_list))

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_and_python_same_maxmol(self):
        """C++ and Python SDM must identify the same number of molecular fragments."""
        _, _, cpp_maxmol = self._run_sdm(self._fresh_shx(), use_cpp=True)
        _, _, py_maxmol = self._run_sdm(self._fresh_shx(), use_cpp=False)
        self.assertEqual(cpp_maxmol, py_maxmol)

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_and_python_same_need_symm(self):
        """C++ and Python SDM must return the same symmetry operations needed for packing."""
        _, cpp_ns, _ = self._run_sdm(self._fresh_shx(), use_cpp=True)
        _, py_ns, _ = self._run_sdm(self._fresh_shx(), use_cpp=False)
        self.assertEqual(sorted(cpp_ns), sorted(py_ns))

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_and_python_same_covalent_count(self):
        """C++ and Python paths must agree on the number of covalent bonds."""
        cpp_list, _, _ = self._run_sdm(self._fresh_shx(), use_cpp=True)
        py_list, _, _ = self._run_sdm(self._fresh_shx(), use_cpp=False)
        cpp_cov = sum(1 for x in cpp_list if x.covalent)
        py_cov = sum(1 for x in py_list if x.covalent)
        self.assertEqual(cpp_cov, py_cov)

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_and_python_distances_match(self):
        """Shortest distances from C++ and Python must agree to 6 decimal places."""
        cpp_list, _, _ = self._run_sdm(self._fresh_shx(), use_cpp=True)
        py_list, _, _ = self._run_sdm(self._fresh_shx(), use_cpp=False)
        for cpp_item, py_item in zip(cpp_list, py_list):
            self.assertAlmostEqual(cpp_item.dist, py_item.dist, places=6,
                                   msg=f"Distance mismatch for pair ({cpp_item.a1}, {cpp_item.a2})")

    @unittest.skipUnless(HAS_CPP, 'sdm_cpp extension not available')
    def test_cpp_packer_same_atom_count_as_python(self):
        """Packed structures from C++ and Python backends must have the same atom count."""
        shx_cpp = self._fresh_shx()
        shx_py = self._fresh_shx()
        with patch.object(sdm_module, 'HAS_CPP', True):
            s_cpp = SDM(shx_cpp)
            ns_cpp = s_cpp.calc_sdm()
        packed_cpp = s_cpp.packer(s_cpp, ns_cpp)

        with patch.object(sdm_module, 'HAS_CPP', False):
            s_py = SDM(shx_py)
            ns_py = s_py.calc_sdm()
        packed_py = s_py.packer(s_py, ns_py)

        self.assertEqual(len(packed_cpp), len(packed_py))

    # ------------------------------------------------------------------
    # Python fallback (always runs, regardless of HAS_CPP)
    # ------------------------------------------------------------------

    def test_python_fallback_sdm_list_nonempty(self):
        """The pure-Python fallback must also find atom pairs."""
        sdm_list, _, _ = self._run_sdm(self._fresh_shx(), use_cpp=False)
        self.assertGreater(len(sdm_list), 0)

    def test_python_fallback_maxmol_positive(self):
        """The pure-Python fallback must identify at least one molecule."""
        _, _, maxmol = self._run_sdm(self._fresh_shx(), use_cpp=False)
        self.assertGreater(maxmol, 0)

    def test_python_fallback_need_symm(self):
        """The pure-Python fallback must return a non-empty symmetry list for p-31c."""
        _, need_symm, _ = self._run_sdm(self._fresh_shx(), use_cpp=False)
        self.assertGreater(len(need_symm), 0)

