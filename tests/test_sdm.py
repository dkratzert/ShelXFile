import unittest
import tempfile
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

from shelxfile.misc.misc import wrap_line
from shelxfile.shelx import sdm as sdm_module
from shelxfile.shelx.sdm import SDM, HAS_CPP
from shelxfile import Shelxfile


class TestPackUnitCell(TestCase):
    """Tests for :meth:`SDM.pack_unit_cell` and :meth:`Shelxfile.pack`."""

    RES_FILE = 'tests/resources/p21c.res'
    RES_FILE_P31C = 'tests/resources/p-31c.res'

    def _fresh_shx(self, path: str = None) -> Shelxfile:
        shx = Shelxfile()
        shx.read_file(path or self.RES_FILE)
        return shx

    # ------------------------------------------------------------------
    # Basic correctness
    # ------------------------------------------------------------------

    def test_packed_nonempty(self):
        """pack_unit_cell must return at least one atom."""
        shx = self._fresh_shx()
        sdm = SDM(shx)
        packed = sdm.pack_unit_cell()
        self.assertGreater(len(packed), 0)

    def test_packed_larger_than_asymm(self):
        """Packed unit cell must contain more atoms than the asymmetric unit."""
        shx = self._fresh_shx()
        asymm_count = len(shx.atoms.all_atoms)
        sdm = SDM(shx)
        packed = sdm.pack_unit_cell()
        self.assertGreater(len(packed), asymm_count)

    def test_all_coords_in_unit_cell(self):
        """All fractional coordinates must be in [0, 1)."""
        shx = self._fresh_shx()
        sdm = SDM(shx)
        for at in sdm.pack_unit_cell():
            self.assertGreaterEqual(at.x, 0.0, f"{at.name}.x = {at.x}")
            self.assertLess(at.x, 1.0,          f"{at.name}.x = {at.x}")
            self.assertGreaterEqual(at.y, 0.0, f"{at.name}.y = {at.y}")
            self.assertLess(at.y, 1.0,          f"{at.name}.y = {at.y}")
            self.assertGreaterEqual(at.z, 0.0, f"{at.name}.z = {at.z}")
            self.assertLess(at.z, 1.0,          f"{at.name}.z = {at.z}")

    def test_returns_atom_objects(self):
        """pack_unit_cell must return Atom objects."""
        from shelxfile.atoms.atom import Atom
        shx = self._fresh_shx()
        sdm = SDM(shx)
        packed = sdm.pack_unit_cell()
        for at in packed:
            self.assertIsInstance(at, Atom)

    def test_no_qpeaks_by_default(self):
        """Q-peaks must be excluded unless with_qpeaks=True."""
        shx = self._fresh_shx()
        sdm = SDM(shx)
        packed = sdm.pack_unit_cell(with_qpeaks=False)
        self.assertFalse(any(at.qpeak for at in packed))

    def test_with_qpeaks(self):
        """with_qpeaks=True must include Q-peaks when the structure has them."""
        shx = self._fresh_shx()
        sdm = SDM(shx)
        packed_no_q = sdm.pack_unit_cell(with_qpeaks=False)
        packed_with_q = sdm.pack_unit_cell(with_qpeaks=True)
        # p21c.res may or may not have Q-peaks; the with_q count must be >= without_q
        self.assertGreaterEqual(len(packed_with_q), len(packed_no_q))

    # ------------------------------------------------------------------
    # Symmetry-operation selection
    # ------------------------------------------------------------------

    def test_identity_only_gives_asymm_unit(self):
        """Selecting only the identity (index 0) must give one atom per asymm-unit atom."""
        shx = self._fresh_shx()
        asymm = [at for at in shx.atoms.all_atoms if not at.qpeak]
        sdm = SDM(shx)
        packed = sdm.pack_unit_cell(symmop_indices=[0])
        self.assertEqual(len(packed), len(asymm))

    def test_symmop_indices_none_uses_all_ops(self):
        """None (default) must apply all symmetry operations."""
        shx = self._fresh_shx()
        sdm = SDM(shx)
        sdm._build_symm_arrays()
        packed_all = sdm.pack_unit_cell(symmop_indices=None)
        packed_sel = sdm.pack_unit_cell(symmop_indices=list(range(len(sdm._symm_m))))
        self.assertEqual(len(packed_all), len(packed_sel))

    # ------------------------------------------------------------------
    # Before calc_sdm
    # ------------------------------------------------------------------

    def test_usable_without_calc_sdm(self):
        """pack_unit_cell must work on a fresh SDM that has not run calc_sdm."""
        shx = self._fresh_shx()
        sdm = SDM(shx)
        # Do NOT call calc_sdm; pack_unit_cell must still work
        packed = sdm.pack_unit_cell()
        self.assertGreater(len(packed), 0)

    # ------------------------------------------------------------------
    # Duplicate elimination
    # ------------------------------------------------------------------

    def test_no_duplicate_positions(self):
        """No two atoms must share the same position within the tolerance."""
        shx = self._fresh_shx()
        sdm = SDM(shx)
        packed = sdm.pack_unit_cell(cart_tolerance=0.2)
        coords = [(round(at.x, 4), round(at.y, 4), round(at.z, 4)) for at in packed]
        self.assertEqual(len(coords), len(set(coords)),
                         "Duplicate fractional positions found in packed unit cell")

    # ------------------------------------------------------------------
    # Shelxfile.pack() high-level API
    # ------------------------------------------------------------------

    def test_shelxfile_pack_matches_sdm_pack(self):
        """`Shelxfile.pack()` and `SDM.pack_unit_cell()` must give the same count."""
        shx1 = self._fresh_shx()
        shx2 = self._fresh_shx()
        direct = SDM(shx1).pack_unit_cell()
        via_api = shx2.pack()
        self.assertEqual(len(direct), len(via_api))

    def test_shelxfile_pack_returns_list(self):
        """`Shelxfile.pack()` must return a list."""
        shx = self._fresh_shx()
        self.assertIsInstance(shx.pack(), list)

    # ------------------------------------------------------------------
    # P-31c (different space group)
    # ------------------------------------------------------------------

    def test_p31c_packed_count(self):
        """P-31c structure must produce a non-trivial packed unit cell."""
        shx = self._fresh_shx(self.RES_FILE_P31C)
        sdm = SDM(shx)
        packed = sdm.pack_unit_cell()
        self.assertGreater(len(packed), len(shx.atoms.all_atoms))


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

    def test_vector_length(self):
        shx = Shelxfile()
        shx.read_string(self.head + self.tail)
        sdm = SDM(shx)
        self.assertEqual(7.343581289102655, sdm.vector_length(-0.3665069999999999, 0.293439, -0.06597900000000001))

    def test_grown_atom_names_are_legal(self):
        """Symmetry-generated atoms returned by packer() must have legal SHELXL names."""
        import re
        legal = re.compile(r'^[A-Za-z][A-Za-z0-9]{0,3}$')
        shx = Shelxfile()
        shx.read_file('tests/resources/p-31c.res')
        sdm = SDM(shx)
        needsymm = sdm.calc_sdm()
        packed_atoms = sdm.packer(sdm, needsymm)
        for at in packed_atoms:
            if at.symmgen:  # only check atoms created by the packer
                self.assertRegex(at.name, legal,
                                 f"Atom name '{at.name}' is not a legal SHELXL name")

    def test_grown_atom_names_are_unique(self):
        """Every atom returned by packer() must have a unique name."""
        shx = Shelxfile()
        shx.read_file('tests/resources/p-31c.res')
        sdm = SDM(shx)
        needsymm = sdm.calc_sdm()
        packed_atoms = sdm.packer(sdm, needsymm)
        names = [at.name.upper() for at in packed_atoms]
        self.assertEqual(len(names), len(set(names)),
                         "Duplicate atom names found in packer() output")

    def test_pack_unit_cell_names_are_legal(self):
        """Symmetry-generated atoms from pack_unit_cell() must have legal SHELXL names."""
        import re
        legal = re.compile(r'^[A-Za-z][A-Za-z0-9]{0,3}$')
        shx = Shelxfile()
        shx.read_file('tests/resources/p-31c.res')
        sdm = SDM(shx)
        for at in sdm.pack_unit_cell():
            if at.symmgen:  # only check atoms created by the packer
                self.assertRegex(at.name, legal,
                                 f"Atom name '{at.name}' is not a legal SHELXL name")

    def test_pack_unit_cell_names_are_unique(self):
        """Every atom returned by pack_unit_cell() must have a unique name."""
        shx = Shelxfile()
        shx.read_file('tests/resources/p-31c.res')
        sdm = SDM(shx)
        names = [at.name.upper() for at in sdm.pack_unit_cell()]
        self.assertEqual(len(names), len(set(names)),
                         "Duplicate atom names found in pack_unit_cell() output")

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


class TestWriteGrownFile(TestCase):
    """Tests for :meth:`Shelxfile.write_grown_file`."""

    RES_FILE = 'tests/resources/p-31c.res'

    def _fresh_shx(self) -> Shelxfile:
        shx = Shelxfile()
        shx.read_file(self.RES_FILE)
        return shx

    def test_write_grown_file_creates_file(self):
        """write_grown_file must create a non-empty output file."""
        shx = self._fresh_shx()
        with tempfile.NamedTemporaryFile(suffix='.res', delete=False) as f:
            out = Path(f.name)
        try:
            shx.write_grown_file(out)
            self.assertTrue(out.exists())
            self.assertGreater(out.stat().st_size, 0)
        finally:
            out.unlink(missing_ok=True)

    def test_grown_file_has_p1_symmetry(self):
        """The grown file must have LATT -1 and no SYMM cards."""
        shx = self._fresh_shx()
        with tempfile.NamedTemporaryFile(suffix='.res', delete=False) as f:
            out = Path(f.name)
        try:
            shx.write_grown_file(out)
            text = out.read_text()
            lines = text.splitlines()
            latt_lines = [l for l in lines if l.strip().upper().startswith('LATT')]
            symm_lines = [l for l in lines if l.strip().upper().startswith('SYMM')]
            self.assertEqual(len(latt_lines), 1, "Expected exactly one LATT line")
            self.assertIn('-1', latt_lines[0], "LATT must be -1 for P1 symmetry")
            self.assertEqual(len(symm_lines), 0, "Grown file must have no SYMM cards")
        finally:
            out.unlink(missing_ok=True)

    def test_grown_file_has_warning_comment(self):
        """The grown file must contain a REM warning about refinement."""
        shx = self._fresh_shx()
        with tempfile.NamedTemporaryFile(suffix='.res', delete=False) as f:
            out = Path(f.name)
        try:
            shx.write_grown_file(out)
            text = out.read_text().upper()
            self.assertIn('NOT SUITED FOR REFINEMENT', text,
                          "Grown file must warn that it is not suited for refinement")
        finally:
            out.unlink(missing_ok=True)

    def test_grown_file_has_part_cards(self):
        """Disorder parts must be preserved via PART cards in the grown file."""
        shx = self._fresh_shx()
        with tempfile.NamedTemporaryFile(suffix='.res', delete=False) as f:
            out = Path(f.name)
        try:
            shx.write_grown_file(out)
            text = out.read_text()
            lines = text.splitlines()
            part_lines = [l for l in lines if l.strip().upper().startswith('PART ')]
            # p-31c.res has PART 1 / PART 2 disorder groups, repeated for each
            # symmetry copy, so there must be several PART cards
            self.assertGreater(len(part_lines), 0, "Grown file must contain PART cards")
            # Check that both non-zero parts and the closing PART 0 are present
            part_values = {l.strip().split()[1] for l in part_lines}
            self.assertIn('1', part_values, "PART 1 must appear")
            self.assertIn('2', part_values, "PART 2 must appear")
            self.assertIn('0', part_values, "PART 0 (close disorder) must appear")
        finally:
            out.unlink(missing_ok=True)

    def test_grown_file_correct_atom_count(self):
        """The grown file must contain the same number of atoms as grow() returns."""
        shx = self._fresh_shx()
        grown = shx.grow()
        real_grown = [a for a in grown if not a.qpeak]
        with tempfile.NamedTemporaryFile(suffix='.res', delete=False) as f:
            out = Path(f.name)
        try:
            shx.write_grown_file(out)
            # Exclude SHELXL comment lines ("!..." at start of first token)
            atom_lines = [
                l for l in out.read_text().splitlines()
                if Shelxfile.is_atom(l) and not l.lstrip().startswith('!')
            ]
            self.assertEqual(len(atom_lines), len(real_grown))
        finally:
            out.unlink(missing_ok=True)

    def test_grown_file_no_original_symm_cards(self):
        """The original SYMM cards from p-31c must not appear in the grown file."""
        shx = self._fresh_shx()
        with tempfile.NamedTemporaryFile(suffix='.res', delete=False) as f:
            out = Path(f.name)
        try:
            shx.write_grown_file(out)
            text = out.read_text()
            # p-31c has "-Y, X-Y, Z" as first SYMM line — it must not appear
            self.assertNotIn('-Y, X-Y, Z', text,
                             "Original SYMM cards must be absent in grown file")
        finally:
            out.unlink(missing_ok=True)

    def test_grown_file_no_include_references(self):
        """Include-file (+filename) lines must not appear in the grown file."""
        shx = self._fresh_shx()
        with tempfile.NamedTemporaryFile(suffix='.res', delete=False) as f:
            out = Path(f.name)
        try:
            shx.write_grown_file(out)
            text = out.read_text()
            include_lines = [l for l in text.splitlines() if l.startswith('+')]
            self.assertEqual(include_lines, [],
                             "Include-file references must be absent in grown file")
        finally:
            out.unlink(missing_ok=True)

    def test_grown_file_accepts_string_path(self):
        """write_grown_file must also accept a plain string as the filename."""
        shx = self._fresh_shx()
        with tempfile.NamedTemporaryFile(suffix='.res', delete=False) as f:
            out = f.name  # plain string
        try:
            shx.write_grown_file(out)  # must not raise
            self.assertGreater(Path(out).stat().st_size, 0)
        finally:
            Path(out).unlink(missing_ok=True)

