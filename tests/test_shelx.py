import tempfile
from pathlib import Path
from unittest import TestCase

from shelxfile import Shelxfile
from shelxfile.atoms.atom import Atom
from shelxfile.shelx.cards import HKLF


class TestShelxFileIsAtom(TestCase):

    def test_is_atom_long_line_true(self):
        """
        """
        is_an_atom = Shelxfile.is_atom(atomline='O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        self.assertEqual(True, is_an_atom)

    def test_is_no_atoms_because_type_is_missing(self):
        is_an_atom = Shelxfile.is_atom(atomline='O1    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        self.assertEqual(False, is_an_atom)

    def test_is_no_atom_because_coordinate_is_missing(self):
        is_an_atom = Shelxfile.is_atom(atomline='O1  4  0.120080    0.494426  11.00000   0.01445 ...')
        self.assertEqual(False, is_an_atom)

    def test_is_no_atom_because_its_a_command(self):
        is_atom = Shelxfile.is_atom("AFIX 123")
        self.assertEqual(False, is_atom)

    def test_is_no_atom_because_its_a_command2(self):
        is_atom = Shelxfile.is_atom("AFIX")
        self.assertEqual(False, is_atom)

    def test_is_an_atom_short(self):
        is_atom = Shelxfile.is_atom('O1    3    0.120080   0.336659   0.494426')
        self.assertEqual(True, is_atom)


class TestShelxfileElementToSfac(TestCase):
    """
    SFAC C H O F Al Ga
    """

    def setUp(self) -> None:
        self.shx = Shelxfile(debug=True)
        self.shx.read_file('tests/resources/p21c.res')

    def test_elem2sfac_oxygen(self):
        self.assertEqual(3, self.shx.elem2sfac('O'))

    def test_elem2sfac_carbon(self):
        self.assertEqual(1, self.shx.elem2sfac('c'))

    def test_elem2sfac_Argon(self):
        self.assertEqual(0, self.shx.elem2sfac('Ar'))


class TestShelxfile(TestCase):

    def setUp(self) -> None:
        self.shx = Shelxfile(debug=True)
        self.shx.read_file('tests/resources/p21c.res')

    def test_sfac2elem_C(self):
        self.assertEqual('C', self.shx.sfac2elem(1))

    def test_sfac2elem_H(self):
        self.assertEqual('H', self.shx.sfac2elem(2))

    def test_sfac2elem_O(self):
        self.assertEqual('O', self.shx.sfac2elem(3))

    def test_sfac2elem_Al(self):
        self.assertEqual('Al', self.shx.sfac2elem(5))

    def test_sfac2elem_not_existent(self):
        self.assertEqual('', self.shx.sfac2elem(8))

    def test_sfac2elem_zero(self):
        self.assertEqual('', self.shx.sfac2elem(0))

    def test_sum_formula(self):
        self.assertEqual('C0.25 H0.5 O0.75 F1 AL1.25 GA1.5', self.shx.sum_formula)

    def test_read_file_to_list(self):
        shx = Shelxfile()
        shx._reslist = ['Foo']
        shx.read_file('tests/resources/p21c.res')
        # Make sure read_file() re-initializes the constructor:
        self.assertEqual(['TITL p21c in P2(1)/c',
                          '    created by SHELXL-2018/3 at 16:18:25 on 03-May-2018',
                          'CELL 0.71073 10.5086 20.9035 20.5072 90 94.13 90',
                          'ZERR 4 0.0003 0.0005 0.0005 0 0.001 0',
                          'LATT 1'],
                         [str(x) for x in shx._reslist[:5]])

    def test_read_string(self):
        self.shx.atoms = 'Foo'
        # Make sure read_string re-initializes the constructor:
        self.shx.read_string(Path('tests/resources/p21c.res').read_text(encoding='latin1'))
        assert self.shx.atoms.number == 148

    def test_read_file(self):
        shx = Shelxfile()
        with self.assertRaises(FileNotFoundError):
            shx.read_file('tests/resources/foobar.res')


class TestShelxfileGoodModel(TestCase):

    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/model_finished/p21c.res')

    def test_atom_representation(self):
        a = self.shx._reslist[39]
        self.assertEqual(
            'C1    1   -0.187504    0.282782    0.527531    11.00000    0.01807    0.02353      0.01797    0.00008   -0.00179    0.00036',
            str(a))

    def test_atom_object(self):
        a = self.shx._reslist[39]
        self.assertEqual(
            'C1    1   -0.187504    0.282782    0.527531    11.00000    0.01807    0.02353      0.01797    0.00008   -0.00179    0.00036',
            a)

    def test_qpeak(self):
        a = self.shx._reslist[-20]
        self.assertEqual('Q1   1   0.0784    0.1310    0.6428   11.00000  0.04      0.10     ', str(a))

    def test_r1(self):
        assert self.shx.R1 == 0.194

    def test_wr2(self):
        assert self.shx.wr2 == 0.51

    def test_data(self):
        assert self.shx.data == 10786

    def test_param(self):
        assert self.shx.parameters == 551

    def test_data_to_param(self):
        assert round(self.shx.dat_to_param, 4) == 19.5753

    def test_num_restraints(self):
        assert self.shx.num_restraints == 0

    def test_peak_hole(self):
        assert self.shx.highest_peak == 6.395
        assert self.shx.deepest_hole == -5.556

    def test_peak_goof(self):
        assert self.shx.goof == 5.982
        assert self.shx.rgoof == 5.982


class TestWriteFile(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/2240189.res')
        self.written = Path('tests/resources/2240189-written.res')
        self.testwritepath = Path('tests/resources/testwrite.res')

    def tearDown(self) -> None:
        self.testwritepath.unlink(missing_ok=True)

    def test_write_shelx_file(self):
        self.shx.write_shelx_file(str(self.testwritepath))
        self.assertEqual(self.written.read_bytes(), self.testwritepath.read_bytes())


# ---------------------------------------------------------------------------
# Tests for the improved add_atom / unused_atom_name interface
# ---------------------------------------------------------------------------

class TestAddAtom(TestCase):
    """Tests for Shelxfile.add_atom() and Shelxfile.unused_atom_name()."""

    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    # ------------------------------------------------------------------
    # Basic insertion
    # ------------------------------------------------------------------

    def test_add_atom_returns_atom_object(self):
        a = self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        self.assertIsInstance(a, Atom)

    def test_add_atom_appears_in_atoms_list(self):
        count_before = len(self.shx.atoms)
        self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        self.assertEqual(count_before + 1, len(self.shx.atoms))

    def test_add_atom_lookup_by_name(self):
        self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        found = self.shx.atoms.get_atom_by_name('C99')
        self.assertIsNotNone(found)
        self.assertEqual('C99', found.name)

    # ------------------------------------------------------------------
    # _reslist position: must NOT overwrite TITL (index 0)
    # ------------------------------------------------------------------

    def test_add_atom_does_not_corrupt_reslist_index_0(self):
        self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        first = str(self.shx._reslist[0])
        self.assertIn('TITL', first.upper())

    def test_add_atom_placed_before_hklf(self):
        self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        atom_idx = self.shx._reslist.index(self.shx.atoms.get_atom_by_name('C99'))
        hklf_idx = next(i for i, x in enumerate(self.shx._reslist) if isinstance(x, HKLF))
        self.assertLess(atom_idx, hklf_idx)

    def test_add_atom_placed_after_last_real_atom(self):
        # The last non-Q-peak atom before adding:
        existing_atoms = [x for x in self.shx._reslist if isinstance(x, Atom) and not x.qpeak]
        last_existing_pos = self.shx._reslist.index(existing_atoms[-1])
        self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        new_atom_pos = self.shx._reslist.index(self.shx.atoms.get_atom_by_name('C99'))
        self.assertGreater(new_atom_pos, last_existing_pos)

    # ------------------------------------------------------------------
    # after= parameter
    # ------------------------------------------------------------------

    def test_add_atom_after_specific_atom(self):
        anchor = self.shx.atoms.get_atom_by_name('C1')
        anchor_pos = self.shx._reslist.index(anchor)
        self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3], after=anchor)
        new_pos = self.shx._reslist.index(self.shx.atoms.get_atom_by_name('C99'))
        self.assertEqual(anchor_pos + 1, new_pos)

    # ------------------------------------------------------------------
    # uvals handling
    # ------------------------------------------------------------------

    def test_add_atom_uvals_default_isotropic(self):
        a = self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        self.assertEqual([0.04, 0.0, 0.0, 0.0, 0.0, 0.0], a.uvals)

    def test_add_atom_uvals_single_value_expanded(self):
        a = self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3], uvals=[0.05])
        self.assertEqual([0.05, 0.0, 0.0, 0.0, 0.0, 0.0], a.uvals)

    def test_add_atom_uvals_anisotropic(self):
        uvals = [0.03, 0.04, 0.05, 0.001, 0.002, 0.003]
        a = self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3], uvals=uvals)
        self.assertEqual(uvals, a.uvals)

    # ------------------------------------------------------------------
    # symmgen flag
    # ------------------------------------------------------------------

    def test_add_atom_not_symmgen(self):
        a = self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        self.assertFalse(a.symmgen)

    # ------------------------------------------------------------------
    # Element auto-registration
    # ------------------------------------------------------------------

    def test_add_atom_known_element_carbon(self):
        a = self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3], element='C')
        self.assertEqual('C', a.element)

    def test_add_atom_new_element_auto_added_to_sfac(self):
        # 'Xe' is not in p21c.res SFAC table
        self.assertFalse(self.shx.sfac_table.has_element('Xe'))
        self.shx.add_atom(name='XE1', coordinates=[0.1, 0.2, 0.3], element='Xe')
        self.assertTrue(self.shx.sfac_table.has_element('Xe'))

    def test_add_atom_new_element_sfac_num_nonzero(self):
        a = self.shx.add_atom(name='XE1', coordinates=[0.1, 0.2, 0.3], element='Xe')
        self.assertGreater(a.sfac_num, 0)
        self.assertEqual('Xe', a.element)

    # ------------------------------------------------------------------
    # Duplicate name validation
    # ------------------------------------------------------------------

    def test_add_atom_duplicate_name_raises(self):
        self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        with self.assertRaises(ValueError):
            self.shx.add_atom(name='C99', coordinates=[0.2, 0.3, 0.4])

    # ------------------------------------------------------------------
    # Cartesian coordinates
    # ------------------------------------------------------------------

    def test_add_atom_cartesian_coords(self):
        # Provide a clearly non-fractional Cartesian point; the atom should land in _reslist
        a = self.shx.add_atom(name='C99', coordinates=[1.0, 2.0, 3.0], coords_are_cartesian=True)
        # Fractional coords must all be in a sane range (not equal to 1.0/2.0/3.0 directly)
        self.assertFalse(a.x == 1.0 and a.y == 2.0 and a.z == 3.0)

    # ------------------------------------------------------------------
    # PART and SOF
    # ------------------------------------------------------------------

    def test_add_atom_with_part(self):
        a = self.shx.add_atom(name='C99', coordinates=[0.5, 0.5, 0.5], part=1, sof=21.0)
        self.assertEqual(1, a.part.n)
        self.assertEqual(21.0, a.sof)

    def test_add_atom_occupancy_parameter(self):
        """occupancy=0.5 with default fvar=1 → sof=10.5"""
        a = self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3], occupancy=0.5)
        self.assertAlmostEqual(10.5, a.sof)

    def test_add_atom_occupancy_with_fvar(self):
        """occupancy=1.0, fvar=2 → sof=21.0"""
        a = self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3],
                              occupancy=1.0, fvar=2)
        self.assertAlmostEqual(21.0, a.sof)

    def test_add_atom_occupancy_takes_priority_over_sof(self):
        """Mixing occupancy and sof raises ValueError."""
        with self.assertRaises(ValueError):
            self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3],
                              occupancy=0.5, sof=10.5)

    def test_add_atom_fvar_with_sof_raises(self):
        """Mixing fvar (non-default) and sof raises ValueError."""
        with self.assertRaises(ValueError):
            self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3],
                              fvar=2, sof=21.0)

    def test_add_atom_occupancy_half_disorder(self):
        """Typical 50:50 disorder: occupancy=0.5 for both parts."""
        a1 = self.shx.add_atom(name='C97', coordinates=[0.1, 0.2, 0.3],
                               occupancy=0.5, fvar=1, part=1)
        a2 = self.shx.add_atom(name='C98', coordinates=[0.1, 0.2, 0.31],
                               occupancy=0.5, fvar=1, part=2)
        self.assertAlmostEqual(10.5, a1.sof)
        self.assertAlmostEqual(10.5, a2.sof)
        self.assertEqual(1, a1.part.n)
        self.assertEqual(2, a2.part.n)

    # ------------------------------------------------------------------
    # Written output is valid
    # ------------------------------------------------------------------

    def test_add_atom_written_to_file(self):
        self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3])
        with tempfile.NamedTemporaryFile(suffix='.res', delete=False) as tmp:
            tmppath = Path(tmp.name)
        try:
            self.shx.write_shelx_file(tmppath)
            content = tmppath.read_text()
            self.assertIn('C99', content)
            self.assertIn('TITL', content)
            self.assertIn('HKLF', content)
        finally:
            tmppath.unlink(missing_ok=True)

    # ------------------------------------------------------------------
    # unused_atom_name
    # ------------------------------------------------------------------

    def test_unused_atom_name_carbon(self):
        name = self.shx.unused_atom_name('C')
        # It must be a valid, non-existing name
        self.assertFalse(self.shx.atoms.has_atom(name))
        self.assertTrue(name.startswith('C'))

    def test_unused_atom_name_new_element(self):
        name = self.shx.unused_atom_name('Xe')
        self.assertEqual('Xe1', name)

    def test_unused_atom_name_max_4_chars(self):
        """Every generated name must fit within SHELX's 4-character limit."""
        for element in ('C', 'Fe', 'Al'):
            name = self.shx.unused_atom_name(element)
            self.assertLessEqual(len(name), 4, f"Name {name!r} exceeds 4 chars")

    def test_unused_atom_name_two_char_element_limit(self):
        """For a 2-char element like 'Fe', max suffix is 99 (Fe99 = 4 chars)."""
        # Fill Fe1..Fe99 artificially by patching nameslist lookup
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        # Add Fe1..Fe99 to the atoms dict cache to simulate exhaustion
        fake_names = {f'FE{n}_0' for n in range(1, 100)}  # FE1_0 … FE99_0
        original_nameslist = shx.atoms.nameslist
        # Monkey-patch nameslist for this one test
        from unittest.mock import patch
        with patch.object(type(shx.atoms), 'nameslist',
                          new_callable=lambda: property(lambda self: original_nameslist + tuple(fake_names))):
            with self.assertRaises(ValueError):
                shx.unused_atom_name('Fe')

    def test_unused_atom_name_returns_unique_name(self):
        name1 = self.shx.unused_atom_name('C')
        self.shx.add_atom(name=name1, coordinates=[0.1, 0.2, 0.3])
        name2 = self.shx.unused_atom_name('C')
        self.assertNotEqual(name1, name2)
        self.assertFalse(self.shx.atoms.has_atom(name2))

    def test_unused_atom_name_combines_with_add_atom(self):
        """Round-trip: get a name then add an atom with it."""
        name = self.shx.unused_atom_name('N')
        a = self.shx.add_atom(name=name, coordinates=[0.3, 0.3, 0.3], element='N')
        self.assertEqual(name, a.name)
        self.assertIsNotNone(self.shx.atoms.get_atom_by_name(name))


