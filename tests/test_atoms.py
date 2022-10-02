from unittest import TestCase

from shelxfile.misc.misc import frac_to_cart
from shelxfile.shelx.shelx import Shelxfile


class TestAtoms(TestCase):

    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_u_eq_of_hydrogen_atom_in_single_afix(self):
        h = self.shx.atoms.get_atom_by_name('H34')
        self.assertEqual('C34', h.pivot.name)
        self.assertEqual(0.02959929, round(h.Uiso, 8))
        self.assertEqual(h.pivot.Uiso * 1.2, h.Uiso)

    def test_u_eq_of_hydrogen_atom_in_afix_137(self):
        h1 = self.shx.atoms.get_atom_by_name('H36A')
        h2 = self.shx.atoms.get_atom_by_name('H36B')
        h3 = self.shx.atoms.get_atom_by_name('H36C')
        self.assertEqual('C36', h1.pivot.name)
        self.assertEqual('C36', h2.pivot.name)
        self.assertEqual('C36', h2.pivot.name)
        self.assertEqual(0.05123489, round(h1.Uiso, 8))
        self.assertEqual(0.05123489, round(h2.Uiso, 8))
        self.assertEqual(0.05123489, round(h3.Uiso, 8))

    def test_uiso_of_pivot_of_H36x(self):
        c36 = self.shx.atoms.get_atom_by_name('C36')
        self.assertEqual(None, c36.pivot)
        self.assertEqual(0.03415659, round(c36.Uiso, 8))
        self.assertEqual(0.03415659 * 1.5, 0.051234885)

    def test_number(self):
        self.assertEqual(148, self.shx.atoms.number)

    def test_has_atom_al1(self):
        self.assertEqual(True, self.shx.atoms.has_atom('Al1'))

    def test_has_atom_al1_0(self):
        self.assertEqual(True, self.shx.atoms.has_atom('Al1_0'))

    def test_has_atom_al2_0(self):
        self.assertEqual(False, self.shx.atoms.has_atom('Al2_0'))

    def test_get_atom_by_name(self):
        self.assertEqual("Atom ID: 73", self.shx.atoms.get_atom_by_name('Al1').__repr__())

    def test_get_all_atomcoordinates(self):
        coordinate = {'O1_4': (0.074835, 0.238436, 0.402457), 'C1_4': (0.028576, 0.234542, 0.337234)}
        self.assertEqual(coordinate['O1_4'], self.shx.atoms.get_all_atomcoordinates()['O1_4'])
        self.assertEqual(coordinate['C1_4'], self.shx.atoms.get_all_atomcoordinates()['C1_4'])

    def test_residues(self):
        self.assertEqual([0, 1, 2, 3, 4], self.shx.atoms.residues)

    def test_q_peaks(self):
        self.assertEqual(['Atom ID: 328', 'Atom ID: 329', 'Atom ID: 330', 'Atom ID: 331', 'Atom ID: 332'],
                         [x.__repr__() for x in self.shx.atoms.q_peaks[:5]])

    def test_distance(self):
        self.assertEqual(2.154399, round(self.shx.atoms.distance('F1_2', 'F2_2'), 6))

    def test_distance2(self):
        self.assertEqual(1.332854, round(self.shx.atoms.distance('C2_2', 'F1_2'), 6))

    def test_angle(self):
        at1 = self.shx.atoms.get_atom_by_name('O1_4')
        at2 = self.shx.atoms.get_atom_by_name('C1_4')
        at3 = self.shx.atoms.get_atom_by_name('C2_4')
        # The angle between the three atoms is ... degree:
        self.assertEqual(109.688123, round(self.shx.atoms.angle(at1, at2, at3), 6))

    def test_torsion_angle(self):
        at1 = self.shx.atoms.get_atom_by_name('O1')
        at2 = self.shx.atoms.get_atom_by_name('C1')
        at3 = self.shx.atoms.get_atom_by_name('C2')
        at4 = self.shx.atoms.get_atom_by_name('F1')
        # The torsion angle between the four atoms is ... degree:
        self.assertEqual(74.095731, round(self.shx.atoms.torsion_angle(at1, at2, at3, at4), 6))

        # Thisis the other way round:
        at4 = self.shx.atoms.get_atom_by_name('F2')
        self.assertEqual(-44.467358, round(self.shx.atoms.torsion_angle(at1, at2, at3, at4), 6))

    def test_atoms_in_class(self):
        atoms = ['O1', 'C1', 'C2', 'F1', 'F2', 'F3', 'C3', 'F4', 'F5', 'F6', 'C4', 'F7', 'F8', 'F9']
        self.assertEqual(atoms, self.shx.atoms.atoms_in_class('CCF3'))

    def test_sfac_num(self):
        at = self.shx.atoms.get_atom_by_name('C1_4')
        self.assertEqual(1, at.sfac_num)

    def test_element_c(self):
        at = self.shx.atoms.get_atom_by_name('C1_4')
        self.assertEqual('C', at.element)

    def test_element_o(self):
        at = self.shx.atoms.get_atom_by_name('C1_4')
        at.element = 'O'
        self.assertEqual('O', at.element)
        self.assertEqual(3, at.sfac_num)

    def test_an(self):
        at = self.shx.atoms.get_atom_by_name('C1_4')
        self.assertEqual(6, at.an)

    def test_delete_show(self):
        # We have an atom with ID = 40
        at = self.shx.atoms.get_atom_by_id(40)
        self.assertEqual('Atom ID: 40', at.__repr__())

    def test_delete_show_list(self):
        # Before and after atom 40 are 38 and 42
        self.assertEqual(['Atom ID: 38', 'Atom ID: 40', 'Atom ID: 42'],
                         [x.__repr__() for x in self.shx.atoms.all_atoms[:3]])

    def test_c1_4(self):
        # The name of atom 40 is C1_4
        self.assertEqual('C1_4', self.shx.atoms.get_atom_by_id(40).fullname)

    def test_delete_atom40(self):
        # get atom 40
        at = self.shx.atoms.get_atom_by_id(40)
        # Now delete atom 40
        at.delete()
        # The resulting atoms list does not contain atom 40
        new_atoms = ['Atom ID: 38', 'Atom ID: 41', 'Atom ID: 43']
        self.assertEqual(new_atoms, [x.__repr__() for x in self.shx.atoms.all_atoms[:3]])

    def test_delete2(self):
        self.shx.atoms.get_atom_by_id(40).delete()
        self.assertEqual(None, self.shx.atoms.get_atom_by_id(40))

    def test_delete3(self):
        self.shx.atoms.get_atom_by_id(40).delete()
        self.assertEqual(None, self.shx.atoms.get_atom_by_name('C1_4'))

    def test_equal(self):
        self.assertEqual(True, self.shx.atoms.get_atom_by_id(40) == self.shx.atoms.get_atom_by_name('C1_4'))

    def test_find_atoms_around(self):
        at = self.shx.atoms.get_atom_by_name('C1_4')
        found = [x.__repr__() for x in at.find_atoms_around(dist=2, only_part=2)]
        shouldfind = ['Atom ID: 38', 'Atom ID: 42', 'Atom ID: 50', 'Atom ID: 58']
        self.assertEqual(shouldfind, found)

    def test_get_coordinates(self):
        c = self.shx.atoms.get_atom_by_name('C1_4').cart_coords
        self.assertEqual([-0.19777464582151, 4.902748697, 6.89776640065678], [round(x, 14) for x in c])

    def test_to_isotropic(self):
        # We have regular u values of atom 40
        self.assertEqual([0.02311, 0.03617, 0.01096, -0.01, 0.00201, 0.00356], self.shx.atoms.get_atom_by_id(40).uvals)
        # We make them isotropic
        self.shx.atoms.get_atom_by_id(40).to_isotropic()
        # And the result is only one atomic u value:
        self.assertEqual([0.04, 0.0, 0.0, 0.0, 0.0, 0.0], self.shx.atoms.get_atom_by_id(40).uvals)

    def test_cart_coords(self):
        self.assertEqual([-0.19777464582151, 4.902748697, 6.89776640065678],
                         [round(x, 14) for x in self.shx.atoms.get_atom_by_id(40).cart_coords])

    def test_cartesian_from_method(self):
        self.assertEqual(
            [round(x, 13) for x in frac_to_cart(self.shx.atoms.get_atom_by_id(40).frac_coords, list(self.shx.cell))],
            [round(x, 13) for x in self.shx.atoms.get_atom_by_id(40).cart_coords])

    def test_frac_coords(self):
        self.assertEqual((0.028576, 0.234542, 0.337234), self.shx.atoms.get_atom_by_id(40).frac_coords)

    def test_frac_x(self):
        self.assertEqual(0.028576, self.shx.atoms.get_atom_by_id(40).x)

    def test_cart_x(self):
        self.assertEqual(-0.1977746458, round(self.shx.atoms.get_atom_by_id(40).xc, 10))

    def test_radius(self):
        self.assertEqual(0.77, self.shx.atoms.get_atom_by_id(40).radius)

    def test_resinum(self):
        self.assertEqual(4, self.shx.atoms.get_atom_by_id(40).resinum)

    def test_resiclass(self):
        self.assertEqual('CCF3', self.shx.atoms.get_atom_by_id(40).resiclass)

    def test_fullname(self):
        self.assertEqual('C1_4', self.shx.atoms.get_atom_by_id(40).fullname)

    def test_part(self):
        self.assertEqual('PART 2 -31', self.shx.atoms.get_atom_by_id(40).part.__repr__())

    def test_part_number(self):
        self.assertEqual(2, self.shx.atoms.get_atom_by_id(40).part.n)

    def test_part_occupancy(self):
        self.assertEqual(0.44236, self.shx.atoms.get_atom_by_id(40).occupancy)

    def test_part_sof(self):
        self.assertEqual(-31.0, self.shx.atoms.get_atom_by_id(40).part.sof)

    def test_get_pivot_atom(self):
        self.assertEqual('C34', self.shx.atoms.all_atoms[45].pivot.name)

    def test_get_pivot_atom_of_fluorine(self):
        # This atom should not have a pivot atom, because it is heavy and not in an AFIX:
        self.assertEqual('F2', self.shx.atoms.all_atoms[4].name)
        self.assertEqual(None, self.shx.atoms.all_atoms[4].pivot)

    def test_hydrogen_atoms(self):
        self.assertEqual('[Atom ID: 134, Atom ID: 141, Atom ID: 148]', str(self.shx.atoms.hydrogen_atoms[:3]))

    def test_riding_atoms(self):
        self.assertEqual('[Atom ID: 134, Atom ID: 141, Atom ID: 148]', str(self.shx.atoms.riding_atoms[:3]))


class TestRidingAtoms(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/sad-final.res')

    def test_hydrogen_atoms(self):
        # Returns the list of hydrogen atoms, regardless of their model
        self.assertEqual('[Atom ID: 71, Atom ID: 75, Atom ID: 76]', str(self.shx.atoms.hydrogen_atoms[3:6]))

    def test_riding_atoms(self):
        # Returns only riding atoms, therefore different IDs than above
        self.assertEqual('[Atom ID: 75, Atom ID: 76, Atom ID: 77]', str(self.shx.atoms.riding_atoms[3:6]))

    def test_number_of_anisotropic_atoms(self):
        self.assertEqual(79, self.shx.atoms.n_anisotropic_atoms)

    def test_number_of_isotropic_atoms(self):
        self.assertEqual(49, self.shx.atoms.n_isotropic_atoms)

    def test_number_of_hydrogen_atoms(self):
        self.assertEqual(49, self.shx.atoms.n_hydrogen_atoms)

    def test_number_of_anisotropic_hydrogen_atoms(self):
        self.assertEqual(0, self.shx.atoms.n_anisotropic_hydrogen_atoms)

    def test_number_of_hydrogen_adtoms_with_constrained_u_values(self):
        self.assertEqual(49, self.shx.atoms.n_hydrogen_atoms_with_constr_u_val)


class TestUisoOfFreeRefine(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/2240189.res')

    def test_free_refined_hydrogen(self):
        h1a = self.shx.atoms.get_atom_by_name('H1A')
        self.assertEqual(0.04654, h1a.Uiso)
        self.assertEqual([0.04654, 0.0, 0.0, 0.0, 0.0, 0.0], h1a.uvals)
        self.assertEqual("O3'", h1a.pivot.name)
