from math import acos, sqrt, degrees
from typing import Union, List, TYPE_CHECKING, Tuple, Iterator

if TYPE_CHECKING:
    from shelxfile import Shelxfile
from shelxfile.atoms.atom import Atom
from shelxfile.misc.dsrmath import atomic_distance, Array


class Atoms():
    """
    All atoms from a SHELXL file with their properties.
    """

    def __init__(self, shx: 'Shelxfile'):
        self.shx = shx
        self.all_atoms: List[Atom] = []

    def append(self, atom: 'Atom') -> None:
        """
        Adds a new atom to the list of atoms. Using append is essential.
        """
        self.all_atoms.append(atom)

    @property
    def nameslist(self) -> Tuple[str]:
        return tuple([at.fullname.upper() for at in self.all_atoms])

    def __repr__(self) -> str:
        if self.all_atoms:
            return '\n'.join([str(x) for x in self.all_atoms])
        else:
            return 'No Atoms in file.'

    def __iter__(self) -> Iterator:
        return iter(x for x in self.all_atoms)

    def __getitem__(self, item: int) -> 'Atom':
        return self.get_atom_by_id(item)

    def __len__(self) -> int:
        return len(self.all_atoms)

    def __delitem__(self, key: int) -> None:
        """
        Delete an atom by its atomid:
        del atoms[4]
        """
        for n, at in enumerate(self.all_atoms):
            if key == at.atomid:
                if self.shx.debug:
                    print("deleting atom", at.fullname)
                del self.all_atoms[n]
                del self.shx._reslist[self.shx._reslist.index(at)]
        # if self.shx.debug:
        #    print('Could not delete atom {}'.format(self.get_atom_by_id(key.atomid).fullname))

    @property
    def atomsdict(self):
        return dict((atom.fullname, atom) for atom in self.all_atoms)

    @property
    def number(self) -> int:
        """
        The number of atoms in the current SHELX file.
        """
        return len(self.all_atoms)

    def get_atom_by_id(self, aid: int) -> Union['Atom', None]:
        """
        Returns the atom objext with atomId id.
        """
        for a in self.all_atoms:
            if aid == a.atomid:
                return a

    def has_atom(self, atom_name: str) -> bool:
        """
        Returns true if shelx file has atom.
        """
        if '_' not in atom_name:
            atom_name += '_0'
        if atom_name.upper() in self.nameslist:
            return True
        else:
            return False

    def get_atom_by_name(self, atom_name: str) -> Union['Atom', None]:
        """
        Returns an Atom object using an atom name with residue number like C1, C1_0, F2_4, etc.
        C1 means atom C1 in residue 0.
        """
        if '_' not in atom_name:
            if atom_name == ">" or atom_name == "<":
                return None
            atom_name = f'{atom_name}_0'
        atom = self.atomsdict.get(atom_name.upper(), None)
        #if not atom and self.shx.debug:
        #    print(f"Atom {atom_name} not found in atom list.")
        return atom

    def get_multi_atnames(self, atom_name, residue_class):
        atoms = []
        if residue_class:
            for num in self.shx.residues.residue_classes[residue_class]:
                if '_' not in atom_name:
                    atom_name += '_0'
                else:
                    atom_name += '_{}'.format(num)
                try:
                    atoms.append(self.atomsdict[atom_name.upper()])
                except KeyError:
                    pass
        else:
            try:
                atoms.append(self.atomsdict[atom_name.upper()])
            except KeyError:
                return None
        return atoms

    def get_all_atomcoordinates(self) -> dict:
        """
        Returns a dictionary {'C1': ['1.123', '0.7456', '3.245'], 'C2_2': ...}
        """
        atdict = {}
        for at in self.all_atoms:
            # if at.qpeak:
            #    atdict[at.name] = at.frac_coords
            # else:
            atdict[at.name.upper() + '_' + str(at.resinum)] = at.frac_coords
        return atdict

    def get_frag_fend_atoms(self) -> list:
        """
        Returns a list of atoms with cartesian coordinates. Atom names and sfac are ignored. They come from AFIX 17x.
        [[0.5316439256202359, 7.037351406500001, 10.112963255220803],
        [-1.7511017452002604, 5.461541059000001, 10.01187984858907]]
        """
        atoms = []
        for at in self.all_atoms:
            if at.frag_atom:
                atoms.append([at.xc, at.yc, at.zc])
        return atoms

    @property
    def hydrogen_atoms(self) -> List[Atom]:
        return [x for x in self.shx.atoms.all_atoms if x.is_hydrogen]

    @property
    def n_hydrogen_atoms(self) -> int:
        return len(self.hydrogen_atoms)

    @property
    def n_anisotropic_atoms(self) -> int:
        return len([x for x in self.all_atoms if sum(x.uvals[1:]) > 0.00001])

    @property
    def n_isotropic_atoms(self) -> int:
        return len([x for x in self.all_atoms if sum(x.uvals[1:]) == 0.0])

    @property
    def n_anisotropic_hydrogen_atoms(self) -> int:
        return len([x for x in self.hydrogen_atoms if sum(x.uvals[1:]) > 0.0001])

    @property
    def n_hydrogen_atoms_with_constr_u_val(self) -> int:
        return len([x for x in self.hydrogen_atoms if x.uvals[0] < -1.0])

    @property
    def riding_atoms(self) -> List[Atom]:
        return [x for x in self.hydrogen_atoms if x.afix]

    @property
    def residues(self) -> list:
        """
        Returns a list of the residue numbers in the shelx file.
        """
        return list(set([x.resinum for x in self.all_atoms]))

    @property
    def q_peaks(self) -> list:
        r"""
        Returns a list of q-peaks in the file.
        """
        return [x for x in self.all_atoms if x.qpeak]

    def distance(self, atom1: str, atom2: str) -> float:
        """
        Calculates the (shortest) distance of two atoms given as text names e.g. C1_3.
        """
        a1 = self.get_atom_by_name(atom1)
        a2 = self.get_atom_by_name(atom2)
        try:
            return atomic_distance([a1.xc, a1.yc, a1.zc], [a2.xc, a2.yc, a2.zc])
        except AttributeError:
            return 0.0

    def angle(self, at1: 'Atom', at2: 'Atom', at3: 'Atom') -> float:
        """
        Calculates the angle between three atoms.
        """
        ac1 = Array(at1.cart_coords)
        ac2 = Array(at2.cart_coords)
        ac3 = Array(at3.cart_coords)
        vec1 = ac2 - ac1
        vec2 = ac2 - ac3
        return vec1.angle(vec2)

    def torsion_angle(self, at1: 'Atom', at2: 'Atom', at3: 'Atom', at4: 'Atom') -> float:
        """
        Calculates the torsion angle (dieder angle) between four atoms.

        From the book of Camelo Giacovazzo:
        For a sequence of four atoms A, B, C, D, the torsion angle w(ABCD) is
        defined as the angle between the normals to the planes ABC and BCD.
        By convention w is positive if the sense of rotation from BA to
        CD, viewed down BC, is clockwise, otherwise it is negative.
        """
        ac1 = Array(at1.cart_coords)
        ac2 = Array(at2.cart_coords)
        ac3 = Array(at3.cart_coords)
        ac4 = Array(at4.cart_coords)
        # Three vectors between four atoms:
        v1 = ac2 - ac1
        v2 = ac3 - ac2
        v3 = ac4 - ac3
        # cross product:
        a = v1.cross(v2)
        b = v2.cross(v3)
        # If direction > 0, angle is positive, else negative:
        direction = v1[0] * v2[1] * v3[2] - v1[2] * v1[1] * v3[0] + v1[2] * v2[0] * v3[1] - v1[0] \
                    * v2[2] * v3[1] + v1[1] * v2[2] * v3[0] - v1[1] * v2[0] * v3[2]
        # angle between plane normals:
        ang = acos((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (
                sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) * sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2])))
        return degrees(ang) if direction > 0 else degrees(-ang)

    def atoms_in_class(self, name: str) -> list:
        """
        Returns a list of atoms in residue class 'name'
        """
        atoms = []
        for x in self.all_atoms:
            if x.resiclass == name and x.name not in atoms:
                atoms.append(x.name)
        return atoms
