from fractions import Fraction
from math import acos, sqrt, degrees
from typing import Union, List, TYPE_CHECKING, Iterator, Dict

import numpy as np

if TYPE_CHECKING:
    from shelxfile import Shelxfile
from shelxfile.atoms.atom import Atom
from shelxfile.atoms.pairs import Bond, SymBond
from shelxfile.misc.dsrmath import atomic_distance, Array
from shelxfile.misc.misc import build_conntable


def _build_symm_label(symm_element, h: int, k: int, l: int) -> str:
    """Build a human-readable symmetry label in SHELXL/CIF style.

    Combines the rotation/translation from *symm_element* with the extra
    integer lattice shifts *h*, *k*, *l* into a string such as
    ``'-x+1, y+1/2, -z'``.

    Parameters
    ----------
    symm_element : SymmetryElement
        The symmetry operation from ``Shelxfile.symmcards``.
    h, k, l : int
        Extra integer translations (negative values shift in the opposite
        direction), corresponding to the lattice vector that places the
        generated atom nearest to the reference atom.
    """
    axes = ['x', 'y', 'z']
    parts = []
    for i, shift in enumerate([h, k, l]):
        total_t = float(symm_element.trans[i]) + shift
        row = symm_element.matrix[i]
        term = ''
        for j in range(3):
            c = float(row[j])
            if c == 0.0:
                continue
            if c == 1.0:
                term += f'+{axes[j]}'
            elif c == -1.0:
                term += f'-{axes[j]}'
            elif c > 0:
                term += f'+{c:g}{axes[j]}'
            else:
                term += f'{c:g}{axes[j]}'
        if total_t:
            frac = Fraction(total_t).limit_denominator(12)
            if frac.denominator == 1:
                t_str = f'{int(frac):+d}'
            else:
                sign = '+' if frac > 0 else ''
                t_str = f'{sign}{frac}'
        else:
            t_str = ''
        parts.append((term + t_str).lstrip('+'))
    return ', '.join(parts)


class Atoms():
    """
    All atoms from a SHELXL file with their properties.
    """

    def __init__(self, shx: 'Shelxfile'):
        self.shx = shx
        self.all_atoms: List[Atom] = []
        self._atomsdict: Dict[str, Atom] = {}

    def append(self, atom: 'Atom') -> None:
        """
        Adds a new atom to the list of atoms. Using append is essential.
        """
        self.all_atoms.append(atom)
        self._atomsdict.clear()

    @property
    def nameslist(self) -> tuple[str, ...]:
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
                self._atomsdict.clear()
        # if self.shx.debug:
        #    print('Could not delete atom {}'.format(self.get_atom_by_id(key.atomid).fullname))

    @property
    def atomsdict(self):
        if not self._atomsdict:
            self._atomsdict = dict((atom.fullname.upper(), atom) for atom in self.all_atoms)
        return self._atomsdict

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
        # if not atom and self.shx.debug:
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

    @property
    def conntable(self) -> tuple:
        """Connectivity table for all atoms in the asymmetric unit.

        Returns a tuple of ``(i, j)`` index pairs (i < j) where atoms *i* and
        *j* are considered covalently bonded based on the sum of their covalent
        radii multiplied by 1.2.  Disorder-part rules and H–H bonds are
        handled automatically.

        The indices correspond to positions in :attr:`all_atoms`.

        Returns
        -------
        tuple of (int, int)
        """
        atoms = self.all_atoms
        if not atoms:
            return ()
        coords = np.array([[at.xc, at.yc, at.zc] for at in atoms], dtype=np.float64)
        types = [at.element for at in atoms]
        parts = [at.part.n for at in atoms]
        radii = np.array([at.radius for at in atoms], dtype=np.float64)
        symmgen = [at.symmgen for at in atoms]
        return build_conntable(coords, types, parts, radii=radii, symmgen=symmgen)

    @property
    def bonds(self) -> List[Bond]:
        """All covalent bonds in the asymmetric unit as a human-readable list.

        Each entry is a :class:`~shelxfile.atoms.pairs.Bond` with ``atom1``,
        ``atom2``, and ``distance`` (Å) attributes.  The list is sorted by
        atom1 name, then atom2 name.

        Bonds are detected via the same covalent-radii heuristic used by
        :attr:`conntable` (H–H pairs, cross-disorder and symmetry-generated
        atoms are excluded).

        Example::

            for bond in shx.atoms.bonds:
                print(bond)
            # C1     – C2      1.5210 Å
            # C1     – H1      1.0900 Å
            # ...

            # Unpack into components:
            atom1, atom2, dist = shx.atoms.bonds[0]

        Returns
        -------
        list of :class:`~shelxfile.atoms.pairs.Bond`
        """
        atoms = self.all_atoms
        result: List[Bond] = []
        for i, j in self.conntable:
            a1, a2 = atoms[i], atoms[j]
            dist = float(np.sqrt(
                (a1.xc - a2.xc) ** 2 +
                (a1.yc - a2.yc) ** 2 +
                (a1.zc - a2.zc) ** 2
            ))
            result.append(Bond(a1, a2, dist))
        result.sort(key=lambda b: (b.atom1.fullname_short, b.atom2.fullname_short))
        return result

    def full_bond_list(self, with_qpeaks: bool = False) -> List[SymBond]:
        """Bond list covering every atom in the asymmetric unit and all its
        crystallographic neighbors, including symmetry-generated ones.

        For each atom *A* in the asymmetric unit the list contains one
        :class:`~shelxfile.atoms.pairs.SymBond` per bonded neighbor:

        * **Plain bond** — neighbor *B* is also in the asymmetric unit (same
          unit cell, identity operation).  Stored once with an empty
          ``symm_label``; *A* is chosen as the atom with the lower index so
          duplicate entries are suppressed.
        * **Symmetry bond** — the bonded image of *B* is generated by applying
          a crystallographic symmetry operation (and possibly a lattice
          translation) to the coordinates of *B* in the asymmetric unit.  The
          label uses SHELXL / CIF notation, e.g. ``'-x, y+1/2, -z+1/2'``.
          Both (A → image-of-B) and (B → image-of-A) appear, mirroring
          SHELXL's ``.lst`` bond table where each atom lists all its neighbors.

        The list is sorted by (atom1 name, atom2 name).

        The method runs :meth:`~shelxfile.shelx.sdm.SDM.calc_sdm` internally;
        it is O(N²·S) where *N* is the number of atoms and *S* is the number
        of symmetry operations.

        Parameters
        ----------
        with_qpeaks : bool
            Include Q-peaks (difference-map peaks) in the bond detection.
            Defaults to ``False``.

        Returns
        -------
        list of :class:`~shelxfile.atoms.pairs.SymBond`

        Example::

            for bond in shx.atoms.full_bond_list():
                print(bond)
            # AL1    – O1      1.7236 Å
            # AL1    – O2      1.7095 Å   #2 [-x, y+1/2, -z+1/2]
            # ...

            # Unpack all four fields:
            a1, a2, dist, label = shx.atoms.full_bond_list()[0]

            # Only symmetry bonds:
            sym_bonds = [b for b in shx.atoms.full_bond_list() if b.is_symmetry_bond]
        """
        # Late import avoids circular dependency: SDM imports Atom but not Atoms.
        from shelxfile.shelx.sdm import SDM

        sdm = SDM(self.shx)
        sdm.calc_sdm()

        seen_plain: set = set()  # deduplicate intra-asym-unit bonds
        result: List[SymBond] = []

        for item in sdm.sdm_list:
            if not item.covalent:
                continue
            if not with_qpeaks and (item.atom1.qpeak or item.atom2.qpeak):
                continue

            # --- compute lattice translation ---
            x1, y1, z1 = item.atom1.x, item.atom1.y, item.atom1.z
            x2, y2, z2 = item.atom2.x, item.atom2.y, item.atom2.z
            n = item.symmetry_number
            symm = self.shx.symmcards[n]
            m = symm.matrix  # 3×3 numpy array
            t = symm.trans  # 3-element array

            # Position of atom1's image under symmetry op n:
            px = m[0, 0] * x1 + m[0, 1] * y1 + m[0, 2] * z1 + t[0]
            py = m[1, 0] * x1 + m[1, 1] * y1 + m[1, 2] * z1 + t[1]
            pz = m[2, 0] * x1 + m[2, 1] * y1 + m[2, 2] * z1 + t[2]

            # Integer lattice shift that places the image nearest to atom2:
            h = round(px - x2)
            k = round(py - y2)
            l = round(pz - z2)

            is_plain = (n == 0 and h == 0 and k == 0 and l == 0)

            if is_plain:
                # Deduplicate: keep only the (lower-index, higher-index) pair.
                key = (min(item.a1, item.a2), max(item.a1, item.a2))
                if key in seen_plain:
                    continue
                seen_plain.add(key)
                # Canonical ordering: lower index atom is atom1.
                if item.a1 < item.a2:
                    a_center, a_neighbor = item.atom1, item.atom2
                else:
                    a_center, a_neighbor = item.atom2, item.atom1
                result.append(SymBond(a_center, a_neighbor, item.dist, ''))
            else:
                # Symmetry bond: atom2 bonds to image of atom1 via symm n + (h,k,l).
                symm_label = _build_symm_label(symm, h, k, l)
                result.append(SymBond(item.atom2, item.atom1, item.dist, symm_label, symm_number=n))

        result.sort(key=lambda b: (b.atom1.fullname_short, b.atom2.fullname_short))
        return result

    def atoms_in_class(self, name: str) -> list:
        """
        Returns a list of atoms in residue class 'name'
        """
        atoms = []
        for x in self.all_atoms:
            if x.resiclass == name and x.name not in atoms:
                atoms.append(x.name)
        return atoms
