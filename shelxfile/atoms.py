from typing import List, Any

from dsrmath import atomic_distance, frac_to_cart
from misc import DEBUG, split_fvar_and_parameter, ParseUnknownParam


class Atoms():
    """
    All atoms from a SHELXL file with their properties.
    """

    def __init__(self, shx):
        self.shx = shx
        self.atoms = []
        self.atomsdict = {}
        self.nameslist = []

    def append(self, atom: 'Atom') -> None:
        """
        Adds a new atom to the list of atoms. Using append is essential.
        """
        self.atoms.append(atom)
        name = atom.name + '_{}'.format(atom.resinum)
        self.atomsdict[name] = atom
        self.nameslist.append(name.upper())

    def __repr__(self):
        if self.atoms:
            return '\n'.join([str(x) for x in self.atoms])
        else:
            return 'No Atoms in file.'

    def __iter__(self):
        for x in self.atoms:
            yield x

    def __getitem__(self, item: int) -> None:
        return self.atoms[item]

    def __len__(self) -> int:
        return len(self.atoms)

    def __delitem__(self, key):
        """
        Delete an atom by its atomid:
        del atoms[4]
        """
        for n, at in enumerate(self.atoms):
            if key == at.atomid:
                if DEBUG:
                    print("deleting atom", at.name)
                del self.atoms[n]
                del self.atomsdict[at.name + '_{}'.format(at.resinum)]
                del self.nameslist[self.nameslist.index(at.fullname.upper())]
                for x in at._line_numbers:
                    del self.shx._reslist[x]

    @property
    def number(self) -> int:
        """
        The number of atoms in the current SHELX file.
        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.number
        148
        """
        return len(self.atoms)

    def get_atom_by_id(self, aid: int) -> 'Atom':
        """
        Returns the atom objext with atomId id.
        """
        for a in self.atoms:
            if aid == a.atomid:
                return a

    def has_atom(self, atom_name: str) -> bool:
        """
        Returns true if shelx file has atom.

        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.has_atom('Al1')
        True
        >>> shx.atoms.has_atom('Al1_0')
        True
        >>> shx.atoms.has_atom('Al2_0')
        False
        """
        if '_' not in atom_name:
            atom_name += '_0'
        if atom_name.upper() in self.nameslist:
            return True
        else:
            return False

    def get_atom_by_name(self, atom_name: str) -> 'Atom' or None:
        """
        Returns an Atom object using an atom name with residue number like C1, C1_0, F2_4, etc.
        C1 means atom C1 in residue 0.
        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.get_atom_by_name('Al1')
        Atom ID: 88
        """
        if '_' not in atom_name:
            atom_name += '_0'
        try:
            at = self.atomsdict[atom_name.upper()]
        except KeyError:
            print("Atom {} not found in atom list.".format(atom_name))
            return None
        return at

    def get_all_atomcoordinates(self) -> dict:
        """
        Returns a dictionary {'C1': ['1.123', '0.7456', '3.245'], 'C2_2': ...}
        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.get_all_atomcoordinates() # doctest: +ELLIPSIS
        {'O1_4': [0.074835, 0.238436, 0.402457], 'C1_4': [0.028576, 0.234542, 0.337234], ...}
        """
        atdict = {}
        for at in self.atoms:
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
        for at in self.atoms:
            if at.frag_atom:
                atoms.append([at.xc, at.yc, at.zc])
        return atoms

    @property
    def residues(self) -> list:
        """
        Returns a list of the residue numbers in the shelx file.
        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.residues
        [0, 1, 2, 3, 4]
        """
        return list({x.resinum for x in self.atoms})

    @property
    def q_peaks(self) -> list:
        r"""
        Returns a list of q-peaks in the file.
        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.q_peaks[:5] # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        [Atom ID: 346, Atom ID: 347, Atom ID: 348, Atom ID: 349, Atom ID: 350]
        """
        return [x for x in self.atoms if x.qpeak]

    def distance(self, atom1: str, atom2: str) -> float:
        """
        Calculates the (shortest) distance of two atoms given as text names e.g. C1_3.
        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> round(shx.atoms.distance('F1_2', 'F2_2'), 6)
        2.154399
        >>> round(shx.atoms.distance('C2_2', 'F1_2'), 6)
        1.332854
        """
        a1 = self.get_atom_by_name(atom1)
        a2 = self.get_atom_by_name(atom2)
        try:
            return atomic_distance([a1.xc, a1.yc, a1.zc], [a2.xc, a2.yc, a2.zc])
        except AttributeError:
            return 0.0

    def atoms_in_class(self, name: str) -> list:
        """
        Returns a list of atoms in residue class 'name'
        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.atoms_in_class('CCF3')
        ['O1', 'C1', 'C2', 'F1', 'F2', 'F3', 'C3', 'F4', 'F5', 'F6', 'C4', 'F7', 'F8', 'F9']
        """
        atoms = []
        for x in self.atoms:
            if x.resiclass == name:
                if x.name not in atoms:
                    atoms.append(x.name)
        return atoms


class Atom():
    """
    An Opbect holding all Properties of a shelxl atom plus some extra information like
    kartesian coordinates and element type.

    :type restraints: List[Restraint]
    """
    #                name    sfac     x         y        z       occ      u11      u12 ...
    _anisatomstr = '{:<4.4s}{:>3}{:>12.6f}{:>12.6f}{:>12.6f}{:>12.5f}{:>11.5f}{:>11.5f}' \
                   ' {:>12.5f}{:>11.5f}{:>11.5f}{:>11.5f}'  # Line wrap is handled during file write.
    #               name    sfac     x         y         z         occ      u11
    _isoatomstr = '{:<5.5s} {:<3}{:>10.6f}  {:>10.6f}  {:>9.6f}  {:>9.5f}  {:>9.5f}'
    _qpeakstr = '{:<5.5s} {:<3}{:>8.4f}  {:>8.4f}  {:>8.4f}  {:>9.5f}  {:<9.2f} {:<9.2f}'
    _fragatomstr = '{:<5.5s} {:>10.6f}  {:>10.6f}  {:>9.6f}'

    def __init__(self, shelx, spline: list, line_nums: list, line_number: int, part: int = 0,
                 afix: int = 0, residict: dict = None, sof: float = 0) -> None:
        # super(Atom, self).__init__(shelx)
        self._line_number = line_number
        self._lines = line_nums
        self.sfac_num = None
        self.name = None  # Name without residue number like "C1"
        self.fullname = None  # Name including residue nimber like "C1_2"
        # Site occupation factor including free variable like 31.0
        self.sof = None
        self.atomid = line_number
        self.shx = shelx
        self.element = None
        # fractional coordinates:
        self.x = None
        self.y = None
        self.z = None
        # cartesian coordinates:
        self.xc = None
        self.yc = None
        self.zc = None
        self.frag_atom = False
        self.restraints = []
        self.previous_non_h = self.shx.non_h
        if self.element not in ['H', 'D']:
            self.shx.non_h = line_number
        else:
            self.shx.non_h = None
        if sof:
            # sof defined from outside e.g. by PART 1 31
            self.sof = float(sof)
        elif len(spline) > 5:
            self.sof = float(spline[5])
        else:
            self.sof = 11.0
        # Only the occupancy of the atom *without* the free variable like 0.5
        fvar, self.occupancy = split_fvar_and_parameter(self.sof)
        self.shx.fvars.set_fvar_usage(fvar)
        self.uvals = [0.04]  # [u11 u12 u13 u21 u22 u23]
        self.resiclass = residict['class']
        if not residict['number']:
            self.resinum = 0  # all other atoms are residue 0
        else:
            self.resinum = residict['number']
        self.chain_id = residict['ID']
        self.part = part
        self.afix = afix
        self.qpeak = False
        self.peak_height = 0.0
        self.cell = shelx.cell
        self.parse_line(spline)
        if self.shx.frag:
            self.afix = int(self.shx.frag[0])  # The FRAG AFIX fit code like 176 (must be greater 16)
        if self.shx.anis:
            self.parse_anis()
        for n, u in enumerate(self.uvals):
            if abs(u) > 4.0:
                fvar, uval = split_fvar_and_parameter(u)
                self.uvals[n] = uval
                self.shx.fvars.set_fvar_usage(fvar)

    def parse_anis(self):
        """
        Parses the ANIS card. It can be either ANIS, ANIS name(s) or ANIS number.
        # TODO: Test if ANIS $CL works and if ANIS_* $C works
        """
        try:
            # ANIS with a number as parameter
            if int(self.shx.anis[1]) > 0:
                self.shx.anis[1] -= 1
                if len(self.uvals) < 6:
                    self.uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]
        except (TypeError, KeyError, ValueError, IndexError):
            # ANIS with a list of atoms
            if len(self.shx.anis) > 1:
                # if '_' in self.shx.anis[0]:
                #    resinum = self.shx.anis[0].upper().split('_')[1]
                for x in self.shx.anis[1:]:
                    if '_' in x:
                        name, resinum = x.upper().split('_')
                    else:
                        name = x.upper()
                        resinum = 0
                    if self.name == name and (int(self.resinum) == int(resinum) or resinum == '*'):
                        self.uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]
                        self.shx.anis.pop()
                        if self.shx.anis == ['ANIS']:
                            # ANIS finished, deactivating again:
                            self.shx.anis = None
                    if x.startswith('$'):
                        if name[1:].upper() == self.element \
                                and (int(self.resinum) == int(resinum) or resinum == '*'):
                            self.uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]
                            self.shx.anis.pop()
                            # TODO: This is a mess. Test and fix all sorts of ANIS possibilities.
                            if self.shx.anis == ['ANIS']:
                                # ANIS finished, deactivating again:
                                self.shx.anis = None
            # ANIS for all atoms
            else:
                if len(self.uvals) < 6:
                    self.uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]

    def parse_line(self, line):
        self.name = line[0][:4]
        self.fullname = self.name + '_{}'.format(self.resinum)
        uvals = [float(x) for x in line[6:12]]
        try:
            x, y, z = [float(x) for x in line[2:5]]
        except ValueError as e:
            if DEBUG:
                print(e, 'Line:', self._line_numbers[-1])
            raise ParseUnknownParam
        if abs(x) > 4:
            fvar, x = split_fvar_and_parameter(x)
            self.shx.fvars.set_fvar_usage(fvar)
        if abs(y) > 4:
            fvar, x = split_fvar_and_parameter(y)
            self.shx.fvars.set_fvar_usage(fvar)
        if abs(z) > 4:
            fvar, x = split_fvar_and_parameter(z)
            self.shx.fvars.set_fvar_usage(fvar)
        self.x = x
        self.y = y
        self.z = z
        self.uvals = uvals
        if len(self.uvals) == 2:
            self.peak_height = uvals.pop()
        if self.shx.end:  # After 'END' can only be Q-peaks!
            self.qpeak = True
        self.sfac_num = int(line[1])
        self.element = self.shx.sfac2elem(self.sfac_num).upper()
        self.xc, self.yc, self.zc = frac_to_cart([self.x, self.y, self.z], self.cell)

    def __iter__(self):
        for x in self.__repr__().split():
            yield x

    def __repr__(self) -> str:
        return 'Atom ID: ' + str(self.atomid)

    def __str__(self) -> str:
        """
        Returns a text line of the Atom with SHELXL syntax.
        :return: SHELX-formated atom string
        """
        if self.afix and self.shx.frag:
            # An atom from a FRAG/FEND instruction
            return Atom._fragatomstr.format(self.name, self.x, self.y, self.z)
        else:
            if len(self.uvals) > 2:
                # anisotropic atom
                try:
                    return Atom._anisatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof,
                                                    *self.uvals)
                except(IndexError):
                    return 'REM Error in U values.'
            else:
                # isotropic atom:
                if self.qpeak:
                    return Atom._qpeakstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, 0.04,
                                                 self.peak_height)
                try:
                    return Atom._isoatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof,
                                                   *self.uvals)
                except(IndexError):
                    return Atom._isoatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, 0.04)

    def resolve_restraints(self):
        for num, r in enumerate(self.shx.restraints):
            for at in r.atoms:
                print(r.residue_number, self.resinum, r.residue_class, self.resiclass, self.name, at)
                if r.residue_number == self.resinum and r.residue_class == self.resiclass and self.name == at:
                    self.restraints.append(r)

    @property
    def _line_numbers(self) -> list:
        # Line numbers (indexes) in the resfile.
        return self._lines

    @_line_numbers.setter
    def _line_numbers(self, value: list):
        self._lines = value

    @property
    def position(self):
        # The position in the res file.
        return self.shx._reslist.index(self)

    @property
    def frac_coords(self):
        return [self.x, self.y, self.z]

    @property
    def cart_coords(self):
        return [self.xc, self.yc, self.zc]

    def delete(self):
        """
        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms[:3]
        [Atom ID: 53, Atom ID: 55, Atom ID: 57]
        >>> shx._reslist[55:58]
        [Atom ID: 55, '', Atom ID: 57]
        >>> del shx.atoms[55]
        >>> shx.atoms[:3]
        [Atom ID: 53, Atom ID: 57, Atom ID: 59]
        >>> shx._reslist[55:58]
        ['', '', Atom ID: 59]
        """
        del self.shx.atoms[self.atomid]

    def to_isotropic(self) -> None:
        """
        Makes the current atom isotropic.
        """
        self.uvals = [0.04]

    '''
    def __eq__(self, other) -> bool:
        """
        Returns True if two atoms are of same name, part and residue.
        """
        if isinstance(other, str):
            return self.__repr__() == other
        if self.atomid == other.atomid:
            return True
        else:
            return False
    '''

    def find_atoms_around(self, dist=1.2, only_part=0) -> list:
        """
        Finds atoms around the current atom.
        >>> from shelxfile.shelx import ShelXFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> at = shx.atoms.get_atom_by_name('Al1')
        >>> found = at.find_atoms_around(2)
        >>> [n.atomid for n in found]
        [90, 92]
        >>> del shx.atoms[17]

        """
        found = []
        for at in self.shx.atoms:
            if atomic_distance([self.x, self.y, self.z], [at.x, at.y, at.z], self.cell) < dist:
                # Not the atom itselv:
                if not self == at:
                    # only in special part and no q-peaks:
                    if at.part == only_part and not at.qpeak:
                        found.append(at)
        return found
