from dsrmath import atomic_distance, frac_to_cart, my_isnumeric
from misc import DEBUG, ParseUnknownParam, split_fvar_and_parameter, chunks, ParseParamError, ParseNumError, \
    ParseOrderError


class Restraint():

    def __init__(self, spline: list, line_nums: list):
        """
        Base class for parsing restraints.
        TODO: resolve ranges like SADI_CCF3 O1 > F9
        Therefore, make method to get atoms of residue CCF3 and residue x<-number
        and between C21 ans C25 for atoms outside residues.
        """
        self.line_numbers = line_nums
        self.residue_class = ''
        self.residue_number = 0
        self.textline = ' '.join(spline)
        self.name = None
        self.atoms = []
        self.atoms_involved = []

    def _parse_line(self, spline, pairs=False):
        self.spline = spline
        if '_' in spline[0]:
            self.name, suffix = spline[0].upper().split('_')
            if any([x.isalpha() for x in suffix]):
                self.residue_class = suffix
            else:
                # TODO: implement _+, _- and _*
                if '*' in suffix:
                    self.residue_number = suffix
                else:
                    self.residue_number = int(suffix)
        else:
            self.name = spline[0].upper()
        # Beware! DEFS changes only the non-defined default values:
        # DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
        if DEFS.active:
            if self.name == 'DFIX':
                self.s = DEFS.sd
            if self.name == 'SAME':
                self.s1 = DEFS.sd
                self.s2 = DEFS.sd * 2
            if self.name == 'SADI':
                self.s = DEFS.sd
            if self.name == 'CHIV':
                self.s = DEFS.sf
            if self.name == 'FLAT':
                self.s = DEFS.sf
            if self.name == 'DELU':
                self.s1 = DEFS.su
                self.s2 = DEFS.su
            if self.name == 'SIMU':
                self.s = DEFS.ss
                self.st = DEFS.ss * 2
        params = []
        atoms = []
        for x in spline[1:]:
            if my_isnumeric(x):
                params.append(float(x))
            else:
                atoms.append(x)
        if pairs:
            return params, chunks(atoms, 2)
        else:
            return params, atoms

    def _paircheck(self):
        if not self.atoms:
            return
        if len(self.atoms[-1]) != 2:
            print('*** Wrong number of numerical parameters ***')
            print('Instruction: {}'.format(self.textline))
            # raise ParseNumError

    def __iter__(self):
        for x in self.textline.split():
            yield x

    def __repr__(self):
        return self.textline

    def __str__(self):
        return self.textline

    def split(self):
        return self.textline.split()


class Command():
    """
    A class to parse all general commands except restraints.
    """

    def __init__(self, spline: list, line_nums: list):
        self.line_numbers = line_nums
        self.residue_class = ''
        self.residue_number = 0
        self.textline = ' '.join(spline)
        self.name = None
        self.atoms = []

    def _parse_line(self, spline, intnums=False):
        """
        :param spline: Splitted shelxl line
        :param intnums: if numerical parameters should be integer
        :return: numerical parameters and words
        """
        if '_' in spline[0]:
            self.name, suffix = spline[0].upper().split('_')
            if any([x.isalpha() for x in suffix]):
                self.residue_class = suffix
            else:
                # TODO: implement _+, _- and _*
                self.residue_number = int(suffix)
        else:
            self.name = spline[0].upper()
        numparams = []
        words = []
        for x in spline[1:]:
            if str.isdigit(x[0]) or x[0] in '+-':
                if intnums:
                    numparams.append(int(x))
                else:
                    numparams.append(float(x))
            else:
                words.append(x)
        return numparams, words

    def split(self):
        return self.textline.split()

    def __str__(self):
        return self.textline

    def __repr__(self):
        return self.textline

class ACTA(Command):
    """
    ACTA 2θfull[#]
    >>> from shelxfile.shelx import ShelXlFile
    >>> shx = ShelXlFile('./tests/p21c.res')
    >>> shx.acta
    
    """
    def __init__(self, shx, spline: list, line_nums: list):
        super(ACTA, self).__init__(spline, line_nums)
        self.twotheta, _ = self._parse_line(spline)
        self.shx = shx

    def remove_acta_card(self):
        self.shx.delete_on_write.update([self.shx._reslist.index(self)])
        return self.textline.strip('\r\n')

    def __repr__(self):
        return "ACTA {}".format(self.twotheta if self.twotheta else '')


class FVAR():
    def __init__(self, number: int = 1, value: float = 0.0):
        """
        FVAR osf[1] free variables
        """
        self.fvar_value = value  # value
        self.number = number  # occurence inside of FVAR instructions
        self.usage = 1

    def __str__(self):
        return str(float(self.fvar_value))

    def __repr__(self):
        return str(float(self.fvar_value))


class FVARs():
    def __init__(self, shx):
        super(FVARs, self).__init__()
        self.fvars = []  # free variables
        self.shx = shx
        self._fvarline = 0

    def __iter__(self):
        """
        Must be defined for __repr__() to work.
        """
        for x in self.fvars:
            yield x

    def __getitem__(self, item: int) -> str:
        # SHELXL counts fvars from 1 to x:
        item = item - 1
        if item < 0:
            raise IndexError("*** Illegal free variable number ***")
        return self.fvars[item].fvar_value

    def __setitem__(self, key, fvar_value):
        self.fvars[key] = fvar_value

    def __len__(self) -> int:
        return len(self.fvars)

    def __str__(self) -> str:
        # returnes FVAR as list of FVAR instructions with seven numbers in one line
        lines = chunks(self.as_stringlist, 7)
        fvars = ['   '.join(i) for i in lines]
        fvars = ['FVAR   ' + i for i in fvars]
        return "\n".join(fvars)

    def __repr__(self):
        return str([x for x in self.fvars])

    @property
    def fvarline(self) -> int:
        self._fvarline = self.shx._reslist.index(self)
        return self._fvarline

    def set_free_variables(self, fvar: int, dummy_fvar: float = 0.5):
        """
        Inserts additional free variables according to the fvar number.
        """
        if fvar > 99:
            print('*** SHELXL allows only 99 free variables! ***')
            raise ParseParamError
        varlen = len(self.fvars)
        difference = (abs(fvar) - varlen)
        if difference > 0:
            for n in range(int(difference)):
                fv = FVAR(varlen + n, dummy_fvar)
                self.fvars.append(fv)

    def append(self, fvar) -> None:
        self.fvars.append(fvar)

    def set_fvar_usage(self, fvarnum: int, times: int = 1) -> None:
        """
        Incerements the usage count of a free variable by times.
        """
        fvarnum = abs(fvarnum)
        if len(self.fvars) >= abs(fvarnum):
            self.fvars[fvarnum - 1].usage += times
        elif fvarnum > 1:
            print('*** Free variable {} is not defined but used! ***'.format(fvarnum))
            # raise Exception

    def get_fvar_usage(self, fvarnum):
        """
        Returns the usage (count) of a certain free variable.
        """
        try:
            usage = self.fvars[fvarnum - 1].usage
        except IndexError:
            return 0
        return usage

    def fvars_used(self):
        """
        Retruns a dictionary with the usage of all free variables.
        """
        used = {}
        for num, fv in enumerate(self.fvars):
            used[num + 1] = fv.usage
        return used

    @property
    def as_stringlist(self):
        return [str(x.fvar_value) for x in self.fvars]


def range_resolver(atoms_range: list, atom_names: list) -> list:
    """
    Resolves the atom names of ranges like "C1 > C5"
    and works for each restraint line separately.
    :param atoms_range: atoms with a range definition
    :param atom_names: names of atoms in the fragment
    >>> r = "C2 > C5".split()
    >>> atlist = 'C1 C2 C3 C4 C5'.split()
    >>> range_resolver(r, atlist)
    ['C2', 'C3', 'C4', 'C5']
    >>> r = "C2_2 > C5_2".split()
    >>> atlist = 'C1_1 C1_2 C2_2 C3_2 C4_2 C5_2'.split()
    >>> range_resolver(r, atlist)
    ['C2_2', 'C3_2', 'C4_2', 'C5_2']
    >>> r = "C2_1 > C5_1".split()
    >>> atlist = 'C1_1 C1_2 C2_2 C3_2 C4_2 C5_2'.split()
    >>> range_resolver(r, atlist) # doctest +ELLIPSIS
    Traceback (most recent call last):
     ...
    ValueError: 'C2_1' is not in list
    """
    # dict with lists of positions of the > or < sign:
    rightleft = {'>': [], '<': []}
    for rl in rightleft:
        for num, i in enumerate(atoms_range):
            i = i.upper()
            if rl == i:
                # fill the dictionary:
                rightleft[rl].append(num)
    for rl in rightleft:
        # for each sign:
        for i in rightleft[rl]:
            # for each position of < or >:
            if rl == '>':
                # forward range
                left = atom_names.index(atoms_range[i - 1]) + 1
                right = atom_names.index(atoms_range[i + 1])
                atoms_range[i:i + 1] = atom_names[left:right]
            else:
                # backward range
                left = atom_names.index(atoms_range[i - 1])
                right = atom_names.index(atoms_range[i + 1]) + 1
                names = atom_names[right:left]
                names.reverse()  # counting backwards
                atoms_range[i:i + 1] = names
    return atoms_range


class REM(Command):
    """
    Parses REM lines
    """

    def __init__(self, spline: list, line_nums: list) -> None:
        super(REM, self).__init__(spline, line_nums)


class BOND(Command):
    """
    BOND atomnames
    """

    def __init__(self, spline: list, line_nums: list) -> None:
        super(BOND, self).__init__(spline, line_nums)
        _, self.atoms = self._parse_line(spline)


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
                del self.atomsdict[at.name+'_{}'.format(at.resinum)]
                del self.nameslist[self.nameslist.index(at.fullname.upper())]
                for x in at.line_numbers:
                    del self.shx._reslist[x]
                #self.shx.delete_on_write.update(at.line_numbers)

    @property
    def number(self) -> int:
        """
        The number of atoms in the current SHELX file.
        >>> from shelxfile.shelx import ShelXlFile
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

        >>> from shelxfile.shelx import ShelXlFile
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
        >>> from shelxfile.shelx import ShelXlFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.get_atom_by_name('Al1')
        ID: 88
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
        >>> from shelxfile.shelx import ShelXlFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.get_all_atomcoordinates() # doctest: +ELLIPSIS
        {'O1_4': [0.074835, 0.238436, 0.402457], 'C1_4': [0.028576, 0.234542, 0.337234], ...}
        """
        atdict = {}
        for at in self.atoms:
            #if at.qpeak:
            #    atdict[at.name] = at.frac_coords
            #else:
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
        >>> from shelxfile.shelx import ShelXlFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.residues
        [0, 1, 2, 3, 4]
        """
        return list({x.resinum for x in self.atoms})

    @property
    def q_peaks(self) -> list:
        r"""
        Returns a list of q-peaks in the file.
        >>> from shelxfile.shelx import ShelXlFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms.q_peaks[:5] # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        [ID: 346, ID: 347, ID: 348, ID: 349, ID: 350]
        """
        return [x for x in self.atoms if x.qpeak]

    def distance(self, atom1: str, atom2: str) -> float:
        """
        Calculates the (shortest) distance of two atoms given as text names e.g. C1_3.
        >>> from shelxfile.shelx import ShelXlFile
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
        >>> from shelxfile.shelx import ShelXlFile
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
    """
    #                name    sfac     x         y        z       occ      u11      u12 ...
    _anisatomstr = '{:<4.4s}{:>3}{:>12.6f}{:>12.6f}{:>12.6f}{:>12.5f}{:>11.5f}{:>11.5f}' \
                  ' {:>12.5f}{:>11.5f}{:>11.5f}{:>11.5f}'
    #               name    sfac     x         y         z         occ      u11
    _isoatomstr = '{:<5.5s} {:<3}{:>10.6f}  {:>10.6f}  {:>9.6f}  {:>9.5f}  {:>9.5f}'
    _fragatomstr = '{:<5.5s} {:>10.6f}  {:>10.6f}  {:>9.6f}'

    def __init__(self, shelx, spline: list, line_nums: list, line_number: int, part: int = 0,
                 afix: int = 0, residict: dict = None, sof: float = 0) -> None:
        #super(Atom, self).__init__(shelx)
        self._line_number = line_number
        self._lines = line_nums
        self.sfac_num = None
        self.name = None      # Name without residue number like "C1"
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
        self.fullname = self.name+'_{}'.format(self.resinum)
        uvals = [float(x) for x in line[6:12]]
        try:
            x, y, z = [float(x) for x in line[2:5]]
        except ValueError as e:
            if DEBUG:
                print(e, 'Line:', self.line_numbers[-1])
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

    def __repr__(self) -> str:
        return 'ID: ' + str(self.atomid)

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
                    return Atom._anisatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, *self.uvals)
                except(IndexError):
                    return 'REM Error in U values.'
            else:
                # isotropic atom
                try:
                    return Atom._isoatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, *self.uvals)
                except(IndexError):
                    return Atom._isoatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, 0.04)

    @property
    def line_numbers(self) -> list:
        return self._lines

    @line_numbers.setter
    def line_numbers(self, value: list):
        self._lines = value

    @property
    def frac_coords(self):
        return [self.x, self.y, self.z]

    @property
    def cart_coords(self):
        return [self.xc, self.yc, self.zc]

    def delete(self):
        """
        >>> from shelxfile.shelx import ShelXlFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.atoms[:3]
        [ID: 53, ID: 55, ID: 57]
        >>> shx._reslist[55:58]
        [ID: 55, '         0.01096   -0.01000    0.00201    0.00356', ID: 57]
        >>> del shx.atoms[55]
        >>> shx.atoms[:3]
        [ID: 53, ID: 57, ID: 59]
        >>> shx._reslist[55:58]
        ['         0.01096   -0.01000    0.00201    0.00356', '         0.01555   -0.00485   -0.00023    0.01102', ID: 59]
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
        >>> from shelxfile.shelx import ShelXlFile
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


class Restraints():
    """
    Base class for the list of restraints.
    """

    def __init__(self):
        self.restraints = []

    def append(self, restr):
        self.restraints.append(restr)

    def __iter__(self):
        for x in self.restraints:
            yield x

    def __getitem__(self, item):
        return self.restraints[item]

    def __repr__(self):
        if self.restraints:
            return "\n".join([str(x) for x in self.restraints])
        else:
            return 'No Restraints in file.'


class DEFS(Restraint):
    """
    DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
    Changes the *default* effective standard deviations for the following
    DFIX, SAME, SADI, CHIV, FLAT, DELU and SIMU restraints.
    """
    # keeps track if DEFS was previously activated:
    active = False
    sd = 0.02
    sf = 0.1
    su = 0.01
    ss = 0.04
    maxsof = 1

    def __init__(self, spline: list, line_nums: list):
        super(DEFS, self).__init__(spline, line_nums)
        DEFS.active = True
        p, _ = self._parse_line(spline)
        if _:
            raise ParseParamError
        if len(p) > 0:
            DEFS.sd = p[0]
        if len(p) > 1:
            DEFS.sf = p[1]
        if len(p) > 2:
            DEFS.su = p[2]
        if len(p) > 3:
            DEFS.ss = p[3]
        if len(p) > 4:
            DEFS.maxsof = p[4]

    @property
    def all(self):
        return DEFS.sd, DEFS.sf, DEFS.su, DEFS.ss, DEFS.maxsof


class NCSY(Restraint):
    """
    NCSY DN sd[0.1] su[0.05] atoms
    """

    def __init__(self, spline: list, line_nums: list):
        super(NCSY, self).__init__(spline, line_nums)
        self.sd = 0.1
        self.su = 0.05
        self.DN = None
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.DN = p[0]
        if len(p) > 1:
            self.sd = p[1]
        if len(p) > 2:
            self.su = p[2]
        if not self.DN:
            raise ParseNumError


class ISOR(Restraint):
    """
    ISOR s[0.1] st[0.2] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(ISOR, self).__init__(spline, line_nums)
        self.s = 0.1
        self.st = 0.2
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        if len(p) > 1:
            self.st = p[1]


class FLAT(Restraint):
    """
    FLAT s[0.1] four or more atoms
    """

    def __init__(self, spline: list, line_nums: list):
        super(FLAT, self).__init__(spline, line_nums)
        self.s = 0.1
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        # TODO: Have to resolve ranges first:
        # if len(self.atoms) < 4:
        #    raise ParseParamError


class BUMP(Restraint):
    """
    BUMP s [0.02]
    """

    def __init__(self, spline, line_nums):
        super(BUMP, self).__init__(spline, line_nums)
        self.s = 0.02
        p, _ = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        if _:
            raise ParseParamError


class DFIX(Restraint):
    """
    DFIX d s[0.02] atom pairs
    """

    def __init__(self, spline, line_nums):
        super(DFIX, self).__init__(spline, line_nums)
        self.s = 0.02
        p, self.atoms = self._parse_line(spline, pairs=True)
        if len(p) > 0:
            self.d = p[0]
        if len(p) > 1:
            self.s = p[1]
        self._paircheck()
        if not self.d:
            raise ParseNumError
        if 0.0001 < self.d <= self.s:  # Raise exception if d is smaller than s
            print('*** WRONG ODER of INSTRUCTIONS. d is smaller than s ***')
            print("{}".format(self.textline))


class DANG(Restraint):
    """
    DANG d s[0.04] atom pairs
    """

    def __init__(self, spline, line_nums):
        super(DANG, self).__init__(spline, line_nums)
        self.s = 0.04
        p, self.atoms = self._parse_line(spline, pairs=True)
        if len(p) > 0:
            self.d = p[0]
        if len(p) > 1:
            self.s = p[1]
        self._paircheck()
        if not self.d:
            raise ParseNumError
        if 0.0001 < self.d <= self.s:  # Raise exception if d is smaller than s
            raise ParseOrderError


class SADI(Restraint):
    """
    SADI s[0.02] pairs of atoms
    """

    def __init__(self, spline, line_nums):
        super(SADI, self).__init__(spline, line_nums)
        self.s = 0.02
        p, self.atoms = self._parse_line(spline, pairs=True)
        if len(p) > 0:
            self.s = p[0]
        self._paircheck()


class SAME(Restraint):
    """
    SAME s1[0.02] s2[0.04] atomnames
    """

    def __init__(self, spline, line_nums):
        super(SAME, self).__init__(spline, line_nums)
        self.s1 = 0.02
        self.s2 = 0.04
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s1 = p[0]
        if len(p) > 1:
            self.s2 = p[1]


class RIGU(Restraint):
    """
    RIGU s1[0.004] s2[0.004] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(RIGU, self).__init__(spline, line_nums)
        self.s1 = 0.004
        self.s2 = 0.004
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s1 = p[0]
        if len(p) > 1:
            self.s2 = p[1]


class SIMU(Restraint):
    """
    SIMU s[0.04] st[0.08] dmax[2.0] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(SIMU, self).__init__(spline, line_nums)
        self.s = 0.04
        self.st = 0.08
        self.dmax = 2.0
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        if len(p) > 1:
            self.st = p[1]
        if len(p) > 2:
            self.dmax = p[2]


class DELU(Restraint):
    """
    DELU s1[0.01] s2[0.01] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(DELU, self).__init__(spline, line_nums)
        self.s1 = 0.01
        self.s2 = 0.01
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s1 = p[0]
        if len(p) > 1:
            self.s2 = p[1]


class CHIV(Restraint):
    """
    CHIV V[0] s[0.1] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(CHIV, self).__init__(spline, line_nums)
        self.s = 0.1
        self.V = 0.0
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.V = p[0]
        if len(p) > 1:
            self.s = p[1]


class EADP(Restraint):
    """
    EADP atomnames
    """

    def __init__(self, spline: list, line_nums: list) -> None:
        super(EADP, self).__init__(spline, line_nums)
        _, self.atoms = self._parse_line(spline, pairs=False)


class EXYZ(Restraint):
    """
    EADP atomnames
    """

    def __init__(self, spline: list, line_nums: list) -> None:
        super(EXYZ, self).__init__(spline, line_nums)
        _, self.atoms = self._parse_line(spline, pairs=False)


class DAMP(Command):
    """
    DAMP damp[0.7] limse[15]
    """
    def __init__(self, spline, line_nums):
        super(DAMP, self).__init__(spline, line_nums)
        values, _ = self._parse_line(spline, intnums=False)
        self.damp, self.limse = 0, 0
        if len(values) > 0:
            self.damp = values[0]
        if len(values) > 1:
            self.damp, self.limse = values

    def __repr__(self) -> str:
        if self.limse == 0:
            return "DAMP  {:,g}".format(self.damp)
        else:
            return "DAMP  {:,g} {:,g}".format(self.damp, self.limse)


class HFIX(Command):
    """
    HFIX mn U[#] d[#] atomnames
    """
    def __init__(self, spline: list, line_nums: list):
        super(HFIX, self).__init__(spline, line_nums)
        self.params, self.atoms = self._parse_line(spline, intnums=True)

    def __repr__(self):
        return "HFIX {} {}".format(" ".join([str(x) for x in self.params]) if self.params else '',
                                   " ".join(self.atoms) if self.atoms else '')


class HKLF(Command):
    """
    HKLF N[0] S[1] r11...r33[1 0 0 0 1 0 0 0 1] sm[1] m[0]
    """
    def __init__(self, spline: list, line_nums: list):
        super(HKLF, self).__init__(spline, line_nums)
        p, _ = self._parse_line(spline)
        self.n = 0
        self.s = 1
        self.matrix = [1, 0, 0, 0, 1, 0, 0, 0, 1]
        self.sm = 1
        self.m = 0
        if len(p) > 0:
            self.n = p[0]
        if len(p) > 1:
            self.s = p[1]
        if len(p) > 10:
            self.matrix = p[3:11]
        if len(p) > 11:
            self.sm = p[12]
        if len(p) > 12:
            self.m = p[13]

    def __repr__(self):
        return "HKLF {:,g} {:,g}  {}  {:,g} {:,g}".format(self.n, self.s, ' '.join([str(i) for i in self.matrix]), self.sm, self.m)


class SUMP(Command):
    """
    SUMP for linear equation eypressions with free variables.
    SUMP c sigma c1 m1 c2 m2 ...
    """

    def __init__(self, spline, line_nums):
        super(SUMP, self).__init__(spline, line_nums)
        p, _ = self._parse_line(spline)
        self.c = p.pop(0)
        self.fvars = {}
        self.sigma = p.pop(0)
        # this is to have integer free variables
        fvars = [int(x) for x in p[1::2]]
        times = [x for x in p[0::2]]
        self.fvars = [[x, y] for x, y in zip(times, fvars)]

    def __getitem__(self, item):
        return self.fvars[item]


class SYMM(Command):
    """
    Container for a symm card.
    """

    def __init__(self, spline, line_nums):
        super(SYMM, self).__init__(spline, line_nums)
        self.symmcards = self._parse_line(spline)

    def _parse_line(self, spline, intnums=False):
        symmcards = []
        line = ''.join(spline[1:])  # removes whitespace
        symmcards.append(line.split(','))
        return symmcards

    def __repr__(self):
        return "\n".join(["SYMM  " + "  ".join(x) for x in self.symmcards])

    def __str__(self):
        return "\n".join(["SYMM  " + "  ".join(x) for x in self.symmcards])


class LSCycles():
    def __init__(self, shx, spline: list, line_number: int = 0):
        """
        L.S. nls[0] nrf[0] nextra[0]
        If nrf is positive, it is the number of these cycles that should be performed before applying ANIS.
        Negative nrf indicates which reflections should be ignored during the refinement but used instead for
        the calculation of free R-factors in the final structure factor summation.
        nextra is the number of additional parameters that were derived from the data when 'squeezing' the
        structure etc.
        """
        self.shx = shx
        self.cgls = False
        self.cycles = 0
        self.nrf = ''
        self.nextra = ''
        self.line_number = line_number  # line number in res file
        try:
            self.cycles = int(spline[1])
        except (IndexError, NameError):
            raise ParseNumError
        try:
            self.nrf = spline[2]
        except IndexError:
            pass
        try:
            self.nextra = spline[3]
        except IndexError:
            pass
        if spline[0].upper() == 'CGLS':
            self.cgls = True

    @property
    def number(self):
        return self.cycles

    def set_refine_cycles(self, number: int):
        """
        Sets the number of refinement cycles for the current res file.
        >>> from shelxfile.shelx import ShelXlFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.cycles.set_refine_cycles(44)
        >>> shx._reslist[shx.cycles.line_number]
        L.S. 44
        """
        self.cycles = number
        #self.shx.reslist[self.shx.reslist.index(self)] = self.text.strip('\r\n')

    @property
    def text(self):
        """
        'CGLS 10 2 '
        """
        return self.__repr__()

    def __repr__(self):
        return '{} {} {} {}'.format('CGLS' if self.cgls else 'L.S.', self.cycles,
                                    self.nrf if self.nrf else '', self.nextra if self.nextra else '').strip()


class SFACTable():
    def __init__(self, shx):
        """
        Holds the information of SFAC instructions. Either with default values and only elements
        SFAC elements
        or as explicit scattering factor in the form of an exponential series, followed by real and
        imaginary dispersion terms, linear absorption coefficient, covalent radius and atomic weight.

        SFAC elements  or  SFAC E a1 b1 a2 b2 a3 b3 a4 b4 c f' f" mu r wt
        """
        self.sfac_table = []
        self.shx = shx
        self.elements_list = []

    def __iter__(self):
        for x in self.sfac_table:
            yield x['element'].capitalize()

    def __repr__(self):
        regular = []
        exponential = []
        for sf in self.sfac_table:
            if not self.is_exp(sf):
                regular.append(sf['element'].capitalize())
            else:
                values = []
                for x in ['element', 'a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'a4', 'b4', 'c',
                          'fprime', 'fdprime', 'mu', 'r', 'wt']:
                    values.append(sf[x])
                exponential.append("SFAC " + "  ".join(values))
        sfac = "SFAC {:<4s}".format("  ".join(regular))
        if exponential:
            return "\n".join(exponential) + "\n" + sfac
        else:
            return sfac

    def __getitem__(self, index: int):
        """
        Returns the n-th element in the sfac table, beginning with 1.
        """
        if index == 0:
            raise IndexError
        if index < 0:
            index = len(self.sfac_table) + index + 1
        return self.sfac_table[index - 1]['element'].capitalize()


    def parse_element_line(self, spline: list):
        """
        Adds a new SFAC card to the list of cards.
        >>> from shelxfile.shelx import ShelXlFile
        >>> shx = ShelXlFile('./tests/p21c.res')
        >>> shx.sfac_table
        SFAC C  H  O  F  Al  Ga
        >>> shx.unit
        UNIT 1  2  3  4  5  6
        >>> shx.sfac_table.add_element('Au')
        >>> shx.sfac_table
        SFAC C  H  O  F  Al  Ga  Au
        >>> shx.unit
        UNIT 1  2  3  4  5  6  1
        """
        if not ''.join(spline[1:]).isalpha():  # joining, because space is not alphabetical
            # Excplicit with all values
            sfdic = {}
            for n, x in enumerate(['element', 'a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'a4', 'b4', 'c',
                                   'fprime', 'fdprime', 'mu', 'r', 'wt']):
                if n == 0:
                    self.elements_list.append(spline[n + 1].upper())
                try:
                    sfdic[x] = spline[n + 1]
                except IndexError:
                    raise ParseNumError()
            self.sfac_table.append(sfdic)
        else:
            # Just the elements
            for x in spline[1:]:
                self.elements_list.append(x.upper())
                self.sfac_table.append({'element': x.upper()})

    def has_element(self, element):
        return element.upper() in self.elements_list

    @staticmethod
    def is_exp(item):
        return 'a1' in item

    def add_element(self, element: str):
        """
        Adds an element to the SFAC list. 
        """
        if self.has_element(element):
            return
        self.elements_list.append(element.upper())
        self.sfac_table.append({'element': element.upper(), 'line_number': None})
        self.shx.unit.add_number(1)

    def remove_element(self, element: str):
        del self.sfac_table[self.shx.elem2sfac(element.upper()) - 1]
        del self.elements_list[self.elements_list.index(element.upper())]


class UNIT(Command):
    """
    UNIT n1 n2 ...
    """
    def __init__(self, spline: list, line_nums: list):
        super(UNIT, self).__init__(spline, line_nums)
        self.values, _ = self._parse_line(spline)

    def add_number(self, number: float):
        self.values.append(number)

    def __iter__(self):
        yield [x for x in self.values]

    def __repr__(self):
        return "UNIT " + "  ".join(["{:,g}".format(x) for x in self.values])

    def __setitem__(self, key, value):
        self.values[key] = value

    def __getitem__(self, item):
        return self.values[item]

    def __add__(self, other):
        self.values.append(other)


class BASF(Command):
    """
    BASF scale factors
    BASF can occour in multiple lines.
    """

    def __init__(self, spline: list, line_numbers: list):
        super(BASF, self).__init__(spline, line_numbers)
        self.scale_factors, _ = self._parse_line(spline)
        del self.atoms

    def __iter__(self):
        yield self.scale_factors


class TWIN(Command):
    """
    TWIN 3x3 matrix [-1 0 0 0 -1 0 0 0 -1] N[2]
    +N     -N  m = |N|
    m-1 to 2m-1
    m-1   (2*abs(m)/2)-1
    """

    def __init__(self, spline: list, line_nums: list):
        super(TWIN, self).__init__(spline, line_nums)
        self.matrix = [-1, 0, 0, 0, -1, 0, 0, 0, -1]
        self.allowed_N = 2
        self.n_value = 2
        if len(spline) > 1:
            p, _ = self._parse_line(spline, intnums=False)
            if len(p) == 9:
                self.matrix = p
            elif len(p) == 10:
                self.matrix = p[:9]
                self.n_value = int(p[9])
            else:
                raise ParseNumError("*** Check TWIN instruction. ***")
        m = abs(self.n_value) / 2
        if self.n_value > 0:
            self.allowed_N = abs(self.n_value) - 1
        else:
            self.allowed_N = (2 * m) - 1


class WGHT(Command):
    """
    The weighting scheme is defined as follows:
    w = q / [ σ²(Fo²) + (a*P)² + b*P + d + e*sin(θ)/$lambda; ]

    WGHT a[0.1] b[0] c[0] d[0] e[0] f[.33333]
    Usually only WGHT a b
    """

    def __init__(self, spline: list, line_nums: list):
        super(WGHT, self).__init__(spline, line_nums)
        self.a = 0.1
        self.b = 0.0
        self.c = 0.0
        self.d = 0.0
        self.e = 0.0
        self.f = 0.33333
        self.line_numbers = line_nums[0]
        p, _ = self._parse_line(spline)
        if len(p) > 0:
            self.a = p[0]
        if len(p) > 1:
            self.b = p[1]
        if len(p) > 2:
            self.c = p[2]
        if len(p) > 3:
            self.d = p[3]
        if len(p) > 4:
            self.e = p[4]
        if len(p) > 5:
            self.f = p[5]

    def __repr__(self):
        wght = 'WGHT {} {}'.format(self.a, self.b)
        # It is very unlikely that someone changes other parameter than a and b:
        if (self.c + self.d + self.e + self.f) != 0.33333:
            wght += ' {} {} {} {}'.format(self.c, self.d, self.e, self.f)
        return wght