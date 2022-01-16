import re
from math import cos, radians, sqrt
from typing import List

from shelxfile.misc.dsrmath import my_isnumeric, SymmetryElement
from shelxfile.misc.misc import chunks, ParseParamError, ParseNumError, \
    ParseOrderError, DEBUG, ParseSyntaxError

"""
SHELXL cards:

ABIN n1 n2
ACTA 2θfull[#]
AFIX mn d[#] sof[11] U[10.08]
ANIS n
ANIS names
ANSC six coefficients
ANSR anres[0.001]
BASF scale factors
BIND atom1 atom2
BIND m n
BLOC n1 n2 atomnames
BOND atomnames
BUMP s [0.02]
CELL λ a b c α β γ
CGLS nls[0] nrf[0] nextra[0]
CHIV V[0] s[0.1] atomnames
CONF atomnames max_d[1.9] max_a[170]
CONN bmax[12] r[#] atomnames or CONN bmax[12]
DAMP damp[0.7] limse[15]
DANG d s[0.04] atom pairs
DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
DELU s1[0.01] s2[0.01] atomnames
DFIX d s[0.02] atom pairs
DISP E f' f"[#] mu[#]
EADP atomnames
END
EQIV $n symmetry operation
EXTI x[0]
EXYZ atomnames
FEND
FLAT s[0.1] four or more atoms
FMAP code[2] axis[#] nl[53]
FRAG code[17] a[1] b[1] c[1] α[90] β[90] γ[90]
FREE atom1 atom2
FVAR osf[1] free variables
GRID sl[#] sa[#] sd[#] dl[#] da[#] dd[#]
HFIX mn U[#] d[#] atomnames
HKLF N[0] S[1] r11...r33[1 0 0 0 1 0 0 0 1] sm[1] m[0]
HTAB dh[2.0]
HTAB donor-atom acceptor-atom
ISOR s[0.1] st[0.2] atomnames
LATT N[1]
LAUE E
LIST m[#] mult[1]
L.S. nls[0] nrf[0] nextra[0]
MERG n[2]
MORE m[1]
MOVE dx[0] dy[0] dz[0] sign[1]
MPLA na atomnames
NCSY DN sd[0.1] su[0.05] atoms
NEUT
OMIT atomnames
OMIT s[-2] 2θ(lim)[180]
OMIT h k l
PART n sof
PLAN npeaks[20] d1[#] d2[#]
PRIG p[#]
REM
RESI class[ ] number[0] alias
RIGU s1[0.004] s2[0.004] atomnames
RTAB codename atomnames
SADI s[0.02] pairs of atoms
SAME s1[0.02] s2[0.04] atomnames
SFAC elements
SFAC E a1 b1 a2 b2 a3 b3 a4 b4 c f' f" mu r wt
SHEL lowres[infinite] highres[0]
SIMU s[0.04] st[0.08] dmax[2.0] atomnames
SIZE dx dy dz
SPEC del[0.2]
STIR sres step[0.01]
SUMP c sigma c1 m1 c2 m2 ...
SWAT g[0] U[2]
SYMM symmetry operation
TEMP T[20]
TITL [ ]
TWIN 3x3 matrix [-1 0 0 0 -1 0 0 0 -1] N[2]
TWST N[0]
UNIT n1 n2 ...
WGHT a[0.1] b[0] c[0] d[0] e[0] f[.33333]
WIGL del[0.2] dU[0.2]
WPDB n[1]
XNPD Umin[-0.001]
ZERR Z esd(a) esd(b) esd(c) esd(α) esd(β) esd(γ)
"""


class Residue():
    def __init__(self):
        self.shx = None
        self._spline = None

    @property
    def residue_number(self):
        if '_' in self._spline[0]:
            _, suffix = self._spline[0].upper().split('_')
            if any([x.isalpha() for x in suffix]):
                self.residue_class = suffix
            else:
                # TODO: implement _+, _- and _*
                if '*' in suffix:
                    return list(self.shx.residues.residue_numbers.keys())
                else:
                    return [int(suffix)]
        return [0]


class Restraint(Residue):

    def __init__(self, shx, spline: list):
        """
        Base class for parsing restraints.
        TODO: resolve ranges like SADI_CCF3 O1 > F9
        """
        super().__init__()
        self.shx = shx
        self.residue_class = ''  # '' is the default class (with residue number 0)
        self.textline = ' '.join(spline)
        self.name = None
        self.atoms = []

    @property
    def index(self):
        return self.shx.index_of(self)

    def _parse_line(self, spline, pairs=False):
        """
        Residues may be referenced by any instruction that allows atom names; the reference takes
        the form of the character '_' followed by either the residue class or number without intervening
        spaces.
        Individual atom names in an instruction may be followed by '_' and a residue number, but not by '_* ' or '_'
        and a residue class. If an atom name is not followed by a residue number, the current residue is
        assumed (unless overridden by a global residue number or class appended to the instruction
        codeword).
        """
        self._spline = spline
        # TODO: check if suffix is a residue class and set instance variable accordingly
        if '_' in spline[0]:
            self.name, suffix = spline[0].upper().split('_')
        else:
            self.name = spline[0].upper()
        self.set_defs_values()
        self.check_if_class_name_fits_to_command()
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

    def check_if_class_name_fits_to_command(self):
        if DEBUG and self.__class__.__name__ != self.name:
            print('*** Trying to parse restraint with wrong class ***')
            raise ParseSyntaxError

    def set_defs_values(self):
        # Beware! DEFS changes only the non-defined default values:
        # DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
        if self.shx.defs:  # and self.shx.defs.active:
            if self.name == 'DFIX':
                self.s = self.shx.defs.sd
            if self.name == 'SAME':
                self.s1 = self.shx.defs.sd
                self.s2 = self.shx.defs.sd * 2
            if self.name == 'SADI':
                self.s = self.shx.defs.sd
            if self.name == 'CHIV':
                self.s = self.shx.defs.sf
            if self.name == 'FLAT':
                self.s = self.shx.defs.sf
            if self.name == 'DELU':
                self.s1 = self.shx.defs.su
                self.s2 = self.shx.defs.su
            if self.name == 'SIMU':
                self.s = self.shx.defs.ss
                self.st = self.shx.defs.ss * 2

    def _paircheck(self):
        if not self.atoms:
            return
        if len(self.atoms[-1]) != 2 and DEBUG:
            print('*** Wrong number of numerical parameters ***')
            print('Instruction: {}'.format(self.textline))
            raise ParseNumError

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

    def __init__(self, shx, spline: list):
        self._shx = shx
        self._spline = spline
        self.residue_class = ''
        self._textline = ' '.join(spline)

    def _parse_line(self, spline, intnums=False):
        """
        :param spline: Splitted shelxl line
        :param intnums: if numerical parameters should be integer
        :return: numerical parameters and words
        """
        if '_' in spline[0]:
            self._card_name, suffix = spline[0].upper().split('_')
        else:
            self._card_name = spline[0].upper()
        numparams = []
        words = []
        for x in spline[1:]:  # all values after SHELX card
            if str.isdigit(x[0]) or x[0] in '+-':
                if intnums:
                    numparams.append(int(x))
                else:
                    numparams.append(float(x))
            else:
                words.append(x)
        return numparams, words

    @property
    def index(self):
        return self._shx.index_of(self)

    @property
    def position(self) -> int:
        return self.index

    def __iter__(self):
        for x in self.__repr__().split():
            yield x

    def split(self):
        return self._textline.split()

    def __str__(self):
        return self._textline

    def __repr__(self):
        return self._textline


class ABIN(Command):

    def __init__(self, shx, spline):
        """
        ABIN n1 n2
        """
        super(ABIN, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        if len(p) > 0:
            self.n1 = p[0]
        if len(p) > 1:
            self.n2 = p[1]


class ANIS(Command):

    def __init__(self, shx, spline: List):
        """
        ANIS n
        ANIS names
        """
        super(ANIS, self).__init__(shx, spline)
        p, self.atoms = self._parse_line(spline)
        self.over_all = True
        if len(p) > 0:
            self.over_all = False
            self.n = p[0]
        if len(self.atoms) > 0:
            self.over_all = False

    def __bool__(self):
        return True


class MPLA(Command):

    def __init__(self, shx, spline: list):
        """
        MPLA na atomnames
        """
        super(MPLA, self).__init__(shx, spline)
        p, self.atoms = self._parse_line(spline, intnums=True)
        if len(p) > 0:
            self.na = p[0]


class MORE(Command):

    def __init__(self, shx, spline: list):
        """
        MORE m[1]
        """
        super(MORE, self).__init__(shx, spline)
        self.m = 1
        p, _ = self._parse_line(spline, intnums=True)
        self.m = p[0]


class CELL(Command):

    def __init__(self, shx, spline: list):
        """
        CELL λ a b c α β γ
        """
        super(CELL, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self._cell_list = []
        if len(p) > 0:
            self.wavelen = float(p[0])
        if len(p) > 6:
            self._cell_list = p[1:]
            self.a = p[1]
            self.b = p[2]
            self.c = p[3]
            self.alpha = p[4]
            self.beta = p[5]
            self.gamma = p[6]
            self.cosal = cos(radians(self.alpha))
            self.cosbe = cos(radians(self.beta))
            self.cosga = cos(radians(self.gamma))
        else:
            raise ParseSyntaxError

    @property
    def volume(self) -> float:
        """
        calculates the volume of a unit cell
        """
        try:
            ca, cb, cg = cos(radians(self.alpha)), cos(radians(self.beta)), cos(radians(self.gamma))
            v = self.a * self.b * self.c * sqrt(1 + 2 * ca * cb * cg - ca ** 2 - cb ** 2 - cg ** 2)
        except AttributeError:
            # No valid celll
            v = 0.0
        return v

    @property
    def V(self) -> float:
        return self.volume

    def __iter__(self):
        return iter(self._cell_list)

    def __getitem__(self, item):
        return self._cell_list[item]


class ZERR(Command):

    def __init__(self, shx, spline: list):
        """
        ZERR Z esd(a) esd(b) esd(c) esd(α) esd(β) esd(γ)
        """
        super(ZERR, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.Z = 1
        if len(p) > 0:
            self.Z = p[0]
        if len(p) > 6:
            self.esd_list = p[1:]
            self.esd_a = p[0]
            self.esd_b = p[1]
            self.esd_c = p[2]
            self.esd_al = p[3]
            self.esd_be = p[4]
            self.esd_ga = p[5]


class AFIX(Command):

    def __init__(self, shx, spline: list):
        """
        AFIX mn d[#] sof[11] U[10.08]
        """
        super(AFIX, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.U = 10.08
        self.sof = 11.0
        if len(p) > 0:
            self.mn = int(p[0])
        if len(p) > 1:
            self.d = p[1]
        if len(p) > 2:
            self.sof = p[2]
        if len(p) > 3:
            self.U = p[3]

    def __bool__(self):
        if self.mn > 0:
            return True
        else:
            return False


class Residues():

    def __init__(self, shx):
        self.shx = shx
        self.all_residues = []
        self.residue_classes = {}  # class: numbers
        # self.residue_numbers = {}  # number: class

    def append(self, resi: 'RESI') -> None:
        """
        Adds a new residues to the list of residues.
        """
        self.all_residues.append(resi)
        # Collect dict with class: numbers
        if resi.residue_class in self.residue_classes:
            self.residue_classes[resi.residue_class].append(resi.residue_number)
        else:
            self.residue_classes[resi.residue_class] = [resi.residue_number]
        """
        # Collect dict with number: classes
        if resi.residue_number in self.residue_numbers:
            if DEBUG:
                print('*** Duplicate residue number {} found! ***'.format(resi.residue_number))
        else:
            self.residue_numbers[resi.residue_number] = resi.residue_class
        """

    @property
    def residue_numbers(self):
        return dict((x.residue_number, x.residue_class) for x in self.shx.residues.all_residues)


class RESI():

    def __init__(self, shx, spline: list):
        """
        RESI class[ ] number[0] alias
        """
        self.shx = shx
        self.residue_class = ''
        self.residue_number = 0
        self.alias = None
        self.chain_id = None
        self.textline = ' '.join(spline)
        if len(spline) < 2 and DEBUG:
            print('*** Wrong RESI definition found! Check your RESI instructions ***')
            raise ParseParamError
        self.get_resi_definition(spline)
        if self.residue_number < -999 or self.residue_number > 9999:
            print('*** Invalid residue number given. ****')
            raise ParseSyntaxError

    def get_resi_definition(self, resi: list) -> tuple:
        """
        RESI class[ ] number[0] alias

        Returns the residue number and class of a string like 'RESI TOL 1'
        or 'RESI 1 TOL'

        Residue names may now begin with a digit.
        They must however contain at least one letter

        Allowed residue numbers is now from -999 to 9999 (2017/1)
        """
        alpha = re.compile('[a-zA-Z]')
        for x in resi:
            if alpha.search(x):
                if ':' in x:
                    # contains ":" thus must be a chain-id+number
                    self.chain_id, self.residue_number = x.split(':')[0], int(x.split(':')[1])
                else:
                    # contains letters, must be a name (class)
                    self.residue_class = x
            else:
                # everything else can only be a number
                if self.residue_number > 0:
                    self.alias = int(x)
                else:
                    try:
                        self.residue_number = int(x)
                    except ValueError:
                        self.residue_number = 0
        return self.residue_class, self.residue_number, self.chain_id, self.alias

    def _parse_line(self, spline, intnums=False):
        """
        :param spline: Splitted shelxl line
        :param intnums: if numerical parameters should be integer
        :return: numerical parameters and words
        """
        if '_' in spline[0]:
            self.card_name, suffix = spline[0].upper().split('_')
            if any([x.isalpha() for x in suffix]):
                self.residue_class = suffix
            else:
                # TODO: implement _+, _- and _*
                self.residue_number = int(suffix)
        else:
            self.card_name = spline[0].upper()
        numparams = []
        words = []
        for x in spline[1:]:  # all values after SHELX card
            if str.isdigit(x[0]) or x[0] in '+-':
                if intnums:
                    numparams.append(int(x))
                else:
                    numparams.append(float(x))
            else:
                words.append(x)
        return numparams, words

    @property
    def index(self):
        return self.shx.index_of(self)

    def __iter__(self):
        for x in self.__repr__().split():
            yield x

    def split(self):
        return self.textline.split()

    def __str__(self):
        return self.textline

    def __repr__(self):
        return self.textline

    def __bool__(self):
        if self.residue_number > 0:
            return True
        else:
            return False


class PART(Command):

    def __init__(self, shx, spline: list):
        """
        PART n sof
        """
        super(PART, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.sof = 11.0
        self.n = 0
        try:
            self.n = int(p[0])
        except(ValueError, IndexError):
            if DEBUG:
                print('*** Wrong PART definition in line {} found! '
                      'Check your PART instructions ***'.format(shx.error_line_num))
                raise
            self.n = 0
        if len(p) > 1:
            self.sof = float(p[1])

    def __bool__(self):
        if self.n > 0:
            return True
        else:
            return False


class XNPD(Command):

    def __init__(self, shx, spline: list):
        """
        XNPD Umin[-0.001]
        """
        super(XNPD, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.Umin = -0.001
        if len(p) > 0:
            self.Umin = p[0]


class SIZE(Command):

    def __init__(self, shx, spline: list):
        """
        SIZE dx dy dz
        """
        super(SIZE, self).__init__(shx, spline)
        self.dx, self.dy, self.dz = None, None, None
        p, _ = self._parse_line(spline)
        if len(p) > 0:
            self.dx = p[0]
        if len(p) > 1:
            self.dy = p[1]
        if len(p) > 2:
            self.dz = p[2]

    def _as_text(self):
        if all([self.dx, self.dy, self.dz]):
            return "SIZE {:,g} {:,g} {:,g}".format(self.dx, self.dy, self.dz)
        else:
            return ""

    def __repr__(self):
        return self._as_text()

    def __str__(self):
        return self._as_text()


class SHEL(Command):

    def __init__(self, shx, spline: list):
        """
        SHEL lowres[infinite] highres[0]
        """
        super(SHEL, self).__init__(shx, spline)
        params, _ = self._parse_line(spline)
        if len(params) > 0:
            self.lowres = params[0]
        if len(params) > 1:
            self.highres = params[1]


class WIGL(Command):

    def __init__(self, shx, spline: list):
        """
        WIGL del[0.2] dU[0.2]
        """
        super(WIGL, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.d = 0.2
        self.dU = 0.2
        if len(p) > 0:
            self.d = p[0]
        if len(p) > 1:
            self.dU = p[1]


class WPDB(Command):

    def __init__(self, shx, spline: list):
        """
        WPDB n[1]
        """
        super(WPDB, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.n = 1
        if len(p) > 0:
            self.n = p[0]


class SPEC(Command):

    def __init__(self, shx, spline: list):
        """
        SPEC d[0.2]
        """
        super(SPEC, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        if len(p) > 0:
            self.d = p[0]


class STIR(Command):

    def __init__(self, shx, spline: list):
        """
        STIR sres step[0.01]
        """
        super(STIR, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.step = 0.01
        if len(p) > 0:
            self.sres = p[0]
        if len(p) > 1:
            self.step = p[1]

    def __repr__(self):
        return "STIR {} {}".format(self.sres if self.sres else '', self.step)

    def __str__(self):
        return "STIR {} {}".format(self.sres if self.sres else '', self.step)


class TWST(Command):

    def __init__(self, shx, spline: list):
        """
        TWST N[0] (N[1] after SHELXL-2018/3)
        Twin component number to be used for the completeness and Friedel completeness statistics.
        """
        super(TWST, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.N = 1
        if len(p) > 0:
            self.N = p[0]


class RTAB(Command):

    def __init__(self, shx, spline: list):
        """
        RTAB codename atomnames
        """
        super(RTAB, self).__init__(shx, spline)
        self.code = spline.pop(1)
        _, self.atoms = self._parse_line(spline)


class PRIG(Command):

    def __init__(self, shx, spline: list):
        """
        PRIG p[#]
        """
        super(PRIG, self).__init__(shx, spline)
        params, _ = self._parse_line(spline)
        if len(params) > 0:
            self.p = params[0]


class PLAN(Command):

    def __init__(self, shx, spline: list):
        """
        PLAN npeaks[20] d1[#] d2[#]
        """
        super(PLAN, self).__init__(shx, spline)
        params, _ = self._parse_line(spline)
        if len(params) > 0:
            self.npeaks = params[0]
        if len(params) > 1:
            self.d1 = params[1]
        if len(params) > 2:
            self.d2 = params[2]


class FRAG(Command):

    def __init__(self, shx, spline: list):
        """
        FRAG code[17] a[1] b[1] c[1] α[90] β[90] γ[90]
        """
        super(FRAG, self).__init__(shx, spline)
        params, _ = self._parse_line(spline)
        self.code = params[0]
        self.cell = params[1:7]


class FREE(Command):

    def __init__(self, shx, spline: list):
        """
        FREE atom1 atom2
        """
        super(FREE, self).__init__(shx, spline)
        _, atoms = self._parse_line(spline)
        try:
            self.atom1 = atoms[0]
            self.atom2 = atoms[1]
        except IndexError:
            raise ParseParamError


class FMAP(Command):
    """
    FMAP code[2] axis[#] nl[53]
    """

    def __init__(self, shx, spline: list):
        super(FMAP, self).__init__(shx, spline)
        params, _ = self._parse_line(spline)
        if len(params) > 0:
            self.code = params[0]
        if len(params) > 1:
            self.axis = params[1]
        if len(params) > 2:
            self.axis = params[2]


class MOVE(Command):

    def __init__(self, shx, spline: list):
        """
        MOVE dx[0] dy[0] dz[0] sign[1]
        """
        super(MOVE, self).__init__(shx, spline)
        params, _ = self._parse_line(spline)
        if len(params) > 2:
            self.dxdydz = params[:3]
        if len(params) > 3:
            self.sign = params[3]


class MERG(Command):

    def __init__(self, shx, spline: list):
        """
        MERG n[2]
        """
        super(MERG, self).__init__(shx, spline)
        self.n = None
        _n, _ = self._parse_line(spline)
        if len(_n) > 0:
            self.n = _n[0]


class HTAB(Command):

    def __init__(self, shx, spline: list):
        """
        HTAB dh[2.0]
        HTAB donor-atom acceptor-atom
        """
        super(HTAB, self).__init__(shx, spline)
        self.dh = None
        dh, atoms = self._parse_line(spline)
        if dh:
            self.dh = dh[0]
        if len(atoms) == 2:
            self.donor = atoms[0]
            self.acceptor = atoms[1]


class GRID(Command):

    def __init__(self, shx, spline: list):
        """
        GRID sl[#] sa[#] sd[#] dl[#] da[#] dd[#]
        """
        super(GRID, self).__init__(shx, spline)
        params, _ = self._parse_line(spline)
        if len(params) > 0:
            self.sl = params[0]
        if len(params) > 1:
            self.sa = params[1]
        if len(params) > 2:
            self.sd = params[2]
        if len(params) > 3:
            self.dl = params[3]
        if len(params) > 4:
            self.da = params[4]
        if len(params) > 5:
            self.dd = params[5]


class ACTA(Command):

    def __init__(self, shx, spline: list):
        """
        ACTA 2θfull[#]
        """
        super(ACTA, self).__init__(shx, spline)
        self.twotheta, _ = self._parse_line(spline)
        self.shx = shx

    def _as_str(self):
        if self.twotheta:
            return "ACTA {:,g}".format(self.twotheta[0])
        else:
            return "ACTA"

    def __repr__(self):
        return self._as_str()

    def __str__(self):
        return self._as_str()


class BLOC(Command):

    def __init__(self, shx, spline: list):
        """
        BLOC n1 n2 atomnames
        """
        super(BLOC, self).__init__(shx, spline)
        params, self.atoms = self._parse_line(spline)
        if len(params) > 1:
            self.n2 = params[1]
        if len(params) > 0:
            self.n1 = params[0]
        self.shx = shx


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
        item = abs(item) - 1
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
    def position(self) -> int:
        return self.shx.index_of(self)

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
        elif fvarnum > 1 and DEBUG:
            print('*** Free variable {} is not defined but used! ***'.format(fvarnum))

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


class CONF(Command):

    def __init__(self, shx, spline: list) -> None:
        """
        CONF atomnames max_d[1.9] max_a[170]
        """
        super(CONF, self).__init__(shx, spline)


class CONN(Command):

    def __init__(self, shx, spline: list) -> None:
        """
        CONN bmax[12] r[#] atomnames or CONN bmax[12]
        """
        super(CONN, self).__init__(shx, spline)


class REM(Command):

    def __init__(self, shx, spline: list) -> None:
        """
        Parses REM lines
        """
        super(REM, self).__init__(shx, spline)


class BIND(Command):

    def __init__(self, shx, spline: list) -> None:
        """
        BIND atom1 atom2
        BIND m n
        """
        super(BIND, self).__init__(shx, spline)
        self.parts, self.atoms = self._parse_line(spline)


class BOND(Command):

    def __init__(self, shx, spline: list) -> None:
        """
        BOND atomnames
        """
        super(BOND, self).__init__(shx, spline)
        _, self.atoms = self._parse_line(spline)


class DISP(Command):
    """
    DISP E f' f"[#] mu[#]
    """

    def __init__(self, shx, spline: list) -> None:
        super(DISP, self).__init__(shx, spline)
        self.element, self.parameter = self._parse_line(spline)


class Restraints():
    """
    Base class for the list of restraints.
    :type _restraints: List[Restraint]
    """

    def __init__(self):
        """
        """
        self._restraints = []

    def append(self, restr):
        self._restraints.append(restr)

    def __iter__(self):
        for x in self._restraints:
            yield x

    def __getitem__(self, item):
        return self._restraints[item]

    def __repr__(self):
        if self._restraints:
            return "\n".join([str(x) for x in self._restraints])
        else:
            return 'No Restraints in file.'

    @staticmethod
    def _resolve_atoms(shx, restr: Restraint):
        atoms = restr.atoms
        for atnum, ap in enumerate(atoms):
            if not isinstance(ap, (list, tuple)):
                # ignores all ranges: ['O1', '>', 'F9']
                try:
                    atoms[atnum] = shx.atoms.get_atom_by_name(ap)
                except TypeError:
                    continue
                continue
            for num, atname in enumerate(ap):
                # without range: [['O1', 'C1'], ...]
                try:
                    atoms[atnum][num] = shx.atoms.get_atom_by_name(atname)
                except TypeError:
                    continue


class DEFS(Restraint):

    def __init__(self, shx, spline: list):
        """
        DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
        Changes the *default* effective standard deviations for the following
        DFIX, SAME, SADI, CHIV, FLAT, DELU and SIMU restraints.
        """
        super(DEFS, self).__init__(shx, spline)
        self.sf = 0.1
        self.su = 0.01
        self.ss = 0.04
        self.maxsof = 1
        self.sd = 0.02
        self.active = True
        p, _ = self._parse_line(spline)
        if _:
            raise ParseParamError
        if len(p) > 0:
            self.sd = p[0]
        if len(p) > 1:
            self.sf = p[1]
        if len(p) > 2:
            self.su = p[2]
        if len(p) > 3:
            self.ss = p[3]
        if len(p) > 4:
            self.maxsof = p[4]

    @property
    def all(self):
        return self.sd, self.sf, self.su, self.ss, self.maxsof


class NCSY(Restraint):
    """
    NCSY DN sd[0.1] su[0.05] atoms
    """

    def __init__(self, shx, spline: list):
        super(NCSY, self).__init__(shx, spline)
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

    def __init__(self, shx, spline: list):
        """
        ISOR s[0.1] st[0.2] atomnames
        """
        super(ISOR, self).__init__(shx, spline)
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

    def __init__(self, shx, spline: list):
        super(FLAT, self).__init__(shx, spline)
        self.s = 0.1
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        # TODO: Have to resolve ranges first:
        # if len(self.atoms) < 4:
        #    raise ParseParamError


class BUMP(Restraint):

    def __init__(self, shx, spline):
        """
        BUMP s [0.02]
        """
        super(BUMP, self).__init__(shx, spline)
        self.s = 0.02
        p, _ = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        if _:
            raise ParseParamError


class DFIX(Restraint):

    def __init__(self, shx, spline):
        """
        DFIX d s[0.02] atom pairs
        """
        super(DFIX, self).__init__(shx, spline)
        self.s = 0.02
        p, self.atoms = self._parse_line(spline, pairs=True)
        if len(p) > 0:
            self.d = p[0]
        if len(p) > 1:
            self.s = p[1]
        self._paircheck()
        if not self.d:
            raise ParseNumError
        if DEBUG and 0.0001 < self.d <= self.s:
            print('*** WRONG ODER of INSTRUCTIONS. d is smaller than s ***')
            print("{}".format(self.textline))


class DANG(Restraint):

    def __init__(self, shx, spline):
        """
        DANG d s[0.04] atom pairs
        """
        super(DANG, self).__init__(shx, spline)
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

    def __init__(self, shx, spline):
        """
        SADI s[0.02] pairs of atoms
        Instructions with only two atoms are ignored by SHELXL: SADI C3 C4
        SADI_3 C3 C4 C3_4 C4_4 creates SADI C3_3 C4_3  C3_4 C4_4
        SADI_CCF3 C3 C4 C3_4 C4_4 creates SADI C3_1 C4_1  C3_2 C4_2  C3_4 C4_4 if there are residues 1, 2 and 4
        """
        super(SADI, self).__init__(shx, spline)
        self.s = 0.02
        p, self.atoms = self._parse_line(spline, pairs=True)
        if len(p) > 0:
            self.s = p[0]
        self._paircheck()


class SAME(Restraint):

    def __init__(self, shx, spline):
        """
        SAME s1[0.02] s2[0.04] atomnames
        """
        super(SAME, self).__init__(shx, spline)
        self.s1 = 0.02
        self.s2 = 0.04
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s1 = p[0]
        if len(p) > 1:
            self.s2 = p[1]


class RIGU(Restraint):

    def __init__(self, shx, spline: list):
        """
        RIGU s1[0.004] s2[0.004] atomnames
        """
        super(RIGU, self).__init__(shx, spline)
        self.s1 = 0.004
        self.s2 = 0.004
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s1 = p[0]
        if len(p) > 1:
            self.s2 = p[1]


class SIMU(Restraint):

    def __init__(self, shx, spline: list):
        """
        SIMU s[0.04] st[0.08] dmax[2.0] atomnames
        """
        super(SIMU, self).__init__(shx, spline)
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

    def __init__(self, shx, spline: list):
        """
        DELU s1[0.01] s2[0.01] atomnames
        """
        super(DELU, self).__init__(shx, spline)
        self.s1 = 0.01
        self.s2 = 0.01
        p, self.atoms = self._parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s1 = p[0]
        if len(p) > 1:
            self.s2 = p[1]


class CHIV(Restraint):

    def __init__(self, shx, spline: list):
        """
        CHIV V[0] s[0.1] atomnames
        """
        super(CHIV, self).__init__(shx, spline)
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

    def __init__(self, shx, spline: list) -> None:
        super(EADP, self).__init__(shx, spline)
        _, self.atoms = self._parse_line(spline, pairs=False)


class EXYZ(Restraint):
    """
    EADP atomnames
    """

    def __init__(self, shx, spline: list) -> None:
        super(EXYZ, self).__init__(shx, spline)
        _, self.atoms = self._parse_line(spline, pairs=False)


class DAMP(Command):
    """
    DAMP damp[0.7] limse[15]
    """

    def __init__(self, shx, spline: list):
        super(DAMP, self).__init__(shx, spline)
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

    def __init__(self, shx, spline: list):
        """
        HFIX mn U[#] d[#] atomnames
        """
        super(HFIX, self).__init__(shx, spline)
        self.params, self.atoms = self._parse_line(spline, intnums=True)

    def __repr__(self):
        return "HFIX {} {}".format(" ".join([str(x) for x in self.params]) if self.params else '',
                                   " ".join(self.atoms) if self.atoms else '')


class HKLF(Command):

    def __init__(self, shx, spline: list):
        """
        HKLF N[0] S[1] r11...r33[1 0 0 0 1 0 0 0 1] sm[1] m[0]
        """
        super(HKLF, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.n = 0
        self.s = 1
        self.matrix = [1, 0, 0, 0, 1, 0, 0, 0, 1]
        self.sm = 1
        self.m = 0
        if len(p) > 0:
            self.n = int(p[0])
        if len(p) > 1:
            self.s = p[1]
        if len(p) > 10:
            self.matrix = p[3:11]
        if len(p) > 11:
            self.sm = p[12]
        if len(p) > 12:
            self.m = p[13]

    def __repr__(self):
        return "HKLF {:,g} {:,g}  {}  {:,g} {:,g}".format(self.n, self.s, ' '.join([str(i) for i in self.matrix]),
                                                          self.sm, self.m)


class SUMP(Command):
    """
    SUMP for linear equation eypressions with free variables.
    SUMP c sigma c1 m1 c2 m2 ...
    """

    def __init__(self, shx, spline: list):
        super(SUMP, self).__init__(shx, spline)
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


class LATT(Command):
    lattdict = {1: [],  # Primitive
                2: [SymmetryElement(['0.5', '0.5', '0.5'])],  # I-centered
                3: [SymmetryElement(['1/3', '2/3', '2/3']),  # Rhombohedral
                    SymmetryElement(['2/3', '1/3', '1/3'])],
                4: [SymmetryElement(['0.0', '0.5', '0.5']),  # F-centered
                    SymmetryElement(['0.5', '0.0', '0.5']),
                    SymmetryElement(['0.5', '0.5', '0.0'])],
                5: [SymmetryElement(['0.0', '0.5', '0.5'])],  # A-centered
                6: [SymmetryElement(['0.5', '0.0', '0.5'])],  # B-centered
                7: [SymmetryElement(['0.5', '0.5', '0.0'])]}  # C-centered

    lattint_to_str = {1: 'P', 2: 'I', 3: 'R', 4: 'F', 5: 'A', 6: 'B', 7: 'C'}

    def __init__(self, shx, spline: list):
        """
        LATT N[1]
        """
        super(LATT, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self.centric = False
        try:
            self.N = int(p[0])
        except ValueError:
            self.N = -1
        self.N_str = self.lattint_to_str[abs(self.N)]
        if self.N > 0:  # centrosymmetric space group:
            self.centric = True
        self.latt_ops = LATT.lattdict[abs(self.N)]


class SYMM(Command):

    def __init__(self, shx, spline: list):
        """
        SYMM symmetry operation
        """
        super(SYMM, self).__init__(shx, spline)
        self.symmcard = self._parse_line(spline)

    def _parse_line(self, spline: list, intnums: bool = False) -> list:
        symmcard = ''.join(spline[1:]).split(',')  # removes whitespace
        return symmcard

    def _as_str(self):
        return "SYMM  " + ", ".join(self.symmcard)

    def __repr__(self):
        return self._as_str()

    def __str__(self) -> str:
        return self._as_str()


class SymmCards():
    """
    Contains the list of SYMM cards
    """

    def __init__(self, shx):
        self.shx = shx
        self._symmcards = [SymmetryElement(['X', 'Y', 'Z'])]
        self.latt_ops = []

    def _as_str(self) -> str:
        return "\n".join([str(x) for x in self._symmcards])

    def __repr__(self) -> str:
        return self._as_str()

    def __str__(self) -> str:
        return self._as_str()

    def __getitem__(self, item):
        return self._symmcards[item]

    def __iter__(self):
        for x in self._symmcards:
            yield x

    def append(self, symm_data: list) -> None:
        """
        Add the content of a Shelxl SYMM command to generate the appropriate SymmetryElement instance.
        :param symm_data: list of strings. eg.['1/2+X', '1/2+Y', '1/2+Z']
        :return: None
        """
        new_symm = SymmetryElement(symm_data)
        self._symmcards.append(new_symm)
        for symm in self.shx.latt.latt_ops:
            latt_symm = new_symm.apply_latt_symm(symm)
            if latt_symm not in self._symmcards:
                self._symmcards.append(latt_symm)
        if self.shx.latt.centric:
            self._symmcards.append(SymmetryElement(symm_data, centric=True))
            for symm in self.shx.latt.latt_ops:
                latt_symm = new_symm.apply_latt_symm(symm)
                latt_symm.centric = True
                self._symmcards.append(latt_symm)

    def set_centric(self, value: bool):
        """
        Defines the instance as representing a centrosymmetric structure. Generates the appropriate SymmetryElement
        instances automatically if called before adding further SYMM commands via self.addSymm().
        """
        self.shx.latt.centric = value
        self._symmcards.append(SymmetryElement(['-X', '-Y', '-Z']))
        self._symmcards[-1].centric = True

    def set_latt_ops(self, lattops: list) -> None:
        """
        Adds lattice operations. If called before adding SYMM commands, the appropriate lattice operations are used
        automatically to generate further SymmetryElements.
        :param lattops: list of SymmetryElement instances.
        """
        self.latt_ops = lattops


class LSCycles(Command):
    def __init__(self, shx, spline: list):
        """
        L.S. nls[0] nrf[0] nextra[0]
        #TODO: support nls parameter
        If nrf is positive, it is the number of these cycles that should be performed before applying ANIS.
        Negative nrf indicates which reflections should be ignored during the refinement but used instead for
        the calculation of free R-factors in the final structure factor summation.
        nextra is the number of additional parameters that were derived from the data when 'squeezing' the
        structure etc.
        """
        super(LSCycles, self).__init__(shx, spline)
        p, _ = self._parse_line(spline)
        self._shx = shx
        self.cgls = False
        self._cycles = 0
        self._nrf = ''
        self._nextra = ''
        try:
            self._cycles = int(p[0])
        except (IndexError, NameError, ValueError):
            raise ParseNumError
        try:
            self._nrf = int(p[1])
        except IndexError:
            pass
        try:
            self._nextra = int(p[2])
        except IndexError:
            pass
        if spline[0].upper() == 'CGLS':
            self.cgls = True

    @property
    def number(self):
        return self._cycles

    @number.setter
    def number(self, n: int) -> None:
        self._cycles = int(n)
        self.__init__(self._shx, self._as_str().split())

    def set_refine_cycles(self, number: int) -> None:
        """
        Sets the number of refinement cycles for the current res file.
        """
        self.number = int(number)

    @property
    def text(self) -> str:
        """
        'CGLS 10 2 '
        """
        return self.__repr__()

    def __iter__(self):
        for x in self.__repr__().split():
            yield x

    def _as_str(self):
        return '{} {} {} {}'.format('CGLS' if self.cgls else 'L.S.', self._cycles,
                                    self._nrf if self._nrf else '', self._nextra if self._nextra else '').strip()

    def __repr__(self):
        return self._as_str()


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
        sftext = ''
        elements = []
        for sf in self.sfac_table:
            if not self.is_exp(sf) and sf['element'].capitalize() not in elements:
                elements.append(sf['element'].capitalize())
            else:
                if elements:
                    sftext = self._extend_sfac_text(elements, sftext)
                    elements = []
                values = []
                for x in ('element', 'a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'a4', 'b4', 'c',
                          'fprime', 'fdprime', 'mu', 'r', 'wt'):
                    values.append(sf[x])
                sftext += "\nSFAC " + "  ".join(values)
        if elements:
            sftext = self._extend_sfac_text(elements, sftext)
        return sftext[1:]

    def _extend_sfac_text(self, elements: List[str], sftext: str) -> str:
        sftext += "\nSFAC " + "  ".join(elements)
        return sftext

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
        self.shx.unit.add_number(1.0)

    def remove_element(self, element: str):
        del self.sfac_table[self.shx.elem2sfac(element.upper()) - 1]
        del self.elements_list[self.elements_list.index(element.upper())]


class UNIT(Command):

    def __init__(self, shx, spline: list):
        """
        UNIT n1 n2 ...
        """
        super(UNIT, self).__init__(shx, spline)
        self.values, _ = self._parse_line(spline)

    def add_number(self, number: float):
        self.values.append(number)

    def __iter__(self):
        yield [x for x in self.values]

    def __repr__(self):
        return "UNIT " + "  ".join(["{:,g}".format(x) for x in self.values])

    def __str__(self):
        return self.__repr__()

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

    def __init__(self, shx, spline):
        super(BASF, self).__init__(shx, spline)
        self.scale_factors, _ = self._parse_line(spline)

    def __iter__(self):
        yield self.scale_factors


class TWIN(Command):

    def __init__(self, shx, spline: list):
        """
        TWIN 3x3 matrix [-1 0 0 0 -1 0 0 0 -1] N[2]
        +N     -N  m = |N|
        m-1 to 2m-1
        m-1   (2*abs(m)/2)-1
        """
        super(TWIN, self).__init__(shx, spline)
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

    def __init__(self, shx, spline: list):
        """
        The weighting scheme is defined as follows:
        w = q / [ σ²(Fo²) + (a*P)² + b*P + d + e*sin(θ)/$lambda; ]

        WGHT a[0.1] b[0] c[0] d[0] e[0] f[.33333]
        Usually only WGHT a b
        """
        super(WGHT, self).__init__(shx, spline)
        self.shx = shx
        self.a = 0.1
        self.b = 0.0
        self.c = 0.0
        self.d = 0.0
        self.e = 0.0
        self.f = 0.33333
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

    def _as_string(self):
        wght = 'WGHT   {} {}'.format(self.a, self.b)
        # It is very unlikely that someone changes other parameter than a and b:
        if (self.c + self.d + self.e + self.f) != 0.33333:
            wght += ' {} {} {} {}'.format(self.c, self.d, self.e, self.f)
        return wght

    def difference(self) -> List[float]:
        """
        Returns a list with the weight differences of the parameters a and b.
        """
        try:
            adiff = abs(self.shx.wght.a - self.shx.wght_suggested.a)
            bdiff = abs(self.shx.wght.b - self.shx.wght_suggested.b)
        except AttributeError:
            print("No suggested weighting scheme found. Unable to proceed.")
            return [0.0, 0.0]
        return [round(adiff, 3), round(bdiff, 3)]

    def __repr__(self):
        return self._as_string()

    def __str__(self):
        return self._as_string()
