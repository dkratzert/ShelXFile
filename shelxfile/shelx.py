# -*- encoding: utf-8 -*-
# möpß
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
"""
This is a full implementation of the SHELXL file syntax. Additionally it is able to edit SHELX properties with Python.
The implementation is Python3-only and supports SHELXL after 2017 (You should not use old versions anyway).
"""
import os
import re
import sys

from refine.shx_refine import ShelxlRefine
from shelxfile.cards import ACTA, FVAR, FVARs, REM, BOND, Restraints, DEFS, NCSY, ISOR, FLAT, \
    BUMP, DFIX, DANG, SADI, SAME, RIGU, SIMU, DELU, CHIV, EADP, EXYZ, DAMP, HFIX, HKLF, SUMP, SYMM, LSCycles, \
    SFACTable, UNIT, BASF, TWIN, WGHT, BLOC, SymmCards, CONN, CONF, BIND, DISP, GRID, HTAB, MERG, FRAG, FREE, FMAP, \
    MOVE, PLAN, PRIG, RTAB, SHEL, SIZE, SPEC, STIR, TWST, WIGL, WPDB, XNPD, ZERR, CELL, LATT, MORE, MPLA, AFIX, PART, \
    RESI, ABIN, ANIS
from shelxfile.atoms import Atoms, Atom
from misc import DEBUG, ParseOrderError, ParseNumError, ParseUnknownParam, \
    split_fvar_and_parameter, flatten, time_this_method, multiline_test, dsr_regex, wrap_line

from math import radians, cos, sin, sqrt

from dsrmath import Matrix

"""
TODO:

- Atoms.add_atom(position=None) method. default for position is after FVAR table
- fit fragment without shelxl
- check if atoms in restraints are also in structure
- restraints involved with an atom should also be part of the atoms properties
------------------------------------------------------------
- deleting atoms should also remove them from restraints
- add remove_hydrogen_atoms(atom) method.
- shx.remove_all_H([list of atoms], or all)
- bond list
- shx.update_weight
- shx.weight_difference
- shx.atoms.angle(at1, at2, at3)
- shx.atoms.tors(at1, at2, at3, at4)
- shx.atom.change_type('xx')
- implement an add_afix(afixm, afixn, atoms, frag_fend=False, position=None, afix_options=None)
  default position is directly behind FVAR or FRAG/FEND if enabled
- will read in lst file after refinement to fill shx.lst_file properties.
- shx.move_in_part([list of atoms])
- shx.move_in_resi([list of atoms])
- shx.sort.residue(3) (low priority)
- shx.lst_file.residuals -> density, r-values, goof, movements, bad reflections
   bad restraints, shelx *** errors
- shx.unused_atom_name('C') -> get a carbon atom with unused number
- shx.sort -> sort file (low priority)
"""

SHX_CARDS = ('TITL', 'CELL', 'ZERR', 'LATT', 'SYMM', 'SFAC', 'UNIT', 'LIST', 'L.S.', 'CGLS',
             'BOND', 'FMAP', 'PLAN', 'TEMP', 'ACTA', 'CONF', 'SIMU', 'RIGU', 'WGHT', 'FVAR',
             'DELU', 'SAME', 'DISP', 'LAUE', 'REM ', 'MORE', 'TIME', 'END ', 'HKLF', 'OMIT',
             'SHEL', 'BASF', 'TWIN', 'EXTI', 'SWAT', 'HOPE', 'MERG', 'SPEC', 'RESI', 'MOVE',
             'ANIS', 'AFIX', 'HFIX', 'FRAG', 'FEND', 'EXYZ', 'EADP', 'EQIV', 'CONN', 'BIND',
             'FREE', 'DFIX', 'BUMP', 'SADI', 'CHIV', 'FLAT', 'DEFS', 'ISOR', 'NCSY', 'SUMP',
             'BLOC', 'DAMP', 'STIR', 'MPLA', 'RTAB', 'HTAB', 'SIZE', 'WPDB', 'GRID', 'MOLE',
             'XNPD', 'REST', 'CHAN', 'FLAP', 'RNUM', 'SOCC', 'PRIG', 'WIGL', 'RANG', 'TANG',
             'ADDA', 'STAG', 'NEUT', 'ABIN', 'ANSC', 'ANSR', 'NOTR', 'TWST', 'PART', 'DANG',
             'BEDE', 'LONE', 'REM', 'END')


class ShelXFile():
    """
    Class for data from a SHELXL res file. Includes Atoms, cards and unit cell.

    :type restraints: List[Restraint]
    """
    delete_on_write = None
    atoms = None
    sump = []
    _r1_regex = re.compile(r'^REM\s+R1\s+=', re.IGNORECASE)
    _wr2_regex = re.compile(r'^REM\s+wR2\s+=', re.IGNORECASE)
    _parameters_regex = re.compile(r'REM\s+\d+\s+parameters\s+refined', re.IGNORECASE)

    def __init__(self: 'ShelXFile', resfile: str):
        """
        Reads the shelx file and extracts information.

        :param resfile: file path
        """
        self.temp_in_Kelvin = 0.0
        self.shelx_max_line_length = 79  # maximum character lenth per line in SHELXL
        self.nohkl = False
        self._a, self._b, self._c, self._alpha, self._beta, self._gamma, self.V = \
            None, None, None, None, None, None, None
        self.cell = None
        self.ansc = None
        self.abin = None
        self.acta = None
        self.fmap = None
        self.xnpd = None
        self.wpdb = None
        self.wigl = None
        self.temp = 20
        self.swat = None
        self.stir = None
        self.spec = None
        self.twst = None
        self.plan = None
        self.prig = None
        self.merg = None
        self.more = None
        self.move = None
        self.defs = None
        self.zerr = None
        self.wght = None
        self.frag = None
        self.twin = None
        self.basf = None
        self.latt = None
        self.anis = None
        self.damp = None
        self.unit = None
        self.R1 = None
        self.wR2 = None
        self.data = None
        self.parameters = None
        self.dat_to_param = None
        self.num_restraints = None
        self.sump = None
        self.end = False
        self.maxsof = 1.0
        self.size = None
        self.htab = None
        self.shel = None
        self.mpla = None
        self.rtab = None
        self.omit = None
        self.hklf = None
        self.grid = None
        self.free = None
        self.titl = ""
        self.exti = 0
        self.eqiv = None
        self.disp = None
        self.conn = None
        self.conv = None
        self.bind = None
        self.ansr = 0.001
        self.bloc = None
        self.afix = None  # AFIX(self, [''])
        self.part = PART(self, ['PART',  '0'])
        self.resi = RESI(self, ['RESI',  '0'])
        self.dsrlines = []
        self.dsrline_nums = []
        self.symmcards = SymmCards(self)
        self.hfixes = []
        self.Z = 1
        self.rem = []
        self.indexes = {}
        self.atoms = Atoms(self)
        self.fvars = FVARs(self)
        self.restraints = Restraints()
        self.sfac_table = SFACTable(self)
        self.delete_on_write = set()
        self.wavelen = None
        self.global_sadi = None
        self.cycles = None
        self.list = 0
        self.theta_full = 0
        self.non_h = None
        self.error_line_num = -1  # Only used to tell the line number during an exception.
        self.restrdict = {}
        self.resfile = os.path.abspath(resfile)
        if DEBUG:
            print('Resfile is:', self.resfile)
        try:
            self._reslist = self.read_file_to_list(self.resfile)
        except UnicodeDecodeError:
            if DEBUG:
                print('*** Unable to read file', self.resfile, '***')
            return
        try:
            self.parse_shx_file()
            pass
        except Exception as e:
            if DEBUG:
                print(e)
                print("*** Syntax error found in file {}, line {} ***".format(self.resfile, self.error_line_num + 1))
            if DEBUG:
                raise
            else:
                print('*** Parser Error ***')
                return
        else:
            self.run_after_parse()

    @time_this_method
    def parse_shx_file(self):
        """
        Extracts the atoms and other information from the res file.

        line is upper() after multiline_test()
        spline is as in .res file.
        """
        lastcard = ''
        sof = 0
        fvarnum = 1
        for line_num, line in enumerate(self._reslist):
            self.error_line_num = line_num  # For exception during parsing.
            list_of_lines = [line_num]  # list of lines where a card appears, e.g. for atoms with two lines
            if line[:1] == ' ' or line == '':
                continue
            if not self.titl and line[:4] == 'TITL':
                # TITL[]  ->  = and ! can be part of the TITL!
                self.titl = line[5:76]
                lastcard = 'TITL'
                continue
            wrapindex = 0
            # This while loop makes wrapped lines look like they are not wrapped. The following lines are then
            # beginning with a space character and thus are ignored. The 'lines' list holds the line nnumbers where
            # 'line' is located ([line_num]) plus the wrapped lines.
            while multiline_test(self._reslist[line_num + wrapindex]):
                # Glue together the two lines wrapped with "=":
                wrapindex += 1
                line = line.rpartition('=')[0] + self._reslist[line_num + wrapindex]
                self.delete_on_write.update([line_num + wrapindex])
                self._reslist[line_num + wrapindex] = ''
                list_of_lines.append(line_num + wrapindex)  # list containing the lines of a multiline command
            # The current line splitted:
            spline = line.split('!')[0].split()  # Ignore comments with "!", see how this performes
            # The current line as string:
            line = line.upper().split('!')[0]  # Ignore comments with "!", see how this performes
            word = line[:4]
            # get RESI:
            if line.startswith(('END', 'HKLF')) and self.resi:
                self.resi.num = 0
                if DEBUG:
                    print('RESI in line {} was not closed'.format(line_num + 1))
                continue
            if line.startswith('RESI'):
                self.resi = RESI(self, spline)
                self.assign_card(self.resi, line_num)
                continue
            # Now collect the PART:
            if line.startswith(('END', 'HKLF')) and self.part:
                self.part.n = 0
                if DEBUG:
                    print('PART in line {} was not closed'.format(line_num + 1))
                continue
            if line.startswith('PART'):
                self.part = PART(self, spline)
                self.assign_card(self.part, line_num)
                continue
            # collect AFIX:
            if line.startswith(('END', 'HKLF')) and self.afix:
                self.afix.mn = 0
                if DEBUG:
                    print('AFIX in line {} was not closed'.format(line_num + 1))
                continue
            elif line.startswith('AFIX'):
                self.afix = AFIX(self, spline)
                self.assign_card(self.afix, line_num)
                continue
            elif self.is_atom(line):
                # A SHELXL atom:
                # F9    4    0.395366   0.177026   0.601546  21.00000   0.03231  ( 0.03248 =
                #            0.03649  -0.00522  -0.01212   0.00157 )
                a = Atom(self, spline, list_of_lines, line_num, part=self.part, afix=self.afix, resi=self.resi, sof=sof)
                self.append_card(self.atoms, a, line_num)
                continue
            elif self.end:
                # Prevents e.g. parsing of second WGHT after END:
                continue
            elif word == 'SADI':
                # SADI s[0.02] pairs of atoms
                # or SADI
                if len(spline) == 1:
                    self.global_sadi = line_num
                self.append_card(self.restraints, SADI(self, spline), line_num)
                continue
            elif word == 'DFIX':
                # DFIX d s[0.02] atom pairs
                self.append_card(self.restraints, DFIX(self, spline), line_num)
                continue
            elif word == 'SIMU':
                # SIMU s[0.04] st[0.08] dmax[2.0] atomnames
                self.append_card(self.restraints, SIMU(self, spline), line_num)
                continue
            elif word == 'DELU':
                # DELU s1[0.01] s2[0.01] atomnames
                self.append_card(self.restraints, DELU(self, spline), line_num)
                continue
            elif word == 'RIGU':
                # RIGU s1[0.004] s2[0.004] atomnames
                self.append_card(self.restraints, RIGU(self, spline), line_num)
                continue
            elif word == 'BASF':
                # BASF scale factors
                self.append_card(self.restraints, BASF(spline, list_of_lines), line_num)
                continue
            elif word == 'HFIX':
                # HFIX mn U[#] d[#] atomnames
                self.append_card(self.hfixes, HFIX(spline, list_of_lines), line_num)
                continue
            elif word == 'DANG':
                # DANG d s[0.04] atom pairs
                self.append_card(self.restraints, DANG(self, spline), line_num)
                continue
            elif word == 'EADP':
                self.append_card(self.restraints, EADP(self, spline), line_num)
                continue
            elif line[:3] == 'REM':
                if dsr_regex.match(line):
                    self.dsrlines.append(" ".join(spline))
                    self.dsrline_nums.extend(list_of_lines)
                self.append_card(self.rem, REM(self, spline), line_num)
                self._get_residuals(spline, line)
                continue
            elif word == 'AFIX':
                pass
            elif word == 'CELL':
                # CELL λ a b c α β γ
                if not lastcard == 'TITL':
                    if DEBUG:
                        print('TITL is missing.')
                    # raise ParseOrderError
                self.cell = CELL(self, spline)
                self.assign_card(self.cell, line_num)
                self._a, self._b, self._c, self._alpha, self._beta, self._gamma = self.cell.cell_list
                self.V = self.vol_unitcell(self._a, self._b, self._c, self._alpha, self._beta, self._gamma)
                self.cell.volume = self.V
                # self.A = self.orthogonal_matrix()
                self.wavelen = self.cell.wavelen
                lastcard = 'CELL'
                continue
            elif word == "ZERR":
                # ZERR Z esd(a) esd(b) esd(c) esd(α) esd(β) esd(γ)
                if not lastcard == 'CELL':
                    if DEBUG:
                        print('*** Invalid SHELX file!')
                    raise ParseOrderError
                if not self.cell:
                    raise ParseOrderError('*** Cell parameters missing! ***')
                if len(spline) >= 8:
                    self.zerr = ZERR(self, spline)
                    self.Z = self.zerr.Z
                    self.assign_card(self.zerr, line_num)
                lastcard = 'ZERR'
                continue
            elif word == "SYMM":
                # SYMM symmetry operation
                #  Being more greedy, because many files do this wrong:
                # if not lastcard == 'ZERR':
                #    raise ParseOrderError
                # if not self.zerr:
                #    raise ParseOrderError
                s = SYMM(self, spline)
                if self.latt.centric:
                    self.symmcards.set_centric(True)
                self.symmcards.append(s.symmcard)
                if s not in self._reslist:
                    self._reslist[line_num] = s
                else:
                    self.delete_on_write.update([line_num])
                    self._reslist[line_num] = ' '
                lastcard = 'SYMM'
                continue
            elif word == 'SFAC':
                # SFAC elements or
                # SFAC E a1 b1 a2 b2 a3 b3 a4 b4 c f' f" mu r wt
                # Being less strict to be able to parse files without cell errors:
                # if not (lastcard == 'LATT' or lastcard == 'ZERR'):
                #    raise ParseOrderError
                # if not self.symmcards:
                #    raise ParseOrderError
                if len(spline) <= 1:
                    continue
                self.sfac_table.parse_element_line(spline)
                if self.sfac_table not in self._reslist:
                    self._reslist[line_num] = self.sfac_table
                else:
                    self.delete_on_write.update([line_num])
                    self._reslist[line_num] = ' '
                lastcard = 'SFAC'
                continue
            elif word == 'UNIT':
                # UNIT n1 n2 ...
                # Number of atoms of each type in the unit-cell, in SFAC order.
                if not lastcard == 'SFAC':
                    raise ParseOrderError
                if self.sfac_table:
                    try:
                        self.unit = self.assign_card(UNIT(self, spline), line_num)
                    except ValueError:
                        if DEBUG:
                            print('*** Non-numeric value in SFAC instruction! ***')
                        raise
                else:
                    raise ParseOrderError
                if len(self.unit.values) != len(self.sfac_table.elements_list):
                    raise ParseNumError
                lastcard = 'UNIT'
                continue
            elif word == "LATT":
                # LATT N[1]
                # 1=P, 2=I, 3=rhombohedral obverse on hexagonal axes, 4=F, 5=A, 6=B, 7=C.
                # negative is non-centrosymmetric
                self.latt = LATT(self, spline)
                self.assign_card(self.latt, line_num)
                if not lastcard == 'ZERR':
                    if DEBUG:
                        print('*** ZERR instruction is missing! ***')
                    # raise ParseOrderError
                continue
            elif word in ['L.S.', 'CGLS']:
                # CGLS nls[0] nrf[0] nextra[0]
                # L.S. nls[0] nrf[0] nextra[0]
                self.cycles = self.assign_card(LSCycles(self, spline, line_num), line_num)
                continue
            elif word == "LIST":
                # LIST m[#] mult[1] (mult is for list 4 only)
                self.list = int(spline[1])
                continue
            elif word == "FVAR":
                # FVAR osf[1] free variables
                for fvvalue in spline[1:]:
                    fvarnum += 1
                    self.append_card(self.fvars, FVAR(fvarnum, float(fvvalue)), line_num)
                    if self.fvars not in self._reslist:
                        self._reslist[line_num] = self.fvars
                    else:
                        self.delete_on_write.update([line_num])
            elif word == 'ANIS':
                # ANIS n or ANIS names
                # Must be before Atom(), to know which atom is anis.
                self.anis = ANIS(self, spline)
                self.assign_card(self.anis, line_num)
                continue
            elif word == 'WGHT':
                # WGHT a[0.1] b[0] c[0] d[0] e[0] f[.33333]
                self.wght = self.assign_card(WGHT(self, spline), line_num)
                continue
            elif word == 'ACTA':
                # ACTA 2θfull[#] -> optional parameter NOHKL
                self.acta = ACTA(self, spline)
                self.assign_card(self.acta, line_num)
                continue
            elif word == 'DAMP':
                # DAMP damp[0.7] limse[15]
                self.damp = DAMP(self, spline)
                self.assign_card(self.damp, line_num)
                continue
            elif word == 'ABIN':
                # ABIN n1 n2
                self.abin = ABIN(self, spline)
                self.assign_card(self.abin, line_num)
                continue
            elif word == 'ANSC':
                # ANSC six coefficients
                if len(spline) == 7:
                    self.ansc = [float(x) for x in spline[:1]]
                continue
            elif word == 'ANSR':
                # ANSR anres[0.001]
                if len(spline) == 2:
                    self.ansr = float(spline[1])
                continue
            elif word == 'BIND':
                # BIND atom1 atom2
                if len(spline) == 3:
                    self.append_card(self.bind , BIND(self, spline), line_num)
                continue
            elif word == 'BLOC':
                # BLOC n1 n2 atomnames
                self.append_card(self.bloc, BLOC(self, spline), line_num)
                continue
            elif word == 'BOND':
                # BOND atomnames
                self.bond =BOND(self, spline) 
                self.assign_card(self.bond, line_num)
                continue
            elif word == 'BUMP':
                # BUMP s [0.02]
                self.append_card(self.restraints, BUMP(self, spline), line_num)
                continue
            elif word == 'CHIV':
                # CHIV V[0] s[0.1] atomnames
                self.append_card(self.restraints, CHIV(self, spline), line_num)
                continue
            elif word == 'CONF':
                # CONF atomnames max_d[1.9] max_a[170]
                self.conf = CONF(self, spline) 
                self.assign_card(self.conf, line_num)
                continue
            elif word == 'CONN':
                # CONN bmax[12] r[#] atomnames or CONN bmax[12]
                # bonded are d < (r1 + r2 + 0.5) Å
                self.conn = CONN(self, spline) 
                self.assign_card(self.conn, line_num)
                continue
            elif word == 'DEFS':
                # DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
                self.defs = DEFS(self, spline)
                self.assign_card(self.defs, line_num)
                continue
            elif word == 'DISP':
                # DISP E f' f"[#] mu[#]
                if not lastcard == 'SFAC':
                    raise ParseOrderError
                self.append_card(self.disp, DISP(self, spline), line_num)
                continue
            elif word == 'EQIV':
                # EQIV $n symmetry operation
                if len(spline) > 1:
                    if spline[1].startswith('$'):
                        self.eqiv.append(spline[1:])
                continue
            elif word == 'EXTI':
                # EXTI x[0]
                self.exti = float(spline[1])
                continue
            elif word == 'EXYZ':
                # EXYZ atomnames
                self.append_card(self.restraints, EXYZ(self, spline), line_num)
                continue
            elif word == 'FRAG':
                # FRAG code[17] a[1] b[1] c[1] α[90] β[90] γ[90]
                if len(spline) == 8:
                    self.frag = FRAG(self, spline)
                    self.assign_card(self.frag, line_num)
                continue
            elif word == 'FEND':
                # FEND (must follow FRAG)
                if not self.frag:
                    raise ParseOrderError
                self.frag = None  # Turns frag mode off.
                continue
            elif word == 'FLAT':
                # FLAT s[0.1] four or more atoms
                self.append_card(self.restraints, FLAT(self, spline), line_num)
                continue
            elif word == 'FREE':
                # FREE atom1 atom2
                free = FREE(self, spline)
                self.free.append(free)
                continue
            elif word == 'GRID':
                # GRID sl[#] sa[#] sd[#] dl[#] da[#] dd[#]
                self.grid = GRID(self, spline)
                self.assign_card(self.grid, line_num)
                continue
            elif word == 'HKLF':
                # HKLF N[0] S[1] r11...r33[1 0 0 0 1 0 0 0 1] sm[1] m[0]
                self.hklf = HKLF(self, spline)
                self.assign_card(self.hklf, line_num)
                continue
            elif line.startswith('END'):
                # END (after HKLF or ends an include file)
                self.end = True
                continue
            elif word == 'HTAB':
                # HTAB dh[2.0]  or  HTAB donor-atom acceptor-atom
                self.htab = HTAB(self, spline)
                self.assign_card(self.htab, line_num)
                continue
            elif word == 'ISOR':
                # ISOR s[0.1] st[0.2] atomnames
                self.append_card(self.restraints, ISOR(self, spline), line_num)
                continue
            elif word == 'LAUE':
                # LAUE E
                # I completely do not understand the LAUE instruction description in the manual!
                continue
            elif word == 'MERG':
                # MERG n[2]
                self.merg = MERG(self, spline)
                self.assign_card(self.merg, line_num)
                continue
            elif word == 'MORE':
                # MORE m[1]
                self.more = MORE(self, spline)
                self.assign_card(self.more, line_num)
                continue
            elif word == 'FMAP':
                # FMAP code[2] axis[#] nl[53]
                self.fmap = FMAP(self, spline)
                self.assign_card(self.fmap, line_num)
                continue
            elif word == 'MOVE':
                # MOVE dx[0] dy[0] dz[0] sign[1]
                self.move = MOVE(self, spline)
                self.assign_card(self.move, line_num)
                continue
            elif word == 'MPLA':
                # MPLA na atomnames
                self.mpla = MPLA(self, spline)
                self.assign_card(self.mpla, line_num)
                continue
            elif word == 'NCSY':
                # NCSY DN sd[0.1] su[0.05] atoms
                self.append_card(self.restraints, NCSY(self, spline), line_num)
                continue
            elif word == 'NEUT':
                # NEUT
                if not lastcard == 'SYMM':
                    raise ParseOrderError
                continue
            elif word == 'OMIT':
                # OMIT atomnames  or  OMIT s[-2] 2θ(lim)[180]  or  OMIT h k l
                self.omit.append(spline[1:])
                continue
            elif word == 'PLAN':
                # PLAN npeaks[20] d1[#] d2[#]
                self.plan = PLAN(self, spline)
                self.assign_card(self.plan, line_num)
                continue
            elif word == 'PRIG':
                # PRIG p[#]
                self.prig = PRIG(self, spline)
                self.assign_card(self.prig, line_num)
                continue
            elif word == 'RTAB':
                # RTAB codename atomnames  -->  codename: e.g. 'omeg' gets tabualted in the lst
                self.append_card(self.rtab, RTAB(self, spline), line_num)
                continue
            elif word == 'SAME':
                # SAME s1[0.02] s2[0.04] atomnames
                self.append_card(self.restraints, SAME(self, spline), line_num)
                continue
            elif word == 'SHEL':
                # SHEL lowres[infinite] highres[0]
                self.shel = SHEL(self, spline)
                self.assign_card(self.shel, line_num)
                continue
            elif word == 'SIZE':
                # SIZE dx dy dz
                self.size = SIZE(self, spline)
                self.assign_card(self.size, line_num)
                continue
            elif word == 'SPEC':
                # SPEC del[0.2]
                if len(spline) > 1:
                    self.spec = SPEC(self, spline)
                    self.assign_card(self.spec, line_num)
                continue
            elif word == 'STIR':
                # STIR sres step[0.01]   -> stepwise improvement in the resolution sres
                self.stir = STIR(self, spline)
                self.assign_card(self.stir, line_num)
                continue
            elif word == 'SUMP':
                # SUMP c sigma c1 m1 c2 m2 ...
                self.append_card(self.sump, SUMP(self, spline), line_num)
                continue
            elif word == 'SWAT':
                # SWAT g[0] U[2]
                self.swat = spline[1:]
                continue
            elif word == 'TEMP':
                # TEMP T[20]  -> in Celsius
                self.temp = float(spline[1])
                self.temp_in_Kelvin = self.temp + 273.15
                continue
            elif word == 'TWIN':
                # TWIN 3x3 matrix [-1 0 0 0 -1 0 0 0 -1] N[2]
                self.twin = TWIN(self, spline)
                self.assign_card(self.twin, line_num)
                continue
            elif word == 'TWST':
                # TWST N[0] (N[1] after SHELXL-2018/3)
                if len(spline) > 1:
                    self.twst = TWST(self, spline)
                    self.assign_card(self.twst, line_num)
                continue
            elif word == 'WIGL':
                # WIGL del[0.2] dU[0.2]
                self.wigl = WIGL(self, spline)
                self.assign_card(self.wigl, line_num)
                continue
            elif word == 'WPDB':
                # WPDB n[1]
                self.wpdb = WPDB(self, spline)
                self.assign_card(self.wpdb, line_num)
                continue
            elif word == 'XNPD':
                # XNPD Umin[-0.001]
                self.xnpd = XNPD(self, spline)
                self.assign_card(self.xnpd, line_num)
                continue
            elif word == 'BEDE':
                # Later...
                continue
            elif word == 'LONE':
                # Later...
                continue
            elif word == 'MOLE':
                # print('*** MOLE is deprecated! Do not use it! ***')
                pass
            elif word == 'HOPE':
                # print('*** HOPE is deprecated! Do not use it! ***')
                pass
            elif line[:1] == '+':
                pass
            else:
                if DEBUG:
                    print(line)
                    raise ParseUnknownParam

    @time_this_method
    def run_after_parse(self):
        """
        Runs all what is left after parsing all lines. E.G. sanity checks.
        # TODO: check if EXYZ and anisotropy fit togeter
        """
        if self.sump:
            for x in self.sump:
                for y in x:
                    self.fvars.set_fvar_usage(int(y[1]))
        for r in self.restraints:
            if r.name == "DFIX" or r.name == "DANG":
                if abs(r.d) > 4:
                    fvar, value = split_fvar_and_parameter(r.d)
                    self.fvars.set_fvar_usage(fvar)
        if self.abin:
            if len(self.abin) > 1:
                self.fvars.set_fvar_usage(self.abin[0])
                self.fvars.set_fvar_usage(self.abin[1])
            else:
                self.fvars.set_fvar_usage(self.abin[0])
        # Check if basf parameters are consistent:
        if self.basf:
            if self.twin:
                basfs = flatten(self.basf)
                if int(self.twin.allowed_N) != len(basfs):
                    if DEBUG:
                        print('*** Invalid TWIN instruction! BASF with wrong number of parameters. ***')
        #for a in self.atoms:
        #    a.resolve_restraints()

    def restore_acta_card(self, acta: str):
        """
        Place ACTA after UNIT
        """
        self.acta = ACTA(self, acta.split())
        self._reslist.insert(self.unit.position + 1, self.acta)

    def orthogonal_matrix(self):
        """
        Converts von fractional to cartesian by .
        Invert the matrix to do the opposite.

        Old tests:
        #>>> import mpmath as mpm
        #>>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
        #>>> coord = (-0.186843,   0.282708,   0.526803)
        #>>> print(mpm.nstr(A*mpm.matrix(coord)))
        [-2.74151]
        [ 5.90959]
        [ 10.7752]
        #>>> cartcoord = mpm.matrix([['-2.74150542399906'], ['5.909586678'], ['10.7752007008937']])
        #>>> print(mpm.nstr(A**-1*cartcoord))
        [-0.186843]
        [ 0.282708]
        [ 0.526803]
        """
        return Matrix([[self._a, self._b * cos(self._gamma), self._c * cos(self._beta)],
                       [0, self._b * sin(self._gamma),
                        (self._c * (cos(self._alpha) - cos(self._beta) * cos(self._gamma)) / sin(self._gamma))],
                       [0, 0, self.V / (self._a * self._b * sin(self._gamma))]])

    @staticmethod
    def vol_unitcell(a, b, c, al, be, ga) -> float:
        """
        calculates the volume of a unit cell
        >>> v = ShelXFile.vol_unitcell(2, 2, 2, 90, 90, 90)
        >>> print(v)
        8.0
        """
        ca, cb, cg = cos(radians(al)), cos(radians(be)), cos(radians(ga))
        v = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca ** 2 - cb ** 2 - cg ** 2)
        return v

    def __repr__(self):
        """
        Represents the shelxl object.
        """
        resl = []
        for num, line in enumerate(self._reslist):
            if num in self.delete_on_write:
                if DEBUG:
                    pass
                    # print('Deleted line {}'.format(num + 1))
                continue
            try:
                if line == '' and self._reslist[num + 1] == '':
                    continue
            except IndexError:
                pass
            # Prevent wrapping long lines with \n breaks by splitting first:
            line = "\n".join([wrap_line(x) for x in str(line).split("\n")])
            resl.append(line)
        return "\n".join(resl)

    def write_shelx_file(self, filename=None, verbose=False):
        if not filename:
            filename = self.resfile
        with open(filename, 'w') as f:
            for num, line in enumerate(self._reslist):
                if num in self.delete_on_write:
                    if DEBUG:
                        pass
                        #print('Deleted line {}'.format(num + 1))
                    continue
                if line == '' and self._reslist[num + 1] == '':
                    continue
                # Prevent wrapping long lines with \n breaks by splitting first:
                line = "\n".join([wrap_line(x) for x in str(line).split("\n")])
                f.write(str(line) + '\n')
        if verbose or DEBUG:
            print('File successfully written to {}'.format(os.path.abspath(filename)))
            return True
        return True

    # @time_this_method
    def read_file_to_list(self, resfile: str) -> list:
        """
        Read in shelx file and returns a list without line endings. +include files are inserted
        also.
        :param resfile: The path to a SHLEL .res or .ins file.
        """
        reslist = []
        #resnodes = ResList()
        includefiles = []
        try:
            with open(resfile, 'r') as f:
                reslist = f.read().splitlines(keepends=False)
                #for ll in reslist:
                    #resnodes.append(ll)
                for n, line in enumerate(reslist):
                    if line.startswith('+'):
                        try:
                            include_filename = line.split()[1]
                            # Detect recoursive file inclusion:
                            if include_filename in includefiles:
                                raise ValueError('*** Recoursive include files detected! ***')
                            includefiles.append(include_filename)
                            newfile = ShelXFile.read_nested_file_to_list(os.path.join(os.path.dirname(resfile),
                                                                                      include_filename))
                            if newfile:
                                for num, l in enumerate(newfile):
                                    lnum = n + 1 + num
                                    # '+filename' include files are not copied to res file,
                                    #  so I have to delete these lines on write.
                                    # '++filename' copies them to the .res file where appropriate
                                    if l.startswith('+') and l[:2] != '++':
                                        self.delete_on_write.update([lnum])
                                    reslist.insert(lnum, l)
                                continue
                        except IndexError:
                            if DEBUG:
                                print('*** CANNOT READ INCLUDE FILE {} ***'.format(line))
                            #del reslist[n]
        except (IOError) as e:
            print(e)
            print('*** CANNOT READ FILE {} ***'.format(resfile))
        return reslist

    @staticmethod
    def read_nested_file_to_list(resfile: str) -> list:
        """
        Read in shelx file and returns a list without line endings.
        :param resfile: The path to a SHLEL .res or .ins file.
        """
        reslist = []
        try:
            with open(os.path.abspath(resfile), 'r') as f:
                reslist = f.read().splitlines(keepends=False)
        except (IOError) as e:
            if DEBUG:
                print(e)
                print('*** CANNOT OPEN NESTED INPUT FILE {} ***'.format(resfile))
            return reslist
        return reslist

    def reload(self):
        """
        Reloads the shelx file and parses it again.
        """
        if DEBUG:
            print('loading file:', self.resfile)
        self.__init__(self.resfile)

    def refine(self, cycles: int = 0) -> bool:
        bc = 0
        acta = ''
        bc = int(self.cycles.cycles)
        try:
            self._reslist[self.index_of(self.acta)] = ' '
        except ValueError:
            pass
        # THis does not work:
        #acta = self.acta.remove_acta_card()
        self.cycles.cycles = cycles
        filen, _ = os.path.splitext(self.resfile)
        self.write_shelx_file(filen + '.ins')
        #shutil.copyfile(filen+'.res', filen+'.ins')
        ref = ShelxlRefine(self, self.resfile)
        ref.run_shelxl()
        self.reload()
        # Does not work:
        #self.restore_acta_card(acta)
        self.cycles.cycles = bc
        self.write_shelx_file(filen + '.res')
        return True

    def append_card(self, obj, card, line_num):
        """
        Appends SHELX card to an object list, e.g. self.restraints and
        assigns the line_num in reslist with the card instance.
        """
        obj.append(card)
        self._reslist[line_num] = card
        return card

    def assign_card(self, card, line_num):
        self._reslist[line_num] = card
        return card

    @staticmethod
    def is_atom(atomline: str) -> bool:
        """
        Returns True is line contains an atom.
        atomline:  'O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...'
        >>> ShelXFile.is_atom(atomline = 'O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        True
        >>> ShelXFile.is_atom(atomline = 'O1    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        False
        >>> ShelXFile.is_atom(atomline = 'O1  4  0.120080    0.494426  11.00000   0.01445 ...')
        True
        >>> ShelXFile.is_atom("AFIX 123")
        False
        >>> ShelXFile.is_atom("AFIX")
        False
        >>> ShelXFile.is_atom('O1    3    0.120080   0.336659   0.494426')
        True
        """
        # no empty line, not in cards and not space at start:
        if atomline[:4].upper() not in SHX_CARDS:  # exclude all non-atom cards
            # Too few parameter for an atom:
            if len(atomline.split()) < 5:
                return False
            # means sfac number is missing:
            if '.' in atomline.split()[1]:
                return False
            return True
        else:
            return False

    def elem2sfac(self, atom_type: str) -> int:
        """
        returns an sfac-number for the element given in "atom_type"
        >>> shx = ShelXFile('../tests/p21c.res')
        >>> shx.elem2sfac('O')
        3
        >>> shx.elem2sfac('c')
        1
        >>> shx.elem2sfac('Ar')

        """
        for num, element in enumerate(self.sfac_table, 1):
            if atom_type.upper() == element.upper():
                return num  # return sfac number

    def sfac2elem(self, sfacnum: int) -> str:
        """
        returns an element and needs an sfac-number
        :param sfacnum: string like '2'
        >>> shx = ShelXFile('../tests/p21c.res')
        >>> shx.sfac2elem(1)
        'C'
        >>> shx.sfac2elem(2)
        'H'
        >>> shx.sfac2elem(3)
        'O'
        >>> shx.sfac2elem(5)
        'Al'
        >>> shx.sfac2elem(8)
        ''
        >>> shx.sfac2elem(0)
        ''
        """
        try:
            elem = self.sfac_table[int(sfacnum)]
        except IndexError:
            return ''
        return elem

    def add_line(self, linenum: int, obj):
        """
        Adds a new SHELX card to the reslist after linenum.
        """
        self._reslist.insert(linenum + 1, obj)

    def replace_line(self, obj, new_line: str):
        """
        Replaces a single line in the res file with new_line.
        """
        self._reslist[self.index_of(obj)] = new_line

    def index_of(self, obj):
        return self._reslist.index(obj)

    def insert_frag_fend_entry(self, dbatoms: list, cell: list):
        """
        Inserts the FRAG ... FEND entry in the res file.
        :param dbatoms:   list of atoms in the database entry
        :param cell:  string with "FRAG 17 cell" from the database entry
        """
        dblist = []
        for line in dbatoms:
            dblist.append("{:4} {:<4} {:>8}  {:>8}  {:>8}".format(*line))
        dblines = ' The following is from DSR:\n'
        dblines = dblines + 'FRAG 17 {} {} {} {} {} {}'.format(*cell) + '\n'
        dblines = dblines + '\n'.join(dblist)
        dblines = dblines + '\nFEND\n'
        # insert the db entry right after FVAR
        self.add_line(self.fvars.position, dblines)

    def _get_residuals(self, spline, line):
        if ShelXFile._r1_regex.match(line):
            self.R1 = float(spline[3])
            try:
                self.R1 = float(spline[3])
            except IndexError:
                if DEBUG:
                    raise
                pass
            try:
                self.data = int(spline[-2])
            except IndexError:
                if DEBUG:
                    raise
                pass
        if ShelXFile._wr2_regex.match(line):
            try:
                self.wR2 = float(spline[3].split(",")[0])
            except IndexError:
                if DEBUG:
                    raise
                pass
        if ShelXFile._parameters_regex.match(line):
            try:
                self.parameters = int(spline[1])
                if self.data and self.parameters:
                    self.dat_to_param = self.data / self.parameters
            except IndexError:
                if DEBUG:
                    raise
                pass
            try:
                self.num_restraints = int(spline[-2])
            except IndexError:
                if DEBUG:
                    raise
                pass



if __name__ == "__main__":
    def runall():
        """
        >>> file = r'p21c.res'
        >>> try:
        >>>     shx = ShelXFile(file)
        >>> except Exception:
        >>>    raise
        """
        pass

    """
    def get_commands():
        url = "http://shelx.uni-goettingen.de/shelxl_html.php"
        response = urlopen('{}/version.txt'.format(url))
        html = response.read().decode('UTF-8')
        #res = BeautifulSoup(html, "html5lib")
        tags = res.findAll("p", {"class": 'instr'})
        for l in tags:
            if l:
                print(str(l).split(">")[1].split("<")[0])
    """
    #get_commands()
    #sys.exit()

    file = r'tests/p21c.res'
    try:
        shx = ShelXFile(file)
    except Exception:
        raise

    print(shx.atoms.distance('Ga1', 'Al1'))
    print(shx.hklf)
    print(shx.sfac_table.is_exp(shx.sfac_table[1]))
    print(shx.sfac_table)
    shx.sfac_table.add_element('Zn')
    print(shx.unit)
    #print(shx.sfac_table.remove_element('In'))
    print(shx.sfac_table)
    print(shx.unit)
    shx.cycles.set_refine_cycles(33)
    shx.write_shelx_file(r'./test.ins')
    print('\n\n')
    print(shx.hklf)
    print('######################')
    sys.exit()
    # for x in shx.atoms:
    #    print(x)
    # shx.reload()
    # for x in shx.restraints:
    #    print(x)
    # for x in shx.rem:
    #    print(x)
    # print(shx.size)
    # for x in shx.sump:
    #    print(x)
    # print(float(shx.temp)+273.15)
    # print(shx.atoms.atoms_in_class('CCF3'))
    #sys.exit()
    files = walkdir(r'/Users/daniel', 'res')
    print('finished')
    for f in files:
        if "dsrsaves" in str(f) or ".olex" in str(f) or 'ED' in str(f) or 'shelXlesaves' in str(f):
            continue
        #path = f.parent
        #file = f.name
        #print(path.joinpath(file))
        #id = id_generator(size=4)
        #copy(str(f), Path(r"d:/Github/testresfiles/").joinpath(id+file))
        #print('copied', str(f.name))
        print(f)
        shx = ShelXlFile(f, debug=False)
        # print(len(shx.atoms), f)
