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
from __future__ import print_function

import os
import re
import sys

from shelxfile.cards import ACTA, FVAR, FVARs, REM, BOND, Atoms, Atom, Restraints, DEFS, NCSY, ISOR, FLAT, \
    BUMP, DFIX, DANG, SADI, SAME, RIGU, SIMU, DELU, CHIV, EADP, EXYZ, DAMP, HFIX, HKLF, SUMP, SYMM, LSCycles, \
    SFACTable, UNIT, BASF, TWIN, WGHT
from shelxfile.misc import DEBUG, ParseOrderError, ParseNumError, ParseUnknownParam, \
    split_fvar_and_parameter, flatten, time_this_method, multiline_test, dsr_regex

ver = sys.version_info
if ver < (3, 4, 0):
    print("You need at least Python 3.4 to run this program!")
    sys.exit()
import textwrap
from math import radians, cos, sin, sqrt

from shelxfile.dsrmath import Matrix

"""
TODO:

- Atoms.add_atom(position=None) method. default for position is after FVAR table
- fit fragment without shelxl
------------------------------------------------------------
- make all __repr__() unwrapped strings and wrap lines during write_res_file()
- add remove_hydrogen_atoms(atom) method.
- shx.remove_all_H([list of atoms], or all)
- check if atoms in restraints are also in structure
- deleting atoms should also remove them from restraints  
- shx.update_weight
- shx.weight_difference
- bond list
- shx.atoms.angle(at1, at2, at3)
- shx.atoms.tors(at1, at2, at3, at4)
- shx.atom.change_type('xx')
- restraints involved with an atom should also be part of the atoms properties
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


class ShelXlFile():
    """
    Class for data from a SHELXL res file. Includes Atoms, cards and unit cell.
    """
    delete_on_write = None
    atoms = None
    sump = []
    r1_regex = re.compile(r'^REM\s+R1\s+=', re.IGNORECASE)
    wr2_regex = re.compile(r'^REM\s+wR2\s+=', re.IGNORECASE)
    parameters_regex = re.compile(r'REM\s+\d+\s+parameters\s+refined', re.IGNORECASE)
    fvars = None
    sfac_table = None
    _reslist = None
    restraints = None
    dsrlines = None
    symmcards = None

    def __init__(self: 'ShelXlFile', resfile: str):
        """
        Reads the shelx file and extracts information.
        #TODO: Exchange reslist items with SHELX objects like SFAC or ATOM. Add Reslist() class.
        #TODO: line number of objects are then inside each object self.sxh.index(self)
        :param resfile: file path
        """
        self.shelx_max_line_length = 79  # maximum character lenth per line in SHELXL
        self.nohkl = False
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.V = None, None, None, None, None, None, None
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
        self.sump = []
        self.end = False
        self.maxsof = 1.0
        self.commands = []
        self.size = {}
        self.htab = []
        self.shel = []
        self.mpla = []
        self.rtab = []
        self.omit = []
        self.hklf = None
        self.grid = []
        self.free = []
        self.titl = ""
        self.exti = 0
        self.eqiv = []
        self.disp = []
        self.conn = []
        self.conv = []
        self.bind = []
        self.ansr = 0.001
        self.bloc = []
        self.cell = []
        self.dsrlines = []
        self.dsrline_nums = []
        self.symmcards = []
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
        part = False
        partnum = 0
        resi = False
        residict = {'class': '', 'number': 0, 'ID': ''}
        sof = 0
        afix = False
        afixnum = 0
        fvarnum = 1
        resinull = re.compile(r'^RESI\s+0')
        partnull = re.compile(r'^PART\s+0')
        afixnull = re.compile(r'^AFIX\s+0')
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
                list_of_lines.append(line_num + wrapindex)  # list containing the lines of a multiline command
            # The current line splitted:
            spline = line.split('!')[0].split()  # Ignore comments with "!", see how this performes
            # The current line as string:
            line = line.upper().split('!')[0]  # Ignore comments with "!", see how this performes
            # RESI class[ ] number[0] alias
            if resinull.match(line) and resi:  # RESI 0
                resi = False
                residict['number'] = 0
                residict['class'] = ''
                residict['ID'] = ''
                continue
            if resinull.match(line) and not resi:
                # A second RESI 0
                continue
            if line.startswith(('END', 'HKLF', 'RESI')) and resi:
                self._reslist.insert(line_num, "RESI 0")
                if DEBUG:
                    print('RESI in line {} was not closed'.format(line_num + 1))
                resi = False
                continue
            if line.startswith('RESI') and not resinull.match(line):
                resi = True
                residict = self.get_resi_definition_dict(spline[1:])
                continue
            # Now collect the part:
            # PART n sof
            if partnull.match(line) and part:  # PART 0
                part = False
                partnum = 0
                sof = 0
                continue
            if partnull.match(line) and not part:
                # A second PART 0
                continue
            if line.startswith(('END', 'HKLF', 'PART')) and part:
                self._reslist.insert(line_num, "PART 0")
                if DEBUG:
                    print('PART in line {} was not closed'.format(line_num + 1))
                part = False
                continue
            if line.startswith('PART') and not partnull.match(line):
                part = True
                partnum, sof = self.get_partnumber(line)
                continue
            # Now collect the AFIXes:
            # AFIX mn d[#] sof[11] U[10.08]
            if afixnull.match(line) and afix:  # AFIX 0
                afix = False
                afixnum = 0
                continue
            if afixnull.match(line) and not afix:
                # A second AFIX 0
                continue
            if line.startswith(('END', 'HKLF', 'AFIX')) and afix:
                self._reslist.insert(line_num, "AFIX 0")
                if DEBUG:
                    print('AFIX in line {} was not closed'.format(line_num + 1))
                afix = False
                continue
            elif line.startswith('AFIX') and not afixnull.match(line):
                # TODO: if afixnum > 17x: compare atoms count in frag and afix
                afix = True
                # TODO: only afixnum and sof are used at the moment. Use also u and d for AFIX():
                afixnum, d, sof, u = self.get_afix_numbers(spline, line_num)
                continue
            elif self.is_atom(line):
                # A SHELXL atom:
                # F9    4    0.395366   0.177026   0.601546  21.00000   0.03231  ( 0.03248 =
                #            0.03649  -0.00522  -0.01212   0.00157 )
                a = Atom(self, spline, list_of_lines, line_num, part=partnum, afix=afixnum, residict=residict, sof=sof)
                self.append_card(self.atoms, a, line_num)
                continue
            elif self.end:
                # Prevents e.g. parsing of second WGHT after END:
                continue
            elif line[:4] == 'SADI':
                # SADI s[0.02] pairs of atoms
                # or SADI
                if len(spline) == 1:
                    self.global_sadi = line_num
                self.append_card(self.restraints, SADI(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'DFIX':
                # DFIX d s[0.02] atom pairs
                self.append_card(self.restraints, DFIX(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'SIMU':
                # SIMU s[0.04] st[0.08] dmax[2.0] atomnames
                self.append_card(self.restraints, SIMU(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'DELU':
                # DELU s1[0.01] s2[0.01] atomnames
                self.append_card(self.restraints, DELU(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'RIGU':
                # RIGU s1[0.004] s2[0.004] atomnames
                self.append_card(self.restraints, RIGU(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'BASF':
                # BASF scale factors
                self.append_card(self.restraints, BASF(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'HFIX':
                # HFIX mn U[#] d[#] atomnames
                self.append_card(self.hfixes, HFIX(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'DANG':
                # DANG d s[0.04] atom pairs
                self.append_card(self.restraints, DANG(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'EADP':
                self.append_card(self.restraints, EADP(spline, list_of_lines), line_num)
                continue
            elif line[:3] == 'REM':
                if dsr_regex.match(line):
                    self.dsrlines.append(" ".join(spline))
                    self.dsrline_nums.extend(list_of_lines)
                self.append_card(self.rem, REM(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'AFIX':
                pass
            elif line[:4] == 'CELL':
                # CELL λ a b c α β γ
                if not lastcard == 'TITL':
                    if DEBUG:
                        print('TITL is missing.')
                    # raise ParseOrderError
                if len(spline) >= 8:
                    self.cell = [float(x) for x in spline[2:8]]
                    self.a, self.b, self.c, self.alpha, self.beta, self.gamma = self.cell
                    self.V = self.vol_unitcell(self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
                    # self.A = self.orthogonal_matrix()
                self.wavelen = float(spline[1])
                lastcard = 'CELL'
                continue
            elif line[:4] == "ZERR":
                # ZERR Z esd(a) esd(b) esd(c) esd(α) esd(β) esd(γ)
                if not lastcard == 'CELL':
                    if DEBUG:
                        print('*** Invalid SHELX file!')
                    raise ParseOrderError
                if not self.cell:
                    raise ParseOrderError('*** Cell parameters missing! ***')
                if len(spline) >= 8:
                    self.Z = spline[1]
                    self.zerr = [float(x) for x in spline[2:8]]
                lastcard = 'ZERR'
                continue
            elif line[:4] == "SYMM":
                # SYMM symmetry operation
                #  Being more greedy, because many files do this wrong:
                # if not lastcard == 'ZERR':
                #    raise ParseOrderError
                # if not self.zerr:
                #    raise ParseOrderError
                self.append_card(self.symmcards, SYMM(spline, list_of_lines), line_num)
                lastcard = 'SYMM'
                continue
            elif line[:4] == 'SFAC':
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
                lastcard = 'SFAC'
                continue
            elif line[:4] == 'UNIT':
                # UNIT n1 n2 ...
                # Number of atoms of each type in the unit-cell, in SFAC order.
                if not lastcard == 'SFAC':
                    raise ParseOrderError
                if self.sfac_table:
                    try:
                        self.unit = self.assign_card(UNIT(spline, list_of_lines), line_num)
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
            elif line[:4] == "LATT":
                # LATT N[1]
                # 1=P, 2=I, 3=rhombohedral obverse on hexagonal axes, 4=F, 5=A, 6=B, 7=C.
                # negative is non-centrosymmetric
                self.latt = spline[1]
                if not lastcard == 'ZERR':
                    if DEBUG:
                        print('*** ZERR instruction is missing! ***')
                    # raise ParseOrderError
                continue
            elif line[:4] in ['L.S.', 'CGLS']:
                # CGLS nls[0] nrf[0] nextra[0]
                # L.S. nls[0] nrf[0] nextra[0]
                self.cycles = self.assign_card(LSCycles(self, spline, line_num), line_num)
                continue
            elif line[:4] == "LIST":
                # LIST m[#] mult[1] (mult is for list 4 only)
                self.list = int(spline[1])
                continue
            elif line[:4] == "FVAR":
                # FVAR osf[1] free variables
                # TODO: assign value
                for fvvalue in spline[1:]:
                    fvarnum += 1
                    self.append_card(self.fvars, FVAR(fvarnum, float(fvvalue)), line_num)
                    if self.fvars not in self._reslist:
                        self._reslist[line_num] = self.fvars
                    else:
                        self.delete_on_write.update([line_num])
            elif line[:4] == 'ANIS':
                # ANIS n or ANIS names
                # Must be before Atom(), to know which atom is anis.
                self.anis = spline
                continue
            elif line[:4] == 'WGHT':
                # WGHT a[0.1] b[0] c[0] d[0] e[0] f[.33333]
                self.wght = self.assign_card(WGHT(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'ACTA':
                # ACTA 2θfull[#] -> optional parameter NOHKL
                self.acta = self.append_card(self.commands, ACTA(self, spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'DAMP':
                # DAMP damp[0.7] limse[15]
                self.damp = self.append_card(self.commands, DAMP(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'ABIN':
                # ABIN n1 n2   ->   Reads h, k, l, A and B from the file name.fab
                self.abin = [float(x) for x in spline[1:]]
                continue
            elif line[:4] == 'ANSC':
                # ANSC six coefficients
                if len(spline) == 7:
                    self.ansc = [float(x) for x in spline[:1]]
                continue
            elif line[:4] == 'ANSR':
                # ANSR anres[0.001]
                if len(spline) == 2:
                    self.ansr = float(spline[1])
                continue
            elif line[:4] == 'BIND':
                # BIND atom1 atom2
                if len(spline) == 3:
                    self.bind.append(spline[1:])
                continue
            elif line[:4] == 'BLOC':
                # BLOC n1 n2 atomnames
                # TODO: Make class that resolves atomnames and cycles
                self.bloc.append(spline[1:])
                continue
            elif line[:4] == 'BOND':
                # BOND atomnames
                self.append_card(self.commands, BOND(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'BUMP':
                # BUMP s [0.02]
                self.append_card(self.restraints, BUMP(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'CHIV':
                # CHIV V[0] s[0.1] atomnames
                self.append_card(self.restraints, CHIV(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'CONF':
                # CONF atomnames max_d[1.9] max_a[170]
                self.conv.append(spline[1:])
                continue
            elif line[:4] == 'CONN':
                # CONN bmax[12] r[#] atomnames or CONN bmax[12]
                # bonded are d < (r1 + r2 + 0.5) Å
                self.conn.append(spline[1:])
                continue
            elif line[:4] == 'DEFS':
                # DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
                self.defs = DEFS(spline, list_of_lines)
                self.assign_card(self.defs, line_num)
                continue
            elif line[:4] == 'DISP':
                # DISP E f' f"[#] mu[#]
                if not lastcard == 'SFAC':
                    raise ParseOrderError
                self.disp.append(spline[1:])
                continue
            elif line[:4] == 'EQIV':
                # EQIV $n symmetry operation
                if len(spline) > 1:
                    if spline[1].startswith('$'):
                        self.eqiv.append(spline[1:])
                continue
            elif line[:4] == 'EXTI':
                # EXTI x[0]
                self.exti = float(spline[1])
                continue
            elif line[:4] == 'EXYZ':
                # EXYZ atomnames
                self.append_card(self.restraints, EXYZ(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'FRAG':
                # FRAG code[17] a[1] b[1] c[1] α[90] β[90] γ[90]
                if len(spline) == 8:
                    self.frag = spline[1:]
                continue
            elif line[:4] == 'FEND':
                # FEND (must follow FRAG)
                if not self.frag:
                    raise ParseOrderError
                self.frag = None  # Turns frag mode off.
                continue
            elif line[:4] == 'FLAT':
                # FLAT s[0.1] four or more atoms
                self.append_card(self.restraints, FLAT(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'FREE':
                # FREE atom1 atom2
                self.free = spline[1:]
                continue
            elif line[:4] == 'GRID':
                # GRID sl[#] sa[#] sd[#] dl[#] da[#] dd[#]
                self.grid = spline[1:]
                continue
            elif line[:4] == 'HKLF':
                # HKLF N[0] S[1] r11...r33[1 0 0 0 1 0 0 0 1] sm[1] m[0]
                self.hklf = HKLF(spline, list_of_lines)
                self.assign_card(self.hklf, line_num)
                continue
            elif line.startswith('END'):
                # END (after HKLF or ends an include file)
                # TODO: run sanity checks after END like checking if EXYZ and
                # anisotropy fit togeter
                self.end = True
                continue
            elif line[:4] == 'HTAB':
                # HTAB dh[2.0]  or  HTAB donor-atom acceptor-atom
                self.htab = spline[1:]
                continue
            elif line[:4] == 'ISOR':
                # ISOR s[0.1] st[0.2] atomnames
                self.append_card(self.restraints, ISOR(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'LAUE':
                # LAUE E
                # I completely do not understand the LAUE instruction description in the manual!
                continue
            elif line[:4] == 'MERG':
                # MERG n[2]
                self.merg = spline[1:]
                continue
            elif line[:4] == 'MORE':
                # MORE m[1]
                self.more = spline[1]
                continue
            elif line[:4] == 'FMAP':
                # FMAP code[2] axis[#] nl[53]
                self.fmap = spline[1:]
                continue
            elif line[:4] == 'MOVE':
                # MOVE dx[0] dy[0] dz[0] sign[1]
                self.move = spline[1:]
                continue
            elif line[:4] == 'MPLA':
                # MPLA na atomnames
                self.mpla.append(spline[1:])
                continue
            elif line[:4] == 'NCSY':
                # NCSY DN sd[0.1] su[0.05] atoms
                self.append_card(self.restraints, NCSY(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'NEUT':
                # NEUT
                if not lastcard == 'SYMM':
                    raise ParseOrderError
                continue
            elif line[:4] == 'OMIT':
                # OMIT atomnames  or  OMIT s[-2] 2θ(lim)[180]  or  OMIT h k l
                self.omit.append(spline[1:])
                continue
            elif line[:4] == 'PLAN':
                # PLAN npeaks[20] d1[#] d2[#]
                self.plan = spline[1:]
                continue
            elif line[:4] == 'PRIG':
                # PRIG p[#]
                self.prig = spline[1:]
                continue
            elif line[:4] == 'RTAB':
                # RTAB codename atomnames  -->  codename: e.g. 'omeg' gets tabualted in the lst
                self.rtab.append(spline[1:])
                continue
            elif line[:4] == 'SAME':
                # SAME s1[0.02] s2[0.04] atomnames
                self.append_card(self.restraints, SAME(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'SHEL':
                # SHEL lowres[infinite] highres[0]
                self.shel = spline[1:]
                continue
            elif line[:4] == 'SIZE':
                # SIZE dx dy dz
                if len(spline) == 4:
                    self.size['x'] = spline[1]
                    self.size['y'] = spline[2]
                    self.size['z'] = spline[3]
                continue
            elif line[:4] == 'SPEC':
                # SPEC del[0.2]
                # This implementation is not enough, but more is maybe only needed for a
                # refinement program:
                if len(spline) > 1:
                    self.spec = spline[1]
                continue
            elif line[:4] == 'STIR':
                # STIR sres step[0.01]   -> stepwise improvement in the resolution sres
                self.stir = spline[1:]
                continue
            elif line[:4] == 'SUMP':
                # SUMP c sigma c1 m1 c2 m2 ...
                self.append_card(self.sump, SUMP(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'SWAT':
                # SWAT g[0] U[2]
                self.swat = spline[1:]
                continue
            elif line[:4] == 'TEMP':
                # TEMP T[20]  -> in Celsius
                self.temp = spline[1]
                continue
            elif line[:4] == 'TWIN':
                # TWIN 3x3 matrix [-1 0 0 0 -1 0 0 0 -1] N[2]
                self.twin = TWIN(spline, list_of_lines)
                self.assign_card(self.twin, line_num)
                continue
            elif line[:4] == 'TWST':
                # TWST N[0] (N[1] after SHELXL-2018/3)
                if len(spline) > 1:
                    self.twst = spline[1]
                continue
            elif line[:4] == 'WIGL':
                # WIGL del[0.2] dU[0.2]
                if len(spline) == 1:
                    self.wigl = True
                if len(spline) > 1:
                    self.wigl = spline[1:]
                continue
            elif line[:4] == 'WPDB':
                # WPDB n[1]
                self.wpdb = spline[1:]
                continue
            elif line[:4] == 'XNPD':
                # XNPD Umin[-0.001]
                self.xnpd = spline[1:]
                continue
            elif line[:4] == 'BEDE':
                # Later...
                continue
            elif line[:4] == 'LONE':
                # Later...
                continue
            elif ShelXlFile.r1_regex.match(line):
                try:
                    self.R1 = float(spline[3])
                except IndexError:
                    pass
                try:
                    self.data = float(spline[-2])
                except IndexError:
                    pass
                continue
            elif ShelXlFile.wr2_regex.match(line):
                try:
                    self.wR2 = float(spline[3])
                except IndexError:
                    pass
                continue
            elif ShelXlFile.parameters_regex.match(line):
                try:
                    self.parameters = float(spline[1])
                    if self.data and self.parameters:
                        self.dat_to_param = self.data / self.parameters
                except IndexError:
                    pass
                try:
                    self.num_restraints = float(spline[-2])
                except IndexError:
                    pass
                continue
            elif line[:4] == 'MOLE':
                # print('*** MOLE is deprecated! Do not use it! ***')
                pass
            elif line[:4] == 'HOPE':
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

    def restore_acta_card(self, acta: str):
        """
        Place ACTA after UNIT
        """
        self.add_line(self.unit.line_numbers[-1] + 1, acta)

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
        return Matrix([[self.a, self.b * cos(self.gamma), self.c * cos(self.beta)],
                       [0, self.b * sin(self.gamma),
                        (self.c * (cos(self.alpha) - cos(self.beta) * cos(self.gamma)) / sin(self.gamma))],
                       [0, 0, self.V / (self.a * self.b * sin(self.gamma))]])

    @staticmethod
    def vol_unitcell(a, b, c, al, be, ga) -> float:
        """
        calculates the volume of a unit cell
        >>> v = ShelXlFile.vol_unitcell(2, 2, 2, 90, 90, 90)
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
            if line == '' and self._reslist[num + 1] == '':
                continue
            line = self.wrap_line(line)
            resl.append(line)
        return "\n".join(resl)

    def wrap_line(self, line: str) -> str:
        # Generally, all Shelx opbjects have no line wrap. I do this now:
        line = textwrap.wrap(line, 78, subsequent_indent='  ', drop_whitespace=False, replace_whitespace=False)
        if len(line) > 1:
            newline = []
            for n, l in enumerate(line):
                if n < len(line) - 1:
                    l += ' =\n'
                newline.append(l)
            line = ' '.join(newline)
        else:
            line = ''.join(line)
        return line

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
                            newfile = ShelXlFile.read_nested_file_to_list(os.path.join(os.path.dirname(resfile),
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
                line = self.wrap_line(str(line))
                f.write(line + '\n')
        if verbose or DEBUG:
            print('File successfully written to {}'.format(os.path.abspath(filename)))
            return True
        return True

    def append_card(self, obj, card, line_num):
        obj.append(card)
        self._reslist[line_num] = card

    def assign_card(self, card, line_num):
        self._reslist[line_num] = card
        return card

    @staticmethod
    def get_resi_definition_dict(resi: list) -> dict:
        """
        Returns the residue number and class of a string like 'RESI TOL 1'
        or 'RESI 1 TOL'

        Residue names may now begin with a digit.
        They must however contain at least one letter

        Allowed residue numbers is now from -999 to 9999 (2017/1)

        TODO: support alias

        :param resi: ['number', 'class']
        :type resi: list or string

        >>> sorted(list(ShelXlFile.get_resi_definition_dict('RESI 1 TOL'.split()[1:])))
        ['ID', 'class', 'number']
        >>> sorted(ShelXlFile.get_resi_definition_dict('RESI 1 TOL'.split()[1:]).items())
        [('ID', None), ('class', 'TOL'), ('number', 1)]
        >>> ShelXlFile.get_resi_definition_dict('RESI 1 TOL'.split()[1:])
        {'class': 'TOL', 'number': 1, 'ID': None}
        >>> ShelXlFile.get_resi_definition_dict('RESI A:100 TOL'.split()[1:])
        {'class': 'TOL', 'number': 100, 'ID': 'A'}
        >>> ShelXlFile.get_resi_definition_dict('RESI -10 TOL'.split()[1:])
        {'class': 'TOL', 'number': -10, 'ID': None}
        >>> ShelXlFile.get_resi_definition_dict('RESI b:-10 TOL'.split()[1:])
        {'class': 'TOL', 'number': -10, 'ID': 'b'}
        """
        resi_dict = {'class': None, 'number': None, 'ID': None}
        for x in resi:
            if re.search('[a-zA-Z]', x):
                if ':' in x:
                    # contains :, must be a chain-id+number
                    resi_dict['ID'], resi_dict['number'] = x.split(':')[0], int(x.split(':')[1])
                else:
                    # contains letters, must be a name
                    resi_dict['class'] = x
            else:
                # everything else can only be a number
                resi_dict['number'] = int(x)
        return resi_dict

    @staticmethod
    def get_partnumber(partstring: str) -> (int, float):
        """
        get the part number from a string like PART 1 oder PART 2 -21

        PART n sof
        partstring: string like 'PART 2 -21'

        >>> ShelXlFile.get_partnumber(partstring='PART 2 -21')
        (2, -21.0)
        """
        part = partstring.upper().split()
        sof = 0
        try:
            partnum = int(part[1])
        except(ValueError, IndexError):
            if DEBUG:
                print('*** Wrong PART definition found! Check your PART instructions ***')
            partnum = 0
        if len(part) > 2:
            sof = float(part[2])
        return partnum, sof

    @staticmethod
    def get_afix_numbers(spline: list, line_num: int) -> ((int, float, float, float), [int, float, float, float]):
        """
        Returns a tuple with the AFIX instructions. afixnum, d, sof, u

        AFIX mn d[#] sof[11] U[10.08]

        >>> ShelXlFile.get_afix_numbers(['AFIX', 137], 1)
        (137, 0.0, 11, 10.08)
        >>> ShelXlFile.get_afix_numbers(['AFIX', '137b'], 1)
        *** Wrong AFIX definition in line 1. Check your AFIX instructions ***
        (0, 0.0, 11, 10.08)
        >>> ShelXlFile.get_afix_numbers(['AFIX', 13, 1.234], 1)
        (13, 1.234, 11, 10.08)
        >>> ShelXlFile.get_afix_numbers(['AFIX', 13, 1.234, -21, 10.05], 1)
        (13, 1.234, -21.0, 10.05)
        """
        d = 0.0
        sof = 11
        u = 10.08
        try:
            afixnum = int(spline[1])
        except(ValueError, IndexError):
            if DEBUG:
                print('*** Wrong AFIX definition in line {}. Check your AFIX instructions ***'.format(line_num))
            afixnum = 0
        if len(spline) > 2:
            d = float(spline[2])
        if len(spline) > 3:
            sof = float(spline[3])
        if len(spline) > 4:
            u = float(spline[4])
        return afixnum, d, sof, u

    @staticmethod
    def is_atom(atomline: str) -> bool:
        """
        Returns True is line contains an atom.
        atomline:  'O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...'
        >>> ShelXlFile.is_atom(atomline = 'O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        True
        >>> ShelXlFile.is_atom(atomline = 'O1    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        False
        >>> ShelXlFile.is_atom(atomline = 'O1  4  0.120080    0.494426  11.00000   0.01445 ...')
        True
        >>> ShelXlFile.is_atom("AFIX 123")
        False
        >>> ShelXlFile.is_atom("AFIX")
        False
        >>> ShelXlFile.is_atom('O1    3    0.120080   0.336659   0.494426')
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
        >>> shx = ShelXlFile('p21c.res')
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
        >>> shx = ShelXlFile('./p21c.res')
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
        print(self.fvars.line_number)
        self.add_line(self.fvars.line_number, dblines)



if __name__ == "__main__":
    def runall():
        """
        >>> file = r'p21c.res'
        >>> try:
        >>>     shx = ShelXlFile(file)
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
        shx = ShelXlFile(file)
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