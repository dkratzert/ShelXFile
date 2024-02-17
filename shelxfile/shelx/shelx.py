# -*- encoding: utf-8 -*-
# möpß
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
__doc__ = """
This is a full implementation of the SHELXL file syntax. Additionally it is able to edit SHELX properties with Python.
The implementation is Python3-only and supports SHELXL after 2017 (You should not use old versions anyway).

The parser is quiet about the most errors unless you enable debug or verbose during initialization. 
The parser will try to read the SHELX file even if it has syntax errors, but if for example, the SFAC and UNIT 
instruction is not consistent it will fail. 
"""

import re
import sys
from contextlib import suppress
from pathlib import Path
from typing import Union, List, Optional

from shelxfile.atoms.atom import Atom
from shelxfile.atoms.atoms import Atoms
from shelxfile.cif.cif_write import CifFile
from shelxfile.misc.dsrmath import Array
from shelxfile.misc.elements import weight_from_symbol
# noinspection PyUnresolvedReferences
from shelxfile.misc.misc import ParseOrderError, ParseNumError, ParseUnknownParam, \
    multiline_test, dsr_regex, wrap_line, ParseSyntaxError
from shelxfile.refine.refine import ShelxlRefine
from shelxfile.shelx.cards import ACTA, FVAR, FVARs, REM, BOND, Restraints, DEFS, NCSY, ISOR, FLAT, \
    BUMP, DFIX, DANG, SADI, SAME, RIGU, SIMU, DELU, CHIV, EADP, EXYZ, DAMP, HFIX, HKLF, SUMP, SYMM, LSCycles, \
    SFACTable, UNIT, BASF, TWIN, WGHT, BLOC, SymmCards, CONN, CONF, BIND, DISP, GRID, HTAB, MERG, FRAG, FREE, FMAP, \
    MOVE, PLAN, PRIG, RTAB, SHEL, SIZE, SPEC, STIR, TWST, WIGL, WPDB, XNPD, ZERR, CELL, LATT, MORE, MPLA, AFIX, PART, \
    RESI, ABIN, ANIS, Residues, SWAT, Command, Restraint
from shelxfile.shelx.sdm import SDM
from shelxfile.version import VERSION

__version__ = VERSION

"""
TODO:
- Handle BEDE & LONE plus their results
- Rotate ellipsoids with kabsch
- killall.Q, killall.C 
- Q-peak printing is wrong: Q1    1   0.9828    1.1159    0.3148   11.00000  0.04      0.00  
- Write out parts and afix in grow mode
- Delete atoms (H) in AFIX -> delete entire afix group
- Is backup file and reccovery from failed refinement working?
------------------------------------------------------------
- add remove_hydrogen_atoms(atom) method.
- shx.remove_all_H([list of atoms], or all)
- bond list
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


class Shelxfile():
    """
    Class for data from a SHELXL res file. Includes Atoms, cards and unit cell.
    """
    _r1_regex = re.compile(r'^REM\s+R1\s+=', re.IGNORECASE)
    _wr2_regex = re.compile(r'^REM\s+wR2\s+=', re.IGNORECASE)
    _parameters_regex = re.compile(r'^REM\s+\d+\s+parameters\s+refined', re.IGNORECASE)
    _diff_peak_regex = re.compile(r'^REM\sHighest\sdifference', re.IGNORECASE)
    _goof_regex = re.compile(r'^REM\swR2\s=\s.*,\sGooF', re.IGNORECASE)
    _spgrp_regex = re.compile(r'^REM\s+\S+\s+in\s+\S+', re.IGNORECASE)

    def __init__(self, verbose: bool = False, debug: bool = False) -> None:
        if debug and verbose:
            raise ValueError("Either 'verbose' or 'debug' allowed, not both.")
        self.debug = False
        self.verbose = False
        if debug:
            self.debug = True
        if verbose:
            self.verbose = True
        # print(f'DEBUG: {self.debug}, VERBOSE: {self.verbose}')
        self.restraint_errors: List[str] = []
        self.temp_in_kelvin: float = 0.0
        self.shelx_max_line_length: int = 79  # maximum character lenth per line in SHELXL
        self.cell: Optional[CELL] = None
        self.ansc: List[float] = []
        self.abin: Optional[ABIN] = None
        self.acta: Optional[ACTA] = None
        self.fmap: Optional[FMAP] = None
        self.xnpd: Optional[XNPD] = None
        self.wpdb: Optional[WPDB] = None
        self.wigl: Optional[WIGL] = None
        self.temp: Union[int, float] = 20
        self.swat: Optional[SWAT] = None
        self.stir: Optional[STIR] = None
        self.spec: Optional[SPEC] = None
        self.twst: Optional[TWST] = None
        self.plan: Optional[PLAN] = None
        self.prig: Optional[PRIG] = None
        self.merg: Optional[MERG] = None
        self.more: Optional[MORE] = None
        self.move: Optional[MOVE] = None
        self.defs: Optional[DEFS] = None
        self.zerr: Optional[ZERR] = None
        self.wght: Optional[WGHT] = None
        self.frag: Optional[FRAG] = None
        self.twin: Optional[TWIN] = None
        self.basf: Optional[BASF] = None
        self.latt: Optional[LATT] = None
        self.anis: Optional[ANIS] = None
        self.damp: Optional[DAMP] = None
        self.unit: Optional[UNIT] = None
        self.size: Optional[SIZE] = None
        self.htab: Optional[HTAB] = None
        self.shel: Optional[SHEL] = None
        self.mpla: Optional[MPLA] = None
        self.hklf: Optional[HKLF] = None
        self.grid: Optional[GRID] = None
        self.conn: Optional[CONN] = None
        self.conf: Optional[CONF] = None
        self.afix: Optional[AFIX] = None
        self.rtab: List[RTAB] = []
        self.omit: List[str] = []
        self.free: List[FREE] = []
        self.eqiv: List[str] = []
        self.bonds: List[BOND] = []
        self.disp: List[BOND] = []
        self.bind: List[BIND] = []
        self.bloc: List[BLOC] = []
        self.part: PART = PART(self, ['PART', '0'])
        self.resi: RESI = RESI(self, ['RESI', '0'])
        self.residues = Residues(self)
        self.dsrlines: List[str] = []
        self.dsrline_nums: List[int] = []
        self.symmcards: SymmCards = SymmCards(self)
        self.hfixes: List[HFIX] = []
        self.sump: List[SUMP] = []
        self.wght_suggested: Optional[WGHT] = None
        self.Z: int = 1
        self.titl: str = ""
        self.exti: float = 0.0
        self.ansr: float = 0.001
        self.rem: List[REM] = []
        self.atoms: Atoms = Atoms(self)
        self.fvars: FVARs = FVARs(self)
        self.restraints: Restraints = Restraints()
        self.sfac_table: SFACTable = SFACTable(self)
        self.cycles: Optional[LSCycles] = None
        self.R1: Optional[float] = None
        self.wr2: Optional[float] = None
        self.goof: Optional[float] = None
        self.rgoof: Optional[float] = None
        self.space_group: Optional[str] = None
        self.data: Optional[int] = None
        self.parameters: Optional[int] = None
        self.dat_to_param: Optional[float] = None
        self.num_restraints: Optional[int] = None
        self.highest_peak: Optional[float] = None
        self.deepest_hole: Optional[float] = None
        self.formula_weight: Optional[float] = None
        self.end: bool = False
        self.maxsof: float = 1.0
        self.delete_on_write: set = set()
        self.wavelength: float = 0.0
        self.global_sadi: Optional[int] = None
        self.list: int = 0
        self.theta_full: float = 0.0
        self.error_line_num: int = -1  # Only used to tell the line number during an exception.
        self.resfile: Optional[Path] = None
        self._reslist: List[Union[str, Command, SFACTable, FVARs, Atom, SYMM]] = []

    def write_shelx_file(self, filename: Union[str, Path, None] = None) -> None:
        if not self._reslist:
            print('*** No file was loaded for writing. ***')
            return None
        if not filename:
            filename = self.resfile
        if isinstance(filename, str):
            filename = Path(filename)
        with open(filename, 'w') as f:
            for num, line in enumerate(self._reslist):
                if num in self.delete_on_write:
                    if self.debug:
                        # print('Deleted line {}'.format(num + 1))
                        pass
                    continue
                if line == '':  # and self._reslist[num + 1] == '':
                    continue
                # Prevent wrapping long lines with \n breaks by splitting first:
                line = "\n".join([wrap_line(x) for x in str(line).split("\n")])
                f.write(str(line) + '\n')
        if self.verbose or self.debug:
            print(f'*** File successfully written to {filename.resolve()} ***')

    def read_file(self, resfile: Union[Path, str]) -> None:
        """
        Read input from a file path.
        """
        self.__init__(debug=self.debug, verbose=self.verbose)
        if isinstance(resfile, str):
            resfile = Path(resfile)
        self.resfile = resfile.resolve()
        if self.debug:
            print(f'Resfile is: {resfile}')
        try:
            self._reslist: List = resfile.read_text().splitlines(keepends=False)
            self._test_if_file_is_valid(resfile)
        except UnicodeDecodeError:
            if self.debug or self.verbose:
                print(f'*** Unable to read file {resfile} ***')
            return
        self._find_included_files()
        self.parse_cards()

    def read_string(self, resfile_string: str):
        """
        Read input as string.
        This will not read files included with "+filename" syntax!
        """
        self.__init__(debug=self.debug, verbose=self.verbose)
        self._reslist = resfile_string.splitlines(keepends=False)
        self.parse_cards()

    def parse_cards(self) -> None:
        try:
            self._parse_cards()
        except Exception as e:
            if self.debug or self.verbose:
                self.show_line_where_error_occured(e)
                if self.debug:
                    raise
            else:
                return
        self.restraint_errors = self._assign_atoms_to_restraints()

    def _assign_atoms_to_restraints(self) -> List[str]:
        warnings = []
        for restraint in self.restraints:
            bad_atoms = []
            for restraint_atom in restraint.atoms:
                if restraint_atom in ('>', '<', '='):
                    continue
                if (restraint.residue_class or sum(restraint.residue_number) > 0) and '_' not in restraint_atom:
                    for num in restraint.residue_number:
                        self.does_atom_exist(f'{restraint_atom}_{num}', bad_atoms, f'{restraint_atom}_{num}')
                elif '_' in restraint_atom:
                    self.does_atom_exist(f'{restraint_atom}', bad_atoms, restraint_atom)
                else:
                    self.does_atom_exist(f'{restraint_atom}_{0}', bad_atoms, restraint_atom)
            if bad_atoms:
                sorted_atoms = list(set(bad_atoms))
                sorted_atoms.sort()
                warnings.append(f'*** Unknown atom{"s" if len(bad_atoms) > 0 else ""} in restraint: {restraint}, '
                                f'line {restraint.index + 1} ***')
                warnings.append(f'*** Atom list has no --> {", ".join(sorted_atoms)} ***')
                bad_atoms.clear()
            if restraint.residue_class and sum(restraint.residue_number) == 0:
                warnings.append(f"*** Restraint '{restraint}', line {restraint.index + 1}, "
                                f"has a residue class, but no residues are defined. ***")
        if self.debug or self.verbose:
            print('\n'.join(warnings))
        return warnings

    def does_atom_exist(self, atom_name: str, bad_atoms: List[str], restraint_atom: str):
        residue_number_is_wildcard = '_' in atom_name and atom_name.split('_')[-1] == '*'
        if atom_name.startswith('$'):
            return None
        elif residue_number_is_wildcard:
            for num in self.residues.residue_numbers.keys():
                residue_atom = f"{atom_name.split('_')[0]}_{num}"
                if not self.atoms.get_atom_by_name(residue_atom):
                    bad_atoms.append(residue_atom)
        else:
            if not self.atoms.get_atom_by_name(atom_name):
                bad_atoms.append(restraint_atom)

    def _test_if_file_is_valid(self, resfile: Path) -> None:
        if len(self._reslist) < 20 and (self.debug or self.verbose):
            print('*** Not a SHELXL file: {} ***'.format(resfile))
            if self.debug:
                sys.exit()

    def show_line_where_error_occured(self, e):
        try:
            print(f'Error near:\n {self._reslist[self.error_line_num]}')
        except IndexError:
            pass
        print(e)
        print(f"*** Syntax error found in file {self.resfile}, line {self.error_line_num + 1} ***")

    def _find_included_files(self) -> None:
        # Tracks the file names of included files in order to find recursive inclusion:
        includefiles = []
        for line_num, line in enumerate(self._reslist):
            if line.startswith('+'):
                try:
                    file_included_in_includefile = self._read_included_file(includefiles, line)
                    if file_included_in_includefile:
                        for line_num_includefile, l in enumerate(file_included_in_includefile):
                            reslist_position = line_num + 1 + line_num_includefile
                            # '+filename' include files are not copied to res file,
                            #  so I have to delete these lines on write.
                            # '++filename' copies them to the .res file where appropriate
                            # I leave this out, because I am not SHELXL:
                            # if l.startswith('+') and l[:2] != '++':
                            #    self.delete_on_write.update([lnum])
                            self._reslist.insert(reslist_position, l)
                        continue
                except IndexError:
                    if self.debug or self.verbose:
                        print(f'*** CANNOT READ INCLUDE FILE {line} ***')
                    # Not sure if this is a good idea: del reslist[n]

    def _read_included_file(self, includefiles: List[str], line: str):
        include_filename: Path = self.resfile.resolve().parent.joinpath(line[1:])
        # Detect recursive file inclusion:
        if include_filename.name in includefiles:
            raise ValueError('*** Recoursive include files detected! ***')
        includefiles.append(include_filename.name)
        try:
            newfile = include_filename.read_text().splitlines(keepends=False)
        except IOError as e:
            if self.debug or self.verbose:
                print(e)
                print(f'*** CANNOT OPEN NESTED INPUT FILE {include_filename} ***')
            return []
        return newfile

    def reload(self) -> None:
        """
        Reloads the shelx file and parses it again.
        """
        if self.debug or self.verbose:
            print(f'*** reloading file: {self.resfile} ***')
        self.read_file(self.resfile.resolve())

    def _parse_cards(self) -> None:
        last_nonhydrogen_atom: Optional[Atom] = None
        lastcard = ''
        fvarnum = 1
        for line_num, line in enumerate(self._reslist):
            self.error_line_num = line_num  # For exception during parsing.
            list_of_lines = [line_num]  # list of lines where a card appears, e.g. for atoms with two lines
            if line.startswith(' ') or line == '':
                continue
            wrapindex = 0
            # This while loop makes wrapped lines look like they are not wrapped. The following lines are then
            # beginning with a space character and thus are ignored. The 'lines' list holds the line nnumbers where
            # 'line' is located ([line_num]) plus the wrapped lines.
            if multiline_test(self._reslist[line_num]):
                multiline = True
            else:
                multiline = False
            while multiline:
                # Glue together the two lines wrapped with "=":
                wrapindex += 1
                line = line.rpartition('=')[0] + self._reslist[line_num + wrapindex]
                # self.delete_on_write.update([line_num + wrapindex])
                list_of_lines.append(line_num + wrapindex)  # list containing the lines of a multiline command
                # Do not activate this, otherwise, the unwrapping stops after two lines.
                if multiline_test(self._reslist[line_num + wrapindex]):
                    multiline = True
                else:
                    multiline = False
                self._reslist[line_num + wrapindex] = ''
            # The current line split:
            spline: list = line.split('!')[0].split()  # Ignore comments with "!", see how this performes
            # The current line as string:
            line = line.upper().split('!')[0]  # Ignore comments with "!", see how this performes
            word = line[:4]
            # get RESI:
            if line.startswith(('END', 'HKLF')) and self.resi:
                self.resi.num = 0
                if self.debug or self.verbose:
                    print('RESI in line {} was not closed'.format(line_num + 1))
                # Do not continue here, otherwise HKLF is not parsed
                # continue
            if word == 'RESI':
                self.resi = self._assign_card(RESI(self, spline), line_num)
                if self.resi.residue_number > 0:
                    self.residues.append(self.resi)
                continue
            # Now collect the PART:
            if line.startswith(('END', 'HKLF')) and self.part:
                self.part.n = 0
                if self.debug or self.verbose:
                    print('PART in line {} was not closed'.format(line_num + 1))
                # Do not continue here, otherwise HKLF is not parsed
                # continue
            if word == 'PART':
                self.part = self._assign_card(PART(self, spline), line_num)
                continue
            # collect AFIX:
            if line.startswith(('END', 'HKLF')) and self.afix:
                self.afix.mn = 0
                if self.debug or self.verbose:
                    print('AFIX in line {} was not closed'.format(line_num + 1))
            elif word == 'AFIX':
                self.afix = self._assign_card(AFIX(self, spline), line_num)
            elif self.is_atom(line):
                # A SHELXL atom:
                # F9    4    0.395366   0.177026   0.601546  21.00000   0.03231  ( 0.03248 =
                #            0.03649  -0.00522  -0.01212   0.00157 )
                a = Atom(self)
                if last_nonhydrogen_atom:
                    a.pivot = last_nonhydrogen_atom
                a.parse_line(spline, list_of_lines, part=self.part, afix=self.afix, resi=self.resi)
                if not a.is_hydrogen:
                    last_nonhydrogen_atom = a
                    a.pivot = None
                self._append_card(self.atoms, a, line_num)
            elif word == 'SADI':
                # SADI s[0.02] pairs of atoms
                # or SADI
                if len(spline) == 1:
                    self.global_sadi = line_num
                self._append_card(self.restraints, SADI(self, spline), line_num)
            elif word == 'DFIX':
                # DFIX d s[0.02] atom pairs
                self._append_card(self.restraints, DFIX(self, spline), line_num)
            elif word == 'SIMU':
                # SIMU s[0.04] st[0.08] dmax[2.0] atomnames
                self._append_card(self.restraints, SIMU(self, spline), line_num)
            elif word == 'DELU':
                # DELU s1[0.01] s2[0.01] atomnames
                self._append_card(self.restraints, DELU(self, spline), line_num)
            elif word == 'RIGU':
                # RIGU s1[0.004] s2[0.004] atomnames
                self._append_card(self.restraints, RIGU(self, spline), line_num)
            elif word == 'BASF':
                # BASF scale factors
                self._assign_card(BASF(self, spline), line_num)
            elif word == 'HFIX':
                # HFIX mn U[#] d[#] atomnames
                self._append_card(self.hfixes, HFIX(self, spline), line_num)
            elif word == 'DANG':
                # DANG d s[0.04] atom pairs
                self._append_card(self.restraints, DANG(self, spline), line_num)
            elif word == 'EADP':
                self._append_card(self.restraints, EADP(self, spline), line_num)
            elif line.startswith('REM'):
                if dsr_regex.match(line):
                    self.dsrlines.append(" ".join(spline))
                    self.dsrline_nums.extend(list_of_lines)
                self._append_card(self.rem, REM(self, spline), line_num)
                self._get_residuals(spline, line)
            elif word == 'CELL':
                # CELL λ a b c α β γ
                if lastcard != 'TITL' and (self.debug or self.verbose):
                    print('*** TITL is missing. ***')
                self.cell = self._assign_card(CELL(self, spline), line_num)
                self.orthogonal_matrix = self.cell.o
                self.wavelength = self.cell.wavelen
                lastcard = 'CELL'
            elif word == "ZERR":
                # ZERR Z esd(a) esd(b) esd(c) esd(α) esd(β) esd(γ)
                if lastcard != 'CELL':
                    if self.debug or self.verbose:
                        print('*** Invalid SHELX file: CELL must occur before ZERR. ***')
                    if self.debug:
                        raise ParseOrderError(debug=shx.debug, verbose=shx.verbose)
                if not self.cell:
                    raise ParseOrderError('*** Cell parameters missing! ***', debug=shx.debug, verbose=shx.verbose)
                if len(spline) >= 8:
                    self.zerr: ZERR = self._assign_card(ZERR(self, spline), line_num)
                    self.Z = self.zerr.Z
                    if self.Z < 1:
                        self.Z = 1
                        if self.verbose or self.debug:
                            print('*** Warning: Z value is zero. ***')
                lastcard = 'ZERR'
            elif word == "LATT":
                # LATT N[1]
                # 1=P, 2=I, 3=rhombohedral obverse on hexagonal axes, 4=F, 5=A, 6=B, 7=C.
                # negative is non-centrosymmetric
                self.latt = self._assign_card(LATT(self, spline), line_num)
                if lastcard != 'ZERR' and (self.verbose or self.debug):
                    print('*** ZERR instruction is missing! ***')
                if self.latt.centric:
                    self.symmcards.set_centric(True)
            elif word == "SYMM":
                # SYMM symmetry operation
                #  Being more greedy, because many files do this wrong:
                # if not lastcard == 'ZERR':
                #    raise ParseOrderError
                # if not self.zerr:
                #    raise ParseOrderError
                s = SYMM(self, spline)
                if not self.latt and (self.debug or self.verbose):
                    print("*** LATT instruction is missing! ***")
                    if self.debug:
                        raise ParseSyntaxError(debug=self.debug, verbose=self.verbose)
                # Have to do this after parsing, because P-1 has no SYMM!
                # if self.latt.centric:
                #    self.symmcards.set_centric(True)
                self.symmcards.append(s.symmcard)
                if s not in self._reslist:
                    self._reslist[line_num] = s
                else:
                    self.delete_on_write.update([line_num])
                    self._reslist[line_num] = ' '
                lastcard = 'SYMM'
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
            elif word == 'UNIT':
                # UNIT n1 n2 ...
                # Number of atoms of each type in the unit-cell, in SFAC order.
                if lastcard != 'SFAC':
                    raise ParseOrderError(debug=shx.debug, verbose=shx.verbose)
                if self.sfac_table:
                    try:
                        self.unit = self._assign_card(UNIT(self, spline), line_num)
                    except ValueError:
                        if self.debug or self.verbose:
                            print('*** Non-numeric value in SFAC instruction! ***')
                        if self.debug:
                            raise
                else:
                    raise ParseOrderError(debug=shx.debug, verbose=shx.verbose)
                if len(self.unit.values) != len(self.sfac_table.elements_list) and (self.debug or self.verbose):
                    print('*** Number of UNIT and SFAC values differ! ***')
                    if self.debug:
                        raise ParseNumError(debug=self.shx.debug, verbose=self.shx.verbose)
                lastcard = 'UNIT'
            elif word in ['L.S.', 'CGLS']:
                # CGLS nls[0] nrf[0] nextra[0]
                # L.S. nls[0] nrf[0] nextra[0]
                self.cycles = self._assign_card(LSCycles(self, spline), line_num)
            elif word == "LIST":
                # LIST m[#] mult[1] (mult is for list 4 only)
                self.list = int(spline[1])
            elif word == "FVAR":
                # FVAR osf[1] free variables
                for fvvalue in spline[1:]:
                    fvarnum += 1
                    self._append_card(self.fvars, FVAR(fvarnum, float(fvvalue)), line_num)
                    if self.fvars not in self._reslist:
                        self._reslist[line_num] = self.fvars
                    else:
                        self.delete_on_write.update([line_num])
            elif word == 'ANIS':
                # ANIS n or ANIS names
                # Must be before Atom(), to know which atom is anis.
                self.anis = self._assign_card(ANIS(self, spline), line_num)
            elif word == 'WGHT':
                # WGHT a[0.1] b[0] c[0] d[0] e[0] f[.33333]
                if self.end:
                    self.wght_suggested = self._assign_card(WGHT(self, spline), line_num)
                    continue
                self.wght = self._assign_card(WGHT(self, spline), line_num)
            elif word == 'ACTA':
                # ACTA 2θfull[#] -> optional parameter NOHKL
                self.acta = self._assign_card(ACTA(self, spline), line_num)
            elif word == 'DAMP':
                # DAMP damp[0.7] limse[15]
                self.damp = self._assign_card(DAMP(self, spline), line_num)
            elif word == 'ABIN':
                # ABIN n1 n2
                self.abin = self._assign_card(ABIN(self, spline), line_num)
            elif word == 'ANSC':
                # ANSC six coefficients
                if len(spline) == 7:
                    self.ansc = [float(x) for x in spline[:1]]
            elif word == 'ANSR':
                # ANSR anres[0.001]
                if len(spline) == 2:
                    self.ansr = float(spline[1])
            elif word == 'BIND':
                # BIND atom1 atom2
                if len(spline) == 3:
                    self._append_card(self.bind, BIND(self, spline), line_num)
            elif word == 'BLOC':
                # BLOC n1 n2 atomnames
                self._append_card(self.bloc, BLOC(self, spline), line_num)
            elif word == 'BOND':
                # BOND atomnames
                self._append_card(self.bonds, BOND(self, spline), line_num)
            elif word == 'BUMP':
                # BUMP s [0.02]
                self._assign_card(BUMP(self, spline), line_num)
            elif word == 'CHIV':
                # CHIV V[0] s[0.1] atomnames
                self._append_card(self.restraints, CHIV(self, spline), line_num)
            elif word == 'CONF':
                # CONF atomnames max_d[1.9] max_a[170]
                self.conf = self._assign_card(CONF(self, spline), line_num)
            elif word == 'CONN':
                # CONN bmax[12] r[#] atomnames or CONN bmax[12]
                # bonded are d < (r1 + r2 + 0.5) Å
                self.conn = self._assign_card(CONN(self, spline), line_num)
            elif word == 'DEFS':
                # DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
                self.defs = self._assign_card(DEFS(self, spline), line_num)
            elif word == 'DISP':
                # DISP E f' f"[#] mu[#]
                if lastcard != 'SFAC':
                    raise ParseOrderError(debug=shx.debug, verbose=shx.verbose)
                self._append_card(self.disp, DISP(self, spline), line_num)
            elif word == 'EQIV':
                # EQIV $n symmetry operation
                # TODO: implement EQUIV class
                if len(spline) > 1 and spline[1].startswith('$'):
                    self.eqiv.append(spline[1:])
            elif word == 'EXTI':
                # EXTI x[0]
                self.exti = float(spline[1])
            elif word == 'EXYZ':
                # EXYZ atomnames
                self._append_card(self.restraints, EXYZ(self, spline), line_num)
            elif word == 'FRAG':
                # FRAG code[17] a[1] b[1] c[1] α[90] β[90] γ[90]
                if len(spline) == 8:
                    self.frag = self._assign_card(FRAG(self, spline), line_num)
            elif word == 'FEND':
                # FEND (must follow FRAG)
                if not self.frag:
                    raise ParseOrderError(debug=shx.debug, verbose=shx.verbose)
                self.frag = None  # Turns frag mode off.
            elif word == 'FLAT':
                # FLAT s[0.1] four or more atoms
                self._append_card(self.restraints, FLAT(self, spline), line_num)
            elif word == 'FREE':
                # FREE atom1 atom2
                self._append_card(self.free, FREE(self, spline), line_num)
            elif word == 'GRID':
                # GRID sl[#] sa[#] sd[#] dl[#] da[#] dd[#]
                self.grid = self._assign_card(GRID(self, spline), line_num)
            elif word == 'HKLF':
                # HKLF N[0] S[1] r11...r33[1 0 0 0 1 0 0 0 1] sm[1] m[0]
                self.hklf = self._assign_card(HKLF(self, spline), line_num)
            elif line.startswith('END'):
                # END (after HKLF or ends an include file)
                self.end = True
            elif word == 'HTAB':
                # HTAB dh[2.0]  or  HTAB donor-atom acceptor-atom
                self.htab = self._assign_card(HTAB(self, spline), line_num)
            elif word == 'ISOR':
                # ISOR s[0.1] st[0.2] atomnames
                self._append_card(self.restraints, ISOR(self, spline), line_num)
            elif word == 'LAUE':
                # LAUE E
                # I completely do not understand the LAUE instruction description in the manual!
                continue
            elif word == 'MERG':
                # MERG n[2]
                self.merg = self._assign_card(MERG(self, spline), line_num)
            elif word == 'MORE':
                # MORE m[1]
                self.more = self._assign_card(MORE(self, spline), line_num)
            elif word == 'FMAP':
                # FMAP code[2] axis[#] nl[53]
                self.fmap = self._assign_card(FMAP(self, spline), line_num)
            elif word == 'MOVE':
                # MOVE dx[0] dy[0] dz[0] sign[1]
                self.move = self._assign_card(MOVE(self, spline), line_num)
            elif word == 'MPLA':
                # MPLA na atomnames
                self.mpla = self._assign_card(MPLA(self, spline), line_num)
            elif word == 'NCSY':
                # NCSY DN sd[0.1] su[0.05] atoms
                self._append_card(self.restraints, NCSY(self, spline), line_num)
            elif word == 'NEUT':
                # NEUT
                # TODO: Implement NEUT class
                if lastcard != 'SYMM':
                    raise ParseOrderError(debug=shx.debug, verbose=shx.verbose)
            elif word == 'OMIT':
                # OMIT atomnames  or  OMIT s[-2] 2θ(lim)[180]  or  OMIT h k l
                # TODO: Implement OMIT class
                self.omit.append(spline[1:])
            elif word == 'PLAN':
                # PLAN npeaks[20] d1[#] d2[#]
                self.plan = self._assign_card(PLAN(self, spline), line_num)
            elif word == 'PRIG':
                # PRIG p[#]
                self.prig = self._assign_card(PRIG(self, spline), line_num)
            elif word == 'RTAB':
                # RTAB codename atomnames  -->  codename: e.g. 'omeg' gets tabualted in the lst
                self._append_card(self.rtab, RTAB(self, spline), line_num)
            elif word == 'SAME':
                # SAME s1[0.02] s2[0.04] atomnames
                self._append_card(self.restraints, SAME(self, spline), line_num)
            elif word == 'SHEL':
                # SHEL lowres[infinite] highres[0]
                self.shel = self._assign_card(SHEL(self, spline), line_num)
            elif word == 'SIZE':
                # SIZE dx dy dz
                self.size = self._assign_card(SIZE(self, spline), line_num)
            elif word == 'SPEC':
                # SPEC del[0.2]
                if len(spline) > 1:
                    self.spec = self._assign_card(SPEC(self, spline), line_num)
            elif word == 'STIR':
                # STIR sres step[0.01]   -> stepwise improvement in the resolution sres
                self.stir = self._assign_card(STIR(self, spline), line_num)
            elif word == 'SUMP':
                # SUMP c sigma c1 m1 c2 m2 ...
                self._append_card(self.sump, SUMP(self, spline), line_num)
            elif word == 'SWAT':
                # SWAT g[0] U[2]
                self.swat = self._assign_card(SWAT(self, spline), line_num)
            elif word == 'TEMP':
                # TEMP T[20]  -> in Celsius
                self.temp = float(spline[1].split('(')[0])
                self.temp_in_kelvin = self.temp + 273.15
            elif word == 'TWIN':
                # TWIN 3x3 matrix [-1 0 0 0 -1 0 0 0 -1] N[2]
                self.twin = self._assign_card(TWIN(self, spline), line_num)
            elif word == 'TWST':
                # TWST N[0] (N[1] after SHELXL-2018/3)
                if len(spline) > 1:
                    self.twst = self._assign_card(TWST(self, spline), line_num)
            elif word == 'WIGL':
                # WIGL del[0.2] dU[0.2]
                self.wigl = self._assign_card(WIGL(self, spline), line_num)
            elif word == 'WPDB':
                # WPDB n[1]
                self.wpdb = self._assign_card(WPDB(self, spline), line_num)
            elif word == 'XNPD':
                # XNPD Umin[-0.001]
                self.xnpd = self._assign_card(XNPD(self, spline), line_num)
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
            elif line.startswith('+'):
                pass
            elif word == 'TITL':
                self.titl = line[5:76]
                lastcard = 'TITL'
            else:
                if not line.strip():
                    continue
                if self.debug or self.verbose:
                    print("Error in line: {} -> {}".format(line_num + 1, line))
                    if self.debug:
                        raise ParseUnknownParam(debug=self.debug, verbose=self.verbose)

    def add_atom(self, name: str = None, coordinates: list = None, element='C', uvals: list = None, part: int = 0,
                 sof: float = 11.0):
        """
        Adds an atom to the ShelxFile.atoms list. If no element is given, carbon atoms are assumed.
        """
        if uvals is None:
            uvals = [0.04]
        part = PART(self, f'PART {part}'.split())
        afix = AFIX(self, 'AFIX 0'.split())
        resi = RESI(self, 'RESI 0'.split())
        a = Atom(self)
        sfac_num = self.elem2sfac(element)
        a.set_atom_parameters(name=name, sfac_num=sfac_num, coords=coordinates,
                              part=part, afix=afix, resi=resi, site_occupation=sof, uvals=uvals)
        self._append_card(self.atoms, a, 0)

    def frac_to_cart(self, coordinates: list) -> Array:
        """
        fractional to cartesian coordinates by applying the orthogonal matrix.
        """
        return self.orthogonal_matrix * Array(coordinates)

    def __repr__(self):
        """
        Represents the shelxl object.
        """
        resl = []
        for num, line in enumerate(self._reslist):
            if num in self.delete_on_write:
                continue
            try:
                if line == '' and self._reslist[num + 1] == '':
                    continue
            except IndexError:
                pass
            # Prevent wrapping long lines with \n breaks by splitting first:
            line = "".join([wrap_line(x) for x in str(line).split("\n")])
            resl.append(line)
        return "\n".join(resl)

    def grow(self, with_qpeaks: bool = False):
        """
        Returns a list of atoms that represent the complete molecules of the structure.
        """
        sdm = SDM(self)
        needsymm = sdm.calc_sdm()
        packed_atoms = sdm.packer(sdm, needsymm, with_qpeaks=with_qpeaks)
        return packed_atoms

    def refine(self, cycles: Union[int, None] = None, backup_before: bool = True) -> bool:
        if self.resfile:
            filen = self.resfile.stem
            if cycles is not None:
                self.cycles.number = cycles
            ref = ShelxlRefine(self, self.resfile)
            ref.remove_acta_card(self.acta)
            self.write_shelx_file(filen + '.ins')
            ref.run_shelxl(backup_before=backup_before)
            self.reload()
            ref.restore_acta_card()
            # self.write_shelx_file(filen + '.res')
            return True
        return False

    def refine_weight_convergence(self, stop_after: int = 10):
        """
        Tries to refine weigting sheme from SHELXL until it converged (self.weight_difference() is zero) or
        stopt_after cycles are reached.
        """
        for _ in range(stop_after):
            difference = self.wght.difference()
            print("Weighting difference = {} {}".format(*difference))
            if self._weight_converged(difference):
                return True
            else:
                self.update_weight()
                self.refine(9)
        print("Maximum number of refinement cycles reached, but no WGHT convergence.")
        return False

    def _weight_converged(self, diff: List[float]) -> bool:
        return diff == [0.0, 0.0]

    def _append_card(self, obj, card, line_num: int) -> Command:
        """
        Appends SHELX card to an object list, e.g. self.restraints and
        assigns the line_num in reslist with the card instance.
        """
        obj.append(card)
        self._reslist[line_num] = card
        return card

    def _assign_card(self, card, line_num: int):
        self._reslist[line_num] = card
        return card

    @staticmethod
    def is_atom(atomline: str) -> bool:
        """
        Returns True is line contains an atom.
        """
        # no empty line, not in cards and not space at start:
        if atomline[:4].upper() not in SHX_CARDS:  # exclude all non-atom cards
            spline: List[str] = atomline.split()
            # Too few parameter for an atom:
            if len(spline) < 5:
                return False
            # means sfac number is missing:
            if '.' in spline[1]:
                return False
            if Shelxfile._coordinates_are_unrealistic(spline):
                return False
            # Exclude lone pairs:
            if len(spline) > 5 and spline[5] == '!':
                return False
            return True
        else:
            return False

    @staticmethod
    def _coordinates_are_unrealistic(spline: List[str]) -> bool:
        return any(float(y) > 4.0 for y in spline[2:5])

    def to_cif(self, filename: str = None, template: Optional[str] = None) -> None:
        """
        Writes a CIF file from the ShelxFile object.
        """
        if not filename:
            filename = self.resfile.stem + '.cif'
        CifFile(self, template).write_cif(Path(filename))

    def elem2sfac(self, atom_type: str) -> int:
        """
        returns an sfac-number for the element given in "atom_type"
        """
        for num, element in enumerate(self.sfac_table, 1):
            if atom_type.capitalize() == element.capitalize():
                return num  # return sfac number
        # Element was not found in sfac table
        return 0

    def sfac2elem(self, sfacnum: int) -> str:
        """
        returns an element and needs an sfac-number
        :param sfacnum: string like '2'
        """
        try:
            elem = self.sfac_table[int(sfacnum)]
        except IndexError:
            return ''
        return elem

    def add_line(self, linenum: int, line: str) -> None:
        """
        Adds a new SHELX card to the reslist after linenum.
        e.g. shx.add_line(shx.unit.position, 'ANIS')
        """
        self._reslist.insert(linenum + 1, line)

    def replace_line(self, obj, new_line: str) -> None:
        """
        Replaces a single line in the res file with new_line.
        """
        self._reslist[self.index_of(obj)] = new_line

    def index_of(self, obj: Union[Atom, Restraint, Command]) -> int:
        return self._reslist.index(obj)

    @property
    def sum_formula(self) -> str:
        """
        The sum formula of the structure in regard to the UNIT instruction.
        """
        formstring = ''
        formula_weight = 0.0
        try:
            val = self.unit.values
            eli = self.sfac_table.elements_list
        except AttributeError:
            return ''
        if len(val) == len(eli):
            for el, num in zip(self.sfac_table.elements_list, self.unit.values):
                try:
                    elcount = num / self.Z
                    formula_weight += elcount * float(weight_from_symbol(el.capitalize()))
                    formstring += f"{el}{elcount :,g} "
                except ZeroDivisionError:
                    return ''
        self.formula_weight = round(formula_weight, 3)
        return formstring.strip()

    def update_weight(self) -> None:
        try:
            self.wght.a = self.wght_suggested.a
            self.wght.b = self.wght_suggested.b
            self.wght.c = self.wght_suggested.c
            self.wght.d = self.wght_suggested.d
            self.wght.e = self.wght_suggested.e
            self.wght.f = self.wght_suggested.f
        except AttributeError:
            return

    def insert_anis(self, atoms: str = '', residue: str = ''):
        """
        Inserts ANIS into a results file for refinement of anisotropic displacement parameters.

        TODO: implement ANIS n
        @param atoms: Specify secific atoms or wildcards of atoms.
        @param residue: Specify a residue like ANIS_ABC with residue='ABC' or ANIS_* (residue='*')
        """
        if atoms:
            self.add_line(self.unit.position, f'ANIS{"_" if residue else ""}{residue} {atoms}')
        else:
            self.add_line(self.unit.position, 'ANIS')

    @property
    def sum_formula_exact(self) -> str:
        """
        The sum formula of the structure with all atom occupancies summed together as string.
        """
        formstring = ''
        sumdict = self.sum_formula_exact_as_dict()
        for el in sumdict:
            formstring += f"{el.capitalize()}{round(sumdict[el], 2):,g} "
        return formstring.strip()

    def sum_formula_exact_as_dict(self) -> dict:
        """
        The sum formula of the structure with all atom occupancies summed together as dictionary.
        """
        sumdict = {}
        for el in self.sfac_table.elements_list:
            for atom in self.atoms:
                if atom.element.upper() == el.upper() and not atom.qpeak:
                    if el in sumdict:
                        sumdict[el] += atom.occupancy
                    else:
                        sumdict[el] = atom.occupancy
            if el not in sumdict:
                sumdict[el] = 0.0
        return sumdict

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

    def _get_residuals(self, spline: List[str], line: str) -> None:
        if Shelxfile._r1_regex.match(line):
            self._get_r1(spline)
        if Shelxfile._wr2_regex.match(line):
            self._get_wr2(spline)
        if Shelxfile._parameters_regex.match(line):
            self._get_params_and_restraints(spline)
        if Shelxfile._diff_peak_regex.match(line):
            self._get_peak_hole(spline)
        if Shelxfile._goof_regex.match(line):
            self._get_goof(spline)
        if Shelxfile._spgrp_regex.match(line):
            self._get_space_group(spline)

    def _get_space_group(self, spline: List[str]) -> None:
        try:
            self.space_group = spline[3]
        except(IndexError, ValueError):
            pass

    def _get_goof(self, spline: List[str]) -> None:
        with suppress(IndexError, ValueError):
            self.goof = float(spline[8].split(',')[0])
            self.rgoof = float(spline[12].split(',')[0])

    def _get_peak_hole(self, spline: List[str]) -> None:
        # REM Highest difference peak  0.407,  deepest hole -0.691,  1-sigma level  0.073
        with suppress(IndexError, ValueError):
            self.highest_peak = float(spline[4].split(",")[0])
            self.deepest_hole = float(spline[7].split(",")[0])

    def _get_params_and_restraints(self, spline: List[str]) -> None:
        with suppress(IndexError):
            self.parameters = int(spline[1])
            if self.data and self.parameters:
                self.dat_to_param = float(self.data) / float(self.parameters)
        with suppress(IndexError, ValueError):
            self.num_restraints = int(spline[-2])

    def _get_wr2(self, spline: List[str]) -> None:
        with suppress(IndexError, ValueError):
            self.wr2 = float(spline[3].split(",")[0])

    def _get_r1(self, spline: List[str]) -> None:
        with suppress(IndexError, ValueError):
            self.R1 = float(spline[3])
        with suppress(IndexError, ValueError):
            self.data = int(spline[-2])


if __name__ == "__main__":
    print(Path('.').resolve())
    # file = r'../shelxfile/tests/resources/p21c.res'
    file = r'tests/resources/p-31c.res'
    shx = Shelxfile()
    shx.read_file(file)
    print(shx.atoms)
    print(shx.sum_formula_exact)
    print(shx.sum_formula)
    print(shx.sum_formula_exact_as_dict())
    print(shx.restraints)
    print(shx.atoms.nameslist)
    sys.exit()

    # noinspection PyUnreachableCode
    """
    #To get all available SHELX commands:
    def get_shelx_commands():
        url = "http://shelx.uni-goettingen.de/shelxl_html.php"
        response = urlopen('{}/version.txt'.format(url))
        html = response.read().decode('UTF-8')
        #res = BeautifulSoup(html, "html5lib")
        tags = res.findAll("p", {"class": 'instr'})
        for l in tags:
            if l:
                print(str(l).split(">")[1].split("<")[0])
    """
