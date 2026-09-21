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
from __future__ import annotations

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
from typing import cast

from shelxfile.atoms.atom import Atom, BedeLoneResultAtom
from shelxfile.atoms.atoms import Atoms
from shelxfile.cif.cif_write import CifFile
from shelxfile.edit.card_meta import CardLifetime
from shelxfile.misc.dsrmath import Array
from shelxfile.misc.elements import weight_from_symbol
# noinspection PyUnresolvedReferences
from shelxfile.misc.misc import ParseOrderError, ParseNumError, ParseUnknownParam, \
    multiline_test, dsr_regex, wrap_line, ParseSyntaxError, cart_to_frac, \
    split_fvar_and_parameter
from shelxfile.refine.refine import ShelxlRefine
from shelxfile.shelx.cards import ACTA, FVAR, FVARs, REM, BOND, Restraints, DEFS, NCSY, ISOR, FLAT, \
    BUMP, DFIX, DANG, SADI, SAME, RIGU, SIMU, DELU, CHIV, EADP, EXYZ, DAMP, HFIX, HKLF, SUMP, SYMM, LSCycles, \
    SFACTable, UNIT, BASF, TWIN, WGHT, BLOC, SymmCards, CONN, CONF, BIND, DISP, GRID, HTAB, MERG, FRAG, FREE, FMAP, \
    MOVE, PLAN, PRIG, RTAB, SHEL, SIZE, SPEC, STIR, TWST, WIGL, WPDB, XNPD, ZERR, CELL, LATT, MORE, MPLA, AFIX, PART, \
    RESI, ABIN, ANIS, Residues, SWAT, Command, Restraint, BEDE, LONE, EQIV
from shelxfile.shelx.sdm import SDM
from shelxfile.version import VERSION

__version__ = VERSION

ResListEntry = str | Command | Restraint | SFACTable | FVARs | Atom | SYMM | BedeLoneResultAtom

#: Fractional coordinates are small numbers, so a much larger value on an
#: atom line usually means the line is not an atom at all.
MAX_PLAIN_COORDINATE = 4.0

#: How large a *decoded* ``10*m + p`` coordinate may be. A fractional
#: coordinate lies in ``[0, 1]``, well inside the ``abs(p) < 5`` the
#: manual allows for the encoding in general.
MAX_DECODED_COORDINATE = 1.0


def _decodes_to_a_coordinate(value: float) -> bool:
    """Whether *value* is a coordinate written in SHELXL's ``10*m + p`` form.

    *"To fix any atom parameter, add 10"*, and more generally *"if any
    atom parameter is given as (10*m + p), where abs(p) is less than 5
    and m is an integer, it is interpreted as p*fvm"*.  A coordinate
    fixed at 0.6666 is therefore written ``10.666600``, and one tied to
    free variable 2 is written ``21.000000`` -- both far outside the
    plain fractional range, both perfectly ordinary atoms.

    The one code deliberately refused is ``m = 1`` with ``abs(p) >= 1``,
    of which ``11.00000`` is the only common case: as a coordinate it
    would mean a site fixed exactly on the cell edge, but it is the
    standard way to write a fixed full occupancy, so on an atom line it
    almost always means a coordinate is missing and the sof has slid
    into its column.  ``m >= 2`` carries no such ambiguity -- there is no
    free variable 1 -- so ``21.000000`` is read as ``1.0 * fv2``.

    :class:`~shelxfile.atoms.atom.Atom` already decodes these; without
    this check :meth:`Shelxfile.is_atom` would reject the line first and
    the atom would be silently kept as raw text -- and its continuation
    line would then be swallowed by the ``=`` of the line above.
    """
    fvar, decoded = split_fvar_and_parameter(value)
    if abs(decoded) > MAX_DECODED_COORDINATE:
        return False
    if abs(decoded) >= MAX_DECODED_COORDINATE and abs(fvar) < 2:
        return False
    return True

"""
TODO:
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

# Matches a trailing '_$n' EQIV symmetry-equivalent atom suffix, e.g. 'H9a_$1'.
_EQIV_ATOM_SUFFIX_RE = re.compile(r'^(.*)_\$(\d+)$')


class Shelxfile:
    """
    Class for data from a SHELXL res file. Includes Atoms, cards and unit cell.
    """
    _r1_regex = re.compile(r'^REM\s+R1\s+=', re.IGNORECASE)
    _wr2_regex = re.compile(r'^REM\s+wR2\s+=', re.IGNORECASE)
    _parameters_regex = re.compile(r'^REM\s+\d+\s+parameters\s+refined', re.IGNORECASE)
    _diff_peak_regex = re.compile(r'^REM\sHighest\sdifference', re.IGNORECASE)
    _goof_regex = re.compile(r'^REM\swR2\s=\s.*,\sGooF', re.IGNORECASE)
    _spgrp_regex = re.compile(r'^REM\s+\S+\s+in\s+\S+', re.IGNORECASE)

    #: How many individually unparsable instruction lines to tolerate
    #: before giving up on a file. A handful signals an oddity worth
    #: reporting; dozens signal that this is not a SHELXL file at all.
    MAX_UNPARSABLE_LINES = 20

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
        self.restraint_errors: list[str] = []
        self.temp_in_kelvin: float = 0.0
        self.shelx_max_line_length: int = 79  # maximum character lenth per line in SHELXL
        self.cell: CELL | None = None
        self.ansc: list[float] = []
        self.abin: ABIN | None = None
        self.acta: ACTA | None = None
        self.fmap: FMAP | None = None
        self.xnpd: XNPD | None = None
        self.wpdb: WPDB | None = None
        self.wigl: WIGL | None = None
        self.temp: int | float = 20
        self.swat: SWAT | None = None
        self.stir: STIR | None = None
        self.spec: SPEC | None = None
        self.twst: TWST | None = None
        self.plan: PLAN | None = None
        self.prig: PRIG | None = None
        self.merg: MERG | None = None
        self.more: MORE | None = None
        self.move: MOVE | None = None
        self.defs: DEFS | None = None
        self.zerr: ZERR | None = None
        self.wght: WGHT | None = None
        self.frag: FRAG | None = None
        self.twin: TWIN | None = None
        self.basf: BASF | None = None
        self.latt: LATT | None = None
        self.anis: ANIS | None = None
        self.damp: DAMP | None = None
        self.unit: UNIT | None = None
        self.size: SIZE | None = None
        self.htab: HTAB | None = None
        self.shel: SHEL | None = None
        self.mpla: MPLA | None = None
        self.hklf: HKLF | None = None
        self.grid: GRID | None = None
        self.conn: CONN | None = None
        self.conf: CONF | None = None
        self.afix: AFIX | None = None
        self.rtab: list[RTAB] = []
        self.omit: list[list[str]] = []
        self.free: list[FREE] = []
        self.eqiv: list[EQIV] = []
        self.bonds: list[BOND] = []
        self.disp: list[DISP] = []
        self.bind: list[BIND] = []
        self.bloc: list[BLOC] = []
        self.part: PART = PART(self, ['PART', '0'])
        self.resi: RESI = RESI(self, ['RESI', '0'])
        self.residues = Residues(self)
        self.dsrlines: list[str] = []
        self.dsrline_nums: list[int] = []
        self.symmcards: SymmCards = SymmCards(self)
        self.hfixes: list[HFIX] = []
        self.sump: list[SUMP] = []
        self.wght_suggested: WGHT | None = None
        self.Z: int = 1
        self.titl: str = ""
        self.exti: float = 0.0
        self.ansr: float = 0.001
        self.rem: list[REM] = []
        self.atoms: Atoms = Atoms(self)
        self.bede_cards: list[BEDE] = []
        self.lone_cards: list[LONE] = []
        self.bede_lone_results: list[BedeLoneResultAtom] = []
        self.fvars: FVARs = FVARs(self)
        self.restraints: Restraints = Restraints()
        self.sfac_table: SFACTable = SFACTable(self)
        self.cycles: LSCycles | None = None
        self.R1: float | None = None
        self.wr2: float | None = None
        self.goof: float | None = None
        self.rgoof: float | None = None
        self.space_group: str | None = None
        self.data: int | None = None
        self.parameters: int | None = None
        self.dat_to_param: float | None = None
        self.num_restraints: int | None = None
        self.highest_peak: float | None = None
        self.deepest_hole: float | None = None
        self.formula_weight: float | None = None
        self.end: bool = False
        self.maxsof: float = 1.0
        self.delete_on_write: set[int] = set()
        self.wavelength: float = 0.0
        self.global_sadi: int | None = None
        self.list: int = 0
        self.theta_full: float = 0.0
        self.error_line_num: int = -1  # Only used to tell the line number during an exception.
        #: Codec the file was decoded with, reused when writing so that an
        #: unedited file comes back byte for byte. See :meth:`_read_text`.
        self.encoding: str = 'utf-8'
        self.resfile: Path | None = None
        self.orthogonal_matrix: Array | None = None
        self._reslist: list[ResListEntry] = []
        #: Bumped on every structural change. Lets the edit layer notice a
        #: stale cache without the model having to know about it.
        self.model_version: int = 0
        #: Instruction lines that could not be parsed, as messages.
        self.parse_errors: list[str] = []
        #: ``_reslist`` indices quarantined after a parse failure. They
        #: stay in the file as text so they still round-trip.
        self._unparsable_lines: set[int] = set()
        #: Pristine copy of the input lines, so parsing can restart.
        self._source_lines: list[ResListEntry] = []
        #: Source line of the card currently being parsed, with '='
        #: continuations glued on. Parse-time scratch used by
        #: :meth:`_tag_lifetime`.
        self._current_raw_line: str = ''

    def touch(self) -> None:
        """Record that the structure changed."""
        self.model_version += 1

    def dumps(self) -> str:
        """
        Returns the current content of the SHELX file as a single string, with
        every line wrapped after 79 characters as SHELXL requires. This is the
        single source of truth used by both :meth:`write_shelx_file` and
        :meth:`__repr__` so the two never drift apart.
        """
        resl = []
        for num, line in enumerate(self._reslist):
            if num in self.delete_on_write:
                continue
            if line == '':  # and self._reslist[num + 1] == '':
                continue
            # Prevent wrapping long lines with \n breaks by splitting first:
            line = "\n".join([wrap_line(x) for x in str(line).split("\n")])
            resl.append(line)
        return "\n".join(resl)

    def write_shelx_file(self, filename: str | Path | None = None) -> None:
        if not self._reslist:
            print('*** No file was loaded for writing. ***')
            return None
        if not filename:
            filename = self.resfile
        if isinstance(filename, str):
            filename = Path(filename)
        filename = cast(Path, filename)
        with open(filename, 'w', encoding=self.encoding) as f:
            f.write(self.dumps() + '\n')
        if self.verbose or self.debug:
            print(f'*** File successfully written to {filename.resolve()} ***')

    def read_file(self, resfile: Path | str) -> None:
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
            text = self._read_text(resfile)
            self._reslist = cast(list[ResListEntry], text.splitlines(keepends=False))
            self._test_if_file_is_valid(resfile)
        except UnicodeDecodeError:
            if self.debug or self.verbose:
                print(f'*** Unable to read file {resfile} ***')
            return
        self._find_included_files()
        self.parse_cards()

    def _read_text(self, resfile: Path) -> str:
        """Decode *resfile*, remembering which encoding worked.

        SHELX files predate Unicode, and non-ASCII bytes turn up in ``REM``
        comments and titles.  Relying on the platform default makes the
        same file decode differently on Windows and Linux, so the encoding
        is pinned here and reused when writing: a file read and written
        unchanged must come back byte for byte.

        ``latin-1`` is the fallback because it maps every byte, so a file
        in some other 8-bit codepage still round-trips exactly even though
        the characters may display oddly.
        """
        for encoding in ('utf-8', 'latin-1'):
            try:
                text = resfile.read_text(encoding=encoding)
            except UnicodeDecodeError:
                continue
            self.encoding = encoding
            return text
        raise UnicodeDecodeError('utf-8', b'', 0, 1, 'undecodable file')

    def read_string(self, resfile_string: str) -> None:
        """
        Read input as string.
        This will not read files included with "+filename" syntax!
        """
        self.__init__(debug=self.debug, verbose=self.verbose)
        self._reslist = cast(list[ResListEntry], resfile_string.splitlines(keepends=False))
        self.parse_cards()

    def parse_cards(self) -> None:
        """Parse the instruction lines, tolerating individual bad ones.

        A single malformed instruction used to cost the whole file: the
        exception propagated out of :meth:`_parse_cards`, was swallowed
        here, and the caller received a model silently missing every card
        and atom below the failure.  In the reference corpus that hit
        7.6 % of files -- one ``TEMP -100.0C`` discarded six ``MPLA``
        cards and the atom list beneath it.

        SHELXL keeps reading past an instruction it cannot use, so the
        offending line is recorded in :attr:`parse_errors`, quarantined,
        and parsing restarts.  The line stays in ``_reslist`` as text, so
        it still round-trips to the output.

        Debug mode raises instead, which is what it is for.
        """
        self._source_lines = list(self._reslist)
        for _ in range(self.MAX_UNPARSABLE_LINES + 1):
            try:
                self._parse_cards()
                break
            except Exception as error:
                if self.debug or self.verbose:
                    self.show_line_where_error_occured(error)
                if self.debug:
                    raise
                failed_at = self.error_line_num
                if failed_at < 0 or failed_at in self._unparsable_lines:
                    return
                self._record_line_error(failed_at, error)
                self._unparsable_lines.add(failed_at)
                self._restart_parsing()
        else:
            return
        self.restraint_errors = self._assign_atoms_to_restraints()

    def _restart_parsing(self) -> None:
        """Reset parse state, keeping the source text and what we learnt."""
        source = self._source_lines
        quarantined = self._unparsable_lines
        errors = self.parse_errors
        resfile = self.resfile
        self.__init__(debug=self.debug, verbose=self.verbose)
        self.resfile = resfile
        self._source_lines = source
        self._unparsable_lines = quarantined
        self.parse_errors = errors
        self._reslist = list(source)

    def _record_line_error(self, line_num: int, error: Exception) -> None:
        """Note that one instruction could not be parsed, and carry on."""
        try:
            text = str(self._reslist[line_num])
        except IndexError:
            text = ''
        self.parse_errors.append(
            f'*** Could not parse line {line_num + 1}: {text.strip()!r} '
            f'({type(error).__name__}: {error}) ***'
        )
        if self.verbose:
            print(self.parse_errors[-1])

    def _assign_atoms_to_restraints(self) -> list[str]:
        warnings: list[str] = []
        populated_resi_nums_by_class: dict[str, set[int]] = {}
        for atom in self.atoms:
            cls = atom.resiclass
            num = atom.resinum
            if cls not in populated_resi_nums_by_class:
                populated_resi_nums_by_class[cls] = set()
            populated_resi_nums_by_class[cls].add(num)
        # Warn about residues that are defined (via RESI) but contain no atoms.
        empty_residues: list[str] = []
        for resi in self.residues.all_residues:
            if resi.residue_number not in populated_resi_nums_by_class.get(resi.residue_class, set()):
                empty_residues.append(f'{resi.residue_class} {resi.residue_number}')
        if empty_residues:
            warnings.append(f'*** Empty residue(s) detected (no atoms): {", ".join(empty_residues)} ***')
        for restraint in self.restraints:
            if not restraint.lifetime:
                # After END: inert output from the previous refinement,
                # which SHELXL never reads back. Its atom names may carry
                # ACTA TABS part suffixes (C1^a) that deliberately do not
                # resolve against the atom list. See CardLifetime.
                continue
            bad_atoms = []
            missing_eqiv = []
            for restraint_atom in restraint.atoms:
                if restraint_atom in ('>', '<', '='):
                    continue
                # A trailing '_$n' is a symmetry reference, not a residue, so
                # the name underneath may still need the instruction's residue
                # scope applied to it: 'RTAB_23 ... O1_$3' means (O1_23)_$3.
                eqiv_match = _EQIV_ATOM_SUFFIX_RE.match(restraint_atom)
                base_name = eqiv_match.group(1) if eqiv_match else restraint_atom
                has_own_residue = '_' in base_name
                if (restraint.residue_class or sum(restraint.residue_number) > 0) and not has_own_residue:
                    populated_nums = populated_resi_nums_by_class.get(restraint.residue_class, set())
                    # A class-scoped instruction is applied once per residue,
                    # and "if some or all of the named atoms cannot be found
                    # for a particular residue, the instruction is simply
                    # ignored for that residue". So a name missing from one
                    # residue is normal; only a name missing from *every*
                    # residue is worth reporting.
                    per_residue: list[str] = []
                    found_somewhere = False
                    for num in restraint.residue_number:
                        # Skip residue numbers of this class that have no atoms at all.
                        # SHELXL silently ignores empty residues for restraints too.
                        if restraint.residue_class and num not in populated_nums:
                            continue
                        probe: list[str] = []
                        if eqiv_match:
                            self.does_atom_exist(restraint_atom, probe, restraint_atom,
                                                 missing_eqiv, residue_scope=num)
                        else:
                            self.does_atom_exist(f'{restraint_atom}_{num}', probe, f'{restraint_atom}_{num}',
                                                 missing_eqiv)
                        if probe:
                            per_residue.extend(probe)
                        else:
                            found_somewhere = True
                    if not found_somewhere:
                        bad_atoms.extend(per_residue)
                elif '_' in restraint_atom:
                    self.does_atom_exist(f'{restraint_atom}', bad_atoms, restraint_atom, missing_eqiv)
                else:
                    self.does_atom_exist(f'{restraint_atom}_{0}', bad_atoms, restraint_atom, missing_eqiv)
            if bad_atoms:
                sorted_atoms = list(set(bad_atoms))
                sorted_atoms.sort()
                warnings.append(f'*** Unknown atom{"s" if len(bad_atoms) > 0 else ""} in restraint: {restraint}, '
                                f'line {restraint.index + 1} ***')
                warnings.append(f'*** Atom list has no --> {", ".join(sorted_atoms)} ***')
                bad_atoms.clear()
            if missing_eqiv:
                sorted_missing = sorted(set(missing_eqiv))
                warnings.append(f'*** Undefined EQIV in restraint: {restraint}, line {restraint.index + 1} ***')
                warnings.append(f'*** No EQIV instruction defines --> {", ".join(sorted_missing)} ***')
                missing_eqiv.clear()
            if restraint.residue_class and sum(restraint.residue_number) == 0:
                warnings.append(f"*** Restraint '{restraint}', line {restraint.index + 1}, "
                                f"has a residue class, but no residues are defined. ***")
        if self.debug or self.verbose:
            print('\n'.join(warnings))
        return warnings

    def does_atom_exist(
        self,
        atom_name: str,
        bad_atoms: list[str],
        restraint_atom: str,
        missing_eqiv: list[str],
        residue_scope: int = 0,
    ) -> None:
        # A trailing '_$n' references a symmetry equivalent atom defined by an EQIV
        # instruction (see EQIV documentation), not a residue number. It has to be
        # stripped off before checking whether the underlying atom really exists.
        eqiv_match = _EQIV_ATOM_SUFFIX_RE.match(atom_name)
        eqiv_id = None
        if eqiv_match:
            base_name, eqiv_num = eqiv_match.groups()
            eqiv_id = f'${eqiv_num}'
            # The manual applies the residue *before* the symmetry operation,
            # so an unqualified name inherits the instruction's residue rather
            # than defaulting to residue 0.
            atom_name = base_name if '_' in base_name else f'{base_name}_{residue_scope}'
        residue_number_is_wildcard = '_' in atom_name and atom_name.split('_')[-1] == '*'
        if atom_name.startswith('$'):
            return None
        bad_atoms_before = len(bad_atoms)
        if residue_number_is_wildcard:
            for num in self.residues.residue_numbers.keys():
                residue_atom = f"{atom_name.split('_')[0]}_{num}"
                if not self.atoms.get_atom_by_name(residue_atom):
                    bad_atoms.append(residue_atom)
        else:
            if not self.atoms.get_atom_by_name(atom_name):
                bad_atoms.append(restraint_atom)
        atom_was_found = len(bad_atoms) == bad_atoms_before
        # Only complain about a missing EQIV definition if the underlying atom itself
        # exists; otherwise the 'unknown atom' warning above already covers it.
        if eqiv_id and atom_was_found and not any(entry.id == eqiv_id for entry in self.eqiv):
            missing_eqiv.append(restraint_atom)

    def _test_if_file_is_valid(self, resfile: Path) -> None:
        if len(self._reslist) < 20 and (self.debug or self.verbose):
            print('*** Not a SHELXL file: {} ***'.format(resfile))
            if self.debug:
                sys.exit()

    def show_line_where_error_occured(self, e: Exception) -> None:
        try:
            print(f'Error near:\n {self._reslist[self.error_line_num]}')
        except IndexError:
            pass
        print(e)
        print(f"*** Syntax error found in file {self.resfile}, line {self.error_line_num + 1} ***")

    def _find_included_files(self) -> None:
        # Tracks the file names of included files in order to find recursive inclusion:
        includefiles: list[str] = []
        for line_num, line in enumerate(self._reslist):
            if not isinstance(line, str):
                continue
            if line.startswith('+'):
                try:
                    file_included_in_includefile = self._read_included_file(includefiles, line)
                    if file_included_in_includefile:
                        for line_num_includefile, include_line in enumerate(file_included_in_includefile):
                            reslist_position = line_num + 1 + line_num_includefile
                            # '+filename' include files are not copied to res file,
                            #  so I have to delete these lines on write.
                            # '++filename' copies them to the .res file where appropriate
                            # I leave this out, because I am not SHELXL:
                            # if include_line.startswith('+') and include_line[:2] != '++':
                            #    self.delete_on_write.update([lnum])
                            self._reslist.insert(reslist_position, include_line)
                        continue
                except IndexError:
                    if self.debug or self.verbose:
                        print(f'*** CANNOT READ INCLUDE FILE {line} ***')
                    # Not sure if this is a good idea: del reslist[n]

    def _read_included_file(self, includefiles: list[str], line: str) -> list[str]:
        include_filename: Path = cast(Path, self.resfile).resolve().parent.joinpath(line[1:])
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
        self.read_file(cast(Path, self.resfile).resolve())

    def _parse_cards(self) -> None:
        last_nonhydrogen_atom: Atom | None = None
        # Reference for riding U values: the previous atom whose own U is not a
        # riding code. Unlike the pivot atom this is not restricted to hydrogen.
        last_unconstrained_atom: Atom | None = None
        lastcard = ''
        fvarnum = 1
        for line_num, line in enumerate(self._reslist):
            if not isinstance(line, str):
                continue
            if line_num in self._unparsable_lines:
                # Quarantined by an earlier pass: kept as text so it still
                # round-trips, but not interpreted.
                continue
            self.error_line_num = line_num  # For exception during parsing.
            list_of_lines = [line_num]  # list of lines where a card appears, e.g. for atoms with two lines
            if line.startswith(' ') or line == '':
                continue
            wrapindex = 0
            # This while loop makes wrapped lines look like they are not wrapped. The following lines are then
            # beginning with a space character and thus are ignored. The 'lines' list holds the line nnumbers where
            # 'line' is located ([line_num]) plus the wrapped lines.
            # The untouched physical lines, kept so post-END output can be
            # echoed exactly as read (see _tag_lifetime).
            raw_physical_lines: list[str] = [line]
            if multiline_test(line):
                multiline = True
            else:
                multiline = False
            while multiline:
                # Glue together the two lines wrapped with "=":
                wrapindex += 1
                wrapped_line = cast(str, self._reslist[line_num + wrapindex])
                raw_physical_lines.append(wrapped_line)
                line = line.rpartition('=')[0] + wrapped_line
                # self.delete_on_write.update([line_num + wrapindex])
                list_of_lines.append(line_num + wrapindex)  # list containing the lines of a multiline command
                # Do not activate this, otherwise, the unwrapping stops after two lines.
                if multiline_test(wrapped_line):
                    multiline = True
                else:
                    multiline = False
                self._reslist[line_num + wrapindex] = ''
            raw_line_with_comment = line  # keep comment for BEDE/LONE parsing before stripping it below
            # Source of this card as the original *physical* lines, '=' breaks
            # and all. _tag_lifetime() echoes these for post-END output.
            # Storing the glued single line instead would be unstable: the
            # continuation's leading spaces end up inside the text and
            # wrap_line() adds a fresh indent on every round-trip.
            # Trailing whitespace is dropped because it carries no meaning
            # and can push a line past wrap_line()'s 79-character limit,
            # which would make the output grow a continuation per round-trip.
            self._current_raw_line = '\n'.join(x.rstrip() for x in raw_physical_lines)
            # The current line split:
            spline: list[str] = line.split('!')[0].split()  # Ignore comments with "!"
            # The current line as string:
            line = line.upper().split('!')[0]  # Ignore comments with "!"
            word = line[:4]
            # get RESI:
            if line.startswith(('END', 'HKLF')) and self.resi:
                setattr(self.resi, 'num', 0)
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
                # Deliberately not an 'elif' chain below: closing an open
                # AFIX must not consume the dispatch, or the very HKLF line
                # that triggered it never gets parsed. The RESI and PART
                # guards above avoid 'continue' for the same reason.
            if word == 'AFIX':
                self.afix = self._assign_card(AFIX(self, spline), line_num)
            elif self._is_bede_lone_result(raw_line_with_comment, spline):
                # A BEDE/LONE bond/lone-pair electron density pseudo-atom result, e.g.:
                # L50     2    0.822405    0.640999    0.461525  !    0.235    0.164  C6
                # Kept out of self.atoms on purpose (see BedeLoneResultAtom docstring).
                post_tokens = raw_line_with_comment.split('!', 1)[1].split()
                result_atom = BedeLoneResultAtom(self)
                result_atom.parse_result_line(spline, post_tokens, list_of_lines)
                self._append_card(self.bede_lone_results, result_atom, line_num)
            elif self.is_atom_spline(word, spline):
                # A SHELXL atom:
                # F9    4    0.395366   0.177026   0.601546  21.00000   0.03231  ( 0.03248 =
                #            0.03649  -0.00522  -0.01212   0.00157 )
                a = Atom(self)
                if last_nonhydrogen_atom:
                    a.pivot = last_nonhydrogen_atom
                a.u_reference = last_unconstrained_atom
                a.parse_line(spline, list_of_lines, part=self.part, afix=cast(AFIX, self.afix), resi=self.resi)
                if not a.is_hydrogen:
                    last_nonhydrogen_atom = a
                    a.pivot = None
                if not a.qpeak and not Atom.is_riding_u(a.uvals):
                    last_unconstrained_atom = a
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
                unit = cast(UNIT, self.unit)
                if len(unit.values) != len(self.sfac_table.elements_list) and (self.debug or self.verbose):
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
                eqiv_card = EQIV(self, spline)
                if eqiv_card.number is not None:
                    if any(existing.number == eqiv_card.number for existing in self.eqiv):
                        if self.debug or self.verbose:
                            print(f'*** Duplicate EQIV {eqiv_card.id}: the same $n may '
                                  f'not appear on two EQIV instructions ***')
                    self._append_card(self.eqiv, eqiv_card, line_num)
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
                # Written with a stray unit now and then ('-100.0C'), and
                # occasionally with an esd in parentheses.
                self.temp = float(
                    spline[1].split('(')[0].rstrip('Cc\u00b0').strip()
                )
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
                # BEDE name1 name2 d a b1 b2 [!BOND! direction]
                self._append_card(self.bede_cards, BEDE(self, raw_line_with_comment.split()), line_num)
            elif word == 'LONE':
                # LONE atomName code a b1 b2 d [angle]
                self._append_card(self.lone_cards, LONE(self, raw_line_with_comment.split()), line_num)
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
                    print(f"Error in line: {line_num + 1} -> {line}")
                    if self.debug:
                        raise ParseUnknownParam(debug=self.debug, verbose=self.verbose)

    def _find_atom_insert_position(self, after: Atom | None = None) -> int:
        """
        Returns the index in ``_reslist`` at which a new atom should be inserted.

        Priority order:
        1. Directly after *after* (if given and present in ``_reslist``).
        2. Directly after the last real (non-Q-peak) :class:`~shelxfile.atoms.atom.Atom`
           in ``_reslist``.
        3. Directly before the ``HKLF`` card.
        4. One position before the end of ``_reslist`` as a last resort.
        """
        if after is not None:
            try:
                return self._reslist.index(after) + 1
            except ValueError:
                pass
        # Find last non-Q-peak atom
        last_atom_pos = None
        for i, item in enumerate(self._reslist):
            if isinstance(item, Atom) and not item.qpeak:
                last_atom_pos = i
        if last_atom_pos is not None:
            return last_atom_pos + 1
        # Fall back to before HKLF
        for i, item in enumerate(self._reslist):
            if isinstance(item, HKLF):
                return i
        return max(len(self._reslist) - 1, 0)

    def unused_atom_name(self, element: str) -> str:
        """
        Returns the first atom name of the given element that is not yet used in residue 0.

        SHELX limits atom labels to **4 characters**, so the numeric suffix is capped
        accordingly: a 1-character element (e.g. ``'C'``) allows up to ``C999``; a
        2-character element (e.g. ``'Fe'``) allows up to ``Fe99``.

        Example::

            shx.unused_atom_name('C')   # -> 'C1' if C1_0 is free, else 'C2', …
            shx.unused_atom_name('Fe')  # -> 'Fe1' if Fe1_0 is free, …

        :param element: Chemical element symbol (case-insensitive), e.g. ``'C'``, ``'Fe'``.
        :returns: An atom name string (without residue suffix) that is not present in the
                  file and is at most 4 characters long.
        :raises ValueError: if no unused name can be found within the allowed range.
        """
        prefix = element.capitalize()
        max_digits = 4 - len(prefix)
        if max_digits < 1:
            raise ValueError(
                f"Element symbol {prefix!r} is {len(prefix)} characters long; "
                f"no room for a numeric suffix within SHELX's 4-character name limit."
            )
        max_num = 10 ** max_digits - 1  # e.g. 999 for 1-char element, 99 for 2-char
        existing = set(self.atoms.nameslist)  # upper-cased "NAME_RESINUM"
        for n in range(1, max_num + 1):
            candidate = f'{prefix}{n}'
            if f'{candidate}_0'.upper() not in existing:
                return candidate
        raise ValueError(
            f"No unused atom name found for element {element!r} "
            f"(tried {prefix}1 … {prefix}{max_num})"
        )

    def add_atom(self,
                 name: str,
                 coordinates: list[float | int],
                 element: str = 'C',
                 uvals: list[float | int] | None = None,
                 part: int = 0,
                 sof: float | None = None,
                 occupancy: float | None = None,
                 fvar: int = 1,
                 resi: int = 0,
                 after: Atom | None = None,
                 coords_are_cartesian: bool = False) -> Atom:
        """
        Add a new atom to the structure and insert it into the ``_reslist`` at the
        correct position so that ``write_shelx_file`` produces a valid file.

        **Occupation encoding** — two mutually exclusive styles are supported:

        * *High-level*: pass ``occupancy`` (0.0–1.0) and optionally ``fvar``
          (free-variable number, default 1).  The raw SHELX site-occupation factor
          is computed automatically as ``fvar * 10 + occupancy``, e.g.
          ``occupancy=0.5, fvar=2`` → ``sof=20.5``.
        * *Raw SHELX*: pass ``sof`` directly, e.g. ``sof=21.0`` (fvar 2,
          occupancy 1.0).

        Passing both ``occupancy`` (or ``fvar``) **and** ``sof`` raises a
        :exc:`ValueError` to prevent ambiguous input.  If neither is given the
        atom is fully occupied on fvar 1 (``sof=11.0``).

        :param name: Atom label, e.g. ``'C1'`` (max 4 characters).  Must not
                     already exist in residue *resi*; a :exc:`ValueError` is raised
                     if it does.
        :param coordinates: Fractional (default) or Cartesian (if
                            *coords_are_cartesian* is ``True``) coordinates
                            ``[x, y, z]``.
        :param element: Chemical element symbol (default ``'C'``).  If the element
                        is not yet in the SFAC table it is added automatically
                        together with a matching UNIT count of 1.
        :param uvals: Displacement parameters.  Pass one value ``[Uiso]`` for
                      isotropic, or six values ``[U11, U22, U33, U23, U13, U12]``
                      for anisotropic.  A single-element list is silently expanded
                      to six values.  Defaults to ``[0.04, 0, 0, 0, 0, 0]``
                      (isotropic, Uiso = 0.04).
        :param part: PART number (default 0 = no disorder part).
        :param sof: Raw SHELXL site-occupation factor, e.g. ``11.0`` (fvar 1 ×
                    1.0 = fully occupied).  Mutually exclusive with *occupancy* /
                    *fvar*; defaults to ``11.0`` when neither *sof* nor *occupancy*
                    is given.
        :param occupancy: Fractional occupancy in the range 0.0–1.0.  When given,
                          the raw ``sof`` is computed as ``fvar * 10 + occupancy``.
                          Mutually exclusive with the raw *sof* parameter.
        :param fvar: Free-variable number used together with *occupancy* (default
                     1).  Mutually exclusive with the raw *sof* parameter.
        :param resi: Residue number (default 0).
        :param after: If given, the new atom is inserted directly after this
                      :class:`~shelxfile.atoms.atom.Atom` in the file.  Otherwise
                      the atom is appended after the last real atom (before
                      ``HKLF``).
        :param coords_are_cartesian: If ``True``, *coordinates* are treated as
                                     Cartesian (Å) and converted to fractional
                                     automatically.
        :returns: The new :class:`~shelxfile.atoms.atom.Atom` object.
        :raises ValueError: if *name* is already used in residue *resi*, or if
                            both *sof* and *occupancy*/*fvar* are provided.
        """
        # --- guard against mixed occupation styles ---
        high_level = occupancy is not None or fvar != 1
        raw_sof = sof is not None
        if high_level and raw_sof:
            raise ValueError(
                "Specify occupation using either 'occupancy'/'fvar' or 'sof', not both."
            )
        # --- validate name uniqueness ---
        full_name = f'{name}_{resi}'
        if self.atoms.has_atom(full_name):
            raise ValueError(f"Atom '{full_name}' already exists in the structure.")
        # --- resolve site-occupation factor ---
        if occupancy is not None:
            sof = fvar * 10 + occupancy
        elif sof is None:
            sof = 11.0  # default: fvar 1, fully occupied
        # --- normalise uvals ---
        if uvals is None:
            uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]
        elif len(uvals) == 1:
            uvals = [uvals[0], 0.0, 0.0, 0.0, 0.0, 0.0]
        # --- convert coordinates if needed ---
        if coords_are_cartesian:
            coordinates = list(cart_to_frac(coordinates, list(cast(CELL, self.cell))))
        # --- auto-register element in SFAC/UNIT ---
        if not self.sfac_table.has_element(element):
            self.sfac_table.add_element(element)
        sfac_num = self.elem2sfac(element)
        # --- build context cards ---
        part_card = PART(self, f'PART {part}'.split())
        afix_card = AFIX(self, 'AFIX 0'.split())
        # Resolve residue: look for an existing RESI card matching *resi*, else create one.
        resi_card = None
        for r in self.residues.all_residues:
            if r.residue_number == resi:
                resi_card = r
                break
        if resi_card is None:
            resi_card = RESI(self, ['RESI', str(resi)])
        # --- construct atom ---
        a = Atom(self)
        a.set_atom_parameters(
            name=name,
            sfac_num=sfac_num,
            coords=coordinates,
            part=part_card,
            afix=afix_card,
            resi=resi_card,
            site_occupation=sof,
            uvals=uvals,
            symmgen=False,
        )
        # --- insert into _reslist at the correct position ---
        insert_pos = self._find_atom_insert_position(after=after)
        self._reslist.insert(insert_pos, a)
        self.atoms.append(a)
        self.atoms._atomsdict.clear()
        return a

    #: Restraint keywords accepted by :meth:`add_restraint`, mapped to their
    #: parsing class. Kept in sync with the dispatch in :meth:`_parse_cards`.
    #: Public so GUI tooling (e.g. a restraint-type combo box) can list the
    #: supported keywords without reaching into a private attribute.
    RESTRAINT_CARD_CLASSES: dict[str, type] = {
        'DEFS': DEFS, 'NCSY': NCSY, 'ISOR': ISOR, 'FLAT': FLAT, 'BUMP': BUMP,
        'DFIX': DFIX, 'DANG': DANG, 'SADI': SADI, 'SAME': SAME, 'RIGU': RIGU,
        'SIMU': SIMU, 'DELU': DELU, 'CHIV': CHIV, 'EADP': EADP, 'EXYZ': EXYZ,
    }

    def _find_restraint_insert_position(self, after: Restraint | None = None) -> int:
        """
        Returns the index in ``_reslist`` at which a new restraint should be
        inserted.

        Priority order:
        1. Directly after *after* (if given and present in ``_reslist``).
        2. Directly after the last existing :class:`~shelxfile.shelx.cards.Restraint`
           in ``_reslist``.
        3. Directly before the first :class:`~shelxfile.atoms.atom.Atom`.
        4. One position before the end of ``_reslist`` as a last resort.
        """
        if after is not None:
            try:
                return self._reslist.index(after) + 1
            except ValueError:
                pass
        last_restraint_pos = None
        for i, item in enumerate(self._reslist):
            if isinstance(item, Restraint):
                last_restraint_pos = i
        if last_restraint_pos is not None:
            return last_restraint_pos + 1
        for i, item in enumerate(self._reslist):
            if isinstance(item, Atom):
                return i
        return max(len(self._reslist) - 1, 0)

    def add_restraint(self, text: str, after: Restraint | None = None) -> Restraint:
        """
        Parse a single restraint instruction line and insert it into the
        structure at the correct position so that ``write_shelx_file``
        produces a valid file.

        :param text: A single SHELXL restraint instruction, e.g.
                     ``'SADI 0.02 C1 C2 C1 C3'`` or ``'DFIX 1.54 C1 C2'``.
                     The residue-class suffix (e.g. ``'SADI_CCF3'``) is
                     supported. Only restraint keywords are accepted; see
                     :data:`RESTRAINT_CARD_CLASSES` for the whitelist.
        :param after: If given, the new restraint is inserted directly after
                      this :class:`~shelxfile.shelx.cards.Restraint` in the
                      file. Otherwise it is appended after the last existing
                      restraint (or before the first atom if there is none).
        :returns: The new :class:`~shelxfile.shelx.cards.Restraint` instance.
        :raises ValueError: if *text* is empty or its keyword is not a
                             supported restraint instruction.
        """
        spline = text.split()
        if not spline:
            raise ValueError('add_restraint() requires a non-empty instruction line.')
        keyword = spline[0].upper().split('_', 1)[0]
        card_cls = self.RESTRAINT_CARD_CLASSES.get(keyword)
        if card_cls is None:
            supported = ', '.join(sorted(self.RESTRAINT_CARD_CLASSES))
            raise ValueError(
                f"Unsupported restraint keyword {spline[0]!r}. Supported: {supported}"
            )
        restraint = card_cls(self, spline)
        insert_pos = self._find_restraint_insert_position(after=after)
        self._reslist.insert(insert_pos, restraint)
        self.restraints.append(restraint)
        return restraint

    def frac_to_cart(self, coordinates: list[float | int]) -> Array:
        """
        fractional to cartesian coordinates by applying the orthogonal matrix.
        """
        return cast(Array, cast(Array, self.orthogonal_matrix) * Array(coordinates))

    def __repr__(self) -> str:
        """
        Represents the shelxl object.
        """
        return self.dumps()

    def grow(self, with_qpeaks: bool = False) -> list[Atom]:
        """
        Returns a list of atoms that represent the complete molecules of the structure.
        """
        sdm = SDM(self)
        needsymm = sdm.calc_sdm()
        packed_atoms = sdm.packer(sdm, needsymm, with_qpeaks=with_qpeaks)
        return packed_atoms

    def _make_grown_names_unique(self, atoms: list[Atom]) -> None:
        """Give *atoms* SHELXL-legal, unique names, in place.

        In memory a symmetry image is labelled after its parent and the
        operation that made it (``O1>>2``), which is readable and says
        where the atom came from.  A ``.res`` file cannot carry that: a
        name is "up to 4 characters, of which the first must be a
        letter", and "the combination of atom name, PART and RESI numbers
        must be unique".

        So on the way out, anything that is not a legal name is replaced
        by a plain sequential one.  Atoms of the asymmetric unit keep
        their own labels, so the file still reads the way it was written.
        """
        from shelxfile.shelx.sdm import _unique_legal_atom_name

        def slot(atom: Atom) -> str:
            resi = atom.resi.residue_number if atom.resi else 0
            part = atom.part.n if atom.part else 0
            return f'{atom.name.upper()}_{resi}_{part}'

        def is_legal(name: str) -> bool:
            return (bool(name) and len(name) <= 4 and name[0].isalpha()
                    and name.isalnum())

        seen: set[str] = set()
        used: set[str] = {a.name.upper() for a in atoms if is_legal(a.name)}
        for atom in atoms:
            if is_legal(atom.name) and slot(atom) not in seen:
                seen.add(slot(atom))
                continue
            element = self.sfac2elem(atom.sfac_num) or 'C'
            atom._name = _unique_legal_atom_name(element, used)
            seen.add(slot(atom))

    def write_grown_file(self, filename: str | Path, with_qpeaks: bool = False) -> None:
        """Write a grown (complete molecule) .res file in P1 symmetry.

        The output file:

        - Completes all molecular fragments by applying crystal symmetry to the
          asymmetric unit (calls :meth:`grow` internally).
        - Sets the space group to P1 (``LATT -1``, no ``SYMM`` cards) so that
          viewers and SHELXL do **not** re-apply symmetry to the already-grown
          atoms — this prevents duplicate atoms and wrong bonds.
        - Preserves disorder parts via ``PART`` cards so the bond graph is
          correct (atoms in different disorder alternatives are not connected).
        - Adds a ``REM`` warning that the file is not suited for refinement.

        Restraints, ``AFIX``, ``HFIX``, ``ACTA``, ``CONF``, ``ANIS`` and other
        instruction cards that reference atoms of the original asymmetric unit
        are omitted from the output.

        :param filename: Output file path (string or :class:`~pathlib.Path`).
        :param with_qpeaks: Include Q-peaks in the grown structure (default
            ``False``).
        """
        if isinstance(filename, str):
            filename = Path(filename)

        grown_atoms = self.grow(with_qpeaks=with_qpeaks)
        real_grown = [a for a in grown_atoms if not a.qpeak] if not with_qpeaks else grown_atoms
        # Grown atoms keep the label of the atom they came from, which is
        # what makes them recognisable in a viewer. A .res file cannot: the
        # combination of name, PART and RESI has to be unique, so the
        # duplicates are relabelled here, at the point where it matters.
        self._make_grown_names_unique(real_grown)

        lines: list[str] = []

        # ── Warning comment ───────────────────────────────────────────────
        lines.append('REM This file was grown from the asymmetric unit to complete molecules.')
        lines.append('REM It has P1 symmetry and is NOT suited for refinement unless you know what you are doing.')

        # ── Iterate the original _reslist, filtering/transforming as needed ──
        atoms_inserted = False
        for num, item in enumerate(self._reslist):
            if num in self.delete_on_write:
                continue
            if item == '':
                continue

            # Skip original atoms — replaced by the grown atom list
            if isinstance(item, Atom):
                continue

            # Skip SYMM cards — P1 has no extra symmetry operations
            if isinstance(item, SYMM):
                continue

            # Skip all restraints — they reference asymmetric-unit atom names
            if isinstance(item, Restraint):
                continue

            # Skip instruction cards that only make sense in the context of
            # the original asymmetric unit
            if isinstance(item, (AFIX, PART, RESI, HFIX, ANIS)):
                continue

            # Skip refinement-specific cards not needed in a visualization file
            if isinstance(item, (ACTA, CONF)):
                continue

            # Skip include-file references (+filename) — they are local to the
            # original file's directory and would not resolve for the grown output
            if isinstance(item, str) and item.startswith('+'):
                continue

            # Replace LATT with LATT -1 (primitive, non-centrosymmetric = P1)
            if isinstance(item, LATT):
                lines.append('LATT -1')
                continue

            # Before the HKLF card: insert grown atoms with PART grouping
            if isinstance(item, HKLF) and not atoms_inserted:
                atoms_inserted = True
                current_part = 0
                for atom in real_grown:
                    atom_part = atom.part.n
                    if atom_part != current_part:
                        lines.append(f'PART {atom_part}')
                        current_part = atom_part
                    lines.append(wrap_line(str(atom)))
                if current_part != 0:
                    lines.append('PART 0')
                lines.append('')

            line_str = str(item)
            if not line_str.strip():
                continue
            line_str = '\n'.join([wrap_line(x) for x in line_str.split('\n')])
            lines.append(line_str)

        with open(filename, 'w') as f:
            f.write('\n'.join(lines) + '\n')

        if self.verbose or self.debug:
            print(f'*** Grown file written to {filename.resolve()} ***')

    def pack(self, with_qpeaks: bool = False) -> list[Atom]:
        """Returns a list of atoms representing the packed unit cell.

        Applies all symmetry operations to the asymmetric unit and folds every
        position back into [0, 1) fractional coordinates, removing duplicates.
        Unlike :meth:`grow`, this does not stitch molecular fragments together —
        it simply fills one unit cell.

        :param with_qpeaks: include Q-peaks (difference-map peaks) in the output.
        :returns: List of :class:`~shelxfile.atoms.atom.Atom` objects.
        """
        sdm = SDM(self)
        return sdm.pack_unit_cell(with_qpeaks=with_qpeaks)

    def refine(self, cycles: int | None = None, backup_before: bool = True) -> bool:
        if self.resfile:
            filen = self.resfile.stem
            if cycles is not None:
                cast(LSCycles, self.cycles).number = cycles
            ref = ShelxlRefine(self, self.resfile)
            ref.remove_acta_card(cast(ACTA, self.acta))
            self.write_shelx_file(filen + '.ins')
            ref.run_shelxl(backup_before=backup_before)
            self.reload()
            ref.restore_acta_card()
            # self.write_shelx_file(filen + '.res')
            return True
        return False

    def refine_weight_convergence(self, stop_after: int = 10) -> bool:
        """
        Tries to refine weigting sheme from SHELXL until it converged (self.weight_difference() is zero) or
        stopt_after cycles are reached.
        """
        for _ in range(stop_after):
            difference = cast(WGHT, self.wght).difference()
            print("Weighting difference = {} {}".format(*difference))
            if self._weight_converged(difference):
                return True
            else:
                self.update_weight()
                self.refine(9)
        print("Maximum number of refinement cycles reached, but no WGHT convergence.")
        return False

    def _weight_converged(self, diff: list[float]) -> bool:
        return diff == [0.0, 0.0]

    def _append_card(self, obj, card, line_num: int) -> Command:
        """
        Appends SHELX card to an object list, e.g. self.restraints and
        assigns the line_num in reslist with the card instance.
        """
        self._tag_lifetime(card, line_num)
        obj.append(card)
        self._reslist[line_num] = card
        return card

    def _assign_card(self, card, line_num: int):
        self._tag_lifetime(card, line_num)
        self._reslist[line_num] = card
        return card

    def _tag_lifetime(self, card, line_num: int) -> None:
        """Mark *card* as live input or inert post-``END`` output.

        Cards after ``END`` keep their **original source line** so they
        round-trip verbatim (plan item D-10).  Parsing normalises
        whitespace via ``' '.join(spline)``, which would silently reflow
        the generated ``SADI`` block that SHELXL wrote out; that block is
        a record of a refinement that already happened and must not be
        rewritten.
        """
        if not self.end:
            card.lifetime = CardLifetime.INPUT
            return
        card.lifetime = CardLifetime.POST_END_OUTPUT
        raw = self._current_raw_line
        if not isinstance(raw, str) or not raw:
            return
        # Restraint renders from .textline, Command from ._textline.
        if hasattr(card, 'textline'):
            card.textline = raw
        if hasattr(card, '_textline'):
            card._textline = raw

    @staticmethod
    def is_atom(atomline: str) -> bool:
        """
        Returns True is line contains an atom.
        """
        # no empty line, not in cards and not space at start:
        if atomline[:4].upper() not in SHX_CARDS:  # exclude all non-atom cards
            spline: list[str] = atomline.split()
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
    def is_atom_spline(word4: str, spline: list[str]) -> bool:
        """
        Fast atom check using a pre-split, pre-uppercased spline.

        Avoids re-splitting the line and uppercasing it again — the caller
        (``_parse_cards``) has already computed both ``word = line[:4]`` and
        ``spline`` before reaching this check.

        :param word4: The first four characters of the line, upper-cased.
        :param spline: The line already split on whitespace (comments stripped).
        """
        if word4 in SHX_CARDS:
            return False
        if len(spline) < 5:
            return False
        if '.' in spline[1]:
            return False
        # Inline coordinate-realism check (avoids genexpr + slice + function call
        # for the overwhelmingly common case of plain fractional coordinates):
        try:
            for token in (spline[2], spline[3], spline[4]):
                value = float(token)
                if value > MAX_PLAIN_COORDINATE and not _decodes_to_a_coordinate(value):
                    return False
        except ValueError:
            return False
        if len(spline) > 5 and spline[5] == '!':
            return False
        return True

    @staticmethod
    def _coordinates_are_unrealistic(spline: list[str]) -> bool:
        return any(float(y) > MAX_PLAIN_COORDINATE
                   and not _decodes_to_a_coordinate(float(y))
                   for y in spline[2:5])

    @staticmethod
    def _is_bede_lone_result(raw_line: str, spline: list[str]) -> bool:
        """
        Detects a BEDE/LONE bond/lone-pair electron-density pseudo-atom
        result line, e.g.::

            L50     2    0.822405    0.640999    0.461525  !    0.235    0.164  C6

        Unlike every real SHELXL atom line (which always has at least 6
        tokens before any comment: name, sfac, x, y, z, occ[, U...]), these
        pseudo-atom result lines have exactly 5 tokens (name, sfac, x, y, z)
        followed by a ``!`` comment carrying the BEDE/LONE ``b1``, ``b2`` and
        owner-atom-name metadata.

        :param raw_line: the line before comment-stripping (still contains '!').
        :param spline: the same line already split on whitespace with any
            comment removed.
        """
        if '!' not in raw_line:
            return False
        if len(spline) != 5:
            return False
        if '.' in spline[1]:
            return False
        try:
            if any(abs(float(v)) > 4.0 for v in spline[2:5]):
                return False
        except ValueError:
            return False
        return True

    def to_cif(self, filename: str | None = None, template: str | None = None) -> None:
        """
        Writes a CIF file from the ShelxFile object.
        """
        if not filename:
            filename = cast(Path, self.resfile).stem + '.cif'
        CifFile(self, template).write_cif(Path(filename))

    def get_bede_for_atom(self, atom_name: str) -> list[BEDE]:
        """
        Returns all BEDE cards where the given atom name occurs as name1 or name2.
        """
        atom_name = atom_name.upper()
        return [b for b in self.bede_cards if atom_name in (b.name1, b.name2)]

    def get_lone_for_atom(self, atom_name: str) -> list[LONE]:
        """
        Returns all LONE cards for the given atom name.
        """
        atom_name = atom_name.upper()
        return [lo for lo in self.lone_cards if lo.name == atom_name]

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

    def replace_line(self, obj: ResListEntry, new_line: str) -> None:
        """
        Replaces a single line in the res file with new_line.
        """
        self._reslist[self.index_of(obj)] = new_line

    def eqiv_by_id(self, eqiv_id: str | int) -> EQIV | None:
        """The ``EQIV`` card with the given ``$n`` id, or ``None``.

        Accepts either ``'$2'`` or ``2``.
        """
        if isinstance(eqiv_id, str):
            number = int(eqiv_id.lstrip('$')) if eqiv_id.lstrip('$').isdigit() else None
        else:
            number = eqiv_id
        if number is None:
            return None
        for card in self.eqiv:
            if card.number == number:
                return card
        return None

    def index_of(self, obj: ResListEntry) -> int:
        """Position of *obj* in ``_reslist``, matched by **identity**.

        ``list.index()`` uses ``__eq__``, and :meth:`Atom.__eq__` compares
        the serialised line.  Two atoms that render identically would
        therefore resolve to the same index and the wrong line would be
        edited or deleted.  SHELXL guarantees that the combination of atom
        name, ``PART`` and ``RESI`` is unique, but the rendered string is
        not, so identity is the only safe key here.

        :raises ValueError: if *obj* is not present.
        """
        for num, item in enumerate(self._reslist):
            if item is obj:
                return num
        # Fall back to equality for plain strings, which have no identity
        # of their own once they have been copied around.
        if isinstance(obj, str):
            return self._reslist.index(obj)
        # Deliberately not %r: repr() of an Atom asks for its atomid, which
        # calls back into index_of() and would recurse forever.
        raise ValueError(f'{type(obj).__name__} instance is not in the reslist')

    def remove_from_reslist(self, obj: ResListEntry) -> int:
        """Remove *obj* from ``_reslist`` by identity and keep bookkeeping sane.

        Every deletion path must go through here.  Besides dropping the
        entry it re-maps :attr:`delete_on_write`, which stores *indices*:
        without the shift, deleting any line would silently suppress the
        wrong line on the next write.

        :returns: the index the object occupied.
        :raises ValueError: if *obj* is not present.
        """
        index = self.index_of(obj)
        del self._reslist[index]
        self._shift_delete_on_write(index)
        self.touch()
        return index

    def _shift_delete_on_write(self, removed_index: int) -> None:
        """Keep :attr:`delete_on_write` pointing at the same lines.

        Indices above *removed_index* move down by one; an index equal to
        it refers to a line that no longer exists and is dropped.
        """
        if not self.delete_on_write:
            return
        self.delete_on_write = {
            num - 1 if num > removed_index else num
            for num in self.delete_on_write
            if num != removed_index
        }

    @property
    def sum_formula(self) -> str:
        """
        The sum formula of the structure in regard to the UNIT instruction.
        """
        formstring = ''
        formula_weight = 0.0
        try:
            unit = cast(UNIT, self.unit)
            val = unit.values
            eli = self.sfac_table.elements_list
        except AttributeError:
            return ''
        if len(val) == len(eli):
            for el, num in zip(self.sfac_table.elements_list, unit.values):
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
            wght = cast(WGHT, self.wght)
            wght_suggested = cast(WGHT, self.wght_suggested)
            wght.a = wght_suggested.a
            wght.b = wght_suggested.b
            wght.c = wght_suggested.c
            wght.d = wght_suggested.d
            wght.e = wght_suggested.e
            wght.f = wght_suggested.f
        except AttributeError:
            return

    def insert_anis(self, atoms: str = '', residue: str = '') -> None:
        """
        Inserts ANIS into a results file for refinement of anisotropic displacement parameters.

        TODO: implement ANIS n
        @param atoms: Specify secific atoms or wildcards of atoms.
        @param residue: Specify a residue like ANIS_ABC with residue='ABC' or ANIS_* (residue='*')
        """
        unit = cast(UNIT, self.unit)
        if atoms:
            self.add_line(unit.position, f'ANIS{"_" if residue else ""}{residue} {atoms}')
        else:
            self.add_line(unit.position, 'ANIS')

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

    def sum_formula_exact_as_dict(self) -> dict[str, float]:
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

    def insert_frag_fend_entry(self, dbatoms: list, cell: list) -> None:
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

    def _get_residuals(self, spline: list[str], line: str) -> None:
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

    def _get_space_group(self, spline: list[str]) -> None:
        try:
            self.space_group = spline[3]
        except(IndexError, ValueError):
            pass

    def _get_goof(self, spline: list[str]) -> None:
        with suppress(IndexError, ValueError):
            self.goof = float(spline[8].split(',')[0])
            self.rgoof = float(spline[12].split(',')[0])

    def _get_peak_hole(self, spline: list[str]) -> None:
        # REM Highest difference peak  0.407,  deepest hole -0.691,  1-sigma level  0.073
        with suppress(IndexError, ValueError):
            self.highest_peak = float(spline[4].split(",")[0])
            self.deepest_hole = float(spline[7].split(",")[0])

    def _get_params_and_restraints(self, spline: list[str]) -> None:
        with suppress(IndexError):
            self.parameters = int(spline[1])
            if self.data and self.parameters:
                self.dat_to_param = float(self.data) / float(self.parameters)
        with suppress(IndexError, ValueError):
            self.num_restraints = int(spline[-2])

    def _get_wr2(self, spline: list[str]) -> None:
        with suppress(IndexError, ValueError):
            self.wr2 = float(spline[3].split(",")[0])

    def _get_r1(self, spline: list[str]) -> None:
        with suppress(IndexError, ValueError):
            self.R1 = float(spline[3])
        with suppress(IndexError, ValueError):
            self.data = int(spline[-2])


if __name__ == "__main__":
    print(Path('.').resolve())
    # file = r'../shelxfile/tests/resources/p21c.res'
    # file = r'D:\_DEV\GitHub\ShelXFile\tests\resources\test_bedelone.res'
    file = r'/Users/daniel/Documents/GitHub/FinalCif/tests/examples/Esser_JW367_0m-finalcif.res'
    shx = Shelxfile(debug=True)
    shx.read_file(file)
    # print(shx.atoms)
    # print(shx.sum_formula_exact)
    # print(shx.sum_formula)
    # print(shx.sum_formula_exact_as_dict())
    # print(shx.restraints)
    # print(shx.atoms.nameslist)
    print(shx)
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
