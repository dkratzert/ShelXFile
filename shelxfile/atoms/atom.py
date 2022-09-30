from contextlib import suppress
from typing import Union, List, Tuple

import numpy
from numpy import trace

with suppress(Exception):
    from shelxfile import Shelxfile
from shelxfile.misc.dsrmath import atomic_distance, Array, Matrix
from shelxfile.misc.elements import get_atomic_number, get_radius_from_element
from shelxfile.misc.misc import split_fvar_and_parameter, DEBUG, ParseSyntaxError, frac_to_cart, ParseUnknownParam, \
    VERBOSE
from shelxfile.shelx.cards import PART, AFIX, RESI, CELL, Restraints


class Atom():
    """
    An object holding all Properties of a shelxl atom plus some extra information like
    kartesian coordinates and element type.

    atomname sfac x y z sof[11] U[0.05] or U11 U22 U33 U23 U13 U12
    """
    #                name    sfac     x         y        z       occ      u11      u22 ...
    _anisatomstr = '{:<5s} {:>2}{:>12.6f}{:>12.6f}{:>12.6f}{:>12.5f}{:>11.5f}{:>11.5f}' \
                   ' {:>12.5f}{:>11.5f}{:>11.5f}{:>11.5f}'  # Line wrap is handled during file write.
    #               name    sfac     x         y         z         occ      u11
    _isoatomstr = '{:<5s} {:<2}{:>10.6f}  {:>10.6f}  {:>9.6f}  {:>9.5f}  {:>9.5f}'
    _qpeakstr = '{:<5s} {:<2}{:>8.4f}  {:>8.4f}  {:>8.4f}  {:>9.5f}  {:<9.2f} {:<9.2f}'
    _fragatomstr = '{:<5s} {:>10.6f}  {:>10.6f}  {:>9.6f}'

    def __init__(self, shx: 'Shelxfile') -> None:
        self.shx = shx
        self._cell: CELL = shx.cell
        self.sfac_num: int = 1
        self.resi: Union[RESI, None] = None
        self.part: PART = PART(shx, ['PART', '0'])
        self.afix: Union[AFIX, None] = None
        self.name = 'name'  # Name without residue number like "C1"
        # Site occupation factor including free variable like 31.0
        self.sof = 11.0
        # fractional coordinates:
        self.x: float = 0.0
        self.y: float = 0.0
        self.z: float = 0.0
        # cartesian coordinates:
        self.xc: float = 0.0
        self.yc: float = 0.0
        self.zc: float = 0.0
        self.qpeak: bool = False
        self.peak_height: float = 0.0
        self.uvals: List[float] = [0.04, 0.0, 0.0, 0.0, 0.0]  # [U] or [U11 U22 U33 U23 U13 U12]
        self.uvals_orig: List[float] = [0.04, 0.0, 0.0, 0.0, 0.0]
        self.frag_atom: bool = False
        self.restraints: List[Restraints] = []
        self._line_numbers = None
        self._occupancy: float = 1.0
        self.molindex: int = 0
        # Indicates if this atom is generated by symmetry:
        self.symmgen: bool = False

    def __eq__(self, other):
        return self.__str__() == other.__str__()

    @property
    def atomid(self) -> int:
        if self.symmgen:
            return 0
        try:
            return self.shx._reslist.index(self)
        except ValueError:
            return 0

    @property
    def fullname(self) -> str:
        return self.name + '_' + str(self.resinum)  # Name including residue nimber like "C1_2"

    @property
    def resiclass(self) -> str:
        return self.resi.residue_class

    @property
    def resinum(self) -> int:
        return self.resi.residue_number

    @property
    def chain_id(self) -> int:
        return self.resi.chain_id

    @property
    def fvar(self) -> int:
        # Be aware: fvar can be negative!
        fvar, _ = split_fvar_and_parameter(self.sof)
        return fvar

    @property
    def occupancy(self) -> float:
        # Only the occupancy of the atom like 0.5 (without the free variable)
        _, occ = split_fvar_and_parameter(self.sof)
        # Fractional occupancy:
        if abs(self.fvar) == 1:
            return occ
        else:
            if occ > 0:
                occ = self._get_positive_occupancy(occ)
            else:
                occ = self._get_negative_occupancy(occ)
        return occ

    def _get_negative_occupancy(self, occ):
        try:
            occ = 1 + (self.shx.fvars[self.fvar] * occ)
        except IndexError:
            occ = 1.0
            if DEBUG:
                raise ParseSyntaxError
            if VERBOSE:
                print(f'*** Could not get occupancy of free variable {self.fvar} ***')
        return occ

    def _get_positive_occupancy(self, occ):
        try:
            occ = self.shx.fvars[self.fvar] * occ
        except IndexError:
            occ = 1.0  # Happens if the self.fvar is not defined
            if DEBUG:
                raise ParseSyntaxError
            if VERBOSE:
                print(f'*** Could not get occupancy of free variable {self.fvar} ***')
        return occ

    @occupancy.setter
    def occupancy(self, occ: float):
        self._occupancy = occ

    @property
    def is_hydrogen(self) -> bool:
        """
        Returns True if the current atom is a hydrogen isotope.
        """
        if self.element in {'H', 'D', 'T'}:
            return True
        else:
            return False

    @property
    def ishydrogen(self) -> bool:
        return self.is_hydrogen

    def set_atom_parameters(self, name: str = 'C', sfac_num: int = 1, coords: List[float] = None, part: PART = None,
                            afix: AFIX = None, resi: RESI = None, site_occupation: float = 11.0,
                            uvals: (list, tuple) = None,
                            symmgen: bool = True):
        """
        Sets atom properties manually if not parsed from a SHELXL file.
        """
        self.name = name
        self.sfac_num = sfac_num
        self.frac_coords = coords
        self.x, self.y, self.z = coords[0], coords[1], coords[2]
        self.xc, self.yc, self.zc = frac_to_cart(self.frac_coords, list(self._cell))
        self.part = part
        self.afix = afix
        self.resi = resi
        self.sof = site_occupation
        self.uvals = uvals
        self.symmgen = symmgen

    def set_uvals(self, uvals: List[float]) -> None:
        """
        Sets u values and checks if a free variable was used.
        """
        self.uvals = uvals
        if uvals[2] == 0.0:  # 0 is Uiso and 1 q-peak hight
            for n, u in enumerate(uvals):
                if abs(u) > 4.0:
                    fvar, uval = split_fvar_and_parameter(u)
                    # self.uvals[n] = uval
                    self.shx.fvars.set_fvar_usage(fvar)
        else:
            if abs(uvals[0]) > 4.0:
                fvar, uval = split_fvar_and_parameter(uvals[0])
                self.shx.fvars.set_fvar_usage(fvar)

    def parse_line(self, atline: List, list_of_lines: List, part: PART, afix: AFIX, resi: RESI) -> None:
        """
        Parses the text line of an atom from SHELXL to initialize the atom parameters.
        """
        uvals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.name = atline[0][:4]  # Atom names are limited to 4 characters
        for n, u in enumerate(atline[6:12]):
            uvals[n] = float(u)
        self.uvals_orig = uvals[:]
        self.set_uvals(uvals)
        self._line_numbers = list_of_lines
        self.part = part
        self.afix = afix
        self.resi = resi
        self._get_part_and_occupation(atline)
        self.x, self.y, self.z = self._get_atom_coordinates(atline)
        self.xc, self.yc, self.zc = self._cell.o * Array(self.frac_coords)
        if abs(self.uvals[1]) > 0.0 and self.uvals[2] == 0.0 and self.shx.hklf:  # qpeaks are always behind hklf
            self.peak_height = uvals[0]
            self.qpeak = True
        if self.shx.end:  # After 'END' can only be Q-peaks!
            self.qpeak = True
        self.sfac_num = int(atline[1])
        self.shx.fvars.set_fvar_usage(self.fvar)
        # U(star) = N * U(cif) * N.T
        # U(cart) = A * U(star) * A.T
        # U(star) = R * U(star) * R^t
        # U(cif) = N^-1 * U(star) * (N^-1).T
        # U(star) = A^-1 * U(cart) * A^-1.T
        self.u_cart = self._cell.o * Array(self.uvals)
        print(self.name, uvals)
        self.Ucif = Matrix(((uvals[0], uvals[1], uvals[3]),
                            (uvals[1], uvals[2], uvals[4]),
                            (uvals[3], uvals[4], uvals[5])))
        self.Ustar = self._cell.N * self.Ucif * self._cell.N.T
        self.Ucart = self._cell.o * self.Ustar * self._cell.o.T
        #print(self.Ustar)
        print(self.Ucart.trace/3, '###1')

    def _get_part_and_occupation(self, atline: List[str]) -> None:
        # TODO: test all variants of PART and AFIX sof combinations:
        if self.part.sof != 11.0:
            if self.afix and self.afix.sof:  # handles position of afix and part:
                if self.afix.index > self.part.index:
                    self.sof = self.afix.sof
            else:
                self.sof = self.part.sof
        elif self.afix and self.afix.sof:
            if self.part.sof != 11.0:
                if self.part.index > self.afix.index:
                    self.sof = self.part.sof
            else:
                self.sof = self.afix.sof
        else:
            self.sof = float(atline[5])

    def _get_atom_coordinates(self, atline: List[str]) -> Tuple[float, float, float]:
        try:
            x, y, z = [float(x) for x in atline[2:5]]
        except ValueError as e:
            if DEBUG or VERBOSE:
                print(e, 'Line:', self._line_numbers[-1])
                raise ParseUnknownParam
            else:
                x, y, z = 0.1, 0.1, 0.1
        if abs(x) > 4:
            fvar, x = split_fvar_and_parameter(x)
            self.shx.fvars.set_fvar_usage(fvar)
        if abs(y) > 4:
            fvar, y = split_fvar_and_parameter(y)
            self.shx.fvars.set_fvar_usage(fvar)
        if abs(z) > 4:
            fvar, z = split_fvar_and_parameter(z)
            self.shx.fvars.set_fvar_usage(fvar)
        return x, y, z

    @property
    def element(self) -> str:
        """
        Chemical element character
        """
        return self.shx.sfac2elem(self.sfac_num).capitalize()

    @property
    def an(self) -> int:
        return get_atomic_number(self.element)

    @element.setter
    def element(self, new_element: str) -> None:
        """
        Sets the element type of atom.
        """
        self.sfac_num = self.shx.elem2sfac(new_element)

    @property
    def radius(self) -> float:
        """
        Returns the atomic covalence radius in angstrom.
        """
        return get_radius_from_element(self.element)

    def __iter__(self):
        for x in self.__repr__().split():
            yield x

    def __repr__(self) -> str:
        return 'Atom ID: {}'.format(self.atomid)

    def __str__(self) -> str:
        """
        Returns a text line of the Atom with SHELXL syntax.
        :return: SHELX-formated atom string
        """
        if self.afix and self.shx.frag:
            # An atom from a FRAG/FEND instruction
            return Atom._fragatomstr.format(self.name, self.x, self.y, self.z)
        else:
            if self.uvals[0] > 0 and self.uvals[1] != 0.0 and self.uvals[2] != 0.0 and self.uvals[3] != 0.0 \
                    and self.uvals[4] != 0.0:
                # anisotropic atom
                try:
                    return Atom._anisatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof,
                                                    *self.uvals)
                except IndexError:
                    return 'REM Error in U values.'
            else:
                # isotropic atom:
                if self.qpeak:
                    return Atom._qpeakstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, 0.04,
                                                 self.peak_height)
                try:
                    return Atom._isoatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof,
                                                   *self.uvals)
                except IndexError:
                    return Atom._isoatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, 0.04)

    @property
    def index(self) -> int:
        # The position in the res file as index number (starting from 0).
        return self.shx.index_of(self)

    @property
    def frac_coords(self) -> tuple:
        return self.x, self.y, self.z

    @frac_coords.setter
    def frac_coords(self, coords: List[float]):
        self.x, self.y, self.z = coords

    @property
    def cart_coords(self) -> Tuple[float, float, float]:
        return self.xc, self.yc, self.zc

    def delete(self):
        """
        Delete atom(s) in the file.
        """
        del self.shx.atoms[self.index]

    def to_isotropic(self) -> None:
        """
        Makes the current atom isotropic.
        """
        self.uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]

    def find_atoms_around(self, dist=1.2, only_part=0) -> List['Atom']:
        """
        Finds atoms around the current atom.
        """
        found = []
        for at in self.shx.atoms:
            if atomic_distance([self.x, self.y, self.z], [at.x, at.y, at.z], self._cell) < dist \
                    and self != at and at.part.n == only_part and not at.qpeak:
                # only in special part and no q-peaks:
                found.append(at)
        return found

    def get_pivot_atom(self) -> Union['Atom', None]:
        """
        Returns the pivot atom (C1) of a riding hydrogen atom e.g. (H1, H2, or H3).
           /H1
        -C1--H2
           \H3
        """
        pivots = self.find_atoms_around(dist=1.2)
        return pivots[0] if pivots and self.afix.mn else None
