"""Structural edits: moving atoms, disorder parts, duplication, naming.

These are the model-level primitives behind the matching
:class:`~shelxfile.edit.ShelxDocument` methods.  They mutate ``_reslist``
directly, which is why they live in the edit layer, and they do **no**
notification or history bookkeeping -- the document wraps every call in a
:meth:`~shelxfile.edit.ShelxDocument.batch`.

Two SHELXL rules shape most of the code here:

* ``PART`` and ``AFIX`` are *brackets*: an instruction applies to every atom
  that follows it until the next one.  Putting atoms into a new part means
  inserting ``PART n`` / ``PART 0`` around them, and those cards must never
  land inside an ``AFIX`` group, or the group's riding atoms would end up in
  a different part than the atom they ride on.
* A riding U value such as ``-1.2`` means *"1.2 times U(eq) of the previous
  atom not constrained this way"*, so the **order** of atoms matters.  A
  duplicated moiety is therefore written as one block in the original file
  order, with each ``AFIX`` card copied in front of the atoms it governed.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, Sequence

from shelxfile.edit.brackets import BracketResolver

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.atoms.atom import Atom

__all__ = [
    'MAX_ATOM_NAME_LENGTH',
    'add_free_variable',
    'afix_closure_problem',
    'assign_part',
    'duplicate_atoms',
    'free_atom_name',
    'move_atom',
    'split_name',
]

#: SHELXL atom names are "up to 4 characters, of which the first must be a
#: letter".
MAX_ATOM_NAME_LENGTH = 4

#: A coordinate change smaller than this is not a move, so a fixed
#: coordinate code survives it.
_COORDINATE_TOLERANCE = 1e-6


# ---------------------------------------------------------------- helpers

def _file_order(shx: Shelxfile, atoms: Iterable[Atom]) -> list[Atom]:
    """*atoms* without duplicates, sorted by their position in the file."""
    unique: dict[int, Atom] = {}
    for atom in atoms:
        unique[id(atom)] = atom
    return sorted(unique.values(), key=shx.index_of)


def _check_real_atoms(atoms: Sequence[Atom]) -> None:
    for atom in atoms:
        if atom.symmgen:
            raise ValueError(f'{atom.name} is a symmetry-generated atom, not a line in the file')
        if atom.qpeak:
            raise ValueError(f'{atom.name} is a Q-peak')


def _is_blank(item: object) -> bool:
    return isinstance(item, str) and not item.strip()


def _is_comment(item: object) -> bool:
    """A ``REM`` line, which is a comment and never part of a bracket."""
    if type(item).__name__ == 'REM':
        return True
    return isinstance(item, str) and item.strip().upper().startswith('REM')


def _is_skippable(item: object) -> bool:
    return _is_blank(item) or _is_comment(item)


def _is_afix(item: object) -> bool:
    return type(item).__name__ == 'AFIX'


def _is_part(item: object) -> bool:
    return type(item).__name__ == 'PART'


def _next_meaningful(reslist: list, index: int) -> int:
    """First index at or after *index* that is not a blank or comment line.

    ``REM`` lines are skipped as well as blanks: SHELXL treats them as
    comments, so a ``REM`` between the last riding atom and its closing
    ``AFIX 0`` does not end the group.  Files written with an embedded
    ``REM <hkl>`` block really do look like that.
    """
    while index < len(reslist) and _is_skippable(reslist[index]):
        index += 1
    return index


def _previous_meaningful(reslist: list, index: int) -> int:
    """Last index at or before *index* that is not a blank or comment line."""
    while index >= 0 and _is_skippable(reslist[index]):
        index -= 1
    return index


def afix_closure_problem(shx: Shelxfile, atoms: Iterable[Atom]) -> str | None:
    """Why *atoms* cannot be moved around as a unit, or ``None``.

    Every ``AFIX`` group touched must be taken whole, together with the
    atom it rides on, and an atom that carries riders must bring them
    along: otherwise a ``PART`` card or a duplicated block would separate
    riders from their pivot.
    """
    selected = {id(atom) for atom in atoms}
    resolver = BracketResolver(shx)
    for scope in resolver.afix_scopes():
        if scope.is_reset or not scope.members:
            continue
        inside = [member for member in scope.members if id(member) in selected]
        pivot = resolver.pivot_of(scope)
        pivot_selected = pivot is not None and id(pivot) in selected
        mn = getattr(scope.card, 'mn', '?')
        if inside and len(inside) != len(scope.members):
            names = ', '.join(member.name for member in scope.members)
            return f'the AFIX {mn} group ({names}) is only partly selected'
        if inside and pivot is not None and not pivot_selected:
            return (f'the AFIX {mn} group of {inside[0].name} is selected '
                    f'without the atom it rides on ({pivot.name})')
        if pivot_selected and not inside:
            return (f'{pivot.name} is selected without its riding AFIX {mn} '
                    f'group ({scope.members[0].name}, ...)')
    return None


# ----------------------------------------------------------------- naming

def _name_is_free(shx: Shelxfile, name: str, resinum: int, taken: set[str]) -> bool:
    if not name or not name[0].isalpha() or len(name) > MAX_ATOM_NAME_LENGTH:
        return False
    if name.upper() in taken:
        return False
    return not shx.atoms.has_atom(f'{name}_{resinum}')


def free_atom_name(shx: Shelxfile, element: str, resinum: int = 0,
                   taken: Iterable[str] = ()) -> str:
    """The first ``<element><number>`` name unused in residue *resinum*.

    :param taken: Further names (case-insensitive) to treat as used, for
        callers that are about to create several atoms at once.
    :raises ValueError: when the 4-character limit leaves no free name.
    """
    prefix = element.capitalize()
    reserved = {name.upper() for name in taken}
    digits = MAX_ATOM_NAME_LENGTH - len(prefix)
    if digits < 1:
        raise ValueError(f'Element symbol {prefix!r} leaves no room for a number')
    for number in range(1, 10 ** digits):
        candidate = f'{prefix}{number}'
        if _name_is_free(shx, candidate, resinum, reserved):
            return candidate
    raise ValueError(f'No unused atom name left for element {element!r}')


def _family_name(shx: Shelxfile, atom: Atom, reserved: set[str]) -> str | None:
    """``H10A`` → ``H10D`` style: same stem, the next free trailing letter.

    This is how SHELXL itself names the second set of a split methyl group
    (``AFIX 127``: ``H..A`` to ``H..F``), and it keeps labels readable where
    appending a suffix would exceed four characters.
    """
    stem = atom.name.rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz')
    if not stem or not stem[0].isalpha():
        stem = atom.name
    if len(stem) >= MAX_ATOM_NAME_LENGTH:
        return None
    for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
        candidate = f'{stem}{letter}'
        if candidate.upper() == atom.name.upper():
            continue
        if _name_is_free(shx, candidate, atom.resinum, reserved):
            return candidate
    return None


def split_name(shx: Shelxfile, atom: Atom, suffix: str,
               taken: Iterable[str] = (), *, keep_own_name: bool = False) -> str:
    """Name for one half of a split atom: its name plus *suffix*.

    ``C1`` becomes ``C1A`` or ``C1B``.  When that would exceed SHELXL's four
    characters or is already used in the atom's residue:

    * with *keep_own_name* the atom's current name is returned (the half
      that stays where it is keeps its label, and the part number tells the
      two halves apart);
    * otherwise the same stem with the next free trailing letter is tried
      (``H10A`` → ``H10D``), and finally a free ``<element><number>`` name.
    """
    reserved = {name.upper() for name in taken}
    candidate = f'{atom.name}{suffix}'
    if _name_is_free(shx, candidate, atom.resinum, reserved):
        return candidate
    if keep_own_name:
        return atom.name
    family = _family_name(shx, atom, reserved)
    if family is not None:
        return family
    return free_atom_name(shx, atom.element, atom.resinum, reserved)


# ----------------------------------------------------------------- moving

def move_atom(shx: Shelxfile, atom: Atom, frac: Sequence[float]) -> list[str]:
    """Put *atom* at fractional coordinates *frac*.

    A coordinate held by a ``10*m + p`` code (fixed, or tied to a free
    variable) that actually changes loses its code: keeping the code would
    re-impose the old constraint on the new value.

    :returns: One message per released coordinate constraint.
    """
    from shelxfile.misc.misc import frac_to_cart_fast

    _check_real_atoms([atom])
    messages: list[str] = []
    old = (atom.x, atom.y, atom.z)
    codes = list(atom._coord_codes)
    for axis, (before, after) in enumerate(zip(old, frac)):
        if codes[axis] and abs(float(after) - before) > _COORDINATE_TOLERANCE:
            messages.append(
                f'{atom.name}: released the constraint on {"xyz"[axis]} '
                f'(code {codes[axis]}) because the coordinate moved'
            )
            codes[axis] = 0
    atom._coord_codes = (codes[0], codes[1], codes[2])
    atom.x, atom.y, atom.z = (float(value) for value in frac)
    atom.xc, atom.yc, atom.zc = frac_to_cart_fast(atom.x, atom.y, atom.z, atom._cell)
    shx.touch()
    return messages


# --------------------------------------------------------- free variables

def add_free_variable(shx: Shelxfile, value: float) -> int:
    """Append a free variable with starting *value*; return its number.

    SHELXL numbers free variables from 1 (the overall scale factor), so the
    new one is ``len(FVAR) + 1``.
    """
    from shelxfile.atoms.atom import Atom as AtomClass
    from shelxfile.shelx.cards import FVAR

    number = len(shx.fvars) + 1
    if number > 99:
        raise ValueError('SHELXL allows only 99 free variables')
    if number == 1:
        # No FVAR yet: the scale factor has to come first.
        shx.fvars.append(FVAR(1, 1.0))
        number = 2
    shx.fvars.append(FVAR(number, float(value)))
    if not any(item is shx.fvars for item in shx._reslist):
        position = next(
            (i for i, item in enumerate(shx._reslist) if isinstance(item, AtomClass)),
            len(shx._reslist),
        )
        shx.insert_into_reslist(position, shx.fvars)
    shx.touch()
    return number


# ------------------------------------------------------------------ parts

def _per_atom(value: float | Sequence[float], count: int, what: str) -> list[float]:
    """*value* repeated *count* times, or checked to be *count* long."""
    if isinstance(value, (int, float)):
        return [float(value)] * count
    values = [float(v) for v in value]
    if len(values) != count:
        raise ValueError(f'{what} must be one number or one per atom')
    return values


def assign_part(shx: Shelxfile, atoms: Iterable[Atom], part: int,
                sof: float | Sequence[float]) -> list:
    """Move *atoms* from ``PART 0`` into ``PART part`` with occupancy *sof*.

    Each contiguous run of the atoms is wrapped in ``PART part`` ...
    ``PART 0``.  A run may span ``AFIX`` cards (riders travel with their
    pivot), and the brackets are placed outside any ``AFIX`` group so a
    group is never split.

    :param sof: The raw SHELXL site-occupation factor written on every atom,
        e.g. ``21.0`` for "free variable 2" or ``-21.0`` for "1 - fv 2" --
        either one value for all atoms or one per atom, parallel to *atoms*
        (atoms on special positions keep their reduced occupancy that way,
        e.g. ``20.5``).
    :returns: The ``PART`` cards that open each run.
    :raises ValueError: for atoms already in a disorder part, symmetry
        images, Q-peaks, ``part == 0`` or a partly selected ``AFIX`` group.
    """
    from shelxfile.atoms.atom import Atom as AtomClass
    from shelxfile.shelx.cards import PART

    if part == 0:
        raise ValueError('Use a non-zero PART number for a disorder part')
    given = list(atoms)
    sofs = _per_atom(sof, len(given), 'sof')
    ordered = _file_order(shx, given)
    if not ordered:
        return []
    _check_real_atoms(ordered)
    for atom in ordered:
        if atom.part.n != 0:
            raise ValueError(f'{atom.name} is already in PART {atom.part.n}')
    problem = afix_closure_problem(shx, ordered)
    if problem:
        raise ValueError(f'Cannot assign a PART: {problem}')

    reslist = shx._reslist
    selected = {id(atom) for atom in ordered}
    runs: list[list[Atom]] = [[ordered[0]]]
    for atom in ordered[1:]:
        previous = runs[-1][-1]
        between = reslist[shx.index_of(previous) + 1:shx.index_of(atom)]
        joined = all(
            _is_skippable(item) or _is_afix(item)
            or (isinstance(item, AtomClass) and id(item) in selected)
            for item in between
        )
        if joined:
            runs[-1].append(atom)
        else:
            runs.append([atom])

    openers = []
    # Work backwards so earlier insertions do not move later positions.
    for run in reversed(runs):
        close_at = shx.index_of(run[-1]) + 1
        following = _next_meaningful(reslist, close_at)
        if following < len(reslist) and _is_afix(reslist[following]) \
                and not reslist[following].mn:
            close_at = following + 1  # after the AFIX 0 that ends the run's group
        open_at = shx.index_of(run[0])
        preceding = _previous_meaningful(reslist, open_at - 1)
        if preceding >= 0 and _is_afix(reslist[preceding]) and reslist[preceding].mn:
            open_at = preceding  # the run starts a rigid group: open before its AFIX

        opener = PART(shx, ['PART', str(part)])
        shx.insert_into_reslist(close_at, PART(shx, ['PART', '0']))
        shx.insert_into_reslist(open_at, opener)
        for atom in run:
            atom.part = opener
        openers.append(opener)

    for atom, value in zip(given, sofs):
        atom.sof = value
    shx.atoms._atomsdict.clear()
    shx.touch()
    return list(reversed(openers))


# ------------------------------------------------------------ duplication

def _copy_insert_position(shx: Shelxfile, last_source: Atom) -> int:
    """Where a duplicated block goes: right after the sources' own part.

    Inserting inside the sources' ``PART`` bracket would end that bracket
    early (the block closes with ``PART 0``), so the block goes after the
    bracket's last member instead.
    """
    reslist = shx._reslist
    anchor = last_source
    if last_source.part.n != 0:
        scope = BracketResolver(shx).scope_of(last_source.part)
        if scope is not None and scope.members:
            anchor = scope.members[-1]
    position = shx.index_of(anchor) + 1
    following = _next_meaningful(reslist, position)
    if following < len(reslist) and _is_afix(reslist[following]) and not reslist[following].mn:
        position = following + 1
    if anchor.part.n != 0:
        following = _next_meaningful(reslist, position)
        if following < len(reslist) and _is_part(reslist[following]) and reslist[following].n == 0:
            position = following + 1
    return position


def duplicate_atoms(
    shx: Shelxfile,
    atoms: Sequence[Atom],
    names: Sequence[str],
    part: int,
    sof: float | Sequence[float],
    frac_coordinates: Sequence[Sequence[float] | None] | None = None,
    uvals: Sequence[Sequence[float] | None] | None = None,
) -> list[Atom]:
    """Append copies of *atoms* as a new ``PART part`` block.

    The copies are written in the sources' file order, in one
    ``PART part`` ... ``PART 0`` block placed after the sources, each
    ``AFIX`` card copied in front of the atoms it governed and closed with
    ``AFIX 0`` -- so riding hydrogens stay riding and riding U codes such as
    ``-1.2`` still refer to the (copied) atom before them.

    ``Atom`` objects are not hashable, so the per-atom arguments are
    sequences parallel to *atoms*.

    :param names: New name per source atom.  Must be valid and unused in the
        source's residue.
    :param sof: Raw site-occupation factor, one for every copy or one per atom.
    :param frac_coordinates: Optional new fractional coordinates per source;
        ``None`` (or no sequence at all) copies the atom in place.
    :param uvals: Optional U values per source (one value for isotropic or
        six); ``None`` copies the source's raw values, codes included.
    :returns: The copies, parallel to *atoms*.
    """
    from shelxfile.atoms.atom import Atom as AtomClass
    from shelxfile.shelx.cards import AFIX, PART

    atoms = list(atoms)
    if not atoms:
        return []
    if len(names) != len(atoms):
        raise ValueError('names must be parallel to atoms')
    if frac_coordinates is not None and len(frac_coordinates) != len(atoms):
        raise ValueError('frac_coordinates must be parallel to atoms')
    if uvals is not None and len(uvals) != len(atoms):
        raise ValueError('uvals must be parallel to atoms')
    if len({id(atom) for atom in atoms}) != len(atoms):
        raise ValueError('an atom is listed twice')
    _check_real_atoms(atoms)
    problem = afix_closure_problem(shx, atoms)
    if problem:
        raise ValueError(f'Cannot duplicate: {problem}')
    reserved: set[str] = set()
    for atom, name in zip(atoms, names):
        if not _name_is_free(shx, name, atom.resinum, reserved):
            raise ValueError(f'Atom name {name!r} is invalid or already used')
        reserved.add(name.upper())

    position_of = {id(atom): i for i, atom in enumerate(atoms)}
    sofs = _per_atom(sof, len(atoms), 'sof')
    frac_list = list(frac_coordinates) if frac_coordinates is not None else [None] * len(atoms)
    u_list = list(uvals) if uvals is not None else [None] * len(atoms)
    ordered = _file_order(shx, atoms)
    reslist = shx._reslist
    part_card = PART(shx, ['PART', str(part)])
    block: list = [part_card]
    copy_of: dict[int, Atom] = {}
    current_afix = None
    for source in ordered:
        i = position_of[id(source)]
        index = shx.index_of(source)
        preceding = _previous_meaningful(reslist, index - 1)
        if preceding >= 0 and _is_afix(reslist[preceding]) and reslist[preceding].mn:
            current_afix = AFIX(shx, list(reslist[preceding]._spline))
            block.append(current_afix)
        elif not (source.afix and source.afix.mn):
            current_afix = None

        coordinates = frac_list[i]
        if coordinates is None:
            coordinates = (source.x, source.y, source.z)
        u_values = list(u_list[i] if u_list[i] is not None else source.uvals)
        if len(u_values) == 1:
            u_values += [0.0] * 5
        copy = AtomClass(shx)
        copy.set_atom_parameters(
            name=names[i], sfac_num=source.sfac_num,
            coords=[float(value) for value in coordinates],
            part=part_card, afix=current_afix, resi=source.resi,
            site_occupation=sofs[i], uvals=u_values, symmgen=False,
        )
        copy.uvals_orig = list(u_values)
        copy_of[id(source)] = copy
        block.append(copy)

        following = _next_meaningful(reslist, index + 1)
        if current_afix is not None and following < len(reslist) \
                and _is_afix(reslist[following]) and not reslist[following].mn:
            block.append(AFIX(shx, ['AFIX', '0']))
            current_afix = None
    if current_afix is not None:
        block.append(AFIX(shx, ['AFIX', '0']))
    block.append(PART(shx, ['PART', '0']))

    # Riding U values and pivots must point at the copies, not the sources.
    for source in ordered:
        copy = copy_of[id(source)]
        for attribute in ('u_reference', 'pivot'):
            reference = getattr(source, attribute)
            if reference is not None:
                reference = copy_of.get(id(reference), reference)
            setattr(copy, attribute, reference)

    position = _copy_insert_position(shx, ordered[-1])
    for offset, item in enumerate(block):
        shx.insert_into_reslist(position + offset, item)
    for source in ordered:
        shx.atoms.append(copy_of[id(source)])
    return [copy_of[id(atom)] for atom in atoms]
