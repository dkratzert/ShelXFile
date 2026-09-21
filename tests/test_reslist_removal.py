"""Tests for identity-based ``_reslist`` removal (plan items A2, A3, A4, D-6e).

These cover the prerequisites that everything else in the
atom/restraint-linking work sits on:

* **A2** ``delete_on_write`` holds *indices* and must be re-mapped when a
  line disappears, or the wrong line gets suppressed on the next write.
* **A3** removal must match by identity, not by ``__eq__``; ``Atom.__eq__``
  compares the rendered line, so look-alike atoms would collide.
* **A4** ``Atoms.__delitem__`` must not mutate the list it enumerates.
* **D-6e** symmetry-generated atoms are not lines in the file and must
  never be deletion targets.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.atoms.atom import Atom
from shelxfile.shelx.cards import Restraint

RES = 'tests/resources/p21c.res'


@pytest.fixture
def shx() -> Shelxfile:
    s = Shelxfile()
    s.read_file(RES)
    return s


# --------------------------------------------------------------- A3: identity

def test_index_of_matches_by_identity_not_equality(shx: Shelxfile) -> None:
    """Two equal-rendering atoms must resolve to their own positions."""
    first = shx.atoms.all_atoms[0]
    clone = Atom(shx)
    clone.__dict__.update(first.__dict__)
    assert clone == first, 'precondition: the clone renders identically'

    position = shx.index_of(first)
    shx._reslist.insert(position + 1, clone)

    assert shx.index_of(first) == position
    assert shx.index_of(clone) == position + 1, (
        'index_of() fell back to __eq__ and returned the first look-alike'
    )


def test_index_of_raises_for_absent_object(shx: Shelxfile) -> None:
    with pytest.raises(ValueError):
        shx.index_of(Atom(shx))


def test_deleting_one_of_two_lookalike_atoms_keeps_the_other(shx: Shelxfile) -> None:
    """A3 - the surviving twin must still be present and addressable."""
    original = shx.atoms.all_atoms[0]
    clone = Atom(shx)
    clone.__dict__.update(original.__dict__)
    position = shx.index_of(original)
    shx._reslist.insert(position + 1, clone)
    shx.atoms.all_atoms.insert(1, clone)
    shx.atoms._atomsdict.clear()

    clone.delete()

    assert any(item is original for item in shx._reslist), 'deleted the wrong twin'
    assert not any(item is clone for item in shx._reslist)


# ------------------------------------------------- A2: delete_on_write shifts

def _count_lines_starting_with(dump: str, name: str) -> int:
    """How many rendered lines begin with *name*.

    ``dumps()`` wraps long lines with ``=`` continuations, so comparing
    whole atom lines is unreliable.  Names are also not unique -- the p21c
    fixture holds five atoms called ``F2`` in different residues -- so the
    meaningful assertion is how the *count* changes, not presence.
    """
    total = 0
    for line in dump.splitlines():
        if line.startswith(' '):
            continue  # continuation line
        parts = line.split()
        if parts and parts[0].upper() == name.upper():
            total += 1
    return total


def test_delete_on_write_indices_follow_removed_lines(shx: Shelxfile) -> None:
    """A2 - suppressing line N must still suppress that same *line* afterwards."""
    atoms = shx.atoms.all_atoms
    victim, marked = atoms[0], atoms[4]
    baseline = _count_lines_starting_with(shx.dumps(), marked.name)

    shx.delete_on_write.add(shx.index_of(marked))
    suppressed = _count_lines_starting_with(shx.dumps(), marked.name)
    assert suppressed == baseline - 1, 'precondition: the mark suppresses one line'

    victim.delete()

    assert _count_lines_starting_with(shx.dumps(), marked.name) == suppressed, (
        'delete_on_write was not re-mapped: the marked line reappeared'
    )


def test_delete_on_write_does_not_suppress_an_unrelated_line(shx: Shelxfile) -> None:
    """The inverse of the above: a *different* line must not get suppressed."""
    atoms = shx.atoms.all_atoms
    victim, marked, neighbour = atoms[0], atoms[4], atoms[5]
    neighbour_before = _count_lines_starting_with(shx.dumps(), neighbour.name)
    shx.delete_on_write.add(shx.index_of(marked))

    victim.delete()

    assert _count_lines_starting_with(shx.dumps(), neighbour.name) == neighbour_before, (
        'the shift over-corrected and suppressed the following line'
    )


def test_delete_on_write_shift_is_exact(shx: Shelxfile) -> None:
    """Pure index arithmetic, independent of rendering."""
    atoms = shx.atoms.all_atoms
    victim = atoms[3]
    position = shx.index_of(victim)
    shx.delete_on_write.update({position - 1, position, position + 1, position + 5})

    victim.delete()

    assert shx.delete_on_write == {position - 1, position, position + 4}, (
        'indices below the removal must stay, the removed one must go, '
        'and those above must move down by exactly one'
    )


def test_delete_on_write_drops_the_index_of_the_removed_line(shx: Shelxfile) -> None:
    victim = shx.atoms.all_atoms[0]
    shx.delete_on_write.add(shx.index_of(victim))

    victim.delete()

    assert shx.delete_on_write == set()


def test_remove_from_reslist_returns_the_freed_index(shx: Shelxfile) -> None:
    atom = shx.atoms.all_atoms[2]
    expected = shx.index_of(atom)
    assert shx.remove_from_reslist(atom) == expected
    assert not any(item is atom for item in shx._reslist)


# ---------------------------------------------- A4: safe deletion while paging

def test_deleting_every_atom_removes_them_all(shx: Shelxfile) -> None:
    """A4 - enumerating while deleting used to skip every other entry."""
    total = len(shx.atoms.all_atoms)
    assert total > 5, 'fixture should have a decent number of atoms'

    for atom in list(shx.atoms.all_atoms):
        atom.delete()

    assert shx.atoms.all_atoms == []
    assert not any(isinstance(item, Atom) for item in shx._reslist)


def test_deleting_consecutive_atoms_leaves_the_rest_intact(shx: Shelxfile) -> None:
    atoms = shx.atoms.all_atoms
    doomed = [atoms[1], atoms[2], atoms[3]]
    survivors = [a.fullname for a in atoms if a not in doomed]

    for atom in doomed:
        atom.delete()

    assert [a.fullname for a in shx.atoms.all_atoms] == survivors


# ------------------------------------------------------- D-6e: symmgen guards

def test_symmetry_generated_atoms_are_not_deletable(shx: Shelxfile) -> None:
    """D-6e - grown atoms are not lines in the file."""
    atom = shx.atoms.all_atoms[0]
    atom.symmgen = True
    before = len(shx._reslist)

    atom.delete()

    assert len(shx._reslist) == before
    assert any(item is atom for item in shx._reslist)


def test_packed_atoms_are_never_removed_by_deletion(shx: Shelxfile) -> None:
    """D-6e - symmetry mates from pack() are not lines in the file."""
    packed = shx.pack()
    symmgen = [a for a in packed if a.symmgen]
    assert symmgen, 'precondition: pack() produced symmetry mates'
    before = len(shx._reslist)

    for atom in symmgen[:5]:
        atom.delete()

    assert len(shx._reslist) == before


# --------------------------------------------- D-9: restraint removal is safe

def _restraint_with_atoms(shx: Shelxfile) -> Restraint:
    """First restraint that actually names atoms.

    The p21c fixture opens with a bare ``DELU``, which under SHELXL rules
    means "all non-hydrogen atoms" (plan item B2) and names nothing.
    """
    for restraint in shx.restraints:
        if [a for a in restraint.atoms if a not in ('>', '<', '=')]:
            return restraint
    raise AssertionError('fixture has no atom-bearing restraint')


def test_deleting_a_restraint_keeps_its_atoms(shx: Shelxfile) -> None:
    """D-9b - card removal must never reach an atom."""
    restraint = _restraint_with_atoms(shx)
    named = [a for a in restraint.atoms if a not in ('>', '<', '=')]
    before = len(shx.atoms.all_atoms)

    restraint.delete()

    assert len(shx.atoms.all_atoms) == before
    assert not any(item is restraint for item in shx._reslist)
    for name in named:
        assert shx.atoms.get_atom_by_name(f'{name}_0') is not None, (
            f'atom {name} was removed along with its restraint'
        )


def test_deleting_a_restraint_removes_it_from_the_collection(shx: Shelxfile) -> None:
    restraint = _restraint_with_atoms(shx)
    restraint.delete()
    assert all(item is not restraint for item in shx.restraints)
    assert not any(isinstance(i, Restraint) and i is restraint for i in shx._reslist)


def test_bare_restraint_with_no_atoms_is_still_deletable(shx: Shelxfile) -> None:
    """B2 - a card authored without atoms is valid and removable as a card."""
    bare = next(r for r in shx.restraints if not r.atoms)
    before = len(shx.atoms.all_atoms)

    bare.delete()

    assert len(shx.atoms.all_atoms) == before
    assert not any(item is bare for item in shx._reslist)
