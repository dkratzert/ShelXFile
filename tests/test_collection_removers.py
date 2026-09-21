"""Symmetric removal helpers for the card collections (plan item C4).

``Restraints`` and ``Residues`` both had ``append()`` but no ``remove()``,
so deletion code reached into ``_restraints`` directly and no one undid
the bookkeeping that ``append()`` performs.

``Residues.append()`` in particular maintains a ``class -> [numbers]``
map that allows duplicates, so removal has to drop exactly one entry and
retire a class once its last residue is gone.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.shelx.cards import RESI, Restraint

RES = 'tests/resources/p21c.res'

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C\n'
    'UNIT 4\n'
)
FOOTER = 'HKLF 4\nEND\n'


@pytest.fixture
def shx() -> Shelxfile:
    s = Shelxfile()
    s.read_file(RES)
    return s


def _restraint_with_atoms(shx: Shelxfile) -> Restraint:
    for restraint in shx.restraints:
        if [a for a in restraint.atoms if a not in ('>', '<', '=')]:
            return restraint
    raise AssertionError('fixture has no atom-bearing restraint')


# ------------------------------------------------------------- Restraints

def test_remove_drops_the_restraint_everywhere(shx: Shelxfile) -> None:
    restraint = _restraint_with_atoms(shx)
    shx.restraints.remove(restraint)
    assert all(item is not restraint for item in shx.restraints)
    assert not any(item is restraint for item in shx._reslist)


def test_remove_keeps_the_referenced_atoms(shx: Shelxfile) -> None:
    """D-9b again, at the collection level."""
    restraint = _restraint_with_atoms(shx)
    before = len(shx.atoms)
    shx.restraints.remove(restraint)
    assert len(shx.atoms) == before


def test_remove_is_idempotent(shx: Shelxfile) -> None:
    restraint = _restraint_with_atoms(shx)
    shx.restraints.remove(restraint)
    shx.restraints.remove(restraint)  # must not raise


def test_restraint_delete_uses_the_collection(shx: Shelxfile) -> None:
    """One implementation, so the two can't drift apart."""
    restraint = _restraint_with_atoms(shx)
    before = len(shx.restraints)
    restraint.delete()
    assert len(shx.restraints) == before - 1


def test_removing_one_of_several_leaves_the_rest(shx: Shelxfile) -> None:
    victims = [r for r in shx.restraints if type(r).__name__ == 'SADI']
    assert len(victims) > 2, 'fixture should have several SADI cards'
    survivors = [r for r in shx.restraints if r is not victims[1]]
    shx.restraints.remove(victims[1])
    assert [r for r in shx.restraints] == survivors


# --------------------------------------------------------------- Residues

def _with_residues() -> Shelxfile:
    shx = Shelxfile()
    shx.read_string(
        HEADER
        + 'RESI 1 TOL\n'
          'C1    1     0.1  0.2  0.3  11.00000  0.05\n'
          'RESI 2 TOL\n'
          'C2    1     0.2  0.3  0.4  11.00000  0.05\n'
          'RESI 3 BEN\n'
          'C3    1     0.3  0.4  0.5  11.00000  0.05\n'
          'RESI 0\n'
        + FOOTER
    )
    return shx


def test_remove_drops_the_residue_everywhere() -> None:
    shx = _with_residues()
    resi = next(r for r in shx.residues.all_residues if r.residue_number == 2)
    shx.residues.remove(resi)
    assert all(r is not resi for r in shx.residues.all_residues)
    assert not any(item is resi for item in shx._reslist)


def test_remove_drops_exactly_one_class_entry() -> None:
    """``TOL`` holds two numbers; removing residue 2 must leave residue 1."""
    shx = _with_residues()
    assert shx.residues.residue_classes['TOL'] == [1, 2]

    resi = next(r for r in shx.residues.all_residues if r.residue_number == 2)
    shx.residues.remove(resi)

    assert shx.residues.residue_classes['TOL'] == [1]


def test_class_is_retired_when_its_last_residue_goes() -> None:
    shx = _with_residues()
    resi = next(r for r in shx.residues.all_residues if r.residue_class == 'BEN')
    shx.residues.remove(resi)
    assert 'BEN' not in shx.residues.residue_classes


def test_duplicate_numbers_lose_only_one_entry() -> None:
    """A class split across several RESI cards may repeat a number."""
    shx = _with_residues()
    extra = RESI(shx, ['RESI', '1', 'TOL'])
    shx.residues.append(extra)
    assert shx.residues.residue_classes['TOL'] == [1, 2, 1]

    shx.residues.remove(extra)

    assert shx.residues.residue_classes['TOL'] == [2, 1], (
        'removal dropped every matching number instead of one'
    )


def test_remove_is_idempotent_for_residues() -> None:
    shx = _with_residues()
    resi = shx.residues.all_residues[0]
    shx.residues.remove(resi)
    shx.residues.remove(resi)  # must not raise


def test_residue_numbers_reflect_removal() -> None:
    shx = _with_residues()
    resi = next(r for r in shx.residues.all_residues if r.residue_number == 3)
    shx.residues.remove(resi)
    assert 3 not in shx.residues.residue_numbers
