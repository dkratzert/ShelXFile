"""Bracket scopes for ``RESI``, ``AFIX`` and ``PART`` (plan items B5, B6, D-4).

The three cards all bracket a run of atoms but end differently, and
``AFIX`` additionally requires a *specific number* of members -- a
constraint the original edit plan did not model at all.

Includes the synthetic fixtures for constructs the reference corpus does
not cover; frequency says nothing about whether the code has to be right.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit.brackets import (
    AFIX_EXPECTED_ATOMS,
    BracketResolver,
    expected_atom_count,
)
from shelxfile.edit.card_meta import AfixDependency

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C H\n'
    'UNIT 8 8\n'
)
FOOTER = 'HKLF 4\nEND\n'


def _resolver(body: str) -> tuple[Shelxfile, BracketResolver]:
    shx = Shelxfile()
    shx.read_string(HEADER + body + FOOTER)
    return shx, BracketResolver(shx)


def _atom(name: str, num: int = 1, sfac: int = 1) -> str:
    return f'{name:<5} {sfac}     0.{num:02d}0  0.{num:02d}5  0.{num:02d}9  11.00000  0.05\n'


# ------------------------------------------------- B5: expected atom counts

@pytest.mark.parametrize('mn, expected', [
    (13, 1),    # tertiary C-H
    (23, 2),    # secondary CH2
    (33, 3),    # CH3
    (43, 1),    # aromatic C-H
    (56, 5),    # pentagon
    (66, 6),    # hexagon
    (83, 1),    # OH
    (93, 2),    # terminal X=CH2
    (105, 10),  # Cp*
    (116, 8),   # naphthalene
    (123, 6),   # disordered methyl
    (137, 3),   # CH3 with torsion from the map
    (163, 1),   # acetylenic C-H
])
def test_geometry_code_sets_a_member_count(mn, expected) -> None:
    assert expected_atom_count(mn) == expected


def test_only_the_leading_digits_pick_the_geometry() -> None:
    """``m`` selects the geometry, ``n`` the constraint."""
    assert expected_atom_count(66) == expected_atom_count(65)


@pytest.mark.parametrize('mn', [0, 3, 9, None])
def test_codes_without_a_geometry_have_no_requirement(mn) -> None:
    assert expected_atom_count(mn) is None


def test_frag_fitted_groups_have_no_fixed_count() -> None:
    """``m`` > 16 fits a FRAG block, whose size only that block knows."""
    assert expected_atom_count(177) is None


def test_every_documented_geometry_is_covered() -> None:
    assert set(AFIX_EXPECTED_ATOMS) == set(range(1, 17))


# ---------------------------------------------------------- B6: scope ends

def test_afix_scope_ends_at_the_next_afix() -> None:
    shx, resolver = _resolver(
        _atom('C1') + 'AFIX 43\n' + _atom('H1', 2, 2) + 'AFIX 0\n' + _atom('C2', 3)
    )
    scopes = resolver.afix_scopes()
    assert [a.name for a in scopes[0].members] == ['H1']


def test_part_scope_ends_at_the_next_part() -> None:
    shx, resolver = _resolver(
        'PART 1\n' + _atom('C1') + 'PART 2\n' + _atom('C2', 2) + 'PART 0\n'
    )
    scopes = resolver.part_scopes()
    assert [a.name for a in scopes[0].members] == ['C1']
    assert [a.name for a in scopes[1].members] == ['C2']


def test_resi_scope_ends_at_the_next_resi() -> None:
    shx, resolver = _resolver(
        'RESI 1 TOL\n' + _atom('C1') + 'RESI 2 TOL\n' + _atom('C2', 2) + 'RESI 0\n'
    )
    scopes = resolver.resi_scopes()
    assert [a.name for a in scopes[0].members] == ['C1']


def test_rotating_group_runs_through_further_n7_cards() -> None:
    """"until the next AFIX with n not equal to 7"."""
    shx, resolver = _resolver(
        _atom('C1')
        + 'AFIX 137\n' + _atom('H1', 2, 2)
        + 'AFIX 137\n' + _atom('H2', 3, 2)
        + 'AFIX 0\n'
    )
    scopes = resolver.afix_scopes()
    rotating = [s for s in scopes if getattr(s.card, 'mn', None) == 137]
    assert len(rotating) == 1, 'a second n=7 card must not open a new scope'
    assert [a.name for a in rotating[0].members] == ['H1', 'H2']


def test_a_non_seven_afix_does_close_a_rotating_group() -> None:
    shx, resolver = _resolver(
        _atom('C1')
        + 'AFIX 137\n' + _atom('H1', 2, 2)
        + 'AFIX 43\n' + _atom('H2', 3, 2)
        + 'AFIX 0\n'
    )
    rotating = [s for s in resolver.afix_scopes()
                if getattr(s.card, 'mn', None) == 137]
    assert [a.name for a in rotating[0].members] == ['H1']


# ------------------------------------------------------- resets are sacred

def test_afix_zero_is_a_reset() -> None:
    shx, resolver = _resolver(_atom('C1') + 'AFIX 43\n' + _atom('H1', 2, 2) + 'AFIX 0\n')
    closers = [s for s in resolver.afix_scopes() if s.is_reset]
    assert closers, 'AFIX 0 must be recognised as a reset'


def test_part_zero_is_a_reset_but_negative_parts_are_not() -> None:
    """``PART -n`` suppresses special-position constraints; it is a real group."""
    shx, resolver = _resolver(
        'PART -1\n' + _atom('C1') + 'PART 0\n' + _atom('C2', 2)
    )
    scopes = resolver.part_scopes()
    assert scopes[0].is_reset is False
    assert scopes[1].is_reset is True


def test_resi_zero_is_a_reset() -> None:
    shx, resolver = _resolver('RESI 1 TOL\n' + _atom('C1') + 'RESI 0\n')
    scopes = resolver.resi_scopes()
    assert scopes[0].is_reset is False
    assert scopes[1].is_reset is True


def test_a_reset_is_never_under_populated() -> None:
    shx, resolver = _resolver(_atom('C1') + 'AFIX 0\n')
    assert all(not s.is_under_populated for s in resolver.afix_scopes() if s.is_reset)


# ------------------------------------------------- B5: under-population

def test_full_group_is_not_under_populated() -> None:
    body = _atom('C1') + 'AFIX 66\n' + ''.join(
        _atom(f'C{i}', i) for i in range(2, 8)
    ) + 'AFIX 0\n'
    shx, resolver = _resolver(body)
    hexagon = next(s for s in resolver.afix_scopes()
                   if getattr(s.card, 'mn', None) == 66)
    assert len(hexagon.members) == 6
    assert hexagon.is_under_populated is False


def test_losing_one_atom_under_populates_a_hexagon() -> None:
    """Five atoms cannot be fitted to a six-atom ring."""
    body = _atom('C1') + 'AFIX 66\n' + ''.join(
        _atom(f'C{i}', i) for i in range(2, 7)
    ) + 'AFIX 0\n'
    shx, resolver = _resolver(body)
    hexagon = next(s for s in resolver.afix_scopes()
                   if getattr(s.card, 'mn', None) == 66)
    assert len(hexagon.members) == 5
    assert hexagon.is_under_populated is True


def test_empty_group_is_reported() -> None:
    shx, resolver = _resolver(_atom('C1') + 'AFIX 43\n' + 'AFIX 0\n')
    riding = next(s for s in resolver.afix_scopes()
                  if getattr(s.card, 'mn', None) == 43)
    assert riding.is_empty
    assert riding.is_under_populated


# ------------------------------------------------------ D-4: outside pivots

def test_riding_group_pivots_on_the_preceding_atom() -> None:
    shx, resolver = _resolver(
        _atom('C1') + 'AFIX 43\n' + _atom('H1', 2, 2) + 'AFIX 0\n'
    )
    riding = next(s for s in resolver.afix_scopes()
                  if getattr(s.card, 'mn', None) == 43)
    assert riding.dependency is AfixDependency.RIDING
    assert resolver.pivot_of(riding).name == 'C1'


def test_rigid_group_has_no_outside_pivot() -> None:
    """Its pivot is the first atom *inside* the bracket."""
    shx, resolver = _resolver(
        _atom('C1') + 'AFIX 66\n' + ''.join(_atom(f'C{i}', i) for i in range(2, 8))
        + 'AFIX 0\n'
    )
    hexagon = next(s for s in resolver.afix_scopes()
                   if getattr(s.card, 'mn', None) == 66)
    assert resolver.pivot_of(hexagon) is None


def test_pivot_search_skips_other_riders() -> None:
    """"The atom on which riding is performed may not itself be a riding atom"."""
    shx, resolver = _resolver(
        _atom('C1')
        + 'AFIX 43\n' + _atom('H1', 2, 2) + 'AFIX 0\n'
        + 'AFIX 43\n' + _atom('H2', 3, 2) + 'AFIX 0\n'
    )
    groups = [s for s in resolver.afix_scopes()
              if getattr(s.card, 'mn', None) == 43]
    assert resolver.pivot_of(groups[1]).name == 'C1', (
        'the second group must not ride on H1'
    )


def test_scopes_pivoted_on_finds_the_orphan_candidates() -> None:
    shx, resolver = _resolver(
        _atom('C1') + 'AFIX 43\n' + _atom('H1', 2, 2) + 'AFIX 0\n'
    )
    c1 = shx.atoms.get_atom_by_name('C1_0')
    assert len(resolver.scopes_pivoted_on(c1)) == 1


def test_an_unrelated_atom_orphans_nothing() -> None:
    shx, resolver = _resolver(
        _atom('C1') + 'AFIX 43\n' + _atom('H1', 2, 2) + 'AFIX 0\n' + _atom('C9', 9)
    )
    c9 = shx.atoms.get_atom_by_name('C9_0')
    assert resolver.scopes_pivoted_on(c9) == []


def test_dependent_atoms_do_not_carry_their_own_count() -> None:
    """``n`` = 5 continues the enclosing rigid group.

    ``AFIX 65`` marks dependent atoms of a hexagon already opened by an
    ``AFIX 66``, so the six-atom requirement belongs to that group, not to
    this continuation.  Treating it as a new group reported 61 false
    under-populations across the reference corpus.
    """
    shx, resolver = _resolver(
        _atom('C1') + 'AFIX 66\n' + ''.join(_atom(f'C{i}', i) for i in range(2, 8))
        + 'AFIX 65\n' + _atom('C8', 8) + 'AFIX 0\n'
    )
    dependent = next(s for s in resolver.afix_scopes()
                     if getattr(s.card, 'mn', None) == 65)
    assert dependent.dependency is AfixDependency.RIGID_DEPENDENT
    assert dependent.expected_count is None
    assert dependent.is_under_populated is False


def test_the_opening_rigid_group_still_carries_the_count() -> None:
    shx, resolver = _resolver(
        _atom('C1') + 'AFIX 66\n' + _atom('C2', 2) + 'AFIX 65\n' + _atom('C3', 3)
        + 'AFIX 0\n'
    )
    hexagon = next(s for s in resolver.afix_scopes()
                   if getattr(s.card, 'mn', None) == 66)
    assert hexagon.expected_count == 6
    assert hexagon.is_under_populated is True


# --------------------------------------------------------------- real file

def test_real_file_scopes_resolve() -> None:
    shx = Shelxfile()
    shx.read_file('tests/resources/p21c.res')
    resolver = BracketResolver(shx)
    assert resolver.afix_scopes() or resolver.part_scopes() or resolver.resi_scopes()
    for scope in resolver.resi_scopes():
        assert scope.card is not None
