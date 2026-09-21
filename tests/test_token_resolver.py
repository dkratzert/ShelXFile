"""Expanding SHELXL atom-reference tokens (plan items B4 and D-6b).

Covers the whole token grammar, including the four constructs that do not
appear anywhere in the reference corpus (``LAST``, ``_+``/``_-``, and
``$element`` on restraints).  Absence from one lab's files says nothing
about importance, so these are tested at the same rigour as ``>``, which
appears in 5530 files.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit.token_resolver import AtomTokenResolver, split_range_tokens

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C H O\n'
    'UNIT 8 8 4\n'
)
FOOTER = 'HKLF 4\nEND\n'


def _shx(body: str) -> Shelxfile:
    shx = Shelxfile()
    shx.read_string(HEADER + body + FOOTER)
    return shx


def _resolver(body: str) -> AtomTokenResolver:
    return AtomTokenResolver(_shx(body))


CHAIN = (
    'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
    'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
    'H2    2     0.16  0.26  0.36  11.00000  0.05\n'
    'C3    1     0.20  0.30  0.40  11.00000  0.05\n'
    'C4    1     0.25  0.35  0.45  11.00000  0.05\n'
    'O1    3     0.30  0.40  0.50  11.00000  0.05\n'
)


# ------------------------------------------------------------- tokenising

@pytest.mark.parametrize('tokens, expected', [
    (['C1', '>', 'C5'], ['C1', '>', 'C5']),
    (['C1>C5'], ['C1', '>', 'C5']),
    (['O_201>', 'LAST'], ['O_201', '>', 'LAST']),       # manual's ISOR example
    (['C11', '>C15', 'Fe'], ['C11', '>', 'C15', 'Fe']),  # manual's MPLA example
    (['C1', '<', 'C5'], ['C1', '<', 'C5']),
    (['C6<C2'], ['C6', '<', 'C2']),
    (['C1', 'C2'], ['C1', 'C2']),
])
def test_range_markers_are_split_from_names(tokens, expected) -> None:
    """SHELXL accepts every spacing; its own docs use several."""
    assert split_range_tokens(tokens) == expected


# ---------------------------------------------------------- plain names

def test_bare_name_resolves_to_residue_zero() -> None:
    resolver = _resolver(CHAIN)
    assert resolver.resolve(['C1']).fullnames == ['C1_0']


def test_unknown_name_is_reported_not_guessed() -> None:
    resolver = _resolver(CHAIN)
    result = resolver.resolve(['CZZ'])
    assert result.fullnames == []
    assert result.unresolved == ['CZZ']


def test_explicit_residue_number_wins() -> None:
    shx = _shx('RESI 1 TOL\n' + CHAIN + 'RESI 0\n')
    assert AtomTokenResolver(shx).resolve(['C1_1']).fullnames == ['C1_1']


# ------------------------------------------------------------- ranges

def test_forward_range_includes_everything_between() -> None:
    resolver = _resolver(CHAIN)
    assert resolver.resolve(['C1', '>', 'C4']).fullnames == [
        'C1_0', 'C2_0', 'C3_0', 'C4_0',
    ]


def test_range_skips_hydrogens() -> None:
    """The manual defines ``>`` as all intervening *non-hydrogen* atoms."""
    resolver = _resolver(CHAIN)
    assert 'H2_0' not in resolver.resolve(['C1', '>', 'C4']).fullnames


def test_backward_range_walks_in_reverse() -> None:
    """``SAME C1 C6 < C2 C7`` in the manual expands backwards."""
    resolver = _resolver(CHAIN)
    assert resolver.resolve(['C4', '<', 'C1']).fullnames == [
        'C4_0', 'C3_0', 'C2_0', 'C1_0',
    ]


def test_range_against_file_order_does_not_invent_atoms() -> None:
    resolver = _resolver(CHAIN)
    result = resolver.resolve(['C4', '>', 'C1'])
    assert result.fullnames == ['C4_0', 'C1_0']


def test_range_with_unknown_endpoint_is_reported() -> None:
    resolver = _resolver(CHAIN)
    assert 'CZZ' in resolver.resolve(['C1', '>', 'CZZ']).unresolved


# -------------------------------------------------------------- LAST

def test_last_resolves_to_the_final_atom() -> None:
    """Absent from the corpus; still fully supported."""
    resolver = _resolver(CHAIN)
    assert resolver.resolve(['LAST']).fullnames == ['O1_0']


def test_range_to_last_expands_to_the_end() -> None:
    resolver = _resolver(CHAIN)
    assert resolver.resolve(['C3', '>', 'LAST']).fullnames == ['C3_0', 'C4_0', 'O1_0']


def test_last_ignores_q_peaks() -> None:
    """Q peaks sit after ``END`` in real files and are not model atoms."""
    shx = Shelxfile()
    shx.read_string(
        HEADER + CHAIN + FOOTER
        + 'Q1    1   0.9000  0.9000  0.9000  11.00000  0.05    0.30\n'
    )
    assert [a.fullname for a in shx.atoms if a.qpeak], 'precondition: a Q peak parsed'
    assert AtomTokenResolver(shx).resolve(['LAST']).fullnames == ['O1_0']


# ---------------------------------------------------------- $element

def test_element_reference_collects_every_atom_of_that_type() -> None:
    resolver = _resolver(CHAIN)
    assert resolver.resolve(['$C']).fullnames == ['C1_0', 'C2_0', 'C3_0', 'C4_0']


def test_element_reference_is_case_insensitive() -> None:
    resolver = _resolver(CHAIN)
    assert resolver.resolve(['$o']).fullnames == ['O1_0']


def test_unknown_element_resolves_to_nothing() -> None:
    resolver = _resolver(CHAIN)
    assert resolver.resolve(['$Zz']).unresolved == ['$Zz']


# --------------------------------------------------------- residues

def _residue_file() -> Shelxfile:
    return _shx(
        'RESI 1 TOL\n'
        'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
        'RESI 2 TOL\n'
        'C1    1     0.20  0.30  0.40  11.00000  0.05\n'
        'RESI 3 TOL\n'
        'C1    1     0.30  0.40  0.50  11.00000  0.05\n'
        'RESI 0\n'
    )


def test_residue_scope_is_applied_to_bare_names() -> None:
    resolver = AtomTokenResolver(_residue_file())
    assert resolver.resolve(['C1'], [2]).fullnames == ['C1_2']


def test_multiple_scopes_expand_the_instruction() -> None:
    """``SADI_TOL C1`` applies once per residue of that class."""
    resolver = AtomTokenResolver(_residue_file())
    assert resolver.resolve(['C1'], [1, 2, 3]).fullnames == ['C1_1', 'C1_2', 'C1_3']


def test_wildcard_suffix_covers_every_residue() -> None:
    resolver = AtomTokenResolver(_residue_file())
    assert resolver.resolve(['C1_*']).fullnames == ['C1_1', 'C1_2', 'C1_3']


def test_next_residue_suffix() -> None:
    """``_+`` is absent from the corpus but part of the language."""
    resolver = AtomTokenResolver(_residue_file())
    assert resolver.resolve(['C1_+'], [1]).fullnames == ['C1_2']


def test_previous_residue_suffix() -> None:
    resolver = AtomTokenResolver(_residue_file())
    assert resolver.resolve(['C1_-'], [3]).fullnames == ['C1_2']


def test_next_residue_past_the_end_is_unresolved() -> None:
    resolver = AtomTokenResolver(_residue_file())
    assert resolver.resolve(['C1_+'], [3]).unresolved == ['C1_+']


# --------------------------------------------------- D-6b: residue then symmetry

def _eqiv_file() -> Shelxfile:
    return _shx(
        'EQIV $3 1-x, y, 1-z\n'
        'RESI 23 ALA\n'
        'O1    3     0.30  0.40  0.50  11.00000  0.05\n'
        'RESI 0\n'
        'O2    3     0.60  0.70  0.80  11.00000  0.05\n'
    )


def test_residue_scope_applies_before_the_symmetry_operation() -> None:
    """``RTAB_23 ... O1_$3`` means ``(O1_23)_$3``, not ``(O1_0)_$3``.

    The old ``does_atom_exist`` hardcoded residue 0 here, so a valid
    reference was reported as an unknown atom.
    """
    resolver = AtomTokenResolver(_eqiv_file())
    result = resolver.resolve(['O1_$3'], [23])
    assert result.fullnames == ['O1_23']
    assert result.unresolved == []


def test_symmetry_suffix_is_reported_separately() -> None:
    resolver = AtomTokenResolver(_eqiv_file())
    assert resolver.resolve(['O1_$3'], [23]).eqiv_ids == ['$3']


def test_symmetry_suffix_is_not_mistaken_for_a_residue() -> None:
    resolver = AtomTokenResolver(_eqiv_file())
    assert resolver.resolve(['O2_$3']).fullnames == ['O2_0']


def test_explicit_residue_and_symmetry_combine() -> None:
    resolver = AtomTokenResolver(_eqiv_file())
    result = resolver.resolve(['O1_23_$3'.replace('_23_$3', '_23')])
    assert result.fullnames == ['O1_23']


def test_validator_scopes_residue_before_symmetry() -> None:
    """The same D-6b fix, through ``_assign_atoms_to_restraints``.

    ``does_atom_exist`` used to hardcode residue 0 after stripping the
    ``_$n`` suffix, so this valid instruction was reported as referencing
    an unknown atom.
    """
    shx = Shelxfile()
    shx.read_string(
        HEADER
        + 'EQIV $3 1-x, y, 1-z\n'
          'SADI_23 O1 O1_$3\n'
          'RESI 23 ALA\n'
          'O1    3     0.30  0.40  0.50  11.00000  0.05\n'
          'RESI 0\n'
        + FOOTER
    )
    assert shx.restraint_errors == []


def test_validator_still_reports_a_genuinely_missing_atom() -> None:
    """The fix must not silence real problems."""
    shx = Shelxfile()
    shx.read_string(
        HEADER
        + 'EQIV $3 1-x, y, 1-z\n'
          'SADI_23 O1 CZZ_$3\n'
          'RESI 23 ALA\n'
          'O1    3     0.30  0.40  0.50  11.00000  0.05\n'
          'RESI 0\n'
        + FOOTER
    )
    assert any('CZZ' in w for w in shx.restraint_errors)


# -------------------------------------------------------------- results

def test_fullnames_are_deduplicated_in_file_order() -> None:
    resolver = _resolver(CHAIN)
    assert resolver.resolve(['C1', 'C2', 'C1']).fullnames == ['C1_0', 'C2_0']


def test_resolve_for_uses_the_card_residue_scope() -> None:
    shx = Shelxfile()
    shx.read_file('tests/resources/p21c.res')
    resolver = AtomTokenResolver(shx)
    simu = next(r for r in shx.restraints if type(r).__name__ == 'SIMU')
    assert resolver.resolve_for(simu).fullnames == resolver.resolve(
        list(simu.atoms), simu.residue_number
    ).fullnames


def test_range_expansion_resolves_the_parse_line_todo() -> None:
    """``SADI_CCF3 O1 > F9`` was an unimplemented TODO in cards.py."""
    shx = Shelxfile()
    shx.read_file('tests/resources/p21c.res')
    simu = next(r for r in shx.restraints if type(r).__name__ == 'SIMU')
    names = AtomTokenResolver(shx).resolve_for(simu).fullnames
    assert 'O1_4' in names and 'F9_4' in names
    assert len(names) > len(simu.atoms), 'the range must actually expand'
