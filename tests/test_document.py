"""The Qt-free editing façade (plan decision D-8).

``ShelxDocument`` is the only object allowed to reach into ``_reslist``.
A view asks it what sits on a line, tells it what the user did, and
renders the result -- it never parses SHELXL or decides card semantics.

Also covers **D-9**: deleting an atom may cascade, but removing a card
never deletes atoms unless the explicitly named destructive method is
used.
"""

from __future__ import annotations

import pytest

from shelxfile.edit import (
    DeletionReport,
    RemovalReason,
    ShelxDocument,
)
from shelxfile.atoms.atom import Atom
from shelxfile.shelx.cards import Restraint

RES = 'tests/resources/p21c.res'


@pytest.fixture
def doc() -> ShelxDocument:
    return ShelxDocument.from_file(RES)


# ---------------------------------------------------------------- rendering

def test_text_matches_the_model_dump(doc: ShelxDocument) -> None:
    """If these diverge, editor cursor positions point at the wrong card."""
    assert doc.text == doc.shelxfile.dumps()


def test_line_count_matches_the_text(doc: ShelxDocument) -> None:
    assert doc.line_count == len(doc.text.splitlines())


def test_from_string_builds_a_usable_document() -> None:
    doc = ShelxDocument.from_string(
        'TITL t\nCELL 0.71073 10 10 10 90 90 90\nZERR 4 0 0 0 0 0 0\n'
        'LATT 1\nSFAC C\nUNIT 4\n'
        'C1    1     0.1  0.2  0.3  11.00000  0.05\n'
        'HKLF 4\nEND\n'
    )
    assert doc.line_of_atom('C1') is not None


# ------------------------------------------------------------ line lookups

def test_item_at_line_round_trips_with_line_of(doc: ShelxDocument) -> None:
    atom = doc.shelxfile.atoms.all_atoms[3]
    line = doc.line_of(atom)
    assert line is not None
    assert doc.item_at_line(line) is atom


def test_item_at_line_returns_none_out_of_range(doc: ShelxDocument) -> None:
    assert doc.item_at_line(-1) is None
    assert doc.item_at_line(doc.line_count + 10) is None


def test_wrapped_cards_map_every_line_to_one_item(doc: ShelxDocument) -> None:
    """Long atom lines wrap; all their lines belong to the same entry."""
    atom = doc.shelxfile.atoms.all_atoms[0]
    line = doc.line_of(atom)
    assert doc.item_at_line(line) is atom
    assert doc.item_at_line(line + 1) is atom, 'fixture atom should wrap'


def test_items_in_lines_deduplicates_wrapped_entries(doc: ShelxDocument) -> None:
    atom = doc.shelxfile.atoms.all_atoms[0]
    line = doc.line_of(atom)
    items = doc.items_in_lines(line, line + 1)
    assert items.count(atom) == 1


def test_items_in_lines_follows_file_order(doc: ShelxDocument) -> None:
    atoms = doc.shelxfile.atoms.all_atoms[:3]
    first = doc.line_of(atoms[0])
    last = doc.line_of(atoms[2])
    found = [i for i in doc.items_in_lines(first, last) if isinstance(i, Atom)]
    assert found == atoms


def test_line_of_atom_finds_by_short_name(doc: ShelxDocument) -> None:
    line = doc.line_of_atom('C1_4')
    assert line is not None
    item = doc.item_at_line(line)
    assert isinstance(item, Atom) and item.fullname_short == 'C1_4'


def test_line_of_unknown_item_is_none(doc: ShelxDocument) -> None:
    assert doc.line_of_atom('ZZ99') is None
    assert doc.line_of(Atom(doc.shelxfile)) is None


# ------------------------------------------------------------- observation

def test_observers_fire_after_an_edit(doc: ShelxDocument) -> None:
    seen: list[ShelxDocument] = []
    doc.subscribe(seen.append)
    doc.delete_atom(doc.shelxfile.atoms.all_atoms[0])
    assert seen == [doc]


def test_observers_are_registered_once(doc: ShelxDocument) -> None:
    seen: list[object] = []
    callback = seen.append
    doc.subscribe(callback)
    doc.subscribe(callback)
    doc.delete_atom(doc.shelxfile.atoms.all_atoms[0])
    assert len(seen) == 1


def test_unsubscribe_stops_notifications(doc: ShelxDocument) -> None:
    seen: list[object] = []
    doc.subscribe(seen.append)
    doc.unsubscribe(seen.append)
    doc.delete_atom(doc.shelxfile.atoms.all_atoms[0])
    assert seen == []


def test_text_is_refreshed_after_an_edit(doc: ShelxDocument) -> None:
    before = doc.text
    doc.delete_atom(doc.shelxfile.atoms.all_atoms[0])
    assert doc.text != before
    assert doc.text == doc.shelxfile.dumps()


# ---------------------------------------------------------------- deletion

def test_delete_atom_reports_what_it_removed(doc: ShelxDocument) -> None:
    atom = doc.shelxfile.atoms.all_atoms[0]
    report = doc.delete_atom(atom)
    assert len(report.atoms) == 1
    assert report.atoms[0].item is atom
    assert report.atoms[0].reason is RemovalReason.REQUESTED
    assert 'atom' in report.summary()


def test_delete_atoms_skips_symmetry_generated(doc: ShelxDocument) -> None:
    packed = doc.shelxfile.pack()
    symmgen = [a for a in packed if a.symmgen][:3]
    assert symmgen, 'precondition: pack() produced symmetry mates'
    report = doc.delete_atoms(symmgen)
    assert report.is_empty


def test_empty_deletion_does_not_notify(doc: ShelxDocument) -> None:
    seen: list[object] = []
    doc.subscribe(seen.append)
    doc.delete_atoms([])
    assert seen == []


# ----------------------------------------- D-2/D-4: AFIX group cascades

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C H\n'
    'UNIT 8 8\n'
)
FOOTER = 'HKLF 4\nEND\n'


def _atom_line(name: str, index: int, sfac: int = 1) -> str:
    return (f'{name:<5} {sfac}     0.{index:02d}0  0.{index:02d}5  '
            f'0.{index:02d}9  11.00000  0.05\n')


def _hexagon_doc() -> ShelxDocument:
    body = _atom_line('C1', 1) + 'AFIX 66\n' + ''.join(
        _atom_line(f'C{i}', i) for i in range(2, 8)
    ) + 'AFIX 0\n'
    return ShelxDocument.from_string(HEADER + body + FOOTER)


def _riding_doc() -> ShelxDocument:
    body = _atom_line('C1', 1) + 'AFIX 43\n' + _atom_line('H1', 2, 2) + 'AFIX 0\n'
    return ShelxDocument.from_string(HEADER + body + FOOTER)


def test_under_populated_afix_takes_its_members_with_it() -> None:
    """D-2: five atoms cannot be fitted to a six-atom ring."""
    doc = _hexagon_doc()
    doc.delete_atom(doc.shelxfile.atoms.get_atom_by_name('C4_0'))
    assert [a.name for a in doc.shelxfile.atoms] == ['C1'], (
        'the whole rigid group should be gone, leaving only the pivot atom'
    )


def test_under_populated_afix_removes_its_cards() -> None:
    doc = _hexagon_doc()
    report = doc.delete_atom(doc.shelxfile.atoms.get_atom_by_name('C4_0'))
    reasons = {c.reason for c in report.cards}
    assert reasons == {RemovalReason.AFIX_UNDER_POPULATED}
    assert 'AFIX 66' not in doc.text


def test_cascade_attributes_every_removal() -> None:
    doc = _hexagon_doc()
    report = doc.delete_atom(doc.shelxfile.atoms.get_atom_by_name('C4_0'))
    requested = [a for a in report.atoms if a.reason is RemovalReason.REQUESTED]
    collateral = [a for a in report.atoms
                  if a.reason is RemovalReason.AFIX_UNDER_POPULATED]
    assert len(requested) == 1
    assert len(collateral) == 5


def test_removing_an_atom_from_a_full_group_keeps_it_when_allowed() -> None:
    """A group with spare members is not invalidated."""
    body = _atom_line('C1', 1) + 'AFIX 43\n' + _atom_line('H1', 2, 2) \
        + _atom_line('H2', 3, 2) + 'AFIX 0\n'
    doc = ShelxDocument.from_string(HEADER + body + FOOTER)
    doc.delete_atom(doc.shelxfile.atoms.get_atom_by_name('H1_0'))
    assert 'AFIX 43' in doc.text
    assert [a.name for a in doc.shelxfile.atoms] == ['C1', 'H2']


def test_deleting_a_pivot_removes_its_riding_group() -> None:
    """D-4: the pivot sits outside the bracket, so the group is orphaned."""
    doc = _riding_doc()
    doc.delete_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'))
    assert [a.name for a in doc.shelxfile.atoms] == []
    assert 'AFIX 43' not in doc.text


def test_deleting_a_carbon_removes_its_hydrogens() -> None:
    """The behaviour a user expects, stated plainly."""
    doc = _riding_doc()
    report = doc.delete_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'))
    removed = {a.item.name for a in report.atoms}
    assert removed == {'C1', 'H1'}
    assert any(a.reason is RemovalReason.AFIX_PIVOT_DELETED for a in report.atoms)


def test_cascade_terminates(doc: ShelxDocument) -> None:
    """A chain of dependent groups must settle, not loop."""
    body = (
        _atom_line('C1', 1)
        + 'AFIX 43\n' + _atom_line('H1', 2, 2) + 'AFIX 0\n'
        + _atom_line('C2', 3)
        + 'AFIX 43\n' + _atom_line('H2', 4, 2) + 'AFIX 0\n'
    )
    chained = ShelxDocument.from_string(HEADER + body + FOOTER)
    report = chained.delete_atom(chained.shelxfile.atoms.get_atom_by_name('C1_0'))
    assert {a.item.name for a in report.atoms} == {'C1', 'H1'}
    assert [a.name for a in chained.shelxfile.atoms] == ['C2', 'H2']


def test_no_atom_is_reported_twice() -> None:
    doc = _hexagon_doc()
    report = doc.delete_atom(doc.shelxfile.atoms.get_atom_by_name('C4_0'))
    names = [a.item.name for a in report.atoms]
    assert len(names) == len(set(names))


def test_resets_are_not_treated_as_groups() -> None:
    """``AFIX 0`` on its own must never trigger a cascade."""
    body = _atom_line('C1', 1) + 'AFIX 0\n' + _atom_line('C2', 2)
    doc = ShelxDocument.from_string(HEADER + body + FOOTER)
    doc.delete_atom(doc.shelxfile.atoms.get_atom_by_name('C2_0'))
    assert [a.name for a in doc.shelxfile.atoms] == ['C1']


# ------------------------------------------------------ D-9: card removal

def _restraint_with_atoms(doc: ShelxDocument) -> Restraint:
    for restraint in doc.shelxfile.restraints:
        if [a for a in restraint.atoms if a not in ('>', '<', '=')]:
            return restraint
    raise AssertionError('fixture has no atom-bearing restraint')


def test_remove_restraint_keeps_the_atoms(doc: ShelxDocument) -> None:
    """D-9b - card removal must never reach an atom."""
    restraint = _restraint_with_atoms(doc)
    before = len(doc.shelxfile.atoms)

    report = doc.remove_restraint(restraint)

    assert len(doc.shelxfile.atoms) == before
    assert report.atoms == []
    assert len(report.cards) == 1


def test_delete_restraint_with_atoms_is_destructive(doc: ShelxDocument) -> None:
    """The opt-in variant, named so it cannot be triggered by accident."""
    restraint = _restraint_with_atoms(doc)
    before = len(doc.shelxfile.atoms)

    report = doc.delete_restraint_with_atoms(restraint)

    assert len(doc.shelxfile.atoms) < before
    assert report.atoms, 'atoms should have been removed'
    assert len(report.cards) == 1


def test_the_safe_method_has_no_destructive_flag() -> None:
    """A boolean switch would let callers destroy atoms unintentionally."""
    import inspect
    params = inspect.signature(ShelxDocument.remove_restraint).parameters
    assert list(params) == ['self', 'restraint']


# ----------------------------------------------------------------- reports

def test_deletion_reports_merge(doc: ShelxDocument) -> None:
    first = doc.delete_atom(doc.shelxfile.atoms.all_atoms[0])
    second = doc.delete_atom(doc.shelxfile.atoms.all_atoms[0])
    first.extend(second)
    assert len(first) == 2


def test_empty_report_summary_is_explicit() -> None:
    assert DeletionReport().summary() == 'nothing removed'
    assert DeletionReport().is_empty


# ------------------------------------------------------------------ output

def test_write_roundtrips_through_the_document(doc: ShelxDocument, tmp_path) -> None:
    doc.delete_atom(doc.shelxfile.atoms.all_atoms[0])
    out = tmp_path / 'out.res'
    doc.write(out)
    assert ShelxDocument.from_file(out).text == doc.text
