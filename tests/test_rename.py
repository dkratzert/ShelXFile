"""Renaming an atom so the instructions follow (user request).

Setting ``atom.name`` on its own updated nothing else: restraints kept
the old name, and the lookup dict went stale -- the *old* name still
resolved while the new one returned ``None``.

Renaming has a subtlety deletion does not.  Deleting an atom that a
residue-class-scoped card still resolves elsewhere is harmless, because
SHELXL "simply ignores" the instruction for that residue.  *Renaming* the
token would instead redirect every residue, so those references are
deliberately left alone and reported.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit import ShelxDocument

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C H\n'
    'UNIT 8 8\n'
)
FOOTER = 'HKLF 4\nEND\n'
ATOMS = ''.join(
    f'C{i:<4} 1     0.{i:02d}0  0.{i:02d}5  0.{i:02d}9  11.00000  0.05\n'
    for i in range(1, 6)
)
RESIDUES = (
    'RESI 1 TOL\n'
    'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
    'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
    'RESI 2 TOL\n'
    'C1    1     0.40  0.50  0.60  11.00000  0.05\n'
    'C2    1     0.45  0.55  0.65  11.00000  0.05\n'
    'RESI 0\n'
)


def _doc(cards: str = '', atoms: str = ATOMS) -> ShelxDocument:
    return ShelxDocument.from_string(HEADER + cards + atoms + FOOTER)


# ------------------------------------------------------------ the basics

def test_the_atom_is_renamed() -> None:
    doc = _doc()
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    assert report.ok
    assert report.new_name == 'CX_0'


def test_lookup_follows_the_new_name() -> None:
    """The dict is keyed on fullname and used to go stale."""
    doc = _doc()
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    assert doc.shelxfile.atoms.get_atom_by_name('CX_0') is not None
    assert doc.shelxfile.atoms.get_atom_by_name('C1_0') is None


def test_the_rendered_file_shows_the_new_name() -> None:
    doc = _doc()
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    assert any(line.startswith('CX') for line in doc.text.splitlines())


def test_renaming_to_the_same_name_is_a_no_op() -> None:
    doc = _doc('SADI C1 C2 C3 C4\n')
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'C1')
    assert report.ok
    assert report.updated == []


# --------------------------------------------------- instructions follow

def test_restraints_are_updated() -> None:
    doc = _doc('SADI C1 C2 C3 C4\n')
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    assert 'SADI CX C2 C3 C4' in doc.text


def test_several_cards_are_updated() -> None:
    doc = _doc('SADI C1 C2 C3 C4\nDFIX 1.5 C1 C2\n')
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    assert len(report.updated) == 2
    assert 'DFIX 1.5 CX C2' in doc.text


def test_the_report_shows_the_old_line() -> None:
    doc = _doc('SADI C1 C2 C3 C4\n')
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    assert report.updated[0].before == 'SADI C1 C2 C3 C4'


def test_named_operand_cards_follow() -> None:
    """``FREE``/``HTAB`` keep operands in fields, not a list."""
    doc = _doc('FREE C1 C2\n')
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    assert doc.shelxfile.free[0].referenced_atoms == ['CX', 'C2']


def test_a_range_endpoint_is_updated() -> None:
    doc = _doc('SIMU C1 > C4\n')
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    assert 'SIMU CX > C4' in doc.text


def test_an_atom_inside_a_range_needs_no_change() -> None:
    """A range is positional, so the middle is never spelled out."""
    doc = _doc('SIMU C1 > C4\n')
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C2_0'), 'CX')
    assert report.updated == []
    assert 'SIMU C1 > C4' in doc.text


def test_unreferenced_atom_renames_cleanly() -> None:
    doc = _doc('SADI C1 C2 C3 C4\n')
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C5_0'), 'CX')
    assert report.ok
    assert report.updated == []


# ----------------------------------- the residue-class case is left alone

def test_a_token_shared_across_residues_is_not_rewritten() -> None:
    """Rewriting it would redirect every residue of the class."""
    doc = _doc('SADI_TOL C1 C2\n', atoms=RESIDUES)
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_2'), 'CX')
    assert report.ok
    assert report.updated == []
    assert 'SADI_TOL C1 C2' in doc.text


def test_the_untouched_reference_is_reported() -> None:
    doc = _doc('SADI_TOL C1 C2\n', atoms=RESIDUES)
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_2'), 'CX')
    assert len(report.skipped) == 1
    assert report.skipped[0].token == 'C1'
    assert 'other residues' in report.skipped[0].reason
    assert 'left unchanged' in report.summary()


def test_the_other_residue_keeps_working() -> None:
    doc = _doc('SADI_TOL C1 C2\n', atoms=RESIDUES)
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_2'), 'CX')
    reparsed = Shelxfile()
    reparsed.read_string(doc.text)
    assert reparsed.atoms.get_atom_by_name('C1_1') is not None
    assert [w for w in reparsed.restraint_errors if 'Unknown atom' in w] == []


def test_the_last_residue_instance_is_rewritten() -> None:
    """With only one residue left, the token is no longer shared."""
    doc = _doc(
        'SADI_TOL C1 C2\n',
        atoms=(
            'RESI 1 TOL\n'
            'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
            'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
            'RESI 0\n'
        ),
    )
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_1'), 'CX')
    assert len(report.updated) == 1
    assert 'SADI_TOL CX C2' in doc.text


# ------------------------------------------------------------ validation

@pytest.mark.parametrize('bad, fragment', [
    ('C2', 'already used'),          # uniqueness
    ('1X', 'must start with a letter'),
    ('TOOLONG', 'longer than 4'),
    ('', 'must not be empty'),
    ('C1_2', 'without a residue number'),
])
def test_invalid_names_are_refused(bad, fragment) -> None:
    doc = _doc()
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), bad)
    assert not report.ok
    assert fragment in report.error


def test_a_refused_rename_changes_nothing() -> None:
    doc = _doc('SADI C1 C2 C3 C4\n')
    before = doc.text
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'C2')
    assert doc.text == before


def test_the_same_name_in_another_residue_is_allowed() -> None:
    """Uniqueness is per (name, PART, RESI), not global."""
    doc = _doc('', atoms=RESIDUES)
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C2_2'), 'C9')
    assert report.ok


def test_symmetry_generated_atoms_cannot_be_renamed() -> None:
    doc = ShelxDocument.from_file('tests/resources/p21c.res')
    packed = [a for a in doc.shelxfile.pack() if a.symmgen]
    assert packed
    report = doc.rename_atom(packed[0], 'CX')
    assert not report.ok


# ------------------------------------------------------------- round-trip

def test_a_rename_survives_a_round_trip() -> None:
    doc = _doc('SADI C1 C2 C3 C4\nDFIX 1.5 C1 C2\n')
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    reparsed = Shelxfile()
    reparsed.read_string(doc.text)
    assert reparsed.atoms.get_atom_by_name('CX_0') is not None
    assert [w for w in reparsed.restraint_errors if 'Unknown atom' in w] == []


def test_renaming_leaves_no_dangling_references() -> None:
    """The point of the exercise."""
    doc = ShelxDocument.from_file('tests/resources/p21c.res')
    before = len([w for w in doc.shelxfile.restraint_errors
                  if 'Unknown atom' in w])
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('O1_4'), 'OX')
    reparsed = Shelxfile()
    reparsed.read_string(doc.text)
    after = len([w for w in reparsed.restraint_errors if 'Unknown atom' in w])
    assert after <= before


def test_observers_are_notified() -> None:
    doc = _doc('SADI C1 C2 C3 C4\n')
    seen: list[object] = []
    doc.subscribe(seen.append)
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1_0'), 'CX')
    assert seen == [doc]
