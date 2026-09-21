"""What a ``SAME`` instruction actually links (plan item B1).

The original edit plan modelled ``SAME`` as carrying "two
residue-correlated fragment lists (template + target), matched
positionally".  The manual says otherwise, and in two different ways:

* plain ``SAME`` / ``SAME_<n>`` is compared against *"the same number of
  atoms which follow the SAME instruction"* -- the counterpart fragment
  is not on the card at all;
* ``SAME_<class>`` *"no longer uses the following atoms but is applied to
  all residues with the name XYZ"*.

The practical consequence, and the reason this matters for deletion:
removing an atom that a ``SAME`` never mentions can still break it.
"""

from __future__ import annotations

from shelxfile import Shelxfile
from shelxfile.edit.card_meta import AtomGrouping
from shelxfile.edit.same_links import SameResolver

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


def _fragments(body: str):
    shx = _shx(body)
    card = next(r for r in shx.restraints if type(r).__name__ == 'SAME')
    return shx, SameResolver(shx).fragments(card)


#: The manual's own layout: template fragment, then SAME, then the target.
TWO_FRAGMENTS = (
    'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
    'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
    'C3    1     0.20  0.30  0.40  11.00000  0.05\n'
    'SAME C1 > C3\n'
    'C1B   1     0.40  0.50  0.60  11.00000  0.05\n'
    'H1    2     0.41  0.51  0.61  11.00000  0.05\n'
    'C2B   1     0.45  0.55  0.65  11.00000  0.05\n'
    'C3B   1     0.50  0.60  0.70  11.00000  0.05\n'
)


# --------------------------------------------- plain SAME: following atoms

def test_counterpart_is_the_run_of_atoms_after_the_card() -> None:
    _, frags = _fragments(TWO_FRAGMENTS)
    assert frags.mode is AtomGrouping.SAME_FOLLOWING_ATOMS
    assert [a.fullname for a in frags.fragments[0]] == ['C1_0', 'C2_0', 'C3_0']
    assert [a.fullname for a in frags.fragments[1]] == ['C1B_0', 'C2B_0', 'C3B_0']


def test_atoms_never_named_on_the_card_are_still_linked() -> None:
    """The whole point of B1: ``C2B`` appears nowhere on the instruction."""
    shx, frags = _fragments(TWO_FRAGMENTS)
    linked = [a.fullname for a in frags.all_atoms]
    assert 'C2B_0' in linked
    card = next(r for r in shx.restraints if type(r).__name__ == 'SAME')
    assert 'C2B' not in card.atoms


def test_counterparts_are_positional() -> None:
    shx, frags = _fragments(TWO_FRAGMENTS)
    c2 = shx.atoms.get_atom_by_name('C2_0')
    assert [a.fullname for a in frags.counterparts(c2)] == ['C2B_0']


def test_counterparts_work_from_either_fragment() -> None:
    shx, frags = _fragments(TWO_FRAGMENTS)
    c2b = shx.atoms.get_atom_by_name('C2B_0')
    assert [a.fullname for a in frags.counterparts(c2b)] == ['C2_0']


def test_hydrogens_are_ignored() -> None:
    """"Since hydrogen atoms are ignored by SAME"."""
    _, frags = _fragments(TWO_FRAGMENTS)
    assert all(not a.is_hydrogen for a in frags.all_atoms)
    assert 'H1_0' not in [a.fullname for a in frags.all_atoms]


def test_only_as_many_following_atoms_as_named_are_taken() -> None:
    _, frags = _fragments(
        TWO_FRAGMENTS + 'O1    3     0.90  0.90  0.90  11.00000  0.05\n'
    )
    assert 'O1_0' not in [a.fullname for a in frags.all_atoms]


def test_position_in_the_file_changes_the_meaning() -> None:
    """"The position of a SAME instruction in the input file is critical"."""
    _, early = _fragments(
        'SAME C1 C2\n'
        'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
        'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
    )
    _, late = _fragments(
        'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
        'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
        'SAME C1 C2\n'
        'C1B   1     0.40  0.50  0.60  11.00000  0.05\n'
        'C2B   1     0.45  0.55  0.65  11.00000  0.05\n'
    )
    assert [a.fullname for a in early.fragments[1]] != \
           [a.fullname for a in late.fragments[1]]


def test_missing_following_atoms_are_reported_as_truncated() -> None:
    _, frags = _fragments(
        'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
        'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
        'C3    1     0.20  0.30  0.40  11.00000  0.05\n'
        'SAME C1 > C3\n'
        'C1B   1     0.40  0.50  0.60  11.00000  0.05\n'
    )
    assert frags.truncated
    assert frags.width == 1


def test_residue_number_suffix_still_uses_following_atoms() -> None:
    _, frags = _fragments(
        'RESI 1 TOL\n'
        'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
        'RESI 0\n'
        'SAME_1 C1\n'
        'C1B   1     0.40  0.50  0.60  11.00000  0.05\n'
    )
    assert frags.mode is AtomGrouping.SAME_FOLLOWING_ATOMS


# ------------------------------------------------- SAME_<class>: residues

THREE_RESIDUES = (
    'SAME_TOL C1 > C3\n'
    'RESI 1 TOL\n'
    'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
    'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
    'C3    1     0.20  0.30  0.40  11.00000  0.05\n'
    'RESI 2 TOL\n'
    'C1    1     0.40  0.50  0.60  11.00000  0.05\n'
    'C2    1     0.45  0.55  0.65  11.00000  0.05\n'
    'C3    1     0.50  0.60  0.70  11.00000  0.05\n'
    'RESI 3 TOL\n'
    'C1    1     0.70  0.80  0.90  11.00000  0.05\n'
    'C2    1     0.75  0.85  0.95  11.00000  0.05\n'
    'C3    1     0.78  0.88  0.98  11.00000  0.05\n'
    'RESI 0\n'
)


def test_residue_class_mode_builds_one_fragment_per_residue() -> None:
    _, frags = _fragments(THREE_RESIDUES)
    assert frags.mode is AtomGrouping.SAME_RESIDUE_CLASS
    assert len(frags.fragments) == 3


def test_residue_class_mode_ignores_following_atoms() -> None:
    """``SAME_XYZ`` "no longer uses the following atoms"."""
    _, frags = _fragments(THREE_RESIDUES)
    residues = {a.resinum for a in frags.all_atoms}
    assert residues == {1, 2, 3}


def test_residue_class_counterparts_span_every_residue() -> None:
    shx, frags = _fragments(THREE_RESIDUES)
    c2 = shx.atoms.get_atom_by_name('C2_1')
    names = sorted(a.fullname for a in frags.counterparts(c2))
    assert names == ['C2_2', 'C2_3']


def test_residue_class_alignment_is_by_position() -> None:
    _, frags = _fragments(THREE_RESIDUES)
    assert [a.fullname for a in frags.column(0)] == ['C1_1', 'C1_2', 'C1_3']


def test_uneven_residues_are_reported_as_truncated() -> None:
    body = THREE_RESIDUES.replace(
        'RESI 3 TOL\n'
        'C1    1     0.70  0.80  0.90  11.00000  0.05\n'
        'C2    1     0.75  0.85  0.95  11.00000  0.05\n'
        'C3    1     0.78  0.88  0.98  11.00000  0.05\n',
        'RESI 3 TOL\n'
        'C1    1     0.70  0.80  0.90  11.00000  0.05\n'
        'C2    1     0.75  0.85  0.95  11.00000  0.05\n',
    )
    _, frags = _fragments(body)
    assert frags.truncated


# ----------------------------------------------------------------- shared

def test_unknown_atom_is_not_reported_as_a_counterpart() -> None:
    shx, frags = _fragments(TWO_FRAGMENTS)
    stranger = shx.atoms.get_atom_by_name('H1_0')
    assert frags.counterparts(stranger) == []


def test_all_atoms_is_deduplicated() -> None:
    _, frags = _fragments(TWO_FRAGMENTS)
    names = [a.fullname for a in frags.all_atoms]
    assert len(names) == len(set(names))


def test_width_is_the_shortest_fragment() -> None:
    _, frags = _fragments(TWO_FRAGMENTS)
    assert frags.width == 3


def test_real_file_same_card_resolves() -> None:
    shx = Shelxfile()
    shx.read_file('tests/resources/p21c.res')
    card = next(r for r in shx.restraints if type(r).__name__ == 'SAME')
    frags = SameResolver(shx).fragments(card)
    assert frags.mode is AtomGrouping.SAME_RESIDUE_CLASS
    assert frags.width > 0
    assert all(not a.is_hydrogen for a in frags.all_atoms)
