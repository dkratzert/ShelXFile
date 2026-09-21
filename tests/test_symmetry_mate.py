"""Provenance for symmetry-generated atoms (plan item D-7a).

``grow()`` and ``pack()`` invent a fresh name for every image they
produce.  That is enough to draw a picture but not to edit one: a viewer
letting the user click a symmetry mate and ask for a bond needs to know
*which* operation produced it, so ShelXFile can write ``BIND C1 C2_$3``
against the right ``EQIV``.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.atoms.symmetry_mate import SymmetryMate

RES = 'tests/resources/p21c.res'


@pytest.fixture(scope='module')
def packed() -> tuple[Shelxfile, list]:
    shx = Shelxfile()
    shx.read_file(RES)
    return shx, shx.pack()


# ------------------------------------------------------------- recording

def test_generated_atoms_carry_provenance(packed) -> None:
    _, atoms = packed
    mates = [a for a in atoms if a.symmgen]
    assert mates
    assert all(a.symm_mate is not None for a in mates)


def test_asymmetric_unit_atoms_carry_none(packed) -> None:
    """A real line in the file is not an image of anything."""
    _, atoms = packed
    assert all(a.symm_mate is None for a in atoms if not a.symmgen)


def test_provenance_names_the_parent(packed) -> None:
    shx, atoms = packed
    mate = next(a for a in atoms if a.symm_mate).symm_mate
    assert mate.parent in shx.atoms.all_atoms


def test_provenance_names_the_operation(packed) -> None:
    _, atoms = packed
    mate = next(a for a in atoms if a.symm_mate).symm_mate
    assert mate.symm_number > 0
    assert 0 <= mate.symm_number < len(_symmcards(atoms))


def _symmcards(atoms):
    return atoms[0].shx.symmcards


def test_grow_records_provenance_including_the_lattice_shift() -> None:
    shx = Shelxfile()
    shx.read_file('tests/resources/014EP-4_a_shelxl.res')
    mates = [a for a in shx.grow() if a.symm_mate]
    assert mates, 'fixture should grow symmetry mates'
    assert any(m.symm_mate.hkl_shift != (0, 0, 0) for m in mates), (
        'a grown structure should include lattice-shifted images'
    )


# ---------------------------------------------------------------- labels

def test_label_is_in_shelxl_form(packed) -> None:
    _, atoms = packed
    label = next(a for a in atoms if a.symm_mate).symm_mate.symm_label
    assert ',' in label
    assert any(axis in label for axis in 'xyz')


def test_identity_has_no_label() -> None:
    shx = Shelxfile()
    shx.read_file(RES)
    mate = SymmetryMate(parent=shx.atoms.all_atoms[0])
    assert mate.is_identity
    assert mate.symm_label == ''


def test_a_shift_changes_the_label() -> None:
    shx = Shelxfile()
    shx.read_file(RES)
    parent = shx.atoms.all_atoms[0]
    plain = SymmetryMate(parent, symm_number=1)
    shifted = SymmetryMate(parent, symm_number=1, hkl_shift=(1, 0, 0))
    assert plain.symm_label != shifted.symm_label


def test_label_matches_the_bond_convention(packed) -> None:
    """Atom clicks and bond clicks must describe an operation alike."""
    from shelxfile.atoms.atoms import _build_symm_label

    shx, atoms = packed
    mate = next(a for a in atoms if a.symm_mate).symm_mate
    expected = _build_symm_label(shx.symmcards[mate.symm_number], *mate.hkl_shift)
    assert mate.symm_label == expected


def test_str_is_readable(packed) -> None:
    _, atoms = packed
    mate = next(a for a in atoms if a.symm_mate).symm_mate
    text = str(mate)
    assert mate.parent.fullname_short in text
    assert '[' in text


def test_an_out_of_range_operation_degrades_quietly() -> None:
    shx = Shelxfile()
    shx.read_file(RES)
    mate = SymmetryMate(parent=shx.atoms.all_atoms[0], symm_number=999)
    assert mate.symm_label == ''


# ------------------------------------------------------------- behaviour

def test_provenance_does_not_make_mates_deletable(packed) -> None:
    """They are still not lines in the file."""
    shx, atoms = packed
    mate = next(a for a in atoms if a.symm_mate)
    before = len(shx._reslist)
    mate.delete()
    assert len(shx._reslist) == before


def test_mates_are_value_objects() -> None:
    """Equality by content, so two descriptions of the same image match.

    Not hashable: :class:`Atom` defines ``__eq__`` without ``__hash__``,
    so it cannot be a member of a hashed value object.
    """
    shx = Shelxfile()
    shx.read_file(RES)
    parent = shx.atoms.all_atoms[0]
    assert SymmetryMate(parent, 2, (1, 0, 0)) == SymmetryMate(parent, 2, (1, 0, 0))
    assert SymmetryMate(parent, 2, (1, 0, 0)) != SymmetryMate(parent, 3, (1, 0, 0))
