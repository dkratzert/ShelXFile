"""Finding or minting the ``EQIV`` a symmetry reference needs (plan item D-7b).

This is the bridge a viewer needs: the user clicks a symmetry image, and
ShelXFile has to name it in an instruction -- ``C2_$3`` -- which means
knowing, or deciding, that ``$3`` is that operation.

Matching is on the parsed operation rather than its text, so a file does
not collect several numbers for the same symmetry.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit.eqiv_factory import EqivFactory, canonical_symmop

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C\n'
    'UNIT 4\n'
)
ATOMS = (
    'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
    'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
)
FOOTER = 'HKLF 4\nEND\n'


def _factory(body: str = '') -> tuple[Shelxfile, EqivFactory]:
    shx = Shelxfile()
    shx.read_string(HEADER + body + ATOMS + FOOTER)
    return shx, EqivFactory(shx)


# ------------------------------------------------------- canonical form

@pytest.mark.parametrize('a, b', [
    ('1-x, y, 1-z', '-x+1, +y, -z+1'),
    ('-X, 0.5+Y, 0.5-Z', '-x, y+1/2, -z+1/2'),
    ('x, y, z', 'X, Y, Z'),
])
def test_equivalent_spellings_match(a, b) -> None:
    """Text differs, operation does not."""
    assert canonical_symmop(a) == canonical_symmop(b)


def test_different_operations_do_not_match() -> None:
    assert canonical_symmop('1-x, y, 1-z') != canonical_symmop('x, 1-y, z')


def test_a_lattice_translation_is_not_a_new_operation() -> None:
    assert canonical_symmop('x, y, z') == canonical_symmop('x+1, y, z')


@pytest.mark.parametrize('bad', ['', 'nonsense', 'x, y'])
def test_unparsable_input_gives_none(bad) -> None:
    assert canonical_symmop(bad) is None


# -------------------------------------------------------------- lookup

def test_an_existing_eqiv_is_found() -> None:
    shx, factory = _factory('EQIV $1 1-x, y, 1-z\n')
    assert factory.find('1-x, y, 1-z') is shx.eqiv[0]


def test_lookup_ignores_how_it_was_written() -> None:
    shx, factory = _factory('EQIV $1 1-x, y, 1-z\n')
    assert factory.find('-x+1, +y, -z+1') is shx.eqiv[0]


def test_an_unknown_operation_is_not_found() -> None:
    _, factory = _factory('EQIV $1 1-x, y, 1-z\n')
    assert factory.find('x, 1-y, z') is None


def test_post_end_definitions_are_not_reused() -> None:
    """That block records a finished refinement; new input must not lean on it."""
    shx = Shelxfile()
    shx.read_string(HEADER + ATOMS + FOOTER + 'EQIV $9 1-x, y, 1-z\n')
    assert EqivFactory(shx).find('1-x, y, 1-z') is None


# ---------------------------------------------------------------- mint

def test_a_missing_operation_is_created() -> None:
    shx, factory = _factory()
    card = factory.eqiv_for('1-x, y, 1-z')
    assert card is not None
    assert card.id == '$1'
    assert 'EQIV $1 1-x, y, 1-z' in shx.dumps()


def test_the_same_operation_is_not_created_twice() -> None:
    shx, factory = _factory()
    first = factory.eqiv_for('1-x, y, 1-z')
    second = factory.eqiv_for('-x+1, +y, -z+1')
    assert first is second
    assert len(shx.eqiv) == 1


def test_numbers_are_handed_out_from_the_unused_ones() -> None:
    shx, factory = _factory('EQIV $1 1-x, y, 1-z\nEQIV $3 x, 1-y, z\n')
    card = factory.eqiv_for('-x, -y, -z')
    assert card.id == '$2', 'the gap should be filled before going past $3'


def test_existing_numbers_are_never_reassigned() -> None:
    shx, factory = _factory('EQIV $1 1-x, y, 1-z\n')
    factory.eqiv_for('x, 1-y, z')
    assert [c.id for c in shx.eqiv] == ['$1', '$2']
    assert 'EQIV $1 1-x, y, 1-z' in shx.dumps()


def test_create_false_does_not_write_anything() -> None:
    shx, factory = _factory()
    assert factory.eqiv_for('1-x, y, 1-z', create=False) is None
    assert shx.eqiv == []


def test_an_unparsable_operation_is_refused() -> None:
    shx, factory = _factory()
    assert factory.eqiv_for('nonsense') is None
    assert shx.eqiv == []


# ------------------------------------------------------------ position

def test_a_new_eqiv_goes_before_the_atoms() -> None:
    """"Such a symmetry operation must be defined before it is used"."""
    shx, factory = _factory()
    factory.eqiv_for('1-x, y, 1-z')
    lines = shx.dumps().splitlines()
    eqiv_line = next(i for i, ln in enumerate(lines) if ln.startswith('EQIV'))
    atom_line = next(i for i, ln in enumerate(lines) if ln.startswith('C1'))
    assert eqiv_line < atom_line


def test_new_definitions_stay_with_the_existing_ones() -> None:
    shx, factory = _factory('EQIV $1 1-x, y, 1-z\n')
    factory.eqiv_for('x, 1-y, z')
    lines = [ln for ln in shx.dumps().splitlines() if ln.startswith('EQIV')]
    positions = [i for i, ln in enumerate(shx.dumps().splitlines())
                 if ln.startswith('EQIV')]
    assert len(lines) == 2
    assert positions[1] == positions[0] + 1


def test_the_result_reparses() -> None:
    shx, factory = _factory()
    factory.eqiv_for('1-x, y, 1-z')
    reread = Shelxfile()
    reread.read_string(shx.dumps())
    assert [c.id for c in reread.eqiv] == ['$1']


# ----------------------------------------------------- naming a mate

def test_an_asymmetric_unit_atom_needs_no_suffix() -> None:
    shx, factory = _factory()
    atom = shx.atoms.get_atom_by_name('C1_0')
    assert factory.symmetry_atom_name(atom) == 'C1'


def test_a_symmetry_image_is_named_with_its_eqiv() -> None:
    shx = Shelxfile()
    shx.read_file('tests/resources/p21c.res')
    factory = EqivFactory(shx)
    mate = next(a for a in shx.pack() if a.symmgen)

    name = factory.symmetry_atom_name(mate)

    assert name is not None
    assert '_$' in name
    assert name.startswith(mate.symm_mate.parent.fullname_short)


def test_naming_reuses_one_eqiv_per_operation() -> None:
    shx = Shelxfile()
    shx.read_file('tests/resources/p21c.res')
    factory = EqivFactory(shx)
    mates = [a for a in shx.pack() if a.symmgen]
    same_op = [m for m in mates
               if m.symm_mate.symm_number == mates[0].symm_mate.symm_number]
    assert len(same_op) > 1

    names = {factory.symmetry_atom_name(m) for m in same_op}
    suffixes = {n.rsplit('_', 1)[1] for n in names}

    assert len(suffixes) == 1, 'one operation deserves exactly one EQIV'


def test_naming_without_creating_returns_none_when_undefined() -> None:
    shx = Shelxfile()
    shx.read_file('tests/resources/p21c.res')
    factory = EqivFactory(shx)
    mate = next(a for a in shx.pack() if a.symmgen)
    assert factory.symmetry_atom_name(mate, create=False) is None
    assert shx.eqiv == []
