"""Structural edits: moving atoms, disorder parts, duplication and naming.

These are the primitives a model builder (Fastmolwidget's disorder drag,
for one) needs to write a split disorder back into a ``.res`` file.  Every
test re-reads the written text, because the only thing that matters in the
end is what SHELXL will see.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit import ShelxDocument

RES = 'tests/resources/jkd77.res'
METHYL = ('C10', 'H10A', 'H10B', 'H10C')


@pytest.fixture
def doc() -> ShelxDocument:
    return ShelxDocument.from_file(RES)


def atom(doc: ShelxDocument, name: str):
    found = doc.shelxfile.atoms.get_atom_by_name(name)
    assert found is not None, name
    return found


def reparse(doc: ShelxDocument) -> Shelxfile:
    shx = Shelxfile(debug=True)
    shx.read_string(doc.text)
    return shx


def lines_between(doc: ShelxDocument, first: str, count: int) -> list[str]:
    """First word of *count* lines starting at *first*, continuations skipped."""
    lines = [line.split()[0] for line in doc.text.splitlines()
             if line.strip() and not line[0].isspace()]
    start = lines.index(first)
    return lines[start:start + count]


# ------------------------------------------------------------ moving

class TestMoveAtoms:
    def test_cartesian_move_round_trips(self, doc):
        c1 = atom(doc, 'C1')
        target = (c1.xc + 0.5, c1.yc - 0.25, c1.zc + 0.1)
        report = doc.move_atoms([(c1, target)])
        assert report.changed == [c1]
        moved = reparse(doc).atoms.get_atom_by_name('C1')
        assert moved.cart_coords == pytest.approx(target, abs=1e-4)

    def test_fractional_move(self, doc):
        c1 = atom(doc, 'C1')
        doc.move_atoms([(c1, (0.1, 0.2, 0.3))], cartesian=False)
        assert reparse(doc).atoms.get_atom_by_name('C1').frac_coords == pytest.approx((0.1, 0.2, 0.3))

    def test_moving_a_fixed_coordinate_releases_it(self):
        doc = ShelxDocument.from_file(RES)
        c1 = atom(doc, 'C1')
        c1._coord_codes = (1, 0, 0)  # x written as 10.xxxx: held fixed
        report = doc.move_atoms([(c1, (c1.x + 0.01, c1.y, c1.z))], cartesian=False)
        assert c1._coord_codes == (0, 0, 0)
        assert report.messages and 'x' in report.messages[0]

    def test_unmoved_fixed_coordinate_keeps_its_code(self):
        doc = ShelxDocument.from_file(RES)
        c1 = atom(doc, 'C1')
        c1._coord_codes = (1, 0, 0)
        doc.move_atoms([(c1, (c1.x, c1.y + 0.01, c1.z))], cartesian=False)
        assert c1._coord_codes == (1, 0, 0)

    def test_is_one_undo_step(self, doc):
        c1, c2 = atom(doc, 'C1'), atom(doc, 'C2')
        doc.move_atoms([(c1, (0.1, 0.1, 0.1)), (c2, (0.2, 0.2, 0.2))], cartesian=False)
        assert doc.history.undo_labels == ['Move atoms']


# ---------------------------------------------------- free variables

def test_add_free_variable_appends_to_fvar(doc):
    number = doc.add_free_variable(0.3)
    assert number == 2
    shx = reparse(doc)
    assert len(shx.fvars) == 2
    assert shx.fvars[2] == pytest.approx(0.3)


# -------------------------------------------------------------- parts

class TestAssignPart:
    def test_brackets_include_the_riding_group(self, doc):
        doc.assign_part([atom(doc, name) for name in METHYL], 1, 10.5)
        assert lines_between(doc, 'PART', 8) == [
            'PART', 'C10', 'AFIX', 'H10A', 'H10B', 'H10C', 'AFIX', 'PART']
        shx = reparse(doc)
        for name in METHYL:
            reread = shx.atoms.get_atom_by_name(name)
            assert reread.part.n == 1
            assert reread.sof == 10.5
        assert shx.atoms.get_atom_by_name('C11').part.n == 0

    def test_partial_afix_group_is_refused(self, doc):
        with pytest.raises(ValueError, match='partly selected'):
            doc.assign_part([atom(doc, 'C10'), atom(doc, 'H10A')], 1, 21.0)
        assert not doc.can_undo

    def test_pivot_without_riders_is_refused(self, doc):
        with pytest.raises(ValueError, match='without its riding'):
            doc.assign_part([atom(doc, 'C10')], 1, 21.0)

    def test_riders_without_pivot_are_refused(self, doc):
        with pytest.raises(ValueError, match='rides on'):
            doc.assign_part([atom(doc, n) for n in ('H10A', 'H10B', 'H10C')], 1, 21.0)

    def test_non_contiguous_atoms_get_one_bracket_per_run(self, doc):
        doc.assign_part([atom(doc, 'N1'), atom(doc, 'C2')], 1, 10.5)
        assert sum(1 for line in doc.text.splitlines() if line.startswith('PART 1')) == 2
        shx = reparse(doc)
        assert shx.atoms.get_atom_by_name('N2').part.n == 0
        assert shx.atoms.get_atom_by_name('C1').part.n == 0

    def test_atom_already_in_a_part_is_refused(self):
        doc = ShelxDocument.from_file('tests/resources/p21c.res')
        with pytest.raises(ValueError, match='already in PART'):
            doc.assign_part([atom(doc, 'O1_4')], 1, 21.0)


# ---------------------------------------------------------- duplicates

class TestDuplicateAtoms:
    def test_copies_the_afix_group_in_a_new_part(self, doc):
        methyl = [atom(doc, name) for name in METHYL]
        names = ['C10B', 'H10D', 'H10E', 'H10F']
        shifted = [(a.xc + 1.0, a.yc, a.zc) for a in methyl]
        doc.add_free_variable(0.5)
        report = doc.duplicate_atoms(methyl, names, 2, -21.0, shifted)
        assert [copy.name for copy in report.added] == names
        assert lines_between(doc, 'C10B', 7) == [
            'C10B', 'AFIX', 'H10D', 'H10E', 'H10F', 'AFIX', 'PART']
        shx = reparse(doc)
        for name, xyz in zip(names, shifted):
            copy = shx.atoms.get_atom_by_name(name)
            assert copy.part.n == 2
            assert copy.sof == -21.0
            assert copy.cart_coords == pytest.approx(xyz, abs=1e-4)
        # Riding U (-1.5) now refers to the copied carbon.
        assert shx.atoms.get_atom_by_name('H10D').afix.mn == 137
        assert shx.atoms.get_atom_by_name('H10D').Uiso == pytest.approx(
            1.5 * shx.atoms.get_atom_by_name('C10B').Uiso)

    def test_isotropic_override(self, doc):
        c10 = atom(doc, 'C10')
        doc.duplicate_atoms([c10, *(atom(doc, n) for n in METHYL[1:])],
                            ['C10B', 'H10D', 'H10E', 'H10F'], 2, 11.0,
                            uvals=[[0.035], None, None, None])
        copy = reparse(doc).atoms.get_atom_by_name('C10B')
        assert copy.is_isotropic
        assert copy.uvals[0] == pytest.approx(0.035)

    def test_block_goes_after_the_sources_part_bracket(self, doc):
        methyl = [atom(doc, name) for name in METHYL]
        with doc.batch('split'):
            doc.add_free_variable(0.5)
            doc.assign_part(methyl, 1, 21.0)
            doc.duplicate_atoms(methyl, ['C10B', 'H10D', 'H10E', 'H10F'], 2, -21.0)
        words = lines_between(doc, 'PART', 17)
        assert words == ['PART', 'C10', 'AFIX', 'H10A', 'H10B', 'H10C', 'AFIX', 'PART',
                         'PART', 'C10B', 'AFIX', 'H10D', 'H10E', 'H10F', 'AFIX', 'PART', 'C11']

    def test_used_name_is_refused(self, doc):
        with pytest.raises(ValueError, match='already used'):
            doc.duplicate_atoms([atom(doc, 'N1')], ['C1'], 2, 11.0)

    def test_names_must_be_parallel(self, doc):
        with pytest.raises(ValueError, match='parallel'):
            doc.duplicate_atoms([atom(doc, 'N1')], [], 2, 11.0)


# -------------------------------------------------------------- naming

class TestNaming:
    def test_suffix(self, doc):
        assert doc.split_name(atom(doc, 'C10'), 'B') == 'C10B'

    def test_too_long_uses_the_next_letter_of_the_family(self, doc):
        assert doc.split_name(atom(doc, 'H10A'), 'B') == 'H10D'

    def test_keep_own_name(self, doc):
        assert doc.split_name(atom(doc, 'H10A'), 'A', keep_own_name=True) == 'H10A'

    def test_taken_names_are_respected(self, doc):
        assert doc.split_name(atom(doc, 'C10'), 'B', taken={'c10b'}) == 'C10A'

    def test_free_atom_name(self, doc):
        assert doc.free_atom_name('N') == 'N3'


# ---------------------------------------------------- U values and EQIV

def test_set_uvals_makes_an_atom_isotropic(doc):
    doc.set_uvals([(atom(doc, 'C1'), [0.035])])
    reread = reparse(doc).atoms.get_atom_by_name('C1')
    assert reread.is_isotropic
    assert reread.uvals[0] == pytest.approx(0.035)


def test_set_uvals_rejects_wrong_length(doc):
    with pytest.raises(ValueError):
        doc.set_uvals([(atom(doc, 'C1'), [0.1, 0.2])])


class TestNameForOperation:
    def test_identity_is_the_plain_name(self, doc):
        assert doc.name_for_operation(atom(doc, 'C1'), 'X, Y, Z') == 'C1'
        assert not doc.can_undo

    def test_creates_and_reuses_an_eqiv(self, doc):
        first = doc.name_for_operation(atom(doc, 'C1'), '-X+1, -Y+1, -Z')
        again = doc.name_for_operation(atom(doc, 'C2'), '-X+1, -Y+1, -Z')
        assert first == 'C1_$1'
        assert again == 'C2_$1'
        assert doc.history.undo_labels == ['Add EQIV']

    def test_lattice_translation_is_a_different_operation(self, doc):
        doc.name_for_operation(atom(doc, 'C1'), '-X+1, -Y+1, -Z')
        assert doc.name_for_operation(atom(doc, 'C1'), '-X, -Y+1, -Z') == 'C1_$2'


def test_per_atom_sof(doc):
    methyl = [atom(doc, name) for name in METHYL]
    doc.assign_part(methyl, 1, [10.5, 11.0, 11.0, 11.0])
    shx = reparse(doc)
    assert shx.atoms.get_atom_by_name('C10').sof == 10.5
    assert shx.atoms.get_atom_by_name('H10A').sof == 11.0
    doc.duplicate_atoms(methyl, ['C10B', 'H10D', 'H10E', 'H10F'], 2, [10.25, 10.5, 10.5, 10.5])
    assert reparse(doc).atoms.get_atom_by_name('C10B').sof == 10.25


def test_header_restraint_goes_before_the_first_atoms_brackets(doc):
    first = [atom(doc, 'N1')]
    doc.assign_part(first, 1, 10.5)  # N1 is the first atom: PART 1 now opens the atom list
    name = doc.name_for_operation(atom(doc, 'C1'), '-X+1, -Y+1, -Z')
    doc.add_restraint(f'DFIX 1.5 N1 {name}', header=True)
    words = lines_between(doc, 'EQIV', 4)
    assert words == ['EQIV', 'DFIX', 'PART', 'N1']


# ------------------------------------------------- symmetry references


def test_unparsable_symmetry_operation_is_refused(doc):
    """A garbage operation must raise, not mint an invalid EQIV card."""
    with pytest.raises(ValueError):
        doc.name_for_operation(atom(doc, 'C1'), 'foo, bar, baz')
    with pytest.raises(ValueError):
        doc.name_for_operation(atom(doc, 'C1'), 'x, x, x')
    assert not doc.shelxfile.eqiv


def test_the_identity_needs_no_eqiv(doc):
    assert doc.name_for_operation(atom(doc, 'C1'), 'X, Y, Z') == 'C1'
    assert not doc.shelxfile.eqiv


def test_operations_one_cell_apart_get_their_own_eqiv(doc):
    """A restraint one cell further along is a different restraint."""
    from shelxfile.edit.eqiv_factory import EqivFactory

    near = doc.name_for_operation(atom(doc, 'C1'), '-X+1, -Y+1, -Z')
    far = doc.name_for_operation(atom(doc, 'C1'), '-X, -Y+1, -Z')
    assert near != far
    assert len(doc.shelxfile.eqiv) == 2

    factory = EqivFactory(doc.shelxfile)
    exact = factory.find('-X, -Y+1, -Z', reduce_translation=False)
    assert f'C1_{exact.id}' == far
    # The default lookup still ignores lattice translations.
    assert factory.find('-X, -Y+1, -Z') is factory.find('-X+1, -Y+1, -Z')


def test_reusing_an_operation_mints_no_second_card(doc):
    first = doc.name_for_operation(atom(doc, 'C1'), '-X+1, -Y+1, -Z')
    again = doc.name_for_operation(atom(doc, 'C2'), '1-X, 1-Y, -Z')
    assert first.split('_')[1] == again.split('_')[1]
    assert len(doc.shelxfile.eqiv) == 1


# ------------------------------------------------ duplicate arguments


@pytest.mark.parametrize('kwargs', [
    {'coordinates': [[0.1, 0.1, 0.1]]},
    {'uvals': [[0.03]]},
])
def test_short_per_atom_arguments_are_refused(doc, kwargs):
    """Every per-atom argument is parallel to *atoms*, or it is a ValueError.

    These used to fall through to an ``IndexError`` deep in the loop, which
    callers that catch ``ValueError`` could not handle.
    """
    sources = [atom(doc, 'C10'), atom(doc, 'H10A')]
    with pytest.raises(ValueError):
        doc.duplicate_atoms(sources, ['C80', 'H80A'], 2, 1.0, **kwargs)


# --------------------------------------------------------- REM comments


def commented(doc: ShelxDocument, after: str, *comments: str) -> ShelxDocument:
    """Re-read *doc* with *comments* inserted after the line starting with *after*."""
    lines = doc.text.splitlines()
    index = next(i for i, line in enumerate(lines) if line.upper().startswith(after.upper()))
    lines[index + 1:index + 1] = list(comments)
    fresh = ShelxDocument.from_file(RES)
    fresh._restore(fresh._history.capture('\n'.join(lines) + '\n', None))
    return fresh


def test_a_comment_before_the_closing_afix_does_not_shrink_the_group(doc):
    """A REM is a comment; it cannot end an AFIX group.

    Files with an embedded ``REM <hkl>`` block really do put comments
    between the last riding atom and its ``AFIX 0``.
    """
    doc = commented(doc, 'H10C', 'REM <hkl>', 'REM jkd77.hkl', 'REM </hkl>')
    methyl = [atom(doc, name) for name in METHYL]

    doc.duplicate_atoms(methyl, ['C80', 'H80A', 'H80B', 'H80C'], 2, 1.0)
    doc.assign_part(methyl, 1, 0.5)

    shx = reparse(doc)
    for name in ('C80', 'H80A', 'H80B', 'H80C'):
        assert shx.atoms.get_atom_by_name(name) is not None, name
    # The copies went after the sources' AFIX 0, not inside the group.
    written = [line.split()[0] for line in doc.text.splitlines()
               if line.strip() and not line[0].isspace()]
    assert written.index('C80') > written.index('REM')


def test_comments_between_selected_atoms_keep_them_in_one_part(doc):
    doc = commented(doc, 'H10A', 'REM a note')
    methyl = [atom(doc, name) for name in METHYL]
    openers = doc.assign_part(methyl, 1, 0.5).added
    assert len([card for card in openers if getattr(card, 'n', 0) == 1]) == 1


def test_a_comment_above_the_first_atom_keeps_restraints_out_of_its_brackets(doc):
    """`header=True` must stay outside every PART/AFIX/RESI bracket.

    A REM between the opening bracket and the first atom is a comment and
    must not stop the search, or the restraint lands inside the bracket --
    where an unqualified atom name means something else.
    """
    lines = doc.text.splitlines()
    plain = doc.shelxfile._header_insert_position()
    opener = next(i for i, line in enumerate(lines)
                  if line.upper().startswith(doc.shelxfile._reslist[plain].name.upper() + ' '))
    lines[opener:opener] = ['PART 1', 'REM a note']
    fresh = ShelxDocument.from_file(RES)
    fresh._restore(fresh._history.capture('\n'.join(lines) + '\n', None))

    position = fresh.shelxfile._header_insert_position()
    assert type(fresh.shelxfile._reslist[position]).__name__ == 'PART'

    fresh.add_restraint('SADI C1 C2', header=True)
    written = [line.split()[0] for line in fresh.text.splitlines()
               if line.strip() and not line[0].isspace()]
    assert written.index('SADI') < written.index('PART')
