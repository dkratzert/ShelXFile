"""Round-trip-safe restraint serialisation (plan item A1).

``Restraint.__str__`` used to return ``textline``, a string frozen at
parse time.  Editing ``restraint.atoms`` therefore changed the model but
never the written file -- which would have made the whole
atom/restraint-linking feature a silent no-op.

The fix is deliberately conservative:

* an **untouched** card is echoed exactly as it was read, so loading and
  writing a file never reflows instructions nobody asked about;
* an **edited** card rebuilds its atom list, but reuses the original card
  name and numeric parameters verbatim.

Reusing the original parameter *strings* matters: ``_set_defs_values``
resolves ``DEFS``-derived defaults into ``self.s`` and friends, so
re-formatting from the parsed floats would silently freeze a value that
the file expects ``DEFS`` to keep supplying.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile

RES = 'tests/resources/p21c.res'


@pytest.fixture
def shx() -> Shelxfile:
    s = Shelxfile()
    s.read_file(RES)
    return s


def _restraint(shx: Shelxfile, kind: str):
    for restraint in shx.restraints:
        if type(restraint).__name__ == kind and restraint.atoms:
            return restraint
    raise AssertionError(f'fixture has no {kind} with atoms')


# ----------------------------------------------------------- untouched cards

def test_untouched_restraint_is_echoed_verbatim(shx: Shelxfile) -> None:
    for restraint in shx.restraints:
        assert str(restraint) == restraint.textline
        assert not restraint.atoms_were_edited


def test_parsing_and_writing_does_not_reflow_restraints(shx: Shelxfile) -> None:
    """No edit means no change anywhere in the restraint block."""
    before = [str(r) for r in shx.restraints]
    shx.dumps()
    assert [str(r) for r in shx.restraints] == before


# --------------------------------------------------------------- edited cards

def test_removing_an_atom_reaches_the_rendered_line(shx: Shelxfile) -> None:
    """The regression that made the original plan a no-op."""
    restraint = _restraint(shx, 'SADI')
    victim = restraint.atoms[1]

    restraint.atoms.remove(victim)

    assert restraint.atoms_were_edited
    rendered = str(restraint)
    assert rendered != restraint.textline
    assert rendered.split()[1:].count(victim) == restraint.atoms.count(victim)


def test_edited_restraint_reaches_dumps(shx: Shelxfile) -> None:
    restraint = _restraint(shx, 'SADI')
    restraint.atoms[:] = ['C1', 'C2']
    assert str(restraint) in shx.dumps()


def test_in_place_mutation_is_detected(shx: Shelxfile) -> None:
    """Users mutate the list directly; a setter-only hook would miss this."""
    restraint = _restraint(shx, 'SADI')
    restraint.atoms.append('C9')
    assert restraint.atoms_were_edited
    assert str(restraint).endswith('C9')


def test_restoring_the_original_atoms_restores_the_original_line(shx: Shelxfile) -> None:
    restraint = _restraint(shx, 'SADI')
    original = list(restraint.atoms)

    restraint.atoms.append('C9')
    assert str(restraint) != restraint.textline

    restraint.atoms[:] = original
    assert not restraint.atoms_were_edited
    assert str(restraint) == restraint.textline


# ------------------------------------------------- what regeneration preserves

def test_residue_suffix_survives_regeneration(shx: Shelxfile) -> None:
    """``SADI_CCF3`` must not degrade to ``SADI``."""
    restraint = _restraint(shx, 'SADI')
    assert '_' in restraint._prefix_tokens[0], 'fixture should use a residue class'
    head = restraint._prefix_tokens[0]

    restraint.atoms.append('C9')

    assert str(restraint).split()[0] == head


def test_numeric_parameters_are_reused_verbatim(shx: Shelxfile) -> None:
    """Not re-formatted from the parsed floats."""
    restraint = _restraint(shx, 'DFIX')
    original_numbers = restraint.textline.split()[1:len(restraint._prefix_tokens)]

    restraint.atoms[:] = ['O1', 'C2']

    assert str(restraint).split()[1:len(restraint._prefix_tokens)] == original_numbers


def test_parameterless_restraint_regenerates(shx: Shelxfile) -> None:
    restraint = _restraint(shx, 'SAME')
    restraint.atoms[:] = ['C1', 'C2']
    assert str(restraint).split()[0] == restraint._prefix_tokens[0]
    assert str(restraint).split()[-2:] == ['C1', 'C2']


def test_range_tokens_survive_regeneration(shx: Shelxfile) -> None:
    """``>`` is an atom token, not a parameter."""
    restraint = next(r for r in shx.restraints if '>' in r.atoms)
    restraint.atoms.append('C9')
    assert '>' in str(restraint).split()


# ------------------------------------------------------------ consistency

def test_iter_and_split_follow_the_rendered_line(shx: Shelxfile) -> None:
    restraint = _restraint(shx, 'SADI')
    restraint.atoms[:] = ['C1', 'C2']
    assert list(restraint) == str(restraint).split()
    assert restraint.split() == str(restraint).split()


def test_edited_restraint_survives_a_roundtrip(shx: Shelxfile, tmp_path) -> None:
    restraint = _restraint(shx, 'SADI')
    restraint.atoms[:] = ['C1', 'C2', 'C1', 'C4']
    expected = str(restraint)

    out = tmp_path / 'edited.res'
    shx.write_shelx_file(out)
    reread = Shelxfile()
    reread.read_file(out)

    assert any(str(r) == expected for r in reread.restraints)


def test_bare_restraint_without_atoms_is_untouched(shx: Shelxfile) -> None:
    """B2 - a card authored with no atoms must not be seen as edited."""
    bare = next(r for r in shx.restraints if not r.atoms)
    assert not bare.atoms_were_edited
    assert str(bare) == bare.textline
