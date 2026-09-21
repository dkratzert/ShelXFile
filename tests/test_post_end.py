"""Content after ``END`` is inert output (plan items B2a and D-10).

SHELXL stops reading at ``END``.  A bare ``SADI`` makes it emit the
derived restraints *after* that point, using ``ACTA TABS`` part suffixes
(``C1^a`` for PART 1, ``C1A^b`` for PART 2) which -- as the generated
comment in real files states -- "has not yet been implemented for input
into SHELXL".

So those lines must be:

* tagged :attr:`CardLifetime.POST_END_OUTPUT`,
* excluded from atom linking and from restraint validation, and
* **preserved verbatim** (D-10): no edit may rewrite, renumber or drop
  them, because they record a refinement that already happened.

Before this was implemented ShelXFile treated them as live restraints and
produced 434 spurious "unknown atom" warnings across the reference
corpus (6.8 % of all warnings), for 477 restraints in 33 files.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit.card_meta import CardLifetime

RES = 'tests/resources/post_end_sadi.res'


@pytest.fixture
def shx() -> Shelxfile:
    s = Shelxfile()
    s.read_file(RES)
    return s


def _post_end_lines(dump: str) -> list[str]:
    lines = dump.splitlines()
    for num, line in enumerate(lines):
        if line.strip().upper().startswith('END'):
            return lines[num:]
    return []


# ------------------------------------------------------------ classification

def test_restraints_before_end_are_live_input(shx: Shelxfile) -> None:
    live = [r for r in shx.restraints if r.lifetime is CardLifetime.INPUT]
    assert live, 'fixture must contain live restraints'
    assert any(type(r).__name__ == 'SAME' for r in live)


def test_restraints_after_end_are_tagged_as_output(shx: Shelxfile) -> None:
    post = [r for r in shx.restraints if r.lifetime is CardLifetime.POST_END_OUTPUT]
    assert len(post) == 2, 'both generated SADI cards must be tagged as output'
    assert all(type(r).__name__ == 'SADI' for r in post)


def test_card_lifetime_is_falsy_for_output() -> None:
    """``if card.lifetime:`` should read as "is this live input?"."""
    assert CardLifetime.INPUT
    assert not CardLifetime.POST_END_OUTPUT


# -------------------------------------------------------------- no false alarms

def test_caret_names_produce_no_unknown_atom_warnings(shx: Shelxfile) -> None:
    """The 6.8 %-of-all-warnings case."""
    noisy = [w for w in shx.restraint_errors if '^' in w]
    assert noisy == [], (
        'post-END output was validated against the atom list: '
        f'{noisy}'
    )


def test_live_restraints_are_still_validated() -> None:
    """Excluding post-END output must not disable validation generally."""
    shx = Shelxfile()
    shx.read_string(
        'TITL t\nCELL 0.71073 10 10 10 90 90 90\nZERR 4 0 0 0 0 0 0\n'
        'LATT 1\nSFAC C\nUNIT 4\n'
        'SADI C1 CZZ\n'
        'C1    1     0.1  0.2  0.3  11.00000  0.05\n'
        'HKLF 4\nEND\n'
    )
    assert any('CZZ' in w for w in shx.restraint_errors), (
        'a genuinely unknown atom before END must still be reported'
    )


# ------------------------------------------------------ D-10: verbatim output

def test_post_end_block_survives_atom_deletion_verbatim(shx: Shelxfile) -> None:
    before = _post_end_lines(shx.dumps())
    assert before, 'fixture must have a post-END block'

    shx.atoms.get_atom_by_name('C1_0').delete()

    assert _post_end_lines(shx.dumps()) == before, (
        'an atom deletion rewrote the post-END block'
    )


def test_post_end_block_survives_restraint_deletion_verbatim(shx: Shelxfile) -> None:
    before = _post_end_lines(shx.dumps())
    live = next(r for r in shx.restraints if r.lifetime is CardLifetime.INPUT
                and r.atoms)

    live.delete()

    assert _post_end_lines(shx.dumps()) == before


def test_post_end_caret_names_are_never_rewritten(shx: Shelxfile) -> None:
    """The ``^a``/``^b`` suffixes must be echoed exactly as read."""
    shx.atoms.get_atom_by_name('C3_0').delete()
    dump = shx.dumps()
    assert 'C1^a C2^a  C3^b C4^b' in dump
    assert 'C1^a C4^b  C2^a C3^b' in dump


def test_post_end_restraints_are_not_removed_by_atom_deletion(shx: Shelxfile) -> None:
    """Even though C3 is named (as ``C3^b``) in a post-END SADI."""
    before = len([r for r in shx.restraints
                  if r.lifetime is CardLifetime.POST_END_OUTPUT])

    for name in ('C1_0', 'C2_0', 'C3_0', 'C4_0'):
        shx.atoms.get_atom_by_name(name).delete()

    after = len([r for r in shx.restraints
                 if r.lifetime is CardLifetime.POST_END_OUTPUT])
    assert after == before


def test_roundtrip_preserves_the_post_end_block(shx: Shelxfile, tmp_path) -> None:
    out = tmp_path / 'rt.res'
    shx.write_shelx_file(out)
    reread = Shelxfile()
    reread.read_file(out)
    assert _post_end_lines(reread.dumps()) == _post_end_lines(shx.dumps())
