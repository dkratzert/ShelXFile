"""One minimal, hand-written sample per atom-referencing SHELXL card.

These are the **T2 synthetic fixtures** of the plan: written by hand rather
than copied from real data, so that every card gets the same scrutiny
whether or not it happens to be common.  Corpus frequency measures
prevalence, not importance.

Every sample and every expectation below comes from the SHELXL instruction
reference at <https://shelx.uni-goettingen.de/shelxl_html.php>; the quoted
phrase in each ``why`` is the sentence the expectation rests on.  They are
deliberately *not* derived from ``MIN_ATOMS``, which is the thing under
test.

The catalog is checked for completeness against the card registries in
``tests/test_card_coverage.py``: adding a new card class to ``shelxfile``
fails the suite until a sample lands here.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto

HEADER = (
    'TITL card coverage fixture\n'
    'CELL 0.71073 10.0 11.0 12.0 90 90 90\n'
    'ZERR 4 0.001 0.001 0.001 0 0 0\n'
    'LATT -1\n'
    'SYMM -X, 1/2+Y, 1/2-Z\n'
    'SFAC C H O\n'
    'UNIT 24 24 4\n'
    'WGHT 0.1\n'
    'FVAR 0.5\n'
)

#: A fixed pool of atoms every sample draws from, so a failing sample
#: differs from its neighbours only in the instruction under test.
ATOMS = (
    'O1    3     0.0100  0.0200  0.0300  11.00000  0.0500\n'
    'C1    1     0.1000  0.2000  0.3000  11.00000  0.0500\n'
    'C2    1     0.1500  0.2500  0.3500  11.00000  0.0500\n'
    'C3    1     0.2000  0.3000  0.4000  11.00000  0.0500\n'
    'C4    1     0.2500  0.3500  0.4500  11.00000  0.0500\n'
    'C5    1     0.3000  0.4000  0.5000  11.00000  0.0500\n'
    'C6    1     0.3500  0.4500  0.5500  11.00000  0.0500\n'
)

FOOTER = 'HKLF 4\nEND\n'


class Expect(Enum):
    """What should become of the sampled card when its atom is deleted."""

    #: The card survives, no longer naming the deleted atom.
    TRIMMED = auto()

    #: The card is removed: too little of it is left to mean anything.
    REMOVED = auto()

    #: The card is untouched -- it never named the atom in the first
    #: place, either because an empty list means *all atoms* or because
    #: the grammar has no atom list at all.
    UNTOUCHED = auto()


@dataclass(frozen=True)
class CardSample:
    """A card, and what deleting one of its atoms should do.

    :param keyword: the SHELXL instruction this sample covers.
    :param instruction: the instruction line, verbatim.
    :param delete: name of the atom to delete in the deletion test.
    :param expect: the intended fate of the card.
    :param why: the manual's justification.  Present so a failing test
        says what was meant, not merely what differed.
    """

    keyword: str
    instruction: str
    delete: str
    expect: Expect
    why: str

    def source(self) -> str:
        """The full ``.res`` file containing this sample."""
        return HEADER + self.instruction + '\n' + ATOMS + FOOTER


def _s(keyword, instruction, delete, expect, why) -> CardSample:
    return CardSample(keyword, instruction, delete, expect, why)


_SAMPLES: list[tuple[str, CardSample]] = [
    # ------------------------------------------------ distance restraints
    ('', _s('DFIX', 'DFIX 1.45 0.02 C1 C2 C3 C4', 'C1', Expect.TRIMMED,
            '"The distance between the first and second named atom, the '
            'third and fourth ..." -- one pair goes, one remains.')),
    ('_pair', _s('DFIX', 'DFIX 1.45 0.02 C1 C2', 'C1', Expect.REMOVED,
                 'A lone pair loses half of itself; no distance is left.')),
    ('', _s('DANG', 'DANG 2.40 0.04 C1 C3 C2 C4', 'C1', Expect.TRIMMED,
            '"interpreted in exactly the same way as DFIX".')),
    ('', _s('SADI', 'SADI 0.02 C1 C2 C3 C4', 'C1', Expect.TRIMMED,
            '"The distances between the first and second named atoms, the '
            'third and fourth ... are restrained to be equal."')),
    ('_pair', _s('SADI', 'SADI 0.02 C1 C2', 'C1', Expect.REMOVED,
                 'One pair cannot be restrained equal to nothing.')),
    ('', _s('CHIV', 'CHIV 0.0 0.1 C1 C2', 'C2', Expect.TRIMMED,
            '"The chiral volumes of the named atoms are restrained" -- '
            'each atom independently, so one suffices.')),
    ('', _s('FLAT', 'FLAT 0.1 C1 C2 C3 C4 C5', 'C5', Expect.TRIMMED,
            '"FLAT s[0.1] four or more atoms" -- four remain.')),
    ('_minimum', _s('FLAT', 'FLAT 0.1 C1 C2 C3 C4', 'C1', Expect.REMOVED,
                    'Three atoms are below the documented minimum of four.')),
    ('', _s('BUMP', 'BUMP 0.02', 'C1', Expect.UNTOUCHED,
            '"BUMP s[0.02]" -- the grammar has no atom list at all.')),
    ('', _s('DEFS', 'DEFS 0.02 0.1 0.01 0.04', 'C1', Expect.UNTOUCHED,
            '"DEFS sd sf su ss maxsof" -- default esds, names nothing.')),

    # -------------------------------------------- displacement restraints
    ('', _s('SIMU', 'SIMU 0.04 0.08 1.7 C1 C2 C3', 'C1', Expect.TRIMMED,
            '"Atoms closer than dmax are restrained ... to have the same '
            'Uij" -- a pairwise relation, so two must remain.')),
    ('_bare', _s('SIMU', 'SIMU 0.04', 'C1', Expect.UNTOUCHED,
                 '"If no atoms are given, all non-hydrogen atoms are '
                 'understood."')),
    ('', _s('DELU', 'DELU 0.01 0.01 C1 C2 C3', 'C1', Expect.TRIMMED,
            '"All bonds ... connecting atoms on the same DELU instruction '
            'are subject to a rigid bond restraint."')),
    ('', _s('RIGU', 'RIGU 0.004 0.004 C1 C2 C3', 'C1', Expect.TRIMMED,
            'Same free-list shape as DELU.')),
    ('', _s('ISOR', 'ISOR 0.1 C1 C2', 'C1', Expect.TRIMMED,
            '"The named atoms are restrained ... so that their Uij '
            'components approximate to isotropic behavior" -- per atom.')),
    ('_bare', _s('ISOR', 'ISOR 0.1', 'C1', Expect.UNTOUCHED,
                 '"If no atoms are given, all non-hydrogen atoms are '
                 'understood." An empty list is not a dead card (B2).')),
    ('', _s('EADP', 'EADP C1 C2 C3', 'C1', Expect.TRIMMED,
            '"The same ... displacement parameters are used for all the '
            'named atoms" -- two can still share.')),
    ('_pair', _s('EADP', 'EADP C1 C2', 'C1', Expect.REMOVED,
                 'One atom cannot share a parameter with itself.')),
    ('', _s('EXYZ', 'EXYZ C1 C2 C3', 'C1', Expect.TRIMMED,
            '"The same x, y and z parameters are used for all the named '
            'atoms."')),

    # --------------------------------------------------- similarity / NCS
    ('', _s('SAME', 'SAME 0.02 0.04 C4 C5 C6', 'C4', Expect.TRIMMED,
            '"The list of atoms ... is compared with the same number of '
            'atoms which follow the SAME instruction."')),
    ('', _s('NCSY', 'NCSY 1 0.1 0.05 C1 C2 C3', 'C1', Expect.TRIMMED,
            '"a SIMU restraint is generated ... for each pair of '
            'equivalent atoms" -- the mate lives in another residue, so a '
            'single named atom still means something.')),

    # ------------------------------------------------ non-restraint cards
    ('', _s('ANIS', 'ANIS C1 C2', 'C1', Expect.TRIMMED,
            '"The named atoms are made anisotropic" -- one at a time.')),
    ('_bare', _s('ANIS', 'ANIS', 'C1', Expect.UNTOUCHED,
                 '"ANIS on its own ... makes all following non-hydrogen '
                 'atoms anisotropic."')),
    ('_count', _s('ANIS', 'ANIS 2', 'C1', Expect.UNTOUCHED,
                  '"ANIS n: the next n isotropic non-hydrogen atoms are '
                  'made anisotropic" -- a count, not a reference.')),
    ('', _s('HFIX', 'HFIX 43 C1 C2', 'C1', Expect.TRIMMED,
            '"HFIX generates AFIX instructions and dummy hydrogen atoms '
            'bonded to the named atoms" -- per pivot.')),
    ('', _s('BOND', 'BOND C1 C2', 'C1', Expect.REMOVED,
            '"bond lengths for all bonds that involve two atoms '
            'referenced on the same BOND instruction."')),
    ('_bare', _s('BOND', 'BOND', 'C1', Expect.UNTOUCHED,
                 '"A BOND instruction with no parameters outputs bond '
                 'lengths ... for all bonds in the connectivity table."')),
    ('', _s('CONF', 'CONF C1 C2 C3 C4 C5', 'C1', Expect.TRIMMED,
            '"The named atoms define a chain of at least four atoms" -- '
            'four remain.')),
    ('_minimum', _s('CONF', 'CONF C1 C2 C3 C4', 'C1', Expect.REMOVED,
                    'Three atoms cannot define the documented chain of at '
                    'least four.')),
    ('_bare', _s('CONF', 'CONF', 'C1', Expect.UNTOUCHED,
                 '"If no atoms are specified, all possible torsion angles '
                 '... are generated from the connectivity array."')),
    ('', _s('CONN', 'CONN 4 C1 C2', 'C1', Expect.TRIMMED,
            '"the defaults may be overridden for the named atoms" -- a '
            'per-atom coordination limit.')),
    ('_bare', _s('CONN', 'CONN 0', 'C1', Expect.UNTOUCHED,
                 '"CONN without atom names changes the default value of '
                 'bmax for all following atoms."')),
    ('', _s('BLOC', 'BLOC 1 C1 C2', 'C1', Expect.TRIMMED,
            '"the x, y and z parameters of the named atoms are refined in '
            'cycle |n1|" -- per atom.')),
    ('_bare', _s('BLOC', 'BLOC 1', 'C1', Expect.UNTOUCHED,
                 '"A BLOC instruction with no atom names applies to all '
                 'atoms in the specified cycles."')),
    ('', _s('MPLA', 'MPLA 4 C1 C2 C3 C4', 'C4', Expect.TRIMMED,
            '"A least-squares plane is calculated through the first na of '
            'the named atoms ... na must be at least 3."')),
    ('_minimum', _s('MPLA', 'MPLA 3 C1 C2 C3', 'C3', Expect.REMOVED,
                    'Two atoms cannot define a plane.')),
    ('_no_na', _s('MPLA', 'MPLA C1 C2 C3 C4', 'C4', Expect.TRIMMED,
                  '"If na is omitted the plane is fitted to all the atoms '
                  'specified."')),
    ('', _s('BIND', 'BIND C1 C2', 'C1', Expect.REMOVED,
            '"BIND atom1 atom2" -- a bond needs both ends.')),
    ('_parts', _s('BIND', 'BIND 1 2', 'C1', Expect.UNTOUCHED,
                  '"BIND m n: atoms in PART m may bond to atoms in PART n" '
                  '-- a different instruction that names no atoms.')),
    ('', _s('FREE', 'FREE C1 C2', 'C1', Expect.REMOVED,
            '"FREE atom1 atom2" -- a bond needs both ends.')),
    ('', _s('HTAB', 'HTAB O1 C1', 'C1', Expect.REMOVED,
            '"The donor atom D and acceptor A should be specified."')),
    ('_bare', _s('HTAB', 'HTAB 2.0', 'C1', Expect.UNTOUCHED,
                 '"HTAB dh[2.0]" -- a search distance, naming nothing.')),
]

#: Keyed by ``KEYWORD`` plus a suffix when a card has several samples.
CARD_SAMPLES: dict[str, CardSample] = {
    sample.keyword + suffix: sample for suffix, sample in _SAMPLES
}

#: Keywords covered, for the completeness gate.
COVERED_KEYWORDS = frozenset(sample.keyword for _, sample in _SAMPLES)
