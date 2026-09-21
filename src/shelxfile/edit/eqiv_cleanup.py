"""Keeping ``EQIV`` definitions and their ``_$n`` references consistent.

A token like ``C1_$3`` depends on two things: the atom ``C1`` and the
``EQIV $3`` that defines the symmetry operation.  Losing either leaves
the other dangling, so the dependency has to be maintained in both
directions.

Three rules from the manual constrain what may be done about it:

* *"Such a symmetry operation must be defined before it is used"* --
  an ``EQIV`` may be removed, but never moved.
* *"The same $n may not appear on two separate EQIV instructions"* --
  ids are unique, so gaps are never closed by renumbering; a surviving
  reference elsewhere would silently start pointing at a different
  operation.
* ``FREE``/``BIND``: *"Only one of the two atoms may be an equivalent
  atom"*, and ``HTAB``: *"Only the acceptor atom may specify a symmetry
  operation"*.

Post-``END`` definitions are left alone entirely.  SHELXL writes the
``HTAB``/``EQIV`` pairs it generates to the end of the file, and that
block is a record of a finished refinement.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from shelxfile.edit.card_meta import CardLifetime
from shelxfile.edit.reports import DeletionReport, RemovalReason
from shelxfile.edit.token_resolver import EQIV_SUFFIX_RE

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.edit.graph import AtomRestraintGraph
    from shelxfile.shelx.cards import EQIV

#: Cards where only one operand may carry a symmetry suffix, and which
#: one that is allowed to be. ``None`` means "either, but not both".
SYMMETRY_ARITY: dict[str, int | None] = {
    'FREE': None,   # "Only one of the two atoms may be an equivalent atom"
    'BIND': None,   # same wording
    'HTAB': 1,      # "Only the acceptor atom may specify a symmetry operation"
}


def symmetry_suffix(token: str) -> str | None:
    """``'$3'`` for ``'C1_$3'``, else ``None``."""
    match = EQIV_SUFFIX_RE.match(token)
    return f'${match.group(2)}' if match else None


def validate_symmetry_arity(card) -> str | None:
    """Check a card against its symmetry-operand limit.

    :returns: a description of the problem, or ``None`` when the card is
        acceptable.
    """
    limit = SYMMETRY_ARITY.get(type(card).__name__)
    if limit is None and type(card).__name__ not in SYMMETRY_ARITY:
        return None
    operands = card.referenced_atoms
    flagged = [i for i, token in enumerate(operands) if symmetry_suffix(token)]
    if len(flagged) < 2 and limit is None:
        return None
    if limit is None:
        return (f'{type(card).__name__} may carry a symmetry operation on only '
                f'one of its two atoms')
    if not flagged:
        return None
    if flagged != [limit]:
        which = 'acceptor' if limit == 1 else 'donor'
        return (f'{type(card).__name__} may only specify a symmetry operation '
                f'on the {which}')
    return None


class EqivCleaner:
    """Removes ``EQIV`` cards that nothing references any more."""

    def __init__(self, shx: Shelxfile, graph: AtomRestraintGraph) -> None:
        self.shx = shx
        self.graph = graph

    def orphans(self) -> list[EQIV]:
        """Live ``EQIV`` cards with no remaining ``_$n`` reference."""
        return self.graph.unreferenced_eqivs()

    def remove_orphans(self, report: DeletionReport | None = None) -> DeletionReport:
        """Drop every unreferenced ``EQIV``.

        Ids are never reassigned afterwards: the numbers that remain keep
        pointing at the operations they always did.
        """
        report = report if report is not None else DeletionReport()
        for card in self.orphans():
            if card.lifetime is not CardLifetime.INPUT:
                continue
            try:
                self.shx.remove_from_reslist(card)
            except ValueError:
                continue
            if card in self.shx.eqiv:
                self.shx.eqiv.remove(card)
            report.add_card(card, RemovalReason.EQIV_UNREFERENCED)
        return report

    def first_use_index(self, card: EQIV) -> int | None:
        """Position of the first card referencing *card*.

        An ``EQIV`` must stay above this; the manual requires it to be
        defined before it is used.
        """
        users = self.graph.cards_for_eqiv(card.id)
        positions = []
        for user in users:
            try:
                positions.append(self.shx.index_of(user))
            except ValueError:
                continue
        return min(positions) if positions else None

    def is_defined_before_use(self, card: EQIV) -> bool:
        """Whether *card* still precedes everything that references it."""
        first_use = self.first_use_index(card)
        if first_use is None:
            return True
        try:
            return self.shx.index_of(card) < first_use
        except ValueError:
            return True
