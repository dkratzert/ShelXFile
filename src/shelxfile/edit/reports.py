"""Structured results of an edit.

Reports are plain data: no Qt, no callbacks.  A view decides how to show
them; :class:`ShelxDocument` only states what happened and why.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from shelxfile.shelx.cards import Command, Restraint
    from shelxfile.atoms.atom import Atom


class RemovalReason(Enum):
    """Why the model dropped something.

    Every removal is attributed, so a caller can explain a cascade rather
    than presenting the user with silently vanished instructions.
    """

    #: The caller asked for this item directly.
    REQUESTED = auto()
    #: A card fell below the minimum atom count for its type.
    BELOW_MIN_ATOMS = auto()
    #: Keeping the card would have left it with no atoms, which for
    #: several instructions means "apply to all non-hydrogen atoms".
    WOULD_BECOME_GLOBAL = auto()
    #: A ``RESI``/``AFIX``/``PART`` bracket lost its last member.
    EMPTY_BRACKET = auto()
    #: An ``AFIX`` group dropped below its required member count.
    AFIX_UNDER_POPULATED = auto()
    #: The pivot atom of a riding/rotating ``AFIX`` group was deleted.
    AFIX_PIVOT_DELETED = auto()
    #: The positional counterpart on a ``SAME`` instruction went away.
    SAME_COUNTERPART = auto()
    #: The last ``_$n`` reference to an ``EQIV`` disappeared.
    EQIV_UNREFERENCED = auto()


@dataclass(frozen=True)
class RemovedItem:
    """One thing that was removed, and the reason for it."""

    item: Atom | Restraint | Command
    reason: RemovalReason
    #: Rendered form as it was before removal, for user-facing messages.
    text: str = ''

    def __str__(self) -> str:
        return f'{self.text or self.item} ({self.reason.name})'


@dataclass
class EditedCard:
    """A card that survived, with its atom list trimmed."""

    card: object
    before: str

    def __str__(self) -> str:
        return f'{self.before} -> {self.card}'


@dataclass
class DeletionReport:
    """Everything a deletion removed, transitively.

    A cascade may remove far more than was asked for, so the report covers
    the whole closure rather than only the direct request.
    """

    atoms: list[RemovedItem] = field(default_factory=list)
    cards: list[RemovedItem] = field(default_factory=list)
    #: Cards kept but rewritten, e.g. a restraint that lost one atom.
    edited: list[EditedCard] = field(default_factory=list)

    def add_atom(self, atom: Atom, reason: RemovalReason) -> None:
        self.atoms.append(RemovedItem(atom, reason, str(atom)))

    def add_card(self, card: Restraint | Command, reason: RemovalReason) -> None:
        self.cards.append(RemovedItem(card, reason, str(card)))

    def add_edit(self, card: Restraint | Command, before: str) -> None:
        self.edited.append(EditedCard(card, before))

    def extend(self, other: DeletionReport) -> None:
        self.atoms.extend(other.atoms)
        self.cards.extend(other.cards)
        self.edited.extend(other.edited)

    @property
    def is_empty(self) -> bool:
        return not self.atoms and not self.cards and not self.edited

    def __len__(self) -> int:
        return len(self.atoms) + len(self.cards) + len(self.edited)

    def summary(self) -> str:
        if self.is_empty:
            return 'nothing removed'
        parts = []
        if self.atoms:
            parts.append(f'{len(self.atoms)} atom(s)')
        if self.cards:
            parts.append(f'{len(self.cards)} card(s)')
        text = 'removed ' + ' and '.join(parts) if parts else ''
        if self.edited:
            trimmed = f'{len(self.edited)} card(s) trimmed'
            text = f'{text}, {trimmed}' if text else trimmed
        return text


@dataclass
class EditReport:
    """Result of an additive or modifying edit."""

    added: list[object] = field(default_factory=list)
    changed: list[object] = field(default_factory=list)
    removed: DeletionReport = field(default_factory=DeletionReport)

    @property
    def is_empty(self) -> bool:
        return not self.added and not self.changed and self.removed.is_empty
