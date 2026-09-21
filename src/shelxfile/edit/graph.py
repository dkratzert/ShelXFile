"""Which cards reference which atoms.

Derived data, rebuilt from the model rather than maintained alongside it.
That is a deliberate trade: a stale graph is worse than no graph, so
instead of hooking every mutation the graph notes the model version it
was built at and rebuilds itself when that moves on.

Three maps, because atom references run in more directions than the
original plan allowed for:

``atom -> cards``
    the obvious one.
``card -> atoms``
    includes atoms a card never names.  A ``SAME`` is compared against
    the run of atoms that follows it, so those belong to it too (see
    :mod:`shelxfile.edit.same_links`).
``eqiv -> cards``
    a ``C1_$3`` reference depends on ``EQIV $3`` as well as on ``C1``.
    Without this, removing the last such reference would leave an orphan
    ``EQIV`` behind.

Cards whose empty atom list means "everything", and inert post-``END``
output, are excluded: see :class:`AtomListSemantics` and
:class:`CardLifetime`.
"""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, Iterator

from shelxfile.edit.card_meta import (
    AtomGrouping,
    AtomReferencingCard,
    CardLifetime,
)
from shelxfile.edit.same_links import SameResolver
from shelxfile.edit.token_resolver import AtomTokenResolver

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.atoms.atom import Atom
    from shelxfile.shelx.cards import EQIV


class AtomRestraintGraph:
    """Links atoms to the cards that reference them.

    :param shx: the model to index.  Nothing is copied; the graph holds
        references into it and rebuilds when the model changes.
    """

    def __init__(self, shx: Shelxfile) -> None:
        self.shx = shx
        self._tokens = AtomTokenResolver(shx)
        self._same = SameResolver(shx)
        self._atom_to_cards: dict[str, list[AtomReferencingCard]] = defaultdict(list)
        self._card_to_atoms: dict[int, list[str]] = {}
        self._eqiv_to_cards: dict[str, list[AtomReferencingCard]] = defaultdict(list)
        self._cards: list[AtomReferencingCard] = []
        self._built_at: int | None = None

    # ------------------------------------------------------------ rebuild

    @property
    def is_stale(self) -> bool:
        return self._built_at != self.shx.model_version

    def ensure_current(self) -> None:
        """Rebuild if the model moved on since the last build."""
        if self.is_stale:
            self.rebuild()

    def rebuild(self) -> None:
        """Re-scan the model from scratch."""
        self._atom_to_cards = defaultdict(list)
        self._card_to_atoms = {}
        self._eqiv_to_cards = defaultdict(list)
        self._cards = []
        for card in self._linked_cards():
            self._cards.append(card)
            names = self._atoms_of(card)
            self._card_to_atoms[id(card)] = names
            for name in names:
                # Keyed upper-case throughout: atom names keep the file's
                # own casing (``H9a_1``), so a case-sensitive key would
                # make the two directions of the map disagree.
                self._atom_to_cards[name.upper()].append(card)
            for eqiv_id in self._eqiv_ids_of(card):
                self._eqiv_to_cards[eqiv_id].append(card)
        self._built_at = self.shx.model_version

    def _linked_cards(self) -> Iterator[AtomReferencingCard]:
        for item in self.shx._reslist:
            if isinstance(item, AtomReferencingCard) and item.is_atom_linked:
                yield item

    def _atoms_of(self, card: AtomReferencingCard) -> list[str]:
        if card.atom_grouping in (AtomGrouping.SAME_FOLLOWING_ATOMS,
                                  AtomGrouping.SAME_RESIDUE_CLASS):
            return [a.fullname for a in self._same.fragments(card).all_atoms]
        return self._resolve(card).fullnames

    def _eqiv_ids_of(self, card: AtomReferencingCard) -> list[str]:
        return self._resolve(card).eqiv_ids

    def _resolve(self, card: AtomReferencingCard):
        scope = getattr(card, 'residue_number', None)
        if isinstance(scope, int):
            scope = [scope]
        return self._tokens.resolve(card.referenced_atoms, scope)

    # ------------------------------------------------------------ queries

    @property
    def cards(self) -> list[AtomReferencingCard]:
        """Every card currently taking part in atom bookkeeping."""
        self.ensure_current()
        return list(self._cards)

    def cards_for_atom(self, atom: Atom | str) -> list[AtomReferencingCard]:
        """Cards affected if *atom* goes away.

        Includes cards that never name it: a ``SAME`` reaches the atoms
        following it in the file.
        """
        self.ensure_current()
        name = atom if isinstance(atom, str) else atom.fullname
        return list(self._atom_to_cards.get(name.upper(), ()))

    def atoms_for_card(self, card: AtomReferencingCard) -> list[str]:
        """Atom fullnames *card* depends on, expanded."""
        self.ensure_current()
        return list(self._card_to_atoms.get(id(card), ()))

    def cards_for_eqiv(self, eqiv_id: str) -> list[AtomReferencingCard]:
        """Cards referencing the given ``$n`` symmetry operation."""
        self.ensure_current()
        return list(self._eqiv_to_cards.get(eqiv_id, ()))

    def unreferenced_eqivs(self) -> list[EQIV]:
        """Live ``EQIV`` cards nothing points at any more.

        Candidates for removal once a cascade has dropped the last
        ``_$n`` reference to them.

        Post-``END`` definitions are excluded.  SHELXL writes the
        ``HTAB``/``EQIV`` pairs it generates *"to the end of the .res
        file"*, and that block is inert output which must be preserved
        verbatim -- its ``EQIV`` cards look unreferenced only because the
        ``HTAB`` cards using them are inert too.
        """
        self.ensure_current()
        return [card for card in self.shx.eqiv
                if card.id
                and card.lifetime is CardLifetime.INPUT
                and not self._eqiv_to_cards.get(card.id)]

    def referencing_atoms(self) -> list[str]:
        """Every atom fullname some card depends on, upper-cased."""
        self.ensure_current()
        return list(self._atom_to_cards)

    def __len__(self) -> int:
        self.ensure_current()
        return len(self._cards)

    def __repr__(self) -> str:
        self.ensure_current()
        return (f'<AtomRestraintGraph {len(self._cards)} cards, '
                f'{len(self._atom_to_cards)} atoms, '
                f'{len(self._eqiv_to_cards)} eqivs>')
