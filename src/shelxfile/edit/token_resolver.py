"""Expanding SHELXL atom-reference tokens into concrete atoms.

An instruction's atom list is not a plain list of names.  The SHELXL
manual allows, in any combination:

===================  ====================================================
``C1``               an atom in the current residue
``C1_2``             atom ``C1`` of residue 2
``C1_$3``            the image of ``C1`` under ``EQIV $3``
``C1_+`` / ``C1_-``  the same atom in the next / previous residue
``C1_*``             ``C1`` in every residue
``O1 > F9``          every non-hydrogen atom from ``O1`` to ``F9``
``C6 < C2``          the same, scanning backwards
``LAST``             the final atom in the file
``$C``               every atom of ``SFAC`` element C
===================  ====================================================

Two details are easy to get wrong and are handled explicitly here:

* Ranges may be written without spaces.  The manual's own examples
  include ``ISOR 0.1 O_201> LAST`` and ``MPLA 5 C11 >C15 Fe``.
* Residue scoping is applied **before** any symmetry operation:
  *"If the instruction codeword refers to a residue, this is applied to
  the named atoms before any symmetry operation specified with `_$n`"*.
  So ``RTAB_23 O..O OG_12 O_$3`` means ``(O_23)_$3``, not ``(O_0)_$3``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.atoms.atom import Atom

#: Trailing ``_$n`` symmetry-equivalent suffix, e.g. ``H9A_$1``.
EQIV_SUFFIX_RE = re.compile(r'^(.*)_\$(\d+)$')

#: Keyword standing for the last atom in the file.
LAST_KEYWORD = 'LAST'

RANGE_TOKENS = ('>', '<')


@dataclass(frozen=True)
class AtomReference:
    """One resolved atom reference.

    :param token: the source token, unchanged.
    :param fullname: ``NAME_RESINUM`` of the atom, or ``None`` when the
        token could not be resolved.
    :param eqiv_id: ``'$3'`` when the reference carried a symmetry
        suffix, else ``None``.
    """

    token: str
    fullname: str | None
    eqiv_id: str | None = None

    @property
    def resolved(self) -> bool:
        return self.fullname is not None


@dataclass
class ResolvedAtoms:
    """Outcome of expanding one instruction's atom list."""

    references: list[AtomReference] = field(default_factory=list)

    @property
    def fullnames(self) -> list[str]:
        """Resolved atom fullnames, in file order, without duplicates."""
        seen: set[str] = set()
        names: list[str] = []
        for ref in self.references:
            if ref.fullname and ref.fullname not in seen:
                seen.add(ref.fullname)
                names.append(ref.fullname)
        return names

    @property
    def eqiv_ids(self) -> list[str]:
        seen: set[str] = set()
        ids: list[str] = []
        for ref in self.references:
            if ref.eqiv_id and ref.eqiv_id not in seen:
                seen.add(ref.eqiv_id)
                ids.append(ref.eqiv_id)
        return ids

    @property
    def unresolved(self) -> list[str]:
        return [r.token for r in self.references if not r.resolved]

    def __len__(self) -> int:
        return len(self.references)


def split_range_tokens(tokens: list[str]) -> list[str]:
    """Separate ``>`` and ``<`` that are glued to a neighbouring name.

    ``['O_201>', 'LAST']`` becomes ``['O_201', '>', 'LAST']`` and
    ``['C11', '>C15']`` becomes ``['C11', '>', 'C15']``.  SHELXL accepts
    every spacing variant, and its own documentation uses several.
    """
    out: list[str] = []
    for token in tokens:
        if token in RANGE_TOKENS:
            out.append(token)
            continue
        rest = token
        while rest:
            for marker in RANGE_TOKENS:
                position = rest.find(marker)
                if position == -1:
                    continue
                if position > 0:
                    out.append(rest[:position])
                out.append(marker)
                rest = rest[position + 1:]
                break
            else:
                out.append(rest)
                rest = ''
    return [t for t in out if t]


class AtomTokenResolver:
    """Expands instruction atom tokens against a :class:`Shelxfile`."""

    def __init__(self, shx: Shelxfile) -> None:
        self.shx = shx

    # ------------------------------------------------------------ public

    def resolve(
        self,
        tokens: list[str],
        residue_numbers: list[int] | None = None,
    ) -> ResolvedAtoms:
        """Expand *tokens* into concrete atom references.

        :param tokens: the raw ``atoms`` list of an instruction.
        :param residue_numbers: residues the instruction applies to,
            from a ``_n`` / ``_class`` / ``_*`` suffix on the card name.
            ``None`` or empty means residue 0.
        """
        scopes = residue_numbers or [0]
        result = ResolvedAtoms()
        for scope in scopes:
            self._resolve_one_scope(split_range_tokens(tokens), scope, result)
        return result

    def resolve_for(self, card) -> ResolvedAtoms:
        """Expand the atom list of *card*, honouring its residue scope.

        Uses ``referenced_atoms`` so cards that keep their operands in
        named fields (``FREE``, ``HTAB``) work the same as those with a
        flat list.
        """
        numbers = getattr(card, 'residue_number', None)
        if isinstance(numbers, int):
            numbers = [numbers]
        tokens = getattr(card, 'referenced_atoms', None)
        if tokens is None:
            tokens = list(getattr(card, 'atoms', []))
        return self.resolve(list(tokens), numbers)

    # ----------------------------------------------------------- helpers

    def _resolve_one_scope(
        self,
        tokens: list[str],
        scope: int,
        result: ResolvedAtoms,
    ) -> None:
        index = 0
        while index < len(tokens):
            token = tokens[index]
            if token in RANGE_TOKENS:
                # A range needs the atom before it and the one after.
                previous = result.references[-1] if result.references else None
                following = tokens[index + 1] if index + 1 < len(tokens) else None
                if previous is not None and following is not None:
                    self._expand_range(previous, following, token, scope, result)
                    index += 2
                    continue
                index += 1
                continue
            result.references.extend(self._resolve_single(token, scope))
            index += 1

    def _resolve_single(self, token: str, scope: int) -> list[AtomReference]:
        if token.startswith('$'):
            return self._resolve_element(token)
        if token.upper() == LAST_KEYWORD:
            last = self._last_atom()
            return [AtomReference(token, last.fullname if last else None)]

        base, eqiv_id = self._split_eqiv(token)
        name, residue = self._split_residue(base, scope)

        if residue == '*':
            return self._resolve_wildcard(token, name, eqiv_id)

        fullname = f'{name}_{residue}'.upper()
        found = self.shx.atoms.get_atom_by_name(fullname)
        return [AtomReference(token, found.fullname if found else None, eqiv_id)]

    def _split_eqiv(self, token: str) -> tuple[str, str | None]:
        """Peel off a trailing ``_$n``.

        Done first so the remainder can be residue-scoped: the manual
        applies the residue before the symmetry operation.
        """
        match = EQIV_SUFFIX_RE.match(token)
        if not match:
            return token, None
        base, number = match.groups()
        return base, f'${number}'

    def _split_residue(self, base: str, scope: int) -> tuple[str, str | int]:
        """Split ``C1_2`` into ``('C1', 2)``, applying *scope* as default."""
        if '_' not in base:
            return base, scope
        name, _, suffix = base.rpartition('_')
        if not name:  # a token that merely starts with '_'
            return base, scope
        if suffix == '*':
            return name, '*'
        if suffix == '+':
            return name, scope + 1
        if suffix == '-':
            return name, scope - 1
        if suffix.isdigit():
            return name, int(suffix)
        # A residue *class* is not valid on an individual atom name, so
        # treat it as part of the name and let the lookup fail.
        return base, scope

    def _resolve_wildcard(
        self,
        token: str,
        name: str,
        eqiv_id: str | None,
    ) -> list[AtomReference]:
        found: list[AtomReference] = []
        for number in self.shx.residues.residue_numbers:
            atom = self.shx.atoms.get_atom_by_name(f'{name}_{number}'.upper())
            if atom is not None:
                found.append(AtomReference(token, atom.fullname, eqiv_id))
        if not found:
            found.append(AtomReference(token, None, eqiv_id))
        return found

    def _resolve_element(self, token: str) -> list[AtomReference]:
        """``$C`` means every atom of that ``SFAC`` element."""
        element = token[1:].upper()
        matches = [
            AtomReference(token, atom.fullname)
            for atom in self.shx.atoms
            if not atom.qpeak and self._element_of(atom).upper() == element
        ]
        return matches or [AtomReference(token, None)]

    def _element_of(self, atom: Atom) -> str:
        try:
            return self.shx.sfac2elem(atom.sfac_num) or ''
        except (AttributeError, IndexError, KeyError):
            return ''

    def _expand_range(
        self,
        start: AtomReference,
        end_token: str,
        marker: str,
        scope: int,
        result: ResolvedAtoms,
    ) -> None:
        """Append every atom between *start* and *end_token*.

        Hydrogens are skipped: the manual defines ``>`` as "all
        intervening **non-hydrogen** atoms".  The endpoints themselves are
        whatever the file says; only what lies between is filtered.
        """
        end_refs = self._resolve_single(end_token, scope)
        end = end_refs[0] if end_refs else None
        if start.fullname is None or end is None or end.fullname is None:
            result.references.extend(end_refs)
            return

        ordered = self.shx.atoms.all_atoms
        try:
            first = next(i for i, a in enumerate(ordered) if a.fullname == start.fullname)
            last = next(i for i, a in enumerate(ordered) if a.fullname == end.fullname)
        except StopIteration:
            result.references.extend(end_refs)
            return

        step = 1 if marker == '>' else -1
        if (step == 1 and last < first) or (step == -1 and last > first):
            # The file order contradicts the direction; emit the endpoint
            # only rather than inventing a span.
            result.references.extend(end_refs)
            return

        for position in range(first + step, last, step):
            atom = ordered[position]
            if atom.is_hydrogen or atom.qpeak:
                continue
            result.references.append(
                AtomReference(marker, atom.fullname, start.eqiv_id)
            )
        result.references.append(AtomReference(end_token, end.fullname, end.eqiv_id))

    def _last_atom(self) -> Atom | None:
        for atom in reversed(self.shx.atoms.all_atoms):
            if not atom.qpeak:
                return atom
        return None
