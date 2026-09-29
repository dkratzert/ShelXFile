"""Snapshot-based undo/redo for :class:`~shelxfile.edit.ShelxDocument`.

Every edit step stores the complete file text as it was *before* the step,
rather than an inverse operation per kind of edit.  That keeps undo correct
for everything the edit layer can do -- cascading deletions, renames that
retarget restraints, ``PART`` bracketing, duplicated ``AFIX`` groups --
without a hand-written inverse for each, and it cannot drift out of step as
new operations are added.

A SHELX file is small text and compresses well, so the history keeps every
step (no depth limit) as zlib-compressed bytes.

A caller that keeps its own state alongside the document (a viewer that
remembers which atoms it split, say) registers a *state provider* on the
document.  Its return value is stored with every snapshot and handed back
by :meth:`EditHistory.undo` / :meth:`EditHistory.redo`, so the caller can
restore its state to match the restored file.

Nothing here imports Qt.
"""

from __future__ import annotations

import zlib
from dataclasses import dataclass, field
from typing import Any

__all__ = ['EditHistory', 'HistoryState', 'UndoResult']


@dataclass(frozen=True)
class HistoryState:
    """One restorable state: the compressed file text and a caller payload.

    :ivar state_id: Identifies the state, so the history can tell whether
        the current state is the one that was last saved.
    """

    compressed_text: bytes
    payload: Any
    state_id: int

    @classmethod
    def capture(cls, text: str, payload: Any, state_id: int) -> HistoryState:
        return cls(zlib.compress(text.encode('utf-8')), payload, state_id)

    @property
    def text(self) -> str:
        return zlib.decompress(self.compressed_text).decode('utf-8')


@dataclass(frozen=True)
class UndoResult:
    """What an undo or redo restored.

    :ivar label: The label of the step that was undone or redone.
    :ivar payload: The caller payload stored with the restored state
        (``None`` when no state provider was registered).
    """

    label: str
    payload: Any


@dataclass
class EditHistory:
    """The undo and redo stacks of one document.

    Each undo entry is ``(label, state before the step)``; each redo entry
    is ``(label, state after the step)``.  A new step clears the redo
    stack, as in every editor.
    """

    _undo: list[tuple[str, HistoryState]] = field(default_factory=list)
    _redo: list[tuple[str, HistoryState]] = field(default_factory=list)
    _current_id: int = 0
    _next_id: int = 1
    _saved_id: int = 0

    # ----------------------------------------------------------- queries

    @property
    def can_undo(self) -> bool:
        return bool(self._undo)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo)

    @property
    def undo_label(self) -> str | None:
        return self._undo[-1][0] if self._undo else None

    @property
    def redo_label(self) -> str | None:
        return self._redo[-1][0] if self._redo else None

    @property
    def undo_labels(self) -> list[str]:
        """Labels of all undoable steps, most recent last."""
        return [label for label, _ in self._undo]

    @property
    def is_modified(self) -> bool:
        """Whether the current state differs from the last saved one."""
        return self._current_id != self._saved_id

    # ---------------------------------------------------------- mutation

    def capture(self, text: str, payload: Any) -> HistoryState:
        """Snapshot the current state (not yet recorded anywhere)."""
        return HistoryState.capture(text, payload, self._current_id)

    def push(self, label: str, before: HistoryState) -> None:
        """Record a finished step whose prior state was *before*."""
        self._undo.append((label, before))
        self._redo.clear()
        self._current_id = self._next_id
        self._next_id += 1

    def pop_undo(self, now: HistoryState) -> tuple[str, HistoryState]:
        """Take the step to undo; *now* becomes its redo entry."""
        label, before = self._undo.pop()
        self._redo.append((label, now))
        self._current_id = before.state_id
        return label, before

    def pop_redo(self, now: HistoryState) -> tuple[str, HistoryState]:
        """Take the step to redo; *now* becomes its undo entry again."""
        label, after = self._redo.pop()
        self._undo.append((label, now))
        self._current_id = after.state_id
        return label, after

    def mark_saved(self) -> None:
        """Declare the current state to be the one on disk."""
        self._saved_id = self._current_id
