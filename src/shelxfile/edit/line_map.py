"""Mapping between rendered text lines and ``_reslist`` entries.

Lives in the edit layer, not the widget: turning the model into text and
back is SHELXL knowledge, and a view should only ask *"what is on line
N?"*.  See plan decision D-8.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

from shelxfile.misc.misc import wrap_line

if TYPE_CHECKING:
    from shelxfile import Shelxfile


class RenderedFile(NamedTuple):
    """The file as text, plus where each line came from.

    :param text: the rendered SHELX file, identical to
        :meth:`Shelxfile.dumps`.
    :param line_to_index: ``line_to_index[i]`` is the ``_reslist`` index
        that produced text line *i* (0-based).  One ``_reslist`` entry
        that wraps across several lines maps all of them to that index.
    """

    text: str
    line_to_index: list[int]


def render(shx: Shelxfile) -> RenderedFile:
    """Render *shx* to text while recording the origin of every line.

    Deliberately mirrors :meth:`Shelxfile.dumps` step for step: if the two
    ever diverge, cursor positions in an editor would point at the wrong
    instruction.
    """
    chunks: list[str] = []
    line_to_index: list[int] = []
    for num, item in enumerate(shx._reslist):
        if num in shx.delete_on_write:
            continue
        if item == '':
            continue
        wrapped = '\n'.join(wrap_line(x) for x in str(item).split('\n'))
        chunks.append(wrapped)
        line_to_index.extend([num] * (wrapped.count('\n') + 1))
    return RenderedFile('\n'.join(chunks), line_to_index)
