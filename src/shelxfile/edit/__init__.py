"""Qt-free editing layer for ShelXFile.

Everything that understands SHELXL semantics and mutates the model lives
here.  Nothing under this package may import Qt: the GUI depends on the
edit layer, never the other way round.

See ``tests/test_layering.py`` for the guards that enforce this.
"""

from shelxfile.edit.card_meta import CardLifetime
from shelxfile.edit.document import ShelxDocument
from shelxfile.edit.line_map import RenderedFile, render
from shelxfile.edit.reports import (
    DeletionReport,
    EditReport,
    RemovalReason,
    RemovedItem,
)

__all__ = [
    'CardLifetime',
    'DeletionReport',
    'EditReport',
    'RemovalReason',
    'RemovedItem',
    'RenderedFile',
    'ShelxDocument',
    'render',
]
