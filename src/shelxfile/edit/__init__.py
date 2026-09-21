"""Qt-free editing layer for ShelXFile.

Everything that understands SHELXL semantics and mutates the model lives
here.  Nothing under this package may import Qt: the GUI depends on the
edit layer, never the other way round.

See ``tests/test_layering.py`` for the guards that enforce this.
"""

from shelxfile.edit.card_meta import (
    AfixDependency,
    AtomGrouping,
    AtomListSemantics,
    CardLifetime,
)
from shelxfile.edit.document import ShelxDocument
from shelxfile.edit.graph import AtomRestraintGraph
from shelxfile.edit.line_map import RenderedFile, render
from shelxfile.edit.reports import (
    DeletionReport,
    EditReport,
    RemovalReason,
    RemovedItem,
)
from shelxfile.edit.same_links import SameFragments, SameResolver
from shelxfile.edit.token_resolver import (
    AtomReference,
    AtomTokenResolver,
    ResolvedAtoms,
)

__all__ = [
    'AfixDependency',
    'AtomGrouping',
    'AtomListSemantics',
    'AtomReference',
    'AtomRestraintGraph',
    'AtomTokenResolver',
    'CardLifetime',
    'DeletionReport',
    'EditReport',
    'RemovalReason',
    'RemovedItem',
    'RenderedFile',
    'ResolvedAtoms',
    'SameFragments',
    'SameResolver',
    'ShelxDocument',
    'render',
]
