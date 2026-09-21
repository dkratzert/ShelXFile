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
    AtomReferencingCard,
    CardLifetime,
)
from shelxfile.edit.cascade import CascadeEngine, CascadePlan
from shelxfile.edit.document import ParseAttempt, ShelxDocument, restraint_keywords
from shelxfile.edit.eqiv_cleanup import EqivCleaner, validate_symmetry_arity
from shelxfile.edit.eqiv_factory import EqivFactory, canonical_symmop
from shelxfile.edit.graph import AtomRestraintGraph
from shelxfile.edit.line_map import RenderedFile, render
from shelxfile.edit.reports import (
    DeletionReport,
    EditReport,
    EditedCard,
    RemovalReason,
    RemovedItem,
    RenameReport,
    SkippedReference,
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
    'AtomReferencingCard',
    'AtomRestraintGraph',
    'AtomTokenResolver',
    'CardLifetime',
    'CascadeEngine',
    'CascadePlan',
    'DeletionReport',
    'EditReport',
    'EditedCard',
    'EqivCleaner',
    'EqivFactory',
    'ParseAttempt',
    'RemovalReason',
    'RemovedItem',
    'RenameReport',
    'RenderedFile',
    'ResolvedAtoms',
    'SameFragments',
    'SameResolver',
    'ShelxDocument',
    'SkippedReference',
    'canonical_symmop',
    'render',
    'restraint_keywords',
    'validate_symmetry_arity',
]
