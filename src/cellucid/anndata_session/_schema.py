"""
The shapes of the state-domain values a session is applied from.

A session file holds every state domain's value at every step
(:mod:`cellucid.session`); applying one to an AnnData reads four of them, as
``cellucid/assets/js/app/state/domains/`` captures them (``README.md``, and
each module's header): ``cells`` (which cells were open), ``field-overlays``
(renames, deletions and the fields a person or a notebook made),
``highlights`` (pages of highlight groups) and ``gene-sets`` (the notebook's
gene sets). Every key set is written once here and is closed on purpose: a
key this version does not know is a value it does not understand, and
applying it partly would put data on an AnnData that nothing describes.
"""

from __future__ import annotations

#: The state domains a session is applied from.
CELLS_DOMAIN = "cells"
FIELD_OVERLAYS_DOMAIN = "field-overlays"
HIGHLIGHTS_DOMAIN = "highlights"
GENE_SETS_DOMAIN = "gene-sets"
#: Presentation, recorded in the provenance (never applied to the values).
EXPRESSION_DOMAIN = "expression"
#: Presentation, read for the colors a user field's column is written with
#: (:mod:`._colors`): each view's colors, and which view is in focus.
COLORS_DOMAIN = "colors"
VIEWS_DOMAIN = "views"

#: ``expression``: ``{normalize, log1p}``, the expression values the genes
#: were shown in (``app/state/domains/expression.js``); both false: as stored.
EXPRESSION_KEYS = frozenset({"normalize", "log1p"})

#: ``cells``: ``{open: null | the subsample's specification (without nObs), genes}``.
CELLS_KEYS = frozenset({"open", "genes"})

#: ``field-overlays``.
OVERLAY_KEYS = frozenset({"renamedFields", "renamedCategories", "deletedFields", "userFields"})
RENAMED_FIELD_KEYS = frozenset({"source", "key", "name"})
RENAMED_CATEGORY_KEYS = frozenset({"source", "field", "category", "name"})
DELETED_FIELD_KEYS = frozenset({"source", "key", "purged"})
#: Every user field; a purged one, and a continuous copy of another field
#: (``sourceField`` set: its values are read from that field), hold only these.
USER_FIELD_KEYS = frozenset({
    "id", "source", "kind", "key", "deleted", "purged", "createdAt", "sourceField", "operation",
})
#: A categorical user field that is not purged.
USER_CATEGORY_KEYS = USER_FIELD_KEYS | {
    "categories", "codes", "sourcePages", "overlapStrategy", "overlapLabel", "intersectionLabels", "uncoveredLabel",
    "palette",
}
#: A continuous user field holding its own values (a column a notebook pushed), not purged.
USER_VALUES_KEYS = USER_FIELD_KEYS | {"values"}
#: A per-cell column: ``cells`` null (every cell, ``values[p]`` is original
#: position ``p``'s) or the ascending original positions it holds a value for.
COLUMN_KEYS = frozenset({"cells", "values"})
#: What a continuous copy's ``sourceField`` may name.
SOURCE_FIELD_KEYS = frozenset({"sourceKey", "sourceIndex", "kind"})
OVERLAP_STRATEGIES = frozenset({"first", "last", "overlap-label", "intersections"})
#: A categorical field's codes: each code type's largest value marks a cell
#: with no category.
CODE_TYPES = {"Uint8Array": 0xFF, "Uint16Array": 0xFFFF, "Uint32Array": 0xFFFF_FFFF}

#: ``highlights``.
HIGHLIGHTS_KEYS = frozenset({"active", "pages"})
PAGE_KEYS = frozenset({"id", "name", "color", "groups"})
GROUP_KEYS = frozenset({"id", "type", "label", "enabled", "cells"})
CATEGORY_GROUP_KEYS = GROUP_KEYS | {"field", "category"}
RANGE_GROUP_KEYS = GROUP_KEYS | {"field", "range"}
GROUP_TYPES = frozenset({"annotation", "category", "combined", "knn", "lasso", "proximity", "range"})
FIELD_REF_KEYS = frozenset({"source", "key"})
RANGE_KEYS = frozenset({"min", "max"})

#: ``gene-sets``.
GENE_SETS_KEYS = frozenset({"sets"})
GENE_SET_KEYS = frozenset({"name", "description", "genes"})
