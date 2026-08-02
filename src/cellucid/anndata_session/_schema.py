"""
The exact shapes a Cellucid session declares.

Every closed key set, chunk profile, and chunk-id prefix the session contract
defines is written once here, so the chunk validator, the highlight planner and
the field planner cannot drift from each other about what a session may
contain. A key set is closed on purpose: an unknown key is a session this
version does not understand, and applying it partially would put columns on an
AnnData that nothing describes.
"""

from __future__ import annotations

import re
from typing import cast

_WIRE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,179}$")
_FIELD_KEY_MAX_LENGTH = 256
_CATEGORY_LABEL_MAX_LENGTH = 256

_HIGHLIGHT_ROOT_KEYS = {"pages", "activePageId"}
_HIGHLIGHT_PAGE_KEYS = {"id", "name", "color", "highlightedGroups"}
_HIGHLIGHT_GROUP_REQUIRED_KEYS = {
    "id",
    "type",
    "label",
    "enabled",
    "cellCount",
}
_HIGHLIGHT_FIELD_GROUP_KEYS = _HIGHLIGHT_GROUP_REQUIRED_KEYS | {
    "fieldKey",
    "fieldIndex",
    "fieldSource",
}
_HIGHLIGHT_CATEGORY_GROUP_KEYS = _HIGHLIGHT_FIELD_GROUP_KEYS | {
    "categoryIndex",
    "categoryName",
}
_HIGHLIGHT_RANGE_GROUP_KEYS = _HIGHLIGHT_FIELD_GROUP_KEYS | {
    "rangeMin",
    "rangeMax",
}
_HIGHLIGHT_GROUP_TYPES = {
    "annotation",
    "category",
    "combined",
    "knn",
    "lasso",
    "proximity",
    "range",
}
_OVERLAY_ROOT_KEYS = {"renames", "deletedFields", "userDefinedFields"}
_FIELD_COMMON_KEYS = {
    "id",
    "source",
    "kind",
    "key",
    "isDeleted",
    "isPurged",
    "sourceField",
    "operation",
    "createdAt",
}
_CATEGORY_FIELD_KEYS = _FIELD_COMMON_KEYS | {
    "categories",
    "codesLength",
    "codesType",
    "centroidsByDim",
    "normalizedDims",
    "sourcePages",
    "overlapStrategy",
    "overlapLabel",
    "intersectionLabels",
    "uncoveredLabel",
}
_CONTINUOUS_FIELD_KEYS = _FIELD_COMMON_KEYS
_OVERLAP_STRATEGIES = {"first", "last", "overlap-label", "intersections"}
_CHUNK_META_KEYS = {
    "id",
    "contributorId",
    "priority",
    "kind",
    "codec",
    "label",
    "datasetDependent",
    "storedBytes",
    "uncompressedBytes",
}
_STATIC_PROFILE_KEYS = {
    "id",
    "contributorId",
    "priority",
    "kind",
    "codec",
    "label",
    "datasetDependent",
}
_CURRENT_GENERIC_STATIC_CHUNK_PROFILES = (
    {
        "id": "core/field-overlays",
        "contributorId": "field-overlays",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Field overlays",
        "datasetDependent": True,
    },
    {
        "id": "core/state",
        "contributorId": "core-state",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Core state",
        "datasetDependent": True,
    },
    {
        "id": "ui/dockable-layout",
        "contributorId": "dockable-layout",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Floating panels",
        "datasetDependent": False,
    },
    {
        "id": "analysis/windows",
        "contributorId": "analysis-windows",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Analysis windows",
        "datasetDependent": True,
    },
    {
        "id": "highlights/meta",
        "contributorId": "highlights-meta",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Highlight metadata",
        "datasetDependent": True,
    },
    {
        "id": "analysis/cache-inventory",
        "contributorId": "analysis-artifacts",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Analysis cache inventory",
        "datasetDependent": True,
    },
    {
        "id": "cinematic/camera",
        "contributorId": "cinematic-camera",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Cinematic camera path",
        "datasetDependent": True,
    },
)
_STATIC_PROFILE_BY_ID = {
    cast(str, profile["id"]): profile for profile in _CURRENT_GENERIC_STATIC_CHUNK_PROFILES
}
_CURRENT_EXACT_CHUNK_CONTRIBUTORS = {
    chunk_id: cast(str, profile["contributorId"])
    for chunk_id, profile in _STATIC_PROFILE_BY_ID.items()
}
_CONTRIBUTOR_ORDER = (
    "field-overlays",
    "user-defined-codes",
    "cinematic-camera",
    "core-state",
    "dockable-layout",
    "analysis-windows",
    "highlights-meta",
    "highlights-cells",
    "analysis-artifacts",
)
_CONTRIBUTOR_INDEX = {
    contributor_id: index for index, contributor_id in enumerate(_CONTRIBUTOR_ORDER)
}
_ANALYSIS_ARTIFACT_PREFIX = "analysis/artifacts/bulk-gene/"
_HIGHLIGHT_CELLS_PREFIX = "highlights/cells/"
_USER_DEFINED_CODES_PREFIX = "user-defined/codes/"
_JAVASCRIPT_URI_COMPONENT_SAFE = "-_.!~*'()"
