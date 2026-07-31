"""Apply one exact Cellucid session contract to :class:`anndata.AnnData`."""

from __future__ import annotations

import copy
import json
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast
from urllib.parse import quote, unquote_to_bytes

import numpy as np
import pandas as pd
from scipy import sparse

from .session_bundle import CellucidSessionBundle
from .session_codecs import decode_delta_uvarint, decode_user_defined_codes

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


@dataclass(frozen=True)
class ApplySummary:
    """Columns materialized by one successful, atomic session application."""

    added_obs_columns: list[str]


@dataclass(frozen=True)
class _ColumnPlan:
    name: str
    values: pd.Series
    metadata: dict[str, Any]


@dataclass(frozen=True)
class _CurrentChunkInventory:
    ids: set[str]
    metadata_by_id: dict[str, dict[str, Any]]
    analysis_artifact_ids: tuple[str, ...]


def _require_exact_bool(value: Any, *, label: str) -> bool:
    if type(value) is not bool:
        raise TypeError(f"{label} must be exactly True or False")
    return value


def _require_nonempty_string(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise TypeError(f"{label} must be a non-empty string")
    return value


def _require_wire_id(value: Any, *, label: str) -> str:
    identifier = _require_nonempty_string(value, label=label)
    if _WIRE_ID_RE.fullmatch(identifier) is None:
        raise ValueError(
            f"{label} must use only ASCII letters, digits, '.', '_', or '-' "
            "and be at most 180 characters"
        )
    return identifier


def _require_nonnegative_int(value: Any, *, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer")
    if value < 0:
        raise ValueError(f"{label} must be non-negative")
    return int(value)


def _require_exact_keys(
    value: Any,
    expected: set[str],
    *,
    label: str,
    optional_keys: set[str] | None = None,
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise TypeError(f"{label} must be an object")
    allowed = expected | (optional_keys or set())
    if not expected <= set(value) or not set(value) <= allowed:
        missing = sorted(expected - set(value))
        unknown = sorted(set(value) - allowed)
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if unknown:
            details.append("unknown " + ", ".join(unknown))
        raise ValueError(f"{label} has invalid fields ({'; '.join(details)})")
    return value


def _expected_chunk_contributor(chunk_id: str) -> str | None:
    exact = _CURRENT_EXACT_CHUNK_CONTRIBUTORS.get(chunk_id)
    if exact is not None:
        return exact
    if chunk_id.startswith(_HIGHLIGHT_CELLS_PREFIX):
        group_id = chunk_id.removeprefix(_HIGHLIGHT_CELLS_PREFIX)
        return (
            "highlights-cells"
            if re.fullmatch(r"highlight_[1-9]\d*", group_id) is not None
            else None
        )
    if chunk_id.startswith(_USER_DEFINED_CODES_PREFIX):
        field_id = chunk_id.removeprefix(_USER_DEFINED_CODES_PREFIX)
        return "user-defined-codes" if _WIRE_ID_RE.fullmatch(field_id) is not None else None
    if chunk_id.startswith(_ANALYSIS_ARTIFACT_PREFIX):
        identity = chunk_id.removeprefix(_ANALYSIS_ARTIFACT_PREFIX).split("/")
        if len(identity) == 3 and all(_is_canonical_uri_component(segment) for segment in identity):
            return "analysis-artifacts"
    return None


def _is_canonical_uri_component(value: str) -> bool:
    if not value:
        return False
    try:
        decoded = unquote_to_bytes(value).decode("utf-8", errors="strict")
    except UnicodeDecodeError:
        return False
    return (
        bool(decoded)
        and quote(
            decoded,
            safe=_JAVASCRIPT_URI_COMPONENT_SAFE,
            encoding="utf-8",
            errors="strict",
        )
        == value
    )


def _validate_current_chunk_inventory(
    bundle: Any,
    raw_chunk_ids: list[str],
) -> _CurrentChunkInventory:
    manifest = _require_exact_keys(
        bundle.manifest,
        {"createdAt", "datasetFingerprint", "chunks"},
        label="bundle.manifest",
    )
    if manifest["datasetFingerprint"] != bundle.dataset_fingerprint:
        raise ValueError(
            "bundle.manifest.datasetFingerprint must exactly match bundle.dataset_fingerprint"
        )
    manifest_chunks = manifest.get("chunks")
    if not isinstance(manifest_chunks, list):
        raise TypeError("bundle.manifest.chunks must be an array")
    if len(manifest_chunks) != len(raw_chunk_ids):
        raise ValueError("bundle.manifest.chunks must exactly match bundle.list_chunk_ids()")

    metadata_by_id: dict[str, dict[str, Any]] = {}
    saw_lazy = False
    last_contributor_index = {"eager": -1, "lazy": -1}
    for index, raw_meta in enumerate(manifest_chunks):
        raw_meta = _require_exact_keys(
            raw_meta,
            _CHUNK_META_KEYS,
            label=f"bundle.manifest.chunks[{index}]",
        )
        chunk_id = _require_nonempty_string(
            raw_meta["id"],
            label=f"bundle.manifest.chunks[{index}].id",
        )
        if chunk_id != raw_chunk_ids[index]:
            raise ValueError("bundle.manifest.chunks must preserve the exact list_chunk_ids order")
        expected_contributor = _expected_chunk_contributor(chunk_id)
        if expected_contributor is None:
            raise ValueError(f"Unknown current session chunk {chunk_id!r}")
        contributor_id = _require_nonempty_string(
            raw_meta["contributorId"],
            label=f"bundle.manifest.chunks[{index}].contributorId",
        )
        if contributor_id != expected_contributor:
            raise ValueError(
                f"Session chunk {chunk_id!r} requires contributor "
                f"{expected_contributor!r}, received {contributor_id!r}"
            )
        priority = raw_meta["priority"]
        if priority not in {"eager", "lazy"}:
            raise ValueError(f"Session chunk {chunk_id!r} priority must be eager or lazy")
        if priority == "lazy":
            saw_lazy = True
        elif saw_lazy:
            raise ValueError("Session eager chunks must precede lazy chunks")
        contributor_index = _CONTRIBUTOR_INDEX[contributor_id]
        if contributor_index < last_contributor_index[priority]:
            raise ValueError(
                f"Session {priority} contributor groups must match their registered order"
            )
        last_contributor_index[priority] = contributor_index

        if raw_meta["kind"] not in {"json", "binary"}:
            raise ValueError(f"Session chunk {chunk_id!r} kind must be json or binary")
        if raw_meta["codec"] not in {"none", "gzip"}:
            raise ValueError(f"Session chunk {chunk_id!r} codec must be none or gzip")
        _require_nonempty_string(
            raw_meta["label"],
            label=f"Session chunk {chunk_id!r} label",
        )
        _require_exact_bool(
            raw_meta["datasetDependent"],
            label=f"Session chunk {chunk_id!r} datasetDependent",
        )
        stored_bytes = _require_nonnegative_int(
            raw_meta["storedBytes"],
            label=f"Session chunk {chunk_id!r} storedBytes",
        )
        uncompressed_bytes = _require_nonnegative_int(
            raw_meta["uncompressedBytes"],
            label=f"Session chunk {chunk_id!r} uncompressedBytes",
        )
        if raw_meta["codec"] == "none" and stored_bytes != uncompressed_bytes:
            raise ValueError(
                f"Session chunk {chunk_id!r} storedBytes and uncompressedBytes "
                "must match for codec 'none'"
            )

        static_profile = _STATIC_PROFILE_BY_ID.get(chunk_id)
        if static_profile is not None:
            for key in _STATIC_PROFILE_KEYS:
                if raw_meta[key] != static_profile[key]:
                    raise ValueError(
                        f"Session chunk {chunk_id!r} must match its exact current profile"
                    )
        elif contributor_id == "user-defined-codes":
            if (
                raw_meta["kind"] != "binary"
                or raw_meta["codec"] != "gzip"
                or raw_meta["datasetDependent"] is not True
            ):
                raise ValueError(
                    f"Session chunk {chunk_id!r} must be binary, gzip, and dataset-dependent"
                )
        elif contributor_id in {"highlights-cells", "analysis-artifacts"} and (
            raw_meta["priority"] != "lazy"
            or raw_meta["kind"] != "binary"
            or raw_meta["codec"] != "gzip"
            or raw_meta["datasetDependent"] is not True
        ):
            raise ValueError(
                f"Session chunk {chunk_id!r} must be lazy, binary, gzip, and dataset-dependent"
            )
        metadata_by_id[chunk_id] = raw_meta

    missing_singletons = [
        cast(str, profile["id"])
        for profile in _CURRENT_GENERIC_STATIC_CHUNK_PROFILES
        if profile["id"] not in metadata_by_id
    ]
    if missing_singletons:
        raise ValueError(
            "Current session requires singleton chunks: " + ", ".join(missing_singletons)
        )

    raw_inventory = bundle.decode_chunk("analysis/cache-inventory")
    inventory = _require_exact_keys(
        raw_inventory,
        {"artifactIds"},
        label="analysis/cache-inventory",
    )
    raw_artifact_ids = inventory["artifactIds"]
    if not isinstance(raw_artifact_ids, list):
        raise TypeError("analysis/cache-inventory.artifactIds must be an array")
    artifact_ids: list[str] = []
    seen_artifact_ids: set[str] = set()
    for index, raw_artifact_id in enumerate(raw_artifact_ids):
        artifact_id = _require_nonempty_string(
            raw_artifact_id,
            label=f"analysis/cache-inventory.artifactIds[{index}]",
        )
        if _expected_chunk_contributor(artifact_id) != "analysis-artifacts":
            raise ValueError(
                "analysis/cache-inventory artifact IDs must use the exact "
                "current bulk-gene identity"
            )
        if artifact_id in seen_artifact_ids:
            raise ValueError(f"analysis/cache-inventory duplicates artifact {artifact_id!r}")
        seen_artifact_ids.add(artifact_id)
        artifact_ids.append(artifact_id)

    manifest_artifact_ids = [
        chunk_id for chunk_id in raw_chunk_ids if chunk_id.startswith(_ANALYSIS_ARTIFACT_PREFIX)
    ]
    if manifest_artifact_ids != artifact_ids:
        raise ValueError("Analysis artifact chunks must exactly match the cache inventory order")
    canonical_inventory_bytes = len(
        json.dumps(
            {"artifactIds": artifact_ids},
            ensure_ascii=False,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    if metadata_by_id["analysis/cache-inventory"]["uncompressedBytes"] != canonical_inventory_bytes:
        raise ValueError(
            "analysis/cache-inventory uncompressedBytes must match its canonical payload"
        )

    return _CurrentChunkInventory(
        ids=set(raw_chunk_ids),
        metadata_by_id=metadata_by_id,
        analysis_artifact_ids=tuple(artifact_ids),
    )


def _require_field_key(value: Any, *, label: str) -> str:
    key = _require_nonempty_string(value, label=label)
    if key.strip() != key:
        raise ValueError(f"{label} cannot have leading or trailing whitespace")
    if len(key) > _FIELD_KEY_MAX_LENGTH:
        raise ValueError(f"{label} exceeds {_FIELD_KEY_MAX_LENGTH} Unicode code points")
    if ":" in key:
        raise ValueError(f"{label} cannot contain ':'")
    return key


def _validate_fingerprint(
    fingerprint: Any,
    adata: Any,
    *,
    expected_dataset_id: str,
) -> dict[str, Any]:
    # cellOrder ties the session to the cell ordering it was saved against. It
    # is accepted but not required, and never re-derived here: this applies a
    # session to an AnnData whose row order the caller controls, and the
    # viewer's digest is over exported coordinates an AnnData need not contain.
    # Enforcing that identity is the viewer's job; accepting both shapes is what
    # lets the viewer add the field without stranding existing session files.
    expected_keys = {"sourceType", "datasetId", "cellCount", "varCount"}
    fp = _require_exact_keys(
        fingerprint,
        expected_keys,
        label="session datasetFingerprint",
        optional_keys={"cellOrder"},
    )
    _require_nonempty_string(fp["sourceType"], label="datasetFingerprint.sourceType")
    dataset_id = _require_nonempty_string(
        fp["datasetId"],
        label="datasetFingerprint.datasetId",
    )
    cell_count = _require_nonnegative_int(
        fp["cellCount"],
        label="datasetFingerprint.cellCount",
    )
    var_count = _require_nonnegative_int(
        fp["varCount"],
        label="datasetFingerprint.varCount",
    )

    mismatches: list[str] = []
    if dataset_id != expected_dataset_id:
        mismatches.append(
            f"datasetId {dataset_id!r} != expected_dataset_id {expected_dataset_id!r}"
        )
    n_obs = getattr(adata, "n_obs", None)
    n_vars = getattr(adata, "n_vars", None)
    if cell_count != n_obs:
        mismatches.append(f"cellCount {cell_count} != adata.n_obs {n_obs}")
    if var_count != n_vars:
        mismatches.append(f"varCount {var_count} != adata.n_vars {n_vars}")
    if mismatches:
        raise ValueError("Dataset fingerprint mismatch: " + "; ".join(mismatches))
    return copy.deepcopy(fp)


def _reserve_column(existing: set[str], name: str) -> None:
    if name in existing:
        raise ValueError(f"Column already exists: {name}")
    existing.add(name)


def _require_optional_string(value: Any, *, label: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string or null")
    return value


def _require_category_list(
    value: Any,
    *,
    label: str,
    nonempty: bool,
) -> list[str | bool | int | float]:
    if not isinstance(value, list):
        raise TypeError(f"{label} must be an array")
    if nonempty and not value:
        raise ValueError(f"{label} must be a non-empty array")
    output: list[str | bool | int | float] = []
    identities: set[tuple[str, Any]] = set()
    for index, item in enumerate(value):
        exact = _require_category_primitive(item, label=f"{label}[{index}]")
        if isinstance(exact, str) and len(exact) > _CATEGORY_LABEL_MAX_LENGTH:
            raise ValueError(
                f"{label}[{index}] exceeds {_CATEGORY_LABEL_MAX_LENGTH} Unicode code points"
            )
        identity: tuple[str, Any]
        if isinstance(exact, bool):
            identity = ("boolean", exact)
        elif isinstance(exact, str):
            identity = ("string", exact)
        else:
            identity = ("number", float(exact))
        if identity in identities:
            raise ValueError(f"{label} must contain unique exact values")
        identities.add(identity)
        output.append(exact)
    return output


def _require_finite_number(value: Any, *, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"{label} must be a finite number")
    if not np.isfinite(value):
        raise ValueError(f"{label} must be a finite number")
    return float(value)


def _require_highlight_identity(value: Any, *, prefix: str, label: str) -> str:
    identity = _require_nonempty_string(value, label=label)
    if re.fullmatch(rf"{re.escape(prefix)}[1-9]\d*", identity) is None:
        raise ValueError(f"{label} must use the current {prefix}<positive integer> identity")
    return identity


def _require_category_primitive(value: Any, *, label: str) -> str | bool | int | float:
    if isinstance(value, str | bool):
        return value
    if isinstance(value, int | float) and not isinstance(value, bool):
        _require_finite_number(value, label=label)
        return value
    raise TypeError(f"{label} must be a string, finite number, or boolean")


def _validate_highlight_group_metadata(group: dict[str, Any], *, label: str) -> None:
    group_type = group.get("type")
    if group_type not in _HIGHLIGHT_GROUP_TYPES:
        raise ValueError(f"{label}.type is unsupported")
    if group_type == "category":
        expected_keys = _HIGHLIGHT_CATEGORY_GROUP_KEYS
    elif group_type == "range":
        expected_keys = _HIGHLIGHT_RANGE_GROUP_KEYS
    else:
        expected_keys = _HIGHLIGHT_GROUP_REQUIRED_KEYS
    _require_exact_keys(group, expected_keys, label=label)
    _require_highlight_identity(group["id"], prefix="highlight_", label=f"{label}.id")
    _require_nonempty_string(group["label"], label=f"{label}.label")
    _require_exact_bool(group["enabled"], label=f"{label}.enabled")
    cell_count = _require_nonnegative_int(group["cellCount"], label=f"{label}.cellCount")
    if cell_count < 1:
        raise ValueError(f"{label}.cellCount must be positive")
    if group_type in {"category", "range"}:
        _require_field_key(group["fieldKey"], label=f"{label}.fieldKey")
        _require_nonnegative_int(group["fieldIndex"], label=f"{label}.fieldIndex")
        if group["fieldSource"] not in {"obs", "var"}:
            raise ValueError(f"{label}.fieldSource must be 'obs' or 'var'")
    if group_type == "category":
        _require_nonnegative_int(
            group["categoryIndex"],
            label=f"{label}.categoryIndex",
        )
        _require_category_primitive(
            group["categoryName"],
            label=f"{label}.categoryName",
        )
    elif group_type == "range":
        range_min = _require_finite_number(group["rangeMin"], label=f"{label}.rangeMin")
        range_max = _require_finite_number(group["rangeMax"], label=f"{label}.rangeMax")
        if range_min > range_max:
            raise ValueError(f"{label}.rangeMin must not exceed rangeMax")


def _plan_highlights(
    bundle: Any,
    chunk_ids: set[str],
    *,
    n_obs: int,
    obs_index: pd.Index,
    prefix: str,
    existing_columns: set[str],
    materialize: bool,
    metadata_by_id: dict[str, dict[str, Any]],
) -> tuple[list[_ColumnPlan], dict[str, Any] | None, Any | None]:
    membership_prefix = "highlights/cells/"
    membership_ids = {chunk_id for chunk_id in chunk_ids if chunk_id.startswith(membership_prefix)}
    if "highlights/meta" not in chunk_ids:
        if membership_ids:
            raise ValueError("Session contains highlight membership chunks without highlights/meta")
        return [], None, None

    raw_meta = bundle.decode_chunk("highlights/meta")
    meta = _require_exact_keys(
        raw_meta,
        _HIGHLIGHT_ROOT_KEYS,
        label="highlights/meta",
    )
    pages = meta["pages"]
    if not isinstance(pages, list):
        raise TypeError("highlights/meta.pages must be an array")
    if not pages:
        raise ValueError("highlights/meta.pages must be a non-empty array")
    active_page_id = _require_highlight_identity(
        meta["activePageId"],
        prefix="page_",
        label="highlights/meta.activePageId",
    )

    plans: list[_ColumnPlan] = []
    referenced_memberships: set[str] = set()
    page_ids: set[str] = set()
    group_ids: set[str] = set()
    ordered_group_ids: list[str] = []
    highlights_uns: dict[str, Any] = {"groups": {}}

    for page_index, raw_page in enumerate(pages):
        page = _require_exact_keys(
            raw_page,
            _HIGHLIGHT_PAGE_KEYS,
            label=f"highlights/meta.pages[{page_index}]",
        )
        page_id = _require_highlight_identity(
            page["id"],
            prefix="page_",
            label=f"highlights/meta.pages[{page_index}].id",
        )
        if page_id in page_ids:
            raise ValueError(f"Duplicate highlight page id: {page_id}")
        page_ids.add(page_id)
        page_name = _require_nonempty_string(
            page["name"],
            label=f"highlights/meta.pages[{page_index}].name",
        )
        page_color = _require_nonempty_string(
            page["color"],
            label=f"highlights/meta.pages[{page_index}].color",
        )
        if re.fullmatch(r"#[0-9A-Fa-f]{6}", page_color) is None:
            raise ValueError(
                f"highlights/meta.pages[{page_index}].color must be a six-digit hex color"
            )
        groups = page["highlightedGroups"]
        if not isinstance(groups, list):
            raise TypeError(
                f"highlights/meta.pages[{page_index}].highlightedGroups must be an array"
            )

        for group_index, raw_group in enumerate(groups):
            if not isinstance(raw_group, dict):
                raise TypeError(
                    f"highlights/meta.pages[{page_index}].highlightedGroups"
                    f"[{group_index}] must be an object"
                )
            group_label = f"highlights/meta.pages[{page_index}].highlightedGroups[{group_index}]"
            _validate_highlight_group_metadata(raw_group, label=group_label)
            group_id = raw_group["id"]
            if group_id in group_ids:
                raise ValueError(f"Duplicate highlight group id: {group_id}")
            group_ids.add(group_id)
            ordered_group_ids.append(group_id)
            membership_chunk_id = f"{membership_prefix}{group_id}"
            cell_count = raw_group["cellCount"]

            if membership_chunk_id not in chunk_ids:
                raise ValueError(
                    f"Highlight group {group_id!r} is missing required chunk "
                    f"{membership_chunk_id!r}"
                )
            membership_meta = metadata_by_id[membership_chunk_id]
            if membership_meta["label"] != f"Highlight cells: {raw_group['label']}":
                raise ValueError(
                    f"{membership_chunk_id} label must match highlight group {group_id!r}"
                )
            raw_membership = bundle.decode_chunk(membership_chunk_id)
            if not isinstance(raw_membership, bytes | bytearray | memoryview):
                raise TypeError(f"{membership_chunk_id} must decode to binary bytes")
            if membership_meta["uncompressedBytes"] != len(raw_membership):
                raise ValueError(f"{membership_chunk_id} uncompressedBytes must match its payload")
            indices = decode_delta_uvarint(
                raw_membership,
                max_count=n_obs,
                max_index=n_obs - 1,
            )
            if len(indices) != cell_count:
                raise ValueError(
                    f"{membership_chunk_id} contains {len(indices)} indices, "
                    f"but metadata declares {cell_count}"
                )
            referenced_memberships.add(membership_chunk_id)

            column_name = f"{prefix}{group_id}"
            if materialize:
                _reserve_column(existing_columns, column_name)
                mask = np.zeros(n_obs, dtype=bool)
                mask[indices] = True
                plans.append(
                    _ColumnPlan(
                        name=column_name,
                        values=pd.Series(mask, index=obs_index, copy=False),
                        metadata={
                            "kind": "highlight",
                            "group_id": group_id,
                            "page_id": page_id,
                        },
                    )
                )
            highlights_uns["groups"][group_id] = {
                "obs_column": column_name if materialize else None,
                "page_id": page_id,
                "page_name": page_name,
                "group": copy.deepcopy(raw_group),
            }

    if active_page_id not in page_ids:
        raise ValueError(
            f"highlights/meta.activePageId {active_page_id!r} is not in the page inventory"
        )
    if membership_ids != referenced_memberships:
        unexpected = sorted(membership_ids - referenced_memberships)
        raise ValueError(
            "Session contains undeclared highlight membership chunks: " + ", ".join(unexpected)
        )
    manifest_membership_order = [
        chunk_id for chunk_id in metadata_by_id if chunk_id.startswith(_HIGHLIGHT_CELLS_PREFIX)
    ]
    expected_membership_order = [
        f"{_HIGHLIGHT_CELLS_PREFIX}{group_id}" for group_id in ordered_group_ids
    ]
    if manifest_membership_order != expected_membership_order:
        raise ValueError("Highlight membership chunks must match the exact metadata group order")
    return plans, highlights_uns, copy.deepcopy(meta)


def _validate_source_field(value: Any, *, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise TypeError(f"{label} must be an object")
    allowed = {"sourceKey", "sourceIndex", "kind"}
    if "sourceKey" not in value or set(value) - allowed:
        raise ValueError(f"{label} has invalid fields")
    source_key = _require_nonempty_string(
        value["sourceKey"],
        label=f"{label}.sourceKey",
    )
    output: dict[str, Any] = {"sourceKey": source_key}
    if "sourceIndex" in value:
        output["sourceIndex"] = _require_nonnegative_int(
            value["sourceIndex"],
            label=f"{label}.sourceIndex",
        )
    if "kind" in value:
        output["kind"] = _require_nonempty_string(
            value["kind"],
            label=f"{label}.kind",
        )
    return output


def _validate_operation(value: Any, *, label: str) -> None:
    if value is None:
        return
    if not isinstance(value, dict):
        raise TypeError(f"{label} must be an object or null")
    _require_nonempty_string(value.get("type"), label=f"{label}.type")


def _numeric_series(
    values: Any,
    *,
    n_obs: int,
    obs_index: pd.Index,
    label: str,
) -> pd.Series:
    array = np.asarray(values.toarray()) if sparse.issparse(values) else np.asarray(values)
    if array.ndim == 2 and array.shape[1] == 1:
        array = array[:, 0]
    if array.ndim != 1 or array.shape[0] != n_obs:
        raise ValueError(f"{label} must contain exactly {n_obs} cell values")
    if (
        not np.issubdtype(array.dtype, np.number)
        or np.issubdtype(array.dtype, np.bool_)
        or np.issubdtype(array.dtype, np.complexfloating)
    ):
        raise TypeError(f"{label} must contain real numeric values")
    if not np.isfinite(array).all():
        raise ValueError(f"{label} must contain only finite values")
    return pd.Series(array.copy(), index=obs_index)


def _materialize_continuous_alias(
    adata: Any,
    *,
    source: str,
    source_field: dict[str, Any],
    label: str,
) -> pd.Series:
    source_key = source_field["sourceKey"]
    n_obs = int(adata.n_obs)
    obs_index = adata.obs_names
    if source == "obs":
        matching = [index for index, key in enumerate(adata.obs.columns) if key == source_key]
        if len(matching) != 1:
            raise ValueError(f"{label} source obs field {source_key!r} must exist exactly once")
        source_index = matching[0]
        values = adata.obs.iloc[:, source_index].to_numpy(copy=True)
    else:
        matching = [index for index, key in enumerate(adata.var_names) if key == source_key]
        if len(matching) != 1:
            raise ValueError(f"{label} source gene {source_key!r} must exist exactly once")
        source_index = matching[0]
        if getattr(adata, "X", None) is None:
            raise ValueError(f"{label} cannot copy gene values because adata.X is absent")
        values = adata.X[:, source_index]

    declared_index = source_field.get("sourceIndex")
    if declared_index is not None and declared_index != source_index:
        raise ValueError(
            f"{label} sourceIndex {declared_index} does not match exact source "
            f"{source_key!r} at index {source_index}"
        )
    return _numeric_series(
        values,
        n_obs=n_obs,
        obs_index=obs_index,
        label=label,
    )


def _validate_category_metadata(
    field: dict[str, Any],
    *,
    label: str,
) -> list[str | bool | int | float]:
    categories = _require_category_list(
        field["categories"],
        label=f"{label}.categories",
        nonempty=True,
    )
    codes_length = _require_nonnegative_int(
        field["codesLength"],
        label=f"{label}.codesLength",
    )
    codes_type = field["codesType"]
    if codes_type not in {"Uint8Array", "Uint16Array"}:
        raise ValueError(f"{label}.codesType must be Uint8Array or Uint16Array")
    if not isinstance(field["centroidsByDim"], dict):
        raise TypeError(f"{label}.centroidsByDim must be an object")
    normalized_dims = field["normalizedDims"]
    if (
        not isinstance(normalized_dims, list)
        or any(
            isinstance(item, bool) or not isinstance(item, int) or item not in {1, 2, 3}
            for item in normalized_dims
        )
        or len(normalized_dims) != len(set(normalized_dims))
    ):
        raise ValueError(f"{label}.normalizedDims must be unique dimensions 1, 2, or 3")
    if not isinstance(field["sourcePages"], list):
        raise TypeError(f"{label}.sourcePages must be an array")
    if field["overlapStrategy"] not in _OVERLAP_STRATEGIES:
        raise ValueError(f"{label}.overlapStrategy is invalid")
    overlap_label = _require_optional_string(
        field["overlapLabel"],
        label=f"{label}.overlapLabel",
    )
    if overlap_label == "":
        raise ValueError(f"{label}.overlapLabel must be non-empty or null")
    if field["intersectionLabels"] is not None and not isinstance(
        field["intersectionLabels"], dict
    ):
        raise TypeError(f"{label}.intersectionLabels must be an object or null")
    uncovered_label = _require_optional_string(
        field["uncoveredLabel"],
        label=f"{label}.uncoveredLabel",
    )
    if uncovered_label == "":
        raise ValueError(f"{label}.uncoveredLabel must be non-empty or null")
    if field["sourceField"] is not None:
        _validate_source_field(field["sourceField"], label=f"{label}.sourceField")
    _validate_operation(field["operation"], label=f"{label}.operation")
    if codes_length > 0 and codes_type == "Uint8Array" and len(categories) > 255:
        raise ValueError(f"{label} has too many categories for Uint8Array codes")
    return categories


def _validate_field_created_at(value: Any, *, label: str) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"{label} must be a finite non-negative number or null")
    if not np.isfinite(value) or value < 0:
        raise ValueError(f"{label} must be a finite non-negative number or null")


def _validate_string_record(value: Any, *, label: str) -> None:
    if not isinstance(value, dict):
        raise TypeError(f"{label} must be an object")
    for key, item in value.items():
        _require_nonempty_string(key, label=f"{label} key")
        _require_nonempty_string(item, label=f"{label} value for {key!r}")


def _validate_unique_string_array(value: Any, *, label: str) -> list[str]:
    if not isinstance(value, list):
        raise TypeError(f"{label} must be an array")
    output: list[str] = []
    seen: set[str] = set()
    for index, item in enumerate(value):
        exact = _require_nonempty_string(item, label=f"{label}[{index}]")
        if exact in seen:
            raise ValueError(f"{label} duplicates {exact!r}")
        seen.add(exact)
        output.append(exact)
    return output


def _plan_user_defined_fields(
    bundle: Any,
    chunk_ids: set[str],
    adata: Any,
    *,
    prefix: str,
    existing_columns: set[str],
    materialize: bool,
    include_deleted: bool,
    metadata_by_id: dict[str, dict[str, Any]],
) -> tuple[list[_ColumnPlan], Any | None]:
    codes_prefix = "user-defined/codes/"
    code_chunk_ids = {chunk_id for chunk_id in chunk_ids if chunk_id.startswith(codes_prefix)}
    if "core/field-overlays" not in chunk_ids:
        if code_chunk_ids:
            raise ValueError(
                "Session contains user-defined code chunks without core/field-overlays"
            )
        return [], None

    raw_overlays = bundle.decode_chunk("core/field-overlays")
    overlays = _require_exact_keys(
        raw_overlays,
        _OVERLAY_ROOT_KEYS,
        label="core/field-overlays",
    )
    renames = _require_exact_keys(
        overlays["renames"],
        {"fields", "categories"},
        label="core/field-overlays.renames",
    )
    _validate_string_record(
        renames["fields"],
        label="core/field-overlays.renames.fields",
    )
    _validate_string_record(
        renames["categories"],
        label="core/field-overlays.renames.categories",
    )
    deleted_fields = _require_exact_keys(
        overlays["deletedFields"],
        {"deleted", "purged"},
        label="core/field-overlays.deletedFields",
    )
    deleted = set(
        _validate_unique_string_array(
            deleted_fields["deleted"],
            label="core/field-overlays.deletedFields.deleted",
        )
    )
    purged = set(
        _validate_unique_string_array(
            deleted_fields["purged"],
            label="core/field-overlays.deletedFields.purged",
        )
    )
    if not purged.issubset(deleted):
        raise ValueError("core/field-overlays.deletedFields.purged must be a subset of deleted")
    fields = overlays["userDefinedFields"]
    if not isinstance(fields, list):
        raise TypeError("core/field-overlays.userDefinedFields must be an array")

    plans: list[_ColumnPlan] = []
    field_ids: set[str] = set()
    referenced_code_chunks: set[str] = set()
    categorical_field_ids: list[str] = []
    n_obs = int(adata.n_obs)
    for index, raw_field in enumerate(fields):
        if not isinstance(raw_field, dict):
            raise TypeError(f"core/field-overlays.userDefinedFields[{index}] must be an object")
        label = f"core/field-overlays.userDefinedFields[{index}]"
        kind = raw_field.get("kind")
        expected_keys = _CATEGORY_FIELD_KEYS if kind == "category" else _CONTINUOUS_FIELD_KEYS
        field = _require_exact_keys(raw_field, expected_keys, label=label)
        if kind not in {"category", "continuous"}:
            raise ValueError(f"{label}.kind must be 'category' or 'continuous'")
        field_id = _require_wire_id(field["id"], label=f"{label}.id")
        if field_id in field_ids:
            raise ValueError(f"Duplicate user-defined field id: {field_id}")
        field_ids.add(field_id)
        source = field["source"]
        if source not in {"obs", "var"}:
            raise ValueError(f"{label}.source must be 'obs' or 'var'")
        field_key = _require_field_key(field["key"], label=f"{label}.key")
        is_deleted = _require_exact_bool(
            field["isDeleted"],
            label=f"{label}.isDeleted",
        )
        is_purged = _require_exact_bool(
            field["isPurged"],
            label=f"{label}.isPurged",
        )
        if is_purged and not is_deleted:
            raise ValueError(f"{label} purged fields must also be deleted")
        _validate_field_created_at(field["createdAt"], label=f"{label}.createdAt")
        should_materialize = materialize and (include_deleted or not is_deleted) and not is_purged

        column_name = f"{prefix}{field_key}"
        if kind == "continuous":
            source_field = _validate_source_field(
                field["sourceField"],
                label=f"{label}.sourceField",
            )
            _validate_operation(field["operation"], label=f"{label}.operation")
            if should_materialize:
                _reserve_column(existing_columns, column_name)
                values = _materialize_continuous_alias(
                    adata,
                    source=source,
                    source_field=source_field,
                    label=f"Continuous field {field_id!r}",
                )
                plans.append(
                    _ColumnPlan(
                        name=column_name,
                        values=values,
                        metadata={
                            "kind": "continuous",
                            "field_id": field_id,
                            "source": source,
                            "source_field": source_field,
                        },
                    )
                )
            continue

        categories = _validate_category_metadata(field, label=label)
        categorical_field_ids.append(field_id)
        codes_length = field["codesLength"]
        codes_type = field["codesType"]
        code_chunk_id = f"{codes_prefix}{field_id}"
        if codes_length != n_obs:
            raise ValueError(
                f"User-defined categorical field {field_id!r} has codesLength "
                f"{codes_length}, expected {n_obs} cell-aligned codes"
            )
        if n_obs == 0:
            if code_chunk_id in chunk_ids:
                raise ValueError(
                    f"Empty user-defined field {field_id!r} must not contain a codes chunk"
                )
            codes = np.empty(0, dtype=np.uint8 if codes_type == "Uint8Array" else np.uint16)
        else:
            if code_chunk_id not in chunk_ids:
                raise ValueError(
                    f"User-defined categorical field {field_id!r} is missing "
                    f"required chunk {code_chunk_id!r}"
                )
            code_meta = metadata_by_id[code_chunk_id]
            if code_meta["label"] != f"User-defined codes: {field_key}":
                raise ValueError(
                    f"{code_chunk_id} label must match user-defined field {field_id!r}"
                )
            raw_codes = bundle.decode_chunk(code_chunk_id)
            if not isinstance(raw_codes, bytes | bytearray | memoryview):
                raise TypeError(f"{code_chunk_id} must decode to binary bytes")
            if code_meta["uncompressedBytes"] != len(raw_codes):
                raise ValueError(f"{code_chunk_id} uncompressedBytes must match its payload")
            codes = decode_user_defined_codes(
                raw_codes,
                expected_length=codes_length,
                expected_codes_type=codes_type,
            )
            referenced_code_chunks.add(code_chunk_id)
            expected_dtype = np.dtype(np.uint8 if codes_type == "Uint8Array" else np.uint16)
            if codes.dtype != expected_dtype:
                raise ValueError(f"{code_chunk_id} dtype {codes.dtype} does not match {codes_type}")
            if len(codes) != codes_length:
                raise ValueError(
                    f"{code_chunk_id} contains {len(codes)} codes, expected {codes_length}"
                )

        category_codes = codes.astype(np.int64, copy=True)
        missing_sentinel = (
            (255 if codes_type == "Uint8Array" else 65_535)
            if field["uncoveredLabel"] is None
            else None
        )
        if missing_sentinel is not None:
            category_codes[category_codes == missing_sentinel] = -1
        invalid = np.flatnonzero((category_codes < -1) | (category_codes >= len(categories)))
        if invalid.size:
            position = int(invalid[0])
            code = int(category_codes[position])
            raise ValueError(
                f"User-defined categorical field {field_id!r} contains code "
                f"{code} at position {position}, but declares {len(categories)} categories"
            )

        if should_materialize:
            _reserve_column(existing_columns, column_name)
            categorical = pd.Categorical.from_codes(
                cast("Sequence[int]", category_codes),
                categories=pd.Index(categories),
                ordered=False,
            )
            plans.append(
                _ColumnPlan(
                    name=column_name,
                    values=pd.Series(categorical, index=adata.obs_names),
                    metadata={
                        "kind": "category",
                        "field_id": field_id,
                        "source": source,
                    },
                )
            )

    if code_chunk_ids != referenced_code_chunks:
        unexpected = sorted(code_chunk_ids - referenced_code_chunks)
        raise ValueError(
            "Session contains undeclared user-defined code chunks: " + ", ".join(unexpected)
        )
    # The AnnData bridge deliberately leaves core/state opaque. It proves
    # completeness and exact field-overlay order inside each priority bucket;
    # the web restore owner proves which active live/snapshot fields are eager.
    for priority in ("eager", "lazy"):
        actual_order = [
            chunk_id.removeprefix(_USER_DEFINED_CODES_PREFIX)
            for chunk_id, meta in metadata_by_id.items()
            if (chunk_id.startswith(_USER_DEFINED_CODES_PREFIX) and meta["priority"] == priority)
        ]
        expected_order = [
            field_id for field_id in categorical_field_ids if field_id in actual_order
        ]
        if actual_order != expected_order:
            raise ValueError(
                f"User-defined {priority} code chunks must match the exact "
                "field-overlay inventory order"
            )
    return plans, copy.deepcopy(overlays)


def _build_uns(
    adata: Any,
    *,
    bundle: Any,
    fingerprint: dict[str, Any],
    expected_dataset_id: str,
    highlight_meta: Any | None,
    highlights_uns: dict[str, Any] | None,
    overlays: Any | None,
    plans: list[_ColumnPlan],
) -> dict[str, Any]:
    raw_uns = getattr(adata, "uns", None)
    if not isinstance(raw_uns, dict):
        raise TypeError("adata.uns must be a dictionary")
    next_uns = copy.deepcopy(raw_uns)
    cellucid = next_uns.get("cellucid")
    if cellucid is None:
        cellucid = {}
        next_uns["cellucid"] = cellucid
    if not isinstance(cellucid, dict):
        raise TypeError("adata.uns['cellucid'] must be a dictionary if present")
    session = cellucid.get("session")
    if session is not None and not isinstance(session, dict):
        raise TypeError("adata.uns['cellucid']['session'] must be a dictionary if present")
    next_session: dict[str, Any] = {
        "manifest": copy.deepcopy(bundle.manifest),
        "dataset_fingerprint": fingerprint,
        "applied": {
            "expected_dataset_id": expected_dataset_id,
            "contract": "exact",
        },
        "materialized_fields": {plan.name: copy.deepcopy(plan.metadata) for plan in plans},
        "chunks": {},
    }
    if highlight_meta is not None:
        next_session["chunks"]["highlights/meta"] = highlight_meta
    if overlays is not None:
        next_session["chunks"]["core/field-overlays"] = overlays
    if highlights_uns is not None:
        next_session["highlights"] = highlights_uns
    cellucid["session"] = next_session
    return next_uns


def apply_cellucid_session_to_anndata(
    bundle: CellucidSessionBundle | str | Path,
    adata: Any,
    *,
    expected_dataset_id: str,
    inplace: bool = False,
    add_highlights: bool = True,
    highlights_prefix: str = "cellucid_highlight__",
    add_user_defined_fields: bool = True,
    user_defined_prefix: str = "",
    include_deleted_user_defined_fields: bool = False,
    store_uns: bool = True,
    return_summary: bool = False,
) -> Any | tuple[Any, ApplySummary]:
    """Apply a validated session atomically to the exact matching AnnData dataset."""
    expected_dataset_id = _require_nonempty_string(
        expected_dataset_id,
        label="expected_dataset_id",
    )
    _require_exact_bool(inplace, label="inplace")
    _require_exact_bool(add_highlights, label="add_highlights")
    _require_exact_bool(add_user_defined_fields, label="add_user_defined_fields")
    _require_exact_bool(
        include_deleted_user_defined_fields,
        label="include_deleted_user_defined_fields",
    )
    _require_exact_bool(store_uns, label="store_uns")
    _require_exact_bool(return_summary, label="return_summary")
    if not isinstance(highlights_prefix, str):
        raise TypeError("highlights_prefix must be a string")
    if not isinstance(user_defined_prefix, str):
        raise TypeError("user_defined_prefix must be a string")

    if isinstance(bundle, str | Path):
        bundle = CellucidSessionBundle(bundle)
    for attribute in (
        "dataset_fingerprint",
        "manifest",
        "list_chunk_ids",
        "decode_chunk",
    ):
        if not hasattr(bundle, attribute):
            raise TypeError(f"bundle must expose {attribute!r}")
    if not hasattr(adata, "obs") or not hasattr(adata, "var") or not hasattr(adata, "uns"):
        raise TypeError("adata must be an AnnData-compatible object")
    is_backed = getattr(adata, "isbacked", None)
    if type(is_backed) is not bool:
        raise TypeError("adata.isbacked must be exactly True or False")
    if is_backed and not inplace:
        raise ValueError(
            "A backed AnnData target cannot be applied with inplace=False; "
            "pass an explicitly materialized AnnData object or choose inplace=True"
        )

    fingerprint = _validate_fingerprint(
        bundle.dataset_fingerprint,
        adata,
        expected_dataset_id=expected_dataset_id,
    )
    raw_chunk_ids = bundle.list_chunk_ids()
    if not isinstance(raw_chunk_ids, list) or any(
        not isinstance(item, str) for item in raw_chunk_ids
    ):
        raise TypeError("bundle.list_chunk_ids() must return a list of strings")
    if len(raw_chunk_ids) != len(set(raw_chunk_ids)):
        raise ValueError("Session bundle contains duplicate chunk ids")
    inventory = _validate_current_chunk_inventory(bundle, raw_chunk_ids)
    chunk_ids = inventory.ids

    existing_columns = set(adata.obs.columns)
    highlight_plans, highlights_uns, highlight_meta = _plan_highlights(
        bundle,
        chunk_ids,
        n_obs=int(adata.n_obs),
        obs_index=adata.obs_names,
        prefix=highlights_prefix,
        existing_columns=existing_columns,
        materialize=add_highlights,
        metadata_by_id=inventory.metadata_by_id,
    )
    field_plans, overlays = _plan_user_defined_fields(
        bundle,
        chunk_ids,
        adata,
        prefix=user_defined_prefix,
        existing_columns=existing_columns,
        materialize=add_user_defined_fields,
        include_deleted=include_deleted_user_defined_fields,
        metadata_by_id=inventory.metadata_by_id,
    )
    plans = [*highlight_plans, *field_plans]

    next_obs = adata.obs.copy(deep=True)
    for plan in plans:
        next_obs[plan.name] = plan.values
    next_uns = (
        _build_uns(
            adata,
            bundle=bundle,
            fingerprint=fingerprint,
            expected_dataset_id=expected_dataset_id,
            highlight_meta=highlight_meta,
            highlights_uns=highlights_uns,
            overlays=overlays,
            plans=plans,
        )
        if store_uns
        else None
    )

    target = adata if inplace else adata.copy()
    target.obs = next_obs
    if next_uns is not None:
        target.uns = next_uns

    summary = ApplySummary(added_obs_columns=[plan.name for plan in plans])
    if return_summary:
        return target, summary
    return target
