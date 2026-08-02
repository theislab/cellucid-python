"""
Validation of the chunk listing a session declares.

The listing is checked against the profiles in :mod:`._schema`: every generic
chunk must carry exactly the metadata its profile declares, every dataset
dependent chunk id must be one this version knows, and the per-chunk
contributor and priority must be the ones the viewer writes. The dynamic id
families -- highlight cell sets, user-defined codes, analysis artifacts -- are
checked for exact URI-component encoding so a chunk id cannot smuggle a path.
"""

from __future__ import annotations

import json
import re
from typing import Any, cast
from urllib.parse import quote, unquote_to_bytes

from ._primitives import (
    _require_exact_bool,
    _require_exact_keys,
    _require_nonempty_string,
    _require_nonnegative_int,
)
from ._records import _CurrentChunkInventory
from ._schema import (
    _ANALYSIS_ARTIFACT_PREFIX,
    _CHUNK_META_KEYS,
    _CONTRIBUTOR_INDEX,
    _CURRENT_EXACT_CHUNK_CONTRIBUTORS,
    _CURRENT_GENERIC_STATIC_CHUNK_PROFILES,
    _HIGHLIGHT_CELLS_PREFIX,
    _JAVASCRIPT_URI_COMPONENT_SAFE,
    _STATIC_PROFILE_BY_ID,
    _STATIC_PROFILE_KEYS,
    _USER_DEFINED_CODES_PREFIX,
    _WIRE_ID_RE,
)


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
