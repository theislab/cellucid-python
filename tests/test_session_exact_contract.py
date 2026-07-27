from __future__ import annotations

import copy
import gzip
import inspect
import json
import struct
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid.anndata_session import apply_cellucid_session_to_anndata
from cellucid.session_bundle import (
    DEFAULT_MAX_UNCOMPRESSED_CHUNK_BYTES,
    SESSION_BUNDLE_MAGIC,
    CellucidSessionBundle,
)
from cellucid.session_codecs import decode_delta_uvarint, decode_user_defined_codes


def _uvarint(value: int) -> bytes:
    output = bytearray()
    remaining = value
    while True:
        byte = remaining & 0x7F
        remaining >>= 7
        if remaining:
            output.append(byte | 0x80)
        else:
            output.append(byte)
            return bytes(output)


def _delta_indices(indices: list[int]) -> bytes:
    output = bytearray(_uvarint(len(indices)))
    previous = 0
    for position, index in enumerate(indices):
        output.extend(_uvarint(index if position == 0 else index - previous))
        previous = index
    return bytes(output)


def _raw_u8_codes(values: list[int]) -> bytes:
    return bytes([0]) + _uvarint(len(values)) + bytes(values)


def _rle_codes(
    encoding: int,
    length: int,
    pairs: list[tuple[int, int]],
) -> bytes:
    output = bytearray([encoding])
    output.extend(_uvarint(length))
    output.extend(_uvarint(len(pairs)))
    for value, run in pairs:
        output.extend(_uvarint(value))
        output.extend(_uvarint(run))
    return bytes(output)


def _chunk_meta(
    chunk_id: str,
    payload: bytes,
    *,
    kind: str = "binary",
    codec: str = "none",
    contributor_id: str = "test",
    priority: str = "eager",
    dataset_dependent: bool = True,
    uncompressed_bytes: int | None = None,
) -> dict[str, Any]:
    return {
        "id": chunk_id,
        "contributorId": contributor_id,
        "priority": priority,
        "kind": kind,
        "codec": codec,
        "label": chunk_id,
        "datasetDependent": dataset_dependent,
        "storedBytes": len(payload),
        "uncompressedBytes": (len(payload) if uncompressed_bytes is None else uncompressed_bytes),
    }


def _write_bundle(
    path: Path,
    entries: list[tuple[dict[str, Any], bytes]],
    *,
    trailing: bytes = b"",
    fingerprint: dict[str, Any] | None = None,
    root_extra: dict[str, Any] | None = None,
) -> Path:
    manifest = {
        "createdAt": "2026-07-25T00:00:00.000Z",
        "datasetFingerprint": fingerprint
        or {
            "sourceType": "anndata",
            "datasetId": "exact",
            "cellCount": 3,
            "varCount": 2,
        },
        "chunks": [meta for meta, _payload in entries],
    }
    if root_extra is not None:
        manifest.update(root_extra)
    manifest_bytes = json.dumps(manifest, separators=(",", ":")).encode("utf-8")
    with path.open("wb") as handle:
        handle.write(SESSION_BUNDLE_MAGIC)
        handle.write(struct.pack("<I", len(manifest_bytes)))
        handle.write(manifest_bytes)
        for _meta, payload in entries:
            handle.write(struct.pack("<I", len(payload)))
            handle.write(payload)
        handle.write(trailing)
    return path


def _adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.array(
            [
                [1.0, 10.0],
                [2.0, 20.0],
                [3.0, 30.0],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(
            {"score": [0.25, 0.5, 0.75]},
            index=pd.Index(["cell-1", "cell-2", "cell-3"], dtype=object),
        ),
        var=pd.DataFrame(index=pd.Index(["gene-a", "gene-b"], dtype=object)),
    )


class _MemoryBundle:
    def __init__(
        self,
        *,
        chunks: dict[str, Any],
        fingerprint: dict[str, Any] | None = None,
        contributors: dict[str, str] | None = None,
        metadata_overrides: dict[str, dict[str, Any]] | None = None,
        chunk_order: list[str] | None = None,
        omit_chunks: set[str] | None = None,
    ):
        defaults: dict[str, Any] = {
            "core/field-overlays": _field_overlays([]),
            "cinematic/camera": {"notDecoded": True},
            "core/state": {"notDecoded": True},
            "ui/dockable-layout": {"notDecoded": True},
            "analysis/windows": {"notDecoded": True},
            "highlights/meta": _empty_highlights(),
            "analysis/cache-inventory": {
                "artifactIds": [
                    chunk_id
                    for chunk_id in chunks
                    if chunk_id.startswith("analysis/artifacts/bulk-gene/")
                ]
            },
        }
        defaults.update(chunks)
        for chunk_id in omit_chunks or set():
            defaults.pop(chunk_id)
        self._chunks = defaults
        self.dataset_fingerprint = fingerprint or {
            "sourceType": "anndata",
            "datasetId": "exact",
            "cellCount": 3,
            "varCount": 2,
        }
        self.manifest = {
            "createdAt": "2026-07-25T00:00:00.000Z",
            "datasetFingerprint": self.dataset_fingerprint,
            "chunks": [],
        }
        metadata = {
            chunk_id: _application_chunk_meta(
                chunk_id,
                payload,
                chunks=defaults,
            )
            for chunk_id, payload in defaults.items()
        }
        if contributors is not None:
            for chunk_id, contributor_id in contributors.items():
                metadata[chunk_id]["contributorId"] = contributor_id
        if metadata_overrides is not None:
            for chunk_id, overrides in metadata_overrides.items():
                metadata[chunk_id].update(overrides)
        if chunk_order is None:
            contributor_order = {
                contributor_id: index
                for index, contributor_id in enumerate(
                    (
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
                )
            }
            chunk_order = sorted(
                defaults,
                key=lambda chunk_id: (
                    0 if metadata[chunk_id]["priority"] == "eager" else 1,
                    contributor_order.get(metadata[chunk_id]["contributorId"], 99),
                    list(defaults).index(chunk_id),
                ),
            )
        if set(chunk_order) != set(defaults) or len(chunk_order) != len(defaults):
            raise AssertionError("chunk_order must contain every fixture chunk exactly once")
        self._chunks = {chunk_id: defaults[chunk_id] for chunk_id in chunk_order}
        self.manifest["chunks"] = [metadata[chunk_id] for chunk_id in chunk_order]

    def list_chunk_ids(self) -> list[str]:
        return list(self._chunks)

    def decode_chunk(self, chunk_id: str) -> Any:
        return self._chunks[chunk_id]


def _application_contributor(chunk_id: str) -> str:
    if chunk_id == "highlights/meta":
        return "highlights-meta"
    if chunk_id.startswith("highlights/cells/"):
        return "highlights-cells"
    if chunk_id == "core/field-overlays":
        return "field-overlays"
    if chunk_id == "cinematic/camera":
        return "cinematic-camera"
    if chunk_id == "core/state":
        return "core-state"
    if chunk_id == "ui/dockable-layout":
        return "dockable-layout"
    if chunk_id == "analysis/windows":
        return "analysis-windows"
    if chunk_id == "analysis/cache-inventory":
        return "analysis-artifacts"
    if chunk_id.startswith("user-defined/codes/"):
        return "user-defined-codes"
    if chunk_id.startswith("analysis/artifacts/bulk-gene/"):
        return "analysis-artifacts"
    return "unknown-contributor"


_STATIC_CHUNK_PROFILES: dict[str, dict[str, Any]] = {
    "core/field-overlays": {
        "contributorId": "field-overlays",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Field overlays",
        "datasetDependent": True,
    },
    "cinematic/camera": {
        "contributorId": "cinematic-camera",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Cinematic camera path",
        "datasetDependent": True,
    },
    "core/state": {
        "contributorId": "core-state",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Core state",
        "datasetDependent": True,
    },
    "ui/dockable-layout": {
        "contributorId": "dockable-layout",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Floating panels",
        "datasetDependent": False,
    },
    "analysis/windows": {
        "contributorId": "analysis-windows",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Analysis windows",
        "datasetDependent": True,
    },
    "highlights/meta": {
        "contributorId": "highlights-meta",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Highlight metadata",
        "datasetDependent": True,
    },
    "analysis/cache-inventory": {
        "contributorId": "analysis-artifacts",
        "priority": "eager",
        "kind": "json",
        "codec": "gzip",
        "label": "Analysis cache inventory",
        "datasetDependent": True,
    },
}


def _empty_highlights() -> dict[str, Any]:
    return {
        "pages": [
            {
                "id": "page_1",
                "name": "Page 1",
                "color": "#2563eb",
                "highlightedGroups": [],
            }
        ],
        "activePageId": "page_1",
    }


def _json_payload_size(payload: Any) -> int:
    try:
        return len(
            json.dumps(
                payload,
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode("utf-8")
        )
    except TypeError:
        return 1


def _application_chunk_meta(
    chunk_id: str,
    payload: Any,
    *,
    chunks: dict[str, Any],
) -> dict[str, Any]:
    profile = _STATIC_CHUNK_PROFILES.get(chunk_id)
    if profile is not None:
        base = dict(profile)
    elif chunk_id.startswith("user-defined/codes/"):
        field_id = chunk_id.removeprefix("user-defined/codes/")
        overlays = chunks["core/field-overlays"]
        field_key = next(
            (field["key"] for field in overlays["userDefinedFields"] if field["id"] == field_id),
            field_id,
        )
        base = {
            "contributorId": "user-defined-codes",
            "priority": "eager",
            "kind": "binary",
            "codec": "gzip",
            "label": f"User-defined codes: {field_key}",
            "datasetDependent": True,
        }
    elif chunk_id.startswith("highlights/cells/"):
        group_id = chunk_id.removeprefix("highlights/cells/")
        highlights = chunks["highlights/meta"]
        group_label = next(
            (
                group["label"]
                for page in highlights["pages"]
                for group in page["highlightedGroups"]
                if group["id"] == group_id
            ),
            group_id,
        )
        base = {
            "contributorId": "highlights-cells",
            "priority": "lazy",
            "kind": "binary",
            "codec": "gzip",
            "label": f"Highlight cells: {group_label}",
            "datasetDependent": True,
        }
    elif chunk_id.startswith("analysis/artifacts/bulk-gene/"):
        base = {
            "contributorId": "analysis-artifacts",
            "priority": "lazy",
            "kind": "binary",
            "codec": "gzip",
            "label": "Analysis cache: unmaterialized",
            "datasetDependent": True,
        }
    else:
        base = {
            "contributorId": _application_contributor(chunk_id),
            "priority": "eager",
            "kind": "binary",
            "codec": "gzip",
            "label": chunk_id,
            "datasetDependent": True,
        }
    uncompressed_bytes = (
        len(payload)
        if isinstance(payload, bytes | bytearray | memoryview)
        else _json_payload_size(payload)
    )
    return {
        "id": chunk_id,
        **base,
        "storedBytes": uncompressed_bytes,
        "uncompressedBytes": uncompressed_bytes,
    }


def _field_overlays(fields: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "renames": {"fields": {}, "categories": {}},
        "deletedFields": {"deleted": [], "purged": []},
        "userDefinedFields": fields,
    }


def _categorical_field(
    *,
    field_id: str = "user_cat_1",
    source: str = "obs",
    key: str = "status",
    categories: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "id": field_id,
        "source": source,
        "kind": "category",
        "key": key,
        "categories": ["A", "B"] if categories is None else categories,
        "isDeleted": False,
        "isPurged": False,
        "codesLength": 3,
        "codesType": "Uint8Array",
        "centroidsByDim": {},
        "normalizedDims": [],
        "sourceField": None,
        "operation": None,
        "sourcePages": [],
        "overlapStrategy": "first",
        "overlapLabel": None,
        "intersectionLabels": None,
        "uncoveredLabel": None,
        "createdAt": 1,
    }


def _continuous_field(
    *,
    field_id: str,
    source: str,
    key: str,
    source_key: str,
) -> dict[str, Any]:
    return {
        "id": field_id,
        "source": source,
        "kind": "continuous",
        "key": key,
        "isDeleted": False,
        "isPurged": False,
        "sourceField": {
            "sourceKey": source_key,
            "sourceIndex": 0 if source == "obs" else 1,
            "kind": "continuous",
        },
        "operation": {
            "type": "copy-field",
            "source": source,
            "sourceKey": source_key,
        },
        "createdAt": 1,
    }


def test_session_apply_exposes_only_the_exact_failure_contract() -> None:
    parameters = inspect.signature(apply_cellucid_session_to_anndata).parameters
    assert "dataset_mismatch" not in parameters
    assert "column_conflict" not in parameters
    assert parameters["expected_dataset_id"].default is inspect.Parameter.empty

    bundle_parameters = inspect.signature(CellucidSessionBundle.apply_to_anndata).parameters
    assert not any(
        parameter.kind is inspect.Parameter.VAR_KEYWORD for parameter in bundle_parameters.values()
    )
    assert bundle_parameters["expected_dataset_id"].default is inspect.Parameter.empty


def test_session_apply_rejects_dataset_mismatch_without_mutation() -> None:
    adata = _adata()
    before_obs = adata.obs.copy(deep=True)
    before_uns = copy.deepcopy(adata.uns)
    bundle = _MemoryBundle(
        chunks={},
        fingerprint={
            "sourceType": "anndata",
            "datasetId": "other",
            "cellCount": 3,
            "varCount": 2,
        },
    )

    with pytest.raises(ValueError, match="datasetId"):
        apply_cellucid_session_to_anndata(
            bundle,
            adata,
            inplace=True,
            expected_dataset_id="exact",
        )

    pd.testing.assert_frame_equal(adata.obs, before_obs)
    assert adata.uns == before_uns


@pytest.mark.parametrize(
    ("chunk_id", "contributor_id"),
    [
        ("unknown/state", "unknown-contributor"),
        ("core/state", "unknown-contributor"),
        ("unknown/state", "core-state"),
        ("highlights/meta", "highlights-cells"),
        ("highlights/cells/highlight_1", "highlights-meta"),
        ("user-defined/codes/user_cat_1", "field-overlays"),
    ],
)
def test_session_apply_rejects_unknown_or_mismatched_chunk_contributors_before_decode(
    chunk_id: str,
    contributor_id: str,
) -> None:
    adata = _adata()
    before_obs = adata.obs.copy(deep=True)
    before_uns = copy.deepcopy(adata.uns)
    bundle = _MemoryBundle(
        chunks={chunk_id: object()},
        contributors={chunk_id: contributor_id},
    )

    def forbidden_decode(_chunk_id: str) -> Any:
        raise AssertionError("invalid session inventory must fail before payload decode")

    bundle.decode_chunk = forbidden_decode  # type: ignore[method-assign]

    with pytest.raises(ValueError, match="chunk|contributor"):
        apply_cellucid_session_to_anndata(
            bundle,
            adata,
            inplace=True,
            expected_dataset_id="exact",
        )

    pd.testing.assert_frame_equal(adata.obs, before_obs)
    assert adata.uns == before_uns


def test_session_apply_accepts_current_nonmaterialized_chunk_inventory() -> None:
    analysis_chunk = "analysis/artifacts/bulk-gene/bulk_genes%3Apage%2F%CE%B1/IL%2F7/page%2F1"
    bundle = _MemoryBundle(
        chunks={
            "core/state": object(),
            "ui/dockable-layout": object(),
            "analysis/windows": object(),
            "cinematic/camera": object(),
            analysis_chunk: object(),
        },
        contributors={
            "core/state": "core-state",
            "ui/dockable-layout": "dockable-layout",
            "analysis/windows": "analysis-windows",
            "cinematic/camera": "cinematic-camera",
            analysis_chunk: "analysis-artifacts",
        },
    )

    decoded_chunk_ids: list[str] = []
    decode_chunk = bundle.decode_chunk

    def guarded_decode(chunk_id: str) -> Any:
        if chunk_id not in {
            "analysis/cache-inventory",
            "core/field-overlays",
            "highlights/meta",
        }:
            raise AssertionError(f"nonmaterialized current chunk {chunk_id!r} must not be decoded")
        decoded_chunk_ids.append(chunk_id)
        return decode_chunk(chunk_id)

    bundle.decode_chunk = guarded_decode  # type: ignore[method-assign]
    output = apply_cellucid_session_to_anndata(
        bundle,
        _adata(),
        expected_dataset_id="exact",
        add_highlights=False,
        add_user_defined_fields=False,
        store_uns=False,
    )

    assert output.n_obs == 3
    assert output.n_vars == 2
    assert decoded_chunk_ids == [
        "analysis/cache-inventory",
        "highlights/meta",
        "core/field-overlays",
    ]


@pytest.mark.parametrize(
    "missing_chunk",
    ["analysis/cache-inventory", "cinematic/camera"],
)
def test_session_apply_requires_current_cache_and_camera_singletons_before_decode(
    missing_chunk: str,
) -> None:
    adata = _adata()
    before_obs = adata.obs.copy(deep=True)
    before_uns = copy.deepcopy(adata.uns)
    bundle = _MemoryBundle(chunks={}, omit_chunks={missing_chunk})

    def forbidden_decode(_chunk_id: str) -> Any:
        raise AssertionError("missing singleton inventory must fail before payload decode")

    bundle.decode_chunk = forbidden_decode  # type: ignore[method-assign]
    with pytest.raises(ValueError, match=missing_chunk):
        apply_cellucid_session_to_anndata(
            bundle,
            adata,
            inplace=True,
            expected_dataset_id="exact",
        )

    pd.testing.assert_frame_equal(adata.obs, before_obs)
    assert adata.uns == before_uns


@pytest.mark.parametrize(
    ("chunk_id", "overrides"),
    [
        ("core/state", {"priority": "lazy"}),
        ("core/state", {"kind": "binary"}),
        ("core/state", {"codec": "none"}),
        ("core/state", {"label": "State"}),
        ("ui/dockable-layout", {"datasetDependent": True}),
        ("analysis/cache-inventory", {"label": "Cache"}),
        ("cinematic/camera", {"label": "Cinematic camera"}),
    ],
)
def test_session_apply_rejects_static_profile_metadata_mutations_before_decode(
    chunk_id: str,
    overrides: dict[str, Any],
) -> None:
    bundle = _MemoryBundle(
        chunks={},
        metadata_overrides={chunk_id: overrides},
    )

    def forbidden_decode(_chunk_id: str) -> Any:
        raise AssertionError("static profile mutations must fail before payload decode")

    bundle.decode_chunk = forbidden_decode  # type: ignore[method-assign]
    with pytest.raises(ValueError, match="profile|eager chunks"):
        apply_cellucid_session_to_anndata(
            bundle,
            _adata(),
            expected_dataset_id="exact",
        )


def test_session_apply_rejects_registered_contributor_reordering_before_decode() -> None:
    current_order = list(_MemoryBundle(chunks={}).list_chunk_ids())
    camera_index = current_order.index("cinematic/camera")
    core_index = current_order.index("core/state")
    current_order[camera_index], current_order[core_index] = (
        current_order[core_index],
        current_order[camera_index],
    )
    bundle = _MemoryBundle(chunks={}, chunk_order=current_order)

    def forbidden_decode(_chunk_id: str) -> Any:
        raise AssertionError("contributor reordering must fail before payload decode")

    bundle.decode_chunk = forbidden_decode  # type: ignore[method-assign]
    with pytest.raises(ValueError, match="registered order"):
        apply_cellucid_session_to_anndata(
            bundle,
            _adata(),
            expected_dataset_id="exact",
        )


def test_session_apply_requires_exact_analysis_inventory_order_without_artifact_decode() -> None:
    artifact_a = "analysis/artifacts/bulk-gene/cache-a/gene-a/page-a"
    artifact_b = "analysis/artifacts/bulk-gene/cache-b/gene-b/page-b"
    bundle = _MemoryBundle(
        chunks={
            "analysis/cache-inventory": {
                "artifactIds": [artifact_b, artifact_a],
            },
            artifact_a: object(),
            artifact_b: object(),
        }
    )
    decode_chunk = bundle.decode_chunk
    decoded_ids: list[str] = []

    def inventory_only_decode(chunk_id: str) -> Any:
        if chunk_id != "analysis/cache-inventory":
            raise AssertionError("analysis artifact payloads must remain nonmaterialized")
        decoded_ids.append(chunk_id)
        return decode_chunk(chunk_id)

    bundle.decode_chunk = inventory_only_decode  # type: ignore[method-assign]
    with pytest.raises(ValueError, match="cache inventory order"):
        apply_cellucid_session_to_anndata(
            bundle,
            _adata(),
            expected_dataset_id="exact",
        )
    assert decoded_ids == ["analysis/cache-inventory"]


def test_session_apply_requires_user_code_order_from_field_overlay_inventory() -> None:
    field_a = _categorical_field(field_id="field-a", key="A")
    field_b = _categorical_field(field_id="field-b", key="B")
    bundle = _MemoryBundle(
        chunks={
            "core/field-overlays": _field_overlays([field_a, field_b]),
            "user-defined/codes/field-b": _raw_u8_codes([0, 1, 0]),
            "user-defined/codes/field-a": _raw_u8_codes([1, 0, 1]),
        }
    )
    with pytest.raises(ValueError, match="field-overlay inventory order"):
        apply_cellucid_session_to_anndata(
            bundle,
            _adata(),
            expected_dataset_id="exact",
            store_uns=False,
        )


@pytest.mark.parametrize(
    ("chunk_id", "overrides"),
    [
        ("user-defined/codes/user_cat_1", {"label": "User-defined codes"}),
        ("user-defined/codes/user_cat_1", {"kind": "json"}),
        ("user-defined/codes/user_cat_1", {"codec": "none"}),
        ("user-defined/codes/user_cat_1", {"datasetDependent": False}),
        ("user-defined/codes/user_cat_1", {"uncompressedBytes": 99}),
        ("highlights/cells/highlight_1", {"label": "Highlight cells"}),
        ("highlights/cells/highlight_1", {"priority": "eager"}),
        ("highlights/cells/highlight_1", {"kind": "json"}),
        ("highlights/cells/highlight_1", {"codec": "none"}),
        ("highlights/cells/highlight_1", {"datasetDependent": False}),
        ("highlights/cells/highlight_1", {"uncompressedBytes": 99}),
    ],
)
def test_session_apply_rejects_materialized_dynamic_metadata_mutations_atomically(
    chunk_id: str,
    overrides: dict[str, Any],
) -> None:
    adata = _adata()
    before_obs = adata.obs.copy(deep=True)
    before_uns = copy.deepcopy(adata.uns)
    bundle = _MemoryBundle(
        chunks={
            "core/field-overlays": _field_overlays([_categorical_field()]),
            "user-defined/codes/user_cat_1": _raw_u8_codes([0, 1, 0]),
            "highlights/meta": {
                "pages": [
                    {
                        "id": "page_1",
                        "name": "Page 1",
                        "color": "#2563eb",
                        "highlightedGroups": [
                            {
                                "id": "highlight_1",
                                "type": "lasso",
                                "label": "Selected",
                                "enabled": True,
                                "cellCount": 1,
                            }
                        ],
                    }
                ],
                "activePageId": "page_1",
            },
            "highlights/cells/highlight_1": _delta_indices([1]),
        },
        metadata_overrides={chunk_id: overrides},
    )

    with pytest.raises((TypeError, ValueError)):
        apply_cellucid_session_to_anndata(
            bundle,
            adata,
            inplace=True,
            expected_dataset_id="exact",
        )

    pd.testing.assert_frame_equal(adata.obs, before_obs)
    assert adata.uns == before_uns


@pytest.mark.parametrize(
    "overrides",
    [
        {"priority": "eager"},
        {"kind": "json"},
        {"codec": "none"},
        {"datasetDependent": False},
    ],
)
def test_session_apply_rejects_analysis_artifact_metadata_without_heavy_decode(
    overrides: dict[str, Any],
) -> None:
    artifact_id = "analysis/artifacts/bulk-gene/cache/gene/page"
    bundle = _MemoryBundle(
        chunks={artifact_id: object()},
        metadata_overrides={artifact_id: overrides},
    )

    def forbidden_decode(_chunk_id: str) -> Any:
        raise AssertionError("invalid analysis metadata must fail before payload decode")

    bundle.decode_chunk = forbidden_decode  # type: ignore[method-assign]
    with pytest.raises(ValueError, match="lazy, binary, gzip"):
        apply_cellucid_session_to_anndata(
            bundle,
            _adata(),
            expected_dataset_id="exact",
        )


@pytest.mark.parametrize(
    "highlights",
    [
        {
            "pages": [],
            "activePageId": "page_1",
        },
        {
            "pages": [
                {
                    "id": "page-1",
                    "name": "Page 1",
                    "color": "#2563eb",
                    "highlightedGroups": [],
                }
            ],
            "activePageId": "page-1",
        },
        {
            "pages": [
                {
                    "id": "page_1",
                    "name": "Page 1",
                    "color": None,
                    "highlightedGroups": [],
                }
            ],
            "activePageId": "page_1",
        },
        {
            "pages": [
                {
                    "id": "page_1",
                    "name": "Page 1",
                    "color": "#2563eb",
                    "highlightedGroups": [
                        {
                            "id": "highlight_1",
                            "type": "lasso",
                            "label": "Empty",
                            "enabled": True,
                            "cellCount": 0,
                        }
                    ],
                }
            ],
            "activePageId": "page_1",
        },
        {
            "pages": [
                {
                    "id": "page_1",
                    "name": "Page 1",
                    "color": "#2563eb",
                    "highlightedGroups": [
                        {
                            "id": "highlight_1",
                            "type": "lasso",
                            "label": "Selected",
                            "enabled": True,
                            "cellCount": 1,
                            "fieldKey": "not-valid-for-lasso",
                        }
                    ],
                }
            ],
            "activePageId": "page_1",
        },
    ],
)
def test_session_apply_rejects_noncurrent_highlight_metadata(
    highlights: dict[str, Any],
) -> None:
    bundle = _MemoryBundle(chunks={"highlights/meta": highlights})
    with pytest.raises((TypeError, ValueError)):
        apply_cellucid_session_to_anndata(
            bundle,
            _adata(),
            expected_dataset_id="exact",
            store_uns=False,
        )


@pytest.mark.parametrize(
    "chunk_id",
    [
        "analysis/artifacts/bulk-gene/cache/gene",
        "analysis/artifacts/bulk-gene/cache/gene/page/extra",
        "analysis/artifacts/bulk-gene/cache/gene/",
        "analysis/artifacts/bulk-gene/cache/gene%2fvalue/page",
        "analysis/artifacts/bulk-gene/cache/gene value/page",
        "analysis/artifacts/bulk-gene/cache/%FF/page",
    ],
)
def test_session_apply_rejects_noncanonical_analysis_artifact_chunk_ids(
    chunk_id: str,
) -> None:
    bundle = _MemoryBundle(
        chunks={chunk_id: object()},
        contributors={chunk_id: "analysis-artifacts"},
    )

    def forbidden_decode(_chunk_id: str) -> Any:
        raise AssertionError("invalid analysis chunk identity must fail before decode")

    bundle.decode_chunk = forbidden_decode  # type: ignore[method-assign]
    with pytest.raises(ValueError, match="Unknown current session chunk"):
        apply_cellucid_session_to_anndata(
            bundle,
            _adata(),
            expected_dataset_id="exact",
        )


def test_noninplace_session_apply_rejects_backed_anndata_before_copy(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "backed-session-target.h5ad"
    _adata().write_h5ad(source)
    backed = ad.read_h5ad(source, backed="r")
    before_obs = backed.obs.copy(deep=True)
    before_uns = copy.deepcopy(backed.uns)

    def forbidden_copy(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("backed AnnData.copy() must not be attempted")

    monkeypatch.setattr(ad.AnnData, "copy", forbidden_copy)
    try:
        with pytest.raises(ValueError, match="backed AnnData.*inplace=False"):
            apply_cellucid_session_to_anndata(
                _MemoryBundle(chunks={}),
                backed,
                expected_dataset_id="exact",
                inplace=False,
            )
        assert backed.isbacked is True
        pd.testing.assert_frame_equal(backed.obs, before_obs)
        assert backed.uns == before_uns
    finally:
        backed.file.close()


def test_inplace_session_apply_keeps_backed_anndata_backed_and_does_not_write_file(
    tmp_path: Path,
) -> None:
    source = tmp_path / "backed-inplace-session-target.h5ad"
    _adata().write_h5ad(source)
    backed = ad.read_h5ad(source, backed="r")
    bundle = _MemoryBundle(
        chunks={
            "highlights/meta": {
                "pages": [
                    {
                        "id": "page_1",
                        "name": "Page 1",
                        "color": "#2563eb",
                        "highlightedGroups": [
                            {
                                "id": "highlight_1",
                                "type": "lasso",
                                "label": "Selected",
                                "enabled": True,
                                "cellCount": 1,
                            }
                        ],
                    }
                ],
                "activePageId": "page_1",
            },
            "highlights/cells/highlight_1": _delta_indices([1]),
        }
    )

    try:
        output = apply_cellucid_session_to_anndata(
            bundle,
            backed,
            expected_dataset_id="exact",
            inplace=True,
        )
        assert output is backed
        assert output.isbacked is True
        assert output.obs["cellucid_highlight__highlight_1"].tolist() == [
            False,
            True,
            False,
        ]
    finally:
        backed.file.close()

    reopened = ad.read_h5ad(source)
    assert "cellucid_highlight__highlight_1" not in reopened.obs
    assert "cellucid" not in reopened.uns


def test_highlight_metadata_rejects_obsolete_active_page_name() -> None:
    bundle = _MemoryBundle(
        chunks={
            "highlights/meta": {
                "pages": [
                    {
                        "id": "page_1",
                        "name": "Page 1",
                        "color": "#2563eb",
                        "highlightedGroups": [],
                    }
                ],
                "activePageId": "page_1",
                "activePageName": "Page 1",
            }
        }
    )

    with pytest.raises(
        ValueError,
        match="highlights/meta.*unknown activePageName",
    ):
        apply_cellucid_session_to_anndata(
            bundle,
            _adata(),
            expected_dataset_id="exact",
            store_uns=False,
        )


def test_session_apply_rejects_column_conflicts_instead_of_renaming() -> None:
    adata = _adata()
    before_obs = adata.obs.copy(deep=True)
    bundle = _MemoryBundle(
        chunks={
            "core/field-overlays": _field_overlays([_categorical_field(key="score")]),
            "user-defined/codes/user_cat_1": _raw_u8_codes([0, 1, 0]),
        }
    )

    with pytest.raises(ValueError, match="Column already exists: score"):
        apply_cellucid_session_to_anndata(
            bundle,
            adata,
            inplace=True,
            expected_dataset_id="exact",
            add_highlights=False,
            store_uns=False,
        )

    pd.testing.assert_frame_equal(adata.obs, before_obs)
    assert "score__2" not in adata.obs


def test_session_apply_validates_every_chunk_before_any_inplace_mutation() -> None:
    adata = _adata()
    adata.uns["existing"] = {"value": 7}
    before_obs = adata.obs.copy(deep=True)
    before_uns = copy.deepcopy(adata.uns)
    highlight_meta = {
        "pages": [
            {
                "id": "page_1",
                "name": "Page 1",
                "color": "#ff0000",
                "highlightedGroups": [
                    {
                        "id": "highlight_1",
                        "type": "lasso",
                        "label": "Selected",
                        "enabled": True,
                        "cellCount": 1,
                    }
                ],
            }
        ],
        "activePageId": "page_1",
    }
    bundle = _MemoryBundle(
        chunks={
            "highlights/meta": highlight_meta,
            "highlights/cells/highlight_1": _delta_indices([1]),
            "core/field-overlays": _field_overlays([_categorical_field(categories=[])]),
            "user-defined/codes/user_cat_1": _raw_u8_codes([0, 0, 0]),
        }
    )

    with pytest.raises(ValueError, match="categories.*non-empty"):
        apply_cellucid_session_to_anndata(
            bundle,
            adata,
            inplace=True,
            expected_dataset_id="exact",
        )

    pd.testing.assert_frame_equal(adata.obs, before_obs)
    assert adata.uns == before_uns
    assert "cellucid_highlight__highlight_1" not in adata.obs


def test_highlight_metadata_requires_its_exact_membership_chunk() -> None:
    bundle = _MemoryBundle(
        chunks={
            "highlights/meta": {
                "pages": [
                    {
                        "id": "page_1",
                        "name": "Page 1",
                        "color": "#2563eb",
                        "highlightedGroups": [
                            {
                                "id": "highlight_1",
                                "type": "lasso",
                                "label": "Selected",
                                "enabled": True,
                                "cellCount": 1,
                            }
                        ],
                    }
                ],
                "activePageId": "page_1",
            }
        }
    )

    with pytest.raises(ValueError, match="highlights/cells/highlight_1"):
        apply_cellucid_session_to_anndata(
            bundle,
            _adata(),
            expected_dataset_id="exact",
            store_uns=False,
        )


def test_user_defined_category_preserves_exact_key_and_is_cell_aligned() -> None:
    bundle = _MemoryBundle(
        chunks={
            "core/field-overlays": _field_overlays(
                [
                    _categorical_field(
                        source="var",
                        key="Exact cell label",
                    )
                ]
            ),
            "user-defined/codes/user_cat_1": _raw_u8_codes([0, 1, 0]),
        }
    )

    output = apply_cellucid_session_to_anndata(
        bundle,
        _adata(),
        expected_dataset_id="exact",
        add_highlights=False,
        store_uns=False,
    )

    assert "Exact cell label" in output.obs
    assert output.obs["Exact cell label"].tolist() == ["A", "B", "A"]
    assert "Exact_cell_label" not in output.obs
    assert "Exact cell label" not in output.var


def test_user_defined_category_materializes_current_boolean_category_values() -> None:
    bundle = _MemoryBundle(
        chunks={
            "core/field-overlays": _field_overlays(
                [_categorical_field(key="Boolean status", categories=[False, True])]
            ),
            "user-defined/codes/user_cat_1": _raw_u8_codes([0, 1, 0]),
        }
    )

    output = apply_cellucid_session_to_anndata(
        bundle,
        _adata(),
        expected_dataset_id="exact",
        add_highlights=False,
        store_uns=False,
    )

    assert output.obs["Boolean status"].tolist() == [False, True, False]


def test_user_defined_category_materializes_current_rle_uint16_codes() -> None:
    field = _categorical_field(key="RLE Uint16 status")
    field["codesType"] = "Uint16Array"
    bundle = _MemoryBundle(
        chunks={
            "core/field-overlays": _field_overlays([field]),
            "user-defined/codes/user_cat_1": _rle_codes(
                3,
                3,
                [(1, 3)],
            ),
        }
    )

    output = apply_cellucid_session_to_anndata(
        bundle,
        _adata(),
        expected_dataset_id="exact",
        add_highlights=False,
        store_uns=False,
    )

    assert output.obs["RLE Uint16 status"].tolist() == ["B", "B", "B"]


def test_user_defined_category_rejects_rle_length_before_allocation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    field = _categorical_field(key="Bounded RLE status")
    field["codesType"] = "Uint16Array"
    bundle = _MemoryBundle(
        chunks={
            "core/field-overlays": _field_overlays([field]),
            "user-defined/codes/user_cat_1": _rle_codes(
                3,
                0xFFFF_FFFF,
                [(1, 0xFFFF_FFFF)],
            ),
        }
    )
    adata = _adata()

    def forbidden_empty(*_args: Any, **_kwargs: Any) -> np.ndarray:
        raise AssertionError("mismatched RLE length must reject before allocation")

    with (
        pytest.raises(ValueError, match="declared length.*expected length"),
        monkeypatch.context() as bounded,
    ):
        bounded.setattr("cellucid.session_codecs.np.empty", forbidden_empty)
        apply_cellucid_session_to_anndata(
            bundle,
            adata,
            expected_dataset_id="exact",
            add_highlights=False,
            store_uns=False,
        )


@pytest.mark.parametrize(
    ("source", "source_key", "expected"),
    [
        ("obs", "score", [0.25, 0.5, 0.75]),
        ("var", "gene-b", [10.0, 20.0, 30.0]),
    ],
)
def test_continuous_aliases_materialize_exact_cell_values(
    source: str,
    source_key: str,
    expected: list[float],
) -> None:
    field_id = f"user_{source}_cont_1"
    bundle = _MemoryBundle(
        chunks={
            "core/field-overlays": _field_overlays(
                [
                    _continuous_field(
                        field_id=field_id,
                        source=source,
                        key=f"{source} copy",
                        source_key=source_key,
                    )
                ]
            )
        }
    )

    output = apply_cellucid_session_to_anndata(
        bundle,
        _adata(),
        expected_dataset_id="exact",
        add_highlights=False,
        store_uns=False,
    )

    assert output.obs[f"{source} copy"].tolist() == expected


def test_bundle_rejects_duplicate_chunk_ids_and_trailing_bytes(
    tmp_path: Path,
) -> None:
    payload = b"one"
    meta = _chunk_meta("chunk", payload)
    duplicate = _write_bundle(
        tmp_path / "duplicate.cellucid-session",
        [(meta, payload), (dict(meta), payload)],
    )
    with pytest.raises(ValueError, match="duplicate chunk id"):
        CellucidSessionBundle(duplicate).list_chunk_ids()

    trailing = _write_bundle(
        tmp_path / "trailing.cellucid-session",
        [(meta, payload)],
        trailing=b"unexpected",
    )
    with pytest.raises(ValueError, match="trailing bytes"):
        CellucidSessionBundle(trailing).list_chunk_ids()


def test_bundle_enforces_web_writer_exact_current_root_manifest(
    tmp_path: Path,
) -> None:
    current = _write_bundle(
        tmp_path / "current.cellucid-session",
        [],
    )
    assert CellucidSessionBundle(current).list_chunk_ids() == []

    for obsolete_field in ("dataSource", "summary"):
        obsolete = _write_bundle(
            tmp_path / f"obsolete-{obsolete_field}.cellucid-session",
            [],
            root_extra={obsolete_field: None},
        )
        with pytest.raises(
            ValueError,
            match=rf"unknown {obsolete_field}",
        ):
            CellucidSessionBundle(obsolete).list_chunk_ids()


def test_bundle_accepts_web_percent_encoded_analysis_chunk_identity(
    tmp_path: Path,
) -> None:
    chunk_id = (
        "analysis/artifacts/bulk-gene/"
        "bulk_genes%3Apage_1%3AIL-7/%CE%B1/page_1"
    )
    raw = b"current analysis artifact"
    compressed = gzip.compress(raw, mtime=0)
    meta = _chunk_meta(
        chunk_id,
        compressed,
        codec="gzip",
        contributor_id="analysis-artifacts",
        priority="lazy",
        uncompressed_bytes=len(raw),
    )
    path = _write_bundle(
        tmp_path / "percent-encoded-analysis.cellucid-session",
        [(meta, compressed)],
    )

    bundle = CellucidSessionBundle(path)
    assert bundle.list_chunk_ids() == [chunk_id]
    assert bundle.decode_chunk(chunk_id) == raw


def test_on_disk_percent_encoded_analysis_inventory_applies_without_artifact_decode(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact_id = (
        "analysis/artifacts/bulk-gene/"
        "bulk_genes%3Apage_1%3AIL-7/%CE%B1/page_1"
    )
    payloads: dict[str, Any] = {
        "core/field-overlays": _field_overlays([]),
        "cinematic/camera": {"opaque": True},
        "core/state": {"opaque": True},
        "ui/dockable-layout": {"opaque": True},
        "analysis/windows": {"opaque": True},
        "highlights/meta": _empty_highlights(),
        "analysis/cache-inventory": {"artifactIds": [artifact_id]},
        artifact_id: b"opaque heavy analysis artifact",
    }
    entries: list[tuple[dict[str, Any], bytes]] = []
    for chunk_id, payload in payloads.items():
        raw = (
            bytes(payload)
            if isinstance(payload, bytes | bytearray | memoryview)
            else json.dumps(
                payload,
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode("utf-8")
        )
        stored = gzip.compress(raw, mtime=0)
        meta = _application_chunk_meta(chunk_id, payload, chunks=payloads)
        meta["storedBytes"] = len(stored)
        meta["uncompressedBytes"] = len(raw)
        entries.append((meta, stored))
    path = _write_bundle(
        tmp_path / "complete-percent-analysis.cellucid-session",
        entries,
    )
    bundle = CellucidSessionBundle(path)
    decode_chunk = bundle.decode_chunk
    decoded_ids: list[str] = []

    def inventory_only_decode(chunk_id: str) -> Any:
        if chunk_id == artifact_id:
            raise AssertionError("heavy analysis artifact must remain opaque")
        decoded_ids.append(chunk_id)
        return decode_chunk(chunk_id)

    monkeypatch.setattr(bundle, "decode_chunk", inventory_only_decode)
    output = apply_cellucid_session_to_anndata(
        bundle,
        _adata(),
        expected_dataset_id="exact",
        add_highlights=False,
        add_user_defined_fields=False,
        store_uns=False,
    )

    assert output.n_obs == 3
    assert decoded_ids == [
        "analysis/cache-inventory",
        "highlights/meta",
        "core/field-overlays",
    ]


@pytest.mark.parametrize(
    "fingerprint",
    [
        {
            "sourceType": "anndata",
            "datasetId": "exact",
            "cellCount": 3,
        },
        {
            "sourceType": "anndata",
            "datasetId": "exact",
            "cellCount": 3,
            "varCount": 2,
            "legacy": True,
        },
    ],
)
def test_bundle_requires_exact_four_key_dataset_fingerprint(
    tmp_path: Path,
    fingerprint: dict[str, Any],
) -> None:
    path = _write_bundle(
        tmp_path / "fingerprint.cellucid-session",
        [],
        fingerprint=fingerprint,
    )
    with pytest.raises(ValueError, match="datasetFingerprint fields"):
        CellucidSessionBundle(path).list_chunk_ids()


def test_bundle_accepts_exact_nullable_web_dataset_identity(
    tmp_path: Path,
) -> None:
    path = _write_bundle(
        tmp_path / "nullable-fingerprint.cellucid-session",
        [],
        fingerprint={
            "sourceType": None,
            "datasetId": None,
            "cellCount": 0,
            "varCount": 0,
        },
    )
    assert CellucidSessionBundle(path).dataset_fingerprint == {
        "sourceType": None,
        "datasetId": None,
        "cellCount": 0,
        "varCount": 0,
    }


def test_bundle_rejects_removed_depends_on_chunk_metadata(
    tmp_path: Path,
) -> None:
    payload = b"payload"
    meta = _chunk_meta("chunk", payload)
    meta["dependsOn"] = []
    path = _write_bundle(
        tmp_path / "depends-on.cellucid-session",
        [(meta, payload)],
    )
    with pytest.raises(ValueError, match="unknown dependsOn"):
        CellucidSessionBundle(path).list_chunk_ids()


def test_bundle_requires_exact_stored_and_uncompressed_lengths(
    tmp_path: Path,
) -> None:
    raw = b'{"value":1}'
    wrong_stored = _chunk_meta("raw", raw, kind="json")
    wrong_stored["storedBytes"] += 1
    stored_path = _write_bundle(
        tmp_path / "stored.cellucid-session",
        [(wrong_stored, raw)],
    )
    with pytest.raises(ValueError, match="storedBytes"):
        CellucidSessionBundle(stored_path).list_chunk_ids()

    compressed = gzip.compress(raw, mtime=0)
    wrong_uncompressed = _chunk_meta(
        "gzip",
        compressed,
        kind="json",
        codec="gzip",
        uncompressed_bytes=len(raw) + 1,
    )
    compressed_path = _write_bundle(
        tmp_path / "uncompressed.cellucid-session",
        [(wrong_uncompressed, compressed)],
    )
    with pytest.raises(ValueError, match="uncompressedBytes"):
        CellucidSessionBundle(compressed_path).decode_chunk("gzip")


def test_bundle_rejects_oversized_gzip_before_decompression(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    compressed = gzip.compress(b"{}", mtime=0)
    meta = _chunk_meta(
        "gzip",
        compressed,
        kind="json",
        codec="gzip",
        uncompressed_bytes=DEFAULT_MAX_UNCOMPRESSED_CHUNK_BYTES + 1,
    )
    path = _write_bundle(
        tmp_path / "oversized.cellucid-session",
        [(meta, compressed)],
    )

    def forbidden_decompress(_payload: bytes) -> bytes:
        raise AssertionError("gzip.decompress must not run")

    monkeypatch.setattr(
        "cellucid.session_bundle._decompress_gzip_exact",
        forbidden_decompress,
    )
    with pytest.raises(ValueError, match="uncompressedBytes.*limit"):
        CellucidSessionBundle(path).decode_chunk("gzip")


def test_delta_codec_rejects_duplicates_and_trailing_bytes() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        decode_delta_uvarint(_delta_indices([1, 1]), max_count=3, max_index=2)

    with pytest.raises(ValueError, match="trailing bytes"):
        decode_delta_uvarint(_delta_indices([1]) + b"\x00", max_count=3, max_index=2)


@pytest.mark.parametrize(
    ("payload", "expected_length", "expected_codes_type", "message"),
    [
        (_raw_u8_codes([0, 1]) + b"\x00", 2, "Uint8Array", "trailing bytes"),
        (
            bytes([2]) + _uvarint(1) + _uvarint(1) + _uvarint(0) + _uvarint(0),
            1,
            "Uint8Array",
            "run",
        ),
        (
            bytes([2]) + _uvarint(1) + _uvarint(1) + _uvarint(0) + _uvarint(2),
            1,
            "Uint8Array",
            "exceeds",
        ),
        (
            bytes([2]) + _uvarint(1) + _uvarint(1) + _uvarint(65_536) + _uvarint(1),
            1,
            "Uint8Array",
            "declared code type",
        ),
        (
            _rle_codes(2, 2, [(1, 1), (1, 1)]),
            2,
            "Uint8Array",
            "distinct",
        ),
    ],
)
def test_user_defined_codecs_reject_non_exact_payloads(
    payload: bytes,
    expected_length: int,
    expected_codes_type: str,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        decode_user_defined_codes(
            payload,
            expected_length=expected_length,
            expected_codes_type=expected_codes_type,
        )


@pytest.mark.parametrize(
    ("payload", "expected_dtype"),
    [
        (_rle_codes(2, 4, [(1, 3), (0, 1)]), np.dtype(np.uint8)),
        (_rle_codes(3, 4, [(256, 3), (1, 1)]), np.dtype(np.uint16)),
    ],
)
def test_user_defined_codecs_preserve_current_rle_dtype(
    payload: bytes,
    expected_dtype: np.dtype[Any],
) -> None:
    decoded = decode_user_defined_codes(
        payload,
        expected_length=4,
        expected_codes_type=(
            "Uint8Array"
            if expected_dtype == np.dtype(np.uint8)
            else "Uint16Array"
        ),
    )

    assert decoded.dtype == expected_dtype
    assert decoded.tolist() == (
        [1, 1, 1, 0]
        if expected_dtype == np.dtype(np.uint8)
        else [256, 256, 256, 1]
    )
