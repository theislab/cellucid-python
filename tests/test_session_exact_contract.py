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
) -> Path:
    manifest = {
        "createdAt": "2026-07-25T00:00:00.000Z",
        "dataSource": None,
        "datasetFingerprint": fingerprint
        or {
            "sourceType": "anndata",
            "datasetId": "exact",
            "cellCount": 3,
            "varCount": 2,
        },
        "summary": None,
        "chunks": [meta for meta, _payload in entries],
    }
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
        var=pd.DataFrame(
            index=pd.Index(["gene-a", "gene-b"], dtype=object)
        ),
    )


class _MemoryBundle:
    def __init__(
        self,
        *,
        chunks: dict[str, Any],
        fingerprint: dict[str, Any] | None = None,
        contributors: dict[str, str] | None = None,
    ):
        self._chunks = chunks
        self.dataset_fingerprint = fingerprint or {
            "sourceType": "anndata",
            "datasetId": "exact",
            "cellCount": 3,
            "varCount": 2,
        }
        self.manifest = {
            "createdAt": "2026-07-25T00:00:00.000Z",
            "datasetFingerprint": self.dataset_fingerprint,
            "chunks": [
                {
                    "id": chunk_id,
                    "contributorId": (
                        contributors[chunk_id]
                        if contributors is not None
                        else _application_contributor(chunk_id)
                    ),
                }
                for chunk_id in chunks
            ],
        }

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
    if chunk_id.startswith("user-defined/codes/"):
        return "user-defined-codes"
    raise AssertionError(f"Test fixture requires an explicit contributor for {chunk_id!r}")


def _field_overlays(fields: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "renames": None,
        "deletedFields": None,
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

    def forbidden_decode(_chunk_id: str) -> Any:
        raise AssertionError("nonmaterialized current chunks must not be decoded")

    bundle.decode_chunk = forbidden_decode  # type: ignore[method-assign]
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
                        "color": None,
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
                "activePageName": "Page 1",
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
        "activePageName": "Page 1",
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
                        "color": None,
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
                "activePageName": "Page 1",
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
    ("payload", "message"),
    [
        (_raw_u8_codes([0, 1]) + b"\x00", "trailing bytes"),
        (bytes([2]) + _uvarint(1) + _uvarint(1) + _uvarint(0) + _uvarint(0), "run"),
        (bytes([2]) + _uvarint(1) + _uvarint(1) + _uvarint(0) + _uvarint(2), "exceeds"),
        (
            bytes([2]) + _uvarint(1) + _uvarint(1) + _uvarint(65_536) + _uvarint(1),
            "Uint16",
        ),
    ],
)
def test_user_defined_codecs_reject_non_exact_payloads(
    payload: bytes,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        decode_user_defined_codes(payload)
