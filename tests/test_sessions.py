import gzip
import json
import struct
from pathlib import Path

import numpy as np
import pytest

from cellucid.anndata_session import apply_cellucid_session_to_anndata
from cellucid.session_bundle import SESSION_BUNDLE_MAGIC, CellucidSessionBundle
from cellucid.session_codecs import decode_delta_uvarint, decode_user_defined_codes, decode_uvarint


def _encode_uvarint(value: int) -> bytes:
    out = bytearray()
    v = int(value)
    while True:
        b = v & 0x7F
        v >>= 7
        if v:
            out.append(b | 0x80)
        else:
            out.append(b)
            return bytes(out)


def _encode_delta_uvarint(indices: list[int]) -> bytes:
    sorted_idx = sorted(int(i) for i in indices)
    out = bytearray()
    out += _encode_uvarint(len(sorted_idx))
    prev = 0
    for i, idx in enumerate(sorted_idx):
        delta = idx if i == 0 else idx - prev
        out += _encode_uvarint(delta)
        prev = idx
    return bytes(out)


def _write_bundle(path: Path, manifest: dict, chunks: list[bytes]) -> None:
    manifest_bytes = json.dumps(manifest).encode("utf-8")
    with path.open("wb") as f:
        f.write(SESSION_BUNDLE_MAGIC)
        f.write(struct.pack("<I", len(manifest_bytes)))
        f.write(manifest_bytes)
        for chunk in chunks:
            f.write(struct.pack("<I", len(chunk)))
            f.write(chunk)


def test_decode_uvarint_roundtrip():
    for value in [0, 1, 2, 127, 128, 129, 16384, 1_000_000]:
        encoded = _encode_uvarint(value)
        decoded, next_off = decode_uvarint(encoded, 0)
        assert decoded == value
        assert next_off == len(encoded)


def test_decode_delta_uvarint():
    indices = [10, 2, 7]
    encoded = _encode_delta_uvarint(indices)
    decoded = decode_delta_uvarint(encoded, max_count=10, max_index=100)
    assert decoded.tolist() == sorted(indices)


def test_decode_user_defined_codes_raw_u8():
    # enc=0, length=5, payload=5 bytes
    payload = bytes([0]) + _encode_uvarint(5) + bytes([0, 1, 2, 1, 0])
    decoded = decode_user_defined_codes(
        payload,
        expected_length=5,
        expected_codes_type="Uint8Array",
    )
    assert decoded.dtype == np.uint8
    assert decoded.tolist() == [0, 1, 2, 1, 0]


def test_apply_session_bundle_to_anndata(tmp_path):
    import anndata as ad
    import pandas as pd

    n_obs = 5
    n_vars = 3
    adata = ad.AnnData(
        X=np.zeros((n_obs, n_vars), dtype=np.float32),
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)]),
    )

    highlights_meta = {
        "pages": [
            {
                "id": "page_1",
                "name": "Page 1",
                "color": "#ff0000",
                "highlightedGroups": [
                    {
                        "id": "highlight_1",
                        "type": "lasso",
                        "label": "Lasso (2 cells)",
                        "enabled": True,
                        "cellCount": 2,
                    }
                ],
            }
        ],
        "activePageId": "page_1",
    }
    highlights_cells = _encode_delta_uvarint([1, 3])

    overlays = {
        "renames": {"fields": {}, "categories": {}},
        "deletedFields": {"deleted": [], "purged": []},
        "userDefinedFields": [
            {
                "id": "udf_1",
                "source": "obs",
                "kind": "category",
                "key": "my_label",
                "categories": ["A", "B", "C"],
                "codesLength": n_obs,
                "codesType": "Uint8Array",
                "isDeleted": False,
                "isPurged": False,
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
        ],
    }
    udf_codes = bytes([0]) + _encode_uvarint(n_obs) + bytes([0, 1, 2, 1, 0])

    camera = {"keyframes": [], "loop": False, "autoplay": False}
    core_state = {"notMaterializedByPython": True}
    dockable_layout = {"notMaterializedByPython": True}
    analysis_windows = {"notMaterializedByPython": True}
    analysis_inventory = {"artifactIds": []}
    raw_payloads = [
        json.dumps(overlays, separators=(",", ":")).encode("utf-8"),
        udf_codes,
        json.dumps(camera, separators=(",", ":")).encode("utf-8"),
        json.dumps(core_state, separators=(",", ":")).encode("utf-8"),
        json.dumps(dockable_layout, separators=(",", ":")).encode("utf-8"),
        json.dumps(analysis_windows, separators=(",", ":")).encode("utf-8"),
        json.dumps(highlights_meta, separators=(",", ":")).encode("utf-8"),
        json.dumps(analysis_inventory, separators=(",", ":")).encode("utf-8"),
        highlights_cells,
    ]
    chunks = [gzip.compress(payload, mtime=0) for payload in raw_payloads]

    manifest = {
        "createdAt": "2025-01-01T00:00:00.000Z",
        "datasetFingerprint": {
            "sourceType": "jupyter",
            "datasetId": "test",
            "cellCount": n_obs,
            "varCount": n_vars,
        },
        "chunks": [
            {
                "id": "core/field-overlays",
                "contributorId": "field-overlays",
                "priority": "eager",
                "kind": "json",
                "codec": "gzip",
                "label": "Field overlays",
                "datasetDependent": True,
                "storedBytes": len(chunks[0]),
                "uncompressedBytes": len(raw_payloads[0]),
            },
            {
                "id": "user-defined/codes/udf_1",
                "contributorId": "user-defined-codes",
                "priority": "eager",
                "kind": "binary",
                "codec": "gzip",
                "label": "User-defined codes: my_label",
                "datasetDependent": True,
                "storedBytes": len(chunks[1]),
                "uncompressedBytes": len(raw_payloads[1]),
            },
            {
                "id": "cinematic/camera",
                "contributorId": "cinematic-camera",
                "priority": "eager",
                "kind": "json",
                "codec": "gzip",
                "label": "Cinematic camera path",
                "datasetDependent": True,
                "storedBytes": len(chunks[2]),
                "uncompressedBytes": len(raw_payloads[2]),
            },
            {
                "id": "core/state",
                "contributorId": "core-state",
                "priority": "eager",
                "kind": "json",
                "codec": "gzip",
                "label": "Core state",
                "datasetDependent": True,
                "storedBytes": len(chunks[3]),
                "uncompressedBytes": len(raw_payloads[3]),
            },
            {
                "id": "ui/dockable-layout",
                "contributorId": "dockable-layout",
                "priority": "eager",
                "kind": "json",
                "codec": "gzip",
                "label": "Floating panels",
                "datasetDependent": False,
                "storedBytes": len(chunks[4]),
                "uncompressedBytes": len(raw_payloads[4]),
            },
            {
                "id": "analysis/windows",
                "contributorId": "analysis-windows",
                "priority": "eager",
                "kind": "json",
                "codec": "gzip",
                "label": "Analysis windows",
                "datasetDependent": True,
                "storedBytes": len(chunks[5]),
                "uncompressedBytes": len(raw_payloads[5]),
            },
            {
                "id": "highlights/meta",
                "contributorId": "highlights-meta",
                "priority": "eager",
                "kind": "json",
                "codec": "gzip",
                "label": "Highlight metadata",
                "datasetDependent": True,
                "storedBytes": len(chunks[6]),
                "uncompressedBytes": len(raw_payloads[6]),
            },
            {
                "id": "analysis/cache-inventory",
                "contributorId": "analysis-artifacts",
                "priority": "eager",
                "kind": "json",
                "codec": "gzip",
                "label": "Analysis cache inventory",
                "datasetDependent": True,
                "storedBytes": len(chunks[7]),
                "uncompressedBytes": len(raw_payloads[7]),
            },
            {
                "id": "highlights/cells/highlight_1",
                "contributorId": "highlights-cells",
                "priority": "lazy",
                "kind": "binary",
                "codec": "gzip",
                "label": "Highlight cells: Lasso (2 cells)",
                "datasetDependent": True,
                "storedBytes": len(chunks[8]),
                "uncompressedBytes": len(raw_payloads[8]),
            },
        ],
    }

    bundle_path = tmp_path / "test.cellucid-session"
    _write_bundle(bundle_path, manifest, chunks)

    bundle = CellucidSessionBundle(bundle_path)
    out = apply_cellucid_session_to_anndata(
        bundle,
        adata,
        inplace=False,
        expected_dataset_id="test",
    )

    assert "cellucid_highlight__highlight_1" in out.obs.columns
    assert out.obs["cellucid_highlight__highlight_1"].tolist() == [False, True, False, True, False]
    assert "my_label" in out.obs.columns
    assert list(out.obs["my_label"].cat.categories) == ["A", "B", "C"]
    assert out.obs["my_label"].tolist() == ["A", "B", "C", "B", "A"]


def test_apply_session_rejects_out_of_range_categorical_codes_without_mutation() -> None:
    import anndata as ad
    import pandas as pd

    class Bundle:
        def __init__(self) -> None:
            self.dataset_fingerprint = {
                "sourceType": "anndata",
                "datasetId": "exact",
                "cellCount": 3,
                "varCount": 1,
            }
            self._chunks = {
                "core/field-overlays": {
                    "renames": {"fields": {}, "categories": {}},
                    "deletedFields": {"deleted": [], "purged": []},
                    "userDefinedFields": [
                        {
                            "id": "field-1",
                            "source": "obs",
                            "kind": "category",
                            "key": "status",
                            "categories": ["A", "B"],
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
                    ],
                },
                "user-defined/codes/field-1": (bytes([0]) + _encode_uvarint(3) + bytes([0, 2, 1])),
                "cinematic/camera": {"notDecoded": True},
                "core/state": {"notDecoded": True},
                "ui/dockable-layout": {"notDecoded": True},
                "analysis/windows": {"notDecoded": True},
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
                },
                "analysis/cache-inventory": {"artifactIds": []},
            }
            profiles = [
                (
                    "core/field-overlays",
                    "field-overlays",
                    "eager",
                    "json",
                    "Field overlays",
                    True,
                ),
                (
                    "user-defined/codes/field-1",
                    "user-defined-codes",
                    "eager",
                    "binary",
                    "User-defined codes: status",
                    True,
                ),
                (
                    "cinematic/camera",
                    "cinematic-camera",
                    "eager",
                    "json",
                    "Cinematic camera path",
                    True,
                ),
                ("core/state", "core-state", "eager", "json", "Core state", True),
                (
                    "ui/dockable-layout",
                    "dockable-layout",
                    "eager",
                    "json",
                    "Floating panels",
                    False,
                ),
                (
                    "analysis/windows",
                    "analysis-windows",
                    "eager",
                    "json",
                    "Analysis windows",
                    True,
                ),
                (
                    "highlights/meta",
                    "highlights-meta",
                    "eager",
                    "json",
                    "Highlight metadata",
                    True,
                ),
                (
                    "analysis/cache-inventory",
                    "analysis-artifacts",
                    "eager",
                    "json",
                    "Analysis cache inventory",
                    True,
                ),
            ]
            manifest_chunks = []
            for (
                chunk_id,
                contributor_id,
                priority,
                kind,
                label,
                dataset_dependent,
            ) in profiles:
                payload = self._chunks[chunk_id]
                payload_bytes = (
                    payload
                    if isinstance(payload, bytes)
                    else json.dumps(payload, separators=(",", ":")).encode("utf-8")
                )
                manifest_chunks.append(
                    {
                        "id": chunk_id,
                        "contributorId": contributor_id,
                        "priority": priority,
                        "kind": kind,
                        "codec": "gzip",
                        "label": label,
                        "datasetDependent": dataset_dependent,
                        "storedBytes": len(payload_bytes),
                        "uncompressedBytes": len(payload_bytes),
                    }
                )
            self.manifest = {
                "createdAt": "2026-07-27T00:00:00.000Z",
                "datasetFingerprint": self.dataset_fingerprint,
                "chunks": manifest_chunks,
            }

        def list_chunk_ids(self) -> list[str]:
            return list(self._chunks)

        def decode_chunk(self, chunk_id: str):
            return self._chunks[chunk_id]

    adata = ad.AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame(index=["cell-1", "cell-2", "cell-3"]),
        var=pd.DataFrame(index=["gene"]),
    )

    with pytest.raises(ValueError, match="field-1.*code 2.*2 categories"):
        apply_cellucid_session_to_anndata(
            Bundle(),
            adata,
            inplace=True,
            expected_dataset_id="exact",
            add_highlights=False,
            store_uns=False,
        )

    assert "status" not in adata.obs
