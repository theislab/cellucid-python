import gzip
import json
import struct
from pathlib import Path

import numpy as np

from cellucid.anndata_session import apply_cellucid_session_to_anndata
from cellucid.session_bundle import SESSION_BUNDLE_MAGIC, CellucidSessionBundle
from cellucid.session_codecs import decode_delta_uvarint, decode_uvarint, decode_user_defined_codes


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
    indices = [10, 2, 7, 7]
    encoded = _encode_delta_uvarint(indices)
    decoded = decode_delta_uvarint(encoded, max_count=10, max_index=100)
    assert decoded.tolist() == sorted(indices)


def test_decode_user_defined_codes_raw_u8():
    # enc=0, length=5, payload=5 bytes
    payload = bytes([0]) + _encode_uvarint(5) + bytes([0, 1, 2, 1, 0])
    decoded = decode_user_defined_codes(payload)
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
                    {"id": "highlight_1", "type": "lasso", "label": "Lasso (2 cells)", "cellCount": 2}
                ],
            }
        ],
        "activePageId": "page_1",
        "activePageName": "Page 1",
    }
    highlights_cells = _encode_delta_uvarint([1, 3])

    overlays = {
        "renames": None,
        "deletedFields": None,
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
            }
        ],
    }
    udf_codes = bytes([0]) + _encode_uvarint(n_obs) + bytes([0, 1, 2, 1, 0])

    chunks = [
        gzip.compress(json.dumps(highlights_meta).encode("utf-8")),
        gzip.compress(highlights_cells),
        gzip.compress(json.dumps(overlays).encode("utf-8")),
        gzip.compress(udf_codes),
    ]

    manifest = {
        "createdAt": "2025-01-01T00:00:00Z",
        "dataSource": None,
        "datasetFingerprint": {"sourceType": "jupyter", "datasetId": "test", "cellCount": n_obs, "varCount": n_vars},
        "summary": None,
        "chunks": [
            {
                "id": "highlights/meta",
                "contributorId": "highlights-meta",
                "priority": "eager",
                "kind": "json",
                "codec": "gzip",
                "label": "Highlight metadata",
                "datasetDependent": True,
                "storedBytes": len(chunks[0]),
                "uncompressedBytes": len(json.dumps(highlights_meta).encode("utf-8")),
            },
            {
                "id": "highlights/cells/highlight_1",
                "contributorId": "highlights-cells",
                "priority": "lazy",
                "kind": "binary",
                "codec": "gzip",
                "label": "Highlight cells",
                "datasetDependent": True,
                "storedBytes": len(chunks[1]),
                "uncompressedBytes": len(highlights_cells),
            },
            {
                "id": "core/field-overlays",
                "contributorId": "field-overlays",
                "priority": "eager",
                "kind": "json",
                "codec": "gzip",
                "label": "Field overlays",
                "datasetDependent": True,
                "storedBytes": len(chunks[2]),
                "uncompressedBytes": len(json.dumps(overlays).encode("utf-8")),
            },
            {
                "id": "user-defined/codes/udf_1",
                "contributorId": "user-defined-codes",
                "priority": "eager",
                "kind": "binary",
                "codec": "gzip",
                "label": "User-defined codes",
                "datasetDependent": True,
                "storedBytes": len(chunks[3]),
                "uncompressedBytes": len(udf_codes),
            },
        ],
    }

    bundle_path = tmp_path / "test.cellucid-session"
    _write_bundle(bundle_path, manifest, chunks)

    bundle = CellucidSessionBundle(bundle_path)
    out = apply_cellucid_session_to_anndata(bundle, adata, inplace=False)

    assert "cellucid_highlight__highlight_1" in out.obs.columns
    assert out.obs["cellucid_highlight__highlight_1"].tolist() == [False, True, False, True, False]
    assert "my_label" in out.obs.columns
    assert list(out.obs["my_label"].cat.categories) == ["A", "B", "C"]
    assert out.obs["my_label"].tolist() == ["A", "B", "C", "B", "A"]
