import json
import struct
from pathlib import Path

import numpy as np
import pytest

from cellucid.anndata_session import apply_cellucid_session_to_anndata
from cellucid.session_bundle import (
    NOT_A_SESSION_FILE_MESSAGE,
    SESSION_FORMAT_VERSION,
    CellucidSessionBundle,
    other_session_format_message,
    session_signature_refusal,
)
from cellucid.session_codecs import (
    decode_cell_set,
    decode_typed_json,
    decode_user_defined_codes,
    decode_uvarint,
    digest_cell_set,
)

FORMAT3 = Path(__file__).parent / "fixtures" / "session_format3"


def _mapping_fingerprint_keys() -> dict:
    """The cell-order record and data mapping every fingerprint carries."""
    return {
        "cellOrder": {
            "embedding": {"key": "X_umap", "columns": [0, 1]},
            "dimension": 2,
            "digest": "0123456789abcdef",
        },
        "dataMapping": {
            "kind": "anndata",
            "embeddings": {
                "enabled": [{"key": "X_umap", "columns": [0, 1]}],
                "defaults": {"1d": None, "2d": {"key": "X_umap", "columns": [0, 1]}, "3d": None},
            },
            "knn": None,
            "geneNames": "index",
            "expression": "X",
            "log1p": {},
            "highlyVariable": None,
            "latent": None,
            "vectorFields": [],
        },
    }


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


def test_decode_uvarint_roundtrip():
    for value in [0, 1, 2, 127, 128, 129, 16384, 1_000_000]:
        encoded = _encode_uvarint(value)
        decoded, next_off = decode_uvarint(encoded, 0)
        assert decoded == value
        assert next_off == len(encoded)


def _cell_set_sorted_delta(indices: list[int]) -> bytes:
    out = bytearray([0]) + _encode_uvarint(len(indices))
    previous = 0
    for position, index in enumerate(indices):
        out += _encode_uvarint(index if position == 0 else index - previous)
        previous = index
    return bytes(out)


def test_cell_set_encodings_decode_exactly():
    assert decode_cell_set(_cell_set_sorted_delta([2, 7, 10]), cell_count=11).tolist() == [2, 7, 10]
    bitset = bytes([1, 0b1000_0101, 0b0000_0010])
    assert decode_cell_set(bitset, cell_count=10).tolist() == [0, 2, 7, 9]
    runs = bytes([2]) + _encode_uvarint(2) + bytes([3, 2, 4, 3])
    assert decode_cell_set(runs, cell_count=20).tolist() == [3, 4, 9, 10, 11]
    with pytest.raises(ValueError, match="strictly increasing"):
        decode_cell_set(_cell_set_sorted_delta([2, 2]), cell_count=11)
    with pytest.raises(ValueError, match="exceeds the dataset cell count"):
        decode_cell_set(_cell_set_sorted_delta([2, 11]), cell_count=11)
    with pytest.raises(ValueError, match="bits are set past"):
        decode_cell_set(bytes([1, 0, 0b0000_0100]), cell_count=10)
    with pytest.raises(ValueError, match="trailing bytes"):
        decode_cell_set(_cell_set_sorted_delta([2]) + b"\x00", cell_count=11)
    with pytest.raises(ValueError, match="encoding 9"):
        decode_cell_set(bytes([9]), cell_count=11)


def test_cell_set_reference_encodings_check_the_category():
    reference = np.array([1, 4, 6, 9], dtype=np.uint32)
    digest = digest_cell_set(reference).to_bytes(4, "little")
    whole = bytes([4]) + _encode_uvarint(4) + digest
    assert decode_cell_set(whole, cell_count=10, reference=reference).tolist() == [1, 4, 6, 9]
    mask = bytes([3]) + _encode_uvarint(4) + digest + bytes([0b1010])
    assert decode_cell_set(mask, cell_count=10, reference=reference).tolist() == [4, 9]
    moved = np.array([1, 4, 6, 8], dtype=np.uint32)
    with pytest.raises(ValueError, match="no longer holds the same cells"):
        decode_cell_set(whole, cell_count=10, reference=moved)
    with pytest.raises(ValueError, match="its cells are required"):
        decode_cell_set(whole, cell_count=10)


def test_user_defined_codes_encodings():
    raw = bytes([0]) + _encode_uvarint(5) + bytes([0, 1, 2, 1, 0])
    decoded = decode_user_defined_codes(raw, expected_length=5, expected_codes_type="Uint8Array")
    assert decoded.dtype == np.uint8
    assert decoded.tolist() == [0, 1, 2, 1, 0]

    values = [0, 70_000, 4_294_967_295, 70_000]
    raw32 = bytes([0]) + _encode_uvarint(4) + np.array(values, dtype="<u4").tobytes()
    decoded = decode_user_defined_codes(raw32, expected_length=4, expected_codes_type="Uint32Array")
    assert decoded.dtype == np.uint32
    assert decoded.tolist() == values

    runs = bytes([1]) + _encode_uvarint(5) + _encode_uvarint(2)
    runs += _encode_uvarint(70_000) + _encode_uvarint(3)
    runs += _encode_uvarint(4_294_967_295) + _encode_uvarint(2)
    decoded = decode_user_defined_codes(runs, expected_length=5, expected_codes_type="Uint32Array")
    assert decoded.tolist() == [70_000] * 3 + [4_294_967_295] * 2

    # 2 bits per code, least significant bits first: 1, 2, 3, 0, then 2.
    packed = bytes([2]) + _encode_uvarint(5) + bytes([2, 0b0011_1001, 0b0000_0010])
    decoded = decode_user_defined_codes(packed, expected_length=5, expected_codes_type="Uint16Array")
    assert decoded.dtype == np.uint16
    assert decoded.tolist() == [1, 2, 3, 0, 2]
    with pytest.raises(ValueError, match="set bits past their length"):
        decode_user_defined_codes(
            bytes([2]) + _encode_uvarint(5) + bytes([2, 0, 0b0000_1000]),
            expected_length=5,
            expected_codes_type="Uint8Array",
        )
    with pytest.raises(ValueError, match="1, 2 or 4 bits"):
        decode_user_defined_codes(
            bytes([2]) + _encode_uvarint(1) + bytes([3, 0]),
            expected_length=1,
            expected_codes_type="Uint8Array",
        )
    with pytest.raises(ValueError, match="wrong length"):
        decode_user_defined_codes(raw32[:-1], expected_length=4, expected_codes_type="Uint32Array")


def test_typed_json_restores_arrays_and_refuses_dangling_references():
    text = json.dumps({"a": {"$typed": 0}, "b": [1, "x"]}).encode()
    array = np.array([1.5, -2.0], dtype="<f4").tobytes()
    payload = struct.pack("<I", len(text)) + text + bytes([5]) + _encode_uvarint(2) + array
    value = decode_typed_json(payload)
    assert value["b"] == [1, "x"]
    assert value["a"].dtype == np.float32
    assert value["a"].tolist() == [1.5, -2.0]
    with pytest.raises(ValueError, match="does not name"):
        decode_typed_json(struct.pack("<I", 2) + b"{}" + bytes([1]) + _encode_uvarint(0))
    with pytest.raises(ValueError, match="does not hold"):
        decode_typed_json(struct.pack("<I", len(text)) + text)


def test_web_app_format3_file_reads_exactly():
    """A file the web app's own encoders wrote (fixtures/session_format3)."""
    import anndata as ad
    import pandas as pd

    expected = json.loads((FORMAT3 / "expected.json").read_text())
    bundle = CellucidSessionBundle(FORMAT3 / "session.cellucid-session")
    assert bundle.dataset_fingerprint["name"] == "Format 3 fixture"
    assert bundle.dataset_source == {"kind": "file", "name": "format3.h5ad", "dataset": "format3"}
    assert bundle.required_fields == {"obs": ["ctype"], "var": []}
    windows = decode_typed_json(bundle.decode_chunk("analysis/windows"))
    assert windows["windows"][0]["result"]["log2fc"].tolist() == [1.5]
    assert windows["sizes"] == {"window": {"differential": {"width": 300, "height": 200}}, "expandedView": {}}
    assert windows["windows"][0]["settings"]["options"]["matrix"] == "mapped"
    docked = decode_typed_json(bundle.decode_chunk("analysis/docked"))
    assert docked["openMode"] == "differential"
    assert docked["modes"]["differential"]["result"]["pValue"].dtype == np.float64

    for name, digest in (("A", 0), ("B", 1)):
        start, stop = ((0, 100), (100, 900))[digest]
        assert digest_cell_set(np.arange(start, stop)) == expected["digests"][name]
    large = np.array([0, 1, 65535, 65536, 561946, 4294967295], dtype=np.uint32)
    assert digest_cell_set(large) == expected["digests"]["large"]

    n_obs = expected["cellCount"]
    adata = ad.AnnData(
        X=np.zeros((n_obs, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"ctype": pd.Categorical(expected["ctype"])},
            index=[f"cell_{i}" for i in range(n_obs)],
        ),
        var=pd.DataFrame(index=["gene_0"]),
        obsm={"X_umap": np.zeros((n_obs, 2), dtype=np.float32)},
    )
    out = apply_cellucid_session_to_anndata(bundle, adata, expected_dataset_id="format3")
    for group_id, cells in expected["groups"].items():
        column = out.obs[f"cellucid_highlight__{group_id}"].to_numpy()
        assert np.flatnonzero(column).tolist() == cells, group_id
    regions = ["R1", "R2", "rest"]
    assert out.obs["regions"].tolist() == [regions[code] for code in expected["regionsCodes"]]

    # The category the reference encodings were saved against has to hold the
    # same cells: one relabelled cell is refused, not silently reselected.
    moved = adata.copy()
    labels = list(expected["ctype"])
    labels[0] = "B"
    moved.obs["ctype"] = pd.Categorical(labels)
    with pytest.raises(ValueError, match="no longer holds the same cells"):
        apply_cellucid_session_to_anndata(bundle, moved, expected_dataset_id="format3")

    missing = adata.copy()
    del missing.obs["ctype"]
    with pytest.raises(ValueError, match="obs column 'ctype'"):
        apply_cellucid_session_to_anndata(bundle, missing, expected_dataset_id="format3")


def test_a_file_without_the_signature_is_refused_with_the_viewer_sentence(tmp_path):
    path = tmp_path / "unversioned.cellucid-session"
    manifest = json.dumps({"createdAt": "2025-01-01T00:00:00.000Z"}).encode()
    for content in (
        b"CELLUCID_SESSION\n" + struct.pack("<I", len(manifest)) + manifest,
        b"NOT_A_SESSION_FILE_AT_ALL\n",
    ):
        path.write_bytes(content)
        with pytest.raises(ValueError) as error:
            CellucidSessionBundle(path).list_chunk_ids()
        assert str(error.value) == NOT_A_SESSION_FILE_MESSAGE


def test_a_session_of_another_format_is_refused_naming_its_version(tmp_path):
    # The web app's sentences (bundle/format.js otherSessionFormatMessage),
    # word for word: a format-2 file's cell-order digest was taken over
    # display coordinates, which no later version can check.
    path = tmp_path / "other.cellucid-session"
    manifest = json.dumps({"createdAt": "2025-01-01T00:00:00.000Z"}).encode()
    assert SESSION_FORMAT_VERSION == 3
    for version, sentence in (
        (
            2,
            "This session was saved by an earlier version of Cellucid (session format 2): "
            "the check that it marks the same cells cannot be verified, so it was not opened. "
            "Re-create it in this version of Cellucid.",
        ),
        (
            4,
            "This session was saved by a newer version of Cellucid (session format 4), "
            "which this version cannot read. Open it in that version of Cellucid.",
        ),
    ):
        path.write_bytes(f"CELLUCID_SESSION/{version}\n".encode() + struct.pack("<I", len(manifest)) + manifest)
        with pytest.raises(ValueError) as error:
            CellucidSessionBundle(path).list_chunk_ids()
        assert str(error.value) == sentence == other_session_format_message(version)
    assert session_signature_refusal(b"CELLUCID_SESSION/x\n") == NOT_A_SESSION_FILE_MESSAGE


def test_apply_session_rejects_out_of_range_categorical_codes_without_mutation() -> None:
    import anndata as ad
    import pandas as pd

    class Bundle:
        def __init__(self) -> None:
            self.dataset_fingerprint = {
                "datasetId": "exact",
                "name": None,
                "cellCount": 3,
                "varCount": 1,
                **_mapping_fingerprint_keys(),
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
                "ui/figure-export": {"notDecoded": True},
                "analysis/windows": {"notDecoded": True},
                "analysis/docked": {"notDecoded": True},
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
                    "ui/figure-export",
                    "figure-export",
                    "eager",
                    "json",
                    "Figure export",
                    False,
                ),
                (
                    "analysis/windows",
                    "analysis-windows",
                    "eager",
                    "binary",
                    "Analysis windows",
                    True,
                ),
                (
                    "analysis/docked",
                    "analysis-docked",
                    "eager",
                    "binary",
                    "Docked analysis",
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
                "source": None,
                "requires": {"obs": [], "var": []},
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
        obsm={"X_umap": np.zeros((3, 2), dtype=np.float32)},
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


def _mock_bundle(*, n_cells: int, user_fields: list, codes: dict, groups: list, cells: dict):
    """A minimal current session: the generic chunks, overlays, codes and highlights given."""
    static = [
        ("core/field-overlays", "field-overlays", "json", "Field overlays", True),
        ("cinematic/camera", "cinematic-camera", "json", "Cinematic camera path", True),
        ("core/state", "core-state", "json", "Core state", True),
        ("ui/dockable-layout", "dockable-layout", "json", "Floating panels", False),
        ("ui/figure-export", "figure-export", "json", "Figure export", False),
        ("analysis/windows", "analysis-windows", "binary", "Analysis windows", True),
        ("analysis/docked", "analysis-docked", "binary", "Docked analysis", True),
        ("highlights/meta", "highlights-meta", "json", "Highlight metadata", True),
    ]
    chunks: dict = {
        "core/field-overlays": {
            "renames": {"fields": {}, "categories": {}},
            "deletedFields": {"deleted": [], "purged": []},
            "userDefinedFields": user_fields,
        },
        "highlights/meta": {
            "pages": [{"id": "page_1", "name": "Page 1", "color": "#2563eb", "highlightedGroups": groups}],
            "activePageId": "page_1",
        },
    }
    entries = []
    for chunk_id, contributor, kind, label, dependent in static[:1]:
        entries.append((chunk_id, contributor, "eager", kind, label, dependent))
    for field_id, (key, payload) in codes.items():
        chunks[f"user-defined/codes/{field_id}"] = payload
        entries.append((f"user-defined/codes/{field_id}", "user-defined-codes", "eager", "binary", f"User-defined codes: {key}", True))
    for chunk_id, contributor, kind, label, dependent in static[1:]:
        chunks.setdefault(chunk_id, {"notDecoded": True})
        entries.append((chunk_id, contributor, "eager", kind, label, dependent))
    for group in groups:
        chunk_id = f"highlights/cells/{group['id']}"
        chunks[chunk_id] = cells[group["id"]]
        entries.append((chunk_id, "highlights-cells", "lazy", "binary", f"Highlight cells: {group['label']}", True))

    class Bundle:
        dataset_fingerprint = {
            "datasetId": "exact",
            "name": None,
            "cellCount": n_cells,
            "varCount": 1,
            **_mapping_fingerprint_keys(),
        }

        def __init__(self) -> None:
            manifest_chunks = []
            for chunk_id, contributor, priority, kind, label, dependent in entries:
                payload = chunks[chunk_id]
                size = len(payload) if isinstance(payload, bytes) else len(json.dumps(payload).encode())
                manifest_chunks.append({
                    "id": chunk_id, "contributorId": contributor, "priority": priority, "kind": kind,
                    "codec": "gzip", "label": label, "datasetDependent": dependent,
                    "storedBytes": size, "uncompressedBytes": size,
                })
            self.manifest = {
                "createdAt": "2026-07-27T00:00:00.000Z",
                "datasetFingerprint": self.dataset_fingerprint,
                "source": None,
                "requires": {"obs": [], "var": []},
                "chunks": manifest_chunks,
            }

        def list_chunk_ids(self) -> list[str]:
            return [entry[0] for entry in entries]

        def decode_chunk(self, chunk_id: str):
            return chunks[chunk_id]

    return Bundle()


def test_a_continuous_alias_of_a_column_with_missing_values_applies_as_the_viewer_copied_it() -> None:
    """NaN is a float obs column's missing value; the viewer copies it as it is."""
    import anndata as ad
    import pandas as pd

    adata = ad.AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame({"score": [1.5, np.nan, -2.0]}, index=["a", "b", "c"]),
        var=pd.DataFrame(index=["gene"]),
        obsm={"X_umap": np.zeros((3, 2), dtype=np.float32)},
    )
    field = {
        "id": "field-1", "source": "obs", "kind": "continuous", "key": "score copy",
        "isDeleted": False, "isPurged": False, "createdAt": 1,
        "sourceField": {"kind": "continuous-obs", "sourceKey": "score", "sourceIndex": 0},
        "operation": {"type": "copy-field", "source": "obs", "sourceKey": "score"},
    }
    out = apply_cellucid_session_to_anndata(
        _mock_bundle(n_cells=3, user_fields=[field], codes={}, groups=[], cells={}),
        adata, expected_dataset_id="exact",
    )
    np.testing.assert_array_equal(out.obs["score copy"].to_numpy(), [1.5, np.nan, -2.0])


def test_the_recorded_session_leaves_out_where_the_viewer_opened_the_dataset() -> None:
    """A link, a server address or a file name never travels in ``uns`` with the object."""
    import anndata as ad
    import pandas as pd

    adata = ad.AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame(index=["a", "b", "c"]),
        var=pd.DataFrame(index=["gene"]),
        obsm={"X_umap": np.zeros((3, 2), dtype=np.float32)},
    )
    bundle = _mock_bundle(n_cells=3, user_fields=[], codes={}, groups=[], cells={})
    bundle.manifest["source"] = {
        "kind": "link",
        "locator": "https://lab.example/private/atlas.h5ad?versionId=3",
        "dataset": "exact",
    }
    out = apply_cellucid_session_to_anndata(bundle, adata, expected_dataset_id="exact")
    stored = out.uns["cellucid"]["session"]["manifest"]
    assert "source" not in stored
    assert set(stored) == set(bundle.manifest) - {"source"}
    assert bundle.manifest["source"]["kind"] == "link", "the bundle itself is unchanged"


def test_field_names_with_a_colon_apply_as_the_viewer_allows_them() -> None:
    """The viewer allows ':' in field names (ATAC peaks, MuData modality columns)."""
    import anndata as ad
    import pandas as pd

    adata = ad.AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame({"rna:cluster": pd.Categorical(["x", "y", "x"])}, index=["a", "b", "c"]),
        var=pd.DataFrame(index=["gene"]),
        obsm={"X_umap": np.zeros((3, 2), dtype=np.float32)},
    )
    field = {
        "id": "field-1", "source": "obs", "kind": "category", "key": "chr1:100-200",
        "categories": ["A", "B"], "isDeleted": False, "isPurged": False, "codesLength": 3,
        "codesType": "Uint8Array", "sourceField": None, "operation": None, "sourcePages": [],
        "overlapStrategy": "first", "overlapLabel": None, "intersectionLabels": None,
        "uncoveredLabel": None, "createdAt": 1,
    }
    group = {
        "id": "highlight_1", "type": "category", "label": "x", "enabled": True, "cellCount": 2,
        "fieldKey": "rna:cluster", "fieldSource": "obs", "categoryName": "x",
    }
    out = apply_cellucid_session_to_anndata(
        _mock_bundle(
            n_cells=3,
            user_fields=[field],
            codes={"field-1": ("chr1:100-200", bytes([0]) + _encode_uvarint(3) + bytes([0, 1, 0]))},
            groups=[group],
            cells={"highlight_1": _cell_set_sorted_delta([0, 2])},
        ),
        adata, expected_dataset_id="exact",
    )
    assert out.obs["chr1:100-200"].tolist() == ["A", "B", "A"]
    assert out.obs["cellucid_highlight__highlight_1"].tolist() == [True, False, True]
