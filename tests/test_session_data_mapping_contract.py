"""A session is applied to an AnnData only through the data mapping it was saved with.

``datasetFingerprint.dataMapping`` names the embeddings, gene names and -- for
an AnnData source -- the neighbour graph, expression matrix and latent space a
session was saved against. Before anything is written, every key it names must
exist on the AnnData, ``varCount`` must equal the size of the gene axis the
mapping reads, and a gene a user-defined field copies is found by the name the
viewer showed: the mapping's gene-name set, made unique with the shared rule,
with values from the mapping's expression matrix.
"""

from __future__ import annotations

import copy
import json
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid.anndata_session import apply_cellucid_session_to_anndata

# The duplicate var names are the point of several cases here.
pytestmark = pytest.mark.filterwarnings("ignore:Variable names are not unique:UserWarning")

_PROFILES = {
    "core/field-overlays": ("field-overlays", "json", "Field overlays", True),
    "cinematic/camera": ("cinematic-camera", "json", "Cinematic camera path", True),
    "core/state": ("core-state", "json", "Core state", True),
    "ui/dockable-layout": ("dockable-layout", "json", "Floating panels", False),
    "analysis/windows": ("analysis-windows", "json", "Analysis windows", True),
    "highlights/meta": ("highlights-meta", "json", "Highlight metadata", True),
    "analysis/cache-inventory": ("analysis-artifacts", "json", "Analysis cache inventory", True),
}


class _Bundle:
    """An in-memory bundle with the current chunk inventory and given overlays."""

    def __init__(self, fingerprint: dict[str, Any], fields: list[dict[str, Any]] | None = None):
        self.dataset_fingerprint = fingerprint
        self._chunks: dict[str, Any] = {
            "core/field-overlays": {
                "renames": {"fields": {}, "categories": {}},
                "deletedFields": {"deleted": [], "purged": []},
                "userDefinedFields": fields or [],
            },
            "cinematic/camera": {"notDecoded": True},
            "core/state": {"notDecoded": True},
            "ui/dockable-layout": {"notDecoded": True},
            "analysis/windows": {"notDecoded": True},
            "highlights/meta": {
                "pages": [{"id": "page_1", "name": "Page 1", "color": "#2563eb", "highlightedGroups": []}],
                "activePageId": "page_1",
            },
            "analysis/cache-inventory": {"artifactIds": []},
        }
        chunks = []
        for chunk_id, (contributor_id, kind, label, dataset_dependent) in _PROFILES.items():
            payload = json.dumps(self._chunks[chunk_id], separators=(",", ":")).encode("utf-8")
            chunks.append(
                {
                    "id": chunk_id,
                    "contributorId": contributor_id,
                    "priority": "eager",
                    "kind": kind,
                    "codec": "gzip",
                    "label": label,
                    "datasetDependent": dataset_dependent,
                    "storedBytes": len(payload),
                    "uncompressedBytes": len(payload),
                }
            )
        self.manifest = {
            "createdAt": "2026-10-18T00:00:00.000Z",
            "datasetFingerprint": fingerprint,
            "chunks": chunks,
        }

    def list_chunk_ids(self) -> list[str]:
        return list(self._chunks)

    def decode_chunk(self, chunk_id: str) -> Any:
        return self._chunks[chunk_id]


def _adata() -> ad.AnnData:
    adata = ad.AnnData(
        X=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["c0", "c1"]),
        var=pd.DataFrame(
            {"gene_symbols": ["CD4", None, "CD8A"]},
            index=pd.Index(["MATR3", "MATR3", "G3"], dtype=object),
        ),
    )
    adata.obsm["X_pca"] = np.array([[0.0, 1.0, 2.0, 3.0], [4.0, 5.0, 6.0, 7.0]])
    adata.obsm["X_umap"] = np.array([[0.0, 1.0], [1.0, 0.0]])
    adata.obsm["velocity_umap_2d"] = np.array([[0.5, 0.5], [0.5, 0.5]])
    adata.obsp["connectivities"] = sparse.csr_matrix(np.array([[0.0, 1.0], [1.0, 0.0]]))
    adata.layers["log1p"] = np.array([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]], dtype=np.float32)
    raw = ad.AnnData(
        X=np.array([[7.0, 8.0, 9.0, 10.0], [11.0, 12.0, 13.0, 14.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["c0", "c1"]),
        var=pd.DataFrame({"symbol": ["A", "B", "", "D"]}, index=["r0", "r1", "r2", "r3"]),
    )
    adata.raw = raw
    return adata


def _anndata_mapping(**changes: Any) -> dict[str, Any]:
    mapping: dict[str, Any] = {
        "kind": "anndata",
        "embeddings": {"1d": None, "2d": {"key": "X_pca", "columns": [3, 1]}, "3d": None},
        "knn": "connectivities",
        "geneNames": "index",
        "expression": "X",
        "latent": "X_pca",
        "vectorFields": [
            {"key": "velocity_umap_2d", "embedding": {"key": "X_umap", "columns": [0, 1]}}
        ],
    }
    mapping.update(changes)
    return mapping


def _fingerprint(mapping: dict[str, Any], *, var_count: int = 3) -> dict[str, Any]:
    embedding = next(ref for ref in mapping["embeddings"].values() if ref is not None)
    dimension = 1 if mapping["embeddings"]["1d"] is embedding else 2
    if mapping["embeddings"]["3d"] is embedding:
        dimension = 3
    return {
        "sourceType": "anndata",
        "datasetId": "exact",
        "cellCount": 2,
        "varCount": var_count,
        "cellOrder": {
            "embedding": copy.deepcopy(embedding),
            "dimension": dimension,
            "digest": "0123456789abcdef",
        },
        "dataMapping": mapping,
    }


def _var_alias(source_key: str, *, key: str = "copied") -> dict[str, Any]:
    return {
        "id": "alias_1",
        "source": "var",
        "kind": "continuous",
        "key": key,
        "isDeleted": False,
        "isPurged": False,
        "sourceField": {"sourceKey": source_key, "kind": "continuous"},
        "operation": {"type": "copy-field", "source": "var", "sourceKey": source_key},
        "createdAt": 1,
    }


def _apply(bundle: _Bundle, adata: ad.AnnData) -> ad.AnnData:
    return apply_cellucid_session_to_anndata(
        bundle,
        adata,
        expected_dataset_id="exact",
        add_highlights=False,
    )


def test_a_matching_anndata_mapping_is_applied_and_recorded() -> None:
    adata = _adata()
    mapping = _anndata_mapping()
    output = _apply(_Bundle(_fingerprint(mapping)), adata)
    recorded = output.uns["cellucid"]["session"]["dataset_fingerprint"]
    assert recorded["dataMapping"] == mapping
    assert recorded["cellOrder"]["embedding"] == {"key": "X_pca", "columns": [3, 1]}


def test_an_export_mapping_checks_only_its_gene_name_column() -> None:
    mapping = {
        "kind": "export",
        "embeddings": {"1d": None, "2d": {"key": "exported_id", "columns": None}, "3d": None},
        "geneNames": "column:gene_symbols",
        "vectorFields": [],
    }
    _apply(_Bundle(_fingerprint(mapping)), _adata())

    absent = copy.deepcopy(mapping)
    absent["geneNames"] = "column:not_there"
    with pytest.raises(ValueError, match="geneNames uses var column 'not_there'"):
        _apply(_Bundle(_fingerprint(absent)), _adata())


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        (
            {"embeddings": {"1d": None, "2d": {"key": "X_tsne", "columns": [0, 1]}, "3d": None}},
            r"embeddings\.2d uses obsm\['X_tsne'\], which this dataset does not have",
        ),
        (
            {"embeddings": {"1d": None, "2d": {"key": "X_pca", "columns": [4, 1]}, "3d": None}},
            r"embeddings\.2d uses obsm\['X_pca'\]: Embedding column 4 does not exist",
        ),
        ({"knn": "distances"}, r"knn uses obsp\['distances'\]"),
        ({"latent": "X_scvi"}, r"latent uses obsm\['X_scvi'\]"),
        ({"expression": "layer:counts"}, r"expression uses layers\['counts'\]"),
        ({"geneNames": "column:symbol"}, r"geneNames uses var column 'symbol'"),
        (
            {
                "vectorFields": [
                    {"key": "velocity_umap_3d", "embedding": {"key": "X_umap", "columns": [0, 1]}}
                ]
            },
            r"vector field obsm\['velocity_umap_3d'\] is not in this dataset",
        ),
    ],
)
def test_a_mapping_naming_what_the_anndata_lacks_is_refused_before_any_write(
    changes: dict[str, Any], message: str
) -> None:
    adata = _adata()
    before = adata.obs.copy(deep=True)
    with pytest.raises(ValueError, match=message):
        _apply(_Bundle(_fingerprint(_anndata_mapping(**changes))), adata)
    pd.testing.assert_frame_equal(adata.obs, before)
    assert "cellucid" not in adata.uns


def test_a_raw_mapping_needs_raw() -> None:
    adata = _adata()
    del adata.raw
    with pytest.raises(ValueError, match="expression uses raw, which this dataset does not have"):
        _apply(_Bundle(_fingerprint(_anndata_mapping(expression="raw"), var_count=4)), adata)


@pytest.mark.parametrize(
    ("expression", "var_count"),
    [("X", 3), ("layer:log1p", 3), ("raw", 4), (None, 0)],
)
def test_var_count_is_the_size_of_the_mapped_gene_axis(expression: str | None, var_count: int) -> None:
    mapping = _anndata_mapping(expression=expression)
    _apply(_Bundle(_fingerprint(mapping, var_count=var_count)), _adata())
    with pytest.raises(ValueError, match=f"varCount {var_count + 1} != "):
        _apply(_Bundle(_fingerprint(mapping, var_count=var_count + 1)), _adata())


def test_a_gene_is_found_by_its_made_unique_name() -> None:
    adata = _adata()
    output = _apply(
        _Bundle(_fingerprint(_anndata_mapping()), [_var_alias("MATR3-1")]),
        adata,
    )
    assert output.obs["copied"].tolist() == [2.0, 5.0]


def test_a_repeated_index_name_names_only_its_first_row() -> None:
    output = _apply(
        _Bundle(_fingerprint(_anndata_mapping()), [_var_alias("MATR3")]),
        _adata(),
    )
    assert output.obs["copied"].tolist() == [1.0, 4.0]


def test_a_gene_is_found_by_the_mapping_gene_name_column() -> None:
    # gene_symbols is CD4, <empty>, CD8A: the empty row is named by its index.
    mapping = _anndata_mapping(geneNames="column:gene_symbols")
    output = _apply(_Bundle(_fingerprint(mapping), [_var_alias("CD8A")]), _adata())
    assert output.obs["copied"].tolist() == [3.0, 6.0]
    output = _apply(_Bundle(_fingerprint(mapping), [_var_alias("MATR3")]), _adata())
    assert output.obs["copied"].tolist() == [2.0, 5.0]
    with pytest.raises(ValueError, match="source gene 'MATR3-1' must exist exactly once"):
        _apply(_Bundle(_fingerprint(mapping), [_var_alias("MATR3-1")]), _adata())


def test_values_are_copied_from_the_mapped_layer() -> None:
    mapping = _anndata_mapping(expression="layer:log1p")
    output = _apply(_Bundle(_fingerprint(mapping), [_var_alias("G3")]), _adata())
    assert output.obs["copied"].tolist() == [30.0, 60.0]


def test_values_and_names_are_read_from_raw() -> None:
    mapping = _anndata_mapping(expression="raw", geneNames="column:symbol")
    output = _apply(
        _Bundle(_fingerprint(mapping, var_count=4), [_var_alias("r2")]),
        _adata(),
    )
    assert output.obs["copied"].tolist() == [9.0, 13.0]


def test_a_mapping_without_expression_copies_no_gene() -> None:
    mapping = _anndata_mapping(expression=None)
    with pytest.raises(ValueError, match="cannot copy gene values: The data mapping reads no expression"):
        _apply(_Bundle(_fingerprint(mapping, var_count=0), [_var_alias("G3")]), _adata())


def test_a_user_defined_category_carrying_centroids_is_refused() -> None:
    field = {
        "id": "cat_1",
        "source": "obs",
        "kind": "category",
        "key": "status",
        "categories": ["A"],
        "isDeleted": False,
        "isPurged": False,
        "codesLength": 2,
        "codesType": "Uint8Array",
        "centroidsByDim": {},
        "sourceField": None,
        "operation": None,
        "sourcePages": [],
        "overlapStrategy": "first",
        "overlapLabel": None,
        "intersectionLabels": None,
        "uncoveredLabel": None,
        "createdAt": 1,
    }
    with pytest.raises(ValueError, match="unknown centroidsByDim"):
        _apply(_Bundle(_fingerprint(_anndata_mapping()), [field]), _adata())
    field.pop("centroidsByDim")
    field["normalizedDims"] = []
    with pytest.raises(ValueError, match="unknown normalizedDims"):
        _apply(_Bundle(_fingerprint(_anndata_mapping()), [field]), _adata())


@pytest.mark.parametrize("missing", ["cellOrder", "dataMapping"])
def test_the_apply_path_requires_both_mapping_keys(missing: str) -> None:
    fingerprint = _fingerprint(_anndata_mapping())
    del fingerprint[missing]
    with pytest.raises(ValueError, match=f"missing {missing}"):
        _apply(_Bundle(fingerprint), _adata())
