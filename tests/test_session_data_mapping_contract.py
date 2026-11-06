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
from pathlib import Path
from typing import Any

import _session_files as files
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid.anndata_session import apply_cellucid_session_to_anndata
from cellucid.anndata_session._apply import cell_order_digest

# The duplicate var names are the point of several cases here.
pytestmark = pytest.mark.filterwarnings("ignore:Variable names are not unique:UserWarning")

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
        "embeddings": {
            "enabled": [
                {"key": "X_pca", "columns": [3, 1]},
                {"key": "X_umap", "columns": [0, 1]},
            ],
            "defaults": {"1d": None, "2d": {"key": "X_pca", "columns": [3, 1]}, "3d": None},
        },
        "knn": "connectivities",
        "geneNames": "index",
        "expression": "X",
        "log1p": {},
        "highlyVariable": None,
        "latent": "X_pca",
        "vectorFields": [
            {"key": "velocity_umap_2d", "embedding": {"key": "X_umap", "columns": [0, 1]}}
        ],
    }
    mapping.update(changes)
    return mapping


def _embeddings(enabled: list[dict[str, Any]], default_2d: dict[str, Any]) -> dict[str, Any]:
    """An embeddings group enabling `enabled` with `default_2d` as the 2D default."""
    return {"enabled": enabled, "defaults": {"1d": None, "2d": default_2d, "3d": None}}


def _digest(embedding: dict[str, Any]) -> str:
    values = _adata().obsm.get(embedding["key"])
    if values is None or max(embedding["columns"]) >= values.shape[1]:
        return "0123456789abcdef"
    return cell_order_digest(np.asarray(values)[:, embedding["columns"]])


def _fingerprint(mapping: dict[str, Any], *, var_count: int = 3) -> dict[str, Any]:
    defaults = mapping["embeddings"]["defaults"]
    slot, embedding = next((slot, ref) for slot, ref in defaults.items() if ref is not None)
    dimension = int(slot[0])
    return {
        "datasetId": "exact",
        "name": None,
        "cellCount": 2,
        "varCount": var_count,
        "cellOrder": {
            "embedding": copy.deepcopy(embedding),
            "dimension": dimension,
            # The cells of `_adata()`, in their order: an apply checks the digest
            # (a placeholder where the mapping names what `_adata()` lacks).
            "digest": _digest(embedding),
        },
        "dataMapping": mapping,
    }


def _var_alias(source_key: str, *, key: str = "copied") -> dict[str, Any]:
    return files.user_copy("alias_1", key, "var", source_key)


class _Bundle:
    """A session whose head holds ``fields`` (user fields), saved on ``fingerprint``."""

    def __init__(self, fingerprint: dict[str, Any], fields: list[dict[str, Any]] | None = None):
        self.fingerprint = fingerprint
        self.fields = fields or []


@pytest.fixture(autouse=True)
def _session_directory(tmp_path: Path) -> None:
    global _DIRECTORY
    _DIRECTORY = tmp_path


_DIRECTORY = Path()


def _apply(bundle: _Bundle, adata: ad.AnnData) -> ad.AnnData:
    path = files.session_file(
        _DIRECTORY / "mapping.cellucid-session",
        fingerprint_record=bundle.fingerprint,
        trackers=[{"id": "t1", "name": "Main", "root": files.state(
            mapping=bundle.fingerprint.get("dataMapping"), field_overlays=files.overlays(user_fields=bundle.fields),
        )}],
    )
    return apply_cellucid_session_to_anndata(adata, path, dry_run=False).adata


def test_a_matching_anndata_mapping_passes_every_check() -> None:
    adata = _adata()
    mapping = _anndata_mapping()
    path = files.session_file(
        _DIRECTORY / "mapping.cellucid-session",
        fingerprint_record=_fingerprint(mapping),
        trackers=[{"id": "t1", "name": "Main", "root": files.state(mapping=mapping)}],
    )
    plan = apply_cellucid_session_to_anndata(adata, path)
    # The opening state holds nothing to write: a plan with no rows.
    assert (plan.step, plan.rows) == ("s0", ())
    assert plan.apply().applied is False


def test_a_session_saved_with_the_removed_export_kind_is_refused() -> None:
    mapping = _anndata_mapping()
    mapping["kind"] = "export"
    with pytest.raises(ValueError, match='kind must be "anndata"'):
        _apply(_Bundle(_fingerprint(mapping)), _adata())


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        (
            {
                "embeddings": _embeddings(
                    [{"key": "X_tsne", "columns": [0, 1]}, {"key": "X_umap", "columns": [0, 1]}],
                    {"key": "X_tsne", "columns": [0, 1]},
                )
            },
            r"embeddings\.enabled\[0\] uses obsm\['X_tsne'\], which this dataset does not have",
        ),
        (
            {
                "embeddings": _embeddings(
                    [{"key": "X_pca", "columns": [4, 1]}, {"key": "X_umap", "columns": [0, 1]}],
                    {"key": "X_pca", "columns": [4, 1]},
                )
            },
            r"embeddings\.enabled\[0\] uses obsm\['X_pca'\]: Embedding column 4 does not exist",
        ),
        (
            # An enabled representation that is no dimension's default is
            # checked like a default: a view may show it.
            {
                "embeddings": _embeddings(
                    [
                        {"key": "X_pca", "columns": [3, 1]},
                        {"key": "X_phate", "columns": [0, 1]},
                        {"key": "X_umap", "columns": [0, 1]},
                    ],
                    {"key": "X_pca", "columns": [3, 1]},
                )
            },
            r"embeddings\.enabled\[1\] uses obsm\['X_phate'\], which this dataset does not have",
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
    with pytest.raises(ValueError, match=f"saved on a dataset of {var_count + 1} genes?, and this object has {var_count}"):
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
    with pytest.raises(ValueError, match="copies the gene 'MATR3-1', which must exist exactly once"):
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
    # Label positions are derived from codes and coordinates whenever an
    # embedding is shown; a value carrying them is not one this version reads.
    field = files.user_category("cat_1", "status", ["A"], [0, 0])
    field["centroidsByDim"] = {}
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
    with pytest.raises(ValueError, match=f"missing: {missing}"):
        _apply(_Bundle(fingerprint), _adata())
