"""``obs['total_counts']`` and ``uns['cellucid']['totalCounts']``: each cell's total count.

The contract (Cellucid-ready AnnData specification, ``totalCounts``): a file
whose ``X`` holds counts records each cell's total, ``"over": "all genes"``
for the caller's column (``prepare(total_counts=...)``), ``"stored genes"``
for the sums of ``X``'s rows as written. cellucid-r's ``test-total-counts.R``
pins the same.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid import prepare, restore_original_order

N_CELLS = 40
N_GENES = 12


def _counts(seed: int = 0) -> sparse.csr_matrix:
    rng = np.random.default_rng(seed)
    dense = rng.poisson(1.5, size=(N_CELLS, N_GENES)).astype(np.float32)
    dense[3] = 0  # a cell with no counts: total 0
    return sparse.csr_matrix(dense)


def _adata(matrix: object = None, **obs: object) -> ad.AnnData:
    rng = np.random.default_rng(1)
    adata = ad.AnnData(
        X=_counts() if matrix is None else matrix,
        obs=pd.DataFrame(obs or None, index=[f"c{i}" for i in range(N_CELLS)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(N_GENES)]),
    )
    adata.obsm["X_umap"] = rng.normal(size=(N_CELLS, 2))
    return adata


def _write(tmp_path: Path, adata: object, name: str = "a.zarr", **kwargs: object) -> ad.AnnData:
    output = "h5ad" if name.endswith(".h5ad") else "zarr"
    path = prepare(tmp_path / name, format=output, data=adata, dataset_id="t", dataset_name="T", **kwargs)
    return restore_original_order(ad.read_h5ad(path) if output == "h5ad" else ad.read_zarr(path))


def _row_sums(matrix: object) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1), dtype=np.float64).ravel()


@pytest.mark.parametrize(
    "matrix",
    [
        _counts(),
        _counts().tocsc(),
        _counts().toarray(),
        _counts().toarray().astype(np.int32),
        _counts().astype(np.float64),
    ],
    ids=["csr", "csc", "dense", "int32", "float64"],
)
def test_counts_are_summed_over_the_stored_genes(tmp_path: Path, matrix: object) -> None:
    back = _write(tmp_path, _adata(matrix))
    assert back.uns["cellucid"]["totalCounts"] == {"field": "total_counts", "over": "stored genes"}
    assert back.obs["total_counts"].dtype == np.float64
    np.testing.assert_array_equal(back.obs["total_counts"].to_numpy(), _row_sums(_counts()))
    assert back.obs["total_counts"]["c3"] == 0


def test_the_sums_are_of_the_genes_stored_after_gene_identifiers(tmp_path: Path) -> None:
    back = _write(tmp_path, _adata(), gene_identifiers=["g0", "g4", "g7"])
    np.testing.assert_array_equal(back.obs["total_counts"].to_numpy(), _row_sums(_counts()[:, [0, 4, 7]]))
    np.testing.assert_array_equal(back.obs["total_counts"].to_numpy(), _row_sums(back.X))


def test_a_backed_matrix_is_summed_as_one_in_memory(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    _adata().write_h5ad(source)
    back = _write(tmp_path, ad.read_h5ad(source, backed="r"))
    np.testing.assert_array_equal(back.obs["total_counts"].to_numpy(), _row_sums(_counts()))


def test_a_matrix_that_holds_no_counts_records_no_totals(tmp_path: Path) -> None:
    logged = _counts()
    logged.data = np.log1p(logged.data)
    back = _write(tmp_path, _adata(logged))
    assert "totalCounts" not in back.uns["cellucid"]
    assert "total_counts" not in back.obs
    negative = _counts().toarray()
    negative[0, 0] = -1
    assert "totalCounts" not in _write(tmp_path, _adata(negative), "n.zarr").uns["cellucid"]


def test_the_callers_column_is_recorded_over_all_genes(tmp_path: Path) -> None:
    whole = _row_sums(_counts()) + 7  # genes left out before prepare
    back = _write(tmp_path, _adata(nCount_RNA=whole), "a.h5ad", total_counts="nCount_RNA")
    assert back.uns["cellucid"]["totalCounts"] == {"field": "nCount_RNA", "over": "all genes"}
    assert "total_counts" not in back.obs
    np.testing.assert_array_equal(back.obs["nCount_RNA"].to_numpy(), whole)
    # Read back and prepared again, it stays named (``cellucid_record``).
    again = _write(tmp_path, back, "again.zarr")
    assert again.uns["cellucid"]["totalCounts"] == {"field": "nCount_RNA", "over": "all genes"}


@pytest.mark.parametrize(
    ("obs", "named", "message"),
    [
        ({"tc": np.log1p(_row_sums(_counts()))}, "tc", "below the sum of the genes stored for"),
        ({"tc": np.where(np.arange(N_CELLS) == 5, np.nan, 100.0)}, "tc", "with no total of 0 or more"),
        ({"tc": ["a"] * N_CELLS}, "tc", "not a column of numbers"),
        ({}, "tc", "which obs does not hold"),
    ],
    ids=["log", "missing", "text", "absent"],
)
def test_a_column_that_is_not_the_totals_is_refused(tmp_path: Path, obs: dict, named: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _write(tmp_path, _adata(**obs), total_counts=named)
    assert not (tmp_path / "a.zarr").exists()


def test_a_named_total_needs_counts_and_a_name(tmp_path: Path) -> None:
    logged = _counts()
    logged.data = np.log1p(logged.data)
    with pytest.raises(ValueError, match="does not hold counts"):
        _write(tmp_path, _adata(logged, tc=np.full(N_CELLS, 100.0)), total_counts="tc")
    with pytest.raises(TypeError, match="names the obs column"):
        _write(tmp_path, _adata(), total_counts=["tc"])
    with pytest.raises(TypeError, match="names the obs column"):
        _write(tmp_path, _adata(), total_counts={"rna": "tc"})


def test_the_datas_own_total_counts_column_is_kept(tmp_path: Path) -> None:
    sums = _row_sums(_counts())
    same = _write(tmp_path, _adata(total_counts=sums.astype(np.float32)), "same.zarr")
    assert same.uns["cellucid"]["totalCounts"] == {"field": "total_counts", "over": "stored genes"}
    assert same.obs["total_counts"].dtype == np.float32
    other = _write(tmp_path, _adata(total_counts=sums + 3), "other.zarr")
    assert "totalCounts" not in other.uns["cellucid"]
    np.testing.assert_array_equal(other.obs["total_counts"].to_numpy(), sums + 3)
    cased = _write(tmp_path, _adata(Total_Counts=sums), "cased.zarr")
    assert "totalCounts" not in cased.uns["cellucid"] and "total_counts" not in cased.obs


def _mudata() -> md.MuData:
    rna = ad.AnnData(X=_counts(), obs=pd.DataFrame(index=[f"c{i}" for i in range(N_CELLS)]))
    rna.var_names = [f"g{i}" for i in range(N_GENES)]
    prot = ad.AnnData(X=np.abs(np.random.default_rng(2).normal(size=(30, 4))), obs=pd.DataFrame(index=[f"c{i}" for i in range(30)]))
    prot.var_names = [f"p{i}" for i in range(4)]
    mdata = md.MuData({"rna": rna, "prot": prot})
    mdata.obsm["X_umap"] = np.random.default_rng(3).normal(size=(N_CELLS, 2))
    return mdata


def test_a_mudata_records_each_counts_modalitys_totals(tmp_path: Path) -> None:
    path = prepare(tmp_path / "m.zarr", format="zarr", data=_mudata(), dataset_id="m", dataset_name="M")
    back = md.read_zarr(str(path))
    assert dict(back.uns["cellucid"]["totalCounts"]) == {"rna": {"field": "rna:total_counts", "over": "stored genes"}}
    # Each stored row's total is its own row's sum, and the cells' totals as given.
    rna = back.mod["rna"]
    np.testing.assert_array_equal(rna.obs["total_counts"].to_numpy(), _row_sums(rna.X))
    by_name = dict(zip(rna.obs_names, rna.obs["total_counts"], strict=True))
    np.testing.assert_array_equal([by_name[f"c{i}"] for i in range(N_CELLS)], _row_sums(_counts()))
    assert "total_counts" not in back.mod["prot"].obs


def test_a_mudata_names_a_root_or_a_modality_column(tmp_path: Path) -> None:
    mdata = _mudata()
    mdata.obs["nCount_RNA"] = _row_sums(_counts()) + 2
    path = prepare(
        tmp_path / "m.h5mu", format="h5mu", data=mdata, dataset_id="m", dataset_name="M",
        total_counts={"rna": "nCount_RNA"},
    )
    assert dict(md.read_h5mu(str(path)).uns["cellucid"]["totalCounts"]) == {"rna": {"field": "nCount_RNA", "over": "all genes"}}
    mdata = _mudata()
    mdata.mod["rna"].obs["n_counts"] = _row_sums(_counts()) + 1
    mdata.update()
    path = prepare(
        tmp_path / "own.zarr", format="zarr", data=mdata, dataset_id="m", dataset_name="M",
        total_counts={"rna": "rna:n_counts"},
    )
    assert dict(md.read_zarr(str(path)).uns["cellucid"]["totalCounts"]) == {"rna": {"field": "rna:n_counts", "over": "all genes"}}
    with pytest.raises(ValueError, match="Modality 'prot'.*does not hold counts"):
        prepare(tmp_path / "x.zarr", format="zarr", data=_mudata(), dataset_id="m", dataset_name="M",
                total_counts={"prot": "prot:x"})
    with pytest.raises(ValueError, match="not modalities written"):
        prepare(tmp_path / "y.zarr", format="zarr", data=_mudata(), dataset_id="m", dataset_name="M",
                total_counts={"atac": "x"})
    with pytest.raises(TypeError, match="With a MuData"):
        prepare(tmp_path / "z.zarr", format="zarr", data=_mudata(), dataset_id="m", dataset_name="M",
                total_counts="total_counts")


def test_a_column_over_the_sources_gene_panel_is_recorded_as_a_source_panel(tmp_path: Path) -> None:
    # The built-in atlases' totals cover their 8,192-gene integration panels: more genes than are
    # stored, not the whole transcriptome.
    panel = _row_sums(_counts()) + 5
    back = _write(tmp_path, _adata(tc=panel), "p.zarr", total_counts="tc", total_counts_over="source panel")
    assert back.uns["cellucid"]["totalCounts"] == {"field": "tc", "over": "source panel"}
    # Read back and prepared again, it keeps what it covers; a column named anew says its own.
    again = _write(tmp_path, back, "again.zarr")
    assert again.uns["cellucid"]["totalCounts"] == {"field": "tc", "over": "source panel"}
    renamed = _write(tmp_path, back, "renamed.zarr", total_counts="tc")
    assert renamed.uns["cellucid"]["totalCounts"] == {"field": "tc", "over": "all genes"}
    # Checked as "all genes" is.
    with pytest.raises(ValueError, match="below the sum of the genes stored"):
        _write(tmp_path, _adata(tc=_row_sums(_counts()) * 0.5), "short.zarr", total_counts="tc", total_counts_over="source panel")
    with pytest.raises(ValueError, match="no total_counts was given"):
        _write(tmp_path, _adata(), "alone.zarr", total_counts_over="source panel")
    with pytest.raises(ValueError, match="total_counts_over must be 'all genes'.*or 'source panel'"):
        _write(tmp_path, _adata(tc=panel), "stored.zarr", total_counts="tc", total_counts_over="stored genes")


def test_a_mudata_says_what_each_named_column_covers(tmp_path: Path) -> None:
    mdata = _mudata()
    mdata.obs["nCount_RNA"] = _row_sums(_counts()) + 2
    path = prepare(
        tmp_path / "m.zarr", format="zarr", data=mdata, dataset_id="m", dataset_name="M",
        total_counts={"rna": "nCount_RNA"}, total_counts_over={"rna": "source panel"},
    )
    back = md.read_zarr(str(path))
    assert dict(back.uns["cellucid"]["totalCounts"]) == {"rna": {"field": "nCount_RNA", "over": "source panel"}}
    again = prepare(tmp_path / "again.zarr", format="zarr", data=back)
    assert dict(md.read_zarr(str(again)).uns["cellucid"]["totalCounts"]) == {"rna": {"field": "nCount_RNA", "over": "source panel"}}
    with pytest.raises(ValueError, match="total_counts_over names \\['prot'\\]"):
        prepare(tmp_path / "x.zarr", format="zarr", data=mdata, dataset_id="m", dataset_name="M",
                total_counts={"rna": "nCount_RNA"}, total_counts_over={"prot": "source panel"})


def test_a_refused_total_is_quoted_as_a_number_not_a_numpy_repr(tmp_path: Path) -> None:
    totals = _row_sums(_counts()) + 1
    totals[5] = -2.0
    with pytest.raises(ValueError, match=r"\(first 'c5': -2\.0\)\.$"):
        _write(tmp_path, _adata(tc=totals), total_counts="tc")
