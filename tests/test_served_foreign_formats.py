"""Foreign formats served in priority order: MuData, Loom, 10x Cell Ranger.

A file Python reads natively but did not write -- a MuData ``.h5mu`` or MuData
Zarr store, a Loom file, a Cell Ranger ``.h5`` matrix or ``outs/`` folder --
is served by path as a Zarr store whose cells are in priority order
(``cellucid._served_store``): what the browser presents from the file itself
(its layout views, mirrored in ``cellucid._served_formats``), with
``obs['cellucid_original_position']`` naming each row's position in the
source and ``uns['cellucid']['priority']`` the seed, so a subsample is the
first rows of every array. These tests read every key of such a store back
(anndata, mudata) and check it against the source in priority order; the web
app's ``tests/served-priority-store.test.mjs`` opens the same stores and the
files themselves in the browser's reader and compares them.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid._anndata_files import classify_anndata_path
from cellucid._served_store import SERVED_PRIORITY_SEED, PriorityStore
from cellucid.anndata_server import AnnDataServer
from cellucid.prepare_data._priority import (
    ORIGINAL_POSITION_COLUMN,
    PREFIX_COUNTS_KEY,
    CellPriority,
    prefix_levels,
)

N_CELLS = 240


def _order(n: int = N_CELLS) -> np.ndarray:
    priority = CellPriority.of(n, SERVED_PRIORITY_SEED)
    assert priority is not None
    return priority.order


def _dump(store: PriorityStore, folder: Path) -> Path:
    """Every key of a served store, written as a Zarr folder (what the browser reads, key by key)."""
    for key in store.stored_keys():
        value = store.get(key)
        assert value is not None, key
        target = folder / key
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(value)
    return folder


# =============================================================================
# Sources
# =============================================================================


def write_loom(path: Path, n_cells: int = N_CELLS, n_genes: int = 30, *, seed: int = 1) -> dict[str, np.ndarray]:
    """A Loom 2 file as loompy writes one: genes x cells, fixed-length text with XML references."""
    rng = np.random.default_rng(seed)
    matrix = (rng.random((n_genes, n_cells)) * (rng.random((n_genes, n_cells)) < 0.4) * 10).astype(np.float32)
    spliced = (matrix * 2).astype(np.float32)
    cell_ids = np.asarray([f"cell&#233;{i}".encode("ascii") for i in range(n_cells)], dtype="S16")
    with h5py.File(path, "w") as file:
        file.attrs["LOOM_SPEC_VERSION"] = "2.0.1"
        file.attrs["title"] = "A loom"
        chunks = (min(8, n_genes), min(16, n_cells))
        file.create_dataset("matrix", data=matrix, chunks=chunks, compression="gzip")
        file.create_dataset("layers/spliced", data=spliced, chunks=chunks)
        col = file.create_group("col_attrs")
        col.create_dataset("CellID", data=cell_ids)
        col.create_dataset("ClusterName", data=np.asarray([f"c{i % 4}" for i in range(n_cells)], dtype="S4"))
        col.create_dataset("Score", data=rng.random(n_cells))
        col.create_dataset("_X", data=rng.random(n_cells).astype(np.float32))
        col.create_dataset("_Y", data=rng.random(n_cells).astype(np.float32))
        col.create_dataset("X_pca", data=rng.random((n_cells, 3)))
        row = file.create_group("row_attrs")
        row.create_dataset("Gene", data=np.asarray([f"G{j}" for j in range(n_genes)], dtype=object), dtype=h5py.string_dtype())
        row.create_dataset("Chromosome", data=np.asarray([f"chr{j % 3}" for j in range(n_genes)], dtype="S5"))
        edges = rng.integers(0, n_cells, size=(400, 2))
        graph = file.create_group("col_graphs/KNN")
        graph.create_dataset("a", data=edges[:, 0])
        graph.create_dataset("b", data=edges[:, 1])
        graph.create_dataset("w", data=rng.random(400))
    return {"matrix": matrix, "spliced": spliced}


def write_cell_ranger(path: Path, n_cells: int = N_CELLS, *, types: tuple[tuple[str, int], ...] = (("Gene Expression", 12), ("Antibody Capture", 4)), seed: int = 2) -> sparse.csr_matrix:
    """A Cell Ranger 3 matrix: CSC over features x barcodes, features of each type together."""
    n_features = sum(count for _, count in types)
    cells = sparse.random(n_cells, n_features, density=0.3, format="csr", random_state=seed, dtype=np.float64)
    cells.data = np.ceil(cells.data * 20)
    cells.sort_indices()
    by_feature = cells.T.tocsc()
    with h5py.File(path, "w") as file:
        file.attrs["filetype"] = "matrix"
        file.attrs["version"] = 2
        matrix = file.create_group("matrix")
        matrix.create_dataset("barcodes", data=np.asarray([f"AAAC{i:06d}-1" for i in range(n_cells)], dtype="S14"))
        matrix.create_dataset("data", data=cells.data.astype(np.int32), chunks=(64,), compression="gzip")
        matrix.create_dataset("indices", data=cells.indices.astype(np.int64), chunks=(64,), compression="gzip")
        matrix.create_dataset("indptr", data=cells.indptr.astype(np.int64))
        matrix.create_dataset("shape", data=np.asarray(by_feature.shape, dtype=np.int32))
        features = matrix.create_group("features")
        names = [f"{kind[:3]}{j}" for kind, count in types for j in range(count)]
        features.create_dataset("id", data=np.asarray([f"ID{j}" for j in range(n_features)], dtype="S8"))
        features.create_dataset("name", data=np.asarray(names, dtype="S8"))
        features.create_dataset("feature_type", data=np.asarray([kind for kind, count in types for _ in range(count)], dtype="S20"))
        features.create_dataset("genome", data=np.asarray(["GRCh38"] * n_features, dtype="S6"))
        features.create_dataset("_all_tag_keys", data=np.asarray([b"genome"]))
    return cells


def write_outs(folder: Path, n_cells: int = N_CELLS, *, types: tuple[tuple[str, int], ...] = (("Gene Expression", 12),)) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(5)
    folder.mkdir(parents=True)
    cells = write_cell_ranger(folder / "filtered_feature_bc_matrix.h5", n_cells, types=types)
    barcodes = [f"AAAC{i:06d}-1" for i in range(n_cells)]
    umap = rng.random((n_cells, 2))
    projection = folder / "analysis" / "umap" / "gene_expression_2_components" / "projection.csv"
    projection.parent.mkdir(parents=True)
    # Missing the last cell: NaN there.
    projection.write_text("Barcode,UMAP-1,UMAP-2\n" + "".join(f"{b},{float(u[0])!r},{float(u[1])!r}\n" for b, u in zip(barcodes[:-1], umap[:-1], strict=True)))
    clusters = rng.integers(1, 12, n_cells)
    clustering = folder / "analysis" / "clustering" / "gene_expression_graphclust" / "clusters.csv"
    clustering.parent.mkdir(parents=True)
    clustering.write_text("Barcode,Cluster\n" + "".join(f"{b},{c}\n" for b, c in zip(barcodes, clusters, strict=True)))
    return cells, umap, clusters


def mudata_object(n_cells: int = N_CELLS, *, seed: int = 4):
    mudata = pytest.importorskip("mudata")
    rng = np.random.default_rng(seed)
    names = [f"cell-{i}" for i in range(n_cells)]
    rna = ad.AnnData(
        X=sparse.random(n_cells, 9, density=0.4, format="csr", random_state=seed, dtype=np.float32),
        obs=pd.DataFrame({"leiden": pd.Categorical(rng.choice(["a", "b"], n_cells))}, index=names),
        var=pd.DataFrame(index=[f"g{j}" for j in range(9)]),
    )
    rna.obsm["X_umap"] = rng.random((n_cells, 2)).astype(np.float32)
    rna.layers["counts"] = sparse.csc_matrix(rna.X) * 3
    # The protein panel measured three cells in four.
    measured = [name for i, name in enumerate(names) if i % 4]
    prot = ad.AnnData(
        X=rng.random((len(measured), 4)).astype(np.float32),
        obs=pd.DataFrame(index=measured),
        var=pd.DataFrame(index=[f"p{j}" for j in range(4)]),
    )
    prot.obsm["X_pca"] = rng.random((len(measured), 3))
    mdata = mudata.MuData({"rna": rna, "prot": prot})
    mdata.obsm["X_wnn_umap"] = rng.random((n_cells, 2))
    mdata.obs["sample"] = pd.Categorical(rng.choice(["s1", "s2"], n_cells))
    return mdata


# =============================================================================
# MuData
# =============================================================================


@pytest.mark.parametrize("container", ["h5mu", "zarr"])
def test_a_mudata_file_is_served_in_priority_order(tmp_path: Path, container: str) -> None:
    mudata = pytest.importorskip("mudata")
    mdata = mudata_object()
    path = tmp_path / ("data.h5mu" if container == "h5mu" else "data.zarr")
    if container == "h5mu":
        mdata.write_h5mu(path)
    else:
        mdata.write_zarr(path)
    source = classify_anndata_path(path)
    store = PriorityStore.of_mudata_file(source)
    try:
        assert store.layout == "mudata"
        served = mudata.read_zarr(_dump(store, tmp_path / "served.zarr"))
    finally:
        store.close()
    order = _order()
    # The shared cells in priority order, their original positions named.
    assert list(served.obs_names) == list(mdata.obs_names[order])
    assert np.array_equal(served.obs[ORIGINAL_POSITION_COLUMN].to_numpy(), order)
    assert list(served.obs["sample"]) == list(mdata.obs["sample"].to_numpy()[order])
    np.testing.assert_array_equal(served.obsm["X_wnn_umap"], mdata.obsm["X_wnn_umap"][order])
    priority = served.uns["cellucid"]["priority"]
    assert int(priority["seed"]) == SERVED_PRIORITY_SEED
    levels = prefix_levels(N_CELLS)
    assert np.array_equal(np.asarray(priority["levels"]), levels)
    for name in ("rna", "prot"):
        modality = served.mod[name]
        source_modality = mdata.mod[name]
        # Each modality's rows follow its cells: a prefix of the cells is a prefix of it.
        rows = [source_modality.obs_names.get_loc(cell) for cell in served.obs_names if cell in source_modality.obs_names]
        assert list(modality.obs_names) == list(source_modality.obs_names[rows])
        assert np.array_equal(modality.obs[ORIGINAL_POSITION_COLUMN].to_numpy(), rows)
        x = modality.X.toarray() if sparse.issparse(modality.X) else modality.X
        expected = source_modality.X.toarray() if sparse.issparse(source_modality.X) else source_modality.X
        np.testing.assert_array_equal(x, expected[rows])
        placed = np.asarray(served.obsmap[name])
        measured = np.flatnonzero(placed)
        assert np.array_equal(placed[measured], np.arange(1, len(rows) + 1))
        assert list(served.obs_names[measured]) == list(modality.obs_names)
    np.testing.assert_array_equal(served.mod["rna"].layers["counts"].toarray(), mdata.mod["rna"].layers["counts"].toarray()[[mdata.mod["rna"].obs_names.get_loc(c) for c in served.obs_names]])
    # A gene-major layer of a modality is a prefix per gene: its X is cell-major (rows in order).
    assert PREFIX_COUNTS_KEY not in served.mod["prot"].varm


def test_a_mudata_file_in_priority_order_is_served_as_it_is(tmp_path: Path) -> None:
    pytest.importorskip("mudata")
    from cellucid import prepare

    path = tmp_path / "ready.h5mu"
    prepare(path, format="h5mu", data=mudata_object(), dataset_id="m", dataset_name="M")
    server = AnnDataServer(path, port=0, quiet=True)
    try:
        assert server._served.priority is None and server._served.kind == "h5mu"
    finally:
        server.stop()


# =============================================================================
# Loom
# =============================================================================


def test_a_loom_file_is_served_as_the_browser_presents_it_in_priority_order(tmp_path: Path) -> None:
    path = tmp_path / "data.loom"
    arrays = write_loom(path)
    store = PriorityStore.of_loom(classify_anndata_path(path))
    try:
        served = ad.read_zarr(_dump(store, tmp_path / "served.zarr"))
    finally:
        store.close()
    order = _order()
    # The obs index from CellID, its XML character references decoded.
    assert list(served.obs_names) == [f"celln{i}".replace("n", "é") for i in order]
    assert np.array_equal(served.obs[ORIGINAL_POSITION_COLUMN].to_numpy(), order)
    assert set(served.obs.columns) == {"ClusterName", "Score", ORIGINAL_POSITION_COLUMN}
    assert list(served.var_names) == [f"G{j}" for j in range(30)]
    assert list(served.var["Chromosome"]) == [f"chr{j % 3}" for j in range(30)]
    np.testing.assert_array_equal(served.X, arrays["matrix"].T[order])
    np.testing.assert_array_equal(served.layers["spliced"], arrays["spliced"].T[order])
    # _X/_Y stacked into one embedding; an n x 3 attribute as it is.
    assert set(served.obsm) == {"X_embedding", "X_pca"}
    with h5py.File(path) as file:
        expected = np.column_stack([file["col_attrs/_X"][()], file["col_attrs/_Y"][()]])
        graph = file["col_graphs/KNN"]
        knn = sparse.coo_matrix((graph["w"][()], (graph["a"][()], graph["b"][()])), shape=(N_CELLS, N_CELLS)).tocsr()
    np.testing.assert_array_equal(served.obsm["X_embedding"], expected[order])
    np.testing.assert_array_equal(served.obsp["KNN"].toarray(), knn.toarray()[np.ix_(order, order)])
    assert served.uns["title"] == "A loom"
    assert int(served.uns["cellucid"]["priority"]["seed"]) == SERVED_PRIORITY_SEED


def test_a_loom_file_without_an_embedding_is_refused_in_the_browser_s_words(tmp_path: Path) -> None:
    path = tmp_path / "bare.loom"
    with h5py.File(path, "w") as file:
        file.create_dataset("matrix", data=np.ones((4, 6), dtype=np.float32))
        file.create_dataset("col_attrs/CellID", data=np.asarray([f"c{i}" for i in range(6)], dtype="S4"))
        file.create_dataset("row_attrs/Gene", data=np.asarray([f"g{i}" for i in range(4)], dtype="S4"))
    with pytest.raises(ValueError, match="has no cell embedding"):
        PriorityStore.of_loom(classify_anndata_path(path))


def test_a_loom_matrix_reads_bounded_spans(tmp_path: Path) -> None:
    """A served chunk of the transposed matrix reads one block of the rows it needs, never the matrix."""
    from cellucid._served_formats import Transposed

    path = tmp_path / "data.loom"
    arrays = write_loom(path, n_cells=500, n_genes=7)
    with h5py.File(path) as file:
        presented = Transposed(file["matrix"])
        rows = np.asarray([400, 3, 250, 3])
        np.testing.assert_array_equal(presented[rows, slice(2, 5)], arrays["matrix"].T[rows, 2:5])
        np.testing.assert_array_equal(presented[10:20, slice(0, 7)], arrays["matrix"].T[10:20])


# =============================================================================
# 10x Cell Ranger
# =============================================================================


def test_an_outs_folder_is_served_with_its_analysis_in_priority_order(tmp_path: Path) -> None:
    folder = tmp_path / "outs"
    cells, umap, clusters = write_outs(folder)
    store = PriorityStore.of_cell_ranger(classify_anndata_path(folder))
    try:
        served = ad.read_zarr(_dump(store, tmp_path / "served.zarr"))
    finally:
        store.close()
    order = _order()
    assert list(served.obs_names) == [f"AAAC{i:06d}-1" for i in order]
    assert np.array_equal(served.obs[ORIGINAL_POSITION_COLUMN].to_numpy(), order)
    # scanpy's names: the var index from features/name, gene_ids, feature_types, the tags.
    assert list(served.var.columns) == ["gene_ids", "feature_types", "genome"]
    assert sparse.issparse(served.X) and served.X.format == "csr"
    np.testing.assert_array_equal(served.X.toarray(), cells.toarray()[order])
    expected = umap.copy()
    expected[-1] = np.nan
    np.testing.assert_array_equal(served.obsm["X_umap"], expected[order])
    labels = sorted({str(c) for c in clusters}, key=int)
    assert list(served.obs["gene_expression_graphclust"].cat.categories) == labels
    assert list(served.obs["gene_expression_graphclust"].astype(str)) == [str(c) for c in clusters[order]]


def test_several_feature_types_are_modalities_each_a_window_of_the_matrix(tmp_path: Path) -> None:
    mudata = pytest.importorskip("mudata")
    folder = tmp_path / "outs"
    cells, _, _ = write_outs(folder, types=(("Gene Expression", 12), ("Antibody Capture", 4)))
    store = PriorityStore.of_cell_ranger(classify_anndata_path(folder))
    try:
        assert store.layout == "mudata"
        contents = store.contents()
        served = mudata.read_zarr(_dump(store, tmp_path / "served.zarr"))
    finally:
        store.close()
    order = _order()
    assert list(served.mod) == ["Gene Expression", "Antibody Capture"]
    assert contents.modalities is not None and set(contents.modalities) == {"Gene Expression", "Antibody Capture"}
    dense = cells.toarray()[order]
    np.testing.assert_array_equal(served.mod["Gene Expression"].X.toarray(), dense[:, :12])
    np.testing.assert_array_equal(served.mod["Antibody Capture"].X.toarray(), dense[:, 12:])
    assert list(served.mod["Antibody Capture"].var_names) == [f"Ant{j}" for j in range(4)]
    for name in served.mod:
        assert np.array_equal(np.asarray(served.obsmap[name]), np.arange(1, N_CELLS + 1))


def test_a_matrix_alone_is_refused_naming_the_outs_folder(tmp_path: Path) -> None:
    path = tmp_path / "filtered_feature_bc_matrix.h5"
    write_cell_ranger(path)
    with pytest.raises(ValueError, match="serve the whole outs/ folder"):
        PriorityStore.of_cell_ranger(classify_anndata_path(path))


def test_unsorted_cells_cannot_be_cut_into_feature_types(tmp_path: Path) -> None:
    folder = tmp_path / "outs"
    write_outs(folder, types=(("Gene Expression", 12), ("Antibody Capture", 4)))
    with h5py.File(folder / "filtered_feature_bc_matrix.h5", "r+") as file:
        indices = file["matrix/indices"][()]
        pointers = file["matrix/indptr"][()]
        row = next(r for r in range(N_CELLS) if pointers[r + 1] - pointers[r] >= 2)
        indices[pointers[row] : pointers[row + 1]] = indices[pointers[row] : pointers[row + 1]][::-1]
        file["matrix/indices"][...] = indices
    with pytest.raises(ValueError, match="unsorted within a cell"):
        PriorityStore.of_cell_ranger(classify_anndata_path(folder))


# =============================================================================
# The server
# =============================================================================


@pytest.mark.parametrize("kind", ["loom", "tenx-outs", "h5mu"])
def test_the_server_serves_a_foreign_file_through_the_store_in_priority_order(tmp_path: Path, kind: str) -> None:
    if kind == "loom":
        path = tmp_path / "data.loom"
        write_loom(path)
    elif kind == "tenx-outs":
        path = tmp_path / "outs"
        write_outs(path)
    else:
        pytest.importorskip("mudata")
        path = tmp_path / "data.h5mu"
        mudata_object().write_h5mu(path)
    server = AnnDataServer(path, port=0, quiet=True, dataset_id="d", dataset_name="D", serve_web_ui=False)
    try:
        store = server._served.priority
        assert store is not None
        expected_kind = "mudata-zarr" if kind == "h5mu" else "zarr"
        assert server._served.kind == expected_kind
        document = server._datasets_document()["datasets"][0]
        assert document["anndata"]["file"] == {"url": "_cellucid/anndata/d.zarr/", "kind": expected_kind, "layout": "folder"}
        # What the store holds is known before the browser opens it (the launch choices are checked).
        assert server.contents is not None and server.contents.n_obs == N_CELLS
        assert server.server_info["n_cells"] == N_CELLS
        status, _, body = server.answer_request("GET", "_cellucid/anndata/d.zarr/.zmetadata", None)
        assert status == 200
        assert f"obs/{ORIGINAL_POSITION_COLUMN}/.zarray" in json.loads(body)["metadata"]
        status, _, body = server.answer_request("GET", f"_cellucid/anndata/d.zarr/obs/{ORIGINAL_POSITION_COLUMN}/0", None)
        assert (status, body) == (200, store.get(f"obs/{ORIGINAL_POSITION_COLUMN}/0"))
    finally:
        server.stop()
    assert not store._finalizer.alive
