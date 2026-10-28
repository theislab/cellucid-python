"""The store a server presents for data the kernel holds, in priority order.

An AnnData object, and an AnnData file not in priority order, are served as
a Zarr store whose cells are in priority order (``cellucid._served_store``):
the layout ``prepare()`` writes, made chunk by chunk from the source. These
tests read every key of such a store and check, against direct computations:

- the rows are the source's in priority order (``obs['cellucid_original_position']``
  ranks ascending under the served seed), every element bit for bit the
  source's at those rows, read back by anndata;
- a gene-major ``X`` holds each gene's rows ascending, and
  ``varm['cellucid_prefix_counts']`` its values below every level;
- a ``k``-cell prefix of every per-cell array lies in its first chunks, so the
  bytes a subsample reads grow with ``k``, not with the dataset;
- the arrays equal those ``prepare()`` writes for the same object;
- a file in priority order is served as it is; the server publishes the store
  as a Zarr folder and answers its keys, the widget path alike.
"""

from __future__ import annotations

import http.client
import json
import shutil
import tempfile
import threading
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import zarr
from _mapped_files import mapped_files
from scipy import sparse
from zarr.core.buffer import default_buffer_prototype

from cellucid._anndata_files import classify_anndata_path
from cellucid._served_store import SERVED_PRIORITY_SEED, PriorityStore, file_order
from cellucid.anndata_server import AnnDataServer
from cellucid.prepare_data._priority import (
    ORIGINAL_POSITION_COLUMN,
    PREFIX_COUNTS_KEY,
    CellPriority,
    prefix_levels,
)
from cellucid.subsample import _ranks


def _adata(n: int = 600, m: int = 9, *, form: str = "csr", seed: int = 3) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    x = sparse.random(n, m, density=0.3, format="csr", dtype=np.float32, random_state=seed)
    # One gene stored at every cell, and one at none.
    x = sparse.lil_matrix(x)
    x[:, 0] = rng.random((n, 1)).astype(np.float32) + 1
    x[:, m - 1] = 0
    x = sparse.csr_matrix(x)
    x.eliminate_zeros()
    adata = ad.AnnData(
        X=x if form == "csr" else (sparse.csc_matrix(x) if form == "csc" else x.toarray()),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(rng.choice(["a", "b", "c"], n)),
                "score": rng.random(n).astype(np.float32),
                "name": [f"n{i}" for i in range(n)],
                "flag": pd.array(rng.choice([True, False, None], n), dtype="boolean"),
            },
            index=pd.Index([f"cell-{i}" for i in range(n)], dtype=object),
        ),
        var=pd.DataFrame(
            {"symbol": [f"S{j}" for j in range(m)]}, index=pd.Index([f"g{j}" for j in range(m)], dtype=object)
        ),
    )
    adata.layers["counts"] = sparse.csc_matrix(x) * 2
    adata.layers["dense"] = x.toarray()
    adata.obsm["X_umap"] = rng.random((n, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((n, 5))
    adata.obsm["sparse_feature"] = sparse.random(n, 4, density=0.2, format="csr", random_state=seed + 1)
    adata.obsm["table"] = pd.DataFrame({"u": rng.random(n)}, index=adata.obs_names)
    graph = sparse.random(n, n, density=0.01, format="csr", dtype=np.float32, random_state=seed + 2).tolil()
    graph.setdiag(0)
    graph = sparse.csr_matrix(graph)
    # Symmetric with no self-loops, as prepare() requires of a neighbour graph.
    adata.obsp["connectivities"] = sparse.csr_matrix(graph + graph.T)
    adata.uns["cellucid"] = {"id": "served", "name": "Served"}
    adata.uns["note"] = "kept"
    adata.raw = adata.copy()
    return adata


def _read(store: PriorityStore) -> tuple[ad.AnnData, dict[str, bytes]]:
    """Every key of a served store, and the AnnData anndata reads from them."""
    keys = {key: store.get(key) for key in store.stored_keys()}
    assert all(value is not None for value in keys.values())
    memory = {key: default_buffer_prototype().buffer.from_bytes(value) for key, value in keys.items()}
    group = zarr.open_group(zarr.storage.MemoryStore(store_dict=memory), mode="r", zarr_format=2)
    return ad.io.read_elem(group), keys  # type: ignore[return-value]


def _same(a: object, b: object) -> None:
    if sparse.issparse(a) or sparse.issparse(b):
        assert sparse.issparse(a) and sparse.issparse(b)
        assert a.shape == b.shape and (a != b).nnz == 0
    else:
        np.testing.assert_array_equal(np.asarray(a), np.asarray(b))


def _check(served: ad.AnnData, source: ad.AnnData) -> np.ndarray:
    """``served`` is ``source`` in priority order; returns the order."""
    n = source.n_obs
    positions = served.obs[ORIGINAL_POSITION_COLUMN].to_numpy()
    assert positions.dtype == np.uint32
    order = positions.astype(np.int64)
    assert np.array_equal(np.sort(order), np.arange(n))
    ranks = _ranks(order, SERVED_PRIORITY_SEED).astype(np.int64)
    assert (np.diff(ranks) > 0).all(), "rows ascend by rank"
    want = source[order].copy()
    assert list(served.obs_names) == list(want.obs_names)
    pd.testing.assert_frame_equal(served.obs.drop(columns=ORIGINAL_POSITION_COLUMN), want.obs, check_categorical=True)
    pd.testing.assert_frame_equal(served.var, want.var)
    _same(served.X, want.X)
    for key in want.layers:
        _same(served.layers[key], want.layers[key])
    for key in ("X_umap", "X_pca", "sparse_feature"):
        _same(served.obsm[key], want.obsm[key])
    pd.testing.assert_frame_equal(served.obsm["table"], want.obsm["table"])
    _same(served.obsp["connectivities"], source.obsp["connectivities"][order][:, order])
    _same(served.raw.X, want.raw.X)
    pd.testing.assert_frame_equal(served.raw.var, want.raw.var)
    priority = served.uns["cellucid"]["priority"]
    assert (int(priority["seed"]), int(priority["specVersion"])) == (SERVED_PRIORITY_SEED, 1)
    np.testing.assert_array_equal(priority["levels"], prefix_levels(n))
    assert (served.uns["cellucid"]["id"], served.uns["note"]) == ("served", "kept")
    return order


def _check_prefix_counts(served: ad.AnnData) -> None:
    x = served.X
    assert sparse.issparse(x) and x.format == "csc"
    levels = np.asarray(served.uns["cellucid"]["priority"]["levels"])
    counts = served.varm[PREFIX_COUNTS_KEY]
    assert counts.dtype == np.uint32 and counts.shape == (x.shape[1], len(levels))
    dense = x.toarray() != 0
    for gene in range(x.shape[1]):
        rows = x.indices[x.indptr[gene] : x.indptr[gene + 1]]
        assert (np.diff(rows) > 0).all(), f"gene {gene}: rows ascending"
        for j, level in enumerate(levels):
            assert counts[gene, j] == dense[:level, gene].sum(), (gene, level)


# =============================================================================
# Objects
# =============================================================================


@pytest.mark.parametrize("form", ["csr", "csc", "dense"])
def test_an_object_is_served_in_priority_order_bit_for_bit(form: str) -> None:
    source = _adata(form=form)
    store = PriorityStore.of_anndata(source)
    try:
        served, keys = _read(store)
        _check(served, source)
        if form == "dense":
            assert isinstance(served.X, np.ndarray)
            assert PREFIX_COUNTS_KEY not in served.varm
        else:
            # In memory, a cell-major matrix is served gene-major (a gene is one range).
            _check_prefix_counts(served)
        # The CSR layer in obsm and the graph stay cell-major; the other layers keep their form.
        assert served.layers["counts"].format == "csc"
        assert served.obsp["connectivities"].format == "csr"
        # Uncompressed, unfiltered numbers: the browser reads them as they are.
        for key, value in keys.items():
            if key.endswith(".zarray"):
                document = json.loads(value)
                assert document["compressor"] is None, key
                assert document["filters"] in (None, [{"id": "vlen-utf8"}]), key
    finally:
        store.close()
    assert source.X.format == "csr" if form == "csr" else True


def test_the_object_is_not_changed_and_its_own_order_names_the_cells() -> None:
    source = _adata()
    before = source.copy()
    store = PriorityStore.of_anndata(source)
    try:
        served, _ = _read(store)
    finally:
        store.close()
    order = served.obs[ORIGINAL_POSITION_COLUMN].to_numpy().astype(np.int64)
    expected = CellPriority.of(source.n_obs, SERVED_PRIORITY_SEED)
    assert expected is not None
    np.testing.assert_array_equal(order, expected.order)
    assert int(served.uns["cellucid"]["priority"]["seed"]) == SERVED_PRIORITY_SEED
    _same(source.X, before.X)
    assert list(source.obs_names) == list(before.obs_names)


@pytest.mark.parametrize("case", ["verified", "verified file", "column only", "column only file", "subset", "not the seed's order"])
def test_a_cell_keeps_the_original_position_prepare_gives_it_on_every_path(tmp_path: Path, case: str) -> None:
    """A cell's original position is the one prepare() gives it, by the one
    rule for a stored order (`stored_order`): an object in a verified stored
    order (read back from a file prepare() wrote, its record kept) is put
    back first and every cell keeps its position -- served as an object, the
    file served as it is, prepared again; a position column with no record,
    a record that no longer describes the cells (a subset, a column not the
    seed's order) name nothing, the cells taken as they are, with the
    warnings prepare() gives. A session saved on one path applies to the same
    cells on another."""
    import cellucid

    source = _adata(n=200, m=6)
    source.obs_names = pd.Index([f"c{i}" for i in range(source.n_obs)], dtype=object)
    path = tmp_path / "ready.h5ad"
    cellucid.prepare(path, format="h5ad", data=source, dataset_id="t", dataset_name="T")
    stored = ad.read_h5ad(path)
    truth = dict(zip(stored.obs_names, stored.obs[ORIGINAL_POSITION_COLUMN].astype(int), strict=True))
    # The file names each cell by its row in the source.
    assert truth == {f"c{i}": i for i in range(source.n_obs)}
    data = stored
    said: tuple[str, ...] = ()
    if case.startswith("column only"):
        del data.uns["cellucid"]["priority"]
        said = ("has no uns['cellucid']['priority']", "is not written")
    elif case == "subset":
        data = stored[:50].copy()
        said = ("written for 200 cells, and there are 50", "is not written")
    elif case == "not the seed's order":
        column = data.obs[ORIGINAL_POSITION_COLUMN].to_numpy().copy()
        column[[0, 1]] = column[[1, 0]]
        data.obs[ORIGINAL_POSITION_COLUMN] = column
        said = ("is not its seed's order", "is not written")
    # Each row its own position, unless the stored order verifies.
    want = truth if case.startswith("verified") else {name: row for row, name in enumerate(data.obs_names)}

    def warned(make: Callable[[], Any]) -> tuple[Any, list[str]]:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            value = make()
        return value, [str(item.message) for item in caught if issubclass(item.category, UserWarning)]

    if case.endswith("file"):
        other = tmp_path / "given.h5ad"
        data.write_h5ad(other)
        source_store = classify_anndata_path(other)
        assert file_order(source_store) == ("priority" if case == "verified file" else "as-given")
        if case == "verified file":
            # Served as it is: the browser reads its record by the same rule.
            return
        store, served_said = warned(lambda: PriorityStore.of_file(source_store))
    else:
        store, served_said = warned(lambda: PriorityStore.of_anndata(data))
    try:
        served, _ = _read(store)
    finally:
        store.close()
    positions = dict(zip(served.obs_names, served.obs[ORIGINAL_POSITION_COLUMN].astype(int), strict=True))
    assert positions == want
    again = tmp_path / "again.h5ad"
    _, prepared_said = warned(lambda: cellucid.prepare(
        again, format="h5ad", data=data, dataset_id="t", dataset_name="T", priority_seed=SERVED_PRIORITY_SEED
    ))
    reprepared = ad.read_h5ad(again)
    assert dict(zip(reprepared.obs_names, reprepared.obs[ORIGINAL_POSITION_COLUMN].astype(int), strict=True)) == want
    assert list(reprepared.obs_names) == list(served.obs_names)
    for part in said:
        assert any(part in text for text in served_said), (part, served_said)
        assert any(part in text for text in prepared_said), (part, prepared_said)
    if not said:
        assert served_said == [] and prepared_said == []
    if case == "verified":
        # Put back in order, row p is the cell the viewer names p.
        restored = cellucid.restore_original_order(ad.read_h5ad(path))
        assert [positions[name] for name in restored.obs_names] == list(range(source.n_obs))


def test_long_genes_span_chunks_and_are_sorted_once() -> None:
    """A gene with more values than a chunk holds is read across chunks, the
    same entries as one sort of the whole gene."""
    n = 70_000
    rng = np.random.default_rng(1)
    column = sparse.csc_matrix(rng.random((n, 2)).astype(np.float32))
    source = ad.AnnData(X=column, obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
    store = PriorityStore.of_anndata(source)
    try:
        served, keys = _read(store)
        assert "X/indices/2" in keys
        order = served.obs[ORIGINAL_POSITION_COLUMN].to_numpy().astype(np.int64)
        _same(served.X, source.X[order])
        _check_prefix_counts(served)
    finally:
        store.close()


# =============================================================================
# Files
# =============================================================================


@pytest.mark.parametrize("layout", ["h5ad-csr", "h5ad-csc", "zarr", "zarr-zip"])
def test_a_file_is_served_in_priority_order_from_the_file(tmp_path: Path, layout: str) -> None:
    source = _adata(form="csc" if layout == "h5ad-csc" else "csr")
    if layout.startswith("h5ad"):
        path = tmp_path / "source.h5ad"
        source.write_h5ad(path, compression="gzip")
    else:
        path = tmp_path / "source.zarr"
        source.write_zarr(path)
        if layout == "zarr-zip":
            shutil.make_archive(str(tmp_path / "source.zarr"), "zip", root_dir=tmp_path, base_dir="source.zarr")
            path = tmp_path / "source.zarr.zip"
    classified = classify_anndata_path(path)
    assert file_order(classified) == "as-given"
    store = PriorityStore.of_file(classified)
    try:
        served, _ = _read(store)
        _check(served, ad.read_h5ad(path) if layout.startswith("h5ad") else source)
        # On disk, each matrix keeps its form: a cell-major X is read by its rows.
        assert served.X.format == ("csc" if layout == "h5ad-csc" else "csr")
        if served.X.format == "csc":
            _check_prefix_counts(served)
        else:
            assert PREFIX_COUNTS_KEY not in served.varm
        # Arrays stored in chunks are decoded once each, as they are first read,
        # into the scratch folder (rows in priority order are one fancy index
        # of the decoded copy); nothing else of the file is written.
        folder = store._folder.folder
        if folder is not None:
            assert {path.name.split("-")[0] for path in folder.iterdir()} <= {"owner", "source"}
    finally:
        store.close()
    assert folder is None or not folder.exists()


def test_a_file_in_priority_order_is_served_as_it_is(tmp_path: Path) -> None:
    from cellucid.prepare_data import prepare

    path = tmp_path / "prepared.zarr"
    prepare(path, format="zarr", data=_adata(), dataset_id="prepared", dataset_name="Prepared")
    assert file_order(classify_anndata_path(path)) == "priority"
    server = AnnDataServer(path, port=0, quiet=True, serve_web_ui=False)
    try:
        assert server._served.priority is None and server._served.store is not None
        assert server._served.store.path == path.resolve()
    finally:
        server.stop()


# =============================================================================
# Bytes a subsample reads
# =============================================================================


def _chunks_of(document: dict, rows: int) -> set[int]:
    """The chunk numbers along the first axis that rows ``[0, rows)`` lie in."""
    size = document["chunks"][0]
    return set(range(-(-rows // size))) if rows else set()


def test_a_prefix_of_every_per_cell_array_lies_in_its_first_chunks() -> None:
    """A subsample of ``k`` cells under the served seed reads, of every per-cell
    array, only the chunks of rows ``[0, k)`` (and of a CSR matrix the
    entries before row ``k``; of a gene, the entries its prefix count names)."""
    n = 3_000
    source = _adata(n=n, m=40)
    store = PriorityStore.of_anndata(source)
    try:
        served, keys = _read(store)
        metadata = json.loads(keys[".zmetadata"])["metadata"]
        levels = np.asarray(served.uns["cellucid"]["priority"]["levels"])
        counts = served.varm[PREFIX_COUNTS_KEY]
        x = served.X
        for k in (30, 300, 3_000):
            level = int(np.searchsorted(levels, k)) if k < n else None
            budget = 0
            for array in ("obs/score", f"obs/{ORIGINAL_POSITION_COLUMN}", "obsm/X_umap"):
                document = metadata[f"{array}/.zarray"]
                row_bytes = np.dtype(document["dtype"]).itemsize * int(np.prod(document["chunks"][1:]))
                budget += len(_chunks_of(document, k)) * document["chunks"][0] * row_bytes
            gene_bytes = 0
            for gene in range(x.shape[1]):
                stored = int(x.indptr[gene + 1] - x.indptr[gene])
                count = stored if level is None or level >= len(levels) else int(counts[gene, level])
                gene_bytes += count * 8
            # Every gene read up to its level: its prefix's values and rows.
            if k == 30:
                small = budget, gene_bytes
            elif k == n:
                whole = budget, gene_bytes
        # The values a 1% subsample reads of the genes are about 1% of every gene's.
        assert small[1] <= 0.02 * whole[1]
    finally:
        store.close()


def test_served_bytes_of_a_prefix_grow_with_k(tmp_path: Path) -> None:
    """Over HTTP: reading the first ``k`` rows of an embedding and of a
    cell-major X (as the viewer's prefix path does) transfers the chunks of
    the prefix alone -- at most one partly used chunk per array -- where the
    same subsample of the file in its own order touches every chunk."""
    n, m = 200_000, 20
    rng = np.random.default_rng(0)
    x = sparse.random(n, m, density=0.25, format="csr", dtype=np.float32, random_state=0)
    source = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
        obsm={"X_umap": rng.random((n, 2)).astype(np.float32)},
    )
    path = tmp_path / "cells.h5ad"
    source.write_h5ad(path)
    server = AnnDataServer(path, host="127.0.0.1", port=0, quiet=True, serve_web_ui=False,
                           dataset_id="bytes", dataset_name="Bytes")
    server.start_background()
    try:
        assert server._server is not None
        host, port = server._server.server_address[:2]
        base = "/_cellucid/anndata/bytes.zarr"

        def get(key: str) -> bytes:
            connection = http.client.HTTPConnection(host, port, timeout=60)
            try:
                connection.request("GET", f"{base}/{key}")
                response = connection.getresponse()
                assert response.status == 200, key
                return response.read()
            finally:
                connection.close()

        metadata = json.loads(get(".zmetadata"))["metadata"]
        pointers = metadata["X/indptr/.zarray"]
        indptr = np.concatenate([
            np.frombuffer(get(f"X/indptr/{i}"), dtype=pointers["dtype"])
            for i in range(-(-(n + 1) // pointers["chunks"][0]))
        ])[: n + 1]

        def prefix_bytes(k: int) -> int:
            total = sum(len(get(f"obsm/X_umap/{chunk}.0")) for chunk in _chunks_of(metadata["obsm/X_umap/.zarray"], k))
            for part in ("indices", "data"):
                total += sum(len(get(f"X/{part}/{chunk}")) for chunk in _chunks_of(metadata[f"X/{part}/.zarray"], int(indptr[k])))
            return total

        k = n // 100
        small, whole = prefix_bytes(k), prefix_bytes(n)
        one_chunk_each = sum(
            np.dtype(metadata[f"{a}/.zarray"]["dtype"]).itemsize * int(np.prod(metadata[f"{a}/.zarray"]["chunks"]))
            for a in ("obsm/X_umap", "X/indices", "X/data")
        )
        prefix_values = k * 2 * 4 + int(indptr[k]) * 8
        assert small <= prefix_values + one_chunk_each
        assert small < whole / 5
        # The same cells of the file in its own order lie in every chunk of its X.
        kept = np.sort(CellPriority.of(n, SERVED_PRIORITY_SEED).order[:k])  # type: ignore[union-attr]
        entries_chunk = metadata["X/indices/.zarray"]["chunks"][0]
        own = np.unique(np.asarray(x.indptr)[kept] // entries_chunk)
        assert len(own) == -(-x.nnz // entries_chunk)
    finally:
        server.stop()


# =============================================================================
# The same arrays as prepare()
# =============================================================================


def test_the_served_arrays_are_the_ones_prepare_writes(tmp_path: Path) -> None:
    """The store a server presents holds the arrays ``prepare()`` writes for
    the same object: the same rows, gene-major entries and prefix counts, as
    cellucid-r's served store holds ``cellucid_prepare()``'s."""
    from cellucid.prepare_data import prepare

    source = _adata(n=2_000, m=12)
    del source.raw
    path = tmp_path / "prepared.zarr"
    prepare(path, format="zarr", data=source, dataset_id="served", dataset_name="Served")
    prepared = zarr.open_group(str(path), mode="r")
    store = PriorityStore.of_anndata(source)
    try:
        _, keys = _read(store)
        memory = {key: default_buffer_prototype().buffer.from_bytes(value) for key, value in keys.items()}
        served = zarr.open_group(zarr.storage.MemoryStore(store_dict=memory), mode="r", zarr_format=2)
        for key in ("X/indptr", "X/indices", "X/data", f"obs/{ORIGINAL_POSITION_COLUMN}",
                    f"varm/{PREFIX_COUNTS_KEY}", "obsm/X_umap", "uns/cellucid/priority/levels"):
            ours, theirs = served[key], prepared[key]
            assert ours.dtype == theirs.dtype, key
            assert ours.shape == theirs.shape, key
            assert ours.chunks == theirs.chunks, key
            np.testing.assert_array_equal(ours[...], theirs[...], err_msg=key)
    finally:
        store.close()


# =============================================================================
# The server
# =============================================================================


def test_the_server_publishes_the_store_and_the_widget_reads_the_same_keys(tmp_path: Path) -> None:
    source = _adata()
    server = AnnDataServer(source, port=0, quiet=True, serve_web_ui=False)
    try:
        document = server._datasets_document()["datasets"][0]
        assert (document["id"], document["anndata"]["file"]) == (
            "served",
            {"url": "_cellucid/anndata/served.zarr/", "kind": "zarr", "layout": "folder"},
        )
        assert (server.server_info["format"], server.server_info["n_cells"]) == ("zarr", source.n_obs)
        assert ORIGINAL_POSITION_COLUMN in server.contents.obs_columns
        store = server._served.priority
        assert store is not None
        for key in (".zmetadata", "X/indptr/0", f"obs/{ORIGINAL_POSITION_COLUMN}/0", f"varm/{PREFIX_COUNTS_KEY}/0.0"):
            status, headers, body = server.answer_request("GET", f"_cellucid/anndata/served.zarr/{key}", None)
            assert (status, body) == (200, store.get(key)), key
            assert headers["ETag"] == store.etag
        status, headers, body = server.answer_request("GET", "_cellucid/anndata/served.zarr/X/indptr/0", "bytes=0-3")
        assert (status, body) == (206, store.get("X/indptr/0")[:4])
        for absent in ("X/indptr/99", "X/nothing/0", "X/indptr/x", "../x", "obs"):
            assert server.answer_request("GET", f"_cellucid/anndata/served.zarr/{absent}", None)[0] == 404, absent
    finally:
        server.stop()
    assert not store._finalizer.alive


def test_the_served_chunks_are_the_ones_prepare_writes(tmp_path: Path) -> None:
    """Every array the served store and ``prepare()`` both hold is chunked alike,
    by the one set of constants (:mod:`cellucid.prepare_data._store`): per-cell
    arrays by ``CELL_CHUNK_ELEMENTS``, sliced expression entries by
    ``SLICED_CHUNK_ELEMENTS``, anything read whole by ``WHOLE_CHUNK_ELEMENTS``
    -- on a dataset whose per-cell arrays span several chunks."""
    from cellucid.prepare_data import prepare
    from cellucid.prepare_data._store import CELL_CHUNK_ELEMENTS

    n = CELL_CHUNK_ELEMENTS + 4_000
    rng = np.random.default_rng(5)
    x = sparse.random(n, 6, density=0.2, format="csr", dtype=np.float32, random_state=5)
    source = ad.AnnData(
        X=x,
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(rng.choice(["a", "b", "c"], n)),
                "score": rng.random(n).astype(np.float32),
                "name": [f"n{i % 11}" for i in range(n)],
            },
            index=pd.Index([f"cell-{i}" for i in range(n)], dtype=object),
        ),
        var=pd.DataFrame(index=pd.Index([f"g{j}" for j in range(6)], dtype=object)),
    )
    source.obsm["X_umap"] = rng.random((n, 2)).astype(np.float32)
    path = tmp_path / "prepared.zarr"
    prepare(path, format="zarr", zarr_format=2, data=source, dataset_id="served", dataset_name="Served")
    prepared = zarr.open_group(str(path), mode="r")
    store = PriorityStore.of_anndata(source)
    try:
        _, keys = _read(store)
        memory = {key: default_buffer_prototype().buffer.from_bytes(value) for key, value in keys.items()}
        served = zarr.open_group(zarr.storage.MemoryStore(store_dict=memory), mode="r", zarr_format=2)
        compared = []
        for area in ("obs", "obsm", "X", "var", "varm"):

            def arrays(group: zarr.Group, prefix: str) -> dict[str, zarr.Array]:
                found = {}
                for name, node in group.members():
                    if isinstance(node, zarr.Array):
                        found[f"{prefix}/{name}"] = node
                    else:
                        found.update(arrays(node, f"{prefix}/{name}"))
                return found

            ours, theirs = arrays(served[area], area), arrays(prepared[area], area)
            for key in sorted(set(ours) & set(theirs)):
                assert ours[key].chunks == theirs[key].chunks, key
                compared.append(key)
        for key in (f"obs/{ORIGINAL_POSITION_COLUMN}", "obs/score", "obs/cell_type/codes", "obsm/X_umap", "X/data", "X/indices"):
            assert key in compared, key
        assert served[f"obs/{ORIGINAL_POSITION_COLUMN}"].chunks == (CELL_CHUNK_ELEMENTS,)
    finally:
        store.close()


@pytest.mark.parametrize("form", ["csr", "csc", "dense"])
def test_an_object_is_served_from_a_snapshot_edits_after_never_mix_in(form: str) -> None:
    """The store holds the object as it was when it was made: an array changed
    in place afterwards (``np.log1p(adata.X.data, out=adata.X.data)``) is not
    seen, in any chunk, so the bytes under the store's one ETag never mix two
    versions."""
    source = _adata(n=900, m=7, form=form, seed=8)
    before = source.copy()
    store = PriorityStore.of_anndata(source)
    try:
        values = source.X.data if sparse.issparse(source.X) else source.X
        np.log1p(values, out=values)
        source.layers["dense"] += 100
        source.obsm["X_umap"][:] = -1
        served, _ = _read(store)
    finally:
        store.close()
    order = served.obs[ORIGINAL_POSITION_COLUMN].to_numpy().astype(np.int64)
    _same(sparse.csr_matrix(served.X), sparse.csr_matrix(before.X)[order])
    _same(served.layers["dense"], before.layers["dense"][order])
    _same(served.obsm["X_umap"], before.obsm["X_umap"][order])


def test_an_object_s_matrix_is_stored_as_prepare_stores_it() -> None:
    """float64 values float32 holds are served as float32, stored zeros
    dropped, duplicate coordinates summed (``prepare()``'s functions); integer
    counts above 2^24 stay float64; a non-finite value is refused as by
    ``prepare()``."""
    n, m = 300, 5
    rows = np.array([0, 0, 3, 3, 7, 9, 299], dtype=np.int32)
    cols = np.array([1, 1, 2, 4, 0, 3, 4], dtype=np.int32)
    vals = np.array([0.5, 0.25, 0.0, 2.0, 1.5, 0.0, 3.0], dtype=np.float64)
    x = sparse.coo_matrix((vals, (rows, cols)), shape=(n, m)).tocsr()
    # tocsr() sums duplicates; put the duplicate and the stored zeros back as a file or scipy op leaves them.
    x = sparse.csr_matrix((vals[np.argsort(rows, kind="stable")], cols[np.argsort(rows, kind="stable")],
                           np.concatenate([[0], np.cumsum(np.bincount(rows, minlength=n))])), shape=(n, m))
    assert not x.has_canonical_format
    adata = ad.AnnData(X=x, obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
    store = PriorityStore.of_anndata(adata)
    try:
        served, _ = _read(store)
    finally:
        store.close()
    assert served.X.dtype == np.float32
    order = served.obs[ORIGINAL_POSITION_COLUMN].to_numpy().astype(np.int64)
    expected = sparse.csr_matrix(x)
    expected.sum_duplicates()
    expected.eliminate_zeros()
    _same(served.X, expected.astype(np.float32)[order])
    assert served.X.nnz == expected.nnz == 4
    big = ad.AnnData(X=sparse.csr_matrix(np.array([[2.0**25, 1.0]] * 3)), obs=pd.DataFrame(index=["a", "b", "c"]))
    store = PriorityStore.of_anndata(big)
    try:
        assert _read(store)[0].X.dtype == np.float64
    finally:
        store.close()
    for bad in (np.nan, np.inf):
        broken = ad.AnnData(X=sparse.csr_matrix(np.array([[bad, 1.0], [0.0, 2.0]])), obs=pd.DataFrame(index=["a", "b"]))
        with pytest.raises(ValueError, match="finite"):
            PriorityStore.of_anndata(broken)


def test_a_killed_process_s_scratch_folder_is_removed_by_the_next_store(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from cellucid import _served_store

    monkeypatch.setattr(_served_store.tempfile, "tempdir", str(tmp_path))
    dead = tmp_path / "cellucid-served-dead"
    dead.mkdir()
    (dead / "owner").write_text("999999999")
    (dead / "X-data.npy").write_bytes(b"x" * 100)
    alive = tmp_path / "cellucid-served-alive"
    alive.mkdir()
    (alive / "owner").write_text(str(_served_store.os.getppid()))
    unnamed = tmp_path / "cellucid-served-unnamed"
    unnamed.mkdir()
    store = PriorityStore.of_anndata(_adata(n=50, m=4))
    try:
        assert not dead.exists()
        assert alive.exists() and unnamed.exists()
        mine = store._folder.folder
        assert mine is not None and (mine / "owner").read_text() == str(_served_store.os.getpid())
    finally:
        store.close()
    assert not mine.exists()


def test_every_file_a_store_maps_is_let_go_however_it_ends(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A store maps files: an object's snapshots and gene-major layouts in its
    scratch folder, a file's arrays decoded there, a file's uncompressed arrays
    where they lie. Closed, dropped without a close, or failing while it is
    made, it lets every mapping go before its folder is removed (Windows
    removes no mapped file, and lets no program replace one), and the source
    file it read can be replaced."""
    from cellucid import _served_store

    monkeypatch.setattr(_served_store.tempfile, "tempdir", str(tmp_path))
    adata = _adata(n=50, m=4)

    def scratch_folders() -> list[Path]:
        return sorted(tmp_path.glob("cellucid-served-*"))

    # Closed.
    store = PriorityStore.of_anndata(adata)
    _read(store)
    folder = store._folder.folder
    assert folder is not None and mapped_files(folder) != []
    store.close()
    assert mapped_files(tmp_path) == [] and not folder.exists()
    with pytest.raises(RuntimeError, match="closed"):
        store.get("X/data/0")

    # Dropped: the store lets go as it is collected, with no close.
    store = PriorityStore.of_anndata(adata)
    _read(store)
    folder = store._folder.folder
    assert folder is not None and mapped_files(folder) != []
    del store
    assert mapped_files(tmp_path) == [] and not folder.exists()

    # A file: its arrays in one piece mapped where they lie, chunked ones decoded.
    source = tmp_path / "source.h5ad"
    adata.write_h5ad(source)
    store = PriorityStore.of_file(classify_anndata_path(source))
    _read(store)
    assert mapped_files(source) != []
    store.close()
    assert mapped_files(tmp_path) == [] and scratch_folders() == []
    replacement = tmp_path / "replacement.h5ad"
    shutil.copy2(source, replacement)
    replacement.replace(source)

    # Failing while it is made, its layout's files mapped in the frames the error passed through.
    def fail(*_args: object) -> None:
        raise RuntimeError("layout failed")

    monkeypatch.setattr(_served_store._GeneMajor, "_finished", fail)
    with pytest.raises(RuntimeError, match="layout failed") as raised:
        PriorityStore.of_anndata(adata)
    # The error is alive (as it is while the store is closed): it holds no mapping.
    assert raised.value.__traceback__ is not None
    assert mapped_files(tmp_path) == [] and scratch_folders() == []


def test_a_dense_file_matrix_in_whole_row_chunks_is_laid_out_once(tmp_path: Path) -> None:
    """A dense X stored in chunks of whole rows is read once, by blocks of rows,
    into a gene-major layout, and every gene after is a slice of it: the file's
    chunks are decoded once, not once a gene."""
    import h5py

    # A strip of 87 genes (2^18 values over 3,000 rows) is narrower than a chunk of 200.
    n, m = 3_000, 200
    rng = np.random.default_rng(4)
    x = rng.random((n, m)).astype(np.float32)
    adata = ad.AnnData(X=x, obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
    path = tmp_path / "dense.h5ad"
    adata.write_h5ad(path)
    with h5py.File(path, "r+") as file:
        del file["X"]
        file.create_dataset("X", data=x, chunks=(100, m), compression="gzip")
        file["X"].attrs["encoding-type"] = "array"
        file["X"].attrs["encoding-version"] = "0.2.0"
    store = PriorityStore.of_file(classify_anndata_path(path))
    try:
        reads = []
        source = store._layouts[0].array
        original = type(source).__getitem__

        class Counting(type(source)):  # type: ignore[misc]
            def __getitem__(self, key):  # type: ignore[no-untyped-def]
                reads.append(key)
                return original(self, key)

        source.__class__ = Counting
        served, _ = _read(store)
        once = len(reads)
        np.testing.assert_array_equal(served.X, x[served.obs[ORIGINAL_POSITION_COLUMN].to_numpy().astype(np.int64)])
        _read(store)
        assert len(reads) == once, "laid out once"
    finally:
        store.close()


@pytest.mark.parametrize("stored", ["chunked gzip", "contiguous", "zarr"])
def test_rows_in_any_order_are_one_take_from_the_decoded_or_mapped_source(
    tmp_path: Path, stored: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A source array is read by rows in any order (repeats included, and a
    column window of a 2-D one) exactly as numpy indexes it: an uncompressed
    HDF5 array in one piece mapped in place (no copy), any other decoded once
    into the scratch folder, each source chunk when a read first needs it.
    Every mapping -- the source's, the decoded copies' -- is let go on close,
    so the files can be removed and replaced on every system."""
    import h5py

    from cellucid._served_store import _Rows, _Scratch, _Vector

    # h5py names the default driver "windows" on Windows (HDF5's Windows
    # driver is the sec2 driver there): an array is mapped by what its file's
    # driver is, not by that name.
    monkeypatch.setattr(h5py.File, "driver", property(lambda _file: "windows"))
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    rng = np.random.default_rng(3)
    matrix = rng.random((5003, 7)).astype(np.float32)
    vector = rng.integers(0, 1000, 20011).astype(np.int32)
    if stored == "zarr":
        group = zarr.open_group(str(tmp_path / "s.zarr"), mode="w", zarr_format=2)
        group.create_array("m", data=matrix, chunks=(512, 3))
        group.create_array("v", data=vector, chunks=(1000,))
        arrays = (group["m"], group["v"])
    else:
        file = h5py.File(tmp_path / "s.h5", "w")
        options = {"chunks": True, "compression": "gzip"} if stored == "chunked gzip" else {}
        file.create_dataset("m", data=matrix, **options)
        file.create_dataset("v", data=vector, **options)
        arrays = (file["m"], file["v"])
    closers: list = []
    layouts: list = []
    scratch = _Scratch(closers, threading.Event(), layouts)
    rows_of_matrix = _Rows(arrays[0], scratch)
    rows_of_vector = _Vector(arrays[1], scratch)
    picks = rng.integers(0, 5003, 3000)
    picks[:10] = picks[0]
    np.testing.assert_array_equal(rows_of_matrix.take(picks), matrix[picks])
    np.testing.assert_array_equal(rows_of_matrix.take(picks, slice(2, 5)), matrix[picks, 2:5])
    np.testing.assert_array_equal(rows_of_matrix.take(np.zeros(0, dtype=np.int64)), matrix[:0])
    starts = rng.integers(0, 20000, 500)
    ends = starts + rng.integers(0, 11, 500)
    want = np.concatenate([vector[a:b] for a, b in zip(starts, ends, strict=True)])
    np.testing.assert_array_equal(rows_of_vector.ranges(starts, ends), want)
    copied = scratch.folder is not None and any(path.name.startswith("source-") for path in scratch.folder.iterdir())
    assert copied == (stored != "contiguous")
    assert mapped_files(tmp_path) != []
    scratch.release()
    for close in closers:
        close()
    assert mapped_files(tmp_path) == []
    assert scratch.folder is None or not scratch.folder.exists()
    with pytest.raises(RuntimeError, match="closed"):
        rows_of_matrix.take(picks)
