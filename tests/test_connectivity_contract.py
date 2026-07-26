from __future__ import annotations

import gzip
import inspect
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.connectivity_contract import (
    connectivity_index_storage,
    validate_connectivity_edges,
)
from cellucid.prepare_data import prepare


def _graph() -> np.ndarray:
    graph = np.zeros((5, 5), dtype=np.float64)
    for source, destination, weight in [
        (3, 4, 4.25),
        (0, 3, 2.5),
        (1, 3, 3.75),
        (0, 2, 1.125),
    ]:
        graph[source, destination] = weight
        graph[destination, source] = weight
    return graph


def _prepare_kwargs(out_dir, connectivities):
    return {
        "latent_space": np.arange(10, dtype=np.float32).reshape(5, 2),
        "obs": pd.DataFrame(
            {"score": np.arange(5, dtype=np.float32)},
            index=[f"cell-{index}" for index in range(5)],
        ),
        "connectivities": connectivities,
        "X_umap_2d": np.arange(10, dtype=np.float32).reshape(5, 2),
        "out_dir": out_dir,
        "centroid_min_points": 1,
        "dataset_name": "Exact connectivity",
        "dataset_id": "exact-connectivity",
        "obs_categorical_dtype": "uint16",
        "force": True,
    }


def test_prepare_exposes_only_the_canonical_connectivity_paths() -> None:
    parameters = inspect.signature(prepare).parameters
    assert "connectivity_manifest_filename" not in parameters
    assert "connectivity_binary_dirname" not in parameters


@pytest.mark.parametrize("sparse_graph", [False, True])
@pytest.mark.parametrize("graph_dtype", [np.float64, np.bool_])
def test_connectivity_validation_is_exact_canonical_and_non_mutating(
    sparse_graph: bool,
    graph_dtype,
) -> None:
    dense = _graph().astype(graph_dtype)
    if sparse_graph:
        rows, columns = np.nonzero(dense)
        reverse_order = np.arange(rows.size - 1, -1, -1)
        candidate = sparse.coo_matrix(
            (
                dense[rows, columns][reverse_order],
                (
                    rows[reverse_order],
                    columns[reverse_order],
                ),
            ),
            shape=dense.shape,
        )
    else:
        candidate = dense.copy()
    dense_before = candidate.toarray() if sparse_graph else candidate.copy()

    edges = validate_connectivity_edges(candidate, n_cells=5)

    assert edges.sources.tolist() == [0, 0, 1, 3]
    assert edges.destinations.tolist() == [2, 3, 3, 4]
    expected_weights = (
        [1.0, 1.0, 1.0, 1.0]
        if graph_dtype == np.bool_
        else [
            1.125,
            2.5,
            3.75,
            4.25,
        ]
    )
    assert edges.weights.dtype == np.dtype("<f8")
    assert edges.weights.tolist() == expected_weights
    assert edges.n_edges == 4
    assert edges.max_neighbors == 3
    assert edges.index_dtype == "uint16"
    assert edges.index_bytes == 2
    assert edges.sources.dtype == np.dtype("<u2")
    assert edges.destinations.dtype == np.dtype("<u2")
    np.testing.assert_array_equal(
        candidate.toarray() if sparse_graph else candidate,
        dense_before,
    )


@pytest.mark.parametrize(
    ("matrix", "message"),
    [
        (np.array([[0.0, 1.0], [0.0, 0.0]]), "symmetric"),
        (np.array([[0.0, 0.5], [0.75, 0.0]]), "symmetric"),
        (np.array([[0.0, -1.0], [-1.0, 0.0]]), "non-negative"),
        (np.array([[0.0, np.nan], [np.nan, 0.0]]), "finite"),
        (np.array([[1.0, 0.0], [0.0, 0.0]]), "diagonal"),
        (np.zeros((2, 3)), "shape"),
    ],
)
@pytest.mark.parametrize("sparse_graph", [False, True])
def test_invalid_connectivity_is_rejected_without_reinterpretation(
    matrix: np.ndarray,
    message: str,
    sparse_graph: bool,
) -> None:
    candidate = sparse.csr_matrix(matrix) if sparse_graph else matrix
    with pytest.raises(ValueError, match=message):
        validate_connectivity_edges(candidate, n_cells=2)


@pytest.mark.parametrize("sparse_format", ["coo", "csr", "csc", "lil", "bsr"])
def test_sparse_stored_zero_entries_are_rejected(sparse_format: str) -> None:
    row = np.array([0, 1], dtype=np.int32)
    column = np.array([1, 0], dtype=np.int32)
    candidate = sparse.coo_matrix(
        (np.zeros(2, dtype=np.float64), (row, column)),
        shape=(2, 2),
    ).asformat(sparse_format)

    with pytest.raises(ValueError, match="stored zero"):
        validate_connectivity_edges(candidate, n_cells=2)


@pytest.mark.parametrize("sparse_format", ["coo", "csr", "csc", "lil", "bsr"])
def test_duplicate_sparse_coordinates_are_rejected_before_coalescing(
    sparse_format: str,
) -> None:
    data = np.array([0.5, 0.5, 0.5, 0.5], dtype=np.float64)
    indices = np.array([1, 1, 0, 0], dtype=np.int32)
    indptr = np.array([0, 2, 4], dtype=np.int32)
    if sparse_format == "coo":
        candidate = sparse.coo_matrix(
            (
                data,
                (
                    np.array([0, 0, 1, 1], dtype=np.int32),
                    indices,
                ),
            ),
            shape=(2, 2),
        )
    elif sparse_format == "csr":
        candidate = sparse.csr_matrix(
            (data, indices, indptr),
            shape=(2, 2),
        )
    elif sparse_format == "csc":
        candidate = sparse.csc_matrix(
            (data, indices, indptr),
            shape=(2, 2),
        )
    elif sparse_format == "lil":
        candidate = sparse.lil_matrix((2, 2), dtype=np.float64)
        candidate.rows[0] = [1, 1]
        candidate.data[0] = [0.5, 0.5]
        candidate.rows[1] = [0, 0]
        candidate.data[1] = [0.5, 0.5]
    else:
        candidate = sparse.bsr_matrix(
            (
                data.reshape(4, 1, 1),
                indices,
                indptr,
            ),
            shape=(2, 2),
        )

    with pytest.raises(
        ValueError,
        match=r"duplicate sparse coordinates at \(0, 1\)",
    ):
        validate_connectivity_edges(candidate, n_cells=2)


def test_empty_graph_and_cell_index_boundaries_are_exact() -> None:
    empty = validate_connectivity_edges(
        sparse.csr_matrix((4, 4), dtype=np.bool_),
        n_cells=4,
    )
    assert empty.sources.size == 0
    assert empty.destinations.size == 0
    assert empty.weights.size == 0
    assert empty.weights.dtype == np.dtype("<f8")
    assert empty.n_edges == 0
    assert empty.max_neighbors == 0

    for n_cells, dtype_name, index_bytes in [
        (65_535, "uint16", 2),
        (65_536, "uint16", 2),
        (65_537, "uint32", 4),
    ]:
        row = np.array([0, n_cells - 1], dtype=np.int64)
        column = np.array([n_cells - 1, 0], dtype=np.int64)
        graph = sparse.coo_matrix(
            (np.array([0.125, 0.125], dtype=np.float64), (row, column)),
            shape=(n_cells, n_cells),
        )
        edges = validate_connectivity_edges(graph, n_cells=n_cells)
        assert edges.sources.tolist() == [0]
        assert edges.destinations.tolist() == [n_cells - 1]
        assert edges.weights.tolist() == [0.125]
        assert edges.index_dtype == dtype_name
        assert edges.index_bytes == index_bytes
        assert edges.max_neighbors == 1

    assert connectivity_index_storage(4_294_967_296) == (
        np.uint32,
        "uint32",
        4,
    )
    for invalid_n_cells in [0, -1, True, 1.5, 4_294_967_297]:
        with pytest.raises(ValueError, match="positive|uint32"):
            connectivity_index_storage(invalid_n_cells)


@pytest.mark.parametrize("sparse_graph", [False, True])
def test_prepared_writer_emits_the_same_exact_graph_contract(
    tmp_path,
    sparse_graph: bool,
) -> None:
    graph = _graph()
    candidate = sparse.csr_matrix(graph) if sparse_graph else graph
    prepare(**_prepare_kwargs(tmp_path, candidate))

    manifest = json.loads((tmp_path / "connectivity_manifest.json").read_text(encoding="utf-8"))
    identity = json.loads((tmp_path / "dataset_identity.json").read_text(encoding="utf-8"))
    sources = np.fromfile(
        tmp_path / "connectivity" / "edges.src.bin",
        dtype="<u2",
    )
    destinations = np.fromfile(
        tmp_path / "connectivity" / "edges.dst.bin",
        dtype="<u2",
    )
    weights = np.fromfile(
        tmp_path / "connectivity" / "edges.weights.f64.bin",
        dtype="<f8",
    )

    assert sources.tolist() == [0, 0, 1, 3]
    assert destinations.tolist() == [2, 3, 3, 4]
    assert weights.tolist() == [1.125, 2.5, 3.75, 4.25]
    assert list(manifest) == [
        "format",
        "n_cells",
        "n_edges",
        "max_neighbors",
        "index_bytes",
        "index_dtype",
        "sourcesPath",
        "destinationsPath",
        "weightsPath",
        "weight_dtype",
        "weight_bytes",
        "compression",
    ]
    assert manifest == {
        "format": "edge_pairs",
        "n_cells": 5,
        "n_edges": 4,
        "max_neighbors": 3,
        "index_bytes": 2,
        "index_dtype": "uint16",
        "sourcesPath": "connectivity/edges.src.bin",
        "destinationsPath": "connectivity/edges.dst.bin",
        "weightsPath": "connectivity/edges.weights.f64.bin",
        "weight_dtype": "float64",
        "weight_bytes": 8,
        "compression": None,
    }
    assert identity["stats"]["has_connectivity"] is True
    assert identity["stats"]["n_edges"] == 4


@pytest.mark.parametrize(
    ("graph", "expected_weights"),
    [
        (np.zeros((5, 5), dtype=np.float64), []),
        (
            np.array(
                [
                    [0.0, 1.0, 0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                ],
                dtype=np.float64,
            ),
            [1.0, 1.0],
        ),
    ],
)
@pytest.mark.parametrize("compression", [None, 6])
def test_every_present_graph_emits_an_explicit_aligned_weight_payload(
    tmp_path,
    graph: np.ndarray,
    expected_weights: list[float],
    compression: int | None,
) -> None:
    kwargs = _prepare_kwargs(tmp_path, sparse.csr_matrix(graph))
    kwargs["compression"] = compression
    prepare(**kwargs)

    suffix = ".gz" if compression is not None else ""
    weights_path = tmp_path / "connectivity" / f"edges.weights.f64.bin{suffix}"
    weights_bytes = weights_path.read_bytes()
    if compression is not None:
        weights_bytes = gzip.decompress(weights_bytes)
    assert np.frombuffer(weights_bytes, dtype="<f8").tolist() == expected_weights

    manifest = json.loads((tmp_path / "connectivity_manifest.json").read_text(encoding="utf-8"))
    assert manifest["n_edges"] == len(expected_weights)
    assert manifest["weightsPath"] == (f"connectivity/edges.weights.f64.bin{suffix}")
    assert manifest["weight_dtype"] == "float64"
    assert manifest["weight_bytes"] == 8
    for stem in ["edges.src.bin", "edges.dst.bin", "edges.weights.f64.bin"]:
        assert (tmp_path / "connectivity" / f"{stem}{suffix}").is_file()


def test_rejected_connectivity_cannot_mutate_graph_or_identity_publication(
    tmp_path,
) -> None:
    graph_dir = tmp_path / "connectivity"
    graph_dir.mkdir(parents=True)
    sentinels = {
        tmp_path / "connectivity_manifest.json": b"manifest-before",
        graph_dir / "edges.src.bin": b"sources-before",
        graph_dir / "edges.dst.bin": b"destinations-before",
        graph_dir / "edges.weights.f64.bin": b"weights-before",
        tmp_path / "dataset_identity.json": b"identity-before",
    }
    for path, contents in sentinels.items():
        path.write_bytes(contents)

    invalid = _graph()
    invalid[0, 2] = 0.25
    before = invalid.copy()
    with pytest.raises(ValueError, match="symmetric"):
        prepare(**_prepare_kwargs(tmp_path, invalid))

    np.testing.assert_array_equal(invalid, before)
    for path, contents in sentinels.items():
        assert path.read_bytes() == contents


def test_direct_anndata_producer_uses_the_same_canonical_sequence() -> None:
    graph = _graph()
    adata = ad.AnnData(
        X=np.ones((5, 1), dtype=np.float32),
        obs=pd.DataFrame(index=[f"cell-{index}" for index in range(5)]),
        var=pd.DataFrame(index=["gene"]),
    )
    adata.obsm["X_umap_2d"] = np.arange(10, dtype=np.float32).reshape(5, 2)
    adata.obsp["connectivities"] = sparse.csr_matrix(graph)

    adapter = AnnDataAdapter(
        adata,
        dataset_name="Direct exact connectivity",
        dataset_id="direct-exact-connectivity",
    )
    (
        source_bytes,
        destination_bytes,
        weight_bytes,
        n_edges,
        max_neighbors,
    ) = adapter.get_connectivity_edges()

    assert np.frombuffer(source_bytes, dtype="<u2").tolist() == [0, 0, 1, 3]
    assert np.frombuffer(destination_bytes, dtype="<u2").tolist() == [2, 3, 3, 4]
    assert np.frombuffer(weight_bytes, dtype="<f8").tolist() == [
        1.125,
        2.5,
        3.75,
        4.25,
    ]
    assert n_edges == 4
    assert max_neighbors == 3

    direct_manifest = adapter.get_connectivity_manifest()
    assert direct_manifest is not None
    assert list(direct_manifest) == [
        "format",
        "n_cells",
        "n_edges",
        "max_neighbors",
        "index_bytes",
        "index_dtype",
        "sourcesPath",
        "destinationsPath",
        "weightsPath",
        "weight_dtype",
        "weight_bytes",
        "compression",
    ]
    assert direct_manifest == {
        "format": "edge_pairs",
        "n_cells": 5,
        "n_edges": 4,
        "max_neighbors": 3,
        "index_bytes": 2,
        "index_dtype": "uint16",
        "sourcesPath": "connectivity/edges.src.bin",
        "destinationsPath": "connectivity/edges.dst.bin",
        "weightsPath": "connectivity/edges.weights.f64.bin",
        "weight_dtype": "float64",
        "weight_bytes": 8,
        "compression": None,
    }


def test_direct_connectivity_compression_flag_is_exactly_boolean() -> None:
    graph = _graph()
    adata = ad.AnnData(
        X=np.ones((5, 1), dtype=np.float32),
        obs=pd.DataFrame(index=[f"cell-{index}" for index in range(5)]),
        var=pd.DataFrame(index=["gene"]),
    )
    adata.obsm["X_umap_2d"] = np.arange(10, dtype=np.float32).reshape(5, 2)
    adata.obsp["connectivities"] = sparse.csr_matrix(graph)
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Exact compression option",
        dataset_id="exact-compression-option",
    )

    for invalid in [None, 0, 1, "gzip"]:
        with pytest.raises(TypeError, match="compress must be exactly"):
            adapter.get_connectivity_edges(compress=invalid)


def test_force_false_rejects_existing_generation_without_mutating_any_file(
    tmp_path,
) -> None:
    first_graph = np.zeros((5, 5), dtype=np.float64)
    first_graph[0, 1] = first_graph[1, 0] = 1
    prepare(**_prepare_kwargs(tmp_path, first_graph))
    before = {
        path.relative_to(tmp_path): path.read_bytes()
        for path in tmp_path.rglob("*")
        if path.is_file()
    }

    second_graph = first_graph.copy()
    second_graph[2, 3] = second_graph[3, 2] = 1
    second_kwargs = _prepare_kwargs(tmp_path, second_graph)
    second_kwargs["force"] = False

    with pytest.raises(FileExistsError, match="non-empty export directory"):
        prepare(**second_kwargs)

    after = {
        path.relative_to(tmp_path): path.read_bytes()
        for path in tmp_path.rglob("*")
        if path.is_file()
    }
    assert after == before
