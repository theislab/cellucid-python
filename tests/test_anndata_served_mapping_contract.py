"""The data mapping an AnnData server serves, as the viewer reads it.

A session saved on a served AnnData object records the mapping it was shown
with: the embeddings, and the graph, gene names, expression matrix and latent
space the server chose. The identity does not carry the last four, so the
server publishes the served mapping as the ``defaults`` of
``_cellucid/mapping/candidates.json``; the viewer reads it there and writes it
into the session fingerprint. These tests pin that the document says what the
server actually serves, and that it is only metadata.
"""

from __future__ import annotations

import http.client
import json
import threading
from collections.abc import Iterator
from contextlib import contextmanager

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid.anndata_adapter import MAPPING_CANDIDATES_PATH, AnnDataAdapter
from cellucid.anndata_mapping import MappingNotFound, parse_mapping_route
from cellucid.anndata_server import AnnDataServer
from cellucid.data_mapping import validate_data_mapping


def _adata() -> ad.AnnData:
    rng = np.random.default_rng(7)
    n = 6
    adata = ad.AnnData(
        X=rng.random((n, 3)).astype(np.float32),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["a", "b", "a", "b", "a", "b"])},
            index=pd.Index([f"c{i}" for i in range(n)], dtype=object),
        ),
        var=pd.DataFrame(
            {"symbol": ["G1", "G2", "G3"]},
            index=pd.Index(["E1", "E2", "E3"], dtype=object),
        ),
    )
    adata.raw = ad.AnnData(
        X=rng.random((n, 4)).astype(np.float32),
        var=pd.DataFrame(index=pd.Index(["R1", "R2", "R3", "R4"], dtype=object)),
    )
    adata.layers["counts"] = rng.integers(0, 5, (n, 3)).astype(np.float32)
    adata.obsm["X_umap_2d"] = rng.random((n, 2)).astype(np.float32)
    adata.obsm["X_umap_3d"] = rng.random((n, 3)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((n, 5)).astype(np.float32)
    adata.obsm["velocity_umap_3d"] = rng.random((n, 3)).astype(np.float32)
    graph = np.zeros((n, n), dtype=np.float32)
    for i in range(n - 1):
        graph[i, i + 1] = graph[i + 1, i] = 1.0
    adata.obsp["connectivities"] = sparse.csr_matrix(graph)
    return adata


def _adapter(**options: object) -> AnnDataAdapter:
    return AnnDataAdapter(
        _adata(),
        dataset_name="Served mapping",
        dataset_id="served-mapping",
        quiet=True,
        **options,  # type: ignore[arg-type]
    )


def test_the_default_served_mapping_names_what_is_served() -> None:
    adapter = _adapter(serve_vector_fields=True)
    try:
        mapping = adapter.get_data_mapping()
        assert mapping == validate_data_mapping(mapping)
        assert mapping == {
            "kind": "anndata",
            # The served slots are the defaults; every other narrow obsm array
            # is enabled beside them. X_pca (5 columns) fills no slot and the
            # vector field is not a representation, so neither is enabled.
            "embeddings": {
                "enabled": [
                    {"key": "X_umap_2d", "columns": [0, 1]},
                    {"key": "X_umap_3d", "columns": [0, 1, 2]},
                ],
                "defaults": {
                    "1d": None,
                    "2d": {"key": "X_umap_2d", "columns": [0, 1]},
                    "3d": {"key": "X_umap_3d", "columns": [0, 1, 2]},
                },
            },
            "knn": None,
            "geneNames": "index",
            "expression": "X",
            "latent": None,
            "vectorFields": [
                {
                    "key": "velocity_umap_3d",
                    "embedding": {"key": "X_umap_3d", "columns": [0, 1, 2]},
                }
            ],
        }
    finally:
        adapter.close()


def test_the_served_mapping_follows_every_serving_choice() -> None:
    adapter = _adapter(
        embedding_keys={"2d": ("X_pca", [3, 1])},
        knn_key="connectivities",
        gene_id_column="symbol",
        expression_source="layer:counts",
    )
    try:
        mapping = adapter.get_data_mapping()
        assert mapping["embeddings"]["defaults"]["2d"] == {"key": "X_pca", "columns": [3, 1]}
        # The served selection of a wide array is enabled, and so is every
        # narrow array beside it; the wide array is not enabled at its first
        # columns, which nothing serves.
        assert mapping["embeddings"]["enabled"] == [
            {"key": "X_pca", "columns": [3, 1]},
            {"key": "X_umap_2d", "columns": [0, 1]},
            {"key": "X_umap_3d", "columns": [0, 1, 2]},
        ]
        assert mapping["knn"] == "connectivities"
        assert mapping["geneNames"] == "column:symbol"
        assert mapping["expression"] == "layer:counts"
        identity = adapter.get_dataset_identity()
        defaults = identity["embeddings"]["defaults"]
        for slot, ref in mapping["embeddings"]["defaults"].items():
            if ref is None:
                assert slot not in defaults
            else:
                assert identity["embeddings"]["items"][defaults[slot]]["source"] == ref
    finally:
        adapter.close()


def test_raw_expression_is_reported_as_raw() -> None:
    adapter = _adapter(expression_source="raw")
    try:
        assert adapter.get_data_mapping()["expression"] == "raw"
    finally:
        adapter.close()


def test_the_catalog_is_frozen_and_returned_as_a_copy() -> None:
    adapter = _adapter()
    try:
        first = adapter.get_mapping_candidates()
        first["defaults"]["knn"] = "tampered"
        second = adapter.get_mapping_candidates()
        assert second["defaults"] == adapter.get_data_mapping()
        assert second["kind"] == "anndata"
        assert second["n_obs"] == 6
        assert [entry["key"] for entry in second["obsm"]] == [
            "X_umap_2d",
            "X_umap_3d",
            "X_pca",
            "velocity_umap_3d",
        ]
        assert second["guessed"] == []
    finally:
        adapter.close()


def test_concurrent_first_requests_build_one_catalog() -> None:
    adapter = _adapter()
    try:
        results: list[dict[str, object]] = []
        barrier = threading.Barrier(8)

        def read() -> None:
            barrier.wait()
            results.append(adapter.get_mapping_candidates())

        threads = [threading.Thread(target=read) for _ in range(8)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        assert len(results) == 8
        assert all(result == results[0] for result in results)
    finally:
        adapter.close()


@contextmanager
def _running_server() -> Iterator[tuple[str, int]]:
    server = AnnDataServer(
        _adata(),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Served mapping",
        dataset_id="served-mapping",
        knn_key="connectivities",
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        yield str(host), int(port)
    finally:
        server.stop()


@pytest.mark.parametrize("method", ["GET", "HEAD"])
def test_the_catalog_route_serves_the_mapping_document(method: str) -> None:
    with _running_server() as (host, port):
        connection = http.client.HTTPConnection(host, port, timeout=5)
        try:
            connection.request(method, f"/{MAPPING_CANDIDATES_PATH}")
            response = connection.getresponse()
            body = response.read()
            assert response.status == 200
            assert response.getheader("Content-Type") == "application/json"
            if method == "GET":
                document = json.loads(body)
                assert document["defaults"]["knn"] == "connectivities"
                assert validate_data_mapping(document["defaults"]) == document["defaults"]
            else:
                assert body == b""
        finally:
            connection.close()


def _obsm_index(adapter: AnnDataAdapter, key: str) -> int:
    keys = [entry["key"] for entry in adapter.get_mapping_candidates()["obsm"]]
    return keys.index(key)


def test_every_enabled_embedding_is_served_on_request() -> None:
    # Every narrow obsm array is enabled, though the identity lists only the
    # served defaults; the viewer reads each one through the mapping route
    # when a view first shows it.
    adata = _adata()
    adata.obsm["X_tsne"] = np.random.default_rng(3).random((adata.n_obs, 2)).astype(np.float32)
    adapter = AnnDataAdapter(adata, dataset_name="Served mapping", dataset_id="served-mapping", quiet=True)
    try:
        enabled = adapter.get_data_mapping()["embeddings"]["enabled"]
        served = {(item.key, item.columns) for item in adapter.embeddings}
        extra = [ref for ref in enabled if (ref["key"], tuple(ref["columns"])) not in served]
        assert extra, "the fixture must enable an embedding the identity does not list"
        for ref in enabled:
            data = adapter.get_mapped_embedding(_obsm_index(adapter, ref["key"]), ref["columns"])
            values = np.frombuffer(data, dtype="<f4").reshape(adapter.n_cells, len(ref["columns"]))
            assert np.isfinite(values).all()
        # A served embedding read through the route is the same bytes.
        item = adapter.embeddings[0]
        assert adapter.get_mapped_embedding(
            _obsm_index(adapter, item.key), list(item.columns)
        ) == adapter.get_embedding_binary(item.id)
        # Columns of a wider array, chosen per axis.
        pca = _obsm_index(adapter, "X_pca")
        assert len(adapter.get_mapped_embedding(pca, [4, 0])) == adapter.n_cells * 8
        # Routes that name nothing plottable.
        for bad in [f"embeddings/{pca}/5.bin", f"embeddings/{pca}/0-0.bin", "embeddings/99/0.bin"]:
            route = parse_mapping_route(bad)
            assert route is not None, bad
            with pytest.raises(MappingNotFound):
                adapter.mapping.payload(route)
        for unparsed in [f"embeddings/{pca}/0-1-2-3.bin", "embeddings/x/0.bin", f"embeddings/{pca}/01.bin"]:
            assert parse_mapping_route(unparsed) is None, unparsed
    finally:
        adapter.close()


@pytest.mark.parametrize("method", ["GET", "HEAD"])
def test_the_mapping_embedding_route_serves_an_enabled_array(method: str) -> None:
    with _running_server() as (host, port):
        connection = http.client.HTTPConnection(host, port, timeout=5)
        try:
            connection.request("GET", f"/{MAPPING_CANDIDATES_PATH}")
            catalog = json.loads(connection.getresponse().read())
            index = [entry["key"] for entry in catalog["obsm"]].index("X_pca")
            connection.request(method, f"/_cellucid/mapping/embeddings/{index}/2-1.bin")
            response = connection.getresponse()
            body = response.read()
            assert response.status == 200
            assert response.getheader("Content-Type") == "application/octet-stream"
            assert int(response.getheader("Content-Length")) == 6 * 2 * 4
            assert body == (b"" if method == "HEAD" else body)
            if method == "GET":
                assert len(body) == 6 * 2 * 4
            connection.request("GET", f"/_cellucid/mapping/embeddings/{index}/9.bin")
            missing = connection.getresponse()
            missing.read()
            assert missing.status == 404
        finally:
            connection.close()


def test_vector_fields_are_bound_by_default_only_when_served():
    """``serve_vector_fields`` decides the defaults; the catalog lists the arrays either way."""
    adata = _adata()
    quiet = AnnDataAdapter(adata, dataset_name="d", dataset_id="d")
    served = AnnDataAdapter(adata, dataset_name="d", dataset_id="d", serve_vector_fields=True)
    try:
        assert quiet.get_mapping_candidates()["defaults"]["vectorFields"] == []
        assert served.get_mapping_candidates()["defaults"]["vectorFields"] == [
            {"key": "velocity_umap_3d", "embedding": {"key": "X_umap_3d", "columns": [0, 1, 2]}}
        ]
        assert "velocity_umap_3d" in [item["key"] for item in quiet.get_mapping_candidates()["obsm"]]
    finally:
        quiet.close()
        served.close()
