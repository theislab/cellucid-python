"""The ``_cellucid/mapping/`` routes of the AnnData server (design §3.1).

Every candidate array is addressed by its position in the frozen candidate
catalog, read only when it is requested, and refused with a reason when its
values cannot be served: 404 for an index the catalog does not have, 422 with
``{"error", "key"}`` for an entry whose values fail the load-time check. HEAD
answers lengths the catalog determines without reading the payload.
"""

from __future__ import annotations

import gzip
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

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.anndata_mapping import (
    AnnDataMapping,
    MappingLoadError,
    MappingNotFound,
    _ByteLRU,
    parse_mapping_route,
)
from cellucid.anndata_server import AnnDataServer
from cellucid.prepare_data._arrays import _normalize_finite_float32_embedding

N = 40


def _adata() -> ad.AnnData:
    rng = np.random.default_rng(7)
    adata = ad.AnnData(
        X=sparse.random(N, 6, density=0.5, random_state=0, format="csr", dtype=np.float32),
        obs=pd.DataFrame(
            {
                "cluster": pd.Categorical(["a", "b"] * (N // 2)),
                "score": np.linspace(0, 1, N, dtype=np.float32),
            },
            index=[f"cell-{i}" for i in range(N)],
        ),
        var=pd.DataFrame(
            {
                "symbol": ["MATR3", "MATR3", None, "CD4", "CD8A", "MS4A1"],
                "dispersion": np.arange(6, dtype=np.float64),
            },
            index=[f"ENSG{i:02d}" for i in range(6)],
        ),
    )
    adata.obsm["X_umap"] = rng.normal(size=(N, 2))
    adata.obsm["X_pca"] = rng.normal(size=(N, 10))
    adata.obsm["velocity_pca"] = rng.normal(size=(N, 10))
    adata.obsm["velocity_umap_2d"] = rng.normal(size=(N, 2))
    broken = rng.normal(size=(N, 2))
    broken[17, 1] = np.nan
    adata.obsm["broken"] = broken
    adata.obsm["X_table"] = pd.DataFrame({"a": np.arange(N)}, index=adata.obs_names)
    graph = sparse.random(N, N, density=0.1, random_state=1, format="csr")
    graph = graph + graph.T
    graph.setdiag(0)
    graph.eliminate_zeros()
    adata.obsp["connectivities"] = graph.tocsr()
    adata.obsp["distances"] = sparse.triu(
        sparse.random(N, N, density=0.1, random_state=2), k=1
    ).tocsr()
    adata.layers["log1p"] = adata.X.copy()
    adata.raw = adata.copy()
    return adata


@contextmanager
def _server(adata: ad.AnnData) -> Iterator[tuple[AnnDataServer, str, int]]:
    server = AnnDataServer(
        adata,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Mapping routes",
        dataset_id="mapping-routes",
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        yield server, host, port
    finally:
        server.stop()


def _request(
    host: str, port: int, method: str, path: str, headers: dict[str, str] | None = None
) -> tuple[int, dict[str, str], bytes]:
    connection = http.client.HTTPConnection(host, port, timeout=10)
    try:
        connection.request(method, f"/_cellucid/mapping/{path}", headers=headers or {})
        response = connection.getresponse()
        return response.status, dict(response.getheaders()), response.read()
    finally:
        connection.close()


@pytest.fixture(scope="module")
def served() -> Iterator[tuple[AnnDataServer, str, int, dict]]:
    with _server(_adata()) as (server, host, port):
        status, _headers, body = _request(host, port, "GET", "candidates.json")
        assert status == 200
        yield server, host, port, json.loads(body)


def _obsm(catalog: dict, key: str) -> int:
    return [item["key"] for item in catalog["obsm"]].index(key)


def _obsp(catalog: dict, key: str) -> int:
    return [item["key"] for item in catalog["obsp"]].index(key)


def test_the_candidate_catalog_lists_every_entry_and_the_served_defaults(served) -> None:
    server, _host, _port, catalog = served
    assert set(catalog) == {
        "kind",
        "n_obs",
        "obsm",
        "obsp",
        "var",
        "X",
        "layers",
        "raw",
        "defaults",
        "guessed",
    }
    assert catalog["n_obs"] == N
    keys = [item["key"] for item in catalog["obsm"]]
    assert keys == list(server.adapter.adata.obsm.keys())
    table = catalog["obsm"][_obsm(catalog, "X_table")]
    assert (table["storage"], table["usable"], table["reason"]) == (
        "dataframe",
        False,
        "obsm table, not an array",
    )
    assert catalog["defaults"]["embeddings"]["defaults"]["2d"] == {"key": "X_umap", "columns": [0, 1]}
    assert catalog["defaults"]["knn"] is None
    assert catalog["defaults"]["geneNames"] == "index"
    assert catalog["defaults"]["expression"] == "X"
    assert catalog["guessed"] == []
    assert [layer["key"] for layer in catalog["layers"]] == ["log1p"]
    assert catalog["raw"]["n_vars"] == 6


def test_a_mapped_embedding_is_the_normalized_float32_of_its_columns(served) -> None:
    server, host, port, catalog = served
    index = _obsm(catalog, "X_pca")
    status, headers, body = _request(host, port, "GET", f"embeddings/{index}/3-1.bin")
    assert status == 200
    assert headers["Content-Type"] == "application/octet-stream"
    expected, *_ = _normalize_finite_float32_embedding(
        np.asarray(server.adapter.adata.obsm["X_pca"])[:, [3, 1]].astype(np.float64),
        label="X_pca",
    )
    assert body == expected.astype("<f4").tobytes()

    status, headers, head_body = _request(host, port, "HEAD", f"embeddings/{index}/3-1.bin")
    assert (status, int(headers["Content-Length"]), head_body) == (200, N * 2 * 4, b"")

    status, headers, compressed = _request(
        host, port, "GET", f"embeddings/{index}/3-1.bin", {"Accept-Encoding": "gzip"}
    )
    assert headers.get("Content-Encoding") == "gzip"
    assert gzip.decompress(compressed) == body


@pytest.mark.parametrize(
    "path",
    [
        "embeddings/99/0-1.bin",
        "embeddings/1/1-1.bin",
        "embeddings/1/0-10.bin",
        "embeddings/01/0-1.bin",
        "embeddings/1/0-1-2-3.bin",
        "knn/9/connectivity_manifest.json",
        "genes/layer-3/index/var_manifest.json",
        "genes/X/column-9/var_manifest.json",
        "genes/X/index/var/6.values.f32",
        "latent/1/obs/1.outliers.f32",
        "unknown.json",
    ],
)
def test_an_unknown_index_or_route_is_404(served, path: str) -> None:
    _server_obj, host, port, catalog = served
    assert _obsm(catalog, "X_pca") == 1
    for method in ("GET", "HEAD"):
        status, _headers, _body = _request(host, port, method, path)
        assert status == 404, (method, path)


def test_load_time_failures_are_422_with_the_reason_and_the_key(served) -> None:
    _server_obj, host, port, catalog = served
    # `reason` is the clause a viewer lists beside the key, in the words a
    # file read in the browser gives; `error` is the full sentence.
    cases = {
        f"embeddings/{_obsm(catalog, 'broken')}/0-1.bin": ("broken", "NaN at cell 17"),
        f"embeddings/{_obsm(catalog, 'X_table')}/0.bin": ("X_table", "obsm table, not an array"),
        f"knn/{_obsp(catalog, 'distances')}/connectivity_manifest.json": (
            "distances",
            "not symmetric — choose a connectivities matrix",
        ),
        "genes/X/column-1/var_manifest.json": (
            "dispersion",
            "floating-point values cannot name genes",
        ),
    }
    for path, (key, reason) in cases.items():
        status, headers, body = _request(host, port, "GET", path)
        assert status == 422, path
        assert headers["Content-Type"] == "application/json"
        payload = json.loads(body)
        assert set(payload) == {"error", "key", "reason"}
        assert payload["key"] == key
        assert payload["reason"] == reason
    status, _headers, body = _request(
        host, port, "GET", f"knn/{_obsp(catalog, 'distances')}/connectivity_manifest.json"
    )
    assert json.loads(body)["error"] == (
        "obsp['distances'] is not symmetric — choose a connectivities matrix."
    )


def test_a_knn_graph_is_read_only_when_requested_and_served_as_edge_pairs(served) -> None:
    server, host, port, catalog = served
    index = _obsp(catalog, "connectivities")
    status, _headers, body = _request(host, port, "GET", f"knn/{index}/connectivity_manifest.json")
    assert status == 200
    manifest = json.loads(body)
    assert manifest["format"] == "edge_pairs"
    assert manifest["sourcesPath"] == "connectivity/edges.src.bin"
    for name, key in (
        ("edges.src.bin", "index_bytes"),
        ("edges.dst.bin", "index_bytes"),
        ("edges.weights.f64.bin", "weight_bytes"),
    ):
        status, headers, _ = _request(host, port, "HEAD", f"knn/{index}/connectivity/{name}")
        assert int(headers["Content-Length"]) == manifest["n_edges"] * manifest[key]
        status, _, payload = _request(host, port, "GET", f"knn/{index}/connectivity/{name}")
        assert len(payload) == manifest["n_edges"] * manifest[key]
    # The served identity still has no graph: the mapping route is separate.
    assert server.adapter.get_dataset_identity()["stats"]["has_connectivity"] is False


def test_gene_routes_use_the_axis_and_the_made_unique_name_set(served) -> None:
    server, host, port, _catalog = served
    status, _headers, body = _request(host, port, "GET", "genes/X/column-0/var_manifest.json")
    manifest = json.loads(body)
    assert [field[1] for field in manifest["fields"]] == [
        "MATR3",
        "MATR3-1",
        "ENSG02",
        "CD4",
        "CD8A",
        "MS4A1",
    ]
    assert manifest["gene_name_sets"] == [
        {
            "source": "symbol",
            "filled_from_index": 1,
            "renamed": [["MATR3", "MATR3-1"]],
            "names": None,
        }
    ]
    assert manifest["_format"] == "compact_v2"
    adata = server.adapter.adata
    for axis, matrix in (("X", adata.X), ("layer-0", adata.layers["log1p"]), ("raw", adata.raw.X)):
        status, headers, payload = _request(host, port, "GET", f"genes/{axis}/index/var/3.values.f32")
        assert status == 200
        expected = np.asarray(matrix[:, 3].toarray()).reshape(-1).astype("<f4").tobytes()
        assert payload == expected
        status, headers, _ = _request(host, port, "HEAD", f"genes/{axis}/index/var/3.values.f32")
        assert int(headers["Content-Length"]) == N * 4


def test_vectors_bind_to_a_mapped_embedding_and_use_its_scale(served) -> None:
    server, host, port, catalog = served
    adata = server.adapter.adata
    pca, velocity = _obsm(catalog, "X_pca"), _obsm(catalog, "velocity_pca")
    status, _headers, body = _request(host, port, "GET", f"vectors/{velocity}/3-1/{pca}/3-1.bin")
    assert status == 200
    _, _center, scale, _range = _normalize_finite_float32_embedding(
        np.asarray(adata.obsm["X_pca"])[:, [3, 1]].astype(np.float64), label="X_pca"
    )
    expected = (np.asarray(adata.obsm["velocity_pca"])[:, [3, 1]] * scale).astype("<f4")
    assert body == expected.tobytes()
    umap, vumap = _obsm(catalog, "X_umap"), _obsm(catalog, "velocity_umap_2d")
    assert _request(host, port, "GET", f"vectors/{vumap}/0-1/{umap}/0-1.bin")[0] == 200
    # Bound at other columns than the embedding's: refused, not reinterpreted.
    assert _request(host, port, "GET", f"vectors/{velocity}/0-1/{pca}/3-1.bin")[0] == 404
    assert _request(host, port, "HEAD", f"vectors/{velocity}/0-1/{pca}/3-1.bin")[0] == 404


def test_latent_outliers_follow_the_served_obs_fields(served) -> None:
    server, host, port, catalog = served
    pca = _obsm(catalog, "X_pca")
    cluster_index = str(server.adapter.get_obs_payload_index("cluster"))
    # Centroids of a mapped embedding are the viewer's to compute from its
    # positions; the server has no route for them.
    assert _request(host, port, "GET", f"embeddings/{pca}/0-1/centroids.json")[0] == 404
    status, headers, body = _request(host, port, "GET", f"latent/{pca}/obs/{cluster_index}.outliers.f32")
    assert status == 200 and len(body) == N * 4
    quantiles = np.frombuffer(body, dtype="<f4")
    assert np.all((quantiles > 0) & (quantiles <= 1))


def test_concurrent_requests_get_identical_bytes(served) -> None:
    server, host, port, catalog = served
    paths = [
        f"embeddings/{_obsm(catalog, 'X_pca')}/{a}-{b}.bin"
        for a, b in ((0, 1), (2, 5), (9, 4))
    ] + [f"genes/raw/index/var/{gene}.values.f32" for gene in range(6)]
    reference = {path: _request(host, port, "GET", path)[2] for path in paths}
    results: list[tuple[str, bytes]] = []
    errors: list[BaseException] = []

    def worker(offset: int) -> None:
        try:
            for path in paths[offset:] + paths[:offset]:
                results.append((path, _request(host, port, "GET", path)[2]))
        except BaseException as error:  # pragma: no cover - surfaced below
            errors.append(error)

    threads = [threading.Thread(target=worker, args=(offset % len(paths),)) for offset in range(12)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=60)
    assert not errors
    assert len(results) == 12 * len(paths)
    assert all(body == reference[path] for path, body in results)


def test_the_service_itself_is_safe_to_call_from_many_threads() -> None:
    """The kernel thread and the server thread share one adapter."""
    adapter = AnnDataAdapter(_adata(), dataset_name="Threads", dataset_id="threads")
    try:
        mapping = adapter.mapping
        index = [item["key"] for item in mapping.candidates["obsm"]].index("X_pca")
        expected = mapping.embedding_bytes(index, (0, 1))
        mapping.close()
        barrier = threading.Barrier(16)
        outputs: list[bytes] = []

        def call() -> None:
            barrier.wait()
            outputs.append(mapping.embedding_bytes(index, (0, 1)))
            outputs.append(mapping.gene_bytes("X", "column-0", 2))

        threads = [threading.Thread(target=call) for _ in range(16)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=60)
        assert outputs.count(expected) == 16
        assert len(outputs) == 32
    finally:
        adapter.close()


def test_the_payload_cache_is_bounded_by_bytes() -> None:
    cache = _ByteLRU(10)
    cache.put("a", b"12345", 5)
    cache.put("b", b"12345", 5)
    assert cache.get("a") == b"12345"
    cache.put("c", b"12345", 5)
    assert cache.get("b") is None and cache.get("a") == b"12345"
    cache.put("huge", b"x" * 20, 20)
    assert cache.get("huge") is None and cache.resident_bytes == 10


def test_route_grammar_addresses_entries_by_position_only() -> None:
    assert parse_mapping_route("embeddings/2/0-1.bin").parts == (2, (0, 1))
    assert parse_mapping_route("vectors/3/1-0/2/1-0.bin").parts == (3, (1, 0), 2, (1, 0))
    assert parse_mapping_route("genes/layer-12/column-3/var/9.values.f32").parts == (
        "layer-12",
        "column-3",
        9,
    )
    for path in (
        "embeddings/X_umap/0-1.bin",
        "genes/layer-counts/index/var_manifest.json",
        "knn/connectivities/connectivity_manifest.json",
        "embeddings/-1/0-1.bin",
    ):
        assert parse_mapping_route(path) is None


def test_adapter_get_mapped_methods_mirror_the_routes() -> None:
    adapter = AnnDataAdapter(_adata(), dataset_name="Methods", dataset_id="methods")
    try:
        catalog = adapter.get_mapping_candidates()
        assert catalog == adapter.mapping.candidates and catalog is not adapter.mapping.candidates
        pca = [item["key"] for item in catalog["obsm"]].index("X_pca")
        assert len(adapter.get_mapped_embedding(pca, [0, 1])) == N * 8
        velocity = [item["key"] for item in catalog["obsm"]].index("velocity_pca")
        assert len(adapter.get_mapped_vectors(velocity, [0, 1], pca, [0, 1])) == N * 8
        graph = [item["key"] for item in catalog["obsp"]].index("connectivities")
        manifest = adapter.get_mapped_connectivity(graph)
        sources, destinations, weights = adapter.get_mapped_connectivity_edges(graph)
        assert len(weights) == manifest["n_edges"] * 8
        assert len(sources) == len(destinations) == manifest["n_edges"] * manifest["index_bytes"]
        assert len(adapter.get_mapped_var_manifest("raw", "column-0")["fields"]) == 6
        assert len(adapter.get_mapped_gene_values("layer-0", "index", 1)) == N * 4
        assert len(adapter.get_mapped_outliers(pca, 0)) == N * 4
        with pytest.raises(MappingNotFound):
            adapter.get_mapped_embedding(99, [0, 1])
        broken = [item["key"] for item in catalog["obsm"]].index("broken")
        with pytest.raises(MappingLoadError) as raised:
            adapter.get_mapped_embedding(broken, [0, 1])
        assert raised.value.key == "broken"
        assert raised.value.reason == "NaN at cell 17"
        assert str(raised.value) == "Embedding 'broken' cannot be plotted: NaN at cell 17."
    finally:
        adapter.close()
    with pytest.raises(RuntimeError):
        adapter.mapping  # noqa: B018


def test_latent_outliers_use_obsm_without_a_cached_copy() -> None:
    # The latent space is read from obsm as it is; there is no byte budget,
    # and only the per-cell quantiles are kept, never a copy of the latent.
    adata = _adata()
    adapter = AnnDataAdapter(adata, dataset_name="Latent", dataset_id="latent")
    try:
        pca = [item["key"] for item in adapter.mapping.candidates["obsm"]].index("X_pca")
        assert len(adapter.get_mapped_outliers(pca, 0)) == N * 4
        assert adapter.mapping.cache_bytes == N * 4
    finally:
        adapter.close()


def test_the_mapping_service_is_one_per_adapter() -> None:
    adapter = AnnDataAdapter(_adata(), dataset_name="One", dataset_id="one")
    try:
        assert isinstance(adapter.mapping, AnnDataMapping)
        assert adapter.mapping is adapter.mapping
    finally:
        adapter.close()


def test_the_catalog_lists_an_array_the_best_guess_refused_as_unusable() -> None:
    # No UMAP: the guess ranks X_umap_old first, loads it, finds the NaN and
    # moves on to X_tsne. The catalog then says why X_umap_old cannot be
    # chosen -- the same verdict and text the browser's catalog carries -- so
    # no viewer offers a choice this server has already refused.
    rng = np.random.default_rng(11)
    adata = ad.AnnData(
        X=np.ones((N, 2), dtype=np.float32),
        obs=pd.DataFrame(index=[f"cell-{i}" for i in range(N)]),
    )
    umap_old = rng.normal(size=(N, 2))
    umap_old[17, 1] = np.nan
    adata.obsm["X_umap_old"] = umap_old
    adata.obsm["X_tsne"] = rng.normal(size=(N, 2))
    adapter = AnnDataAdapter(adata, dataset_name="guess", dataset_id="guess")
    catalog = adapter.get_mapping_candidates()
    by_key = {item["key"]: item for item in catalog["obsm"]}
    assert by_key["X_umap_old"]["usable"] is False
    assert by_key["X_umap_old"]["reason"] == "NaN at cell 17"
    assert by_key["X_tsne"]["usable"] is True
    assert catalog["defaults"]["embeddings"]["defaults"]["2d"] == {"key": "X_tsne", "columns": [0, 1]}
    assert catalog["guessed"] == ["embeddings.2d"]
    index = [item["key"] for item in catalog["obsm"]].index("X_umap_old")
    with pytest.raises(MappingLoadError) as refused:
        adapter.mapping.embedding_bytes(index, (0, 1))
    assert refused.value.key == "X_umap_old"
