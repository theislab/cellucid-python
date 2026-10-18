"""What a served AnnData object publishes under dataset identity version 3.

The server fills one embedding per dimension slot: the default rule
(``X_umap_<d>d``, else a 1- to 3-column ``X_umap``), any slot named explicitly
with ``embedding_keys``, and -- only when neither fills a slot -- the shared
deterministic best guess. Each served embedding is an identity item keyed by
its canonical id and written to ``embeddings/{index}_{dimension}d.bin``.
"""

from __future__ import annotations

import json
import urllib.request
from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any
from unittest import mock

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.anndata_server import AnnDataServer

N_CELLS = 6


def _adata(**obsm: Any) -> ad.AnnData:
    rng = np.random.default_rng(23)
    return ad.AnnData(
        X=rng.random((N_CELLS, 3), dtype=np.float32),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(["a", "b"] * (N_CELLS // 2))},
            index=[f"cell-{index}" for index in range(N_CELLS)],
        ),
        var=pd.DataFrame(
            {"symbol": ["MATR3", "MATR3", None]},
            index=["ENSG01", "ENSG02", "ENSG03"],
        ),
        obsm=obsm,
    )


def _adapter(adata: ad.AnnData, **kwargs: Any) -> AnnDataAdapter:
    return AnnDataAdapter(
        adata,
        dataset_name="Identity v3",
        dataset_id="identity-v3",
        centroid_min_points=1,
        **kwargs,
    )


def _normalized_float32(values: np.ndarray) -> np.ndarray:
    """The one normalization the format defines, rounded to float32 once."""
    source = np.asarray(values, dtype=np.float64)
    mins = source.min(axis=0)
    maxs = source.max(axis=0)
    scale = 2.0 / float((maxs - mins).max())
    return ((source - (mins + maxs) / 2) * scale).astype(np.float32)


def _ring(n: int, width: int, *, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).normal(size=(n, width))


# ---------------------------------------------------------------------------
# The default rule, byte for byte
# ---------------------------------------------------------------------------


def test_dimensional_umap_keys_publish_one_item_each_with_the_same_bytes() -> None:
    umap2 = _ring(N_CELLS, 2, seed=1)
    umap3 = _ring(N_CELLS, 3, seed=2)
    adapter = _adapter(_adata(X_umap_2d=umap2, X_umap_3d=umap3, X_tsne=_ring(N_CELLS, 2, seed=3)))
    try:
        embeddings = adapter.get_dataset_identity()["embeddings"]
        assert embeddings == {
            "default_dimension": 3,
            "defaults": {"2d": "X_umap_2d", "3d": "X_umap_3d"},
            "items": {
                "X_umap_2d": {
                    "label": "X_umap_2d",
                    "dimension": 2,
                    "file": "embeddings/0_2d.bin",
                    "source": {"key": "X_umap_2d", "columns": [0, 1]},
                },
                "X_umap_3d": {
                    "label": "X_umap_3d",
                    "dimension": 3,
                    "file": "embeddings/1_3d.bin",
                    "source": {"key": "X_umap_3d", "columns": [0, 1, 2]},
                },
            },
        }
        # The coordinates are exactly the ones version 2 served: only the
        # path changed, so a saved cell-order digest still matches.
        assert adapter.get_embedding_binary("X_umap_2d") == _normalized_float32(umap2).astype(
            "<f4"
        ).tobytes()
        assert adapter.get_embedding_binary("X_umap_3d") == _normalized_float32(umap3).astype(
            "<f4"
        ).tobytes()
    finally:
        adapter.close()


def test_a_bare_umap_fills_the_slot_of_its_width_and_keeps_its_bytes() -> None:
    umap = _ring(N_CELLS, 2, seed=4)
    adapter = _adapter(_adata(X_umap=umap))
    try:
        embeddings = adapter.get_dataset_identity()["embeddings"]
        assert embeddings["defaults"] == {"2d": "X_umap"}
        assert embeddings["items"]["X_umap"]["file"] == "embeddings/0_2d.bin"
        assert adapter.get_embedding_binary("X_umap") == _normalized_float32(umap).tobytes()
    finally:
        adapter.close()


def test_a_dimensional_umap_of_the_wrong_width_is_refused_not_sliced() -> None:
    with pytest.raises(ValueError, match=r"'X_umap_2d' must have exactly 2 columns; got 3"):
        _adapter(_adata(X_umap_2d=_ring(N_CELLS, 3, seed=5)))


def test_categorical_centroids_are_keyed_by_served_embedding_id() -> None:
    adapter = _adapter(
        _adata(X_umap_2d=_ring(N_CELLS, 2, seed=6), X_umap_3d=_ring(N_CELLS, 3, seed=7))
    )
    try:
        manifest = adapter.get_obs_manifest()
        assert manifest["_format"] == "compact_v2"
        centroids = manifest["_categoricalFields"][0][5]
        assert list(centroids) == ["X_umap_2d", "X_umap_3d"]
        assert all(len(entry["position"]) == 2 for entry in centroids["X_umap_2d"])
        assert all(len(entry["position"]) == 3 for entry in centroids["X_umap_3d"])
    finally:
        adapter.close()


# ---------------------------------------------------------------------------
# Explicit slots
# ---------------------------------------------------------------------------


def test_embedding_keys_override_their_slot_and_leave_the_others_to_the_default() -> None:
    pca = _ring(N_CELLS, 6, seed=8)
    adapter = _adapter(
        _adata(X_umap_2d=_ring(N_CELLS, 2, seed=9), X_umap_3d=_ring(N_CELLS, 3, seed=10), X_pca=pca),
        embedding_keys={"2d": ("X_pca", [3, 1])},
    )
    try:
        embeddings = adapter.get_dataset_identity()["embeddings"]
        assert embeddings["defaults"] == {"2d": "X_pca[3,1]", "3d": "X_umap_3d"}
        assert embeddings["items"]["X_pca[3,1]"] == {
            "label": "X_pca[3,1]",
            "dimension": 2,
            "file": "embeddings/0_2d.bin",
            "source": {"key": "X_pca", "columns": [3, 1]},
        }
        assert adapter.get_embedding_binary("X_pca[3,1]") == _normalized_float32(
            pca[:, [3, 1]]
        ).tobytes()
    finally:
        adapter.close()


def test_a_bare_key_reads_the_leading_columns_and_a_whole_array_keeps_its_key() -> None:
    adapter = _adapter(
        _adata(X_tsne=_ring(N_CELLS, 2, seed=11), X_pca=_ring(N_CELLS, 5, seed=12)),
        embedding_keys={"2d": "X_tsne", "3d": "X_pca"},
    )
    try:
        assert adapter.get_dataset_identity()["embeddings"]["defaults"] == {
            "2d": "X_tsne",
            "3d": "X_pca[0,1,2]",
        }
    finally:
        adapter.close()


@pytest.mark.parametrize(
    ("embedding_keys", "error", "match"),
    [
        ("X_tsne", TypeError, "embedding_keys must map"),
        ({"4d": "X_tsne"}, ValueError, "'4d' must be '1d', '2d', or '3d'"),
        ({"2d": 7}, TypeError, "obsm key or a \\(key, columns\\) pair"),
        ({"2d": ["X_pca", [0, 1]]}, TypeError, "obsm key or a \\(key, columns\\) pair"),
        ({"2d": ("X_pca", "0,1")}, TypeError, "columns must be a sequence"),
        ({"2d": ""}, ValueError, "non-empty obsm key"),
        ({"2d": "X_missing"}, ValueError, "obsm\\['X_missing'\\], which this object does not have"),
        ({"2d": ("X_pca", [0, 5])}, ValueError, "column 5 does not exist"),
        ({"2d": ("X_pca", [1, 1])}, ValueError, "Each axis needs a different column"),
        ({"2d": ("X_pca", [0, 1, 2])}, ValueError, "exactly 2 columns"),
        ({"2d": "X_line"}, ValueError, "at least 2 columns"),
        ({"2d": "X_sparse"}, ValueError, "sparse obsm arrays cannot be plotted"),
        ({"2d": "X_table"}, ValueError, "obsm table, not an array"),
    ],
)
def test_embedding_keys_are_validated_before_anything_is_read(
    embedding_keys: object,
    error: type[Exception],
    match: str,
) -> None:
    adata = _adata(
        X_umap=_ring(N_CELLS, 2, seed=13),
        X_tsne=_ring(N_CELLS, 2, seed=14),
        X_pca=_ring(N_CELLS, 5, seed=15),
        X_line=np.arange(N_CELLS, dtype=np.float32),
        X_sparse=sparse.csr_matrix(_ring(N_CELLS, 2, seed=16)),
    )
    adata.obsm["X_table"] = pd.DataFrame({"a": np.arange(N_CELLS)}, index=adata.obs_names)
    with pytest.raises(error, match=match):
        _adapter(adata, embedding_keys=embedding_keys)


def test_an_explicit_slot_that_cannot_be_plotted_is_reported_not_skipped() -> None:
    values = _ring(N_CELLS, 2, seed=17)
    values[2, 1] = np.nan
    with pytest.raises(ValueError, match="X_tsne"):
        _adapter(_adata(X_umap=_ring(N_CELLS, 2, seed=18), X_tsne=values), embedding_keys={"2d": "X_tsne"})


# ---------------------------------------------------------------------------
# The best guess
# ---------------------------------------------------------------------------


def test_the_best_guess_skips_a_candidate_that_fails_to_load(capsys) -> None:
    tsne = _ring(N_CELLS, 2, seed=19)
    tsne[4, 0] = np.nan
    spatial = _ring(N_CELLS, 2, seed=20)
    adapter = _adapter(_adata(X_tsne=tsne, spatial=spatial), quiet=False)
    try:
        assert adapter.get_dataset_identity()["embeddings"]["defaults"] == {"2d": "spatial"}
        assert adapter.get_embedding_binary("spatial") == _normalized_float32(spatial).tobytes()
        assert "Embeddings: 2D from obsm['spatial'] (guessed)" in capsys.readouterr().out
    finally:
        adapter.close()


def test_the_best_guess_takes_a_wide_array_at_its_first_two_columns() -> None:
    pca = _ring(N_CELLS, 30, seed=21)
    adapter = _adapter(_adata(X_pca=pca, X_scvi=_ring(N_CELLS, 10, seed=22)), latent_key="X_scvi")
    try:
        embeddings = adapter.get_dataset_identity()["embeddings"]
        # X_scvi is the latent space, so it is not a candidate at all.
        assert embeddings["defaults"] == {"2d": "X_pca[0,1]"}
        assert embeddings["items"]["X_pca[0,1]"]["source"] == {"key": "X_pca", "columns": [0, 1]}
    finally:
        adapter.close()


def test_nothing_plottable_names_every_obsm_entry_and_why() -> None:
    spatial = _ring(N_CELLS, 2, seed=23)
    spatial[3, 0] = np.inf
    adata = _adata(spatial=spatial, velocity_umap_2d=_ring(N_CELLS, 2, seed=24))
    adata.obsm["X_table"] = pd.DataFrame({"a": np.arange(N_CELLS)}, index=adata.obs_names)
    adata.obsm["X_sparse"] = sparse.csr_matrix(_ring(N_CELLS, 2, seed=25))

    with pytest.raises(ValueError) as error:
        _adapter(adata)

    assert str(error.value) == (
        "Cellucid found nothing to plot. It needs an obsm array with one row per cell (6) "
        "and at least one numeric column. obsm holds: spatial (6 × 2, Infinity at cell 3), "
        "velocity_umap_2d (6 × 2, a vector field), X_table (table), X_sparse (sparse 6 × 2)."
    )


def test_an_empty_obsm_says_it_holds_nothing() -> None:
    with pytest.raises(ValueError, match=r"obsm holds: nothing\.$"):
        _adapter(_adata())


# ---------------------------------------------------------------------------
# Vector fields follow the default binding rule
# ---------------------------------------------------------------------------


def test_umap_vectors_bind_to_umap_slots_and_are_routed_by_index_pair() -> None:
    adata = _adata(
        X_umap_2d=_ring(N_CELLS, 2, seed=26),
        X_umap_3d=_ring(N_CELLS, 3, seed=27),
        velocity_umap_2d=_ring(N_CELLS, 2, seed=28),
        velocity_umap_3d=_ring(N_CELLS, 3, seed=29),
    )
    adapter = _adapter(adata, serve_vector_fields=True)
    try:
        assert adapter.get_dataset_identity()["vector_fields"] == {
            "default_field": "velocity_umap",
            "fields": {
                "velocity_umap": {
                    "label": "velocity_umap",
                    "files": {
                        "X_umap_2d": "vectors/0_0.bin",
                        "X_umap_3d": "vectors/0_1.bin",
                    },
                }
            },
        }
        assert adapter.get_vector_field_route("0_1.bin") == ("velocity_umap", "X_umap_3d")
        assert adapter.get_vector_field_route("0_2.bin") is None
        assert adapter.get_vector_field_route("1_0.bin") is None
        assert adapter.get_vector_field_route("0_1d.bin") is None
    finally:
        adapter.close()


def test_a_umap_vector_does_not_bind_to_a_slot_showing_another_layout() -> None:
    adata = _adata(
        X_umap=_ring(N_CELLS, 2, seed=30),
        X_tsne=_ring(N_CELLS, 2, seed=31),
        velocity_umap_2d=_ring(N_CELLS, 2, seed=32),
    )
    adapter = _adapter(adata, serve_vector_fields=True, embedding_keys={"2d": "X_tsne"})
    try:
        assert "vector_fields" not in adapter.get_dataset_identity()
    finally:
        adapter.close()


# ---------------------------------------------------------------------------
# The neighbour graph is named, read when named, and never symmetrized
# ---------------------------------------------------------------------------


def _graphs() -> dict[str, sparse.csr_matrix]:
    symmetric = np.zeros((N_CELLS, N_CELLS))
    for index in range(N_CELLS - 1):
        symmetric[index, index + 1] = symmetric[index + 1, index] = 0.5
    directed = np.zeros((N_CELLS, N_CELLS))
    for index in range(N_CELLS - 1):
        directed[index, index + 1] = 1.0
    return {
        "connectivities": sparse.csr_matrix(symmetric),
        "distances": sparse.csr_matrix(directed),
    }


def test_a_non_symmetric_graph_is_refused_with_the_reason() -> None:
    adata = _adata(X_umap=_ring(N_CELLS, 2, seed=33))
    for key, graph in _graphs().items():
        adata.obsp[key] = graph
    with pytest.raises(ValueError) as error:
        _adapter(adata, knn_key="distances")
    assert str(error.value).startswith(
        "obsp['distances'] is not symmetric — choose a connectivities matrix."
    )


def test_a_named_graph_is_served_and_an_absent_one_is_an_error() -> None:
    adata = _adata(X_umap=_ring(N_CELLS, 2, seed=34))
    for key, graph in _graphs().items():
        adata.obsp[key] = graph
    adapter = _adapter(adata, knn_key="connectivities")
    try:
        assert adapter.get_dataset_identity()["stats"]["n_edges"] == N_CELLS - 1
    finally:
        adapter.close()
    with pytest.raises(ValueError, match=r"knn_key names obsp\['knn'\].*obsp holds"):
        _adapter(adata, knn_key="knn")
    with pytest.raises(TypeError, match="knn_key"):
        _adapter(adata, knn_key=True)


# ---------------------------------------------------------------------------
# Expression source and gene names
# ---------------------------------------------------------------------------


def test_duplicate_and_empty_gene_names_are_made_unique_on_the_served_axis() -> None:
    adata = _adata(X_umap=_ring(N_CELLS, 2, seed=35))
    adapter = _adapter(adata, gene_id_column="symbol")
    try:
        manifest = adapter.get_var_manifest()
        assert manifest["fields"] == [[0, "MATR3"], [1, "MATR3-1"], [2, "ENSG03"]]
        assert manifest["gene_name_sets"] == [
            {
                "source": "symbol",
                "filled_from_index": 1,
                "renamed": [["MATR3", "MATR3-1"]],
                "names": None,
            }
        ]
        second = np.frombuffer(adapter.get_gene_expression("MATR3-1"), dtype="<f4")
        np.testing.assert_array_equal(second, np.asarray(adata.X[:, 1], dtype=np.float32))
        assert adapter.get_gene_id_for_payload_index("1") == "MATR3-1"
    finally:
        adapter.close()


def test_raw_expression_reads_raw_x_and_names_genes_from_raw_var() -> None:
    adata = _adata(X_umap=_ring(N_CELLS, 2, seed=36))
    raw = ad.AnnData(
        X=np.arange(N_CELLS * 4, dtype=np.float32).reshape(N_CELLS, 4),
        obs=adata.obs[[]],
        var=pd.DataFrame({"name": ["A", "B", "C", "D"]}, index=["r1", "r2", "r3", "r4"]),
    )
    adata.raw = raw
    adapter = _adapter(adata, expression_source="raw", gene_id_column="name")
    try:
        assert adapter.n_genes == 4
        assert adapter.get_dataset_identity()["stats"]["n_genes"] == 4
        assert [field[1] for field in adapter.get_var_manifest()["fields"]] == ["A", "B", "C", "D"]
        np.testing.assert_array_equal(
            np.frombuffer(adapter.get_gene_expression("C"), dtype="<f4"),
            raw.X[:, 2],
        )
    finally:
        adapter.close()


def test_a_layer_is_read_by_name_and_an_absent_one_is_named() -> None:
    adata = _adata(X_umap=_ring(N_CELLS, 2, seed=37))
    adata.layers["counts"] = np.full((N_CELLS, 3), 7.0, dtype=np.float32)
    adapter = _adapter(adata, expression_source="layer:counts")
    try:
        np.testing.assert_array_equal(
            np.frombuffer(adapter.get_gene_expression("ENSG01"), dtype="<f4"),
            np.full(N_CELLS, 7.0, dtype=np.float32),
        )
    finally:
        adapter.close()
    with pytest.raises(ValueError, match=r"'layer:log1p' names no layer.*\['counts'\]"):
        _adapter(adata, expression_source="layer:log1p")
    with pytest.raises(ValueError, match="raw"):
        _adapter(adata, expression_source="raw")


@pytest.mark.parametrize("value", ["x", "layers:counts", "layer:", "", None, 1])
def test_expression_source_is_one_exact_spelling(value: object) -> None:
    with pytest.raises(ValueError, match="expression_source must be 'X', 'raw', or 'layer:<name>'"):
        _adapter(_adata(X_umap=_ring(N_CELLS, 2, seed=38)), expression_source=value)


# ---------------------------------------------------------------------------
# Over HTTP, from a backed file
# ---------------------------------------------------------------------------


@contextmanager
def _running(adata_or_path: Any, **kwargs: Any) -> Iterator[tuple[AnnDataServer, str]]:
    server = AnnDataServer(
        adata_or_path,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Identity v3",
        dataset_id="identity-v3",
        centroid_min_points=1,
        **kwargs,
    )
    server.start_background()
    try:
        assert server._server is not None
        host, port = server._server.server_address
        yield server, f"http://{host}:{port}"
    finally:
        server.stop()


def _head_length(url: str) -> int:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=5) as response:
        return int(response.headers["Content-Length"])


def test_a_backed_file_serves_explicit_slots_and_vector_routes_with_exact_lengths(
    tmp_path,
) -> None:
    adata = _adata(
        X_umap_2d=_ring(N_CELLS, 2, seed=39),
        X_pca=_ring(N_CELLS, 5, seed=40),
        velocity_umap_2d=_ring(N_CELLS, 2, seed=41),
    )
    path = tmp_path / "fixture.h5ad"
    adata.write_h5ad(path)

    with _running(
        path,
        embedding_keys={"3d": ("X_pca", [4, 0, 2])},
        serve_vector_fields=True,
        gene_id_column="symbol",
    ) as (server, base):
        with urllib.request.urlopen(f"{base}/dataset_identity.json", timeout=5) as response:
            identity = json.loads(response.read())
        assert identity["version"] == 3
        assert identity["embeddings"]["defaults"] == {"2d": "X_umap_2d", "3d": "X_pca[4,0,2]"}
        assert _head_length(f"{base}/embeddings/0_2d.bin") == N_CELLS * 2 * 4
        assert _head_length(f"{base}/embeddings/1_3d.bin") == N_CELLS * 3 * 4
        assert _head_length(f"{base}/vectors/0_0.bin") == N_CELLS * 2 * 4
        with urllib.request.urlopen(f"{base}/embeddings/1_3d.bin", timeout=5) as response:
            assert response.read() == _normalized_float32(
                np.asarray(adata.obsm["X_pca"])[:, [4, 0, 2]]
            ).tobytes()
        with urllib.request.urlopen(f"{base}/var_manifest.json", timeout=5) as response:
            assert [field[1] for field in json.loads(response.read())["fields"]] == [
                "MATR3",
                "MATR3-1",
                "ENSG03",
            ]
        for missing in ("/embeddings/2_3d.bin", "/embeddings/1_2d.bin", "/vectors/0_1.bin"):
            with pytest.raises(urllib.error.HTTPError) as error:
                urllib.request.urlopen(f"{base}{missing}", timeout=5)
            assert error.value.code == 404


def test_the_dataset_report_names_each_slot_and_the_graph_flag(capsys) -> None:
    adata = _adata(X_umap=_ring(N_CELLS, 2, seed=42), X_pca=_ring(N_CELLS, 5, seed=43))
    adata.obsp["connectivities"] = _graphs()["connectivities"]
    server = AnnDataServer(
        adata,
        host="127.0.0.1",
        port=0,
        quiet=False,
        serve_web_ui=False,
        dataset_name="Report",
        dataset_id="report",
        embedding_keys={"3d": "X_pca"},
    )
    try:
        printed = capsys.readouterr().out
        assert "Embeddings: 2D X_umap, 3D X_pca[0,1,2]" in printed
        assert "pass knn_key='connectivities', or --knn connectivities" in printed
    finally:
        server.stop()


# ---------------------------------------------------------------------------
# The notebook entry points hand every mapping argument to the server
# ---------------------------------------------------------------------------


def test_the_notebook_viewer_forwards_every_mapping_argument() -> None:
    from cellucid.jupyter import _anndata

    with (
        mock.patch("cellucid.anndata_server.AnnDataServer") as server_class,
        mock.patch.object(_anndata.AnnDataViewer, "__init__", return_value=None),
    ):
        viewer = _anndata.AnnDataViewer.__new__(_anndata.AnnDataViewer)
        viewer.data = object()
        viewer.port = 0
        viewer.web_source_url = "https://www.cellucid.com"
        viewer.web_cache_dir = None
        viewer.allowed_hosts = frozenset()
        viewer._start_server(
            latent_key=None,
            gene_id_column="symbol",
            normalize_embeddings=True,
            centroid_outlier_quantile=0.95,
            centroid_min_points=10,
            dataset_name="Notebook",
            dataset_id="notebook",
            obs_keys=None,
            embedding_keys={"2d": ("X_pca", [0, 1])},
            knn_key="connectivities",
            expression_source="layer:counts",
            vector_field_default=None,
            serve_vector_fields=False,
        )
    kwargs = server_class.call_args.kwargs
    assert kwargs["embedding_keys"] == {"2d": ("X_pca", [0, 1])}
    assert kwargs["knn_key"] == "connectivities"
    assert kwargs["expression_source"] == "layer:counts"
    assert kwargs["gene_id_column"] == "symbol"
    assert "serve_connectivity" not in kwargs


def test_show_anndata_forwards_every_mapping_argument() -> None:
    from cellucid.jupyter import _anndata

    with mock.patch.object(_anndata, "AnnDataViewer") as viewer_class:
        _anndata.show_anndata(
            object(),
            dataset_name="Notebook",
            dataset_id="notebook",
            embedding_keys={"1d": "X_line"},
            knn_key="connectivities",
            expression_source="raw",
        )
    kwargs = viewer_class.call_args.kwargs
    assert kwargs["embedding_keys"] == {"1d": "X_line"}
    assert kwargs["knn_key"] == "connectivities"
    assert kwargs["expression_source"] == "raw"
