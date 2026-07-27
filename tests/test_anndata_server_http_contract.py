from __future__ import annotations

import gzip
import http.client
import json
from collections.abc import Iterator
from contextlib import ExitStack, contextmanager
from unittest import mock

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.anndata_server import AnnDataServer, _accepts_gzip


def _minimal_adata(*, connectivity: bool = True) -> ad.AnnData:
    adata = ad.AnnData(
        X=np.array([[1.0], [2.0]], dtype=np.float32),
        obs=pd.DataFrame(
            {"score": np.array([0.25, 0.75], dtype=np.float32)},
            index=pd.Index(["cell-1", "cell-2"], dtype=object),
        ),
        var=pd.DataFrame(index=pd.Index(["Gene_A"], dtype=object)),
    )
    adata.obsm["X_umap_2d"] = np.array(
        [[0.0, 1.0], [1.0, 0.0]],
        dtype=np.float32,
    )
    if connectivity:
        adata.obsp["connectivities"] = sparse.csr_matrix(
            np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.float32)
        )
    return adata


@contextmanager
def _running_server(
    adata: ad.AnnData,
) -> Iterator[tuple[AnnDataServer, str, int]]:
    server = AnnDataServer(
        adata,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="HTTP contract",
        dataset_id="http-contract",
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        yield server, host, port
    finally:
        server.stop()


def _request(
    host: str,
    port: int,
    method: str,
    path: str,
    *,
    headers: dict[str, str] | None = None,
) -> tuple[int, dict[str, str], bytes]:
    connection = http.client.HTTPConnection(host, port, timeout=5)
    try:
        connection.request(method, path, headers=headers or {})
        response = connection.getresponse()
        body = response.read()
        return response.status, dict(response.getheaders()), body
    finally:
        connection.close()


def test_file_identity_never_exposes_or_aliases_the_private_source_path(
    tmp_path,
) -> None:
    source_path = tmp_path / "private_source.h5ad"
    _minimal_adata().write_h5ad(source_path)

    server = AnnDataServer(
        source_path,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Private source",
        dataset_id="private-source",
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        status, _headers, body = _request(
            host,
            port,
            "GET",
            "/dataset_identity.json",
        )
        assert status == 200
        identity = json.loads(body)
        encoded = json.dumps(identity)
        assert str(source_path) not in encoded
        assert str(tmp_path) not in encoded
        assert identity["source"] == {
            "name": "H5AD file",
        }
        assert not hasattr(server.adapter, "_source_path")
    finally:
        server.stop()


def test_direct_server_preserves_biological_names_and_uses_exact_wire_components() -> None:
    adata = ad.AnnData(
        X=np.array(
            [
                [1.25, 2.5],
                [3.75, 5.0],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(
            {"cell type": pd.Categorical(["ductal", "endocrine"])},
            index=pd.Index(["cell-1", "cell-2"], dtype=object),
        ),
        var=pd.DataFrame(
            index=pd.Index(["Gt(ROSA)26Sor", "Gene A"], dtype=object),
        ),
    )
    adata.obsm["X_umap_2d"] = np.array(
        [[0.0, 1.0], [1.0, 0.0]],
        dtype=np.float32,
    )

    with _running_server(adata) as (_server, host, port):
        status, _headers, body = _request(
            host,
            port,
            "GET",
            "/var_manifest.json",
        )
        assert status == 200
        assert json.loads(body)["fields"] == [
            ["Gt(ROSA)26Sor"],
            ["Gene A"],
        ]

        status, _headers, body = _request(
            host,
            port,
            "GET",
            "/var/Gt_ROSA_26Sor.values.f32",
        )
        assert status == 200
        assert np.frombuffer(body, dtype="<f4").tolist() == [1.25, 3.75]

        status, _headers, body = _request(
            host,
            port,
            "GET",
            "/obs/cell_type.codes.u8",
        )
        assert status == 200
        assert np.frombuffer(body, dtype=np.uint8).tolist() == [0, 1]


@pytest.mark.parametrize(
    ("matrix", "cause"),
    [
        (
            np.array([[0.0, np.nan], [np.nan, 0.0]], dtype=np.float32),
            "finite",
        ),
        (
            np.array([[0.0, 1.0], [0.0, 0.0]], dtype=np.float32),
            "symmetric",
        ),
        (
            np.array([[1.0, 0.0], [0.0, 0.0]], dtype=np.float32),
            "diagonal",
        ),
        (
            np.array([[0.0, 0.5], [0.75, 0.0]], dtype=np.float32),
            "symmetric",
        ),
        (
            np.array([[0.0, -1.0], [-1.0, 0.0]], dtype=np.float32),
            "non-negative",
        ),
    ],
)
def test_corrupt_connectivity_aborts_adapter_initialization(
    matrix: np.ndarray,
    cause: str,
) -> None:
    adata = _minimal_adata(connectivity=False)
    adata.obsp["connectivities"] = sparse.csr_matrix(matrix)

    with pytest.raises(
        ValueError,
        match="direct server was not initialized",
    ) as error:
        AnnDataAdapter(
            adata,
            dataset_name="Corrupt connectivity",
            dataset_id="corrupt-connectivity",
        )

    assert error.value.__cause__ is not None
    assert cause in str(error.value.__cause__)


def test_identity_and_manifest_are_one_validated_connectivity_capability() -> None:
    connected = AnnDataAdapter(
        _minimal_adata(),
        dataset_name="Connected",
        dataset_id="connected",
    )
    identity = connected.get_dataset_identity()
    manifest = connected.get_connectivity_manifest()
    assert manifest is not None
    assert identity["stats"]["has_connectivity"] is True
    assert identity["stats"]["n_edges"] == manifest["n_edges"] == 1
    assert manifest["max_neighbors"] == 1
    assert manifest["weightsPath"] == "connectivity/edges.weights.f64.bin"
    assert manifest["weight_dtype"] == "float64"
    assert manifest["weight_bytes"] == 8
    assert np.frombuffer(
        connected.get_connectivity_edges()[2],
        dtype="<f8",
    ).tolist() == [1.0]

    disconnected = AnnDataAdapter(
        _minimal_adata(connectivity=False),
        dataset_name="Disconnected",
        dataset_id="disconnected",
    )
    identity = disconnected.get_dataset_identity()
    assert identity["stats"]["has_connectivity"] is False
    assert identity["stats"]["n_edges"] is None
    assert disconnected.get_connectivity_manifest() is None


def test_present_empty_connectivity_remains_a_three_payload_capability() -> None:
    adata = _minimal_adata(connectivity=False)
    adata.obsp["connectivities"] = sparse.csr_matrix(
        (adata.n_obs, adata.n_obs),
        dtype=np.float64,
    )

    with _running_server(adata) as (server, host, port):
        identity = server.adapter.get_dataset_identity()
        manifest = server.adapter.get_connectivity_manifest()
        assert identity["stats"]["has_connectivity"] is True
        assert identity["stats"]["n_edges"] == 0
        assert manifest is not None
        assert manifest["n_edges"] == 0

        for path in [
            "/connectivity/edges.src.bin",
            "/connectivity/edges.dst.bin",
            "/connectivity/edges.weights.f64.bin",
        ]:
            status, headers, body = _request(host, port, "GET", path)
            assert status == 200
            assert int(headers["Content-Length"]) == 0
            assert body == b""


def test_server_initialization_rejects_an_advertised_missing_manifest() -> None:
    with (
        mock.patch.object(
            AnnDataAdapter,
            "get_connectivity_manifest",
            return_value=None,
        ),
        pytest.raises(
            ValueError,
            match="advertises connectivity without a connectivity manifest",
        ),
    ):
        AnnDataServer(
            _minimal_adata(),
            host="127.0.0.1",
            port=0,
            quiet=True,
            serve_web_ui=False,
            dataset_name="Missing manifest",
            dataset_id="missing-manifest",
        )


def test_head_uses_only_prevalidated_metadata_and_never_materializes_payloads() -> None:
    server = AnnDataServer(
        _minimal_adata(),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="HEAD contract",
        dataset_id="head-contract",
    )
    blocked_methods = [
        "get_dataset_identity",
        "get_obs_manifest",
        "get_var_manifest",
        "get_connectivity_manifest",
        "get_points_binary",
        "get_obs_continuous_values",
        "get_gene_expression",
        "get_connectivity_edges",
    ]

    with ExitStack() as stack:
        for method_name in blocked_methods:
            stack.enter_context(
                mock.patch.object(
                    server.adapter,
                    method_name,
                    side_effect=AssertionError(f"HEAD materialized through {method_name}"),
                )
            )

        server.start_background()
        assert server._server is not None
        host, port = server._server.server_address
        try:
            expected_lengths = {
                **{path: len(body) for path, body in server._metadata_bodies.items()},
                **server._payload_lengths,
            }
            for path, expected_length in expected_lengths.items():
                status, headers, body = _request(
                    host,
                    port,
                    "HEAD",
                    f"/{path}",
                    headers={"Accept-Encoding": "gzip"},
                )
                assert status == 200, path
                assert body == b"", path
                assert int(headers["Content-Length"]) == expected_length, path
                assert "Content-Encoding" not in headers, path
        finally:
            server.stop()


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        ("gzip", True),
        ("br, gzip;q=0.5", True),
        ("GZip; q=1.000", True),
        ("gzip;q=0", False),
        ("gzip;q=0, *;q=1", False),
        ("x-gzip", False),
        ("br, *;q=0.2", True),
        ("gzip;q=1.1", False),
        ("gzip;q=1; q=0", False),
    ],
)
def test_accept_encoding_parser_obeys_exact_tokens_and_quality(
    header: str,
    expected: bool,
) -> None:
    assert _accepts_gzip(header) is expected


def test_get_honors_gzip_quality_zero_and_positive_quality() -> None:
    with _running_server(_minimal_adata()) as (server, host, port):
        expected = server.adapter.get_points_binary(2)

        status, headers, body = _request(
            host,
            port,
            "GET",
            "/points_2d.bin",
            headers={"Accept-Encoding": "br, gzip;q=0"},
        )
        assert status == 200
        assert "Content-Encoding" not in headers
        assert body == expected

        status, headers, body = _request(
            host,
            port,
            "GET",
            "/points_2d.bin",
            headers={"Accept-Encoding": "br, gzip;q=0.5"},
        )
        assert status == 200
        assert headers["Content-Encoding"] == "gzip"
        assert gzip.decompress(body) == expected


def test_connectivity_weights_have_exact_head_get_and_gzip_semantics() -> None:
    adata = _minimal_adata(connectivity=False)
    adata.obsp["connectivities"] = sparse.csr_matrix(
        np.array([[0.0, 2.75], [2.75, 0.0]], dtype=np.float64)
    )
    with _running_server(adata) as (server, host, port):
        (
            _sources,
            _destinations,
            expected_weights,
            n_edges,
            _max_neighbors,
        ) = server.adapter.get_connectivity_edges()
        assert n_edges == 1
        assert np.frombuffer(expected_weights, dtype="<f8").tolist() == [2.75]

        status, headers, body = _request(
            host,
            port,
            "HEAD",
            "/connectivity/edges.weights.f64.bin",
            headers={"Accept-Encoding": "gzip"},
        )
        assert status == 200
        assert body == b""
        assert int(headers["Content-Length"]) == 8
        assert "Content-Encoding" not in headers

        status, headers, body = _request(
            host,
            port,
            "GET",
            "/connectivity/edges.weights.f64.bin",
            headers={"Accept-Encoding": "gzip;q=0"},
        )
        assert status == 200
        assert "Content-Encoding" not in headers
        assert int(headers["Content-Length"]) == len(expected_weights)
        assert body == expected_weights

        status, headers, body = _request(
            host,
            port,
            "GET",
            "/connectivity/edges.weights.f64.bin",
            headers={"Accept-Encoding": "gzip"},
        )
        assert status == 200
        assert headers["Content-Encoding"] == "gzip"
        assert int(headers["Content-Length"]) == len(body)
        assert gzip.decompress(body) == expected_weights


def test_file_backed_weighted_h5ad_serves_all_aligned_graph_payloads(
    tmp_path,
) -> None:
    adata = _minimal_adata(connectivity=False)
    adata.obsp["connectivities"] = sparse.csr_matrix(
        np.array([[0.0, 2.75], [2.75, 0.0]], dtype=np.float64)
    )
    source_path = tmp_path / "weighted.h5ad"
    adata.write_h5ad(source_path)

    server = AnnDataServer(
        source_path,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Weighted file-backed H5AD",
        dataset_id="weighted-file-backed-h5ad",
    )
    assert server.adapter.is_backed is True
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        status, _headers, manifest_body = _request(
            host,
            port,
            "GET",
            "/connectivity_manifest.json",
        )
        assert status == 200
        manifest = json.loads(manifest_body)
        assert manifest["sourcesPath"] == "connectivity/edges.src.bin"
        assert manifest["destinationsPath"] == "connectivity/edges.dst.bin"
        assert manifest["weightsPath"] == "connectivity/edges.weights.f64.bin"
        assert manifest["n_edges"] == 1

        expected_payloads = {
            "/connectivity/edges.src.bin": np.array([0], dtype="<u2").tobytes(),
            "/connectivity/edges.dst.bin": np.array([1], dtype="<u2").tobytes(),
            "/connectivity/edges.weights.f64.bin": np.array([2.75], dtype="<f8").tobytes(),
        }
        for path, expected in expected_payloads.items():
            status, headers, body = _request(
                host,
                port,
                "GET",
                path,
                headers={"Accept-Encoding": "gzip"},
            )
            assert status == 200
            assert headers["Content-Encoding"] == "gzip"
            assert int(headers["Content-Length"]) == len(body)
            assert gzip.decompress(body) == expected
    finally:
        server.stop()


def test_real_weighted_zarr_serves_the_same_float64_payload(
    tmp_path,
) -> None:
    pytest.importorskip("zarr")
    adata = _minimal_adata(connectivity=False)
    adata.obsp["connectivities"] = sparse.csr_matrix(
        np.array([[0.0, 6.125], [6.125, 0.0]], dtype=np.float64)
    )
    source_path = tmp_path / "weighted.zarr"
    with pytest.warns(UserWarning, match="Writing zarr v2 data"):
        adata.write_zarr(source_path)

    server = AnnDataServer(
        source_path,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Weighted Zarr",
        dataset_id="weighted-zarr",
    )
    assert server.adapter.is_backed is False
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        status, headers, body = _request(
            host,
            port,
            "GET",
            "/connectivity/edges.weights.f64.bin",
        )
        assert status == 200
        assert int(headers["Content-Length"]) == 8
        assert np.frombuffer(body, dtype="<f8").tolist() == [6.125]
    finally:
        server.stop()


@pytest.mark.parametrize(
    "path",
    [
        "/connectivity/edges.weights.bin",
        "/connectivity/edges.weights.f32.bin",
        "/connectivity/edges.weights.f64.bin.gz",
        "/connectivity/edges.weights.f64.bin/extra",
    ],
)
def test_connectivity_rejects_every_non_current_weight_route(path: str) -> None:
    with _running_server(_minimal_adata()) as (_server, host, port):
        status, _headers, _body = _request(host, port, "GET", path)
        assert status == 404


def test_server_publishes_the_exact_vector_manifest_and_payload() -> None:
    adata = _minimal_adata()
    vectors = np.array(
        [[0.25, -0.5], [0.75, 1.0]],
        dtype=np.float32,
    )
    adata.obsm["velocity_umap_2d"] = vectors

    with _running_server(adata) as (server, host, port):
        status, _headers, body = _request(
            host,
            port,
            "GET",
            "/dataset_identity.json",
        )
        assert status == 200
        identity = json.loads(body)
        assert identity["vector_fields"] == {
            "default_field": "velocity_umap",
            "fields": {
                "velocity_umap": {
                    "label": "velocity_umap",
                    "basis": "umap",
                    "available_dimensions": [2],
                    "default_dimension": 2,
                    "files": {
                        "2d": "vectors/velocity_umap_2d.bin",
                    },
                }
            },
        }

        expected = server.adapter.get_vector_field_binary(
            "velocity_umap",
            2,
        )
        status, headers, body = _request(
            host,
            port,
            "GET",
            "/vectors/velocity_umap_2d.bin",
            headers={"Accept-Encoding": "gzip"},
        )
        assert status == 200
        assert headers["Content-Encoding"] == "gzip"
        assert gzip.decompress(body) == expected


@pytest.mark.parametrize(
    "path",
    [
        "/points_2d.bin.extra.bin",
        "/points_2d.bin.gz",
        "/points_2d.bin/extra",
        "/points_22d.bin",
    ],
)
def test_points_route_rejects_every_non_exact_suffix(path: str) -> None:
    server = AnnDataServer(
        _minimal_adata(),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Route contract",
        dataset_id="route-contract",
    )
    with mock.patch.object(
        server.adapter,
        "get_points_binary",
        side_effect=AssertionError("invalid route reached point payload"),
    ):
        server.start_background()
        assert server._server is not None
        host, port = server._server.server_address
        try:
            status, _headers, _body = _request(host, port, "GET", path)
            assert status == 404
        finally:
            server.stop()
