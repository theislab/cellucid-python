from __future__ import annotations

import json
import socket
import urllib.request
from http.server import BaseHTTPRequestHandler
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


def _prepared_dataset(
    root: Path,
    *,
    dataset_id: str = "bind-contract",
    dataset_name: str = "Bind contract",
) -> Path:
    root.mkdir()
    (root / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": dataset_id,
                "name": dataset_name,
                "description": "",
                "stats": {
                    "n_cells": 1,
                    "n_genes": 0,
                    "n_obs_fields": 0,
                    "n_categorical_fields": 0,
                    "n_continuous_fields": 0,
                    "has_connectivity": False,
                    "n_edges": None,
                },
                "embeddings": {
                    "available_dimensions": [2],
                    "default_dimension": 2,
                    "files": {"2d": "points_2d.bin"},
                },
                "obs_fields": [],
                "export_settings": {
                    "compression": None,
                    "var_quantization": None,
                    "obs_continuous_quantization": None,
                    "obs_categorical_dtype": "uint16",
                },
            }
        ),
        encoding="utf-8",
    )
    (root / "obs_manifest.json").write_text(
        json.dumps(
            {
                "_format": "compact_v1",
                "n_points": 1,
                "centroid_outlier_quantile": None,
                "latent_key": "latent_space",
                "compression": None,
                "_obsSchemas": {},
                "_continuousFields": [],
                "_categoricalFields": [],
            }
        ),
        encoding="utf-8",
    )
    (root / "points_2d.bin").write_bytes(b"\x00" * 8)
    return root


def _anndata() -> AnnData:
    return AnnData(
        X=np.array([[1.0], [2.0]], dtype=np.float32),
        obs=pd.DataFrame(
            {"label": pd.Categorical(["A", "B"])},
            index=["cell-1", "cell-2"],
        ),
        var=pd.DataFrame(index=["gene-1"]),
        obsm={
            "X_umap_2d": np.array(
                [[0.0, 0.0], [1.0, 1.0]],
                dtype=np.float32,
            )
        },
    )


def test_http_server_uses_the_platforms_nonsharing_port_policy() -> None:
    from cellucid._server_base import CellucidHTTPServer

    exclusive_address_use = getattr(socket, "SO_EXCLUSIVEADDRUSE", None)
    server = CellucidHTTPServer(
        ("127.0.0.1", 0),
        BaseHTTPRequestHandler,
    )
    try:
        if exclusive_address_use is None:
            assert CellucidHTTPServer.allow_reuse_address is True
            assert server.socket.getsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR) != 0
        else:
            assert CellucidHTTPServer.allow_reuse_address is False
            assert (
                server.socket.getsockopt(socket.SOL_SOCKET, exclusive_address_use) == 1
            )
            assert server.socket.getsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR) == 0
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as contender:
            contender.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            with pytest.raises(OSError):
                contender.bind(server.server_address)
    finally:
        server.server_close()


def test_exported_server_publishes_the_exact_os_assigned_port(
    tmp_path: Path,
) -> None:
    from cellucid.server import CellucidServer

    server = CellucidServer(
        _prepared_dataset(tmp_path / "dataset"),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
    )
    server.start_background()
    try:
        assert server._server is not None
        assert server.port == server._server.server_address[1]
        assert server.port > 0
        assert server.url == f"http://127.0.0.1:{server.port}"
        assert server.viewer_url == (f"http://127.0.0.1:{server.port}/?source=remote")
        assert server.server_info["port"] == server.port
    finally:
        server.stop()


def test_exported_server_viewer_url_opens_the_exact_served_catalog(
    tmp_path: Path,
) -> None:
    from cellucid.server import CellucidServer

    root = tmp_path / "catalog"
    root.mkdir()
    _prepared_dataset(root / "alpha", dataset_id="alpha", dataset_name="Alpha")
    _prepared_dataset(root / "beta", dataset_id="beta", dataset_name="Beta")
    server = CellucidServer(
        root,
        host="127.0.0.1",
        port=0,
        open_browser=False,
        quiet=True,
        serve_web_ui=False,
    )
    server.start_background()
    try:
        assert server.viewer_url == f"{server.url}/?source=remote"
        assert "dataset=" not in server.viewer_url
        with urllib.request.urlopen(
            f"{server.url}/_cellucid/datasets",
            timeout=5,
        ) as response:
            payload = json.loads(response.read())
        assert payload == {
            "datasets": [
                {"id": "alpha", "path": "/alpha/", "name": "Alpha"},
                {"id": "beta", "path": "/beta/", "name": "Beta"},
            ]
        }
    finally:
        server.stop()


@pytest.mark.parametrize(
    ("dataset_id", "dataset_name", "message"),
    [
        ("data/set", "Declared name", "dataset_id"),
        ("CON", "Declared name", "dataset_id"),
        ("trailing.", "Declared name", "dataset_id"),
        ("ümlaut", "Declared name", "dataset_id"),
        ("declared-id", " padded", "dataset_name"),
        ("declared-id", "Declared\x7fname", "dataset_name"),
        ("declared-id", "Declared\x85name", "dataset_name"),
    ],
)
def test_exported_server_rejects_nonportable_manifest_identity_before_binding(
    tmp_path: Path,
    dataset_id: str,
    dataset_name: str,
    message: str,
) -> None:
    from cellucid.server import CellucidServer

    with pytest.raises(ValueError, match=message):
        CellucidServer(
            _prepared_dataset(
                tmp_path / "dataset",
                dataset_id=dataset_id,
                dataset_name=dataset_name,
            ),
            host="127.0.0.1",
            port=0,
            quiet=True,
            serve_web_ui=False,
        )


def test_anndata_server_publishes_the_exact_os_assigned_port() -> None:
    from cellucid.anndata_server import AnnDataServer

    server = AnnDataServer(
        _anndata(),
        dataset_name="Bind contract",
        dataset_id="bind-contract",
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
    )
    server.start_background()
    try:
        assert server._server is not None
        assert server.port == server._server.server_address[1]
        assert server.port > 0
        assert server.url == f"http://127.0.0.1:{server.port}"
        assert server.server_info["port"] == server.port
    finally:
        server.stop()


def test_anndata_server_rejects_identity_before_source_classification() -> None:
    from cellucid.anndata_server import AnnDataServer

    with (
        mock.patch(
            "cellucid.anndata_server._classify_anndata_path",
            side_effect=AssertionError("AnnData source was classified"),
        ) as classify,
        pytest.raises(ValueError, match="dataset_id"),
    ):
        AnnDataServer(
            "must-not-be-inspected.h5ad",
            dataset_name="Declared name",
            dataset_id="data/set",
            host="127.0.0.1",
            port=0,
            quiet=True,
            serve_web_ui=False,
        )

    classify.assert_not_called()


def test_explicit_occupied_port_is_rejected_without_port_scanning(
    tmp_path: Path,
) -> None:
    from cellucid.server import CellucidServer

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as occupied:
        occupied.bind(("127.0.0.1", 0))
        occupied.listen(1)
        port = occupied.getsockname()[1]
        server = CellucidServer(
            _prepared_dataset(tmp_path / "dataset"),
            host="127.0.0.1",
            port=port,
            quiet=True,
            serve_web_ui=False,
        )
        with pytest.raises(OSError):
            server.start_background()
        assert server._server is None
        assert server.port == port
        assert server.is_running() is False


def test_anndata_explicit_occupied_port_is_rejected_without_port_scanning() -> None:
    from cellucid.anndata_server import AnnDataServer

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as occupied:
        occupied.bind(("127.0.0.1", 0))
        occupied.listen(1)
        port = occupied.getsockname()[1]
        server = AnnDataServer(
            _anndata(),
            dataset_name="Bind contract",
            dataset_id="bind-contract",
            host="127.0.0.1",
            port=port,
            quiet=True,
            serve_web_ui=False,
        )
        with pytest.raises(OSError):
            server.start_background()
        assert server._server is None
        assert server.port == port
        assert server.is_running() is False


@pytest.mark.parametrize("port", [True, -1, 65536, 8765.0, "8765"])
def test_server_port_rejects_coercion(port: object, tmp_path: Path) -> None:
    from cellucid.server import CellucidServer

    with pytest.raises((TypeError, ValueError), match="port"):
        CellucidServer(
            _prepared_dataset(tmp_path / f"dataset-{type(port).__name__}"),
            port=port,  # type: ignore[arg-type]
            quiet=True,
            serve_web_ui=False,
        )


def test_jupyter_none_port_requests_one_os_assigned_bind(
    monkeypatch,
    tmp_path: Path,
) -> None:
    import cellucid.jupyter as jupyter
    import cellucid.jupyter._exported as jupyter_exported

    seen_ports: list[int] = []

    class DummyServer:
        def __init__(self, **kwargs):
            seen_ports.append(kwargs["port"])
            self.host = kwargs["host"]
            self.port = 43123

        @property
        def url(self) -> str:
            return f"http://{self.host}:{self.port}"

        def start_background(self) -> None:
            return None

        def is_running(self) -> bool:
            return True

        def stop(self) -> None:
            return None

    # ``CellucidViewer._start_server`` resolves ``CellucidServer`` in the module
    # that defines the viewer, so the patch has to name that module.
    monkeypatch.setattr(jupyter_exported, "CellucidServer", DummyServer)
    viewer = jupyter.CellucidViewer(
        data_dir=tmp_path,
        port=None,
        auto_open=False,
    )
    try:
        assert seen_ports == [0]
        assert viewer.port == 43123
    finally:
        viewer.stop()
