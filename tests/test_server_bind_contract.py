from __future__ import annotations

import socket
from http.server import BaseHTTPRequestHandler
from unittest import mock

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


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
            "cellucid.anndata_server.classify_anndata_path",
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
def test_server_port_rejects_coercion(port: object) -> None:
    from cellucid.anndata_server import AnnDataServer

    with pytest.raises((TypeError, ValueError), match="port"):
        AnnDataServer(
            _anndata(),
            dataset_name="Bind contract",
            dataset_id="bind-contract",
            port=port,  # type: ignore[arg-type]
            quiet=True,
            serve_web_ui=False,
        )


def test_jupyter_none_port_requests_8765_and_the_next_free_port(monkeypatch) -> None:
    """8765 is the one local port Google sign-in accepts; a taken one is passed over (as in R)."""
    import cellucid.anndata_server as anndata_server
    import cellucid.jupyter as jupyter

    seen_ports: list[int] = []
    servers: list[object] = []

    class DummyServer:
        def __init__(self, **kwargs):
            seen_ports.append(kwargs["port"])
            servers.append(self)
            self.host = kwargs["host"]
            self.port = 43123
            self.dataset_id = "demo"
            self.server_info = {"n_cells": 3}

        @property
        def url(self) -> str:
            return f"http://{self.host}:{self.port}"

        def start_background(self) -> None:
            return None

        def is_running(self) -> bool:
            return True

        def stop(self) -> None:
            return None

    # ``AnnDataViewer._start_server`` imports ``AnnDataServer`` from its module
    # when it runs, so the patch names that module.
    monkeypatch.setattr(anndata_server, "AnnDataServer", DummyServer)
    viewer = jupyter.AnnDataViewer(
        "data.h5ad",
        port=None,
        auto_open=False,
    )
    try:
        assert seen_ports == [8765]
        assert servers[0].next_free_port is True
        assert viewer.port == 43123
    finally:
        viewer.stop()
