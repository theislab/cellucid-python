from __future__ import annotations

import json
import socket
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellucid.anndata_server import AnnDataServer
from cellucid.server import CellucidServer


def _prepared_dataset(root: Path) -> Path:
    root.mkdir()
    (root / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": "lifecycle-contract",
                "name": "Lifecycle contract",
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
        X=np.array([[1.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["gene-1"]),
        obsm={
            "X_umap_2d": np.array(
                [[0.0, 0.0]],
                dtype=np.float32,
            )
        },
    )


def _server(tmp_path: Path, **kwargs: object) -> CellucidServer:
    return CellucidServer(
        _prepared_dataset(tmp_path / "dataset"),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        **kwargs,
    )


def _assert_socket_released(port: int) -> None:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", port))


@pytest.mark.parametrize("argument", ["open_browser", "quiet"])
@pytest.mark.parametrize("value", [0, 1, None, "false", np.bool_(False)])
def test_server_boolean_arguments_reject_truthiness(
    argument: str,
    value: object,
    tmp_path: Path,
) -> None:
    prepared_kwargs = {
        "data_dir": tmp_path / "missing",
        "port": 0,
        "open_browser": False,
        "quiet": True,
        "serve_web_ui": False,
    }
    prepared_kwargs[argument] = value
    with pytest.raises(TypeError, match=argument):
        CellucidServer(**prepared_kwargs)  # type: ignore[arg-type]

    anndata_kwargs = {
        "data": _anndata(),
        "port": 0,
        "open_browser": False,
        "quiet": True,
        "serve_web_ui": False,
        "dataset_name": "Boolean contract",
        "dataset_id": "boolean-contract",
    }
    anndata_kwargs[argument] = value
    with pytest.raises(TypeError, match=argument):
        AnnDataServer(**anndata_kwargs)  # type: ignore[arg-type]


def test_prepared_server_stop_is_complete_and_single_use(
    tmp_path: Path,
) -> None:
    server = _server(tmp_path)
    server.start_background()
    port = server.port
    server.stop()

    assert server.is_running() is False
    assert server._server is None
    assert server._thread is None
    _assert_socket_released(port)
    with pytest.raises(RuntimeError, match="not running"):
        _ = server.url
    with pytest.raises(RuntimeError, match="not running"):
        _ = server.viewer_url
    with pytest.raises(RuntimeError, match="single-use"):
        server.start_background()


def test_prepared_server_stop_before_start_consumes_instance(
    tmp_path: Path,
) -> None:
    server = _server(tmp_path)
    server.stop()

    with pytest.raises(RuntimeError, match="single-use"):
        server.start_background()
    with pytest.raises(RuntimeError, match="not running"):
        _ = server.url


def test_prepared_server_thread_start_failure_releases_bound_socket(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    server = _server(tmp_path)
    original = RuntimeError("exact prepared thread start failure")

    def fail_start(_thread: object) -> None:
        raise original

    monkeypatch.setattr("cellucid.server.threading.Thread.start", fail_start)
    with pytest.raises(RuntimeError, match="exact prepared thread") as raised:
        server.start_background()

    assert raised.value is original
    assert server._server is None
    assert server.is_running() is False
    _assert_socket_released(server.port)
    with pytest.raises(RuntimeError, match="single-use"):
        server.start_background()


def test_prepared_server_browser_failure_rolls_back_live_thread(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    server = _server(tmp_path, open_browser=True)
    original = RuntimeError("exact prepared browser failure")

    def fail_browser(_url: str) -> bool:
        raise original

    monkeypatch.setattr("cellucid.server.webbrowser.open", fail_browser)
    with pytest.raises(RuntimeError, match="exact prepared browser") as raised:
        server.start_background()

    assert raised.value is original
    assert server._server is None
    assert server._thread is None
    assert server.is_running() is False
    _assert_socket_released(server.port)


def test_prepared_server_false_browser_result_is_a_start_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    server = _server(tmp_path, open_browser=True)
    monkeypatch.setattr(
        "cellucid.server.webbrowser.open",
        lambda _url: False,
    )

    with pytest.raises(RuntimeError, match="Could not open the browser"):
        server.start_background()

    assert server._server is None
    assert server._thread is None
    assert server.is_running() is False
    _assert_socket_released(server.port)


def test_prepared_background_failure_is_retained_for_wait(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    server = _server(tmp_path)
    original = RuntimeError("exact prepared serve failure")

    def fail_serve(_server: object) -> None:
        raise original

    monkeypatch.setattr(
        "cellucid.server.HTTPServer.serve_forever",
        fail_serve,
    )
    server.start_background()
    with pytest.raises(RuntimeError, match="exact prepared serve") as raised:
        server.wait()

    assert raised.value is original
    assert server._server is None
    assert server.is_running() is False
    _assert_socket_released(server.port)
