from __future__ import annotations

import json
import socket
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


def _prepared_dataset(root: Path) -> Path:
    root.mkdir()
    (root / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": "bind-contract",
                "name": "Bind contract",
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
        assert server.server_info["port"] == server.port
    finally:
        server.stop()


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

    monkeypatch.setattr(jupyter, "CellucidServer", DummyServer)
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
