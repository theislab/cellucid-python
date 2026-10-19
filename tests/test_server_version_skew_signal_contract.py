"""The signal a newer viewer uses to name an older server.

``pip install cellucid`` fetches the newest web app from www.cellucid.com, so a
viewer that reads export format 3 will meet servers released before it. Such a
server cannot be taught anything new: whatever the viewer learns from it has to
be something every released server already says. Two facts are that signal:

* ``dataset_identity.json`` ``version`` -- 3 from every server that serves
  format 3, 2 (with ``points_<d>d.bin``) from 0.9.2 and earlier;
* ``GET /_cellucid/info`` ``version`` -- the server's package version, a string
  both the AnnData server and the export server have answered since before
  format 3 existed.

A viewer that reads an identity whose ``version`` it does not support fetches
``/_cellucid/info`` from the same server. A JSON object with a string
``version`` means a Cellucid server, and the viewer says "This Cellucid server
({version}) is older than the viewer; update it with pip install -U cellucid."
Anything else -- a 404 or a non-JSON body from a static host -- means a
prepared export written by an older release, to be re-exported. Both facts are
pinned here so neither can drift.
"""

from __future__ import annotations

import http.client
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import cellucid
from cellucid import prepare
from cellucid.anndata_server import AnnDataServer
from cellucid.server import CellucidServer


def _get_json(host: str, port: int, path: str) -> dict:
    connection = http.client.HTTPConnection(host, port, timeout=10)
    try:
        connection.request("GET", path)
        response = connection.getresponse()
        assert response.status == 200, path
        return json.loads(response.read())
    finally:
        connection.close()


def _check(host: str, port: int, identity_path: str) -> None:
    identity = _get_json(host, port, identity_path)
    assert type(identity["version"]) is int and identity["version"] == 3
    info = _get_json(host, port, "/_cellucid/info")
    assert info["version"] == cellucid.__version__
    health = _get_json(host, port, "/_cellucid/health")
    assert health["version"] == cellucid.__version__


def test_the_anndata_server_states_its_identity_version_and_package_version() -> None:
    adata = ad.AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame(index=["a", "b", "c"]),
        var=pd.DataFrame(index=["g"]),
    )
    adata.obsm["X_umap"] = np.array([[0.0, 0.0], [1.0, 2.0], [2.0, 1.0]])
    server = AnnDataServer(
        adata,
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Skew",
        dataset_id="skew",
    )
    server.start_background()
    try:
        assert server._server is not None
        host, port = server._server.server_address
        _check(host, port, "/dataset_identity.json")
    finally:
        server.stop()


def test_the_export_server_states_its_identity_version_and_package_version(
    tmp_path: Path,
) -> None:
    out = tmp_path / "export"
    prepare(
        latent_space=np.array([[0.0], [1.0], [2.0]]),
        obs=pd.DataFrame({"x": [1.0, 2.0, 3.0]}),
        out_dir=out,
        obs_categorical_dtype="uint8",
        dataset_name="Skew",
        dataset_id="skew",
        embeddings={"X_umap": np.array([[0.0, 0.0], [1.0, 2.0], [2.0, 1.0]])},
    )
    server = CellucidServer(out, port=0, quiet=True, serve_web_ui=False)
    server.start_background()
    try:
        assert server._server is not None
        host, port = server._server.server_address
        _check(host, port, "/dataset_identity.json")
    finally:
        server.stop()
