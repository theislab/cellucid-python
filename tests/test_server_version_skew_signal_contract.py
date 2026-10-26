"""The signal a newer viewer uses to name an older server.

``pip install cellucid`` fetches the newest web app from www.cellucid.com, so a
viewer meets servers released before it. Such a server cannot be taught
anything new: whatever the viewer learns from it has to be something every
released server already says. ``GET /_cellucid/info`` and
``GET /_cellucid/health`` answer the server's package version as a string, and
that is the signal. It is pinned here so it cannot drift.
"""

from __future__ import annotations

import http.client
import json

import anndata as ad
import numpy as np
import pandas as pd

import cellucid
from cellucid.anndata_server import AnnDataServer


def _get_json(host: str, port: int, path: str) -> dict:
    connection = http.client.HTTPConnection(host, port, timeout=10)
    try:
        connection.request("GET", path)
        response = connection.getresponse()
        assert response.status == 200, path
        return json.loads(response.read())
    finally:
        connection.close()


def test_the_anndata_server_states_its_package_version() -> None:
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
        # The server publishes no identity document: the viewer reads the
        # store itself. Its package version is the signal.
        assert _get_json(host, port, "/_cellucid/info")["version"] == cellucid.__version__
        assert _get_json(host, port, "/_cellucid/health")["version"] == cellucid.__version__
    finally:
        server.stop()
