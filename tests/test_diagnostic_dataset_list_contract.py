"""The diagnostic dataset list must accept exactly what the server serves.

``debug_connection()`` fetches ``/_cellucid/datasets`` and validates it with
:func:`cellucid.jupyter._wire._require_debug_dataset_list` before probing every
listed store. The validator and the server hold two separate declarations of
the same entry shape, so these tests tie them together in both directions:
what the server emits is accepted, and anything else is refused.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from _serverless_viewer import serverless_viewer

from cellucid import prepare
from cellucid.anndata_server import AnnDataServer
from cellucid.jupyter._wire import _require_debug_dataset_list


def _served_entry(tmp_path: Path) -> dict[str, Any]:
    path = prepare(
        tmp_path / "alpha.h5ad",
        format="h5ad",
        embeddings={"X_umap": np.zeros((3, 2))},
        dataset_id="alpha",
        dataset_name="Alpha",
        gene_expression=np.ones((3, 1)),
        var=pd.DataFrame(index=["gene-1"]),
    )
    server = AnnDataServer(path, quiet=True, serve_web_ui=False)
    document = server._datasets_document()
    assert len(document["datasets"]) == 1
    return document


def test_diagnostics_accept_the_exact_catalog_the_server_serves(tmp_path: Path) -> None:
    """The served entry is the input; no hand-written stand-in can drift from it."""
    document = _served_entry(tmp_path)
    assert _require_debug_dataset_list(document) == document["datasets"]


def test_debug_connection_probes_every_listed_store(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = serverless_viewer(
        port=8765,
        client_server_url="http://127.0.0.1:8765",
        web_cache_dir=tmp_path / "cache",
    )
    document = _served_entry(tmp_path)
    responses: dict[str, object] = {
        "http://127.0.0.1:8765/_cellucid/health": {"status": "ok"},
        "http://127.0.0.1:8765/_cellucid/info": {"version": "0.9.5"},
        "http://127.0.0.1:8765/_cellucid/datasets": document,
    }
    probed: list[str] = []

    class Response:
        def __init__(self, payload: object) -> None:
            self._payload = json.dumps(payload).encode("utf-8")
            self.headers: dict[str, str] = {}

        def __enter__(self) -> Response:
            return self

        def __exit__(self, *_: object) -> None:
            return None

        def read(self) -> bytes:
            return self._payload

    def urlopen(request: Any, *, timeout: float) -> Response:
        url = request.full_url if hasattr(request, "full_url") else request
        if url not in responses:
            raise OSError(f"unmapped diagnostic URL: {url}")
        return Response(responses[url])

    def probe(url: str) -> dict[str, Any]:
        probed.append(url)
        return {"url": url}

    monkeypatch.setattr("urllib.request.urlopen", urlopen)
    monkeypatch.setattr("cellucid.jupyter._viewer._probe_served_store", probe)
    report = viewer.debug_connection()

    assert "server_datasets_error" not in report
    assert report["server_datasets"] == document["datasets"]
    store_url = f"http://127.0.0.1:8765/{document['datasets'][0]['anndata']['file']['url']}"
    assert probed == [store_url]
    assert report["anndata_store_probes"] == {"alpha": {"url": store_url}}
    assert "dataset_identity_probes" not in report


@pytest.mark.parametrize(
    ("entry", "match"),
    [
        pytest.param(
            {"id": "alpha", "path": "/alpha/", "name": "Alpha"},
            "exactly id, name, and anndata",
            id="a-directory-entry",
        ),
        pytest.param(
            {"id": "alpha", "name": "Alpha", "anndata": {"file": {"url": "x.h5ad"}, "defaults": {}}},
            "anndata must name one served store",
            id="a-store-outside-the-route",
        ),
        pytest.param(
            {"id": "alpha", "name": "Alpha", "anndata": {"file": {"url": "_cellucid/anndata/a.h5ad"}}},
            "anndata must name one served store",
            id="no-defaults",
        ),
        pytest.param(
            {
                "id": "alpha",
                "name": "Alpha",
                "anndata": {"file": {"url": "_cellucid/anndata/a.h5ad"}, "defaults": {}},
                "extra": 1,
            },
            "exactly id, name, and anndata",
            id="unknown-extra-key",
        ),
    ],
)
def test_diagnostics_reject_every_other_entry(entry: dict[str, object], match: str) -> None:
    with pytest.raises((TypeError, ValueError), match=match):
        _require_debug_dataset_list({"datasets": [entry]})
