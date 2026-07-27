from __future__ import annotations

import inspect
import json
from types import SimpleNamespace
from typing import Any

import pytest

import cellucid.jupyter as jupyter
from cellucid.anndata_server import AnnDataRequestHandler
from cellucid.server import CORSRequestHandler


def _undisplayed_viewer() -> jupyter.BaseViewer:
    viewer = object.__new__(jupyter.BaseViewer)
    viewer._displayed = False
    viewer._viewer_id = "viewer-1"
    viewer._viewer_token = "secret-1"
    return viewer


@pytest.mark.parametrize(
    "message",
    [
        None,
        [],
        {},
        {"type": ""},
        {"type": "bad type"},
        {"type": "ping", "viewerId": "forged"},
        {"type": "ping", "viewerToken": "forged"},
        {"type": "custom", "value": float("nan")},
        {"type": "custom", "value": 7},
    ],
)
def test_send_message_rejects_noncanonical_commands_before_display(
    message: object,
) -> None:
    viewer = _undisplayed_viewer()
    with pytest.raises((TypeError, ValueError), match="message|type|routing|JSON"):
        viewer.send_message(message)  # type: ignore[arg-type]


def test_python_command_helpers_reject_coercion_before_publication(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = _undisplayed_viewer()
    published: list[dict[str, Any]] = []
    monkeypatch.setattr(viewer, "send_message", published.append)

    with pytest.raises(TypeError, match="cell_indices"):
        viewer.highlight_cells([0, True])
    with pytest.raises(ValueError, match="non-empty"):
        viewer.highlight_cells([])
    with pytest.raises(ValueError, match="duplicate"):
        viewer.highlight_cells([0, 0])
    with pytest.raises(ValueError, match="color"):
        viewer.highlight_cells([0], color="red")
    with pytest.raises(TypeError, match="field"):
        viewer.set_color_by(7)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="visible"):
        viewer.set_visibility([0], visible=1)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="cell_indices"):
        viewer.set_visibility("0")  # type: ignore[arg-type]
    assert published == []

    viewer.highlight_cells([0, 2], color="#00ff00")
    viewer.set_color_by("cell_type")
    viewer.set_visibility(None, visible=False)
    viewer.reset_view()
    assert published == [
        {
            "type": "highlight",
            "cells": [0, 2],
            "color": "#00ff00",
        },
        {"type": "setColorBy", "field": "cell_type"},
        {"type": "setVisibility", "cells": None, "visible": False},
        {"type": "resetCamera"},
    ]


def test_generated_parent_bridge_throws_instead_of_fabricating_delivery() -> None:
    viewer = _undisplayed_viewer()
    viewer._displayed = True
    viewer.height = 600
    viewer._client_server_url = "https://viewer.example"

    html = viewer._generate_viewer_html()
    assert "Cellucid viewer iframe is unavailable" in html
    assert "console.warn('Cellucid viewer not found" not in html


def test_prepared_health_rejects_missing_server_version() -> None:
    handler = object.__new__(CORSRequestHandler)
    handler.server_info = {}
    handler.datasets = []
    handler.serve_web_ui = False
    handler.path = "/_cellucid/health"
    handler.send_json = lambda *_args, **_kwargs: pytest.fail("health must not fabricate a version")

    with pytest.raises(KeyError):
        handler._handle_read_request(head_only=False)


def test_anndata_health_rejects_missing_server_fields() -> None:
    handler = object.__new__(AnnDataRequestHandler)
    handler.server_info = {}
    handler.adapter = SimpleNamespace(n_cells=2, n_genes=1)
    handler.metadata_bodies = {}
    handler.payload_lengths = {}
    handler.serve_web_ui = False
    handler.path = "/_cellucid/health"
    handler.headers = {}
    responses: list[tuple[int, str]] = []
    handler.send_json = lambda *_args, **_kwargs: pytest.fail("health must not fabricate fields")
    handler.send_error_response = lambda status, message: responses.append((status, message))

    handler.do_GET()
    assert responses == [(500, "Internal server error")]


def test_valid_message_serialization_is_strict_json(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = _undisplayed_viewer()
    viewer._displayed = True
    captured: list[Any] = []
    monkeypatch.setattr("IPython.display.display", captured.append)

    viewer.send_message({"type": "ping", "requestId": "request-1"})

    assert len(captured) == 1
    source = captured[0].data
    assert (
        json.dumps(
            {"type": "ping", "requestId": "request-1"},
            separators=(",", ":"),
        )
        in source
    )
    assert "Cellucid viewer command target is unavailable" in source


def test_debug_connection_exposes_only_current_diagnostics() -> None:
    parameters = inspect.signature(jupyter.BaseViewer.debug_connection).parameters
    assert list(parameters) == ["self", "timeout"]


def test_debug_connection_probes_every_declared_dataset_by_exact_id_and_path(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = jupyter.BaseViewer(
        port=8765,
        client_server_url="http://127.0.0.1:8765",
        web_cache_dir=tmp_path / "web-cache",
    )
    requested: list[str] = []

    class Response:
        def __init__(self, payload: object) -> None:
            self._body = payload if type(payload) is bytes else json.dumps(payload).encode("utf-8")
            self.headers = {"Content-Type": "application/json"}

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self) -> bytes:
            return self._body

    responses: dict[str, object] = {
        "http://127.0.0.1:8765/_cellucid/health": {"status": "ok"},
        "http://127.0.0.1:8765/_cellucid/info": {"version": "0.9.1"},
        "http://127.0.0.1:8765/_cellucid/datasets": {
            "datasets": [
                {"id": "alpha", "path": "/alpha/", "name": "Alpha"},
                {"id": "beta", "path": "/beta/", "name": "Beta"},
            ]
        },
        "http://127.0.0.1:8765/alpha/dataset_identity.json": {
            "id": "alpha",
            "name": "Alpha",
        },
        "http://127.0.0.1:8765/beta/dataset_identity.json": {
            "id": "wrong-id",
            "name": "Beta",
        },
    }

    def urlopen(request, *, timeout: float):
        url = request.full_url if hasattr(request, "full_url") else request
        requested.append(url)
        if url not in responses:
            raise OSError(f"unmapped diagnostic URL: {url}")
        assert timeout == 2
        return Response(responses[url])

    monkeypatch.setattr("urllib.request.urlopen", urlopen)
    report = viewer.debug_connection()

    assert report["server_datasets"] == [
        {"id": "alpha", "path": "/alpha/", "name": "Alpha"},
        {"id": "beta", "path": "/beta/", "name": "Beta"},
    ]
    assert report["dataset_identity_probes"]["alpha"] == {
        "path": "/alpha/",
        "url": "http://127.0.0.1:8765/alpha/dataset_identity.json",
        "identity": {"id": "alpha", "name": "Alpha"},
    }
    assert report["dataset_identity_probes"]["beta"] == {
        "path": "/beta/",
        "url": "http://127.0.0.1:8765/beta/dataset_identity.json",
        "error": "dataset_identity.json id does not match declared dataset 'beta'",
    }
    assert "dataset_identity" not in report
    assert "dataset_identity_url" not in report
    assert "dataset_identity_error" not in report
    assert requested.count("http://127.0.0.1:8765/alpha/dataset_identity.json") == 1
    assert requested.count("http://127.0.0.1:8765/beta/dataset_identity.json") == 1
