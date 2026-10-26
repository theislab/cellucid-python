from __future__ import annotations

import inspect
import json
from typing import Any

import numpy as np
import pytest
from _serverless_viewer import serverless_viewer

import cellucid._viewer_channel as channel
import cellucid.jupyter as jupyter
import cellucid.jupyter._viewer as jupyter_viewer
import cellucid.jupyter._wire as wire


def _undisplayed_viewer() -> jupyter.AnnDataViewer:
    viewer = object.__new__(jupyter.AnnDataViewer)
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

    with pytest.raises(ValueError, match="non-empty"):
        viewer.set_visibility([], visible=False)
    with pytest.raises(TypeError, match="cell_indices"):
        viewer.highlight_cells(np.array([0.0, 1.0]))
    with pytest.raises(TypeError, match="cell_indices"):
        viewer.highlight_cells(np.array([True, False]))
    with pytest.raises(TypeError, match="non-negative"):
        viewer.highlight_cells(np.array([3, -1]))
    with pytest.raises(ValueError, match="below 2\\*\\*32"):
        viewer.highlight_cells(np.array([2**32]))
    assert published == []

    viewer.highlight_cells([2, 0], color="#00ff00")
    viewer.highlight_cells(np.array([7, 5], dtype=np.int64), color="#0000ff")
    viewer.set_color_by("cell_type")
    viewer.set_visibility(None, visible=False)
    viewer.set_visibility(np.array([4, 1]), visible=True)
    viewer.reset_view()
    # The cells are checked once and handed on ascending, as uint32.
    cells = [entry.pop("cells") for entry in published if "cells" in entry]
    assert [None if value is None else (value.dtype.str, value.tolist()) for value in cells] == [
        ("<u4", [0, 2]),
        ("<u4", [5, 7]),
        None,
        ("<u4", [1, 4]),
    ]
    assert published == [
        {"type": "highlight", "color": "#00ff00"},
        {"type": "highlight", "color": "#0000ff"},
        {"type": "setColorBy", "field": "cell_type"},
        {"type": "setVisibility", "visible": False},
        {"type": "setVisibility", "visible": True},
        {"type": "resetCamera"},
    ]


def test_the_notebook_output_is_one_iframe_and_no_script() -> None:
    # The viewer opens its own channel; the output carries no code, so it is
    # the same in every front end, trusted or not.
    viewer = _undisplayed_viewer()
    viewer.height = 600
    viewer._client_server_url = "https://viewer.example"

    html = viewer._generate_viewer_html()
    assert "<script" not in html
    assert "postMessage" not in html
    assert html.count("<iframe") == 1
    assert 'src="https://viewer.example/?viewer=true&amp;viewerId=viewer-1&amp;' in html


def test_the_viewer_page_address_uses_the_shared_protocols_keys() -> None:
    """The page address and the host's name are the ones all three repositories pin."""
    import json
    from pathlib import Path
    from urllib.parse import parse_qsl, urlsplit

    protocol = json.loads(
        (Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text(encoding="utf-8")
    )
    viewer = _undisplayed_viewer()
    viewer._client_server_url = "https://viewer.example"
    query = parse_qsl(urlsplit(viewer.viewer_url).query)
    assert [key for key, _ in query] == protocol["page"]["parameters"]
    assert query[0] == (protocol["page"]["mode"], "true")
    assert set(protocol["runtimes"]) == {"python", "jupyter", "r"}


def test_valid_message_serialization_is_strict_json() -> None:
    viewer = _undisplayed_viewer()
    viewer._message_routing_registered = True
    channel.register_event_callback("viewer-1", "secret-1", lambda _event: None)
    try:
        assert viewer.send_message({"type": "ping", "requestId": "request-1"}) is False
        registration = channel._event_callbacks["viewer-1"]
        assert list(registration.outbox) == [
            json.dumps({"type": "ping", "requestId": "request-1"}, separators=(",", ":"))
        ]
    finally:
        channel.unregister_event_callback("viewer-1", "secret-1")


def test_command_error_names_only_commands_a_viewer_can_send() -> None:
    """The failure channel and the command channel share one command list.

    A `command_error` may only name a command Python can actually send, and
    every such command must have a public method to name in the notice --
    otherwise a viewer could either report a command that does not exist or
    report one Cellucid cannot describe.
    """
    assert set(jupyter_viewer._COMMAND_PYTHON_CALLERS) == set(wire._COMMAND_TYPES)
    for command, caller in jupyter_viewer._COMMAND_PYTHON_CALLERS.items():
        assert caller.startswith("viewer.") and caller.endswith(")"), command
        method = caller[len("viewer.") : caller.index("(")]
        assert callable(getattr(jupyter.AnnDataViewer, method)), command


def test_debug_connection_exposes_only_current_diagnostics() -> None:
    parameters = inspect.signature(jupyter.AnnDataViewer.debug_connection).parameters
    assert list(parameters) == ["self", "timeout"]


def test_debug_connection_reports_the_commands_the_viewer_rejected(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = serverless_viewer(
        port=8765,
        client_server_url="http://127.0.0.1:8765",
        web_cache_dir=tmp_path / "web-cache",
    )
    viewer.register_hook("command_error", lambda _event: None)

    def unreachable(*_args: Any, **_kwargs: Any):
        raise OSError("diagnostic probes are out of scope for this test")

    monkeypatch.setattr("urllib.request.urlopen", unreachable)

    limit = jupyter_viewer._DEBUG_COMMAND_ERROR_LIMIT
    for index in range(limit + 2):
        viewer._handle_frontend_message(
            {
                "type": "command_error",
                "command": "setColorBy",
                "reason": f"failure {index}",
            }
        )
    viewer._handle_frontend_message(
        {
            "type": "ready",
            "n_cells": 3,
            "dimensions": 2,
        }
    )

    recent = viewer.debug_connection()["recent_events"]
    assert recent["count"] == limit + 3
    assert recent["by_type"] == {"command_error": limit + 2, "ready": 1}
    # The most recent failures, oldest first, bounded so one key stays readable.
    assert recent["command_errors"] == [
        {"command": "setColorBy", "reason": f"failure {index}"}
        for index in range(2, limit + 2)
    ]
