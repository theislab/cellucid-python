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


def test_actions_reject_coercion_before_publication(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = serverless_viewer()
    published: list[dict[str, Any]] = []
    monkeypatch.setattr(viewer, "send_message", published.append)

    with pytest.raises(TypeError, match="not a mix"):
        viewer.highlight_cells(cells=[0, True])
    with pytest.raises(ValueError, match="at least 1 cell"):
        viewer.highlight_cells(cells=[])
    with pytest.raises(ValueError, match="more than once"):
        viewer.highlight_cells(cells=[0, 0])
    with pytest.raises(TypeError, match="not a mix"):
        viewer.cells_visibility_set(cells=[0, 1.5], visible=True)  # type: ignore[list-item]
    with pytest.raises(TypeError, match="one string"):
        viewer.cells_visibility_set(cells="0", visible=True)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="rows, a boolean mask or obs names"):
        viewer.highlight_cells(cells=np.array([0.0, 1.0]))
    with pytest.raises(ValueError, match="non-negative"):
        viewer.highlight_cells(cells=np.array([3, -1]))
    with pytest.raises(ValueError, match="beyond 2\\*\\*32"):
        viewer.highlight_cells(cells=np.array([2**32]))
    with pytest.raises(TypeError, match="has no parameter 'colour'"):
        viewer.act("highlight.cells", cells=[0], colour="#00ff00")
    with pytest.raises(TypeError, match="needs 'field'"):
        viewer.act("color.by")
    with pytest.raises(ValueError, match="did you mean 'color.by'"):
        viewer.act("colour.by", field=None)
    # Nothing was sent, and nothing waits for an answer.
    assert published == []
    assert viewer._action_waiters == {}

    viewer.highlight_cells(cells=[2, 0], color="#00ff00")
    viewer.highlight_cells(cells=np.array([7, 5], dtype=np.int64), color="#0000ff", label="Cycling")
    viewer.color_by(field={"source": "obs", "key": "cell_type"})
    viewer.cells_visibility_set(cells=None, visible=False)
    viewer.cells_visibility_set(cells=np.array([4, 1]), visible=True)
    viewer.camera_reset(view="live")
    # Cells go as uint32 original positions, in the order given; the rest as JSON.
    params = [entry["params"] for entry in published]
    cells = [entry.pop("cells") for entry in params if isinstance(entry.get("cells"), np.ndarray)]
    assert [(value.dtype.str, value.tolist()) for value in cells] == [
        ("<u4", [2, 0]),
        ("<u4", [7, 5]),
        ("<u4", [4, 1]),
    ]
    assert [(entry["type"], entry["action"]) for entry in published] == [
        ("action", "highlight.cells"),
        ("action", "highlight.cells"),
        ("action", "color.by"),
        ("action", "cells.visibility.set"),
        ("action", "cells.visibility.set"),
        ("action", "camera.reset"),
    ]
    assert params == [
        {"color": "#00ff00", "tool": "annotation"},
        {"color": "#0000ff", "label": "Cycling", "tool": "annotation"},
        {"field": {"source": "obs", "key": "cell_type"}},
        {"cells": None, "visible": False},
        {"visible": True},
        {"view": "live"},
    ]
    assert len({entry["requestId"] for entry in published}) == 6
    assert len(viewer._action_waiters) == 6


def test_a_browser_tab_viewer_shows_a_link_to_its_tab_and_no_frame() -> None:
    # Inside a notebook cell the viewer is the widget (cellucid.jupyter.Widget);
    # this viewer runs in a tab of its own, and its output only links to it.
    viewer = _undisplayed_viewer()
    viewer._client_server_url = "https://viewer.example"

    html = viewer._repr_html_()
    assert "<script" not in html
    assert "<iframe" not in html
    assert 'href="https://viewer.example/?viewer=true&amp;viewerId=viewer-1&amp;' in html
    assert 'target="_blank"' in html


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
        if command.startswith("agent_"):
            # An agent app's commands: sent by `cellucid mcp`, whose failures
            # reach the agent that sent them, never a notice in the notebook.
            assert caller == "cellucid mcp", command
            continue
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
                "command": "action",
                "requestId": f"request-{index}",
                "reason": f"failure {index}",
            }
        )
    viewer._handle_frontend_message(
        {
            "type": "ready",
            "n_cells": 3,
            "n_obs": 3,
            "dimensions": 2,
        }
    )

    recent = viewer.debug_connection()["recent_events"]
    assert recent["count"] == limit + 3
    assert recent["by_type"] == {"command_error": limit + 2, "ready": 1}
    # The most recent failures, oldest first, bounded so one key stays readable.
    assert recent["command_errors"] == [
        {"command": "action", "requestId": f"request-{index}", "reason": f"failure {index}"}
        for index in range(2, limit + 2)
    ]
