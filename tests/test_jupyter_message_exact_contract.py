from __future__ import annotations

import inspect
import json
from typing import Any

import pytest
from _serverless_viewer import serverless_viewer

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


def test_command_error_names_only_commands_a_viewer_can_send() -> None:
    """The failure channel and the command channel share one command list.

    A `command_error` may only name a command Python can actually send, and
    every such command must have a public method to name in the notice --
    otherwise a viewer could either report a command that does not exist or
    report one Cellucid cannot describe.
    """
    assert set(jupyter_viewer._COMMAND_PYTHON_CALLERS) == set(wire._COMMAND_MESSAGE_FIELDS)
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
                "viewerId": viewer._viewer_id,
                "command": "setColorBy",
                "reason": f"failure {index}",
            }
        )
    viewer._handle_frontend_message(
        {
            "type": "ready",
            "viewerId": viewer._viewer_id,
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
