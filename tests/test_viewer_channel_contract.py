from __future__ import annotations

import inspect
import json
import struct
import tempfile
import threading
import time
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from _serverless_viewer import serverless_viewer
from _viewer_channel_client import (
    BINARY,
    TEXT,
    FakeViewer,
    HandshakeRefused,
    binary_event,
)

import cellucid._viewer_channel as channel
import cellucid._websocket as websocket_module
import cellucid.jupyter as jupyter
import cellucid.jupyter._viewer as jupyter_viewer
import cellucid.jupyter._wire as wire
from cellucid._server_base import CellucidThreadingHTTPServer
from cellucid._websocket import unmask_in_place
from cellucid.anndata_server import AnnDataServer
from cellucid.session_bundle import SESSION_BUNDLE_MAGIC

VIEWER_ID = "security-contract-viewer"
VIEWER_TOKEN = "security-contract-viewer-token"
REQUEST_ID = "security-contract-request"

#: The digest of ``tests/fixtures/viewer_protocol.json``: the web app's
#: ``tests/fixtures/viewer-protocol.json`` byte for byte. The web app and
#: cellucid-r pin the same one.
VIEWER_PROTOCOL_SHA256 = "642c04e1d81b574e72a5d8b3738858b3447175c2c395500d3556c8198279a9e5"


def _known_inbound_events(viewer_id: str) -> list[tuple[str, dict[str, object], dict[str, object]]]:
    return [
        (
            "selection",
            {
                "type": "selection",
                "cells": bytearray(struct.pack("<2I", 0, 2)),
                "source": "lasso",
            },
            {"cells": [0, 2], "source": "lasso"},
        ),
        (
            "hover",
            {
                "type": "hover",
                "cell": 2,
                "position": {"x": 0.25, "y": -0.5, "z": 1.0},
            },
            {"cell": 2, "position": {"x": 0.25, "y": -0.5, "z": 1.0}},
        ),
        (
            "click",
            {
                "type": "click",
                "cell": 1,
                "button": 0,
                "shift": False,
                "ctrl": True,
            },
            {"cell": 1, "button": 0, "shift": False, "ctrl": True},
        ),
        (
            "ready",
            {
                "type": "ready",
                "n_cells": 3,
                "dimensions": 2,
            },
            {"n_cells": 3, "dimensions": 2},
        ),
        (
            "pong",
            {
                "type": "pong",
                "requestId": "ping-request",
                "t": 1_722_000_000_000,
            },
            {"requestId": "ping-request", "t": 1_722_000_000_000},
        ),
        (
            "debug_snapshot",
            {
                "type": "debug_snapshot",
                "requestId": "debug-request",
                "ts": "2026-07-26T12:00:00.000Z",
                "locationHref": "https://viewer.example/?viewer=true",
                "origin": "https://viewer.example",
                "serverUrl": "https://notebook.example/cellucid",
                "connected": True,
                "userAgent": "Synthetic Browser",
            },
            {
                "requestId": "debug-request",
                "ts": "2026-07-26T12:00:00.000Z",
                "locationHref": "https://viewer.example/?viewer=true",
                "origin": "https://viewer.example",
                "serverUrl": "https://notebook.example/cellucid",
                "connected": True,
                "userAgent": "Synthetic Browser",
            },
        ),
        (
            "session_bundle",
            {
                "type": "session_bundle",
                "requestId": "session-request",
                "status": "ok",
                "bytes": 1024,
                "path": "/tmp/session bundle.cellucid-session",
            },
            {
                "requestId": "session-request",
                "status": "ok",
                "bytes": 1024,
                "path": "/tmp/session bundle.cellucid-session",
            },
        ),
        (
            "command_error",
            {
                "type": "command_error",
                "command": "setColorBy",
                "reason": "Field CD8A is not loaded",
            },
            {"command": "setColorBy", "reason": "Field CD8A is not loaded"},
        ),
        (
            "data_mapping",
            {
                "type": "data_mapping",
                "requestId": None,
                "origin": "default",
                "mapping": {
                    "kind": "anndata",
                    "embeddings": {"enabled": [{"key": "X_umap", "columns": [0, 1]}], "defaults": {"1d": None, "2d": {"key": "X_umap", "columns": [0, 1]}, "3d": None}},
                    "knn": None,
                    "geneNames": "index",
                    "expression": "X",
                    "log1p": {},
                    "highlyVariable": None,
                    "latent": None,
                    "vectorFields": [],
                },
            },
            {
                "requestId": None,
                "origin": "default",
                "mapping": {
                    "kind": "anndata",
                    "embeddings": {"enabled": [{"key": "X_umap", "columns": [0, 1]}], "defaults": {"1d": None, "2d": {"key": "X_umap", "columns": [0, 1]}, "3d": None}},
                    "knn": None,
                    "geneNames": "index",
                    "expression": "X",
                    "log1p": {},
                    "highlyVariable": None,
                    "latent": None,
                    "vectorFields": [],
                },
            },
        ),
    ]


def _minimal_adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.array([[1.0], [2.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["gene-1"]),
        obsm={
            "X_umap_2d": np.array(
                [[0.0, 0.0], [1.0, 1.0]],
                dtype=np.float32,
            )
        },
    )


@contextmanager
def _anndata_server() -> Iterator[AnnDataServer]:
    server = AnnDataServer(
        _minimal_adata(),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Jupyter security contract",
        dataset_id="jupyter-security-contract",
    )
    server.start_background()
    assert server._server is not None
    try:
        yield server
    finally:
        server.stop()


@contextmanager
def _running_server() -> Iterator[tuple[str, int]]:
    with _anndata_server() as server:
        assert server._server is not None
        host, port = server._server.server_address
        yield host, port


# The exact current fingerprint shape: a session names the embedding its cell
# order was digested over and the data mapping it was saved with.
_FINGERPRINT = {
    "datasetId": None,
    "name": None,
    "cellCount": 0,
    "varCount": 0,
    "cellOrder": {
        "embedding": {"key": "X_umap", "columns": [0, 1]},
        "dimension": 2,
        "digest": "0123456789abcdef",
    },
    "dataMapping": {
        "kind": "anndata",
        "embeddings": {
            "enabled": [{"key": "X_umap", "columns": [0, 1]}],
            "defaults": {"1d": None, "2d": {"key": "X_umap", "columns": [0, 1]}, "3d": None},
        },
        "knn": None,
        "geneNames": "index",
        "expression": "X",
        "log1p": {},
        "highlyVariable": None,
        "latent": None,
        "vectorFields": [],
    },
}


def _empty_session_bundle() -> bytes:
    manifest = {
        "createdAt": "2026-07-26T00:00:00.000Z",
        "datasetFingerprint": _FINGERPRINT,
        "source": None,
        "requires": {"obs": [], "var": []},
        "chunks": [],
    }
    manifest_bytes = json.dumps(
        manifest,
        separators=(",", ":"),
    ).encode("utf-8")
    return (
        SESSION_BUNDLE_MAGIC + struct.pack("<I", len(manifest_bytes)) + manifest_bytes
    )


def _corrupt_gzip_session_bundle() -> bytes:
    stored = b"not-a-gzip-stream"
    manifest = {
        "createdAt": "2026-07-26T00:00:00.000Z",
        "datasetFingerprint": _FINGERPRINT,
        "source": None,
        "requires": {"obs": [], "var": []},
        "chunks": [
            {
                "id": "corrupt",
                "contributorId": "corrupt",
                "priority": "eager",
                "kind": "binary",
                "codec": "gzip",
                "label": "Corrupt payload",
                "datasetDependent": False,
                "storedBytes": len(stored),
                "uncompressedBytes": 4,
            }
        ],
    }
    manifest_bytes = json.dumps(
        manifest,
        separators=(",", ":"),
    ).encode("utf-8")
    return (
        SESSION_BUNDLE_MAGIC
        + struct.pack("<I", len(manifest_bytes))
        + manifest_bytes
        + struct.pack("<I", len(stored))
        + stored
    )


@pytest.fixture(autouse=True)
def _restore_routing_registries() -> Iterator[None]:
    event_snapshot = dict(channel._event_callbacks)
    pending_snapshot = dict(channel._pending_session_bundle_requests)
    yield
    channel._event_callbacks.clear()
    channel._event_callbacks.update(event_snapshot)
    channel._pending_session_bundle_requests.clear()
    channel._pending_session_bundle_requests.update(pending_snapshot)


def test_routing_functions_expose_only_the_authenticated_current_contract() -> None:
    assert list(inspect.signature(channel.register_event_callback).parameters) == [
        "viewer_id",
        "viewer_token",
        "callback",
    ]
    assert list(inspect.signature(channel.unregister_event_callback).parameters) == [
        "viewer_id",
        "viewer_token",
    ]
    assert list(inspect.signature(channel.route_event).parameters) == [
        "viewer_id",
        "viewer_token",
        "event",
    ]
    assert list(inspect.signature(channel.register_session_bundle_request).parameters) == [
        "viewer_id",
        "viewer_token",
        "request_id",
        "ttl_seconds",
    ]
    assert list(inspect.signature(channel.cancel_session_bundle_request).parameters) == [
        "viewer_id",
        "viewer_token",
        "request_id",
    ]


@pytest.mark.parametrize(
    "ttl_seconds",
    [True, "60", 0, -1.0, float("nan"), float("inf")],
)
def test_pending_request_registration_rejects_ttl_coercion_or_repair(
    ttl_seconds: object,
) -> None:
    channel.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        lambda _event: None,
    )
    before = dict(channel._pending_session_bundle_requests)

    with pytest.raises((TypeError, ValueError), match="ttl_seconds"):
        channel.register_session_bundle_request(
            VIEWER_ID,
            VIEWER_TOKEN,
            REQUEST_ID,
            ttl_seconds=ttl_seconds,  # type: ignore[arg-type]
        )

    assert channel._pending_session_bundle_requests == before


def test_a_request_without_a_deadline_stays_answerable_until_cancelled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """``get_session_bundle(timeout=None)`` waits forever, so its request must too."""
    channel.register_event_callback(VIEWER_ID, VIEWER_TOKEN, lambda _event: None)
    channel.register_session_bundle_request(VIEWER_ID, VIEWER_TOKEN, REQUEST_ID, ttl_seconds=None)
    now = time.monotonic()
    # A day later another registration sweeps expired requests; this one is kept.
    monkeypatch.setattr(channel.time, "monotonic", lambda: now + 86_400.0)
    channel.register_session_bundle_request(VIEWER_ID, VIEWER_TOKEN, "other-request", ttl_seconds=60.0)
    assert (VIEWER_ID, REQUEST_ID) in channel._pending_session_bundle_requests
    channel.cancel_session_bundle_request(VIEWER_ID, VIEWER_TOKEN, REQUEST_ID)
    assert (VIEWER_ID, REQUEST_ID) not in channel._pending_session_bundle_requests


def test_pending_request_registration_rejects_wrong_token_and_duplicate() -> None:
    channel.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        lambda _event: None,
    )

    with pytest.raises(PermissionError):
        channel.register_session_bundle_request(
            VIEWER_ID,
            "wrong-viewer-token",
            REQUEST_ID,
        )
    assert channel._pending_session_bundle_requests == {}

    channel.register_session_bundle_request(
        VIEWER_ID,
        VIEWER_TOKEN,
        REQUEST_ID,
    )
    pending = channel._pending_session_bundle_requests[(VIEWER_ID, REQUEST_ID)]
    with pytest.raises(RuntimeError, match="already pending"):
        channel.register_session_bundle_request(
            VIEWER_ID,
            VIEWER_TOKEN,
            REQUEST_ID,
        )
    assert channel._pending_session_bundle_requests == {(VIEWER_ID, REQUEST_ID): pending}


def test_a_failing_viewer_hook_never_costs_the_notebook_the_event() -> None:
    # One buggy hook used to drop the event from viewer.state and
    # wait_for_event, so get_session_bundle() timed out behind it. The event
    # is recorded first; every hook runs; the failure is raised afterwards for
    # the channel to log.
    viewer = serverless_viewer()
    ran: list[str] = []

    def fail(_event: dict[str, object]) -> None:
        ran.append("fail")
        raise RuntimeError("hook failed")

    viewer.register_hook("message", fail)
    viewer.register_hook("selection", lambda _event: ran.append("selection"))
    waited: dict[str, object] = {}
    waiter = threading.Thread(
        target=lambda: waited.update(payload=viewer.wait_for_event("selection", timeout=10))
    )
    waiter.start()
    deadline = time.monotonic() + 10
    while not viewer._waiters and time.monotonic() < deadline:
        time.sleep(0.001)
    assert viewer._waiters, "wait_for_event registered its waiter"
    with pytest.raises(RuntimeError, match="hook failed"):
        viewer._handle_frontend_message(
            {
                "type": "selection",
                "cells": struct.pack("<I", 0),
                "source": "lasso",
            }
        )
    waiter.join(timeout=10)

    assert ran == ["fail", "selection"]
    assert viewer.state.last_event_type == "selection"
    assert viewer.state.selection == viewer.state.last_event
    assert waited["payload"] == viewer.state.selection
    assert viewer._event_seq == 1
    assert [event for _seq, event, _data in viewer._recent_events] == ["selection"]

    # Several failures are raised together; the event is still recorded.
    viewer.register_hook("selection", fail)
    with pytest.raises(ExceptionGroup) as caught:
        viewer._handle_frontend_message(
            {
                "type": "selection",
                "cells": struct.pack("<I", 1),
                "source": "lasso",
            }
        )
    assert [str(error) for error in caught.value.exceptions] == ["hook failed", "hook failed"]
    assert viewer._event_seq == 2


def test_every_known_inbound_event_preserves_its_exact_payload() -> None:
    viewer = serverless_viewer()
    delivered: list[dict[str, object]] = []
    viewer.register_hook("message", delivered.append)

    known_events = _known_inbound_events(viewer._viewer_id)
    for event_name, message, expected_payload in known_events:
        viewer._handle_frontend_message(message)
        assert delivered[-1] == {"event": event_name, **expected_payload}
        assert viewer.state.last_event_type == event_name
        assert viewer.state.last_event == expected_payload

    assert len(delivered) == len(known_events)


@pytest.mark.parametrize("event_index", range(len(_known_inbound_events("probe"))))
def test_inbound_event_rejects_every_undeclared_property_before_callback_dispatch(
    event_index: int,
) -> None:
    viewer = serverless_viewer()
    delivered: list[dict[str, object]] = []
    viewer.register_hook("message", delivered.append)
    _event_name, message, _expected_payload = _known_inbound_events(viewer._viewer_id)[event_index]
    message["undeclared"] = True

    with pytest.raises(ValueError, match="unknown.*undeclared"):
        viewer._handle_frontend_message(message)

    assert delivered == []
    assert viewer.state.last_event is None
    assert viewer._event_seq == 0
    assert list(viewer._recent_events) == []


@pytest.mark.parametrize(
    "message",
    [
        {
            "type": "custom",
            "value": 1,
        },
        {
            "type": "selection",
            "cells": struct.pack("<I", 0),
        },
        # Events carry no credentials: the channel is the viewer's.
        {
            "type": "selection",
            "viewerId": "a-viewer",
            "cells": struct.pack("<I", 0),
            "source": "lasso",
        },
        {
            "type": "selection",
            "viewerToken": "must-not-reach-the-viewer",
            "cells": struct.pack("<I", 0),
            "source": "lasso",
        },
    ],
)
def test_inbound_event_rejects_unknown_types_missing_fields_and_routing_extras(
    message: dict[str, object],
) -> None:
    viewer = serverless_viewer()
    delivered: list[dict[str, object]] = []
    viewer.register_hook("message", delivered.append)
    candidate = dict(message)

    with pytest.raises(channel.ViewerProtocolError, match="event|fields|unknown"):
        viewer._handle_frontend_message(candidate)

    assert delivered == []
    assert viewer.state.last_event is None


def _command_error(viewer_id: str, **overrides: object) -> dict[str, object]:
    record: dict[str, object] = {
        "type": "command_error",
        "command": "setColorBy",
        "reason": "Jupyter set-color-by field must be exact non-empty text",
    }
    record.update(overrides)
    return record


def test_command_error_reaches_every_notebook_surface(
    capsys: pytest.CaptureFixture[str],
) -> None:
    viewer = serverless_viewer()
    handled: list[dict[str, object]] = []
    viewer.register_hook("command_error", handled.append)

    viewer._handle_frontend_message(_command_error(viewer._viewer_id))

    payload = {
        "command": "setColorBy",
        "reason": "Jupyter set-color-by field must be exact non-empty text",
    }
    assert handled == [payload]
    assert viewer.state.command_error == payload
    assert viewer.state.last_event_type == "command_error"
    assert [entry[1:] for entry in viewer._recent_events] == [
        ("command_error", payload)
    ]
    # A notebook that took the event over is not also reported to.
    assert capsys.readouterr().err == ""


def test_command_error_wakes_a_notebook_blocked_on_wait_for_event() -> None:
    viewer = serverless_viewer()
    viewer.register_hook("command_error", lambda _event: None)
    received: list[dict[str, Any]] = []

    def wait() -> None:
        received.append(viewer.wait_for_event("command_error", timeout=10.0))

    waiter = threading.Thread(target=wait)
    waiter.start()
    try:
        # `wait_for_event` returns the *next* matching event, so publish until
        # the waiter has entered the wait rather than racing it once.
        deadline = time.monotonic() + 10.0
        while waiter.is_alive() and time.monotonic() < deadline:
            viewer._handle_frontend_message(_command_error(viewer._viewer_id))
            waiter.join(timeout=0.05)
    finally:
        waiter.join(timeout=10.0)

    assert not waiter.is_alive()
    assert received == [
        {
            "command": "setColorBy",
            "reason": "Jupyter set-color-by field must be exact non-empty text",
        }
    ]


def test_unhandled_command_error_is_printed_where_the_notebook_looks(
    capsys: pytest.CaptureFixture[str],
) -> None:
    viewer = serverless_viewer()
    observed: list[dict[str, object]] = []
    viewer.register_hook("message", observed.append)

    viewer._handle_frontend_message(_command_error(viewer._viewer_id))
    viewer._handle_frontend_message(
        _command_error(
            viewer._viewer_id,
            command="highlight",
            reason="Jupyter highlight cells must be a non-empty array",
        )
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    # Named by the Python call that sent it, not by the wire type, and repeated
    # for every failure rather than deduplicated.
    assert captured.err == (
        "[cellucid] viewer.set_color_by(...) did not take effect in the viewer: "
        "Jupyter set-color-by field must be exact non-empty text\n"
        "[cellucid] The viewer is unchanged. Handle this in code with "
        "@viewer.on_command_error.\n"
        "[cellucid] viewer.highlight_cells(...) did not take effect in the viewer: "
        "Jupyter highlight cells must be a non-empty array\n"
        "[cellucid] The viewer is unchanged. Handle this in code with "
        "@viewer.on_command_error.\n"
    )
    assert [entry["event"] for entry in observed] == [
        "command_error",
        "command_error",
    ]


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"command": "notACommand"}, "one current Jupyter command"),
        ({"command": "command_error"}, "one current Jupyter command"),
        ({"command": ""}, "command_error command"),
        ({"command": None}, "command_error command"),
        ({"command": "setColorBy "}, "command_error command"),
        ({"reason": ""}, "command_error reason"),
        ({"reason": None}, "command_error reason"),
        ({"reason": "line one\nline two"}, "command_error reason"),
        # This text is printed to a terminal, so escape sequences must not
        # survive the wire.
        ({"reason": "\x1b[31mred"}, "command_error reason"),
        ({"reason": "carriage\rreturn"}, "command_error reason"),
        ({"reason": "x" * 501}, "at most 500 characters"),
        ({"requestId": "unexpected"}, "unknown requestId"),
    ],
)
def test_command_error_rejects_unnamed_commands_and_unbounded_reasons(
    overrides: dict[str, object],
    match: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    viewer = serverless_viewer()
    delivered: list[dict[str, object]] = []
    viewer.register_hook("message", delivered.append)

    with pytest.raises((TypeError, ValueError), match=match):
        viewer._handle_frontend_message(
            _command_error(viewer._viewer_id, **overrides)
        )

    assert delivered == []
    assert viewer.state.command_error is None
    assert list(viewer._recent_events) == []
    assert capsys.readouterr().err == ""


def test_command_error_accepts_the_longest_declared_reason() -> None:
    viewer = serverless_viewer()
    handled: list[dict[str, object]] = []
    viewer.register_hook("command_error", handled.append)

    viewer._handle_frontend_message(
        _command_error(viewer._viewer_id, reason="x" * 500)
    )

    assert handled == [{"command": "setColorBy", "reason": "x" * 500}]


def test_the_published_events_are_exactly_the_events_a_viewer_may_send() -> None:
    """``/_cellucid/protocol`` publishes what the validator accepts, no more.

    ``_known_inbound_events`` is the independent list here: every entry in it is
    driven through ``_handle_frontend_message`` by
    ``test_every_known_inbound_event_preserves_its_exact_payload``, so tying the
    published set to it means the route can only name types this installation
    demonstrably accepts, and can only omit one by leaving a type untested.
    """
    published = wire._protocol_capability_document()["events"]

    assert set(published) == {name for name, _message, _payload in _known_inbound_events(VIEWER_ID)}


def test_no_published_event_is_refused_by_the_validator_as_an_unknown_type() -> None:
    for event_type in wire._protocol_capability_document()["events"]:
        with pytest.raises((TypeError, ValueError)) as rejection:
            wire._require_inbound_jupyter_event({"type": event_type})
        # The stripped-down probe must fail on its missing fields, never on the
        # type: a type rejection here would mean the route advertised an event
        # the channel then closes on as a protocol violation.
        assert "Unknown inbound viewer event type" not in str(rejection.value)

    with pytest.raises(ValueError, match="Unknown inbound viewer event type"):
        wire._require_inbound_jupyter_event({"type": "command_success"})


def test_the_published_commands_are_exactly_the_commands_python_can_send() -> None:
    published = wire._protocol_capability_document()["commands"]

    # `_COMMAND_PYTHON_CALLERS` is declared in a different module, and every key
    # of it names a real `AnnDataViewer` method, so a command published here is one
    # a notebook can actually issue.
    assert set(published) == set(jupyter_viewer._COMMAND_PYTHON_CALLERS)


def test_every_published_command_may_be_named_by_a_command_error() -> None:
    viewer = serverless_viewer()
    handled: list[dict[str, object]] = []
    viewer.register_hook("command_error", handled.append)
    published = wire._protocol_capability_document()["commands"]

    for command in published:
        viewer._handle_frontend_message(
            _command_error(viewer._viewer_id, command=command)
        )

    assert [entry["command"] for entry in handled] == list(published)

    with pytest.raises(ValueError, match="one current Jupyter command"):
        viewer._handle_frontend_message(
            _command_error(viewer._viewer_id, command="setPointSize")
        )


def _routing_snapshot() -> tuple[dict[str, object], set[object]]:
    return dict(channel._event_callbacks), set(jupyter_viewer._active_viewers)


def test_server_start_failure_rolls_back_viewer_registration_and_server_resources(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    instances: list[object] = []

    class FailingServer:
        def __init__(self, **_kwargs: object):
            self.port = 43123
            self.stop_calls = 0
            instances.append(self)

        def start_background(self) -> None:
            raise RuntimeError("server start failed")

        def stop(self) -> None:
            self.stop_calls += 1

    before = _routing_snapshot()
    import cellucid.anndata_server as anndata_server

    monkeypatch.setattr(anndata_server, "AnnDataServer", FailingServer)

    with pytest.raises(RuntimeError, match="server start failed"):
        jupyter.AnnDataViewer(
            _minimal_adata(),
            auto_open=False,
            dataset_name="Failed viewer",
            dataset_id="failed-viewer",
        )

    assert len(instances) == 1
    assert instances[0].stop_calls == 1
    assert _routing_snapshot() == before


def test_display_failure_rolls_back_started_server_and_viewer_registries(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    instances: list[object] = []

    class StartedServer:
        def __init__(self, **_kwargs: object):
            self.port = 43123
            self.stop_calls = 0
            instances.append(self)

        @property
        def url(self) -> str:
            return f"http://127.0.0.1:{self.port}"

        def start_background(self) -> None:
            return None

        def stop(self) -> None:
            self.stop_calls += 1

    def fail_display(_viewer: object) -> None:
        raise RuntimeError("display failed")

    # Each patch names the module where the caller resolves the symbol:
    # ``AnnDataServer`` in the server module the viewer imports it from,
    # ``_detect_jupyter_context`` in the viewer module whose ``__init__`` calls it.
    import cellucid.anndata_server as anndata_server

    monkeypatch.setattr(anndata_server, "AnnDataServer", StartedServer)
    monkeypatch.setattr(
        jupyter_viewer,
        "_detect_jupyter_context",
        lambda: {"in_jupyter": True},
    )
    monkeypatch.setattr(jupyter.AnnDataViewer, "display", fail_display)
    before = _routing_snapshot()

    with pytest.raises(RuntimeError, match="display failed"):
        jupyter.AnnDataViewer(
            _minimal_adata(),
            auto_open=True,
            dataset_name="Failed viewer",
            dataset_id="failed-viewer",
        )

    assert len(instances) == 1
    assert instances[0].stop_calls == 1
    assert _routing_snapshot() == before


def test_a_waited_for_event_survives_any_burst_of_later_events() -> None:
    # The event is handed to the waiter as it is recorded, so the 10,000
    # hover events behind it cannot push it out of a bounded log first.
    viewer = serverless_viewer()
    received: list[dict[str, Any]] = []

    def wait() -> None:
        received.append(viewer.wait_for_event("selection", timeout=10.0))

    waiter = threading.Thread(target=wait)
    waiter.start()
    deadline = time.monotonic() + 5.0
    while not viewer._waiters and time.monotonic() < deadline:
        time.sleep(0.01)
    assert viewer._waiters
    # The waiter is not woken until the whole burst is recorded.
    condition = viewer._event_cv
    notify_all = condition.notify_all
    condition.notify_all = lambda: None  # type: ignore[method-assign]
    viewer._record_event("selection", {"cells": [3], "source": "lasso"})
    for index in range(10_000):
        viewer._record_event("hover", {"cell": index})
    condition.notify_all = notify_all  # type: ignore[method-assign]
    with condition:
        condition.notify_all()
    waiter.join(timeout=10.0)
    assert received == [{"cells": [3], "source": "lasso"}]


# =============================================================================
# The WebSocket: handshake, authentication, commands, events, failures
# =============================================================================


def _wait_until(predicate: Any, timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while not predicate():
        if time.monotonic() > deadline:
            raise AssertionError("condition not reached in time")
        time.sleep(0.005)


@contextmanager
def _channel_viewer() -> Iterator[tuple[jupyter.AnnDataViewer, int]]:
    """A registered notebook viewer and a running data server, no browser."""
    viewer = serverless_viewer()
    viewer._activate()
    with _running_server() as (_host, port):
        try:
            yield viewer, port
        finally:
            viewer.stop()


def _connect(port: int, viewer: jupyter.AnnDataViewer, **kwargs: Any) -> FakeViewer:
    socket_ = FakeViewer(port, **kwargs)
    assert socket_.hello(viewer._viewer_id, viewer._viewer_token) == {"type": "welcome"}
    return socket_


@pytest.mark.parametrize(
    ("headers", "status"),
    [
        ({"Upgrade": None}, 426),
        ({"Sec-WebSocket-Version": "8"}, 426),
        ({"Sec-WebSocket-Key": "c2hvcnQ="}, 400),
    ],
)
def test_the_channel_refuses_a_malformed_upgrade(
    headers: dict[str, str | None], status: int
) -> None:
    import http.client

    request_headers = {
        "Upgrade": "websocket",
        "Connection": "Upgrade",
        "Sec-WebSocket-Version": "13",
        "Sec-WebSocket-Key": "dGhlIHNhbXBsZSBub25jZQ==",
        "Origin": "http://127.0.0.1",
    }
    for name, value in headers.items():
        if value is None:
            request_headers.pop(name)
        else:
            request_headers[name] = value
    with _running_server() as (host, port):
        connection = http.client.HTTPConnection(host, port, timeout=5)
        connection.request("GET", "/_cellucid/ws", headers=request_headers)
        assert connection.getresponse().status == status
        connection.close()


@pytest.mark.parametrize(
    ("origin", "allowed"),
    [
        ("http://127.0.0.1:9999", True),
        ("http://localhost:3000", True),
        ("http://[::1]:8888", True),
        ("https://evil.example", False),
        ("null", False),
        ("", False),
    ],
)
def test_only_this_servers_own_and_loopback_pages_open_the_channel(
    origin: str, allowed: bool
) -> None:
    with _channel_viewer() as (viewer, port):
        if allowed:
            _connect(port, viewer, origin=origin).close()
            return
        with pytest.raises(HandshakeRefused) as refusal:
            FakeViewer(port, origin=origin or "")
        assert refusal.value.status == 403


def test_origin_rule_accepts_the_same_host_and_declared_proxy_names() -> None:
    allowed = frozenset({"hub.example.org"})
    assert channel.origin_is_allowed("https://hub.example.org", "127.0.0.1:8765", allowed)
    assert channel.origin_is_allowed("https://lab.example.org", "lab.example.org", frozenset())
    assert not channel.origin_is_allowed("https://evil.example", "lab.example.org", allowed)
    assert not channel.origin_is_allowed("ftp://127.0.0.1", "127.0.0.1:8765", frozenset())
    assert not channel.origin_is_allowed(None, "127.0.0.1:8765", frozenset())


def test_a_foreign_host_header_is_refused_before_the_upgrade() -> None:
    with _channel_viewer() as (_viewer, port):
        with pytest.raises(HandshakeRefused) as refusal:
            FakeViewer(port, host="attacker.example")
        assert refusal.value.status == 421


@pytest.mark.parametrize(
    "hello",
    [
        {"type": "hello", "viewerId": "unknown", "viewerToken": "x" * 32},
        {"type": "hello", "viewerId": "REPLACE", "viewerToken": "wrong-token"},
        {"type": "hello", "viewerId": "REPLACE", "viewerToken": ""},
    ],
)
def test_unknown_or_wrong_credentials_close_the_channel_unauthorized(
    hello: dict[str, str],
) -> None:
    with _channel_viewer() as (viewer, port):
        socket_ = FakeViewer(port)
        if hello["viewerId"] == "REPLACE":
            hello = {**hello, "viewerId": viewer._viewer_id}
        socket_.send_json(hello)
        assert socket_.receive_close() == (4401, "Unknown viewer, or a wrong viewer token")
        assert not channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token)


@pytest.mark.parametrize(
    "first",
    [
        {"type": "hello", "viewerId": "x"},
        {"type": "hello", "viewerId": "x", "viewerToken": "y", "extra": 1},
        {"type": "click", "cell": 0},
    ],
)
def test_the_first_frame_must_be_exactly_hello(first: dict[str, object]) -> None:
    with _channel_viewer() as (_viewer, port):
        socket_ = FakeViewer(port)
        socket_.send_json(first)
        code, _reason = socket_.receive_close()
        assert code == 1008


def test_commands_sent_before_the_viewer_connects_arrive_in_order_once() -> None:
    with _channel_viewer() as (viewer, port):
        assert viewer.send_message({"type": "setColorBy", "field": "leiden"}) is False
        viewer.highlight_cells([2, 0], color="#00ff00")
        viewer.set_visibility(None, visible=False)
        viewer.reset_view()
        socket_ = _connect(port, viewer)
        assert socket_.receive_json() == {"type": "setColorBy", "field": "leiden"}
        # Commands that name cells are binary frames, their cells ascending.
        assert socket_.receive_command() == {"type": "highlight", "color": "#00ff00", "cells": [0, 2]}
        assert socket_.receive_command() == {"type": "setVisibility", "visible": False, "all": True, "cells": []}
        assert socket_.receive_json() == {"type": "resetCamera"}
        assert viewer.send_message({"type": "clearHighlights"}) is True
        assert socket_.receive_json() == {"type": "clearHighlights"}
        # A reconnect gets only what was sent while it was away.
        socket_.sock.close()
        _wait_until(lambda: not channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token))
        viewer.reset_view()
        again = _connect(port, viewer)
        assert again.receive_json() == {"type": "resetCamera"}
        assert viewer.send_message({"type": "clearHighlights"}) is True
        assert again.receive_json() == {"type": "clearHighlights"}


def test_a_newer_socket_of_the_same_viewer_replaces_the_open_one() -> None:
    with _channel_viewer() as (viewer, port):
        first = _connect(port, viewer)
        second = _connect(port, viewer)
        assert first.receive_close() == (4409, "Replaced by a newer connection of this viewer")
        viewer.reset_view()
        assert second.receive_json() == {"type": "resetCamera"}


def test_events_reach_the_hooks_in_order_and_the_state() -> None:
    with _channel_viewer() as (viewer, port):
        hovers: list[dict[str, Any]] = []
        clicks: list[dict[str, Any]] = []
        viewer.register_hook("hover", hovers.append)
        viewer.register_hook("click", clicks.append)
        socket_ = _connect(port, viewer)
        for cell in range(50):
            socket_.send_json({"type": "hover", "cell": cell, "position": None})
        socket_.send_json({"type": "click", "cell": 7, "button": 0, "shift": False, "ctrl": True})
        socket_.send_json({"type": "ready", "n_cells": 2, "dimensions": 2})
        assert viewer.wait_for_ready(timeout=5) == {"n_cells": 2, "dimensions": 2}
        assert [event["cell"] for event in hovers] == list(range(50))
        assert clicks == [{"cell": 7, "button": 0, "shift": False, "ctrl": True}]


def test_a_fragmented_event_with_a_ping_between_its_frames_is_one_event() -> None:
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        payload = json.dumps({"type": "ready", "n_cells": 2, "dimensions": 3}).encode()
        socket_.sock.sendall(
            socket_.frame(TEXT, payload[:7], fin=False)
            + socket_.frame(0x9, b"are you there")
            + socket_.frame(0x0, payload[7:15], fin=False)
            + socket_.frame(0x0, payload[15:])
        )
        opcode, pong = socket_.receive()
        assert (opcode, pong) == (0xA, b"are you there")
        assert viewer.wait_for_ready(timeout=5) == {"n_cells": 2, "dimensions": 3}


@pytest.mark.parametrize(
    ("frames", "code"),
    [
        # An event the viewer may only send as a binary frame.
        ([(TEXT, json.dumps({"type": "selection", "cells": [0], "source": "lasso"}).encode())], 1008),
        # An event this installation does not declare.
        ([(TEXT, json.dumps({"type": "zoom", "level": 2}).encode())], 1008),
        # An event that names its viewer: events carry no credentials.
        ([(TEXT, json.dumps({"type": "ready", "viewerId": "x", "n_cells": 1, "dimensions": 2}).encode())], 1008),
        ([(TEXT, b"{not json")], 1008),
        ([(TEXT, b"\xff\xfe")], 1007),
        ([(TEXT, b"x" * (channel.MAX_EVENT_TEXT_BYTES + 1))], 1009),
        # A continuation frame outside a message.
        ([(0x0, b"{}")], 1002),
        # A binary event whose header is not one of the two.
        ([(BINARY, binary_event({"type": "hover"}, b""))], 1008),
        # A selection body that is not ascending uint32.
        ([(BINARY, binary_event({"type": "selection", "source": "lasso"}, struct.pack("<2I", 3, 1)))], 1008),
        ([(BINARY, binary_event({"type": "selection", "source": "lasso"}, b"\x01\x02\x03"))], 1008),
        ([(BINARY, binary_event({"type": "selection", "source": "lasso"}, b""))], 1008),
    ],
)
def test_a_frame_that_breaks_the_protocol_closes_the_channel_without_a_hook(
    frames: list[tuple[int, bytes]], code: int
) -> None:
    with _channel_viewer() as (viewer, port):
        delivered: list[dict[str, Any]] = []
        viewer.register_hook("message", delivered.append)
        socket_ = _connect(port, viewer)
        for opcode, payload in frames:
            # The server refuses an oversized frame from its header and closes,
            # but reads the rest of it: the whole frame is sent, never reset.
            socket_.send_frame(opcode, payload)
        assert socket_.receive_close()[0] == code
        assert delivered == []
        _wait_until(lambda: not channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token))


def _open_connections(server: CellucidThreadingHTTPServer) -> int:
    with server._connections_lock:
        return len(server._connections)


@contextmanager
def _viewer_on_a_server() -> Iterator[tuple[jupyter.AnnDataViewer, AnnDataServer, CellucidThreadingHTTPServer, int]]:
    viewer = serverless_viewer()
    viewer._activate()
    try:
        with _anndata_server() as server:
            http_server = server._server
            assert http_server is not None
            yield viewer, server, http_server, http_server.server_address[1]
    finally:
        viewer.stop()


def _close_the_channel(
    port: int, viewer: jupyter.AnnDataViewer, closed_by: str
) -> tuple[FakeViewer, int, list[FakeViewer]]:
    """Open a channel and have the server close it: from the connection's own
    thread (a protocol violation) or from another one (a newer connection)."""
    socket_ = _connect(port, viewer)
    if closed_by == "its own thread":
        socket_.send_json({"type": "zoom", "level": 2})
        return socket_, 1008, []
    return socket_, 4409, [_connect(port, viewer)]


@pytest.mark.parametrize("closed_by", ["its own thread", "another thread"])
def test_a_closed_channel_ends_when_the_page_closes_its_side(
    monkeypatch: pytest.MonkeyPatch, closed_by: str
) -> None:
    # Far longer than the test runs: the connection must end on the page's
    # own close, not on the bound.
    monkeypatch.setattr(websocket_module, "CLOSE_DRAIN_SECONDS", 120.0)
    with _viewer_on_a_server() as (viewer, _server, http_server, port):
        socket_, code, others = _close_the_channel(port, viewer, closed_by)
        assert socket_.receive_close()[0] == code
        # The server half-closed and is reading what the page still sends.
        time.sleep(0.2)
        assert _open_connections(http_server) == 1 + len(others)
        socket_.answer_close()
        _wait_until(lambda: _open_connections(http_server) == len(others))


@pytest.mark.parametrize("closed_by", ["its own thread", "another thread"])
def test_a_page_that_never_closes_its_side_is_cut_off_after_the_drain(
    monkeypatch: pytest.MonkeyPatch, closed_by: str
) -> None:
    monkeypatch.setattr(websocket_module, "CLOSE_DRAIN_SECONDS", 0.2)
    with _viewer_on_a_server() as (viewer, _server, http_server, port):
        socket_, code, others = _close_the_channel(port, viewer, closed_by)
        assert socket_.receive_close()[0] == code
        # The page keeps its socket open and says nothing.
        _wait_until(lambda: _open_connections(http_server) == len(others))


def test_stopping_the_server_closes_an_open_channel_going_away() -> None:
    with _viewer_on_a_server() as (viewer, server, http_server, port):
        socket_ = _connect(port, viewer)
        server.stop()
        assert socket_.receive_close() == (1001, "The server stopped")
        socket_.answer_close()
        _wait_until(lambda: _open_connections(http_server) == 0)


def test_an_unmasked_frame_is_a_protocol_error() -> None:
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        socket_.sock.sendall(struct.pack("!BB", 0x81, 2) + b"{}")
        assert socket_.receive_close()[0] == 1002


def test_a_socket_closed_mid_frame_frees_the_viewer_for_its_next_connection() -> None:
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        frame = socket_.frame(BINARY, binary_event({"type": "selection", "source": "lasso"}, b"\0" * 4000))
        socket_.sock.sendall(frame[:1000])
        socket_.sock.close()
        _wait_until(lambda: not channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token))
        assert viewer.state.selection is None
        viewer.reset_view()
        again = _connect(port, viewer)
        assert again.receive_json() == {"type": "resetCamera"}


def test_a_failing_notebook_hook_keeps_the_channel_open() -> None:
    with _channel_viewer() as (viewer, port):
        def fail(_event: dict[str, Any]) -> None:
            raise RuntimeError("the notebook's hook failed")

        viewer.register_hook("click", fail)
        socket_ = _connect(port, viewer)
        socket_.send_json({"type": "click", "cell": 1, "button": 0, "shift": False, "ctrl": False})
        socket_.send_json({"type": "ready", "n_cells": 2, "dimensions": 2})
        assert viewer.wait_for_ready(timeout=5) == {"n_cells": 2, "dimensions": 2}
        assert channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token)


def test_a_large_selection_arrives_whole_as_one_binary_frame() -> None:
    cells = np.arange(0, 6_000_000, 3, dtype="<u4")  # two million cells, 8 MB
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        started = time.perf_counter()
        socket_.send_binary(binary_event({"type": "selection", "source": "lasso"}, cells.tobytes()))
        event = viewer.wait_for_event("selection", timeout=30)
        elapsed = time.perf_counter() - started
        assert event["source"] == "lasso"
        assert len(event["cells"]) == cells.size
        assert event["cells"][:3] == [0, 3, 6]
        assert event["cells"][-1] == int(cells[-1])
        assert viewer.state.selection is not None
        print(f"[viewer channel] {cells.size:,} cells delivered in {elapsed * 1000:.0f} ms")


def test_a_selection_fragmented_across_frames_is_read_as_one_message() -> None:
    body = binary_event({"type": "selection", "source": "knn"}, struct.pack("<3I", 1, 5, 9))
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        socket_.sock.sendall(
            socket_.frame(BINARY, body[:3], fin=False)
            + socket_.frame(0x9, b"")
            + socket_.frame(0x0, body[3:11], fin=False)
            + socket_.frame(0x0, body[11:])
        )
        assert socket_.receive()[0] == 0xA
        _wait_until(lambda: viewer.state.selection is not None)
        assert viewer.state.selection == {"cells": [1, 5, 9], "source": "knn"}


def test_get_session_bundle_streams_the_frame_to_a_file_once() -> None:
    bundle = _empty_session_bundle()
    with _channel_viewer() as (viewer, port):
        viewer._displayed = True
        viewer.state.ready = {"n_cells": 2, "dimensions": 2}
        socket_ = _connect(port, viewer)
        result: list[Any] = []
        worker = threading.Thread(target=lambda: result.append(viewer.get_session_bundle(timeout=10)))
        worker.start()
        request = socket_.receive_json()
        assert request["type"] == "requestSessionBundle"
        socket_.send_binary(binary_event({"type": "session_bundle", "requestId": request["requestId"]}, bundle))
        worker.join(timeout=10)
        assert result and result[0].path.read_bytes() == bundle
        # The request was consumed: a second copy is dropped and the channel stays.
        socket_.send_binary(binary_event({"type": "session_bundle", "requestId": request["requestId"]}, bundle))
        socket_.send_json({"type": "pong", "requestId": "after", "t": 1})
        _wait_until(lambda: viewer.state.last_event_type == "pong")
        assert channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token)
        result[0].path.unlink()


def test_an_invalid_session_bundle_fails_the_wait_at_once_and_leaves_no_file() -> None:
    with _channel_viewer() as (viewer, port):
        viewer._displayed = True
        viewer.state.ready = {"n_cells": 2, "dimensions": 2}
        socket_ = _connect(port, viewer)
        before = set(Path(tempfile.gettempdir()).glob("cellucid-session-*"))
        failure: list[BaseException] = []

        def capture() -> None:
            try:
                viewer.get_session_bundle(timeout=10)
            except BaseException as error:
                failure.append(error)

        worker = threading.Thread(target=capture)
        worker.start()
        request = socket_.receive_json()
        socket_.send_binary(
            binary_event(
                {"type": "session_bundle", "requestId": request["requestId"]},
                _corrupt_gzip_session_bundle(),
            )
        )
        worker.join(timeout=10)
        assert failure and "did not capture the session bundle" in str(failure[0])
        assert "invalid session bundle" in str(failure[0])
        assert set(Path(tempfile.gettempdir()).glob("cellucid-session-*")) == before
        assert channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token)


def test_a_session_bundle_cut_off_by_the_network_leaves_no_file() -> None:
    with _channel_viewer() as (viewer, port):
        viewer._displayed = True
        viewer.state.ready = {"n_cells": 2, "dimensions": 2}
        socket_ = _connect(port, viewer)
        before = set(Path(tempfile.gettempdir()).glob("cellucid-session-*"))
        failure: list[BaseException] = []

        def capture() -> None:
            try:
                viewer.get_session_bundle(timeout=None)
            except BaseException as error:
                failure.append(error)

        worker = threading.Thread(target=capture)
        worker.start()
        request = socket_.receive_json()
        frame = socket_.frame(
            BINARY,
            binary_event(
                {"type": "session_bundle", "requestId": request["requestId"]},
                SESSION_BUNDLE_MAGIC + b"\0" * 200_000,
            ),
        )
        socket_.sock.sendall(frame[: len(frame) // 2])
        # The bundle is being written to its file when the network drops.
        _wait_until(lambda: set(Path(tempfile.gettempdir()).glob("cellucid-session-*")) != before)
        socket_.sock.close()
        _wait_until(lambda: not channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token))
        _wait_until(lambda: set(Path(tempfile.gettempdir()).glob("cellucid-session-*")) == before)
        # The wait ends with the viewer, at once, never with a bundle.
        viewer.stop()
        worker.join(timeout=5)
        assert not worker.is_alive()
        assert len(failure) == 1 and isinstance(failure[0], RuntimeError)
        assert "the viewer was stopped" in str(failure[0])


def test_a_frame_sent_along_with_the_handshake_is_read() -> None:
    """Bytes the HTTP reader took in past the upgrade request are the channel's
    first, not lost when the channel reads the socket itself."""
    import base64
    import os
    import socket

    with _channel_viewer() as (viewer, port):
        key = base64.b64encode(os.urandom(16)).decode("ascii")
        request = (
            "GET /_cellucid/ws HTTP/1.1\r\n"
            f"Host: 127.0.0.1:{port}\r\n"
            "Upgrade: websocket\r\nConnection: Upgrade\r\n"
            f"Sec-WebSocket-Key: {key}\r\nSec-WebSocket-Version: 13\r\n"
            f"Origin: http://127.0.0.1:{port}\r\n\r\n"
        ).encode("ascii")
        hello = json.dumps({"type": "hello", "viewerId": viewer._viewer_id, "viewerToken": viewer._viewer_token})
        client = object.__new__(FakeViewer)
        client.sock = socket.create_connection(("127.0.0.1", port), timeout=5)
        client._buffer = b""
        client.sock.sendall(request + FakeViewer.frame(TEXT, hello.encode("utf-8")))
        assert client._read_until(b"\r\n\r\n").startswith(b"HTTP/1.1 101 ")
        assert client.receive_json() == {"type": "welcome"}
        client.close()


def test_a_viewer_that_reconnects_during_a_session_request_still_delivers_it() -> None:
    bundle = _empty_session_bundle()
    with _channel_viewer() as (viewer, port):
        viewer._displayed = True
        viewer.state.ready = {"n_cells": 2, "dimensions": 2}
        first = _connect(port, viewer)
        result: list[Any] = []
        worker = threading.Thread(target=lambda: result.append(viewer.get_session_bundle(timeout=10)))
        worker.start()
        request = first.receive_json()
        # The network drops before the bundle is sent; the page reconnects and
        # sends it on the new socket.
        first.sock.close()
        second = _connect(port, viewer)
        second.send_binary(binary_event({"type": "session_bundle", "requestId": request["requestId"]}, bundle))
        worker.join(timeout=10)
        assert result and result[0].path.read_bytes() == bundle
        result[0].path.unlink()


def test_stopping_the_viewer_freezes_an_open_page_and_closes_its_channel() -> None:
    viewer = serverless_viewer()
    viewer._activate()
    with _running_server() as (_host, port):
        socket_ = _connect(port, viewer)
        viewer.stop()
        assert socket_.receive_json() == {"type": "freeze"}
        assert socket_.receive_close() == (1001, "The viewer was stopped")
        # A stopped viewer's credentials no longer open a channel.
        late = FakeViewer(port)
        late.send_json({"type": "hello", "viewerId": viewer._viewer_id, "viewerToken": viewer._viewer_token})
        assert late.receive_close()[0] == 4401
    with pytest.raises(RuntimeError, match="stopped"):
        viewer.reset_view()


def test_the_server_pings_an_idle_channel(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(channel, "KEEPALIVE_SECONDS", 0.05)
    with _channel_viewer() as (viewer, port):
        socket_ = FakeViewer(port)
        socket_.hello(viewer._viewer_id, viewer._viewer_token)
        first, second = socket_._read_exact(2)
        assert first & 0x0F == 0x9


@pytest.mark.parametrize("offset", [0, 1, 2, 3, 5])
@pytest.mark.parametrize("length", [0, 1, 3, 4, 7, 64, 1001])
def test_unmasking_in_place_matches_the_rfc_rule(offset: int, length: int) -> None:
    rng = np.random.default_rng(length * 10 + offset)
    payload = rng.integers(0, 256, size=length + 1, dtype=np.uint8).tobytes()
    mask = bytes(rng.integers(0, 256, size=4, dtype=np.uint8))
    expected = bytes(byte ^ mask[(offset + index) % 4] for index, byte in enumerate(payload[1:]))
    # Slice from 1 so the buffer is misaligned as often as aligned.
    buffer = bytearray(payload)
    with memoryview(buffer) as whole, whole[1:] as view:
        unmask_in_place(view, mask, offset)
    assert bytes(buffer[1:]) == expected


def test_the_python_declarations_are_the_shared_viewer_protocol() -> None:
    fixture = json.loads((Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text())
    assert fixture["schema"] == 1
    assert fixture["channel"]["path"] == channel.CHANNEL_ROUTE == wire._CHANNEL_ROUTE
    assert fixture["channel"]["closeCodes"]["unauthorized"] == channel.CLOSE_UNAUTHORIZED
    assert fixture["channel"]["closeCodes"]["replaced"] == channel.CLOSE_REPLACED
    assert wire._protocol_capability_document() == fixture["protocolDocument"]
    assert {name: frozenset(keys) for name, keys in fixture["commands"].items()} == (
        wire._COMMAND_MESSAGE_FIELDS
    )
    assert {name: frozenset(keys) for name, keys in fixture["binaryCommands"].items()} == (
        wire._BINARY_COMMAND_HEADER_FIELDS
    )
    assert {name: frozenset(keys) for name, keys in fixture["binaryEvents"].items()} == (
        channel.BINARY_EVENT_HEADER_FIELDS
    )
    text_events = {
        name: keys
        for name, keys in wire._INBOUND_EVENT_FIELDS.items()
        if name not in channel.BINARY_EVENT_HEADER_FIELDS
    }
    assert {name: frozenset(keys) for name, keys in fixture["events"].items()} == text_events
    assert list(fixture["dataMappingOrigins"]) == sorted(
        wire._DATA_MAPPING_ORIGINS, key=fixture["dataMappingOrigins"].index
    )
    assert set(fixture["dataMappingOrigins"]) == set(wire._DATA_MAPPING_ORIGINS)
    # Every example command is one Python would send, byte for byte the same
    # once serialized; every example event is one the notebook accepts.
    for command in fixture["examples"]["commands"]:
        serialized = wire._require_frontend_message(command)
        assert isinstance(serialized, str)
        assert json.loads(serialized) == command
    # Every example binary command is the frame Python sends for it, byte for
    # byte, from its cells as a list and as a NumPy array in any order.
    binary_types = set()
    for example in fixture["examples"]["binaryCommands"]:
        header = example["header"]
        binary_types.add(header["type"])
        every = header.get("all") is True
        for cells in (example["cells"], np.array(example["cells"][::-1], dtype=np.int64)):
            message = {key: value for key, value in header.items() if key != "all"}
            message["cells"] = None if every else cells
            assert bytes(wire._require_frontend_message(message)).hex() == example["hex"]
    assert binary_types == set(wire._BINARY_COMMAND_HEADER_FIELDS)
    for event in fixture["examples"]["events"]:
        wire._require_inbound_jupyter_event(event)
    for example in fixture["examples"]["binaryEvents"]:
        frame = bytes.fromhex(example["hex"])
        header_length = int.from_bytes(frame[:4], "little")
        assert json.loads(frame[4 : 4 + header_length]) == example["header"]
        body = frame[4 + header_length :]
        if example["header"]["type"] == "selection":
            assert wire._decode_selection_cells(body) == example["cells"]
        else:
            assert body.decode() == example["bodyText"]


def test_the_shared_viewer_protocol_is_pinned() -> None:
    import hashlib

    data = (Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_bytes()
    assert hashlib.sha256(data).hexdigest() == VIEWER_PROTOCOL_SHA256
