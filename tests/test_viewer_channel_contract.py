from __future__ import annotations

import inspect
import json
import re
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
from _stand_in_server import StandInServer
from _viewer_channel_client import (
    BINARY,
    TEXT,
    FakeViewer,
    HandshakeRefused,
    binary_event,
)
from _web_root import assert_shared_file

import cellucid._viewer_channel as channel
import cellucid._websocket as websocket_module
import cellucid.jupyter as jupyter
import cellucid.jupyter._viewer as jupyter_viewer
import cellucid.jupyter._wire as wire
from cellucid._server_base import CellucidThreadingHTTPServer
from cellucid._websocket import unmask_in_place
from cellucid.anndata_server import AnnDataServer
from cellucid.session import SESSION_SIGNATURE

VIEWER_ID = "security-contract-viewer"
VIEWER_TOKEN = "security-contract-viewer-token"
REQUEST_ID = "security-contract-request"

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
                "n_obs": 5,
                "dimensions": 2,
            },
            {"n_cells": 3, "n_obs": 5, "dimensions": 2},
        ),
        (
            # Every cell visible: no cells, no values, an empty body.
            "result",
            {
                "type": "result",
                "epoch": "0123456789abcdef",
                "revision": 4,
                "kind": "visibility",
                "id": "visible",
                "nObs": 5,
                "meta": {},
                "cells": None,
                "values": None,
                "body": bytearray(),
            },
            {
                "epoch": "0123456789abcdef",
                "revision": 4,
                "kind": "visibility",
                "id": "visible",
                "nObs": 5,
                "meta": {},
                "cells": None,
                "values": None,
            },
        ),
        (
            # An action's answer: its outcome is the binary frame's body, JSON.
            "action_result",
            {
                "type": "action_result",
                "requestId": "action-request",
                "step": {"id": "s3", "label": "Colored by leiden"},
                "unrecorded": None,
                "body": bytearray(b'{"groups":2}'),
            },
            {
                "requestId": "action-request",
                "step": {"id": "s3", "label": "Colored by leiden"},
                "unrecorded": None,
                "outcome": {"groups": 2},
            },
        ),
        (
            "action_result",
            {"type": "action_result", "requestId": "read-request", "step": None, "unrecorded": None, "body": bytearray(b"null")},
            {"requestId": "read-request", "step": None, "unrecorded": None, "outcome": None},
        ),
        (
            # A change that made no step, and why: it was so already.
            "action_result",
            {"type": "action_result", "requestId": "same-request", "step": None, "unrecorded": "unchanged", "body": bytearray(b"null")},
            {"requestId": "same-request", "step": None, "unrecorded": "unchanged", "outcome": None},
        ),
        (
            "result_removed",
            {"type": "result_removed", "epoch": "0123456789abcdef", "revision": 5, "kind": "highlight", "id": "highlight_1"},
            {"epoch": "0123456789abcdef", "revision": 5, "kind": "highlight", "id": "highlight_1"},
        ),
        (
            "results",
            {
                "type": "results",
                "epoch": "0123456789abcdef",
                "dataset": "demo",
                "revision": 5,
                "nObs": 5,
                "items": [{"kind": "visibility", "id": "visible", "revision": 4}],
            },
            {
                "epoch": "0123456789abcdef",
                "dataset": "demo",
                "revision": 5,
                "nObs": 5,
                "items": [{"kind": "visibility", "id": "visible", "revision": 4}],
            },
        ),
        (
            "color_by",
            {"type": "color_by", "field": "leiden", "source": "obs", "symbol": None, "cells_read": None},
            {"field": "leiden", "source": "obs", "symbol": None, "cells_read": None},
        ),
        (
            "color_by",
            # A gene read for the first 40,000 open cells only (a gene sample).
            {"type": "color_by", "field": "ENSG00000168484", "source": "var", "symbol": "SFTPC", "cells_read": 40_000},
            {"field": "ENSG00000168484", "source": "var", "symbol": "SFTPC", "cells_read": 40_000},
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
                "command": "action",
                "requestId": "action-request",
                "reason": "There is no field \"CD8A\" in obs",
            },
            {"command": "action", "requestId": "action-request", "reason": "There is no field \"CD8A\" in obs"},
        ),
        (
            "data_mapping",
            {
                "type": "data_mapping",
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
        (
            "history_changed",
            {"type": "history_changed", "kind": "step", "tracker": "t-1", "head": "s4", "step": "s4", "branch": "b0"},
            {"kind": "step", "tracker": "t-1", "head": "s4", "step": "s4", "branch": "b0"},
        ),
        (
            "session_loaded",
            {"type": "session_loaded", "requestId": "load-request", "tracker": {"id": "t-1", "name": "Shared"}, "head": "s9"},
            {"requestId": "load-request", "tracker": {"id": "t-1", "name": "Shared"}, "head": "s9"},
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


#: A session file of format 4: the shared golden history.
_SESSION_FILE = Path(__file__).parent / "fixtures" / "history_golden" / "history_golden.cellucid-session"


def _session_file() -> bytes:
    return _SESSION_FILE.read_bytes()


def _damaged_session_file() -> bytes:
    """The fixture with one byte of its last chunk flipped: framed as a
    session, and refused when its chunks are read."""
    data = bytearray(_session_file())
    data[-3] ^= 0xFF
    return bytes(data)


def _deeply_nested_session_file(depth: int = 5_000) -> bytes:
    """The fixture with one eager value replaced by lists nested ``depth``
    deep: JSON, nested past what the reader follows (Python's recursion, in
    its JSON parser or in the reader's own walk, wherever it ends first; no
    Cellucid writes one). The text is written as text, never parsed here."""
    from _session_files import frame, payload, store, unframe

    manifest, chunks = unframe(_session_file())
    index = next(i for i, meta in enumerate(manifest["chunks"]) if meta["kind"] == "values")
    meta = manifest["chunks"][index]
    values = json.loads(payload(meta, chunks[index]))
    values[next(iter(values))] = "DEEP"
    text = json.dumps(values).replace('"DEEP"', "[" * depth + "]" * depth, 1)
    kept = {key: value for key, value in meta.items() if key not in ("codec", "storedBytes", "uncompressedBytes")}
    manifest["chunks"][index], chunks[index] = store(kept, text.encode("utf-8"))
    return frame(manifest, chunks)


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
        # Hears the page holding the channel let it go (the agent calls it
        # was running then fail at once).
        "on_page_gone",
        # The cells of the dataset the viewer shows: the bound of every cell
        # and binary body its events name.
        "cell_count",
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
    """``request.result(timeout=None)`` waits forever, so its request must too."""
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
    # A failing hook must not drop the event from viewer.state and
    # wait_for_event (a session request would time out behind it): the event
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


@pytest.mark.core
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
        "command": "action",
        "requestId": "action-request",
        "reason": 'There is no field "CD8A" in obs',
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
        "command": "action",
        "requestId": "action-request",
        "reason": 'There is no field "CD8A" in obs',
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
            "command": "action",
            "requestId": "action-request",
            "reason": 'There is no field "CD8A" in obs',
        }
    ]


def test_unhandled_command_error_is_printed_where_the_notebook_looks(
    capsys: pytest.CaptureFixture[str],
) -> None:
    viewer = serverless_viewer()
    observed: list[dict[str, object]] = []
    viewer.register_hook("message", observed.append)

    # An action this viewer sent is named by its own method; one it did not
    # send (or no longer awaits) by the generic call.
    request = viewer._expect_action("highlight.cells")
    viewer._handle_frontend_message(
        _command_error(viewer._viewer_id, requestId=request.request_id, reason="There is no page \"Notebook\"")
    )
    viewer._handle_frontend_message(_command_error(viewer._viewer_id))
    viewer._handle_frontend_message(
        _command_error(viewer._viewer_id, command="requestSessionBundle", requestId="r", reason="no dataset")
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    # Named by the Python call that sent it, not by the wire type, and repeated
    # for every failure rather than deduplicated.
    assert captured.err == (
        "[cellucid] viewer.highlight_cells(...) did not take effect in the viewer: "
        "There is no page \"Notebook\"\n"
        "[cellucid] The viewer is unchanged. Handle this in code with "
        "@viewer.on_command_error.\n"
        "[cellucid] viewer.act(...) did not take effect in the viewer: "
        'There is no field "CD8A" in obs\n'
        "[cellucid] The viewer is unchanged. Handle this in code with "
        "@viewer.on_command_error.\n"
        "[cellucid] viewer.request_session() did not take effect in the viewer: no dataset\n"
        "[cellucid] The viewer is unchanged. Handle this in code with "
        "@viewer.on_command_error.\n"
    )
    assert [entry["event"] for entry in observed] == [
        "command_error",
        "command_error",
        "command_error",
    ]
    # The failed request is answered, and forgotten by the viewer.
    assert request.done and viewer._action_waiters == {} and viewer._action_types == {}


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"command": "notACommand"}, "one current Jupyter command"),
        ({"command": "command_error"}, "one current Jupyter command"),
        ({"command": ""}, "command_error command"),
        ({"command": None}, "command_error command"),
        ({"command": "action "}, "command_error command"),
        # A state change is an action, never a command of its own.
        ({"command": "setColorBy"}, "one current Jupyter command"),
        ({"requestId": ""}, "command_error requestId"),
        ({"requestId": 7}, "command_error requestId"),
        ({"reason": ""}, "command_error reason"),
        ({"reason": None}, "command_error reason"),
        ({"reason": "line one\nline two"}, "command_error reason"),
        # This text is printed to a terminal, so escape sequences must not
        # survive the wire.
        ({"reason": "\x1b[31mred"}, "command_error reason"),
        ({"reason": "carriage\rreturn"}, "command_error reason"),
        ({"reason": "x" * 4001}, "at most 4000 characters"),
        ({"extra": "unexpected"}, "unknown extra"),
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
        _command_error(viewer._viewer_id, reason="x" * 4000)
    )

    assert handled == [{"command": "action", "requestId": "action-request", "reason": "x" * 4000}]


def test_the_published_events_are_exactly_the_events_a_viewer_may_send() -> None:
    """``/_cellucid/protocol`` publishes what the validator accepts, no more.

    ``_known_inbound_events`` is the independent list here: every entry in it is
    driven through ``_handle_frontend_message`` by
    ``test_every_known_inbound_event_preserves_its_exact_payload``, so tying the
    published set to it means the route can only name types this installation
    demonstrably accepts, and can only omit one by leaving a type untested.
    """
    published = wire._protocol_capability_document()["events"]

    assert set(published) == {name for name, _message, _payload in _known_inbound_events(VIEWER_ID)} | set(
        _known_agent_events()
    )


def _known_agent_events() -> dict[str, dict[str, object]]:
    """One of each agent event (an agent app's answers, ``cellucid mcp``)."""
    return {
        "agent_tools": {"type": "agent_tools", "requestId": "a1", "body": bytearray(b'[{"name": "color_by"}]')},
        "agent_result": {"type": "agent_result", "requestId": "a2", "body": bytearray(b'{"content": [], "isError": false}')},
        "agent_asking": {"type": "agent_asking", "requestId": "a2", "message": "Claude Code (local) asks to color.", "detail": None},
        "agent_answered": {"type": "agent_answered", "requestId": "a2", "allowed": True},
    }


def test_every_agent_event_reaches_the_agent_port_and_never_the_notebook() -> None:
    viewer = serverless_viewer()
    delivered: list[dict[str, object]] = []
    viewer.register_hook("message", delivered.append)
    port: list[tuple[str, dict[str, object]]] = []
    viewer._agent_port.deliver = lambda event, data: port.append((event, data))  # type: ignore[method-assign]
    for message in _known_agent_events().values():
        viewer._handle_frontend_message(dict(message))
    assert [event for event, _data in port] == list(_known_agent_events())
    assert port[0][1] == {"requestId": "a1", "tools": [{"name": "color_by"}]}
    assert port[1][1] == {"requestId": "a2", "result": {"content": [], "isError": False}}
    assert delivered == []
    assert viewer.state.last_event is None
    for message in _known_agent_events().values():
        with pytest.raises(ValueError, match="unknown.*undeclared"):
            viewer._handle_frontend_message({**message, "undeclared": True})


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

    failed: list[dict[str, object]] = []
    viewer._agent_port.command_failed = failed.append  # type: ignore[method-assign]
    for command in published:
        viewer._handle_frontend_message(
            _command_error(viewer._viewer_id, command=command)
        )

    # An agent app's command fails to the agent that sent it (`cellucid mcp`), never to a notebook hook.
    assert [entry["command"] for entry in handled] == [command for command in published if not command.startswith("agent_")]
    assert [entry["command"] for entry in failed] == [command for command in published if command.startswith("agent_")]

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

    class StartedServer(StandInServer):
        def __init__(self, **kwargs: object):
            super().__init__(**{**kwargs, "port": 43123})
            instances.append(self)

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


def test_a_socket_that_never_says_hello_is_closed_and_holds_no_thread(monkeypatch: pytest.MonkeyPatch) -> None:
    # cellucid-r closes it after .viewer_hello_seconds; before, Python held
    # the connection's thread for as long as the peer kept it open.
    monkeypatch.setattr(channel, "HELLO_SECONDS", 0.3)
    monkeypatch.setattr(websocket_module, "CLOSE_DRAIN_SECONDS", 0.2)
    with _viewer_on_a_server() as (_viewer, _server, http_server, port):
        socket_ = FakeViewer(port)
        assert socket_.receive_close() == (1008, "No hello frame was received")
        _wait_until(lambda: _open_connections(http_server) == 0)


def test_a_binary_first_frame_is_refused_from_its_header_unread(monkeypatch: pytest.MonkeyPatch) -> None:
    # A frame declaring a gibibyte, of which a few bytes come: refused at
    # once, before the socket said who it is, never read to its end.
    monkeypatch.setattr(websocket_module, "CLOSE_DRAIN_SECONDS", 0.2)
    with _viewer_on_a_server() as (_viewer, _server, http_server, port):
        socket_ = FakeViewer(port)
        socket_.sock.sendall(struct.pack("!BBQ", 0x80 | BINARY, 0x80 | 127, 1 << 30) + b"\0" * 4 + b"x" * 64)
        assert socket_.receive_close() == (1008, "The first frame must be the text frame 'hello'")
        _wait_until(lambda: _open_connections(http_server) == 0)


def test_a_hello_nested_past_the_interpreter_closes_the_channel_as_a_protocol_violation() -> None:
    # Within the text frame limit, and deeper than Python's JSON parser
    # follows: refused as malformed (before: RecursionError, a traceback on
    # the kernel's stderr for each such socket, and a 1001 close).
    with _channel_viewer() as (_viewer, port):
        socket_ = FakeViewer(port)
        socket_.send_frame(TEXT, b"[" * 200_000)
        code, reason = socket_.receive_close()
        assert (code, reason) == (1008, "The hello frame is nested more deeply than this package reads")


def test_a_binary_header_length_allocates_only_what_arrives() -> None:
    # Four bytes naming a 4 GiB header, then a few bytes and the frame's end:
    # the buffer is never allocated (and zeroed) at the size named before a
    # byte of it arrives.
    import tracemalloc

    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        tracemalloc.start()
        try:
            socket_.send_binary(struct.pack("<I", 0xFFFF_FFFF) + b'{"type"')
            assert socket_.receive_close() == (1008, "Binary message ended early")
            _current, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        assert peak < 16 * 1024 * 1024


def test_a_body_frame_length_allocates_only_what_arrives() -> None:
    # A frame header naming 256 MiB, of which an event header and a few body
    # bytes arrive before the page stops: the body is never allocated (and
    # zeroed) at the frame's length before a byte of it comes.
    import socket as socket_module
    import tracemalloc

    claimed = 256 * 1024 * 1024
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        payload = binary_event({"type": "agent_result", "requestId": "r1"}, b'{"content":')
        tracemalloc.start()
        try:
            socket_.sock.sendall(struct.pack("!BBQ", 0x80 | BINARY, 0x80 | 127, claimed) + b"\0\0\0\0" + payload)
            time.sleep(0.3)
            socket_.sock.shutdown(socket_module.SHUT_WR)
            _wait_until(lambda: not channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token))
            _current, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        assert peak < 16 * 1024 * 1024


def _ready_viewer_socket(viewer: jupyter.AnnDataViewer, port: int, n_obs: int) -> FakeViewer:
    socket_ = _connect(port, viewer)
    socket_.send_json({"type": "ready", "n_cells": n_obs, "n_obs": n_obs, "dimensions": 2})
    _wait_until(lambda: viewer.state.ready is not None and viewer.state.ready["n_obs"] == n_obs)
    return socket_


def test_cells_and_bodies_past_the_dataset_are_refused_before_they_are_read() -> None:
    # What a page names is bounded by the dataset the viewer shows
    # (cellucid-r's .check_event_cells): a selection of more cells than it
    # has, a cell beyond it, a result for a dataset of another size, a hover
    # beyond it. Before, each was delivered (a selection's obs names then
    # failed with an IndexError in the notebook) or read whole.
    result_header = {
        "type": "result", "epoch": "0" * 16, "revision": 1, "kind": "highlight", "id": "g1", "nObs": 5,
        "meta": {}, "cells": 1, "values": None,
    }
    cases = [
        (binary_event({"type": "selection", "source": "lasso"}, struct.pack("<4I", 0, 1, 2, 3)),
         (1009, "Binary message longer than the 12 bytes it can hold")),
        (binary_event({"type": "selection", "source": "lasso"}, struct.pack("<2I", 0, 7)),
         (1008, "A selection names cell 7, and the dataset has 3")),
        (binary_event(result_header, struct.pack("<I", 0)),
         (1008, "A result names a dataset of 5 cells, and this one has 3")),
        (binary_event({**result_header, "nObs": 3, "cells": 9}, struct.pack("<I", 0)),
         (1008, "A result's cells must be a count of at most its nObs, or null")),
    ]
    with _channel_viewer() as (viewer, port):
        for payload, close in cases:
            socket_ = _ready_viewer_socket(viewer, port, 3)
            socket_.send_binary(payload)
            assert socket_.receive_close() == close
        socket_ = _ready_viewer_socket(viewer, port, 3)
        socket_.send_json({"type": "hover", "cell": 9, "position": None})
        assert socket_.receive_close() == (1008, "A hover names cell 9, and the dataset has 3")
        assert viewer.state.selection is None
        # Within the dataset, the same events arrive.
        socket_ = _ready_viewer_socket(viewer, port, 3)
        socket_.send_binary(binary_event({"type": "selection", "source": "lasso"}, struct.pack("<3I", 0, 1, 2)))
        _wait_until(lambda: viewer.state.selection is not None)
        assert viewer.state.selection == {"cells": [0, 1, 2], "source": "lasso"}


def test_the_result_body_bound_names_the_widest_value_the_wire_carries() -> None:
    assert max(dtype.itemsize for dtype in wire._VALUE_DTYPES.values()) == channel._WIDEST_VALUE_BYTES
    assert wire._SELECTION_CELL_DTYPE.itemsize == channel._CELL_INDEX_BYTES


def test_an_event_nested_past_the_interpreter_is_refused_as_malformed() -> None:
    deep: Any = 0
    for _ in range(5_000):
        deep = [deep]
    with pytest.raises(ValueError, match="nested more deeply than this package reads"):
        wire._require_inbound_jupyter_event({"type": "click", "cell": deep, "button": 0, "shift": False, "ctrl": False})
    body = ("[" * 5_000 + "]" * 5_000).encode("utf-8")
    with pytest.raises(ValueError, match="nested more deeply|not UTF-8 JSON"):
        wire._require_inbound_jupyter_event({"type": "agent_result", "requestId": "a", "body": body})


@pytest.mark.core
def test_commands_sent_before_the_viewer_connects_arrive_in_order_once() -> None:
    with _channel_viewer() as (viewer, port):
        assert viewer.send_message({"type": "freeze"}) is False
        colored = viewer.color_by(field={"source": "obs", "key": "leiden"})
        highlighted = viewer.highlight_cells(cells=[2, 0], color="#00ff00")
        hidden = viewer.cells_visibility_set(cells=None, visible=False)
        socket_ = _connect(port, viewer)
        assert socket_.receive_json() == {"type": "freeze"}
        assert socket_.receive_json() == {
            "type": "action", "requestId": colored.request_id, "action": "color.by",
            "params": {"field": {"source": "obs", "key": "leiden"}}, "turn": None, "turnLabel": None,
        }
        # An action whose parameters hold cells is a binary frame: its cells
        # as uint32, in the order given.
        assert socket_.receive_command() == {
            "type": "action", "requestId": highlighted.request_id, "action": "highlight.cells",
            "params": {"color": "#00ff00", "tool": "annotation"}, "turn": None, "turnLabel": None,
            "arrays": [{"path": ["cells"], "type": "uint32", "length": 2, "elements": [2, 0]}],
        }
        assert socket_.receive_json() == {
            "type": "action", "requestId": hidden.request_id, "action": "cells.visibility.set",
            "params": {"cells": None, "visible": False}, "turn": None, "turnLabel": None,
        }
        cleared = viewer.highlight_clear_all()
        assert socket_.receive_json() == {
            "type": "action", "requestId": cleared.request_id, "action": "highlight.clear-all", "params": {},
            "turn": None, "turnLabel": None,
        }
        # A reconnect gets only what was sent while it was away.
        socket_.sock.close()
        _wait_until(lambda: not channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token))
        reset = viewer.camera_reset(view="live")
        again = _connect(port, viewer)
        assert again.receive_json()["requestId"] == reset.request_id
        assert viewer.send_message({"type": "freeze"}) is True
        assert again.receive_json() == {"type": "freeze"}


def test_a_newer_socket_of_the_same_viewer_replaces_the_open_one() -> None:
    with _channel_viewer() as (viewer, port):
        first = _connect(port, viewer)
        second = _connect(port, viewer)
        assert first.receive_close() == (4409, "Replaced by a newer connection of this viewer")
        assert viewer.send_message({"type": "freeze"}) is True
        assert second.receive_json() == {"type": "freeze"}


@pytest.mark.core
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
        socket_.send_json({"type": "ready", "n_cells": 2, "n_obs": 2, "dimensions": 2})
        assert viewer.wait_for_ready(timeout=5) == {"n_cells": 2, "n_obs": 2, "dimensions": 2}
        assert [event["cell"] for event in hovers] == list(range(50))
        assert clicks == [{"cell": 7, "button": 0, "shift": False, "ctrl": True}]


def test_a_fragmented_event_with_a_ping_between_its_frames_is_one_event() -> None:
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        payload = json.dumps({"type": "ready", "n_cells": 2, "n_obs": 2, "dimensions": 3}).encode()
        socket_.sock.sendall(
            socket_.frame(TEXT, payload[:7], fin=False)
            + socket_.frame(0x9, b"are you there")
            + socket_.frame(0x0, payload[7:15], fin=False)
            + socket_.frame(0x0, payload[15:])
        )
        opcode, pong = socket_.receive()
        assert (opcode, pong) == (0xA, b"are you there")
        assert viewer.wait_for_ready(timeout=5) == {"n_cells": 2, "n_obs": 2, "dimensions": 3}


@pytest.mark.parametrize(
    ("frames", "code"),
    [
        # An event the viewer may only send as a binary frame.
        ([(TEXT, json.dumps({"type": "selection", "cells": [0], "source": "lasso"}).encode())], 1008),
        # An event this installation does not declare.
        ([(TEXT, json.dumps({"type": "zoom", "level": 2}).encode())], 1008),
        # An event that names its viewer: events carry no credentials.
        ([(TEXT, json.dumps({"type": "ready", "viewerId": "x", "n_cells": 1, "n_obs": 1, "dimensions": 2}).encode())], 1008),
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
        reset = viewer.camera_reset(view="live")
        again = _connect(port, viewer)
        assert again.receive_json()["requestId"] == reset.request_id


def test_a_failing_notebook_hook_keeps_the_channel_open() -> None:
    with _channel_viewer() as (viewer, port):
        def fail(_event: dict[str, Any]) -> None:
            raise RuntimeError("the notebook's hook failed")

        viewer.register_hook("click", fail)
        socket_ = _connect(port, viewer)
        socket_.send_json({"type": "click", "cell": 1, "button": 0, "shift": False, "ctrl": False})
        socket_.send_json({"type": "ready", "n_cells": 2, "n_obs": 2, "dimensions": 2})
        assert viewer.wait_for_ready(timeout=5) == {"n_cells": 2, "n_obs": 2, "dimensions": 2}
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


def test_a_session_request_streams_the_frame_to_a_file_once() -> None:
    bundle = _session_file()
    with _channel_viewer() as (viewer, port):
        viewer._displayed = True
        viewer.state.ready = {"n_cells": 2, "n_obs": 2, "dimensions": 2}
        socket_ = _connect(port, viewer)
        result: list[Any] = []
        worker = threading.Thread(target=lambda: result.append(viewer.request_session().result(timeout=10)))
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
        viewer.state.ready = {"n_cells": 2, "n_obs": 2, "dimensions": 2}
        socket_ = _connect(port, viewer)
        before = set(Path(tempfile.gettempdir()).glob("cellucid-session-*"))
        failure: list[BaseException] = []

        def capture() -> None:
            try:
                viewer.request_session().result(timeout=10)
            except BaseException as error:
                failure.append(error)

        worker = threading.Thread(target=capture)
        worker.start()
        request = socket_.receive_json()
        socket_.send_binary(
            binary_event(
                {"type": "session_bundle", "requestId": request["requestId"]},
                _damaged_session_file(),
            )
        )
        worker.join(timeout=10)
        assert failure and "did not capture the session" in str(failure[0])
        assert "invalid session bundle" in str(failure[0])
        assert set(Path(tempfile.gettempdir()).glob("cellucid-session-*")) == before
        assert channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token)


def test_a_session_nested_past_the_reader_is_refused_in_a_sentence_and_leaves_no_file() -> None:
    # Before: the interpreter's RecursionError escaped the channel's thread,
    # the person's session stayed in the temporary folder, and the wait ran
    # to its timeout.
    with _channel_viewer() as (viewer, port):
        viewer._displayed = True
        viewer.state.ready = {"n_cells": 2, "n_obs": 2, "dimensions": 2}
        socket_ = _connect(port, viewer)
        before = set(Path(tempfile.gettempdir()).glob("cellucid-session-*"))
        failure: list[BaseException] = []

        def capture() -> None:
            try:
                viewer.request_session().result(timeout=10)
            except BaseException as error:
                failure.append(error)

        worker = threading.Thread(target=capture)
        worker.start()
        request = socket_.receive_json()
        started = time.monotonic()
        socket_.send_binary(
            binary_event(
                {"type": "session_bundle", "requestId": request["requestId"]},
                _deeply_nested_session_file(),
            )
        )
        worker.join(timeout=10)
        assert time.monotonic() - started < 5
        assert failure and "nested more deeply than this package reads" in str(failure[0])
        assert set(Path(tempfile.gettempdir()).glob("cellucid-session-*")) == before
        assert channel.viewer_is_connected(viewer._viewer_id, viewer._viewer_token)


def test_a_session_bundle_whose_read_fails_unexpectedly_leaves_no_file(monkeypatch: pytest.MonkeyPatch) -> None:
    # The written part holds the person's session: whatever ends its read
    # (before: anything but the expected errors left it in the temp folder).
    import cellucid.session as session_module

    def broken(_path: Path) -> None:
        raise RuntimeError("a defect while reading")

    monkeypatch.setattr(session_module, "CellucidSession", broken)
    with _channel_viewer() as (viewer, port):
        viewer._displayed = True
        viewer.state.ready = {"n_cells": 2, "n_obs": 2, "dimensions": 2}
        socket_ = _connect(port, viewer)
        before = set(Path(tempfile.gettempdir()).glob("cellucid-session-*"))
        # The channel ended with the defect: the wait runs to its timeout.
        worker = threading.Thread(target=lambda: pytest.raises(TimeoutError, viewer.request_session().result, timeout=3))
        worker.start()
        request = socket_.receive_json()
        socket_.send_binary(binary_event({"type": "session_bundle", "requestId": request["requestId"]}, _session_file()))
        assert socket_.receive_close()[0] == 1001
        worker.join(timeout=10)
        assert set(Path(tempfile.gettempdir()).glob("cellucid-session-*")) == before


def test_a_session_file_nested_past_the_reader_is_refused_in_a_sentence(tmp_path: Path) -> None:
    from cellucid.session import TOO_DEEP_SESSION_MESSAGE, CellucidSession

    path = tmp_path / "deep.cellucid-session"
    path.write_bytes(_deeply_nested_session_file())
    with pytest.raises(ValueError) as refusal:
        CellucidSession(path)
    assert str(refusal.value) == TOO_DEEP_SESSION_MESSAGE


def test_a_session_bundle_cut_off_by_the_network_leaves_no_file() -> None:
    with _channel_viewer() as (viewer, port):
        viewer._displayed = True
        viewer.state.ready = {"n_cells": 2, "n_obs": 2, "dimensions": 2}
        socket_ = _connect(port, viewer)
        before = set(Path(tempfile.gettempdir()).glob("cellucid-session-*"))
        failure: list[BaseException] = []

        def capture() -> None:
            try:
                viewer.request_session().result(timeout=None)
            except BaseException as error:
                failure.append(error)

        worker = threading.Thread(target=capture)
        worker.start()
        request = socket_.receive_json()
        frame = socket_.frame(
            BINARY,
            binary_event(
                {"type": "session_bundle", "requestId": request["requestId"]},
                SESSION_SIGNATURE + b"\0" * 200_000,
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
    bundle = _session_file()
    with _channel_viewer() as (viewer, port):
        viewer._displayed = True
        viewer.state.ready = {"n_cells": 2, "n_obs": 2, "dimensions": 2}
        first = _connect(port, viewer)
        result: list[Any] = []
        worker = threading.Thread(target=lambda: result.append(viewer.request_session().result(timeout=10)))
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
        viewer.camera_reset(view="live")


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


def _example_params(example: dict[str, Any]) -> dict[str, Any]:
    """A binary action example's parameters with each array in place (as a
    NumPy array of its type, appended in order), as a notebook gives them."""
    import copy

    params = copy.deepcopy(example["header"]["params"])
    for entry in example["arrays"]:
        node: Any = params
        for key in entry["path"][:-1]:
            node = node[key]
        node[entry["path"][-1]] = np.array(entry["elements"], dtype=wire._VALUE_DTYPES[entry["type"]])
    return params


@pytest.mark.core
def test_the_python_declarations_are_the_shared_viewer_protocol() -> None:
    fixture = json.loads((Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text(encoding="utf-8"))
    assert fixture["schema"] == 5
    assert fixture["channel"]["path"] == channel.CHANNEL_ROUTE == wire._CHANNEL_ROUTE
    assert fixture["channel"]["closeCodes"]["unauthorized"] == channel.CLOSE_UNAUTHORIZED
    assert fixture["channel"]["closeCodes"]["replaced"] == channel.CLOSE_REPLACED
    assert wire._protocol_capability_document() == fixture["protocolDocument"]
    assert {kind: frozenset(keys) for kind, keys in fixture["results"]["kinds"].items()} == wire._RESULT_META_FIELDS
    assert set(fixture["results"]["highlightTools"]) == wire._HIGHLIGHT_TOOLS
    assert set(fixture["results"]["origins"]) == wire._RESULT_ORIGINS
    assert {name: wire._VALUE_DTYPES[name].itemsize for name in fixture["results"]["valueTypes"]} == (
        fixture["results"]["valueTypes"]
    )
    assert {name: frozenset(keys) for name, keys in fixture["commands"].items()} == (
        wire._COMMAND_MESSAGE_FIELDS
    )
    assert {name: frozenset(keys) for name, keys in fixture["binaryCommands"].items()} == (
        wire._BINARY_COMMAND_HEADER_FIELDS
    )
    assert tuple(fixture["action"]["arrayFields"]) == wire._ACTION_ARRAY_FIELDS
    assert frozenset(fixture["action"]["stepFields"]) == wire._ACTION_STEP_FIELDS
    assert fixture["action"]["typePattern"] == wire._ACTION_TYPE_PATTERN.pattern
    assert tuple(fixture["history"]["changeKinds"]) == wire._HISTORY_CHANGE_KINDS
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
    assert fixture["commandErrorReasonLimit"] < wire._COMMAND_ERROR_REASON_LIMIT
    # Every example command is one Python would send: an action without
    # arrays as JSON text.
    for command in fixture["examples"]["commands"]:
        serialized = wire._require_frontend_message(command)
        assert isinstance(serialized, str)
        assert json.loads(serialized) == command
    # Every example binary action is the frame Python sends for it, byte for
    # byte, from its parameters with NumPy arrays in place.
    for example in fixture["examples"]["binaryCommands"]:
        header = example["header"]
        if header["type"] == "load_session":
            message = {
                "type": "load_session", "requestId": header["requestId"], "name": header["name"],
                "file": example["bodyText"].encode("utf-8"),
            }
        else:
            message = {
                "type": "action",
                "requestId": header["requestId"],
                "action": header["action"],
                "params": _example_params(example),
                "turn": header["turn"],
                "turnLabel": header["turnLabel"],
            }
        assert bytes(wire._require_frontend_message(message)).hex() == example["hex"]
    for event in fixture["examples"]["events"]:
        wire._require_inbound_jupyter_event(event)
    binary_types = set()
    for example in fixture["examples"]["binaryEvents"]:
        frame = bytes.fromhex(example["hex"])
        header_length = int.from_bytes(frame[:4], "little")
        assert json.loads(frame[4 : 4 + header_length]) == example["header"]
        body = frame[4 + header_length :]
        binary_types.add(example["header"]["type"])
        if example["header"]["type"] == "selection":
            assert wire._decode_selection_cells(body) == example["cells"]
        elif example["header"]["type"] == "result":
            _kind, payload = wire._require_inbound_jupyter_event({**example["header"], "body": bytearray(body)})
            assert (None if payload["cells"] is None else payload["cells"].tolist()) == example["cells"]
            assert (None if payload["values"] is None else payload["values"].tolist()) == example["values"]
        elif example["header"]["type"] == "action_result":
            _kind, payload = wire._require_inbound_jupyter_event({**example["header"], "body": bytearray(body)})
            assert payload == {
                "requestId": example["header"]["requestId"],
                "step": example["header"]["step"],
                "unrecorded": example["header"]["unrecorded"],
                "outcome": example["outcome"],
            }
        elif example["header"]["type"] in {"agent_tools", "agent_result"}:
            # An agent app's tools, or a call's MCP result: the body is UTF-8 JSON.
            _kind, payload = wire._require_inbound_jupyter_event({**example["header"], "body": bytearray(body)})
            key = "tools" if example["header"]["type"] == "agent_tools" else "result"
            assert payload == {"requestId": example["header"]["requestId"], key: example["body"]}
        else:
            assert body.decode() == example["bodyText"]
    assert binary_types == set(channel.BINARY_EVENT_HEADER_FIELDS)


@pytest.mark.parametrize(
    ("params", "match"),
    [
        ({"list": [np.array([1], dtype=np.uint32)]}, "inside a list"),
        ({"cells": np.array([1], dtype=np.int32)}, "uint8, uint16, uint32, float32 or float64"),
        ({"cells": np.zeros((2, 2), dtype=np.uint32)}, "one-dimensional"),
        ({"n": float("nan")}, "finite"),
        ({1: "x"}, "keys must be native strings"),
    ],
)
def test_an_action_the_viewer_cannot_read_is_refused_before_it_is_sent(params: dict[Any, Any], match: str) -> None:
    with pytest.raises((TypeError, ValueError), match=match):
        wire._require_frontend_message(
            {"type": "action", "requestId": "r", "action": "x.y", "params": params, "turn": None, "turnLabel": None}
        )
    for action in ("color", "Color.by", "color.by.", ""):
        with pytest.raises(ValueError, match="dotted action type"):
            wire._require_frontend_message(
                {"type": "action", "requestId": "r", "action": action, "params": {}, "turn": None, "turnLabel": None}
            )
    with pytest.raises(ValueError, match="names its turn and the turn's label together"):
        wire._require_frontend_message(
            {"type": "action", "requestId": "r", "action": "x.y", "params": {}, "turn": "t", "turnLabel": None}
        )
    with pytest.raises(ValueError, match="without a folder"):
        wire._require_frontend_message({"type": "load_session", "requestId": "r", "name": "a/b.cellucid-session", "file": b"x"})
    with pytest.raises(TypeError, match="session file's bytes"):
        wire._require_frontend_message({"type": "load_session", "requestId": "r", "name": "b.cellucid-session", "file": b""})
    for retired in ("setColorBy", "highlight", "setField", "setDataMapping", "resetCamera"):
        with pytest.raises(ValueError, match="Unknown viewer command type"):
            wire._require_frontend_message({"type": retired})


def test_the_shared_viewer_protocol_is_the_web_apps_byte_for_byte() -> None:
    assert_shared_file("tests/fixtures/viewer_protocol.json", "tests/fixtures/viewer-protocol.json")


# =============================================================================
# The history: notebook turns, sessions pushed in, history events
# =============================================================================


def test_actions_inside_a_turn_carry_it_and_turns_do_not_nest() -> None:
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        with viewer.turn("  Normalize   and cluster ") as turn_id:
            assert re.fullmatch(r"python-[0-9a-f]{12}", turn_id)
            colored = viewer.color_by(field={"source": "obs", "key": "leiden"})
            highlighted = viewer.highlight_cells(cells=[1], color="#00ff00")
            nested = pytest.raises(RuntimeError, match='A turn \\("Normalize and cluster"\\) is open already; turns do not nest.')
            with nested, viewer.turn("Inner"):
                pass
        after = viewer.history_undo()
        first = socket_.receive_command()
        assert (first["requestId"], first["turn"], first["turnLabel"]) == (colored.request_id, turn_id, "Normalize and cluster")
        second = socket_.receive_command()
        assert (second["requestId"], second["turn"], second["turnLabel"]) == (highlighted.request_id, turn_id, "Normalize and cluster")
        third = socket_.receive_command()
        assert (third["action"], third["turn"], third["turnLabel"]) == ("history.undo", None, None)
        assert third["requestId"] == after.request_id
        with pytest.raises(ValueError, match="A turn needs a label"), viewer.turn("   "):
            pass


def test_a_session_pushed_in_is_one_binary_frame_answered_by_session_loaded() -> None:
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        request = viewer.load_session(_SESSION_FILE)
        command = socket_.receive_command()
        assert command == {"type": "load_session", "requestId": request.request_id, "name": _SESSION_FILE.name, "file": _session_file()}
        assert not request.done
        socket_.send_json({
            "type": "session_loaded", "requestId": request.request_id,
            "tracker": {"id": "t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10", "name": "Annotation"}, "head": "s13",
        })
        assert request.result(timeout=5) == {"tracker": {"id": "t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10", "name": "Annotation"}, "head": "s13"}
        refused = viewer.load_session(_SESSION_FILE)
        socket_.receive_command()
        socket_.send_json({
            "type": "command_error", "command": "load_session", "requestId": refused.request_id,
            "reason": "This session was saved on another dataset.",
        })
        with pytest.raises(RuntimeError, match=f"The viewer did not load {_SESSION_FILE.name}: This session was saved on another dataset."):
            refused.result(timeout=5)
        with pytest.raises(FileNotFoundError):
            viewer.load_session(_SESSION_FILE.with_name("absent.cellucid-session"))


def test_every_history_change_is_kept_and_reaches_the_hooks() -> None:
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        seen: list[dict[str, Any]] = []
        viewer.history.on_change(seen.append)
        assert viewer.history.last_change is None and repr(viewer.history) == "LiveHistory(no change told yet)"
        change = {"type": "history_changed", "kind": "step", "tracker": "t-1", "head": "s4", "step": "s4", "branch": "b0"}
        socket_.send_json(change)
        moved = {**change, "kind": "moved", "head": "s2", "step": None, "branch": None}
        socket_.send_json(moved)
        _wait_until(lambda: len(seen) == 2)
        assert seen == [{key: value for key, value in event.items() if key != "type"} for event in (change, moved)]
        history = viewer.history
        assert (history.tracker, history.head, history.branch) == ("t-1", "s2", None)
        assert viewer.state.history == seen[-1]
        # A kind the protocol does not have breaks it.
        socket_.send_json({**change, "kind": "renamed"})
        code, _reason = socket_.receive_close()
        assert code == 1008


def test_the_historys_actions_are_the_generated_methods() -> None:
    """The history is changed through the generated history_* methods: the
    live history only tells what changed, it does not wrap the actions."""
    with _channel_viewer() as (viewer, port):
        socket_ = _connect(port, viewer)
        assert not any(callable(getattr(viewer.history, name, None)) for name in ("undo", "jump", "note", "bookmark"))
        requests = [
            (viewer.history_jump(step="s3"), "history.jump", {"step": "s3"}),
            (viewer.history_note_add(target={"step": "s3"}, text="Looks right."), "history.note.add",
             {"target": {"step": "s3"}, "text": "Looks right.", "from_reply": False}),
            (viewer.history_turn_undo(turn="python-1", only_own=True), "history.turn.undo", {"turn": "python-1", "only_own": True}),
        ]
        for request, action, params in requests:
            command = socket_.receive_command()
            assert (command["requestId"], command["action"], command["params"]) == (request.request_id, action, params)
