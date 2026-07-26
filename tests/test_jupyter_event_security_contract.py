from __future__ import annotations

import http.client
import inspect
import json
import struct
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import cellucid._server_base as server_base
import cellucid.jupyter as jupyter
from cellucid.anndata_server import AnnDataServer

VIEWER_ID = "security-contract-viewer"
VIEWER_TOKEN = "security-contract-viewer-token"
REQUEST_ID = "security-contract-request"


def _known_inbound_events(viewer_id: str) -> list[tuple[str, dict[str, object], dict[str, object]]]:
    return [
        (
            "selection",
            {
                "type": "selection",
                "viewerId": viewer_id,
                "cells": [0, 2],
                "source": "lasso",
            },
            {"cells": [0, 2], "source": "lasso"},
        ),
        (
            "hover",
            {
                "type": "hover",
                "viewerId": viewer_id,
                "cell": 2,
                "position": {"x": 0.25, "y": -0.5, "z": 1.0},
            },
            {"cell": 2, "position": {"x": 0.25, "y": -0.5, "z": 1.0}},
        ),
        (
            "click",
            {
                "type": "click",
                "viewerId": viewer_id,
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
                "viewerId": viewer_id,
                "n_cells": 3,
                "dimensions": 2,
            },
            {"n_cells": 3, "dimensions": 2},
        ),
        (
            "pong",
            {
                "type": "pong",
                "viewerId": viewer_id,
                "requestId": "ping-request",
                "t": 1_722_000_000_000,
            },
            {"requestId": "ping-request", "t": 1_722_000_000_000},
        ),
        (
            "debug_snapshot",
            {
                "type": "debug_snapshot",
                "viewerId": viewer_id,
                "requestId": "debug-request",
                "ts": "2026-07-26T12:00:00.000Z",
                "locationHref": "https://viewer.example/?jupyter=true",
                "origin": "https://viewer.example",
                "serverUrl": "https://notebook.example/cellucid",
                "connected": True,
                "parentOrigin": "https://notebook.example",
                "userAgent": "Synthetic Browser",
            },
            {
                "requestId": "debug-request",
                "ts": "2026-07-26T12:00:00.000Z",
                "locationHref": "https://viewer.example/?jupyter=true",
                "origin": "https://viewer.example",
                "serverUrl": "https://notebook.example/cellucid",
                "connected": True,
                "parentOrigin": "https://notebook.example",
                "userAgent": "Synthetic Browser",
            },
        ),
        (
            "session_bundle",
            {
                "type": "session_bundle",
                "viewerId": viewer_id,
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
def _running_server() -> Iterator[tuple[str, int]]:
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
    host, port = server._server.server_address
    try:
        yield host, port
    finally:
        server.stop()


def _request(
    host: str,
    port: int,
    path: str,
    *,
    body: bytes | None = None,
    headers: dict[str, str] | None = None,
) -> tuple[int, dict[str, Any] | bytes]:
    connection = http.client.HTTPConnection(host, port, timeout=5)
    try:
        connection.request("POST", path, body=body, headers=headers or {})
        response = connection.getresponse()
        response_body = response.read()
        content_type = response.getheader("Content-Type")
        if content_type == "application/json":
            return response.status, json.loads(response_body)
        return response.status, response_body
    finally:
        connection.close()


def _event_body(*, viewer_token: str | None, event_type: str = "selection") -> bytes:
    event: dict[str, object] = {
        "type": event_type,
        "viewerId": VIEWER_ID,
        "cells": [0],
    }
    if viewer_token is not None:
        event["viewerToken"] = viewer_token
    return json.dumps(event).encode("utf-8")


def _event_request(
    host: str,
    port: int,
    *,
    viewer_token: str | None,
) -> tuple[int, dict[str, Any] | bytes]:
    return _request(
        host,
        port,
        "/_cellucid/events",
        body=_event_body(viewer_token=viewer_token),
        headers={"Content-Type": "application/json"},
    )


def _empty_session_bundle() -> bytes:
    manifest = {
        "createdAt": "2026-07-26T00:00:00.000Z",
        "dataSource": None,
        "datasetFingerprint": None,
        "summary": None,
        "chunks": [],
    }
    manifest_bytes = json.dumps(
        manifest,
        separators=(",", ":"),
    ).encode("utf-8")
    return (
        server_base.SESSION_BUNDLE_MAGIC + struct.pack("<I", len(manifest_bytes)) + manifest_bytes
    )


def _corrupt_gzip_session_bundle() -> bytes:
    stored = b"not-a-gzip-stream"
    manifest = {
        "createdAt": "2026-07-26T00:00:00.000Z",
        "dataSource": None,
        "datasetFingerprint": None,
        "summary": None,
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
        server_base.SESSION_BUNDLE_MAGIC
        + struct.pack("<I", len(manifest_bytes))
        + manifest_bytes
        + struct.pack("<I", len(stored))
        + stored
    )


def _session_path(*, viewer_token: str | None = VIEWER_TOKEN) -> str:
    pieces = [
        f"viewerId={VIEWER_ID}",
        f"requestId={REQUEST_ID}",
    ]
    if viewer_token is not None:
        pieces.append(f"viewerToken={viewer_token}")
    return "/_cellucid/session_bundle?" + "&".join(pieces)


def _session_request(
    host: str,
    port: int,
    *,
    viewer_token: str | None = VIEWER_TOKEN,
    body: bytes | None = None,
    headers: dict[str, str] | None = None,
) -> tuple[int, dict[str, Any] | bytes]:
    return _request(
        host,
        port,
        _session_path(viewer_token=viewer_token),
        body=_empty_session_bundle() if body is None else body,
        headers={"Content-Type": "application/octet-stream"} if headers is None else headers,
    )


@pytest.fixture(autouse=True)
def _restore_routing_registries() -> Iterator[None]:
    event_snapshot = dict(server_base._event_callbacks)
    pending_snapshot = dict(server_base._pending_session_bundle_requests)
    yield
    server_base._event_callbacks.clear()
    server_base._event_callbacks.update(event_snapshot)
    server_base._pending_session_bundle_requests.clear()
    server_base._pending_session_bundle_requests.update(pending_snapshot)


def test_routing_functions_expose_only_the_authenticated_current_contract() -> None:
    assert list(inspect.signature(server_base.register_event_callback).parameters) == [
        "viewer_id",
        "viewer_token",
        "callback",
    ]
    assert list(inspect.signature(server_base.unregister_event_callback).parameters) == [
        "viewer_id",
        "viewer_token",
    ]
    assert list(inspect.signature(server_base.route_event).parameters) == [
        "viewer_id",
        "viewer_token",
        "event",
    ]
    assert list(inspect.signature(server_base.register_session_bundle_request).parameters) == [
        "viewer_id",
        "viewer_token",
        "request_id",
        "ttl_seconds",
    ]
    assert list(inspect.signature(server_base.cancel_session_bundle_request).parameters) == [
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
    server_base.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        lambda _event: None,
    )
    before = dict(server_base._pending_session_bundle_requests)

    with pytest.raises((TypeError, ValueError), match="ttl_seconds"):
        server_base.register_session_bundle_request(
            VIEWER_ID,
            VIEWER_TOKEN,
            REQUEST_ID,
            ttl_seconds=ttl_seconds,  # type: ignore[arg-type]
        )

    assert server_base._pending_session_bundle_requests == before


def test_pending_request_registration_rejects_wrong_token_and_duplicate() -> None:
    server_base.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        lambda _event: None,
    )

    with pytest.raises(PermissionError):
        server_base.register_session_bundle_request(
            VIEWER_ID,
            "wrong-viewer-token",
            REQUEST_ID,
        )
    assert server_base._pending_session_bundle_requests == {}

    server_base.register_session_bundle_request(
        VIEWER_ID,
        VIEWER_TOKEN,
        REQUEST_ID,
    )
    pending = server_base._pending_session_bundle_requests[(VIEWER_ID, REQUEST_ID)]
    with pytest.raises(RuntimeError, match="already pending"):
        server_base.register_session_bundle_request(
            VIEWER_ID,
            VIEWER_TOKEN,
            REQUEST_ID,
        )
    assert server_base._pending_session_bundle_requests == {(VIEWER_ID, REQUEST_ID): pending}


@pytest.mark.parametrize("viewer_token", [None, "wrong-viewer-token"])
def test_event_post_rejects_absent_or_wrong_viewer_token(
    viewer_token: str | None,
) -> None:
    deliveries: list[dict[str, object]] = []
    server_base.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        deliveries.append,
    )

    with _running_server() as (host, port):
        status, response = _event_request(
            host,
            port,
            viewer_token=viewer_token,
        )

    assert status == 403
    assert response == b"Invalid viewer credentials"
    assert deliveries == []


def test_event_post_authenticates_with_compare_digest_and_delivers_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    deliveries: list[dict[str, object]] = []
    comparisons: list[tuple[str, str]] = []
    real_compare_digest = server_base.hmac.compare_digest

    def compare_digest(left: str, right: str) -> bool:
        comparisons.append((left, right))
        return real_compare_digest(left, right)

    monkeypatch.setattr(server_base.hmac, "compare_digest", compare_digest)
    server_base.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        deliveries.append,
    )

    with _running_server() as (host, port):
        status, response = _event_request(
            host,
            port,
            viewer_token=None,
        )
        assert status == 403
        assert response == b"Invalid viewer credentials"

        status, response = _event_request(
            host,
            port,
            viewer_token="wrong-viewer-token",
        )
        assert status == 403
        assert response == b"Invalid viewer credentials"

        status, response = _event_request(
            host,
            port,
            viewer_token=VIEWER_TOKEN,
        )

    assert status == 200
    assert response == {"status": "ok", "delivered": True}
    assert deliveries == [
        {
            "type": "selection",
            "viewerId": VIEWER_ID,
            "cells": [0],
        }
    ]
    assert comparisons == [
        (VIEWER_TOKEN, ""),
        (VIEWER_TOKEN, "wrong-viewer-token"),
        (VIEWER_TOKEN, VIEWER_TOKEN),
    ]


def test_event_post_reports_callback_absence_and_failure_without_false_success() -> None:
    calls = 0

    def fail(_event: dict[str, object]) -> None:
        nonlocal calls
        calls += 1
        raise RuntimeError("callback failed")

    with _running_server() as (host, port):
        status, response = _event_request(
            host,
            port,
            viewer_token=VIEWER_TOKEN,
        )
        assert status == 404
        assert response == b"Viewer is not registered"

        server_base.register_event_callback(
            VIEWER_ID,
            VIEWER_TOKEN,
            fail,
        )
        status, response = _event_request(
            host,
            port,
            viewer_token=VIEWER_TOKEN,
        )

    assert status == 500
    assert response == b"Viewer callback failed"
    assert calls == 1


def test_failed_viewer_hook_does_not_publish_partial_event_state() -> None:
    viewer = jupyter.BaseViewer()

    def fail(_event: dict[str, object]) -> None:
        raise RuntimeError("hook failed")

    viewer.register_hook("selection", fail)
    with pytest.raises(RuntimeError, match="hook failed"):
        viewer._handle_frontend_message(
            {
                "type": "selection",
                "viewerId": viewer._viewer_id,
                "cells": [0],
                "source": "lasso",
            }
        )

    assert viewer.state.last_event is None
    assert viewer.state.selection is None
    assert viewer._event_seq == 0
    assert list(viewer._recent_events) == []


def test_every_known_inbound_event_preserves_its_exact_payload() -> None:
    viewer = jupyter.BaseViewer()
    delivered: list[dict[str, object]] = []
    viewer.register_hook("message", delivered.append)

    known_events = _known_inbound_events(viewer._viewer_id)
    for event_name, message, expected_payload in known_events:
        viewer._handle_frontend_message(message)
        assert delivered[-1] == {"event": event_name, **expected_payload}
        assert viewer.state.last_event_type == event_name
        assert viewer.state.last_event == expected_payload

    assert len(delivered) == len(known_events)


@pytest.mark.parametrize("event_index", range(7))
def test_inbound_event_rejects_every_undeclared_property_before_callback_dispatch(
    event_index: int,
) -> None:
    viewer = jupyter.BaseViewer()
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
            "viewerId": "replace-at-runtime",
            "value": 1,
        },
        {
            "type": "selection",
            "viewerId": "replace-at-runtime",
            "cells": [0],
        },
        {
            "type": "selection",
            "viewerId": "replace-at-runtime",
            "viewerToken": "must-not-reach-the-viewer",
            "cells": [0],
            "source": "lasso",
        },
    ],
)
def test_inbound_event_rejects_unknown_types_missing_fields_and_routing_extras(
    message: dict[str, object],
) -> None:
    viewer = jupyter.BaseViewer()
    delivered: list[dict[str, object]] = []
    viewer.register_hook("message", delivered.append)
    candidate = dict(message)
    candidate["viewerId"] = viewer._viewer_id

    with pytest.raises((TypeError, ValueError), match="event|fields|unknown"):
        viewer._handle_frontend_message(candidate)

    assert delivered == []
    assert viewer.state.last_event is None


def test_event_http_path_rejects_undeclared_property_without_hook_delivery() -> None:
    viewer = jupyter.BaseViewer()
    delivered: list[dict[str, object]] = []
    viewer.register_hook("selection", delivered.append)
    viewer._activate()
    body = json.dumps(
        {
            "type": "selection",
            "viewerId": viewer._viewer_id,
            "viewerToken": viewer._viewer_token,
            "cells": [0],
            "source": "lasso",
            "undeclared": True,
        }
    ).encode("utf-8")

    try:
        with _running_server() as (host, port):
            status, response = _request(
                host,
                port,
                "/_cellucid/events",
                body=body,
                headers={"Content-Type": "application/json"},
            )
        assert status == 500
        assert response == b"Viewer callback failed"
        assert delivered == []
        assert viewer.state.last_event is None
    finally:
        viewer.stop()


@pytest.mark.parametrize(
    ("body", "headers", "expected_status"),
    [
        (b"", {"Content-Type": "application/octet-stream"}, 411),
        (
            _empty_session_bundle(),
            {"Content-Type": "text/plain"},
            415,
        ),
        (
            server_base.SESSION_BUNDLE_MAGIC,
            {"Content-Type": "application/octet-stream"},
            400,
        ),
        (
            b"not-a-session-bundle",
            {"Content-Type": "application/octet-stream"},
            400,
        ),
        (
            _corrupt_gzip_session_bundle(),
            {"Content-Type": "application/octet-stream"},
            400,
        ),
    ],
)
def test_invalid_session_body_does_not_consume_the_pending_request(
    body: bytes,
    headers: dict[str, str],
    expected_status: int,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    deliveries: list[dict[str, object]] = []
    real_named_temporary_file = tempfile.NamedTemporaryFile

    def named_temporary_file(*args: Any, **kwargs: Any):
        kwargs["dir"] = tmp_path
        return real_named_temporary_file(*args, **kwargs)

    monkeypatch.setattr(
        server_base.tempfile,
        "NamedTemporaryFile",
        named_temporary_file,
    )
    server_base.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        deliveries.append,
    )
    server_base.register_session_bundle_request(
        VIEWER_ID,
        VIEWER_TOKEN,
        REQUEST_ID,
        ttl_seconds=60.0,
    )

    with _running_server() as (host, port):
        status, _response = _session_request(
            host,
            port,
            body=body,
            headers=headers,
        )
        assert status == expected_status
        assert deliveries == []
        assert list(tmp_path.iterdir()) == []

        status, response = _session_request(host, port)

    assert status == 200
    assert response == {"status": "ok", "bytes": len(_empty_session_bundle())}
    assert len(deliveries) == 1
    delivered_path = Path(str(deliveries[0]["path"]))
    assert delivered_path.is_file()
    delivered_path.unlink()


@pytest.mark.parametrize("viewer_token", [None, "wrong-viewer-token"])
def test_session_upload_is_bound_to_the_pending_viewer_token(
    viewer_token: str | None,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    deliveries: list[dict[str, object]] = []
    real_named_temporary_file = tempfile.NamedTemporaryFile

    def named_temporary_file(*args: Any, **kwargs: Any):
        kwargs["dir"] = tmp_path
        return real_named_temporary_file(*args, **kwargs)

    monkeypatch.setattr(
        server_base.tempfile,
        "NamedTemporaryFile",
        named_temporary_file,
    )
    server_base.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        deliveries.append,
    )
    server_base.register_session_bundle_request(
        VIEWER_ID,
        VIEWER_TOKEN,
        REQUEST_ID,
        ttl_seconds=60.0,
    )

    with _running_server() as (host, port):
        status, response = _session_request(
            host,
            port,
            viewer_token=viewer_token,
        )
        assert status == 403
        assert response == b"Invalid viewer credentials"
        assert deliveries == []
        assert list(tmp_path.iterdir()) == []

        status, _response = _session_request(host, port)

    assert status == 200
    delivered_path = Path(str(deliveries[0]["path"]))
    delivered_path.unlink()


def test_session_upload_delivers_once_and_consumes_only_the_valid_request(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    deliveries: list[dict[str, object]] = []
    real_named_temporary_file = tempfile.NamedTemporaryFile

    def named_temporary_file(*args: Any, **kwargs: Any):
        kwargs["dir"] = tmp_path
        return real_named_temporary_file(*args, **kwargs)

    monkeypatch.setattr(
        server_base.tempfile,
        "NamedTemporaryFile",
        named_temporary_file,
    )
    server_base.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        deliveries.append,
    )
    server_base.register_session_bundle_request(
        VIEWER_ID,
        VIEWER_TOKEN,
        REQUEST_ID,
        ttl_seconds=60.0,
    )

    with _running_server() as (host, port):
        status, response = _session_request(host, port)
        assert status == 200
        assert response == {"status": "ok", "bytes": len(_empty_session_bundle())}

        status, response = _session_request(host, port)
        assert status == 404
        assert response == b"No pending session bundle request"

    assert len(deliveries) == 1
    assert deliveries[0]["status"] == "ok"
    assert deliveries[0]["requestId"] == REQUEST_ID
    assert deliveries[0]["bytes"] == len(_empty_session_bundle())
    delivered_path = Path(str(deliveries[0]["path"]))
    assert delivered_path.read_bytes() == _empty_session_bundle()
    delivered_path.unlink()


@pytest.mark.parametrize("callback_state", ["absent", "failure"])
def test_session_callback_absence_or_failure_is_non_success_and_removes_temp_file(
    callback_state: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = 0
    real_named_temporary_file = tempfile.NamedTemporaryFile

    def named_temporary_file(*args: Any, **kwargs: Any):
        kwargs["dir"] = tmp_path
        return real_named_temporary_file(*args, **kwargs)

    def callback(_event: dict[str, object]) -> None:
        nonlocal calls
        calls += 1
        if callback_state == "failure":
            raise RuntimeError("session callback failed")

    monkeypatch.setattr(
        server_base.tempfile,
        "NamedTemporaryFile",
        named_temporary_file,
    )
    server_base.register_event_callback(
        VIEWER_ID,
        VIEWER_TOKEN,
        callback,
    )
    server_base.register_session_bundle_request(
        VIEWER_ID,
        VIEWER_TOKEN,
        REQUEST_ID,
        ttl_seconds=60.0,
    )
    if callback_state == "absent":
        del server_base._event_callbacks[VIEWER_ID]

    with _running_server() as (host, port):
        status, _response = _session_request(host, port)
        repeat_status, repeat_response = _session_request(host, port)

    assert status == (410 if callback_state == "absent" else 500)
    assert repeat_status == 404
    assert repeat_response == b"No pending session bundle request"
    assert calls == (0 if callback_state == "absent" else 1)
    assert list(tmp_path.iterdir()) == []


def _routing_snapshot() -> tuple[dict[str, object], set[object]]:
    return dict(server_base._event_callbacks), set(jupyter._active_viewers)


def test_missing_prepared_data_does_not_change_viewer_registries(
    tmp_path: Path,
) -> None:
    before = _routing_snapshot()

    with pytest.raises(FileNotFoundError, match="Data directory"):
        jupyter.CellucidViewer(
            tmp_path / "missing",
            auto_open=False,
        )

    assert _routing_snapshot() == before


@pytest.mark.parametrize("viewer_type", ["prepared", "anndata"])
def test_server_start_failure_rolls_back_viewer_registration_and_server_resources(
    viewer_type: str,
    tmp_path: Path,
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
    if viewer_type == "prepared":
        monkeypatch.setattr(jupyter, "CellucidServer", FailingServer)

        def constructor() -> object:
            return jupyter.CellucidViewer(
                tmp_path,
                auto_open=False,
            )

    else:
        import cellucid.anndata_server as anndata_server

        monkeypatch.setattr(anndata_server, "AnnDataServer", FailingServer)

        def constructor() -> object:
            return jupyter.AnnDataViewer(
                _minimal_adata(),
                auto_open=False,
                dataset_name="Failed viewer",
                dataset_id="failed-viewer",
            )

    with pytest.raises(RuntimeError, match="server start failed"):
        constructor()

    assert len(instances) == 1
    assert instances[0].stop_calls == 1
    assert _routing_snapshot() == before


def test_display_failure_rolls_back_started_server_and_viewer_registries(
    tmp_path: Path,
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

    monkeypatch.setattr(jupyter, "CellucidServer", StartedServer)
    monkeypatch.setattr(
        jupyter,
        "_detect_jupyter_context",
        lambda: {"in_jupyter": True},
    )
    monkeypatch.setattr(jupyter.CellucidViewer, "display", fail_display)
    before = _routing_snapshot()

    with pytest.raises(RuntimeError, match="display failed"):
        jupyter.CellucidViewer(
            tmp_path,
            auto_open=True,
        )

    assert len(instances) == 1
    assert instances[0].stop_calls == 1
    assert _routing_snapshot() == before
