"""
The viewer channel: one authenticated WebSocket per notebook viewer.

A viewer page opens ``<serverUrl>/_cellucid/ws`` and says ``hello`` with the
viewer id and token from its address. The server checks them against the
registration :class:`cellucid.jupyter.AnnDataViewer` made, answers
``welcome``, and from then on the socket is that viewer's:

- **commands** (Python → viewer) are text frames, and a command that names
  cells (a highlight, a visibility change) a binary frame. A command sent
  while no socket is open waits in the viewer's outbox and is sent, in order,
  the moment one is (:func:`send_viewer_command`); nothing is sent twice.
- **events** (viewer → Python) are text frames, and a selection, a result or a
  session bundle a binary frame (``uint32`` little-endian header length, JSON header,
  body). Every event is delivered to the viewer's callback on the socket's own
  thread, one at a time and in order.

The same channel also runs over a notebook widget's comm
(:mod:`cellucid.jupyter._widget`): the widget relays each frame between the
page it embeds and the kernel, and :func:`serve_viewer_channel` reads it
through the same :class:`ChannelSocket` interface a WebSocket offers, so
authentication, ordering, replacement, binary events and session bundles are
one implementation for both.

A newer ``hello`` for the same viewer (its page reloaded, or it reconnected
after the network dropped) replaces the open socket, which is closed with
``4409``. Stopping the viewer, or the server, closes its socket with ``1001``.
Every close is the orderly one of :meth:`cellucid._websocket.WebSocket.close`:
the page reads the code and reason it was sent, never a reset.

The rules are the cross-language ones in ``viewer-protocol.json`` (the web
app's ``data/viewer-protocol.js``); :mod:`cellucid.jupyter._wire` checks each
message's fields.

This module also keeps the pending ``requestSessionBundle`` requests: a
session bundle frame is accepted only for a request Python made and has not
cancelled.
"""

from __future__ import annotations

import collections
import hmac
import io
import ipaddress
import json
import logging
import math
import tempfile
import threading
import time
import zlib
from collections.abc import Callable
from http import HTTPStatus
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol, cast
from urllib.parse import urlsplit

from ._server_base import CellucidThreadingHTTPServer
from ._websocket import (
    CLOSE_GOING_AWAY,
    CLOSE_POLICY_VIOLATION,
    WebSocket,
    WebSocketClosed,
    WebSocketProtocolError,
    accept_key,
    require_client_key,
)
from .session_bundle import NOT_A_SESSION_FILE_MESSAGE, SESSION_BUNDLE_MAGIC

if TYPE_CHECKING:
    import socket
    from http.client import HTTPMessage
    from http.server import HTTPServer
    from typing import BinaryIO

logger = logging.getLogger("cellucid.viewer_channel")

#: Where the channel is, relative to the server's base address.
CHANNEL_ROUTE = "_cellucid/ws"

#: Seconds between the server's pings on an idle socket. A proxy between the
#: notebook page and this server (jupyter-server-proxy behind nginx, a
#: JupyterHub ingress) closes a connection that carries nothing for a while;
#: a ping every twenty seconds keeps it open, and the browser answers it.
KEEPALIVE_SECONDS = 20.0

#: Longest text frame a viewer sends. Every text event is a small record
#: (hover, click, a data mapping); the large ones are binary frames.
MAX_EVENT_TEXT_BYTES = 1024 * 1024

CLOSE_UNAUTHORIZED = 4401
CLOSE_REPLACED = 4409

#: Event types that arrive only as binary frames.
BINARY_EVENT_HEADER_FIELDS: dict[str, frozenset[str]] = {
    "selection": frozenset({"type", "source"}),
    "session_bundle": frozenset({"type", "requestId"}),
    "result": frozenset({"type", "requestId", "epoch", "revision", "kind", "id", "nObs", "meta", "cells", "values"}),
}

_HELLO_FIELDS = frozenset({"type", "viewerId", "viewerToken"})
_WELCOME = json.dumps({"type": "welcome"}, separators=(",", ":"))

# Read size for a session bundle streamed to disk.
_BUNDLE_CHUNK_BYTES = 1024 * 1024


class ChannelMessage(Protocol):
    """One binary frame being read in order, as a stream: a WebSocket message
    (:class:`cellucid._websocket.BinaryMessage`) or a widget's frame."""

    done: bool

    def readinto(self, buffer: memoryview) -> int: ...
    def read_exact(self, size: int) -> bytearray: ...
    def read_rest(self, initial: int = 0) -> bytearray: ...
    def drain(self) -> None: ...


class ChannelSocket(Protocol):
    """What one viewer channel runs over: the viewer's WebSocket
    (:class:`cellucid._websocket.WebSocket`), or a notebook widget's comm.

    :meth:`receive` blocks for the next frame and raises
    :class:`cellucid._websocket.WebSocketClosed` once the socket closed; the
    send and close methods may be called from any thread.
    """

    _closed: threading.Event

    @property
    def closed(self) -> bool: ...

    def receive(self) -> str | ChannelMessage: ...
    def send_text(self, text: str) -> None: ...
    def send_binary(self, payload: bytes | bytearray) -> None: ...
    def ping(self) -> None: ...
    def close(self, code: int, reason: str = "") -> None: ...
    def finish(self) -> None: ...


class _ViewerNotRegisteredError(LookupError):
    """The supplied viewer id has no live registration."""


class _ViewerAuthenticationError(PermissionError):
    """The supplied viewer token does not authenticate the viewer."""


class _SessionBundleRequestNotFoundError(LookupError):
    """The supplied session request is absent or expired."""


class ViewerProtocolError(ValueError):
    """A viewer broke the channel's rules; its socket is closed with 1008.

    A delivery callback raises it (from a field check) to tell the channel the
    event itself was wrong, as opposed to a notebook hook that failed.
    """


def _require_routing_string(value: object, *, label: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string")
    if not value or value != value.strip() or any(character.isspace() for character in value):
        raise ValueError(f"{label} must be a non-empty string without whitespace")
    return value


# =============================================================================
# Registrations
# =============================================================================


class _Registration:
    __slots__ = (
        "callback",
        "channel_lock",
        "delivery_lock",
        "outbox",
        "socket",
        "viewer_token",
    )

    def __init__(self, viewer_token: str, callback: Callable[[dict], None]) -> None:
        self.viewer_token = viewer_token
        self.callback = callback
        self.delivery_lock = threading.RLock()
        # Guards `socket` and `outbox` together, so a command is either sent
        # on the open socket or queued, and a new socket is handed every
        # queued command before any later one.
        self.channel_lock = threading.Lock()
        self.socket: ChannelSocket | None = None
        self.outbox: collections.deque[str | bytearray] = collections.deque()


_DUMMY_VIEWER_TOKEN = "\0" * 32
_event_callbacks_lock = threading.Lock()
_event_callbacks: dict[str, _Registration] = {}


def _require_event_registration_locked(viewer_id: str, viewer_token: str) -> _Registration:
    registration = _event_callbacks.get(viewer_id)
    expected_token = registration.viewer_token if registration is not None else _DUMMY_VIEWER_TOKEN
    authenticated = hmac.compare_digest(expected_token, viewer_token)
    if registration is None:
        raise _ViewerNotRegisteredError
    if not authenticated:
        raise _ViewerAuthenticationError
    return registration


def _require_registration(viewer_id: str, viewer_token: str) -> _Registration:
    with _event_callbacks_lock:
        return _require_event_registration_locked(viewer_id, viewer_token)


def register_event_callback(
    viewer_id: str,
    viewer_token: str,
    callback: Callable[[dict], None],
) -> None:
    """Register one token-authenticated viewer and the callback its events reach."""
    _require_routing_string(viewer_id, label="viewer_id")
    _require_routing_string(viewer_token, label="viewer_token")
    if not callable(callback):
        raise TypeError("callback must be callable")
    registration = _Registration(viewer_token, callback)
    with _event_callbacks_lock:
        if viewer_id in _event_callbacks:
            raise RuntimeError("An event callback is already registered for this viewer_id")
        _event_callbacks[viewer_id] = registration


def unregister_event_callback(viewer_id: str, viewer_token: str) -> bool:
    """Unregister only the viewer authenticated by the token, and close its socket."""
    _require_routing_string(viewer_id, label="viewer_id")
    _require_routing_string(viewer_token, label="viewer_token")
    with _event_callbacks_lock:
        registration = _event_callbacks.get(viewer_id)
        if registration is None:
            hmac.compare_digest(_DUMMY_VIEWER_TOKEN, viewer_token)
            return False
        if not hmac.compare_digest(registration.viewer_token, viewer_token):
            raise _ViewerAuthenticationError

    with registration.delivery_lock:
        with _event_callbacks_lock:
            if _event_callbacks.get(viewer_id) is not registration:
                return False
            del _event_callbacks[viewer_id]

        with registration.channel_lock:
            socket_ = registration.socket
            registration.socket = None
            registration.outbox.clear()
        if socket_ is not None:
            socket_.close(CLOSE_GOING_AWAY, "The viewer was stopped")

        with _pending_session_bundle_lock:
            matching_requests = [
                key
                for key, request in _pending_session_bundle_requests.items()
                if key[0] == viewer_id and hmac.compare_digest(request.viewer_token, viewer_token)
            ]
            for key in matching_requests:
                del _pending_session_bundle_requests[key]
        return True


def route_event(viewer_id: str, viewer_token: str, event: dict) -> None:
    """Authenticate and deliver one event exactly once."""
    _require_routing_string(viewer_id, label="viewer_id")
    if not isinstance(viewer_token, str):
        raise TypeError("viewer_token must be a string")
    if not isinstance(event, dict):
        raise TypeError("event must be a dictionary")

    registration = _require_registration(viewer_id, viewer_token)
    with registration.delivery_lock:
        with _event_callbacks_lock:
            if _event_callbacks.get(viewer_id) is not registration:
                raise _ViewerNotRegisteredError
        registration.callback(event)


def _send_command(websocket: ChannelSocket, message: str | bytearray) -> None:
    if isinstance(message, str):
        websocket.send_text(message)
    else:
        websocket.send_binary(message)


def send_viewer_command(
    viewer_id: str,
    viewer_token: str,
    message: str | bytearray,
    *,
    queue: bool = True,
) -> bool:
    """Send one serialized command to the viewer: JSON text, or the bytes of a
    binary command frame (:func:`cellucid.jupyter._wire._require_frontend_message`).

    Returns True when it was written to the viewer's open socket. Without one,
    the command is queued for the next socket the viewer opens and False is
    returned; with ``queue=False`` it is dropped instead (a ``freeze`` sent
    while stopping is meant for an open page only).
    """
    if type(message) not in (str, bytearray):
        raise TypeError("message must be serialized text or the bytes of a binary command frame")
    registration = _require_registration(viewer_id, viewer_token)
    with registration.channel_lock:
        socket_ = registration.socket
        if socket_ is not None and not socket_.closed:
            try:
                _send_command(socket_, message)
                return True
            except WebSocketClosed:
                # The frame did not go out whole, so the viewer never read it:
                # queueing it for the next socket cannot deliver it twice.
                registration.socket = None
        if queue:
            registration.outbox.append(message)
        return False


def viewer_is_connected(viewer_id: str, viewer_token: str) -> bool:
    """Whether the viewer has an open, authenticated socket now."""
    registration = _require_registration(viewer_id, viewer_token)
    with registration.channel_lock:
        return registration.socket is not None and not registration.socket.closed


# =============================================================================
# Pending session bundle requests
# =============================================================================


class _PendingSessionBundleRequest:
    __slots__ = ("expires_at", "viewer_token")

    def __init__(self, viewer_token: str, expires_at: float):
        self.viewer_token = viewer_token
        self.expires_at = expires_at


_pending_session_bundle_lock = threading.Lock()
_pending_session_bundle_requests: dict[tuple[str, str], _PendingSessionBundleRequest] = {}


def register_session_bundle_request(
    viewer_id: str,
    viewer_token: str,
    request_id: str,
    ttl_seconds: float | None = 60.0,
) -> None:
    """Register one authenticated, one-shot session bundle request.

    ``ttl_seconds`` is how long the request stays answerable; None keeps it
    until it is answered or cancelled, for a caller that waits without a
    deadline (``request.result(timeout=None)``).
    """
    _require_routing_string(viewer_id, label="viewer_id")
    _require_routing_string(viewer_token, label="viewer_token")
    _require_routing_string(request_id, label="request_id")
    if ttl_seconds is not None:
        if isinstance(ttl_seconds, bool) or not isinstance(ttl_seconds, int | float):
            raise TypeError("ttl_seconds must be a number or None")
        if not math.isfinite(ttl_seconds) or ttl_seconds <= 0:
            raise ValueError("ttl_seconds must be finite and greater than zero")

    now = time.monotonic()
    expires_at = math.inf if ttl_seconds is None else now + ttl_seconds
    key = (viewer_id, request_id)
    with _event_callbacks_lock:
        _require_event_registration_locked(viewer_id, viewer_token)
        with _pending_session_bundle_lock:
            expired = [
                pending_key
                for pending_key, request in _pending_session_bundle_requests.items()
                if request.expires_at <= now
            ]
            for pending_key in expired:
                del _pending_session_bundle_requests[pending_key]
            if key in _pending_session_bundle_requests:
                raise RuntimeError(
                    "A session bundle request with this viewer_id and request_id is already pending"
                )
            _pending_session_bundle_requests[key] = _PendingSessionBundleRequest(
                viewer_token=viewer_token,
                expires_at=expires_at,
            )


def _session_bundle_request(
    viewer_id: str,
    viewer_token: str,
    request_id: str,
    *,
    consume: bool,
) -> None:
    """Authenticate one live pending request, consuming it when asked."""
    now = time.monotonic()
    key = (viewer_id, request_id)
    with _pending_session_bundle_lock:
        request = _pending_session_bundle_requests.get(key)
        expected_token = request.viewer_token if request is not None else _DUMMY_VIEWER_TOKEN
        authenticated = hmac.compare_digest(expected_token, viewer_token)
        if request is None:
            raise _SessionBundleRequestNotFoundError
        if request.expires_at <= now:
            del _pending_session_bundle_requests[key]
            raise _SessionBundleRequestNotFoundError
        if not authenticated:
            raise _ViewerAuthenticationError
        if consume:
            del _pending_session_bundle_requests[key]


def cancel_session_bundle_request(viewer_id: str, viewer_token: str, request_id: str) -> bool:
    """Cancel one pending request without accepting alternate credentials."""
    _require_routing_string(viewer_id, label="viewer_id")
    _require_routing_string(viewer_token, label="viewer_token")
    _require_routing_string(request_id, label="request_id")
    now = time.monotonic()
    key = (viewer_id, request_id)
    with _pending_session_bundle_lock:
        request = _pending_session_bundle_requests.get(key)
        if request is None:
            return False
        if request.expires_at <= now:
            del _pending_session_bundle_requests[key]
            return False
        if not hmac.compare_digest(request.viewer_token, viewer_token):
            raise _ViewerAuthenticationError
        del _pending_session_bundle_requests[key]
        return True


# =============================================================================
# The WebSocket upgrade
# =============================================================================


class _UpgradeHandler(Protocol):
    """The parts of the HTTP handler the upgrade reads and writes."""

    def _request_headers(self) -> HTTPMessage: ...
    def _request_reader(self) -> BinaryIO: ...
    def _response_writer(self) -> BinaryIO: ...
    def _allowed_host_names(self) -> frozenset[str]: ...
    def _request_server(self) -> HTTPServer: ...
    def send_error_response(self, code: int, message: str) -> None: ...

    close_connection: bool
    connection: socket.socket
    command: str


def _single_header(headers: HTTPMessage, name: str) -> str | None:
    values = headers.get_all(name, [])
    return values[0] if len(values) == 1 else None


def _host_name(authority: str) -> str | None:
    """The lower-case host of ``host[:port]``, brackets removed from IPv6."""
    try:
        parsed = urlsplit(f"//{authority}")
        host = parsed.hostname
        _ = parsed.port
    except ValueError:
        return None
    return host.lower() if host else None


def _is_loopback_name(host: str) -> bool:
    if host == "localhost" or host.endswith(".localhost"):
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


def origin_is_allowed(origin: str | None, host: str | None, allowed_hosts: frozenset[str]) -> bool:
    """Whether a page at ``origin`` may open the channel.

    A browser sends ``Origin`` with every WebSocket handshake and no CORS
    check applies to WebSockets, so this is the cross-site check: the page
    must be this server's own (same host as the request's ``Host``, which
    the ``Host`` check already confined to this server's names), on a
    loopback name, or on a name the operator declared (``allowed_hosts``,
    the notebook proxy's).
    """
    if origin is None:
        return False
    parsed = urlsplit(origin)
    if parsed.scheme not in ("http", "https") or not parsed.hostname:
        return False
    origin_host = parsed.hostname.lower()
    if _is_loopback_name(origin_host) or origin_host in allowed_hosts:
        return True
    return host is not None and _host_name(host) == origin_host


def accept_viewer_channel(handler: _UpgradeHandler) -> None:
    """Upgrade ``GET /_cellucid/ws`` to the viewer channel and serve it to its end."""
    handler.close_connection = True
    headers = handler._request_headers()
    upgrade = (_single_header(headers, "Upgrade") or "").lower()
    connection = (_single_header(headers, "Connection") or "").lower()
    if handler.command != "GET" or upgrade != "websocket" or "upgrade" not in {
        token.strip() for token in connection.split(",")
    }:
        handler.send_error_response(
            HTTPStatus.UPGRADE_REQUIRED,
            "The viewer channel is a WebSocket: open it with new WebSocket(...)",
        )
        return
    if _single_header(headers, "Sec-WebSocket-Version") != "13":
        handler.send_error_response(HTTPStatus.UPGRADE_REQUIRED, "WebSocket version 13 is required")
        return
    try:
        client_key = require_client_key(_single_header(headers, "Sec-WebSocket-Key"))
    except ValueError as error:
        handler.send_error_response(HTTPStatus.BAD_REQUEST, str(error))
        return
    if not origin_is_allowed(
        _single_header(headers, "Origin"),
        _single_header(headers, "Host"),
        handler._allowed_host_names(),
    ):
        logger.warning("Refused a viewer channel from a page this server does not serve")
        handler.send_error_response(
            HTTPStatus.FORBIDDEN, "The viewer channel accepts only this server's own pages"
        )
        return

    # Written whole and by hand: the handler's own header path adds CORS and
    # caching headers, which mean nothing on a protocol switch.
    handler._response_writer().write(
        (
            "HTTP/1.1 101 Switching Protocols\r\n"
            "Upgrade: websocket\r\n"
            "Connection: Upgrade\r\n"
            f"Sec-WebSocket-Accept: {accept_key(client_key)}\r\n\r\n"
        ).encode("ascii")
    )

    sock = handler.connection
    # The channel reads the socket itself from here on. What the HTTP reader
    # already holds past the request is taken without waiting: with the
    # socket non-blocking, an empty buffer peeks as empty instead of reading.
    reader = cast("io.BufferedReader", handler._request_reader())
    sock.setblocking(False)
    try:
        received = reader.peek(1)
    except BlockingIOError:
        # Raised only from a read of the socket, made when nothing is held.
        received = b""
    reader.read(len(received))
    # The HTTP handler's idle timeout does not apply: a viewer channel is
    # idle between interactions, and the keep-alive pings bound a dead peer.
    sock.settimeout(None)
    websocket = WebSocket(
        sock,
        handler._response_writer(),
        received=received,
        max_text_bytes=MAX_EVENT_TEXT_BYTES,
    )
    server = handler._request_server()
    if isinstance(server, CellucidThreadingHTTPServer):
        server.close_with_server(sock, lambda: websocket.close(CLOSE_GOING_AWAY, "The server stopped"))
    serve_viewer_channel(websocket)


def serve_viewer_channel(websocket: ChannelSocket) -> None:
    """Run one viewer channel on ``websocket`` until it ends: its ``hello``,
    then every event, delivered on this thread in order."""
    registration: _Registration | None = None
    try:
        viewer_id, viewer_token, registration = _authenticate(websocket)
        if registration is None:
            return
        _attach(registration, websocket)
        threading.Thread(
            target=_keep_alive,
            args=(websocket,),
            name="cellucid-viewer-channel-keepalive",
            daemon=True,
        ).start()
        while True:
            message = websocket.receive()
            event: dict[str, Any] | None
            if isinstance(message, str):
                event = _parse_text_event(message)
            else:
                event = _receive_binary_event(message, viewer_id, viewer_token)
            if event is None:
                continue
            if websocket.closed:
                # Closed from another thread (replaced, stopped) while this
                # event was read: it is not the viewer's any more.
                _discard_session_file(event)
                return
            if not _deliver(websocket, viewer_id, viewer_token, event):
                return
    except WebSocketProtocolError as error:
        logger.warning("Closing a viewer channel: %s", error.reason)
        websocket.close(error.code, error.reason)
    except ViewerProtocolError as error:
        logger.warning("Closing a viewer channel: %s", error)
        websocket.close(CLOSE_POLICY_VIOLATION, str(error))
    except WebSocketClosed:
        pass
    finally:
        websocket.close(CLOSE_GOING_AWAY, "")
        if registration is not None:
            with registration.channel_lock:
                if registration.socket is websocket:
                    registration.socket = None
        # Read what the viewer still sends until it closes its side, so the
        # socket is not reset under the close frame (RFC 6455 section 7.1.1).
        websocket.finish()


def _authenticate(websocket: ChannelSocket) -> tuple[str, str, _Registration | None]:
    message = websocket.receive()
    if not isinstance(message, str):
        message.drain()
        raise ViewerProtocolError("The first frame must be the text frame 'hello'")
    hello = _json_object(message, label="hello")
    if set(hello) != _HELLO_FIELDS or hello["type"] != "hello":
        raise ViewerProtocolError("The first frame must be exactly {type: 'hello', viewerId, viewerToken}")
    viewer_id = hello["viewerId"]
    viewer_token = hello["viewerToken"]
    try:
        _require_routing_string(viewer_id, label="viewerId")
        _require_routing_string(viewer_token, label="viewerToken")
        registration = _require_registration(viewer_id, viewer_token)
    except (TypeError, ValueError, _ViewerNotRegisteredError, _ViewerAuthenticationError):
        logger.warning("Refused a viewer channel with unknown viewer credentials")
        websocket.close(CLOSE_UNAUTHORIZED, "Unknown viewer, or a wrong viewer token")
        return "", "", None
    return viewer_id, viewer_token, registration


def _attach(registration: _Registration, websocket: ChannelSocket) -> None:
    """Make ``websocket`` the viewer's, flush its outbox, retire the old socket."""
    with registration.channel_lock:
        previous = registration.socket
        registration.socket = websocket
        websocket.send_text(_WELCOME)
        while registration.outbox:
            # Sent before it leaves the outbox: a socket that dies here leaves
            # the command queued for the next one, and never sent twice whole.
            _send_command(websocket, registration.outbox[0])
            registration.outbox.popleft()
    if previous is not None and previous is not websocket:
        previous.close(CLOSE_REPLACED, "Replaced by a newer connection of this viewer")


def _keep_alive(websocket: ChannelSocket) -> None:
    while not websocket._closed.wait(KEEPALIVE_SECONDS):
        try:
            websocket.ping()
        except WebSocketClosed:
            return


def _json_object(text: str, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(text)
    except json.JSONDecodeError as error:
        raise ViewerProtocolError(f"The {label} frame is not JSON") from error
    if type(value) is not dict or type(value.get("type")) is not str:
        raise ViewerProtocolError(f"The {label} frame must be a JSON object with a type")
    return value


def _parse_text_event(text: str) -> dict[str, Any]:
    event = _json_object(text, label="event")
    if event["type"] in BINARY_EVENT_HEADER_FIELDS:
        raise ViewerProtocolError(f"A {event['type']} event is sent as a binary frame")
    return event


def _receive_binary_event(
    message: ChannelMessage, viewer_id: str, viewer_token: str
) -> dict[str, Any] | None:
    """Read one binary event; None when it was received and needs no delivery."""
    header_length = int.from_bytes(message.read_exact(4), "little")
    try:
        header_text = message.read_exact(header_length).decode("utf-8")
    except UnicodeDecodeError as error:
        raise ViewerProtocolError("A binary event header is not UTF-8") from error
    header = _json_object(header_text, label="binary event header")
    expected = BINARY_EVENT_HEADER_FIELDS.get(header["type"])
    if expected is None or set(header) != expected:
        raise ViewerProtocolError(
            "A binary event header must be exactly {type: 'selection', source}, "
            "{type: 'session_bundle', requestId} or a result's header"
        )
    if header["type"] == "result":
        # Its cells and values, read whole: `_wire._decode_result` checks them
        # against the header.
        return {**header, "body": message.read_rest()}
    if header["type"] == "selection":
        # One buffer, sized from the frame and filled in place: the cells are
        # read and unmasked without a second copy.
        cells = message.read_rest()
        return {"type": "selection", "cells": cells, "source": header["source"]}
    return _receive_session_bundle(message, viewer_id, viewer_token, header["requestId"])


def _receive_session_bundle(
    message: ChannelMessage, viewer_id: str, viewer_token: str, request_id: object
) -> dict[str, Any] | None:
    """Stream one session bundle to a temporary file and return its event."""
    if type(request_id) is not str:
        raise ViewerProtocolError("A session bundle requestId must be text")
    try:
        _session_bundle_request(viewer_id, viewer_token, request_id, consume=False)
    except (_SessionBundleRequestNotFoundError, _ViewerAuthenticationError):
        # Python stopped waiting (its timeout passed, or it was cancelled):
        # the bundle is read and dropped, and the channel stays open.
        logger.warning("Dropped a session bundle nobody is waiting for (request %s)", request_id)
        message.drain()
        return None

    tmp_path: Path | None = None
    written = 0
    try:
        magic = message.read_exact(len(SESSION_BUNDLE_MAGIC))
        if magic != SESSION_BUNDLE_MAGIC:
            raise ValueError(NOT_A_SESSION_FILE_MESSAGE)
        with tempfile.NamedTemporaryFile(
            delete=False, prefix="cellucid-session-", suffix=".cellucid-session"
        ) as tmp:
            tmp_path = Path(tmp.name)
            tmp.write(magic)
            written = len(magic)
            buffer = memoryview(bytearray(_BUNDLE_CHUNK_BYTES))
            while count := message.readinto(buffer):
                tmp.write(buffer[:count])
                written += count
        from .session_bundle import CellucidSessionBundle

        for _chunk in CellucidSessionBundle(tmp_path).iter_chunks(decode=True):
            pass
        _session_bundle_request(viewer_id, viewer_token, request_id, consume=True)
    except WebSocketClosed:
        # The channel ended mid-bundle (closed, or the network dropped): the
        # part written is nobody's.
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)
        raise
    except (ValueError, zlib.error, OSError, _SessionBundleRequestNotFoundError) as error:
        if not message.done:
            message.drain()
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)
        if isinstance(error, _SessionBundleRequestNotFoundError):
            logger.warning("Dropped a session bundle nobody is waiting for (request %s)", request_id)
            return None
        reason = f"The viewer sent an invalid session bundle: {error}"
        if isinstance(error, OSError):
            logger.exception("Failed to receive a session bundle")
            reason = f"The session bundle could not be saved: {error}"
        # Python is waiting for this request: tell it now, in the form its
        # wait already accepts, instead of letting it time out.
        return {"type": "command_error", "command": "requestSessionBundle", "reason": reason[:500]}
    return {
        "type": "session_bundle",
        "requestId": request_id,
        "status": "ok",
        "bytes": written,
        "path": str(tmp_path),
    }


def _deliver(websocket: ChannelSocket, viewer_id: str, viewer_token: str, event: dict[str, Any]) -> bool:
    """Deliver one event; False when the channel must end."""
    try:
        route_event(viewer_id, viewer_token, event)
    except (_ViewerNotRegisteredError, _ViewerAuthenticationError):
        websocket.close(CLOSE_GOING_AWAY, "The viewer was stopped")
        _discard_session_file(event)
        return False
    except ViewerProtocolError:
        _discard_session_file(event)
        raise
    except Exception:
        # A notebook hook failed. That is the notebook's error, reported in its
        # log; the viewer did nothing wrong, so the channel stays open.
        logger.exception("A viewer event hook failed")
    return True


def _discard_session_file(event: dict[str, Any]) -> None:
    if event.get("type") == "session_bundle" and isinstance(event.get("path"), str):
        Path(event["path"]).unlink(missing_ok=True)
