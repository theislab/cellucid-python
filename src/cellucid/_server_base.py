"""
Shared building blocks of the Cellucid HTTP server.

The HTTP server class, byte ranges, CORS and ``Host`` checks, the notebook
event and session-bundle routing, the web-build cache paths and the console
banner that :class:`cellucid.anndata_server.AnnDataServer` is built from.

It's an internal module (prefixed with _) and not part of the public API.
"""

from __future__ import annotations

import hmac
import ipaddress
import json
import logging
import math
import re
import shutil
import socket
import stat
import sys
import tempfile
import threading
import time
import zlib
from collections.abc import Callable, Sequence
from html.parser import HTMLParser
from http import HTTPStatus
from http.server import HTTPServer as _StandardHTTPServer
from pathlib import Path
from socketserver import ThreadingMixIn
from typing import TYPE_CHECKING, Any, BinaryIO, NoReturn, TypeGuard, cast

from ._console import console_print
from .session_bundle import NOT_A_SESSION_FILE_MESSAGE as _NOT_A_SESSION_FILE_MESSAGE
from .session_bundle import SESSION_BUNDLE_MAGIC as _SESSION_BUNDLE_MAGIC

if TYPE_CHECKING:
    from http.client import HTTPMessage
    from http.server import BaseHTTPRequestHandler, HTTPServer

logger = logging.getLogger("cellucid.server")

DEFAULT_PORT = 8765
DEFAULT_HOST = "127.0.0.1"
CELLUCID_WEB_URL = "https://www.cellucid.com"
WEB_ASSET_INVENTORY_FILENAME = "cellucid-web-assets.json"

_EXCLUSIVE_ADDRESS_USE: int | None = getattr(socket, "SO_EXCLUSIVEADDRUSE", None)


class CellucidHTTPServer(_StandardHTTPServer):
    """Own a listening port deterministically on every supported platform.

    The standard-library ``HTTPServer`` enables ``SO_REUSEADDR``. On Unix that
    permits a prompt restart without allowing two live TCP listeners to share
    an address. Winsock gives the option different semantics: another process
    can bind the same live address, after which connection routing is
    indeterminate. Windows instead provides ``SO_EXCLUSIVEADDRUSE`` for server
    sockets that must reject an occupied port and prevent later port hijacking.
    """

    allow_reuse_address = _EXCLUSIVE_ADDRESS_USE is None
    allow_reuse_port = False

    def server_bind(self) -> None:
        if _EXCLUSIVE_ADDRESS_USE is not None:
            self.socket.setsockopt(
                socket.SOL_SOCKET,
                _EXCLUSIVE_ADDRESS_USE,
                1,
            )
        super().server_bind()

    def handle_error(self, request: Any, client_address: Any) -> None:
        """Report a request failure, except a client that simply went away.

        A browser that closes a tab, reloads, or navigates away mid-download
        aborts the connection, and the write in progress raises. That is the
        client's decision, not a server fault, and the larger the payload the
        more often it happens -- a points file for millions of cells is many
        megabytes, and any reload during it lands here. Reporting it as an
        unexpected error printed a traceback that reads like a crash.
        """
        if isinstance(sys.exc_info()[1], ConnectionError):
            logger.debug("Client closed the connection from %s", client_address)
            return
        super().handle_error(request, client_address)


class CellucidThreadingHTTPServer(ThreadingMixIn, CellucidHTTPServer):
    """A :class:`CellucidHTTPServer` that answers every connection in its own thread.

    A browser reading a file in place keeps several range requests in flight
    and reuses its connections (HTTP/1.1 keep-alive); one thread would answer
    them one connection at a time. Handler threads are daemons, so an idle
    kept-alive connection never holds up shutdown.
    """

    daemon_threads = True
    block_on_close = False
    # A page opens its six connections to a host at once (and a notebook's
    # proxies more); the standard backlog of five would make one wait for a
    # connection retry, so the platform's own maximum is asked for.
    request_queue_size = socket.SOMAXCONN


_BYTE_RANGE_PATTERN = re.compile(r"bytes=(\d{0,20})-(\d{0,20})")


def parse_byte_range(value: str | None, size: int) -> tuple[int, int] | None:
    """Resolve one ``Range`` header against a resource of ``size`` bytes.

    Returns the inclusive ``(start, end)`` interval, or None when there is no
    header. One interval is accepted: ``bytes=a-b``, ``bytes=a-`` or the suffix
    ``bytes=-n``. Anything else, or an interval outside the resource, raises
    ``ValueError`` and is answered 416.
    """
    if value is None:
        return None
    if type(size) is not int or size < 0:
        raise ValueError("Resource size must be a non-negative integer.")
    match = _BYTE_RANGE_PATTERN.fullmatch(value)
    if match is None or (not match.group(1) and not match.group(2)):
        raise ValueError("Range must contain one exact byte interval.")
    if size == 0:
        raise ValueError("An empty resource has no satisfiable byte range.")
    if match.group(1):
        start = int(match.group(1))
        end = int(match.group(2)) if match.group(2) else size - 1
        if start >= size or end < start:
            raise ValueError("Range is outside the resource.")
        return start, min(end, size - 1)
    suffix_length = int(match.group(2))
    if suffix_length <= 0:
        raise ValueError("Range suffix length must be positive.")
    return max(0, size - suffix_length), size - 1


# Session bundle upload protocol (used for Jupyter "no download" capture): the
# upload is checked against the one signature the session reader accepts.
SESSION_BUNDLE_MAGIC = _SESSION_BUNDLE_MAGIC


def _session_bundle_capacity() -> int:
    """The largest session upload this machine can hold: the free space of the
    temporary directory the upload streams into. The body is never held in
    memory, and every chunk inside it is bounded by its declared size when it
    is read, so free disk is the only limit on the bundle as a whole."""
    return shutil.disk_usage(tempfile.gettempdir()).free

# Selection delivery (``POST /_cellucid/selection``): the body is the selected
# cell indices as little-endian uint32, strictly ascending, so no body can hold
# more than every uint32 index once. That is the wire format's own bound, not a
# size policy: any selection the viewer can make fits, and the server's memory
# is the only other limit (a body it cannot hold is refused with 413).
MAX_SELECTION_BODY_BYTES = 4 * 2**32
#: Longest JSON event record the viewer posts to ``/_cellucid/events``. Every
#: event there is a small record (hover, click, mapping change); selections and
#: session bundles, the only large payloads, have their own routes.
MAX_EVENT_BODY_BYTES = 1024 * 1024
# Event types that only the server may build from a dedicated binary route.
# ``POST /_cellucid/events`` refuses them, so the JSON endpoint can neither
# forge one nor carry a selection past its 1 MiB body limit.
_BINARY_ROUTED_EVENT_PATHS = {
    "selection": "/_cellucid/selection",
    "session_bundle": "/_cellucid/session_bundle",
}
# Longest a refused request's unread body is drained before the socket closes.
_UNREAD_BODY_DISCARD_SECONDS = 5.0

_WEB_BUILD_META_NAME = "cellucid-web-build-id"
_WEB_BUILD_ID_MAX_LEN = 200


class _WebBuildMetaParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.build_ids: list[str] = []

    def handle_starttag(
        self,
        tag: str,
        attrs: list[tuple[str, str | None]],
    ) -> None:
        if tag != "meta":
            return
        names = [value for name, value in attrs if name == "name"]
        if _WEB_BUILD_META_NAME not in names:
            return
        contents = [value for name, value in attrs if name == "content"]
        if names != [_WEB_BUILD_META_NAME] or len(contents) != 1:
            raise ValueError(
                "cellucid-web-build-id meta must contain one name and one content attribute"
            )
        content = contents[0]
        if content is None:
            raise ValueError("cellucid-web-build-id meta content must have a value")
        self.build_ids.append(content)


def _extract_web_build_id(index_html_bytes: bytes) -> str:
    """Require and return the single exact build stamp in ``index.html``."""
    if not isinstance(index_html_bytes, bytes):
        raise TypeError("index.html content must be bytes")
    try:
        text = index_html_bytes.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("index.html must be valid UTF-8") from error

    parser = _WebBuildMetaParser()
    parser.feed(text)
    parser.close()
    if len(parser.build_ids) != 1:
        raise ValueError("index.html must contain exactly one cellucid-web-build-id meta element")
    build_id = parser.build_ids[0]
    if not build_id or build_id != build_id.strip() or len(build_id) > _WEB_BUILD_ID_MAX_LEN:
        raise ValueError("index.html contains an invalid cellucid-web-build-id")
    return build_id


def _require_web_cache_directory_or_missing(cache_dir: Path) -> bool:
    """Require an exact directory or a missing path without following the leaf."""
    try:
        metadata = cache_dir.lstat()
    except FileNotFoundError:
        return False
    is_reparse_point = bool(
        getattr(metadata, "st_file_attributes", 0) & stat.FILE_ATTRIBUTE_REPARSE_POINT
    )
    if stat.S_ISLNK(metadata.st_mode) or is_reparse_point:
        raise ValueError(
            f"Web cache path must not be a symbolic link or reparse point: {cache_dir}"
        )
    if not stat.S_ISDIR(metadata.st_mode):
        raise NotADirectoryError(f"Web cache path must be a directory: {cache_dir}")
    return True


def _remove_web_cache_dir(cache_dir: Path) -> None:
    """Remove one exact cache directory and propagate every failure."""
    if not _require_web_cache_directory_or_missing(cache_dir):
        return
    shutil.rmtree(cache_dir)


def _web_cache_dir() -> Path:
    """Return the process-independent default cache location."""
    return Path(tempfile.gettempdir()) / "cellucid-web-cache"


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
    ttl_seconds: float = 60.0,
) -> None:
    """Register one authenticated, one-shot session bundle upload request."""
    _require_routing_string(viewer_id, label="viewer_id")
    _require_routing_string(viewer_token, label="viewer_token")
    _require_routing_string(request_id, label="request_id")
    if isinstance(ttl_seconds, bool) or not isinstance(ttl_seconds, int | float):
        raise TypeError("ttl_seconds must be a number")
    if not math.isfinite(ttl_seconds) or ttl_seconds <= 0:
        raise ValueError("ttl_seconds must be finite and greater than zero")

    now = time.monotonic()
    expires_at = now + ttl_seconds
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


def _require_session_bundle_request(
    viewer_id: str,
    viewer_token: str,
    request_id: str,
) -> None:
    """Authenticate one live pending request without consuming it."""
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


def _consume_session_bundle_request(
    viewer_id: str,
    viewer_token: str,
    request_id: str,
) -> None:
    """Atomically authenticate and consume one live pending request."""
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
        del _pending_session_bundle_requests[key]


def cancel_session_bundle_request(
    viewer_id: str,
    viewer_token: str,
    request_id: str,
) -> bool:
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
# CLI Output Helpers
# =============================================================================
# These functions provide consistent, professional CLI output across all servers.


def print_step(step: int, total: int, title: str):
    """Print a step header like [1/4] Title..."""
    console_print(f"\n[{step}/{total}] {title}...")


def print_detail(label: str, value: str):
    """Print an indented detail line."""
    console_print(f"      {label}: {value}")


def print_success(message: str = "Done"):
    """Print a success indicator with checkmark."""
    console_print(f"      ✓ {message}")


def print_server_banner(
    url: str,
    viewer_url: str,
    *,
    network_urls: Sequence[tuple[str, str]] = (),
):
    """Print the final server ready banner.

    ``network_urls`` names the addresses a *different* machine can reach, which
    a server bound to every interface has and a loopback server does not. The
    bind wildcard itself is never printed as a URL: ``0.0.0.0`` is a statement
    about which interfaces accept connections, and no browser can open it.
    """
    console_print(f"\n{'═' * 60}")
    console_print("  CELLUCID SERVER RUNNING")
    console_print(f"{'═' * 60}")
    console_print()
    console_print(f"  Local URL:    {url}")
    console_print(f"  Viewer URL:   {viewer_url}")
    if network_urls:
        console_print()
        console_print("  Bound to every network interface. From another machine:")
        for label, value in network_urls:
            console_print(f"  {label}{value}")
    console_print()
    console_print("  Press Ctrl+C to stop")
    console_print()
    console_print(f"{'═' * 60}")


class _ViewerNotRegisteredError(LookupError):
    """The supplied viewer id has no live callback registration."""


class _ViewerAuthenticationError(PermissionError):
    """The supplied viewer token does not authenticate the viewer."""


class _SessionBundleRequestNotFoundError(LookupError):
    """The supplied session request is absent or expired."""


class _RequestContractError(ValueError):
    """An HTTP request does not satisfy the exact endpoint contract."""

    def __init__(self, status: int, message: str):
        super().__init__(message)
        self.status = status


class _EventCallbackRegistration:
    __slots__ = ("callback", "delivery_lock", "viewer_token")

    def __init__(
        self,
        viewer_token: str,
        callback: Callable[[dict], None],
    ):
        self.viewer_token = viewer_token
        self.callback = callback
        self.delivery_lock = threading.RLock()


_DUMMY_VIEWER_TOKEN = "\0" * 32
_event_callbacks_lock = threading.Lock()
_event_callbacks: dict[str, _EventCallbackRegistration] = {}


def _require_routing_string(value: object, *, label: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string")
    if not value or value != value.strip() or any(character.isspace() for character in value):
        raise ValueError(f"{label} must be a non-empty string without whitespace")
    return value


def _require_event_registration_locked(
    viewer_id: str,
    viewer_token: str,
) -> _EventCallbackRegistration:
    registration = _event_callbacks.get(viewer_id)
    expected_token = registration.viewer_token if registration is not None else _DUMMY_VIEWER_TOKEN
    authenticated = hmac.compare_digest(expected_token, viewer_token)
    if registration is None:
        raise _ViewerNotRegisteredError
    if not authenticated:
        raise _ViewerAuthenticationError
    return registration


def register_event_callback(
    viewer_id: str,
    viewer_token: str,
    callback: Callable[[dict], None],
) -> None:
    """Register one token-authenticated callback for a viewer."""
    _require_routing_string(viewer_id, label="viewer_id")
    _require_routing_string(viewer_token, label="viewer_token")
    if not callable(callback):
        raise TypeError("callback must be callable")
    registration = _EventCallbackRegistration(viewer_token, callback)
    with _event_callbacks_lock:
        if viewer_id in _event_callbacks:
            raise RuntimeError("An event callback is already registered for this viewer_id")
        _event_callbacks[viewer_id] = registration


def unregister_event_callback(viewer_id: str, viewer_token: str) -> bool:
    """Unregister only the callback authenticated by the supplied token."""
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
            current = _event_callbacks.get(viewer_id)
            if current is not registration:
                return False
            del _event_callbacks[viewer_id]

        with _pending_session_bundle_lock:
            matching_requests = [
                key
                for key, request in _pending_session_bundle_requests.items()
                if key[0] == viewer_id and hmac.compare_digest(request.viewer_token, viewer_token)
            ]
            for key in matching_requests:
                del _pending_session_bundle_requests[key]
        return True


def route_event(
    viewer_id: str,
    viewer_token: str,
    event: dict,
) -> None:
    """Authenticate and deliver one event exactly once."""
    _require_routing_string(viewer_id, label="viewer_id")
    if not isinstance(viewer_token, str):
        raise TypeError("viewer_token must be a string")
    if not isinstance(event, dict):
        raise TypeError("event must be a dictionary")

    with _event_callbacks_lock:
        registration = _require_event_registration_locked(
            viewer_id,
            viewer_token,
        )

    with registration.delivery_lock:
        with _event_callbacks_lock:
            if _event_callbacks.get(viewer_id) is not registration:
                raise _ViewerNotRegisteredError
        registration.callback(event)


def require_server_port(port: object) -> int:
    """Validate one exact TCP bind port.

    Port ``0`` is the explicit request for one operating-system-assigned port.
    Every nonzero port is bound exactly as supplied.
    """
    if type(port) is not int:
        raise TypeError("port must be an integer")
    if port < 0 or port > 65535:
        raise ValueError("port must be between 0 and 65535")
    return port


def _require_exact_content_type(headers: HTTPMessage, expected: str) -> None:
    values = headers.get_all("Content-Type", [])
    if values != [expected]:
        raise _RequestContractError(
            HTTPStatus.UNSUPPORTED_MEDIA_TYPE,
            f"Content-Type must be exactly {expected}",
        )


def _require_content_length(
    headers: HTTPMessage,
    *,
    maximum: int,
) -> int:
    if headers.get_all("Transfer-Encoding", []):
        raise _RequestContractError(
            HTTPStatus.BAD_REQUEST,
            "Transfer-Encoding is not supported",
        )
    values = headers.get_all("Content-Length", [])
    if not values:
        raise _RequestContractError(
            HTTPStatus.LENGTH_REQUIRED,
            "Missing Content-Length",
        )
    if len(values) != 1:
        raise _RequestContractError(
            HTTPStatus.BAD_REQUEST,
            "Content-Length must appear exactly once",
        )
    raw_value = values[0]
    if (
        not raw_value
        or not raw_value.isascii()
        or not raw_value.isdigit()
        or (len(raw_value) > 1 and raw_value.startswith("0"))
    ):
        raise _RequestContractError(
            HTTPStatus.BAD_REQUEST,
            "Content-Length must be a canonical decimal integer",
        )
    value = int(raw_value)
    if value <= 0:
        raise _RequestContractError(
            HTTPStatus.LENGTH_REQUIRED,
            "Content-Length must be greater than zero",
        )
    if value > maximum:
        raise _RequestContractError(
            HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
            f"Request body exceeds the {maximum}-byte limit",
        )
    return value


def _require_viewer_credentials(viewer_id: str, viewer_token: str) -> None:
    """Authenticate one viewer without delivering anything to it.

    A route that accepts a large body checks the credentials first, so an
    unauthenticated caller is refused before the server reads that body.
    :func:`route_event` authenticates again when it delivers.
    """
    with _event_callbacks_lock:
        _require_event_registration_locked(viewer_id, viewer_token)


def _parse_exact_query(query: str, fields: frozenset[str], *, label: str) -> dict[str, str]:
    """Parse one request query that may name only ``fields``, each at most once."""
    from urllib.parse import parse_qsl

    try:
        pairs = parse_qsl(query, keep_blank_values=True, strict_parsing=True)
    except ValueError as error:
        raise _RequestContractError(
            HTTPStatus.BAD_REQUEST,
            f"{label} query is malformed",
        ) from error
    values: dict[str, str] = {}
    for key, value in pairs:
        if key not in fields or key in values:
            raise _RequestContractError(
                HTTPStatus.BAD_REQUEST,
                f"{label} query must contain each declared field exactly once",
            )
        values[key] = value
    return values


def _is_exact_route_text(value: object) -> TypeGuard[str]:
    return (
        isinstance(value, str)
        and bool(value)
        and value == value.strip()
        and not any(character.isspace() for character in value)
    )


def _read_exact_request_body(read: Callable[[int], bytes], content_length: int) -> bytes:
    try:
        body = read(content_length)
    except MemoryError:
        raise _RequestContractError(
            HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
            f"The server could not allocate {content_length} bytes for the request body",
        ) from None
    if len(body) != content_length:
        raise _RequestContractError(
            HTTPStatus.BAD_REQUEST,
            "Request body is shorter than Content-Length",
        )
    return body


# The Cellucid servers only ever speak cleartext HTTP, so an authority that
# omits its port names port 80.
_HTTP_DEFAULT_PORT = 80
_LOOPBACK_HOST_NAME = "localhost"


def _parse_authority_port(text: str) -> int | None:
    """Return one exact decimal authority port, or ``None`` when malformed."""
    if not text or not text.isascii() or not text.isdigit():
        return None
    if len(text) > 1 and text.startswith("0"):
        return None
    port = int(text)
    if port < 1 or port > 65535:
        return None
    return port


def _parse_host_authority(value: object) -> tuple[str, int | None] | None:
    """Split one ``Host`` header into ``(host, port)``, or ``None`` when malformed.

    ``host`` is returned unbracketed so an IPv6 literal can be classified
    directly. A bare ``split(":")`` cannot be used here because an IP-literal
    authority is bracketed and contains colons of its own.
    """
    if not isinstance(value, str) or not value or not value.isascii():
        return None
    if any(ord(character) <= 32 or ord(character) == 127 for character in value):
        return None
    if "@" in value or "/" in value:
        return None

    if value.startswith("["):
        closing = value.find("]")
        if closing < 0:
            return None
        host = value[1:closing]
        # A scope identifier never names a browser-reachable origin.
        if "%" in host:
            return None
        try:
            ipaddress.IPv6Address(host)
        except ValueError:
            return None
        remainder = value[closing + 1 :]
        if not remainder:
            return host, None
        if not remainder.startswith(":"):
            return None
        port_text = remainder[1:]
    else:
        if "[" in value or "]" in value:
            return None
        host, separator, port_text = value.partition(":")
        if not host:
            return None
        if not separator:
            return host, None

    port = _parse_authority_port(port_text)
    if port is None:
        return None
    return host, port


# A DNS name is at most 253 characters over dot-separated labels of at most 63.
_MAXIMUM_HOST_NAME_LENGTH = 253
_MAXIMUM_HOST_LABEL_LENGTH = 63
_HOST_LABEL_CHARACTERS = frozenset("abcdefghijklmnopqrstuvwxyz0123456789-")


def _is_dns_host_name(lowered: str) -> bool:
    """Report whether one lowercase string is a plain LDH DNS name.

    Only letters, digits and hyphens are accepted, because those are the
    characters a browser can put in an authority without percent-encoding.
    """
    if not lowered or len(lowered) > _MAXIMUM_HOST_NAME_LENGTH:
        return False
    for label in lowered.split("."):
        if not label or len(label) > _MAXIMUM_HOST_LABEL_LENGTH:
            return False
        if label.startswith("-") or label.endswith("-"):
            return False
        if not _HOST_LABEL_CHARACTERS.issuperset(label):
            return False
    return True


def _normalize_host_name(host: str) -> str | None:
    """Return the exact comparison form of one authority host, or ``None``.

    Both an allowlist entry and an incoming ``Host`` are reduced through this
    one function, so ``[0:0:0:0:0:0:0:1]`` and ``::1`` — and ``Hub.Example.Org``
    and ``hub.example.org`` — cannot be two different answers to the same
    question.
    """
    try:
        return ipaddress.ip_address(host).compressed
    except ValueError:
        pass
    lowered = host.lower()
    return lowered if _is_dns_host_name(lowered) else None


def _require_allowed_host_entry(entry: str, *, label: str) -> str:
    """Return one validated allowlist entry in comparison form.

    Every shape this guard does not implement is refused by name rather than
    ignored: an entry that is silently dropped, or silently widened, is the
    original rebinding hole with a configuration file in front of it.
    """

    def refuse(reason: str) -> NoReturn:
        raise ValueError(f"{label} entry {entry!r} {reason}")

    if not entry:
        refuse("is empty; give one host name such as 'hub.example.org'")
    if not entry.isascii():
        refuse("is not ASCII; use the Punycode (IDNA) form of an international name")
    if any(character.isspace() for character in entry):
        refuse("contains whitespace; give exactly one host name per entry")
    if any(ord(character) <= 32 or ord(character) == 127 for character in entry):
        refuse("contains a control character; give one plain host name")
    if "*" in entry:
        refuse(
            "uses a wildcard; a wildcard re-admits DNS rebinding under any name it "
            "covers, so every name must be written out in full"
        )
    if "/" in entry:
        refuse("looks like a URL; give the host name alone, with no scheme and no path")
    if "@" in entry:
        refuse("carries credentials; give the host name alone")
    if "[" in entry or "]" in entry:
        refuse("is bracketed; write an IPv6 address bare, for example '::1'")
    if "%" in entry:
        refuse("carries a zone identifier, which never names a browser-reachable host")

    try:
        return ipaddress.ip_address(entry).compressed
    except ValueError:
        pass
    if ":" in entry:
        refuse(
            "carries a port; write the host name alone, because an allowed name is "
            "matched on whatever port the proxy in front of it publishes"
        )
    lowered = entry.lower()
    if lowered.endswith("."):
        refuse("ends with a dot; write the name without its root label")
    if not _is_dns_host_name(lowered):
        refuse(
            "is not a host name; use dot-separated labels of letters, digits and "
            "hyphens, or one bare IP address"
        )
    return lowered


def require_allowed_hosts(
    value: object,
    *,
    label: str = "allowed_hosts",
) -> frozenset[str]:
    """Validate the extra ``Host`` names one server is allowed to answer to.

    Cellucid binds loopback and answers only its own loopback authority, which
    is what stops DNS rebinding. A reverse proxy — ``jupyter-server-proxy`` in
    particular — forwards the front-end's ``Host`` verbatim, and that ``Host``
    is indistinguishable from a rebound one: both arrive on loopback naming a
    foreign authority, and a rebound page is same-origin so it can forge any
    ``X-Forwarded-*`` header it likes. The proxy therefore cannot be detected,
    only declared, which is what this list is.

    Every entry is one bare host name: a DNS name or an IP address literal. No
    port, no scheme, no path, no wildcard. An empty list — the default — leaves
    the loopback-only rule exactly as it is.
    """
    if value is None:
        return frozenset()
    if isinstance(value, str | bytes):
        raise TypeError(
            f"{label} must be a sequence of host names, not one string; "
            f"pass [{value!r}] to allow exactly that name"
        )
    if not isinstance(value, Sequence):
        raise TypeError(f"{label} must be a list or tuple of host names")
    names: set[str] = set()
    for entry in value:
        if type(entry) is not str:
            raise TypeError(f"{label} entries must be strings, got {type(entry).__name__}")
        names.add(_require_allowed_host_entry(entry, label=label))
    return frozenset(names)


def _host_authority_is_loopback(host: str) -> bool:
    """Report whether one authority host names only the loopback interface.

    IP literals are accepted whenever they are loopback: DNS cannot rebind a
    literal, so such an authority can only come from a page whose own origin is
    that literal. ``localhost`` is the one name accepted, because RFC 6761
    reserves it for the loopback interface; every other name is refused, which
    is exactly what defeats DNS rebinding.
    """
    try:
        address = ipaddress.ip_address(host)
    except ValueError:
        return host.lower() == _LOOPBACK_HOST_NAME
    return address.is_loopback


def _requires_loopback_host_header(bound_host: object) -> bool:
    """Report whether one bound address must serve only loopback authorities.

    Enforcement is skipped only when the operator positively bound a
    non-loopback address (``--host 0.0.0.0`` and friends). That is an explicit
    request for network exposure, and the legitimate names of an exposed
    deployment are unknowable from here; anything undetermined is enforced.
    """
    if not isinstance(bound_host, str):
        return True
    try:
        address = ipaddress.ip_address(bound_host)
    except ValueError:
        return True
    return address.is_loopback


def bind_host_is_wildcard(bound_host: object) -> bool:
    """Report whether one bind address names every interface rather than one.

    ``0.0.0.0``, ``::``, and the empty string all tell the operating system to
    accept connections on every interface. None of them is an address a client
    can connect *to*, so none of them may ever be rendered inside a URL.
    """
    if not isinstance(bound_host, str):
        return False
    if not bound_host:
        return True
    try:
        return ipaddress.ip_address(bound_host).is_unspecified
    except ValueError:
        return False


def format_http_origin(host: str, port: int) -> str:
    """Return one http origin, bracketing an IPv6 literal as a URL requires."""
    try:
        is_ipv6 = ipaddress.ip_address(host).version == 6
    except ValueError:
        is_ipv6 = False
    return f"http://[{host}]:{port}" if is_ipv6 else f"http://{host}:{port}"


def loopback_origin_for_bind(bound_host: str, port: int) -> str:
    """Return the origin a browser on the serving machine itself can open.

    A wildcard bind accepts loopback too, so the machine running the server --
    and anything tunnelling into it, which is how a remote session reaches it --
    always has the loopback origin available.
    """
    try:
        is_ipv6 = ipaddress.ip_address(bound_host).version == 6
    except ValueError:
        is_ipv6 = False
    return format_http_origin("::1" if is_ipv6 else "127.0.0.1", port)


def machine_origin(port: int) -> str | None:
    """Return this machine's own named origin, or None when it has no usable one.

    This is the address another machine on the same network uses. A hostname
    that resolves to loopback names nothing outside this machine, and is
    reported as no origin at all rather than as a URL that cannot travel.
    """
    try:
        host_name = socket.gethostname()
    except OSError:
        return None
    normalized = _normalize_host_name(host_name)
    if normalized is None or _host_authority_is_loopback(normalized):
        return None
    return format_http_origin(normalized, port)


class CORSMixin:
    """
    Mixin class providing CORS headers for HTTP handlers.

    Use in conjunction with SimpleHTTPRequestHandler or similar.
    """

    # Override in subclass: whether to allow caching
    allow_caching: bool = True

    # Body bytes the current request's handler has consumed.
    _request_body_bytes_read: int = 0

    if TYPE_CHECKING:
        allowed_hosts: frozenset[str]
        close_connection: bool
        command: str
        path: str
        serve_web_ui: bool
        web_cache_dir: Path

        def end_headers(self) -> None: ...

        def send_header(self, keyword: str, value: object) -> None: ...

        def send_response(
            self,
            code: int | HTTPStatus,
            message: str | None = None,
        ) -> None: ...

    def parse_request(self) -> bool:
        """Refuse foreign ``Host`` authorities before any route is dispatched.

        A page on any origin can publish a short-lived DNS record that
        re-resolves its own hostname to ``127.0.0.1`` (DNS rebinding). The
        browser then considers this server same-origin, so no CORS check runs
        and the response body — the user's private single-cell data — is handed
        to that page. Validating ``Host`` is the only defence, and this is the
        one boundary every request crosses: ``handle_one_request()`` dispatches
        to ``do_GET``/``do_HEAD``/``do_POST``/``do_OPTIONS`` only after
        ``parse_request()`` returns ``True``.
        """
        base = cast("BaseHTTPRequestHandler", super())
        if not base.parse_request():
            return False
        self._request_body_bytes_read = 0
        if self._host_header_is_authorized():
            return True
        # Requests are otherwise logged at DEBUG, so a refusal would be silent.
        # The offending value goes only to the operator's log, never into the
        # response, and %r keeps a folded header out of the log's line framing.
        logger.warning(
            "Refused a request whose Host header does not name this server: %r",
            self._request_headers().get_all("Host", []),
        )
        self.send_error_response(
            HTTPStatus.MISDIRECTED_REQUEST,
            "Host header does not name this server",
        )
        self._discard_unread_request_body()
        return False

    def _host_header_is_authorized(self) -> bool:
        """Report whether ``Host`` names this server's own or a declared authority."""
        server = self._request_server()
        allowed = self._allowed_host_names()
        # Without a declaration this is exactly the loopback-only rule: a
        # deliberate non-loopback bind is an explicit request for network
        # exposure whose legitimate names cannot be known from here. Declaring
        # names is that missing knowledge, so it turns the check on everywhere.
        if not allowed and not _requires_loopback_host_header(server.server_address[0]):
            return True

        values = self._request_headers().get_all("Host", [])
        if len(values) != 1:
            return False
        authority = _parse_host_authority(values[0])
        if authority is None:
            return False
        host, port = authority
        if allowed:
            normalized = _normalize_host_name(host)
            # A declared name is matched on any port: the operator's proxy
            # publishes a front-end port that has nothing to do with the
            # loopback port bound here, and the port in an attacker's Host is
            # attacker-chosen anyway, so pinning it would only misconfigure.
            if normalized is not None and normalized in allowed:
                return True
        # ``server_port`` is read back from the bound socket, so a server that
        # asked for port 0 is checked against the port it was actually given.
        if (port if port is not None else _HTTP_DEFAULT_PORT) != server.server_port:
            return False
        return _host_authority_is_loopback(host)

    def _allowed_host_names(self) -> frozenset[str]:
        """Return the operator-declared extra authorities, empty unless declared.

        A handler built without a validated declaration must land on the
        strictest rule, never the weakest, so anything that is not the exact
        frozenset ``require_allowed_hosts`` produces counts as no declaration.
        """
        names = getattr(self, "allowed_hosts", None)
        return names if type(names) is frozenset else frozenset()

    def _request_headers(self) -> HTTPMessage:
        return cast(
            "HTTPMessage",
            object.__getattribute__(self, "headers"),
        )

    def _request_server(self) -> HTTPServer:
        return cast(
            "HTTPServer",
            object.__getattribute__(self, "server"),
        )

    def _request_reader(self) -> BinaryIO:
        return cast(
            BinaryIO,
            object.__getattribute__(self, "rfile"),
        )

    def _response_writer(self) -> BinaryIO:
        return cast(
            BinaryIO,
            object.__getattribute__(self, "wfile"),
        )

    def _read_request_body(self, size: int) -> bytes:
        """Read up to ``size`` request-body bytes, counting what is consumed."""
        chunk = self._request_reader().read(size)
        self._request_body_bytes_read += len(chunk)
        return chunk

    def _discard_unread_request_body(self) -> None:
        """Drain the body a refused request declared, so its socket closes cleanly.

        A socket closed with received bytes still unread is reset rather than
        shut down, and a client that has not yet read the response loses it:
        Windows reports WinError 10053, POSIX a reset or a broken pipe while the
        client is still sending. Every early POST refusal -- an expired session
        request, bad credentials, a wrong Content-Type, an oversized body, an
        unknown path, a foreign ``Host`` -- answers before reading the body, so
        without this the client saw a dropped connection instead of the status
        it was sent. Only a body framed by one valid ``Content-Length`` within
        the largest accepted upload is drained, and for a bounded time.
        """
        self.close_connection = True
        try:
            declared = _require_content_length(
                self._request_headers(),
                maximum=_session_bundle_capacity(),
            )
        except _RequestContractError:
            return
        remaining = declared - self._request_body_bytes_read
        deadline = time.monotonic() + _UNREAD_BODY_DISCARD_SECONDS
        connection = cast(socket.socket, object.__getattribute__(self, "connection"))
        try:
            while remaining > 0:
                time_left = deadline - time.monotonic()
                if time_left <= 0:
                    return
                connection.settimeout(time_left)
                chunk = self._read_request_body(min(64 * 1024, remaining))
                if not chunk:
                    return
                remaining -= len(chunk)
        except OSError:
            # The client went away or stalled; the close that follows is all
            # that is left to do.
            return

    def _get_allowed_origin(self) -> str | None:
        # `BaseHTTPRequestHandler` only populates `.headers` after successfully
        # parsing a request. If parsing fails, `send_error()` may call our
        # overridden `end_headers()` while `.headers` is still unset.
        headers = cast("HTTPMessage | None", getattr(self, "headers", None))
        if headers is None:
            return None

        origin = headers.get("Origin")
        if not isinstance(origin, str) or not origin:
            return None

        from urllib.parse import urlparse

        try:
            parsed = urlparse(origin)
        except Exception:
            return None

        if parsed.scheme not in ("http", "https"):
            return None

        host = parsed.hostname
        if not host:
            return None

        if host in ("localhost", "127.0.0.1"):
            return origin

        # Allow the canonical hosted viewer origin.
        if origin == CELLUCID_WEB_URL:
            return origin

        # Allow https subdomains under cellucid.com (optional, but safe).
        if parsed.scheme == "https" and (host == "cellucid.com" or host.endswith(".cellucid.com")):
            return origin

        return None

    def add_cors_headers(self):
        """Add CORS headers to the response."""
        allowed_origin = self._get_allowed_origin()
        if allowed_origin:
            self.send_header("Access-Control-Allow-Origin", allowed_origin)
            self.send_header("Vary", "Origin")
        self.send_header("Access-Control-Allow-Methods", "GET, HEAD, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type, Range")
        # Every header a cross-origin reader is expected to act on has to be
        # listed here, because a header that is sent but not exposed reads as
        # absent. `Content-Encoding` was sent and not exposed, so the browser
        # could not tell a gzipped body from a plain one and sized its progress
        # from the compressed length.
        self.send_header(
            "Access-Control-Expose-Headers",
            "Accept-Ranges, Content-Encoding, Content-Length, Content-Range, ETag",
        )
        if self.allow_caching:
            self.send_header("Cache-Control", "public, max-age=3600")
        else:
            self.send_header("Cache-Control", "no-cache")

    def do_OPTIONS(self):
        """Handle CORS preflight requests."""
        self.send_response(HTTPStatus.NO_CONTENT)
        # Note: CORS headers are added by end_headers() override in subclasses
        self.end_headers()

    def send_json(self, data: dict, head_only: bool = False):
        """Send a JSON response with CORS headers."""
        body = json.dumps(data).encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", len(body))
        # Note: CORS headers are added by end_headers() override in subclasses
        self.end_headers()
        if not head_only:
            self._response_writer().write(body)

    def send_error_response(self, code: int, message: str):
        """Send an error response with CORS headers."""
        self.send_response(code)
        self.send_header("Content-Type", "text/plain")
        body = message.encode("utf-8")
        self.send_header("Content-Length", len(body))
        # Note: CORS headers are added by end_headers() override in subclasses
        self.end_headers()
        if getattr(self, "command", "") != "HEAD":
            self._response_writer().write(body)

    def send_range_not_satisfiable(self, size: int) -> None:
        """Answer 416 for a ``Range`` header :func:`parse_byte_range` refused."""
        body = b"Requested byte range is not satisfiable"
        self.send_response(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Content-Range", f"bytes */{size}")
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        if getattr(self, "command", "") != "HEAD":
            self._response_writer().write(body)

    def _discard_session_upload(
        self,
        path: Path,
        *,
        failure_message: str,
    ) -> bool:
        """Remove one failed upload or emit an explicit cleanup failure."""
        try:
            path.unlink()
        except OSError:
            logger.exception(failure_message)
            self.send_error_response(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                failure_message,
            )
            return False
        return True

    def do_POST(self) -> None:
        """Route a POST to its endpoint, then drain whatever body it left unread."""
        try:
            if not (
                self.handle_event_post()
                or self.handle_selection_post()
                or self.handle_session_bundle_post()
            ):
                self.send_error_response(
                    HTTPStatus.NOT_FOUND,
                    f"POST not supported for path: {self.path}",
                )
        finally:
            self._discard_unread_request_body()

    def handle_event_post(self) -> bool:
        """
        Handle POST request to /_cellucid/events endpoint.

        This enables Frontend → Python communication for hooks in ALL environments.
        The frontend POSTs events here, and they get routed to the appropriate viewer.

        Returns True if the request was handled, False if path doesn't match.
        """
        if self.path != "/_cellucid/events":
            return False

        try:
            headers = self._request_headers()
            _require_exact_content_type(headers, "application/json")
            content_length = _require_content_length(
                headers,
                maximum=MAX_EVENT_BODY_BYTES,
            )
            body = _read_exact_request_body(
                self._read_request_body,
                content_length,
            )
            try:
                event = json.loads(body.decode("utf-8"))
            except (UnicodeDecodeError, json.JSONDecodeError) as error:
                raise _RequestContractError(
                    HTTPStatus.BAD_REQUEST,
                    "Event body must be valid UTF-8 JSON",
                ) from error
            if not isinstance(event, dict):
                raise _RequestContractError(
                    HTTPStatus.BAD_REQUEST,
                    "Event body must be a JSON object",
                )
            viewer_id = event.get("viewerId")
            if (
                not isinstance(viewer_id, str)
                or not viewer_id
                or viewer_id != viewer_id.strip()
                or any(character.isspace() for character in viewer_id)
            ):
                raise _RequestContractError(
                    HTTPStatus.BAD_REQUEST,
                    "viewerId must be a non-empty string without whitespace",
                )
            event_type = event.get("type")
            if (
                not isinstance(event_type, str)
                or not event_type
                or event_type != event_type.strip()
                or any(character.isspace() for character in event_type)
            ):
                raise _RequestContractError(
                    HTTPStatus.BAD_REQUEST,
                    "type must be a non-empty string without whitespace",
                )
            binary_route = _BINARY_ROUTED_EVENT_PATHS.get(event_type)
            if binary_route is not None:
                raise _RequestContractError(
                    HTTPStatus.BAD_REQUEST,
                    f"{event_type} events are not accepted here; the viewer posts "
                    f"them to {binary_route}",
                )
            supplied_token = event.get("viewerToken")
            viewer_token = supplied_token if isinstance(supplied_token, str) else ""
            routed_event = {key: value for key, value in event.items() if key != "viewerToken"}
        except _RequestContractError as error:
            self.send_error_response(error.status, str(error))
            return True

        try:
            route_event(viewer_id, viewer_token, routed_event)
        except _ViewerNotRegisteredError:
            self.send_error_response(
                HTTPStatus.NOT_FOUND,
                "Viewer is not registered",
            )
        except _ViewerAuthenticationError:
            self.send_error_response(
                HTTPStatus.FORBIDDEN,
                "Invalid viewer credentials",
            )
        except Exception:
            logger.exception("Viewer event callback failed")
            self.send_error_response(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                "Viewer callback failed",
            )
        else:
            self.send_json({"status": "ok", "delivered": True})

        return True

    def handle_selection_post(self) -> bool:
        """
        Handle POST request to /_cellucid/selection endpoint.

        The viewer delivers a selection here rather than as a JSON event, so a
        selection of any size the viewer can make is delivered (the body format
        itself bounds it at every ``uint32`` index once,
        :data:`MAX_SELECTION_BODY_BYTES`), not the event endpoint's 1 MiB:

        - the query names ``viewerId``, ``viewerToken`` and ``source``, each once;
        - the credentials are checked before the body is read;
        - the body (``application/octet-stream``) is the selected cell indices as
          little-endian ``uint32``, strictly ascending, one or more;
        - the server routes one ``selection`` event whose ``cells`` field is
          those bytes, and the viewer's validator decodes them into the
          ``list[int]`` its hooks receive.

        Returns True if the request was handled, False if path doesn't match.
        """
        from urllib.parse import urlsplit

        parsed = urlsplit(self.path)
        if parsed.path != "/_cellucid/selection":
            return False

        try:
            values = _parse_exact_query(
                parsed.query,
                frozenset({"viewerId", "viewerToken", "source"}),
                label="Selection delivery",
            )
            viewer_id = values.get("viewerId")
            source = values.get("source")
            if not _is_exact_route_text(viewer_id) or not _is_exact_route_text(source):
                raise _RequestContractError(
                    HTTPStatus.BAD_REQUEST,
                    "viewerId and source must be non-empty strings without whitespace",
                )
            viewer_token = values.get("viewerToken", "")
        except _RequestContractError as error:
            self.send_error_response(error.status, str(error))
            return True

        try:
            _require_viewer_credentials(viewer_id, viewer_token)
        except _ViewerNotRegisteredError:
            self.send_error_response(
                HTTPStatus.NOT_FOUND,
                "Viewer is not registered",
            )
            return True
        except _ViewerAuthenticationError:
            self.send_error_response(
                HTTPStatus.FORBIDDEN,
                "Invalid viewer credentials",
            )
            return True

        # Deferred like every `_wire` import in the servers: reaching it
        # initializes the `cellucid.jupyter` package. A registered viewer means
        # that package is already loaded.
        from .jupyter._wire import _require_selection_cell_array

        try:
            headers = self._request_headers()
            _require_exact_content_type(headers, "application/octet-stream")
            content_length = _require_content_length(
                headers,
                maximum=MAX_SELECTION_BODY_BYTES,
            )
            body = _read_exact_request_body(
                self._read_request_body,
                content_length,
            )
            try:
                _require_selection_cell_array(body)
            except ValueError as error:
                raise _RequestContractError(HTTPStatus.BAD_REQUEST, str(error)) from error
        except _RequestContractError as error:
            self.send_error_response(error.status, str(error))
            return True

        try:
            route_event(
                viewer_id,
                viewer_token,
                {
                    "type": "selection",
                    "viewerId": viewer_id,
                    "cells": body,
                    "source": source,
                },
            )
        except _ViewerNotRegisteredError:
            self.send_error_response(
                HTTPStatus.NOT_FOUND,
                "Viewer is not registered",
            )
        except _ViewerAuthenticationError:
            self.send_error_response(
                HTTPStatus.FORBIDDEN,
                "Invalid viewer credentials",
            )
        except Exception:
            logger.exception("Viewer selection callback failed")
            self.send_error_response(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                "Viewer callback failed",
            )
        else:
            self.send_json({"status": "ok", "delivered": True})

        return True

    def handle_web_asset_get(self, url_path: str, head_only: bool = False) -> bool:
        """Serve one verified inventory object from the active local generation."""
        if not self.serve_web_ui:
            return False
        if not isinstance(url_path, str):
            raise TypeError("Web request path must be a string")
        if url_path in {"/", "/index.html"}:
            asset_path = "index.html"
        else:
            if (
                not url_path.startswith("/")
                or "\\" in url_path
                or "?" in url_path
                or "#" in url_path
                or any(character.isspace() for character in url_path)
            ):
                return False
            asset_path = url_path[1:]
            if any(part in {"", ".", ".."} for part in asset_path.split("/")):
                return False
        if asset_path == WEB_ASSET_INVENTORY_FILENAME:
            return False

        from .web_cache import read_cached_web_asset

        result = read_cached_web_asset(self.web_cache_dir, asset_path)
        if result is None:
            return False
        data, content_type = result
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        if not head_only:
            self._response_writer().write(data)
        return True

    def handle_session_bundle_post(self) -> bool:
        """
        Handle POST request to /_cellucid/session_bundle endpoint.

        This supports the Jupyter "no-download" workflow:
        - Python sends a `requestSessionBundle` command to the viewer (postMessage)
        - The frontend uploads the raw `.cellucid-session` bytes here
        - The server validates + streams to a temp file (keeps memory bounded)
        - The server routes a `session_bundle` event back to Python (hooks/state)

        Returns True if the request was handled, False if path doesn't match.
        """
        from urllib.parse import urlsplit

        parsed = urlsplit(self.path)
        if parsed.path != "/_cellucid/session_bundle":
            return False

        try:
            values = _parse_exact_query(
                parsed.query,
                frozenset({"viewerId", "viewerToken", "requestId"}),
                label="Session upload",
            )
        except _RequestContractError as error:
            self.send_error_response(error.status, str(error))
            return True
        viewer_id = values.get("viewerId")
        request_id = values.get("requestId")
        if not _is_exact_route_text(viewer_id) or not _is_exact_route_text(request_id):
            self.send_error_response(
                HTTPStatus.BAD_REQUEST,
                "viewerId and requestId must be non-empty strings without whitespace",
            )
            return True
        viewer_token = values.get("viewerToken", "")

        try:
            _require_session_bundle_request(
                viewer_id,
                viewer_token,
                request_id,
            )
        except _SessionBundleRequestNotFoundError:
            self.send_error_response(
                HTTPStatus.NOT_FOUND,
                "No pending session bundle request",
            )
            return True
        except _ViewerAuthenticationError:
            self.send_error_response(
                HTTPStatus.FORBIDDEN,
                "Invalid viewer credentials",
            )
            return True

        try:
            headers = self._request_headers()
            _require_exact_content_type(
                headers,
                "application/octet-stream",
            )
            content_length = _require_content_length(
                headers,
                maximum=_session_bundle_capacity(),
            )
        except _RequestContractError as error:
            self.send_error_response(error.status, str(error))
            return True

        tmp_path: Path | None = None
        bytes_written = 0

        try:
            magic_length = len(SESSION_BUNDLE_MAGIC)
            if content_length < magic_length:
                raise _RequestContractError(
                    HTTPStatus.BAD_REQUEST,
                    _NOT_A_SESSION_FILE_MESSAGE,
                )
            magic = _read_exact_request_body(
                self._read_request_body,
                magic_length,
            )
            if magic != SESSION_BUNDLE_MAGIC:
                raise _RequestContractError(
                    HTTPStatus.BAD_REQUEST,
                    _NOT_A_SESSION_FILE_MESSAGE,
                )

            with tempfile.NamedTemporaryFile(
                delete=False,
                prefix="cellucid-session-",
                suffix=".cellucid-session",
            ) as tmp:
                tmp_path = Path(tmp.name)
                tmp.write(magic)
                bytes_written += len(magic)

                remaining = content_length - magic_length
                while remaining > 0:
                    chunk = self._read_request_body(min(64 * 1024, remaining))
                    if not chunk:
                        raise _RequestContractError(
                            HTTPStatus.BAD_REQUEST,
                            "Session bundle body is shorter than Content-Length",
                        )
                    tmp.write(chunk)
                    bytes_written += len(chunk)
                    remaining -= len(chunk)

            from .session_bundle import CellucidSessionBundle

            bundle = CellucidSessionBundle(tmp_path)
            for _chunk in bundle.iter_chunks(decode=True):
                pass
        except _RequestContractError as error:
            if tmp_path is not None and not self._discard_session_upload(
                tmp_path,
                failure_message="Failed to discard invalid session bundle",
            ):
                return True
            self.send_error_response(error.status, str(error))
            return True
        except (ValueError, zlib.error) as error:
            if tmp_path is not None and not self._discard_session_upload(
                tmp_path,
                failure_message="Failed to discard invalid session bundle",
            ):
                return True
            self.send_error_response(
                HTTPStatus.BAD_REQUEST,
                f"Invalid session bundle: {error}",
            )
            return True
        except OSError:
            logger.exception("Failed to receive a session bundle")
            if tmp_path is not None and not self._discard_session_upload(
                tmp_path,
                failure_message=(
                    "Session bundle receipt failed and its temporary file could not be removed"
                ),
            ):
                return True
            self.send_error_response(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                "Failed to receive session bundle",
            )
            return True

        assert tmp_path is not None
        try:
            _consume_session_bundle_request(
                viewer_id,
                viewer_token,
                request_id,
            )
        except _SessionBundleRequestNotFoundError:
            if not self._discard_session_upload(
                tmp_path,
                failure_message=(
                    "The pending session request disappeared and its temporary "
                    "file could not be removed"
                ),
            ):
                return True
            self.send_error_response(
                HTTPStatus.NOT_FOUND,
                "No pending session bundle request",
            )
            return True
        except _ViewerAuthenticationError:
            if not self._discard_session_upload(
                tmp_path,
                failure_message=(
                    "Session authentication failed and its temporary file could not be removed"
                ),
            ):
                return True
            self.send_error_response(
                HTTPStatus.FORBIDDEN,
                "Invalid viewer credentials",
            )
            return True

        try:
            route_event(
                viewer_id,
                viewer_token,
                {
                    "type": "session_bundle",
                    "viewerId": viewer_id,
                    "requestId": request_id,
                    "status": "ok",
                    "bytes": bytes_written,
                    "path": str(tmp_path),
                },
            )
        except _ViewerNotRegisteredError:
            if not self._discard_session_upload(
                tmp_path,
                failure_message=(
                    "The viewer callback disappeared and its temporary session "
                    "file could not be removed"
                ),
            ):
                return True
            self.send_error_response(
                HTTPStatus.GONE,
                "Viewer callback is no longer registered",
            )
            return True
        except _ViewerAuthenticationError:
            if not self._discard_session_upload(
                tmp_path,
                failure_message=(
                    "Session callback authentication failed and its temporary "
                    "file could not be removed"
                ),
            ):
                return True
            self.send_error_response(
                HTTPStatus.FORBIDDEN,
                "Invalid viewer credentials",
            )
            return True
        except Exception:
            logger.exception("Session bundle callback failed")
            if not self._discard_session_upload(
                tmp_path,
                failure_message=(
                    "The viewer callback failed and its temporary session file could not be removed"
                ),
            ):
                return True
            self.send_error_response(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                "Viewer callback failed",
            )
            return True

        self.send_json({"status": "ok", "bytes": bytes_written})

        return True
