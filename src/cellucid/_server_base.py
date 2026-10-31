"""
Shared building blocks of the Cellucid HTTP server.

The HTTP server class, byte ranges, CORS and ``Host`` checks, the notebook
event and session-bundle routing, the web-build cache paths and the console
banner that :class:`cellucid.anndata_server.AnnDataServer` is built from.

It's an internal module (prefixed with _) and not part of the public API.
"""

from __future__ import annotations

import contextlib
import ipaddress
import json
import logging
import os
import re
import shutil
import socket
import stat
import sys
import threading
import time
from collections.abc import Callable, Sequence
from html.parser import HTMLParser
from http import HTTPStatus
from http.server import HTTPServer as _StandardHTTPServer
from pathlib import Path
from socketserver import ThreadingMixIn
from typing import TYPE_CHECKING, Any, BinaryIO, NoReturn, cast

from ._console import console_print

if TYPE_CHECKING:
    from http.client import HTTPMessage
    from http.server import BaseHTTPRequestHandler, HTTPServer

logger = logging.getLogger("cellucid.server")

DEFAULT_PORT = 8765
DEFAULT_HOST = "127.0.0.1"
CELLUCID_WEB_URL = "https://www.cellucid.com"
# The origins the hosted app is served from, which may read a server's data
# across origins (CORS): the canonical one and the bare apex, nothing else.
HOSTED_APP_ORIGINS = frozenset({CELLUCID_WEB_URL, "https://cellucid.com"})
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

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        # Every connection being answered, so closing the server closes them;
        # one that switched protocols (the viewer channel) with the close its
        # protocol has.
        self._connections: set[socket.socket] = set()
        self._protocol_closes: dict[socket.socket, Callable[[], None]] = {}
        self._connections_lock = threading.Lock()
        super().__init__(*args, **kwargs)

    def process_request(self, request: Any, client_address: Any) -> None:
        with self._connections_lock:
            self._connections.add(request)
        super().process_request(request, client_address)

    def close_with_server(self, connection: socket.socket, close: Callable[[], None]) -> None:
        """Have :meth:`server_close` end ``connection`` with ``close``.

        ``close`` must not block: it starts the protocol's own close, and the
        connection's thread finishes it.
        """
        with self._connections_lock:
            if connection in self._connections:
                self._protocol_closes[connection] = close

    def shutdown_request(self, request: Any) -> None:
        with self._connections_lock:
            self._connections.discard(request)
            self._protocol_closes.pop(request, None)
        super().shutdown_request(request)

    def server_close(self) -> None:
        """Close the listening socket and every connection still open.

        A browser keeps its connections to a host alive between requests, and
        a viewer's channel stays open for as long as its page does. Closing
        only the listening socket left each of those answered by its handler
        thread: a server started next on the same port saw the page read the
        stopped server's dataset documents over the kept connections, and
        each idle thread lived until its timeout. Shutting each connection
        down wakes its thread from a blocked read, and it ends.

        A viewer channel is closed the WebSocket way instead (1001, then the
        orderly close of :meth:`cellucid._websocket.WebSocket.close`), so its
        page reads why the channel ended; its thread ends after the page closes
        its side, or after the bounded drain.
        """
        super().server_close()
        with self._connections_lock:
            connections = list(self._connections)
            protocol_closes = dict(self._protocol_closes)
            self._connections.clear()
            self._protocol_closes.clear()
        for connection in connections:
            close = protocol_closes.get(connection)
            if close is not None:
                close()
                continue
            # A connection the client or its own thread closed already is gone.
            with contextlib.suppress(OSError):
                connection.shutdown(socket.SHUT_RDWR)


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
    """Return this user's default web cache directory, ``<user cache>/cellucid/web``.

    The user cache is ``%LOCALAPPDATA%`` on Windows, ``~/Library/Caches`` on
    macOS, and ``$XDG_CACHE_HOME`` (an absolute path) or ``~/.cache``
    elsewhere: per user, as the R package's
    ``tools::R_user_dir("cellucid", "cache")/web``. A shared temporary folder
    would be one directory for every user of a machine: on a shared Linux host
    the first user's copy (private to them) could not be replaced by anyone
    else, so no other user's viewer would start.
    """
    if sys.platform == "win32":
        local = os.environ.get("LOCALAPPDATA", "")
        base = Path(local) if local and Path(local).is_absolute() else Path.home() / "AppData" / "Local"
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    else:
        xdg = os.environ.get("XDG_CACHE_HOME", "")
        base = Path(xdg) if xdg and Path(xdg).is_absolute() else Path.home() / ".cache"
    return base / "cellucid" / "web"


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
        client is still sending. A refusal -- a foreign ``Host``, a method this
        server does not answer -- is sent before any body is read, so without
        this the client saw a dropped connection instead of the status it was
        sent. Only a body framed by one canonical ``Content-Length`` is
        drained, and for a bounded time.
        """
        self.close_connection = True
        values: list[str] = self._request_headers().get_all("Content-Length") or []
        if (
            len(values) != 1
            or not values[0].isascii()
            or not values[0].isdigit()
            or self._request_headers().get_all("Transfer-Encoding")
        ):
            return
        declared = int(values[0])
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

        # The hosted app, exactly: no other page on a cellucid.com name (the
        # documentation site runs scripts of its own host) may read the data.
        if origin in HOSTED_APP_ORIGINS:
            return origin

        return None

    def add_cors_headers(self):
        """Add CORS headers to the response."""
        allowed_origin = self._get_allowed_origin()
        if allowed_origin:
            self.send_header("Access-Control-Allow-Origin", allowed_origin)
            self.send_header("Vary", "Origin")
        self.send_header("Access-Control-Allow-Methods", "GET, HEAD, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type, Range")
        # Every header a cross-origin reader is expected to act on has to be
        # listed here, because a header that is sent but not exposed reads as
        # absent. `Content-Encoding` was sent and not exposed, so the browser
        # could not tell a gzipped body from a plain one and sized its progress
        # from the compressed length.
        self.send_header(
            "Access-Control-Expose-Headers",
            "Accept-Ranges, Content-Encoding, Content-Length, Content-Range, ETag, Retry-After",
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

    def do_POST(self) -> None:
        """Refuse every POST: the server's routes are reads and the viewer channel."""
        try:
            self.send_error_response(
                HTTPStatus.METHOD_NOT_ALLOWED,
                "This server answers GET and HEAD; the viewer talks to it on its "
                "channel (/_cellucid/ws)",
            )
        finally:
            self._discard_unread_request_body()

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
