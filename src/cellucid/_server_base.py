"""
Shared base classes for Cellucid HTTP servers.

This module contains the common functionality used by both:
- CellucidServer (serves pre-exported data)
- AnnDataServer (serves AnnData directly)

It's an internal module (prefixed with _) and not part of the public API.
"""

from __future__ import annotations

import json
import logging
import os
import posixpath
import re
import shutil
import socket
import threading
import time
import tempfile
from functools import partial
from http import HTTPStatus
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from typing import Any, Callable

logger = logging.getLogger("cellucid.server")

DEFAULT_PORT = 8765
DEFAULT_HOST = "127.0.0.1"
CELLUCID_WEB_URL = "https://www.cellucid.com"

# Session bundle upload protocol (used for Jupyter "no download" capture).
SESSION_BUNDLE_MAGIC = b"CELLUCID_SESSION\n"
MAX_SESSION_BUNDLE_BYTES = 512 * 1024 * 1024  # 512MB hard cap (zip-bomb + memory guard)

# When the local web UI assets aren't available, the Python server can
# transparently proxy the hosted web assets from `CELLUCID_WEB_URL` so the
# viewer runs on the same origin as the dataset server (avoids mixed-content).
_WEB_PROXY_ROOT_FILES = {
    "/favicon.ico",
    "/site.webmanifest",
    "/browserconfig.xml",
}

_WEB_PROXY_BUILD_META_NAME = "cellucid-web-build-id"
_WEB_PROXY_BUILD_ID_MAX_LEN = 200


def _extract_web_build_id(index_html_bytes: bytes) -> str | None:
    """
    Extract the build/version stamp from the hosted viewer's `index.html`.

    The web app is expected to include:
        <meta name="cellucid-web-build-id" content="...">
    """
    try:
        text = index_html_bytes.decode("utf-8", errors="ignore")
    except Exception:
        return None

    patterns = [
        rf'<meta[^>]*\bname=["\']{re.escape(_WEB_PROXY_BUILD_META_NAME)}["\'][^>]*\bcontent=["\']([^"\']+)["\']',
        rf'<meta[^>]*\bcontent=["\']([^"\']+)["\'][^>]*\bname=["\']{re.escape(_WEB_PROXY_BUILD_META_NAME)}["\']',
    ]
    for pat in patterns:
        m = re.search(pat, text, flags=re.IGNORECASE)
        if not m:
            continue
        build_id = (m.group(1) or "").strip()
        if not build_id:
            continue
        if len(build_id) > _WEB_PROXY_BUILD_ID_MAX_LEN:
            build_id = build_id[:_WEB_PROXY_BUILD_ID_MAX_LEN]
        return build_id

    return None


def _purge_web_proxy_cache_dir(cache_dir: Path) -> None:
    """Best-effort: delete all cached web proxy files under `cache_dir`."""
    if not cache_dir.exists():
        return
    try:
        for child in cache_dir.iterdir():
            try:
                if child.is_dir():
                    shutil.rmtree(child)
                else:
                    child.unlink(missing_ok=True)
            except Exception:
                # Keep going; cache is best-effort.
                pass
    except Exception:
        # Cache is best-effort.
        pass


def _build_web_proxy_error_html(*, reason: str, cache_dir: Path) -> str:
    cache_hint = str(cache_dir)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Cellucid: Viewer UI unavailable</title>
  <style>
    body {{ font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; padding: 24px; }}
    code {{ background: #f3f4f6; padding: 2px 6px; border-radius: 4px; }}
    .box {{ max-width: 900px; border: 1px solid #e5e7eb; border-radius: 10px; padding: 16px 18px; }}
    h1 {{ font-size: 18px; margin: 0 0 10px; }}
    p {{ margin: 8px 0; line-height: 1.4; }}
    ul {{ margin: 8px 0 0 18px; }}
  </style>
</head>
<body>
  <div class="box">
    <h1>Cellucid viewer UI could not be loaded</h1>
    <p><strong>Reason:</strong> {reason}</p>
    <p>Cellucid is running in <code>hosted-asset proxy</code> mode, which downloads the viewer UI from <code>{CELLUCID_WEB_URL}</code> and caches it locally.</p>
    <p><strong>What to do:</strong></p>
    <ul>
      <li>Ensure the runtime can reach <code>{CELLUCID_WEB_URL}</code> over HTTPS (internet / proxy / firewall).</li>
      <li>If you are offline, run once while online to populate the cache, then retry.</li>
      <li>If the cache directory is not writable/persistent, set <code>CELLUCID_WEB_PROXY_CACHE_DIR</code> to a writable path (current: <code>{cache_hint}</code>).</li>
    </ul>
  </div>
</body>
</html>
"""


def _web_proxy_cache_dir() -> Path:
    env = os.getenv("CELLUCID_WEB_PROXY_CACHE_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return Path(tempfile.gettempdir()) / "cellucid-web-cache"


def _normalize_web_path(url_path: str) -> str | None:
    """Normalize and validate a URL path (POSIX) for proxy/cache usage."""
    if not url_path:
        return None
    if not url_path.startswith("/"):
        url_path = "/" + url_path
    if any(ch.isspace() for ch in url_path):
        return None
    if "*" in url_path or "`" in url_path:
        return None

    norm = posixpath.normpath(url_path)
    if not norm.startswith("/"):
        return None
    # Reject traversal.
    parts = [p for p in norm.split("/") if p]
    if any(p == ".." for p in parts):
        return None
    return norm


def _is_proxyable_web_path(url_path: str) -> bool:
    if url_path == "/" or url_path == "/index.html":
        return True
    if url_path.startswith("/assets/"):
        return True
    if url_path in _WEB_PROXY_ROOT_FILES:
        return True
    return False


class _PendingSessionBundleRequest:
    __slots__ = ("expires_at",)

    def __init__(self, expires_at: float):
        self.expires_at = expires_at


_pending_session_bundle_lock = threading.Lock()
_pending_session_bundle_requests: dict[tuple[str, str], _PendingSessionBundleRequest] = {}


def register_session_bundle_request(viewer_id: str, request_id: str, ttl_seconds: float = 60.0) -> None:
    """Register an expected (viewerId, requestId) pair for a session bundle upload."""
    now = time.monotonic()
    expires_at = now + max(1.0, float(ttl_seconds))
    key = (viewer_id, request_id)
    with _pending_session_bundle_lock:
        # Prune expired entries opportunistically.
        expired = [k for k, v in _pending_session_bundle_requests.items() if v.expires_at <= now]
        for k in expired:
            _pending_session_bundle_requests.pop(k, None)
        _pending_session_bundle_requests[key] = _PendingSessionBundleRequest(expires_at=expires_at)


def _consume_session_bundle_request(viewer_id: str, request_id: str) -> bool:
    """Consume a pending request if present and not expired."""
    now = time.monotonic()
    key = (viewer_id, request_id)
    with _pending_session_bundle_lock:
        req = _pending_session_bundle_requests.get(key)
        if req is None:
            return False
        if req.expires_at <= now:
            _pending_session_bundle_requests.pop(key, None)
            return False
        _pending_session_bundle_requests.pop(key, None)
        return True


# =============================================================================
# CLI Output Helpers
# =============================================================================
# These functions provide consistent, professional CLI output across all servers.

def print_step(step: int, total: int, title: str):
    """Print a step header like [1/4] Title..."""
    print(f"\n[{step}/{total}] {title}...")


def print_detail(label: str, value: str):
    """Print an indented detail line."""
    print(f"      {label}: {value}")


def print_success(message: str = "Done"):
    """Print a success indicator with checkmark."""
    print(f"      ✓ {message}")


def print_server_banner(url: str, viewer_url: str):
    """Print the final server ready banner."""
    print(f"\n{'═' * 60}")
    print("  CELLUCID SERVER RUNNING")
    print(f"{'═' * 60}")
    print()
    print(f"  Local URL:    {url}")
    print(f"  Viewer URL:   {viewer_url}")
    print()
    print("  Press Ctrl+C to stop")
    print()
    print(f"{'═' * 60}")


# Global registry for event callbacks (viewer_id -> callback function)
# This allows HTTP event handlers to route events to the correct viewer
_event_callbacks: dict[str, Callable[[dict], None]] = {}


def register_event_callback(viewer_id: str, callback: Callable[[dict], None]):
    """Register a callback to receive events for a viewer."""
    _event_callbacks[viewer_id] = callback


def unregister_event_callback(viewer_id: str):
    """Unregister a viewer's event callback."""
    _event_callbacks.pop(viewer_id, None)


def route_event(viewer_id: str, event: dict) -> bool:
    """
    Route an event to the appropriate viewer.

    Returns True if the event was delivered, False otherwise.
    """
    callback = _event_callbacks.get(viewer_id)
    if callback:
        try:
            callback(event)
            return True
        except Exception as e:
            logger.error(f"Error in event callback for viewer {viewer_id}: {e}")
    return False


def find_free_port(host: str = DEFAULT_HOST, start_port: int = DEFAULT_PORT, max_attempts: int = 100) -> int:
    """
    Find a free port starting from start_port.

    Parameters
    ----------
    host : str
        Host to bind to.
    start_port : int
        Port to start searching from.
    max_attempts : int
        Maximum number of ports to try.

    Returns
    -------
    int
        A free port number.

    Raises
    ------
    RuntimeError
        If no free port is found within max_attempts.
    """
    for port in range(start_port, start_port + max_attempts):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind((host, port))
                return port
        except OSError:
            continue
    raise RuntimeError(f"Could not find a free port starting from {start_port}")


def ensure_port_available(host: str, port: int, quiet: bool = False) -> int:
    """
    Ensure the specified port is available, finding a new one if needed.

    Parameters
    ----------
    host : str
        Host to bind to.
    port : int
        Preferred port.
    quiet : bool
        If True, suppress log messages.

    Returns
    -------
    int
        The available port (may be different from input if port was in use).
    """
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind((host, port))
            return port
    except OSError:
        new_port = find_free_port(host, port + 1)
        if not quiet:
            logger.info(f"Port {port} in use, using {new_port}")
        return new_port


class CORSMixin:
    """
    Mixin class providing CORS headers for HTTP handlers.

    Use in conjunction with SimpleHTTPRequestHandler or similar.
    """

    # Override in subclass: whether to allow caching
    allow_caching: bool = True

    def _get_allowed_origin(self) -> str | None:
        # `BaseHTTPRequestHandler` only populates `.headers` after successfully
        # parsing a request. If parsing fails, `send_error()` may call our
        # overridden `end_headers()` while `.headers` is still unset.
        headers = getattr(self, "headers", None)
        if headers is None:
            return None

        origin = headers.get("Origin")
        if not origin:
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
        self.send_header("Access-Control-Expose-Headers", "Content-Length, Content-Range")
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
            self.wfile.write(body)

    def send_binary(
        self,
        data: bytes,
        content_type: str = "application/octet-stream",
        head_only: bool = False,
        compressed: bool = False,
    ):
        """Send binary data with CORS headers."""
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", len(data))
        if compressed:
            self.send_header("Content-Encoding", "gzip")
        # Note: CORS headers are added by end_headers() override in subclasses
        self.end_headers()
        if not head_only:
            self.wfile.write(data)

    def send_error_response(self, code: int, message: str):
        """Send an error response with CORS headers."""
        self.send_response(code)
        self.send_header("Content-Type", "text/plain")
        body = message.encode("utf-8")
        self.send_header("Content-Length", len(body))
        # Note: CORS headers are added by end_headers() override in subclasses
        self.end_headers()
        if getattr(self, "command", "") != "HEAD":
            self.wfile.write(body)

    def handle_event_post(self) -> bool:
        """
        Handle POST request to /_cellucid/events endpoint.

        This enables Frontend → Python communication for hooks in ALL environments.
        The frontend POSTs events here, and they get routed to the appropriate viewer.

        Returns True if the request was handled, False if path doesn't match.
        """
        from urllib.parse import urlparse, unquote

        parsed = urlparse(self.path)
        path = unquote(parsed.path)

        if path != "/_cellucid/events":
            return False

        try:
            # Read the POST body
            content_length = int(self.headers.get("Content-Length", 0))
            if content_length > 1024 * 1024:  # 1MB limit
                self.send_error_response(413, "Request too large")
                return True

            body = self.rfile.read(content_length)
            event = json.loads(body.decode("utf-8"))

            # Route to the appropriate viewer
            viewer_id = event.get("viewerId")
            if not viewer_id:
                self.send_json({"status": "error", "message": "Missing viewerId"})
                return True

            delivered = route_event(viewer_id, event)
            self.send_json({
                "status": "ok" if delivered else "not_found",
                "delivered": delivered
            })

        except json.JSONDecodeError as e:
            self.send_json({"status": "error", "message": f"Invalid JSON: {e}"})
        except Exception:
            logger.exception("Error handling event POST")
            self.send_error_response(500, "Internal server error")

        return True

    def handle_web_proxy_get(self, url_path: str, head_only: bool = False) -> bool:
        """
        Best-effort proxy for hosted web UI assets (index + /assets/*) from CELLUCID_WEB_URL.

        By serving the viewer UI from the same origin as the dataset server, we
        avoid mixed-content and cross-origin issues in notebook/webview contexts
        (Jupyter/JupyterLab/Colab/VSCode).

        Returns True if the request was handled, False otherwise.
        """
        # Only proxy when explicitly enabled.
        if not getattr(self, "web_proxy", False):
            return False

        normalized = _normalize_web_path(url_path)
        if not normalized:
            return False
        if not _is_proxyable_web_path(normalized):
            return False

        if normalized == "/":
            normalized = "/index.html"

        cache_dir = _web_proxy_cache_dir()
        rel = normalized.lstrip("/")
        cache_path = (cache_dir / rel).resolve()
        try:
            cache_path.relative_to(cache_dir.resolve())
        except Exception:
            return False

        remote_url = f"{CELLUCID_WEB_URL}{normalized}"

        def _serve_bytes(data: bytes, content_type: str) -> bool:
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            if not head_only:
                self.wfile.write(data)
            return True

        def _guess_content_type() -> str:
            content_type = None
            guess_type = getattr(self, "guess_type", None)
            if callable(guess_type):
                try:
                    content_type = guess_type(str(cache_path))
                except Exception:
                    content_type = None
            if not content_type:
                content_type = (
                    "text/html; charset=utf-8"
                    if normalized.endswith(".html")
                    else "application/octet-stream"
                )
            return content_type

        def _send_proxy_error(reason: str) -> bool:
            # For the main page, return an HTML explanation so beginners know what to do.
            if normalized == "/index.html":
                body = _build_web_proxy_error_html(reason=reason, cache_dir=cache_dir).encode("utf-8")
                self.send_response(HTTPStatus.SERVICE_UNAVAILABLE)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                if not head_only:
                    self.wfile.write(body)
                return True

            # For JS/CSS/assets requests, return a short text error (the index.html already explains).
            msg = f"Cellucid web proxy fetch failed for {remote_url}: {reason}\n".encode("utf-8")
            self.send_response(HTTPStatus.SERVICE_UNAVAILABLE)
            self.send_header("Content-Type", "text/plain; charset=utf-8")
            self.send_header("Content-Length", str(len(msg)))
            self.end_headers()
            if not head_only:
                self.wfile.write(msg)
            return True

        # ------------------------------------------------------------------
        # Cache invalidation via build stamp (index.html only)
        # ------------------------------------------------------------------
        cached_index_bytes: bytes | None = None
        cached_build_id: str | None = None
        if normalized == "/index.html":
            try:
                if cache_path.is_file():
                    cached_index_bytes = cache_path.read_bytes()
                    cached_build_id = _extract_web_build_id(cached_index_bytes)
            except Exception:
                cached_index_bytes = None
                cached_build_id = None

        # ------------------------------------------------------------------
        # Serve from cache if present (fast path; may be refreshed below)
        # ------------------------------------------------------------------
        if normalized != "/index.html":
            try:
                if cache_path.is_file():
                    return _serve_bytes(cache_path.read_bytes(), _guess_content_type())
            except Exception:
                # Fall through to refetch.
                pass

        # ------------------------------------------------------------------
        # Fetch from hosted origin (and cache to disk)
        # ------------------------------------------------------------------
        try:
            import urllib.request

            req = urllib.request.Request(
                remote_url,
                headers={"User-Agent": "cellucid-python (web-proxy)"},
            )
            with urllib.request.urlopen(req, timeout=10) as resp:
                status = getattr(resp, "status", 200)
                if int(status) != 200:
                    return _send_proxy_error(f"HTTP {status}")
                data = resp.read()
                content_type = resp.headers.get("Content-Type") or _guess_content_type()

        except Exception:
            # If we have a cached index.html, we can still serve it offline.
            if normalized == "/index.html" and cached_index_bytes is not None:
                return _serve_bytes(cached_index_bytes, "text/html; charset=utf-8")

            logger.debug("Web proxy fetch failed for %s", remote_url, exc_info=True)
            return _send_proxy_error("Network error (no cached copy available)")

        # If this is index.html, compare build IDs and purge cache if the website changed.
        if normalized == "/index.html":
            remote_build_id = _extract_web_build_id(data)
            if remote_build_id and cached_index_bytes is not None and remote_build_id != cached_build_id:
                _purge_web_proxy_cache_dir(cache_dir)
                # Ensure this response also refreshes the cached index.html.
                try:
                    cache_path = (cache_dir / rel).resolve()
                    cache_path.relative_to(cache_dir.resolve())
                except Exception:
                    pass

        # Persist to cache (best-effort).
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_bytes(data)
        except Exception:
            pass

        return _serve_bytes(data, content_type)

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
        from urllib.parse import urlparse, unquote, parse_qs

        parsed = urlparse(self.path)
        path = unquote(parsed.path)
        if path != "/_cellucid/session_bundle":
            return False

        qs = parse_qs(parsed.query or "")
        viewer_id = (qs.get("viewerId") or [None])[0]
        request_id = (qs.get("requestId") or [None])[0]

        if not viewer_id or not request_id:
            self.send_json({"status": "error", "message": "Missing viewerId/requestId"})
            return True

        if not _consume_session_bundle_request(str(viewer_id), str(request_id)):
            self.send_error_response(404, "No pending session bundle request (viewerId/requestId)")
            return True

        content_length_raw = self.headers.get("Content-Length")
        try:
            content_length = int(content_length_raw or "0")
        except ValueError:
            content_length = 0

        if content_length <= 0:
            self.send_error_response(411, "Missing Content-Length")
            return True

        if content_length > MAX_SESSION_BUNDLE_BYTES:
            self.send_error_response(413, f"Session bundle too large (>{MAX_SESSION_BUNDLE_BYTES} bytes)")
            return True

        tmp_path: str | None = None
        bytes_written = 0

        try:
            with tempfile.NamedTemporaryFile(
                delete=False,
                prefix="cellucid-session-",
                suffix=".cellucid-session",
            ) as tmp:
                tmp_path = tmp.name

                # Validate MAGIC header up front to reject non-session uploads quickly.
                magic_len = len(SESSION_BUNDLE_MAGIC)
                if content_length < magic_len:
                    raise ValueError("Session bundle truncated (missing MAGIC header)")

                magic = self.rfile.read(magic_len)
                if magic != SESSION_BUNDLE_MAGIC:
                    raise ValueError("Invalid session bundle MAGIC header")

                tmp.write(magic)
                bytes_written += len(magic)

                remaining = content_length - magic_len
                while remaining > 0:
                    chunk = self.rfile.read(min(64 * 1024, remaining))
                    if not chunk:
                        raise ValueError("Session bundle truncated while reading body")
                    tmp.write(chunk)
                    bytes_written += len(chunk)
                    remaining -= len(chunk)

            # Route an event back to Python so `viewer.get_session_bundle()` can unblock.
            route_event(
                str(viewer_id),
                {
                    "type": "session_bundle",
                    "viewerId": str(viewer_id),
                    "requestId": str(request_id),
                    "status": "ok",
                    "bytes": bytes_written,
                    "path": tmp_path,
                },
            )

            self.send_json({"status": "ok", "bytes": bytes_written})

        except Exception as e:
            logger.exception("Error handling session bundle upload")
            if tmp_path:
                try:
                    Path(tmp_path).unlink(missing_ok=True)  # py3.8+: ok in 3.10+
                except Exception:
                    pass

            # Best-effort: notify Python so callers don't just hang on a timeout.
            try:
                route_event(
                    str(viewer_id),
                    {
                        "type": "session_bundle",
                        "viewerId": str(viewer_id),
                        "requestId": str(request_id),
                        "status": "error",
                        "error": str(e),
                    },
                )
            except Exception:
                pass

            self.send_error_response(500, f"Failed to receive session bundle: {e}")

        return True
