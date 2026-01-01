"""
Cellucid Jupyter Integration

Provides seamless integration with Jupyter notebooks for visualizing
cellucid datasets directly in notebook cells.

Usage:
    from cellucid.jupyter import CellucidViewer, show, show_anndata

    # Quick visualization
    viewer = show_anndata(adata)

    # Control the viewer from Python
    viewer.highlight_cells([100, 200, 300], color="#ff0000")
    viewer.set_color_by("cell_type")

    # React to user interactions
    @viewer.on_selection
    def handle_selection(event):
        print(f"Selected {len(event['cells'])} cells")

Communication (bidirectional, works in all environments):
    Python → Frontend:
        - viewer.highlight_cells()
        - viewer.set_color_by()
        - viewer.set_visibility()
        - viewer.clear_highlights()
        - viewer.reset_view()
        - viewer.send_message()

    Frontend → Python:
        - @viewer.on_selection
        - @viewer.on_hover
        - @viewer.on_click
        - @viewer.on_ready
        - @viewer.on_message

    Both directions work in Jupyter, JupyterLab, Google Colab, and VSCode.
    Frontend → Python uses HTTP POST to the local data server.
"""

from __future__ import annotations

import json
import logging
import os
import secrets
import threading
import time
from urllib.parse import urlparse
import weakref
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

if TYPE_CHECKING:
    import anndata

from .server import CellucidServer, DEFAULT_PORT
from ._server_base import (
    CELLUCID_WEB_URL,
    _extract_web_build_id,
    _web_proxy_cache_dir,
    register_event_callback,
    register_session_bundle_request,
    unregister_event_callback,
)

logger = logging.getLogger("cellucid.jupyter")


# =============================================================================
# HOOKS SYSTEM
# =============================================================================

# Type alias for hook callbacks
HookCallback = Callable[[dict[str, Any]], None]


@dataclass
class HookRegistry:
    """
    Manages event hooks for viewer interactions.

    This class provides a clean interface for registering, unregistering,
    and triggering callbacks for various viewer events.

    Supported events:
        - selection: User selected cells (via lasso, click, etc.)
        - hover: User hovering over a cell
        - click: User clicked on a cell
        - ready: Viewer finished initial load
        - message: Raw message (catches all events)

    Example:
        hooks = HookRegistry()

        @hooks.on('selection')
        def my_handler(event):
            print(event['cells'])

        # Or without decorator:
        hooks.register('selection', my_handler)

        # Trigger (called internally when frontend sends message):
        hooks.trigger('selection', {'cells': [1, 2, 3]})
    """

    _callbacks: dict[str, list[HookCallback]] = field(default_factory=dict)

    def register(self, event: str, callback: HookCallback) -> HookCallback:
        """
        Register a callback for an event.

        Args:
            event: Event name ('selection', 'hover', 'click', 'ready', 'message')
            callback: Function to call when event fires. Receives event dict.

        Returns:
            The callback (for decorator usage)
        """
        if event not in self._callbacks:
            self._callbacks[event] = []
        self._callbacks[event].append(callback)
        return callback

    def unregister(self, event: str, callback: HookCallback) -> bool:
        """
        Remove a callback from an event.

        Args:
            event: Event name
            callback: The callback function to remove

        Returns:
            True if callback was found and removed, False otherwise
        """
        if event in self._callbacks:
            try:
                self._callbacks[event].remove(callback)
                return True
            except ValueError:
                pass
        return False

    def clear(self, event: str | None = None):
        """
        Clear callbacks for an event, or all events if event is None.

        Args:
            event: Event name to clear, or None to clear all
        """
        if event is None:
            self._callbacks.clear()
        elif event in self._callbacks:
            self._callbacks[event].clear()

    def trigger(self, event: str, data: dict[str, Any]):
        """
        Trigger all callbacks for an event.

        Args:
            event: Event name to trigger
            data: Event data dict passed to callbacks
        """
        # Always trigger 'message' handlers for any event
        for cb in self._callbacks.get('message', []):
            try:
                cb({'event': event, **data})
            except Exception as e:
                logger.error(f"Error in message hook callback: {e}")

        # Trigger specific event handlers
        for cb in self._callbacks.get(event, []):
            try:
                cb(data)
            except Exception as e:
                logger.error(f"Error in {event} hook callback: {e}")

    def on(self, event: str) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator for registering event callbacks.

        Example:
            @hooks.on('selection')
            def handle_selection(event):
                print(event['cells'])
        """
        def decorator(callback: HookCallback) -> HookCallback:
            self.register(event, callback)
            return callback
        return decorator


@dataclass
class ViewerState:
    """
    Thread-safe (read-mostly) snapshot of the latest viewer → Python events.

    This is intentionally small: it stores *the latest* payload per event type.
    For awaiting new events, use `viewer.wait_for_event(...)`.
    """

    ready: dict[str, Any] | None = None
    selection: dict[str, Any] | None = None
    hover: dict[str, Any] | None = None
    click: dict[str, Any] | None = None

    last_event_type: str | None = None
    last_event: dict[str, Any] | None = None
    last_updated_at: float | None = None

# Track active viewers for cleanup
_active_viewers: weakref.WeakSet = weakref.WeakSet()


def _detect_jupyter_context() -> dict:
    """Detect the Jupyter environment context."""
    context = {
        "in_jupyter": False,
        "notebook_type": None,  # 'jupyter', 'jupyterlab', 'colab', 'vscode'
        "kernel_id": None,
        "can_iframe": True,
        "preferred_display": "iframe",
    }

    try:
        from IPython import get_ipython
        ipython = get_ipython()
        if ipython is None:
            return context

        # Check if we're in a notebook environment
        if hasattr(ipython, "kernel"):
            context["in_jupyter"] = True

            # Try to detect specific environment
            if "google.colab" in str(type(ipython)):
                context["notebook_type"] = "colab"
            elif "VSCODE_PID" in os.environ:
                context["notebook_type"] = "vscode"
            elif "JPY_PARENT_PID" in os.environ:
                context["notebook_type"] = "jupyterlab"
            else:
                context["notebook_type"] = "jupyter"

            # Get kernel ID if available
            if hasattr(ipython.kernel, "session"):
                context["kernel_id"] = ipython.kernel.session.session

    except ImportError:
        pass

    return context


class BaseViewer:
    """
    Base class for Cellucid viewers in Jupyter notebooks.

    Provides common functionality for both CellucidViewer (pre-exported data)
    and AnnDataViewer (direct AnnData visualization).

    Event Hooks:
        Register callbacks to respond to user interactions in the viewer:

        @viewer.on_selection
        def handle(event):
            print(event['cells'])  # list of selected cell indices

        @viewer.on_hover
        def handle(event):
            print(event['cell'])  # hovered cell index

        @viewer.on_click
        def handle(event):
            print(event['cell'])  # clicked cell index

        @viewer.on_ready
        def handle(event):
            print("Viewer loaded!")

        @viewer.on_message
        def handle(event):
            print(event)  # raw message dict
    """

    def __init__(
        self,
        port: int | None = None,
        height: int = 600,
    ):
        """Initialize common viewer properties."""
        self.port = port or self._find_free_port()
        self.height = height

        self._server = None  # Subclass sets the appropriate server type
        self._viewer_id = secrets.token_hex(8)
        self._viewer_token = secrets.token_hex(16)
        self._context = _detect_jupyter_context()
        self._displayed = False
        self._client_server_url_cache: str | None = None
        self._client_server_port_cache: int | None = None

        # Initialize hooks system
        self._hooks = HookRegistry()

        # Latest-state snapshot + event waiting primitives.
        self.state = ViewerState()
        self._state_lock = threading.Lock()
        self._event_cv = threading.Condition(self._state_lock)
        self._event_seq = 0
        self._recent_events: deque[tuple[int, str, dict[str, Any]]] = deque(maxlen=512)

        # Register this viewer for receiving messages
        _register_viewer_for_messages(self)

    # =========================================================================
    # HOOK DECORATORS
    # =========================================================================

    @property
    def on_selection(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a selection handler.

        Called when user selects cells (lasso, shift-click, etc.)

        Event data:
            - cells: list[int] - indices of selected cells
            - source: str - selection method ('lasso', 'knn', 'proximity', 'annotation', ...)

        Example:
            @viewer.on_selection
            def handle(event):
                selected = adata[event['cells']]
                sc.pl.violin(selected, 'gene1')
        """
        return self._hooks.on('selection')

    @property
    def on_hover(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a hover handler.

        Called when user hovers over a cell.

        Event data:
            - cell: int | None - cell index (None if not hovering over cell)
            - position: dict - {x, y, z} world coordinates

        Example:
            @viewer.on_hover
            def handle(event):
                if event['cell'] is not None:
                    print(f"Cell {event['cell']}: {adata.obs.iloc[event['cell']]}")
        """
        return self._hooks.on('hover')

    @property
    def on_click(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a click handler.

        Called when user clicks on a cell.

        Event data:
            - cell: int - clicked cell index
            - button: int - mouse button (0=left, 1=middle, 2=right)
            - shift: bool - shift key held
            - ctrl: bool - ctrl/cmd key held

        Example:
            @viewer.on_click
            def handle(event):
                print(f"Clicked cell {event['cell']}")
        """
        return self._hooks.on('click')

    @property
    def on_ready(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a ready handler.

        Called when viewer finishes loading and is ready for interaction.

        Event data:
            - n_cells: int - number of cells in dataset
            - dimensions: int - embedding dimensionality (1, 2, or 3)

        Example:
            @viewer.on_ready
            def handle(event):
                print(f"Loaded {event['n_cells']} cells")
        """
        return self._hooks.on('ready')

    @property
    def on_message(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a raw message handler.

        Called for ALL messages from the viewer. Use for debugging
        or handling custom message types.

        Event data:
            - event: str - the event type
            - ... other fields depend on event type

        Example:
            @viewer.on_message
            def handle(event):
                print(f"Got message: {event}")
        """
        return self._hooks.on('message')

    # =========================================================================
    # HOOK MANAGEMENT METHODS
    # =========================================================================

    def register_hook(self, event: str, callback: HookCallback) -> HookCallback:
        """
        Register a callback for an event (non-decorator style).

        Args:
            event: Event name ('selection', 'hover', 'click', 'ready', 'message')
            callback: Function to call when event fires

        Returns:
            The callback function

        Example:
            def my_handler(event):
                print(event['cells'])

            viewer.register_hook('selection', my_handler)
        """
        return self._hooks.register(event, callback)

    def unregister_hook(self, event: str, callback: HookCallback) -> bool:
        """
        Remove a callback from an event.

        Args:
            event: Event name
            callback: The callback to remove

        Returns:
            True if found and removed, False otherwise
        """
        return self._hooks.unregister(event, callback)

    def clear_hooks(self, event: str | None = None):
        """
        Clear all callbacks for an event, or all events.

        Args:
            event: Event name to clear, or None to clear all
        """
        self._hooks.clear(event)

    def _handle_frontend_message(self, message: dict):
        """
        Handle a message received from the frontend.

        This is called by the message routing system when the frontend
        sends data back to Python.

        Args:
            message: Message dict from frontend
        """
        msg_type = message.get('type', '')

        # Map frontend message types to hook events
        event_map = {
            'selection': 'selection',
            'hover': 'hover',
            'click': 'click',
            'ready': 'ready',
        }

        event = event_map.get(msg_type, msg_type)

        # Remove internal fields before passing to callbacks
        data = {k: v for k, v in message.items() if k not in ('type', 'viewerId')}

        self._record_event(event, data)
        self._hooks.trigger(event, data)

    def _record_event(self, event: str, data: dict[str, Any]):
        """Update `viewer.state` and notify any `wait_for_event(...)` callers."""
        now = time.monotonic()
        with self._event_cv:
            self._event_seq += 1
            seq = self._event_seq
            self._recent_events.append((seq, event, data))

            # Update "latest" snapshot.
            self.state.last_event_type = event
            self.state.last_event = data
            self.state.last_updated_at = now
            if event == "ready":
                self.state.ready = data
            elif event == "selection":
                self.state.selection = data
            elif event == "hover":
                self.state.hover = data
            elif event == "click":
                self.state.click = data

            self._event_cv.notify_all()

    def wait_for_event(
        self,
        event: str,
        timeout: float | None = 30.0,
        *,
        predicate: Callable[[dict[str, Any]], bool] | None = None,
    ) -> dict[str, Any]:
        """
        Block until the next event of the given type arrives.

        Parameters
        ----------
        event
            Event type (e.g. "ready", "selection", "session_bundle").
        timeout
            Seconds to wait. None means wait forever.
        predicate
            Optional filter on the event payload.
        """
        deadline = None if timeout is None else (time.monotonic() + float(timeout))
        with self._event_cv:
            start_seq = self._event_seq
            while True:
                for seq, ev_type, payload in self._recent_events:
                    if seq <= start_seq:
                        continue
                    if ev_type != event:
                        continue
                    if predicate is not None and not predicate(payload):
                        continue
                    return payload

                if deadline is not None and time.monotonic() >= deadline:
                    raise TimeoutError(f"Timed out waiting for event '{event}'")

                remaining = None if deadline is None else max(0.0, deadline - time.monotonic())
                self._event_cv.wait(timeout=remaining)

    def wait_for_ready(self, timeout: float | None = 30.0) -> dict[str, Any]:
        """Convenience: wait for the viewer's first `ready` event."""
        if self.state.ready is not None:
            return self.state.ready
        return self.wait_for_event("ready", timeout=timeout)

    def get_session_bundle(self, timeout: float | None = 60.0):
        """
        Request the current `.cellucid-session` bundle and return it as an object.

        This is the "no browser download" workflow: Python triggers a bundle export
        in the frontend, the frontend uploads bytes back to the local server, and
        Python receives a handle to the resulting temp file.
        """
        from .session_bundle import CellucidSessionBundle

        if self._context.get("in_jupyter") and not self._displayed:
            # Best-effort: make the viewer visible if the user calls this directly.
            self.display()

        # Ensure the frontend has finished wiring session bundle export.
        # We treat `timeout` as an overall deadline for ready+upload.
        deadline = None if timeout is None else (time.monotonic() + float(timeout))
        if self.state.ready is None:
            ready_timeout = None if deadline is None else max(0.0, deadline - time.monotonic())
            self.wait_for_ready(timeout=ready_timeout)

        remaining = None if deadline is None else max(0.0, deadline - time.monotonic())
        if remaining is not None and remaining <= 0:
            raise TimeoutError("Timed out before requesting the session bundle (viewer not ready)")

        request_id = secrets.token_hex(16)
        register_session_bundle_request(
            self._viewer_id,
            request_id,
            ttl_seconds=(3600.0 if remaining is None else remaining),
        )
        self.send_message({"type": "requestSessionBundle", "requestId": request_id})

        event = self.wait_for_event(
            "session_bundle",
            timeout=remaining,
            predicate=lambda e: e.get("requestId") == request_id,
        )

        if event.get("status") != "ok":
            raise RuntimeError(event.get("error") or "Failed to capture session bundle")

        path = event.get("path")
        if not path:
            raise RuntimeError("Session bundle upload succeeded but no temp path was returned")

        return CellucidSessionBundle(Path(path))

    def apply_session_to_anndata(
        self,
        adata: "Any",
        *,
        inplace: bool = False,
        timeout: float | None = 60.0,
        cleanup_bundle: bool = True,
        **kwargs,
    ):
        """
        Convenience wrapper: capture a session bundle and apply it to an AnnData.

        Notes
        -----
        - If you want to keep the artifact, call `bundle = viewer.get_session_bundle()`
          and `bundle.save(...)` before applying.
        - By default, the temporary bundle file created by the server is deleted
          after applying (best-effort).
        """
        bundle = self.get_session_bundle(timeout=timeout)
        try:
            if "expected_dataset_id" not in kwargs:
                server = self._server
                adapter = getattr(server, "adapter", None)
                if adapter is not None and hasattr(adapter, "get_dataset_identity"):
                    try:
                        ident = adapter.get_dataset_identity()
                        if isinstance(ident, dict):
                            dataset_id = ident.get("id")
                            if isinstance(dataset_id, str) and dataset_id:
                                kwargs["expected_dataset_id"] = dataset_id
                    except Exception:
                        pass
            return bundle.apply_to_anndata(adata, inplace=inplace, **kwargs)
        finally:
            if cleanup_bundle:
                try:
                    bundle.path.unlink(missing_ok=True)
                except Exception:
                    pass

    def debug_connection(self, timeout: float | None = 5.0, *, max_console_events: int = 50) -> dict[str, Any]:
        """
        Return a structured connectivity/debug report for this viewer.

        Includes:
        - server health/info probes
        - Python→Frontend ping/pong roundtrip (postMessage + HTTP events)
        - frontend debug snapshot (URL/origin/userAgent as seen by the iframe)
        - recent frontend console warnings/errors forwarded to Python
        """
        report: dict[str, Any] = {
            "viewer_id": self._viewer_id,
            "viewer_url": self.viewer_url,
            "server_url": self.server_url,
            "displayed": self._displayed,
            "notebook_context": dict(self._context),
            "server_running": bool(self._server and getattr(self._server, "is_running", lambda: True)()),
            "web_ui": {
                "proxy_cache_dir": str(_web_proxy_cache_dir()),
                "proxy_source": CELLUCID_WEB_URL,
            },
            "state": {
                "ready": self.state.ready,
                "last_event_type": self.state.last_event_type,
                "last_updated_at": self.state.last_updated_at,
            },
        }

        # Jupyter proxy support (recommended for HTTPS/remote notebooks).
        try:
            import jupyter_server_proxy  # type: ignore

            report["jupyter_server_proxy"] = {"installed": True, "module": getattr(jupyter_server_proxy, "__name__", None)}
        except Exception as e:
            report["jupyter_server_proxy"] = {"installed": False, "error": str(e)}

        try:
            report["client_server_url"] = self._get_client_server_url().rstrip("/")
        except Exception as e:
            report["client_server_url_error"] = str(e)

        # Web proxy cache introspection (helps debug offline / stale-cache issues).
        try:
            cache_dir = _web_proxy_cache_dir()
            index_path = cache_dir / "index.html"
            cache_info: dict[str, Any] = {
                "cache_dir_exists": cache_dir.exists(),
                "index_html_exists": index_path.is_file(),
            }
            if index_path.is_file():
                data = index_path.read_bytes()
                cache_info["index_html_bytes"] = len(data)
                cache_info["index_html_build_id"] = _extract_web_build_id(data)
                try:
                    cache_info["index_html_mtime"] = index_path.stat().st_mtime
                except Exception:
                    pass
            report["web_ui"]["cache"] = cache_info
        except Exception as e:
            report["web_ui"]["cache_error"] = str(e)

        # Prefetch marker (written by `cellucid.web_cache.ensure_web_ui_cached()`).
        try:
            import json as _json

            cache_dir = _web_proxy_cache_dir()
            meta_path = cache_dir / ".cellucid-web-prefetch.json"
            if meta_path.is_file():
                report["web_ui"]["prefetch"] = _json.loads(meta_path.read_text("utf-8"))
        except Exception as e:
            report["web_ui"]["prefetch_error"] = str(e)

        # ------------------------------------------------------------------
        # Server probes (no browser needed)
        # ------------------------------------------------------------------
        try:
            import urllib.request

            with urllib.request.urlopen(f"{self.server_url}/_cellucid/health", timeout=2) as f:
                report["server_health"] = json.loads(f.read().decode("utf-8"))
        except Exception as e:
            report["server_health_error"] = str(e)

        try:
            import urllib.request

            with urllib.request.urlopen(f"{self.server_url}/_cellucid/info", timeout=2) as f:
                report["server_info"] = json.loads(f.read().decode("utf-8"))
        except Exception as e:
            report["server_info_error"] = str(e)

        datasets: list[dict[str, Any]] | None = None
        try:
            import urllib.request

            with urllib.request.urlopen(f"{self.server_url}/_cellucid/datasets", timeout=2) as f:
                payload = json.loads(f.read().decode("utf-8"))
                raw = payload.get("datasets") if isinstance(payload, dict) else None
                if isinstance(raw, list):
                    datasets = [d for d in raw if isinstance(d, dict)]
                report["server_datasets"] = datasets
        except Exception as e:
            report["server_datasets_error"] = str(e)

        # Best-effort: fetch a dataset_identity.json using the first dataset entry.
        if datasets:
            try:
                import urllib.request

                first = datasets[0]
                rel = first.get("path") if isinstance(first, dict) else None
                if isinstance(rel, str) and rel:
                    base = f"{self.server_url.rstrip('/')}{rel}"
                    if not base.endswith("/"):
                        base += "/"
                    url = f"{base}dataset_identity.json"
                    with urllib.request.urlopen(url, timeout=2) as f:
                        report["dataset_identity_url"] = url
                        report["dataset_identity"] = json.loads(f.read().decode("utf-8"))
            except Exception as e:
                report["dataset_identity_error"] = str(e)

        # Best-effort: verify that the hosted-asset proxy is serving an index.html.
        try:
            import urllib.request

            req = urllib.request.Request(
                f"{self.server_url}/index.html",
                headers={"User-Agent": "cellucid-python (debug_connection)"},
            )
            with urllib.request.urlopen(req, timeout=2) as f:
                data = f.read()
                report["viewer_index_probe"] = {
                    "bytes": len(data),
                    "content_type": f.headers.get("Content-Type"),
                    "build_id": _extract_web_build_id(data),
                }
        except Exception as e:
            report["viewer_index_probe_error"] = str(e)

        # ------------------------------------------------------------------
        # Recent event summary (helps debug “stuck” hooks)
        # ------------------------------------------------------------------
        try:
            with self._event_cv:
                by_type: dict[str, int] = {}
                for _seq, ev_type, _payload in self._recent_events:
                    by_type[ev_type] = by_type.get(ev_type, 0) + 1
                report["recent_events"] = {
                    "count": len(self._recent_events),
                    "by_type": by_type,
                }
        except Exception as e:
            report["recent_events_error"] = str(e)

        # ------------------------------------------------------------------
        # Frontend roundtrip probe (requires the viewer iframe to be alive)
        # ------------------------------------------------------------------
        if not self._displayed:
            report["frontend_roundtrip"] = {"ok": False, "error": "Viewer not displayed (call viewer.display())"}
            report["frontend_debug_snapshot"] = {"ok": False, "error": "Viewer not displayed (call viewer.display())"}
        else:
            req_id = secrets.token_hex(8)
            try:
                self.send_message({"type": "ping", "requestId": req_id})
                pong = self.wait_for_event(
                    "pong",
                    timeout=timeout,
                    predicate=lambda e: e.get("requestId") == req_id,
                )
                report["frontend_roundtrip"] = {"ok": True, "pong": pong}
            except Exception as e:
                report["frontend_roundtrip"] = {"ok": False, "error": str(e)}

            snap_id = secrets.token_hex(8)
            try:
                self.send_message({"type": "debug_snapshot", "requestId": snap_id})
                snap = self.wait_for_event(
                    "debug_snapshot",
                    timeout=timeout,
                    predicate=lambda e: e.get("requestId") == snap_id,
                )
                report["frontend_debug_snapshot"] = {"ok": True, "snapshot": snap}
            except Exception as e:
                report["frontend_debug_snapshot"] = {"ok": False, "error": str(e)}

        # ------------------------------------------------------------------
        # Recent frontend console warnings/errors (best-effort)
        # ------------------------------------------------------------------
        console_events: list[dict[str, Any]] = []
        with self._event_cv:
            for _seq, ev_type, payload in reversed(self._recent_events):
                if ev_type != "console":
                    continue
                if isinstance(payload, dict):
                    console_events.append(payload)
                if len(console_events) >= max_console_events:
                    break
        console_events.reverse()
        report["frontend_console"] = console_events

        return report

    # =========================================================================
    # WEB UI CACHE UTILITIES
    # =========================================================================

    def clear_web_cache(self) -> Path:
        """
        Clear the hosted-asset web UI cache.

        This forces the next viewer/server run to re-download the web UI assets.
        """
        from .web_cache import clear_web_cache

        return clear_web_cache()

    def ensure_web_ui_cached(self, *, force: bool = False, show_progress: bool = True):
        """
        Best-effort: prefetch the viewer UI assets into the local cache.

        This makes first-load UX clearer (progress bar) and reduces later
        surprises from lazily-loaded assets when offline.
        """
        from .web_cache import ensure_web_ui_cached

        return ensure_web_ui_cached(force=force, show_progress=show_progress)

    # =========================================================================
    # DISPLAY & SERVER
    # =========================================================================

    def _find_free_port(self, start: int = DEFAULT_PORT) -> int:
        """Find an available port."""
        import socket
        for port in range(start, start + 100):
            try:
                with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                    s.bind(("127.0.0.1", port))
                    return port
            except OSError:
                continue
        raise RuntimeError(f"Could not find a free port starting from {start}")

    @property
    def server_url(self) -> str:
        """Get the data server URL."""
        if self._server is not None and getattr(self._server, "url", None):
            return str(self._server.url)
        return f"http://127.0.0.1:{self.port}"

    def _get_client_server_url(self) -> str:
        """
        Get the URL the browser should use to reach the data server.

        Most notebook environments can reach the server directly at `server_url`
        (loopback). Some environments (notably Google Colab) run the kernel on a
        remote VM; in that case, the browser must use Colab's HTTPS port proxy.
        """
        override = os.getenv("CELLUCID_CLIENT_SERVER_URL")
        if override:
            return override.rstrip("/")

        try:
            parsed = urlparse(self.server_url)
            port = int(parsed.port or self.port)
        except Exception:
            port = int(self.port)

        if self._context.get("in_jupyter") and self._context.get("notebook_type") == "colab":
            if self._client_server_port_cache == port and self._client_server_url_cache:
                return self._client_server_url_cache
            try:
                from google.colab.output import eval_js  # type: ignore

                proxy_url = eval_js(f"google.colab.kernel.proxyPort({port})")
                if isinstance(proxy_url, str) and proxy_url:
                    proxy_url = proxy_url.rstrip("/")
                    self._client_server_port_cache = port
                    self._client_server_url_cache = proxy_url
                    return proxy_url
            except Exception:
                # Fall back to direct loopback URL.
                pass

        return self.server_url.rstrip("/")

    @property
    def viewer_origin(self) -> str:
        """Origin (scheme+host+port) for postMessage targetOrigin."""
        parsed = urlparse(self.viewer_url)
        return f"{parsed.scheme}://{parsed.netloc}"

    @property
    def viewer_url(self) -> str:
        """Get the full viewer URL. Subclasses can override to add extra params."""
        from urllib.parse import urlencode

        # We always serve the UI from the same server that serves the dataset
        # (hosted-asset proxy), avoiding mixed-content.
        base = self._get_client_server_url().rstrip("/")
        query = urlencode(
            {
                "jupyter": "true",
                "viewerId": self._viewer_id,
                "viewerToken": self._viewer_token,
            }
        )
        return f"{base}/?{query}"

    def _get_pre_display_html(self) -> str | None:
        """Override to add HTML before the viewer (e.g., warnings)."""
        return None

    def display(self):
        """Display the viewer in a notebook cell."""
        if not self._context["in_jupyter"]:
            print(f"Not in Jupyter environment. Open manually: {self.viewer_url}")
            return

        from IPython.display import display, HTML

        # Ensure the viewer UI assets are available locally (one-time download).
        # This prints a progress bar (tqdm) on first run, which is important UX
        # for beginners in notebook contexts.
        try:
            summary = self.ensure_web_ui_cached(force=False, show_progress=True)
            if getattr(summary, "downloaded_files", 0) or getattr(summary, "errors", None):
                print(
                    f"Cellucid web UI cache: {summary.downloaded_files} file(s), "
                    f"{summary.downloaded_bytes} byte(s) downloaded to {summary.cache_dir}"
                )
                for err in (summary.errors or [])[:5]:
                    print(f"  - {err}")
        except Exception as e:
            # The iframe will still show a more detailed error page if the UI
            # can't be fetched and no cached copy exists; keep notebook output clean.
            logger.warning("Failed to prefetch Cellucid web UI assets: %s", e)

        # Show any pre-display HTML (e.g., warnings)
        pre_html = self._get_pre_display_html()
        if pre_html:
            display(HTML(pre_html))

        # Use HTML with wrapper div for all notebook types
        # This allows us to reference the iframe for postMessage communication
        html = self._generate_viewer_html()
        display(HTML(html))
        self._displayed = True

    def _generate_viewer_html(self) -> str:
        """Generate HTML for embedding the viewer with message passing support."""
        # Important: the initial iframe `src` is set by the script below.
        # This allows us to choose the best URL per notebook frontend:
        # - direct loopback in local classic/JupyterLab
        # - Jupyter Server Proxy when the notebook runs on HTTPS/remote
        # - Colab's HTTPS port proxy (computed on Python side)
        direct_src = self.viewer_url
        port = int(self.port)

        return f"""
        <div id="cellucid-viewer-{self._viewer_id}" style="width:100%; height:{self.height}px;">
            <iframe
                id="cellucid-iframe-{self._viewer_id}"
                src="about:blank"
                width="100%"
                height="100%"
                frameborder="0"
                allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope"
                allowfullscreen>
            </iframe>
        </div>
        <script>
        (function() {{
            var viewerId = '{self._viewer_id}';
            var viewerToken = '{self._viewer_token}';
            var directSrc = {json.dumps(direct_src)};
            var port = {port};

            // We use '*' because notebook frontends can proxy/transform iframe
            // origins (e.g. Colab), and we authenticate messages with viewerToken.
            var targetOrigin = '*';

            function getNotebookBaseUrl() {{
                try {{
                    var baseUrl =
                        (document.body && document.body.dataset && document.body.dataset.baseUrl) ||
                        (document.body && document.body.getAttribute && document.body.getAttribute('data-base-url')) ||
                        '/';
                    if (!baseUrl) baseUrl = '/';
                    if (!baseUrl.endsWith('/')) baseUrl = baseUrl + '/';
                    return baseUrl;
                }} catch (e) {{
                    return '/';
                }}
            }}

            function buildProxySrc() {{
                // Requires jupyter-server-proxy (common in JupyterHub and many local installs).
                var baseUrl = getNotebookBaseUrl();
                if (!(window.location.protocol === 'http:' || window.location.protocol === 'https:')) {{
                    return null;
                }}
                if (!baseUrl.startsWith('/')) {{
                    return null;
                }}
                var q = '';
                try {{
                    var idx = directSrc.indexOf('?');
                    q = idx >= 0 ? directSrc.slice(idx + 1) : '';
                }} catch (e) {{
                    q = '';
                }}
                var root = window.location.origin + baseUrl + 'proxy/' + port + '/';
                return q ? (root + '?' + q) : root;
            }}

            async function probeHealth(url) {{
                try {{
                    var u = new URL(url);
                    u.search = '';
                    u.hash = '';
                    if (!u.pathname.endsWith('/')) u.pathname = u.pathname + '/';
                    u.pathname = u.pathname + '_cellucid/health';
                    var res = await fetch(u.toString(), {{ cache: 'no-store' }});
                    return !!(res && res.ok);
                }} catch (e) {{
                    return false;
                }}
            }}

            // Set up message passing infrastructure
            window.cellucidViewers = window.cellucidViewers || {{}};
            window.cellucidViewers[viewerId] = {{
                sendMessage: function(msg) {{
                    var iframe = document.getElementById('cellucid-iframe-' + viewerId);
                    if (iframe && iframe.contentWindow) {{
                        iframe.contentWindow.postMessage(
                            Object.assign({{}}, msg, {{ viewerId: viewerId, viewerToken: viewerToken }}),
                            targetOrigin
                        );
                    }}
                }}
            }};

            // Choose an iframe URL that avoids HTTPS→HTTP mixed-content blocking.
            (async function() {{
                var iframe = document.getElementById('cellucid-iframe-' + viewerId);
                if (!iframe) return;

                // If Python already provided an HTTPS-capable URL (e.g. Colab proxy),
                // use it directly.
                var isHttp = (typeof directSrc === 'string') && directSrc.startsWith('http://');
                var notebookIsHttps = (window.location && window.location.protocol === 'https:');
                var directIsLoopback = false;
                try {{
                    directIsLoopback = /^http:\\/\\/(127\\.0\\.0\\.1|localhost)(:|\\/)/.test(directSrc);
                }} catch (e) {{
                    directIsLoopback = false;
                }}
                var hostIsLoopback = false;
                try {{
                    var hn = window.location && window.location.hostname;
                    hostIsLoopback = (hn === '127.0.0.1' || hn === 'localhost');
                }} catch (e) {{
                    hostIsLoopback = false;
                }}
                // Prefer Jupyter Server Proxy when direct loopback is unlikely to work:
                // - HTTPS notebook blocks HTTP loopback (mixed content)
                // - Remote notebooks cannot reach kernel loopback directly
                var shouldPreferProxy = directIsLoopback && (notebookIsHttps || !hostIsLoopback);

                var proxySrc = buildProxySrc();
                if (proxySrc && (shouldPreferProxy || !(window.location && window.location.protocol === 'file:'))) {{
                    var ok = await probeHealth(proxySrc);
                    if (ok) {{
                        iframe.src = proxySrc;
                        return;
                    }}
                    if (shouldPreferProxy) {{
                        iframe.srcdoc = [
                          '<div style="font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; padding: 18px;">',
                          '<h3 style="margin: 0 0 8px; font-size: 16px;">Cellucid: notebook proxy required</h3>',
                          '<p style="margin: 0 0 10px; line-height: 1.35;">This notebook is served from a secure/remote origin, so the viewer cannot load an HTTP loopback server directly.</p>',
                          '<p style="margin: 0 0 10px; line-height: 1.35;"><strong>Fix:</strong> install/enable <code>jupyter-server-proxy</code> (recommended) or set <code>CELLUCID_CLIENT_SERVER_URL</code> to a browser-reachable HTTPS URL for the Cellucid server.</p>',
                          '<p style="margin: 0; line-height: 1.35;">In Python, run <code>viewer.debug_connection()</code> for a detailed report.</p>',
                          '</div>'
                        ].join('');
                        return;
                    }}
                }}

                iframe.src = directSrc;
            }})();

            // Note: Frontend → Python communication uses HTTP POST to /_cellucid/events.
            // The postMessage listener is not used for event routing.
        }})();
        </script>
        """

    def send_message(self, message: dict):
        """
        Send a message to the viewer iframe.

        This is the low-level API for sending commands to the frontend.
        Messages are sent via postMessage to the embedded viewer iframe.

        Args:
            message: Dict with 'type' key and message-specific data

        Note:
            The viewer must be displayed first (via display()).

        Example:
            viewer.send_message({'type': 'highlight', 'cells': [1,2,3], 'color': '#ff0000'})
        """
        if not self._displayed:
            logger.warning("Viewer not yet displayed. Call display() first.")
            return

        from IPython.display import display, Javascript

        js = f"""
        (function() {{
            var viewer = window.cellucidViewers && window.cellucidViewers['{self._viewer_id}'];
            if (viewer) {{
                viewer.sendMessage({json.dumps(message)});
            }} else {{
                console.warn('Cellucid viewer not found: {self._viewer_id}');
            }}
        }})();
        """
        display(Javascript(js))

    # =========================================================================
    # PYTHON → FRONTEND CONVENIENCE METHODS
    # =========================================================================
    # These methods send commands TO the viewer (opposite direction from hooks).
    # They work in all environments (Jupyter, JupyterLab, Colab, VSCode).

    def highlight_cells(self, cell_indices: list[int], color: str = "#ff0000"):
        """
        Highlight specific cells in the viewer.

        Args:
            cell_indices: List of cell indices to highlight
            color: Hex color string (default: red)

        Example:
            viewer.highlight_cells([100, 200, 300], color="#00ff00")
        """
        self.send_message({"type": "highlight", "cells": cell_indices, "color": color})

    def clear_highlights(self):
        """Clear all cell highlights in the viewer."""
        self.send_message({"type": "clearHighlights"})

    def set_color_by(self, field: str):
        """
        Set the coloring field in the viewer.

        Args:
            field: Name of the obs column to color by

        Example:
            viewer.set_color_by("cell_type")
        """
        self.send_message({"type": "setColorBy", "field": field})

    def set_visibility(self, cell_indices: list[int] | None = None, visible: bool = True):
        """
        Set visibility of specific cells.

        Args:
            cell_indices: List of cell indices (None for all cells)
            visible: Whether cells should be visible

        Example:
            viewer.set_visibility([0, 1, 2], visible=False)  # Hide cells
        """
        self.send_message({"type": "setVisibility", "cells": cell_indices, "visible": visible})

    def reset_view(self):
        """Reset the camera to the default view."""
        self.send_message({"type": "resetCamera"})

    # =========================================================================
    # LIFECYCLE
    # =========================================================================

    def stop(self):
        """
        Stop the data server and cleanup all resources.

        This method:
        1. Stops the underlying data server (which closes the adapter for AnnData)
        2. Removes the viewer from the active viewers set
        3. Clears all hooks
        4. Logs the cleanup

        Safe to call multiple times.
        """
        # Best-effort: freeze the frontend before the server disappears so the
        # notebook output stays visually reproducible.
        if self._displayed:
            try:
                self.send_message({"type": "freeze"})
            except Exception:
                pass

        if self._server:
            try:
                self._server.stop()
            except Exception as e:
                logger.warning(f"Error stopping server: {e}")
            self._server = None

        # Clear hooks
        self._hooks.clear()

        # Unregister from message routing
        _unregister_viewer_for_messages(self)

        # Remove from active viewers set to prevent double cleanup
        _active_viewers.discard(self)

        logger.info(f"{self.__class__.__name__} stopped")

    def __del__(self):
        """Cleanup on garbage collection."""
        try:
            self.stop()
        except Exception:
            # Ignore errors during GC - object may be partially destroyed
            pass

    def _repr_html_(self) -> str:
        """HTML representation for Jupyter display."""
        if self._context["in_jupyter"]:
            return f'<a href="{self.viewer_url}" target="_blank">Open Cellucid Viewer</a>'
        return repr(self)


# =============================================================================
# MESSAGE ROUTING (Frontend → Python)
# =============================================================================
#
# Bidirectional communication works in ALL environments via HTTP POST:
#
# 1. The viewer iframe POSTs events to the local data server
# 2. The data server (running on localhost) receives the POST at /_cellucid/events
# 3. The event is routed to the appropriate viewer via _handle_frontend_message
# 4. Hooks fire (on_selection, on_hover, on_click, etc.)
#
# This works because:
# - The data server has CORS headers allowing requests from cellucid.com
# - HTTP POST works universally (no special Jupyter/Colab APIs needed)
# - The viewer already knows the server URL (passed in URL params)
#
# =============================================================================

def _register_viewer_for_messages(viewer: BaseViewer):
    """
    Register a viewer to receive messages from frontend.

    This sets up the HTTP-based event routing which works in ALL environments.
    The frontend POSTs events to /_cellucid/events on the data server,
    which routes them to the appropriate viewer.
    """
    viewer_ref = weakref.ref(viewer)

    def _deliver(event: dict) -> None:
        resolved = viewer_ref()
        if resolved is None:
            return
        resolved._handle_frontend_message(event)

    # Register HTTP event callback - this works in ALL environments.
    # The data server routes POSTed events here.
    register_event_callback(viewer._viewer_id, _deliver)
    logger.debug(f"Registered HTTP event callback for viewer {viewer._viewer_id}")


def _unregister_viewer_for_messages(viewer: BaseViewer):
    """Unregister a viewer from message routing."""
    unregister_event_callback(viewer._viewer_id)


class CellucidViewer(BaseViewer):
    """
    Interactive cellucid viewer for Jupyter notebooks (pre-exported data).

    Embeds the cellucid web viewer in a notebook cell, connected to a
    local data server. Supports bidirectional communication via hooks.

    Example:
        >>> viewer = CellucidViewer("/path/to/dataset")
        >>> viewer.display()
        >>>
        >>> @viewer.on_selection
        ... def handle_selection(event):
        ...     print(f"Selected {len(event['cells'])} cells")
        >>>
        >>> # Low-level message API still available:
        >>> viewer.send_message({'type': 'highlight', 'cells': [1,2,3]})
    """

    def __init__(
        self,
        data_dir: str | Path,
        port: int | None = None,
        height: int = 600,
        auto_open: bool = True,
    ):
        """
        Initialize the viewer.

        Args:
            data_dir: Path to the cellucid dataset directory
            port: Port for the data server (auto-selected if None)
            height: Height of the embedded viewer in pixels
            auto_open: Automatically display when created
        """
        super().__init__(port=port, height=height)

        self.data_dir = Path(data_dir).resolve()

        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {self.data_dir}")

        # Start the server
        self._start_server()

        # Register for cleanup
        _active_viewers.add(self)

        if auto_open and self._context["in_jupyter"]:
            self.display()

    def _start_server(self):
        """Start the background data server."""
        self._server = CellucidServer(
            data_dir=self.data_dir,
            port=self.port,
            host="127.0.0.1",
            open_browser=False,
            quiet=True,
        )
        self._server.start_background()
        try:
            self.port = int(getattr(self._server, "port", self.port))
        except Exception:
            pass
        logger.info(f"Started cellucid server at {self._server.url}")

    def __repr__(self) -> str:
        status = "running" if self._server and self._server.is_running() else "stopped"
        return f"CellucidViewer('{self.data_dir}', port={self.port}, status={status})"


def show(
    data_dir: str | Path,
    height: int = 600,
) -> CellucidViewer:
    """
    Quick function to display a cellucid dataset in a notebook.

    Args:
        data_dir: Path to the dataset directory
        height: Height of the viewer in pixels

    Returns:
        CellucidViewer instance for interaction via hooks.

    Example:
        >>> from cellucid.jupyter import show
        >>> viewer = show("/path/to/my_dataset")
        >>>
        >>> @viewer.on_selection
        ... def handle(event):
        ...     print(f"Selected {len(event['cells'])} cells")
    """
    return CellucidViewer(
        data_dir=data_dir,
        height=height,
        auto_open=True,
    )


def cleanup_all():
    """Stop all active viewers and their servers."""
    for viewer in list(_active_viewers):
        try:
            viewer.stop()
        except Exception:
            pass


# Register cleanup on interpreter exit
import atexit
atexit.register(cleanup_all)


# =============================================================================
# ANNDATA DIRECT VISUALIZATION
# =============================================================================


class AnnDataViewer(BaseViewer):
    """
    Interactive viewer for AnnData objects in Jupyter notebooks.

    This viewer serves AnnData directly without requiring prepare.
    It's more convenient for interactive exploration but slower than using
    pre-exported data. Supports bidirectional communication via hooks.

    Supports:
    - In-memory AnnData objects
    - h5ad files (HDF5-based, with lazy loading via backed mode)
    - zarr stores (directory-based, inherently lazy-loaded)

    Example:
        >>> viewer = AnnDataViewer(adata)
        >>> viewer.display()
        >>>
        >>> @viewer.on_selection
        ... def analyze_selection(event):
        ...     subset = adata[event['cells']]
        ...     sc.pl.violin(subset, ['gene1', 'gene2'])
        >>>
        >>> # From h5ad file with lazy loading
        >>> viewer = AnnDataViewer("/path/to/data.h5ad")
    """

    def __init__(
        self,
        data: "str | Path | anndata.AnnData",
        port: int | None = None,
        height: int = 600,
        auto_open: bool = True,
        **adapter_kwargs,
    ):
        """
        Initialize the AnnData viewer.

        Args:
            data: AnnData object or path to h5ad file or zarr directory.
            port: Port for the data server (auto-selected if None).
            height: Height of the embedded viewer in pixels.
            auto_open: Automatically display when created.
            **adapter_kwargs: Additional arguments passed to AnnDataAdapter.
        """
        super().__init__(port=port, height=height)

        self.data = data

        # Start the AnnData server
        self._start_server(**adapter_kwargs)

        # Register for cleanup
        _active_viewers.add(self)

        if auto_open and self._context["in_jupyter"]:
            self.display()

    def _start_server(self, **adapter_kwargs):
        """Start the AnnData server."""
        from .anndata_server import AnnDataServer

        self._server = AnnDataServer(
            data=self.data,
            port=self.port,
            host="127.0.0.1",
            open_browser=False,
            quiet=True,
            **adapter_kwargs,
        )
        self._server.start_background()
        try:
            self.port = int(getattr(self._server, "port", self.port))
        except Exception:
            pass
        logger.info(f"Started AnnData server at {self._server.url}")

    @property
    def viewer_url(self) -> str:
        """Get the full viewer URL with anndata flag."""
        from urllib.parse import urlencode

        base = self._get_client_server_url().rstrip("/")
        query = urlencode(
            {
                "jupyter": "true",
                "viewerId": self._viewer_id,
                "viewerToken": self._viewer_token,
                "anndata": "true",
            }
        )
        return f"{base}/?{query}"

    def _get_pre_display_html(self) -> str | None:
        """Show warning about AnnData mode being slower."""
        return """
        <div style="background: #fff3cd; border: 1px solid #ffc107; border-radius: 4px;
                    padding: 10px; margin-bottom: 10px; color: #856404;">
            <strong>Note:</strong> Loading data directly from AnnData. This is slower than
            using <code>prepare</code>. For production use, consider exporting
            your data first.
        </div>
        """

    def __repr__(self) -> str:
        status = "running" if self._server and self._server.is_running() else "stopped"
        return f"AnnDataViewer(port={self.port}, status={status})"


def show_anndata(
    data: "str | Path | anndata.AnnData",
    height: int = 600,
    **kwargs,
) -> AnnDataViewer:
    """
    Quickly display an AnnData object, h5ad file, or zarr store in a notebook.

    This is the easiest way to visualize AnnData directly without running
    prepare first. However, it's slower than using pre-exported data.

    Args:
        data: AnnData object, path to h5ad file, or path to zarr directory.
        height: Height of the viewer in pixels.
        **kwargs: Additional arguments passed to AnnDataAdapter
            - latent_key: Key in obsm for latent space (auto-detected)
            - gene_id_column: Column in var for gene IDs ("index" by default)
            - normalize_embeddings: Normalize to [-1,1] (default True)
            - dataset_name: Human-readable name

    Returns:
        AnnDataViewer instance for interaction via hooks.

    Supported formats:
        - In-memory AnnData objects
        - .h5ad files (HDF5-based, lazy loading via backed mode)
        - .zarr directories (directory-based, inherently lazy-loaded)

    Example:
        >>> from cellucid import show_anndata
        >>> viewer = show_anndata(adata)
        >>>
        >>> @viewer.on_selection
        ... def handle(event):
        ...     subset = adata[event['cells']]
        ...     sc.pl.violin(subset, 'gene1')

        >>> # With custom options
        >>> viewer = show_anndata(adata, latent_key="X_pca", height=800)

    Note:
        For production use or sharing, consider using prepare
        to create optimized binary files, then use show() to display them.
    """
    return AnnDataViewer(
        data=data,
        height=height,
        auto_open=True,
        **kwargs,
    )
