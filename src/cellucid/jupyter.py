"""
Cellucid Jupyter Integration

Provides seamless integration with Jupyter notebooks for visualizing
cellucid datasets directly in notebook cells.

Usage:
    from cellucid.jupyter import CellucidViewer, show, show_anndata

    # Quick visualization
    viewer = show_anndata(
        adata,
        dataset_name="Example",
        dataset_id="example",
    )

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

import atexit
import html
import json
import logging
import math
import os
import secrets
import threading
import time
import weakref
from collections import deque
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol
from urllib.parse import urlparse

if TYPE_CHECKING:
    import anndata

from ._server_base import (
    CELLUCID_WEB_URL,
    _extract_web_build_id,
    _web_cache_dir,
    cancel_session_bundle_request,
    register_event_callback,
    register_session_bundle_request,
    require_server_port,
    unregister_event_callback,
)
from .server import CellucidServer

logger = logging.getLogger("cellucid.jupyter")


class _ViewerServer(Protocol):
    """Exact lifecycle surface shared by notebook-owned data servers."""

    port: int

    @property
    def url(self) -> str: ...

    def start_background(self) -> None: ...

    def stop(self) -> None: ...

    def is_running(self) -> bool: ...


def _require_client_server_url(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("client_server_url must be a string")
    if not value or value != value.strip() or any(character.isspace() for character in value):
        raise ValueError("client_server_url must be a non-empty URL without surrounding whitespace")
    if value.endswith("/"):
        raise ValueError("client_server_url must not end with '/'")
    parsed = urlparse(value)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        raise ValueError("client_server_url must be an absolute HTTP or HTTPS URL")
    try:
        if parsed.hostname is None:
            raise ValueError("client_server_url must contain a hostname")
        _ = parsed.port
    except ValueError as error:
        raise ValueError("client_server_url contains an invalid host or port") from error
    if parsed.username is not None or parsed.password is not None:
        raise ValueError("client_server_url must not contain credentials")
    if parsed.query or parsed.fragment:
        raise ValueError("client_server_url must not contain a query or fragment")
    return value


def _require_exact_message_text(
    value: object,
    *,
    label: str,
    allow_whitespace: bool = False,
) -> str:
    if (
        type(value) is not str
        or not value
        or value != value.strip()
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
        or (not allow_whitespace and any(character.isspace() for character in value))
    ):
        raise TypeError(f"{label} must be exact non-empty text")
    return value


def _require_exact_json_value(
    value: object,
    *,
    label: str,
    ancestors: set[int] | None = None,
) -> None:
    if value is None or type(value) in {str, bool, int}:
        return
    if type(value) is float:
        if not math.isfinite(value):
            raise ValueError(f"{label} numbers must be finite JSON numbers")
        return
    if ancestors is None:
        ancestors = set()
    identity = id(value)
    if identity in ancestors:
        raise ValueError(f"{label} must not contain JSON cycles")
    ancestors.add(identity)
    try:
        if type(value) is list:
            for index, entry in enumerate(value):
                _require_exact_json_value(
                    entry,
                    label=f"{label}[{index}]",
                    ancestors=ancestors,
                )
            return
        if type(value) is dict:
            for key, entry in value.items():
                if type(key) is not str:
                    raise TypeError(f"{label} object keys must be native strings")
                _require_exact_json_value(
                    entry,
                    label=f"{label}.{key}",
                    ancestors=ancestors,
                )
            return
        raise TypeError(f"{label} must contain only exact JSON values")
    finally:
        ancestors.remove(identity)


def _require_cell_indices(value: object, *, label: str = "cell_indices") -> list[int]:
    if type(value) is not list:
        raise TypeError(f"{label} must be a native list")
    seen: set[int] = set()
    for index, cell_index in enumerate(value):
        if type(cell_index) is not int or cell_index < 0:
            raise TypeError(f"{label}[{index}] must be a non-negative native integer")
        if cell_index in seen:
            raise ValueError(f"{label} must not contain duplicate cell indices")
        seen.add(cell_index)
    return value


def _require_frontend_message(message: object) -> str:
    if type(message) is not dict:
        raise TypeError("message must be a native dictionary")
    _require_exact_json_value(message, label="message")
    message_type = _require_exact_message_text(
        message.get("type"),
        label="message type",
    )
    if "viewerId" in message or "viewerToken" in message:
        raise TypeError("message must not supply viewer routing credentials")

    expected_keys: dict[str, set[str]] = {
        "ping": {"type", "requestId"},
        "debug_snapshot": {"type", "requestId"},
        "requestSessionBundle": {"type", "requestId"},
        "highlight": {"type", "cells", "color"},
        "clearHighlights": {"type"},
        "setColorBy": {"type", "field"},
        "setVisibility": {"type", "cells", "visible"},
        "resetCamera": {"type"},
        "freeze": {"type"},
    }
    expected = expected_keys.get(message_type)
    if expected is None:
        raise ValueError(f"Unknown current Jupyter command type {message_type!r}")
    if set(message) != expected:
        fields = ", ".join(sorted(expected_keys[message_type]))
        raise TypeError(f"{message_type} message must contain exactly {fields}")

    if message_type in {"ping", "debug_snapshot", "requestSessionBundle"}:
        _require_exact_message_text(
            message["requestId"],
            label=f"{message_type} requestId",
        )
    elif message_type == "highlight":
        highlight_cells = _require_cell_indices(
            message["cells"],
            label="highlight cells",
        )
        if not highlight_cells:
            raise ValueError("highlight cells must be a non-empty list")
        color = message["color"]
        if (
            type(color) is not str
            or len(color) != 7
            or color[0] != "#"
            or any(character not in "0123456789abcdefABCDEF" for character in color[1:])
        ):
            raise ValueError("highlight color must be an exact six-digit hex color")
    elif message_type == "setColorBy":
        _require_exact_message_text(
            message["field"],
            label="setColorBy field",
            allow_whitespace=True,
        )
    elif message_type == "setVisibility":
        cells = message["cells"]
        if cells is not None:
            _require_cell_indices(cells, label="setVisibility cells")
        if type(message["visible"]) is not bool:
            raise TypeError("setVisibility visible must be exactly True or False")

    return json.dumps(
        message,
        ensure_ascii=True,
        allow_nan=False,
        separators=(",", ":"),
    )


def _require_nonnegative_event_integer(value: object, *, label: str) -> int:
    if type(value) is not int:
        raise TypeError(f"{label} must be a native integer")
    if value < 0:
        raise ValueError(f"{label} must be non-negative")
    return value


def _require_exact_event_fields(
    message: dict[str, Any],
    expected: set[str],
    *,
    event_type: str,
) -> None:
    missing = sorted(expected - set(message))
    unknown = sorted(set(message) - expected)
    if missing or unknown:
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if unknown:
            details.append("unknown " + ", ".join(unknown))
        raise ValueError(
            f"Inbound Jupyter {event_type} event has invalid fields ({'; '.join(details)})"
        )


def _require_inbound_jupyter_event(
    message: object,
    *,
    expected_viewer_id: str,
) -> tuple[str, dict[str, Any]]:
    if type(message) is not dict:
        raise TypeError("Inbound Jupyter event must be a native dictionary")
    _require_exact_json_value(message, label="Inbound Jupyter event")
    event_type = _require_exact_message_text(
        message.get("type"),
        label="Inbound Jupyter event type",
    )
    viewer_id = _require_exact_message_text(
        message.get("viewerId"),
        label="Inbound Jupyter event viewerId",
    )
    if viewer_id != expected_viewer_id:
        raise ValueError("Inbound Jupyter event viewerId does not match this viewer")

    fields_by_type: dict[str, set[str]] = {
        "selection": {"type", "viewerId", "cells", "source"},
        "hover": {"type", "viewerId", "cell", "position"},
        "click": {
            "type",
            "viewerId",
            "cell",
            "button",
            "shift",
            "ctrl",
        },
        "ready": {"type", "viewerId", "n_cells", "dimensions"},
        "pong": {"type", "viewerId", "requestId", "t"},
        "debug_snapshot": {
            "type",
            "viewerId",
            "requestId",
            "ts",
            "locationHref",
            "origin",
            "serverUrl",
            "connected",
            "parentOrigin",
            "userAgent",
        },
        "session_bundle": {
            "type",
            "viewerId",
            "requestId",
            "status",
            "bytes",
            "path",
        },
    }
    expected_fields = fields_by_type.get(event_type)
    if expected_fields is None:
        raise ValueError(f"Unknown inbound Jupyter event type {event_type!r}")
    _require_exact_event_fields(
        message,
        expected_fields,
        event_type=event_type,
    )

    if event_type == "selection":
        _require_cell_indices(message["cells"], label="selection cells")
        _require_exact_message_text(
            message["source"],
            label="selection source",
        )
    elif event_type == "hover":
        cell = message["cell"]
        if cell is not None:
            _require_nonnegative_event_integer(cell, label="hover cell")
        position = message["position"]
        if position is not None:
            if type(position) is not dict:
                raise TypeError("hover position must be a native dictionary or None")
            _require_exact_event_fields(
                position,
                {"x", "y", "z"},
                event_type="hover position",
            )
            for axis in ("x", "y", "z"):
                coordinate = position[axis]
                if type(coordinate) not in {int, float} or not math.isfinite(coordinate):
                    raise TypeError(f"hover position {axis} must be a finite number")
    elif event_type == "click":
        _require_nonnegative_event_integer(message["cell"], label="click cell")
        button = message["button"]
        if type(button) is not int or button not in {0, 1, 2}:
            raise TypeError("click button must be the native integer 0, 1, or 2")
        if type(message["shift"]) is not bool or type(message["ctrl"]) is not bool:
            raise TypeError("click shift and ctrl must be native booleans")
    elif event_type == "ready":
        _require_nonnegative_event_integer(message["n_cells"], label="ready n_cells")
        if type(message["dimensions"]) is not int or message["dimensions"] not in {
            1,
            2,
            3,
        }:
            raise TypeError("ready dimensions must be the native integer 1, 2, or 3")
    elif event_type == "pong":
        _require_exact_message_text(
            message["requestId"],
            label="pong requestId",
        )
        _require_nonnegative_event_integer(message["t"], label="pong t")
    elif event_type == "debug_snapshot":
        for key in (
            "requestId",
            "ts",
            "locationHref",
            "origin",
            "serverUrl",
            "parentOrigin",
        ):
            _require_exact_message_text(
                message[key],
                label=f"debug_snapshot {key}",
                allow_whitespace=True,
            )
        if type(message["connected"]) is not bool:
            raise TypeError("debug_snapshot connected must be a native boolean")
        user_agent = message["userAgent"]
        if user_agent is not None:
            _require_exact_message_text(
                user_agent,
                label="debug_snapshot userAgent",
                allow_whitespace=True,
            )
    else:
        if message["status"] != "ok":
            raise ValueError("session_bundle status must equal 'ok'")
        _require_exact_message_text(
            message["requestId"],
            label="session_bundle requestId",
        )
        _require_nonnegative_event_integer(
            message["bytes"],
            label="session_bundle bytes",
        )
        _require_exact_message_text(
            message["path"],
            label="session_bundle path",
            allow_whitespace=True,
        )

    return (
        event_type,
        {key: value for key, value in message.items() if key not in {"type", "viewerId"}},
    )


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
        message_callbacks = list(self._callbacks.get("message", ()))
        for callback in message_callbacks:
            callback({"event": event, **data})

        if event == "message":
            return

        event_callbacks = list(self._callbacks.get(event, ()))
        for callback in event_callbacks:
            callback(data)

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
        *,
        client_server_url: str | None = None,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
    ):
        """Initialize common viewer properties."""
        self.port = 0 if port is None else require_server_port(port)
        self.height = height
        self._client_server_url = (
            _require_client_server_url(client_server_url) if client_server_url is not None else None
        )
        self.web_source_url = web_source_url
        self.web_cache_dir = (
            Path(web_cache_dir).expanduser().resolve()
            if web_cache_dir is not None
            else _web_cache_dir()
        )

        self._server: _ViewerServer | None = None
        self._viewer_id = secrets.token_hex(8)
        self._viewer_token = secrets.token_hex(16)
        self._context = _detect_jupyter_context()
        self._displayed = False
        self._message_routing_registered = False
        self._routing_finalizer: weakref.finalize | None = None

        # Initialize hooks system
        self._hooks = HookRegistry()

        # Latest-state snapshot + event waiting primitives.
        self.state = ViewerState()
        self._state_lock = threading.Lock()
        self._event_cv = threading.Condition(self._state_lock)
        self._event_seq = 0
        self._recent_events: deque[tuple[int, str, dict[str, Any]]] = deque(maxlen=512)

    def _activate(self) -> None:
        """Publish this fully started viewer to the routing registries."""
        _register_viewer_for_messages(self)
        try:
            _active_viewers.add(self)
        except BaseException:
            _unregister_viewer_for_messages(self)
            raise

    def _rollback_failed_construction(self) -> None:
        """Undo every resource acquired during an unsuccessful constructor."""
        failures: list[BaseException] = []
        server = self._server
        if server is not None:
            try:
                server.stop()
            except BaseException as error:
                failures.append(error)
            else:
                self._server = None
        try:
            _unregister_viewer_for_messages(self)
        except BaseException as error:
            failures.append(error)
        _active_viewers.discard(self)
        self._hooks.clear()
        if failures:
            details = "; ".join(f"{type(error).__name__}: {error}" for error in failures)
            raise RuntimeError(
                f"Failed viewer construction could not be rolled back: {details}"
            ) from failures[0]

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
        return self._hooks.on("selection")

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
        return self._hooks.on("hover")

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
        return self._hooks.on("click")

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
        return self._hooks.on("ready")

    @property
    def on_message(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a raw message handler.

        Called for every validated current-schema message from the viewer.
        Use it to observe the complete event stream while debugging.

        Event data:
            - event: str - the event type
            - ... other fields depend on event type

        Example:
            @viewer.on_message
            def handle(event):
                print(f"Got message: {event}")
        """
        return self._hooks.on("message")

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

    def _handle_frontend_message(self, message: dict[str, Any]) -> None:
        """
        Handle a message received from the frontend.

        This is called by the message routing system when the frontend
        sends data back to Python.

        Args:
            message: Message dict from frontend
        """
        event, data = _require_inbound_jupyter_event(
            message,
            expected_viewer_id=self._viewer_id,
        )
        self._hooks.trigger(event, data)
        self._record_event(event, data)

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
            self._viewer_token,
            request_id,
            ttl_seconds=(3600.0 if remaining is None else remaining),
        )
        try:
            self.send_message(
                {
                    "type": "requestSessionBundle",
                    "requestId": request_id,
                }
            )
            event = self.wait_for_event(
                "session_bundle",
                timeout=remaining,
                predicate=lambda e: e.get("requestId") == request_id,
            )
        finally:
            cancel_session_bundle_request(
                self._viewer_id,
                self._viewer_token,
                request_id,
            )

        if event.get("status") != "ok":
            raise RuntimeError(event.get("error") or "Failed to capture session bundle")

        path = event.get("path")
        if not path:
            raise RuntimeError("Session bundle upload succeeded but no temp path was returned")

        return CellucidSessionBundle(Path(path))

    def apply_session_to_anndata(
        self,
        adata: Any,
        *,
        expected_dataset_id: str | None = None,
        inplace: bool = False,
        timeout: float | None = 60.0,
        cleanup_bundle: bool = True,
        add_highlights: bool = True,
        highlights_prefix: str = "cellucid_highlight__",
        add_user_defined_fields: bool = True,
        user_defined_prefix: str = "",
        include_deleted_user_defined_fields: bool = False,
        store_uns: bool = True,
        return_summary: bool = False,
    ):
        """
        Convenience wrapper: capture a session bundle and apply it to an AnnData.

        Notes
        -----
        - If you want to keep the artifact, call `bundle = viewer.get_session_bundle()`
          and `bundle.save(...)` before applying.
        - By default, the temporary bundle file created by the server is deleted
          after applying.
        """
        if type(cleanup_bundle) is not bool:
            raise TypeError("cleanup_bundle must be exactly True or False")
        bundle = self.get_session_bundle(timeout=timeout)
        try:
            resolved_dataset_id = expected_dataset_id
            if resolved_dataset_id is None:
                server = self._server
                adapter = getattr(server, "adapter", None)
                if adapter is None or not hasattr(adapter, "get_dataset_identity"):
                    raise ValueError(
                        "expected_dataset_id is required when the viewer is not "
                        "served directly from AnnData"
                    )
                identity = adapter.get_dataset_identity()
                if not isinstance(identity, dict):
                    raise TypeError("AnnData adapter identity must be a dictionary.")
                resolved_dataset_id = identity.get("id")
                if not isinstance(resolved_dataset_id, str) or not resolved_dataset_id:
                    raise ValueError("AnnData adapter identity must contain a non-empty id.")
            return bundle.apply_to_anndata(
                adata,
                expected_dataset_id=resolved_dataset_id,
                inplace=inplace,
                add_highlights=add_highlights,
                highlights_prefix=highlights_prefix,
                add_user_defined_fields=add_user_defined_fields,
                user_defined_prefix=user_defined_prefix,
                include_deleted_user_defined_fields=include_deleted_user_defined_fields,
                store_uns=store_uns,
                return_summary=return_summary,
            )
        finally:
            if cleanup_bundle:
                bundle.path.unlink(missing_ok=True)

    def debug_connection(self, timeout: float | None = 5.0) -> dict[str, Any]:
        """
        Return a structured connectivity/debug report for this viewer.

        Includes:
        - server health/info probes
        - Python→Frontend ping/pong roundtrip (postMessage + HTTP events)
        - frontend debug snapshot (URL/origin/userAgent as seen by the iframe)
        """
        report: dict[str, Any] = {
            "viewer_id": self._viewer_id,
            "viewer_url": self.viewer_url,
            "server_url": self.server_url,
            "displayed": self._displayed,
            "notebook_context": dict(self._context),
            "server_running": bool(
                self._server and getattr(self._server, "is_running", lambda: True)()
            ),
            "web_ui": {
                "cache_dir": str(self.web_cache_dir),
                "source": self.web_source_url,
            },
            "state": {
                "ready": self.state.ready,
                "last_event_type": self.state.last_event_type,
                "last_updated_at": self.state.last_updated_at,
            },
        }

        try:
            report["client_server_url"] = self._get_client_server_url()
        except Exception as e:
            report["client_server_url_error"] = str(e)

        # Exact cache introspection.
        try:
            from .web_cache import verify_web_ui_cache

            cache_dir = self.web_cache_dir
            index_path = cache_dir / "index.html"
            cache_info: dict[str, Any] = {
                "cache_dir_exists": cache_dir.exists(),
                "index_html_exists": index_path.is_file(),
            }
            if index_path.is_file():
                data = index_path.read_bytes()
                cache_info["index_html_bytes"] = len(data)
                cache_info["index_html_build_id"] = _extract_web_build_id(data)
                cache_info["index_html_mtime"] = index_path.stat().st_mtime
                inventory = verify_web_ui_cache(
                    cache_dir,
                    expected_source_url=self.web_source_url,
                )
                cache_info["verified_build_id"] = inventory.build_id
                cache_info["verified_asset_count"] = len(inventory.assets)
            report["web_ui"]["cache"] = cache_info
        except Exception as e:
            report["web_ui"]["cache_error"] = str(e)

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

        # Probe dataset identity using the first declared dataset entry.
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

        # Inspect the locally served, inventory-verified index.
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
            report["frontend_roundtrip"] = {
                "ok": False,
                "error": "Viewer not displayed (call viewer.display())",
            }
            report["frontend_debug_snapshot"] = {
                "ok": False,
                "error": "Viewer not displayed (call viewer.display())",
            }
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

        return report

    # =========================================================================
    # WEB UI CACHE UTILITIES
    # =========================================================================

    def clear_web_cache(self) -> Path:
        """Clear this viewer's selected web UI cache."""
        from .web_cache import clear_web_cache

        return clear_web_cache(cache_dir=self.web_cache_dir)

    def ensure_web_ui_cached(self, *, force: bool = True, show_progress: bool = True):
        """Establish one complete verified web build in this viewer's cache."""
        from .web_cache import ensure_web_ui_cached

        return ensure_web_ui_cached(
            cache_dir=self.web_cache_dir,
            source_url=self.web_source_url,
            force=force,
            show_progress=show_progress,
        )

    # =========================================================================
    # DISPLAY & SERVER
    # =========================================================================

    @property
    def server_url(self) -> str:
        """Get the data server URL."""
        if self._server is not None and getattr(self._server, "url", None):
            return str(self._server.url)
        return f"http://127.0.0.1:{self.port}"

    def _get_client_server_url(self) -> str:
        """Return the one browser URL selected for this viewer."""
        if self._client_server_url is not None:
            return self._client_server_url
        if self._server is None:
            raise RuntimeError("The Cellucid data server has not started")
        return _require_client_server_url(str(self._server.url))

    @property
    def viewer_origin(self) -> str:
        """Origin (scheme+host+port) for postMessage targetOrigin."""
        parsed = urlparse(self.viewer_url)
        return f"{parsed.scheme}://{parsed.netloc}"

    @property
    def viewer_url(self) -> str:
        """Get the full viewer URL. Subclasses can override to add extra params."""
        from urllib.parse import urlencode

        base = self._get_client_server_url()
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

        from IPython.display import HTML, display

        from .web_cache import verify_web_ui_cache

        verify_web_ui_cache(
            self.web_cache_dir,
            expected_source_url=self.web_source_url,
        )

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
        viewer_src = self.viewer_url
        target_origin = self.viewer_origin

        return f"""
        <div id="cellucid-viewer-{self._viewer_id}" style="width:100%; height:{self.height}px;">
            <iframe
                id="cellucid-iframe-{self._viewer_id}"
                src="{html.escape(viewer_src, quote=True)}"
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
            var targetOrigin = {json.dumps(target_origin)};

            // Set up message passing infrastructure
            window.cellucidViewers = window.cellucidViewers || {{}};
            window.cellucidViewers[viewerId] = {{
                sendMessage: function(msg) {{
                    var iframe = document.getElementById('cellucid-iframe-' + viewerId);
                    if (!iframe || !iframe.contentWindow) {{
                        throw new Error('Cellucid viewer iframe is unavailable: ' + viewerId);
                    }}
                    iframe.contentWindow.postMessage(
                        Object.assign({{}}, msg, {{ viewerId: viewerId, viewerToken: viewerToken }}),
                        targetOrigin
                    );
                }}
            }};

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
        serialized_message = _require_frontend_message(message)
        if not self._displayed:
            raise RuntimeError("Viewer must be displayed before sending a message.")

        from IPython.display import Javascript, display

        js = f"""
        (function() {{
            var viewer = window.cellucidViewers && window.cellucidViewers['{self._viewer_id}'];
            if (!viewer) {{
                throw new Error('Cellucid viewer command target is unavailable: {self._viewer_id}');
            }}
            viewer.sendMessage({serialized_message});
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
        exact_indices = _require_cell_indices(cell_indices)
        if not exact_indices:
            raise ValueError("cell_indices must be a non-empty list")
        if (
            type(color) is not str
            or len(color) != 7
            or color[0] != "#"
            or any(character not in "0123456789abcdefABCDEF" for character in color[1:])
        ):
            raise ValueError("color must be an exact six-digit hex color")
        self.send_message(
            {
                "type": "highlight",
                "cells": exact_indices,
                "color": color,
            }
        )

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
        exact_field = _require_exact_message_text(
            field,
            label="field",
            allow_whitespace=True,
        )
        self.send_message({"type": "setColorBy", "field": exact_field})

    def set_visibility(self, cell_indices: list[int] | None = None, visible: bool = True):
        """
        Set visibility of specific cells.

        Args:
            cell_indices: List of cell indices (None for all cells)
            visible: Whether cells should be visible

        Example:
            viewer.set_visibility([0, 1, 2], visible=False)  # Hide cells
        """
        if type(visible) is not bool:
            raise TypeError("visible must be exactly True or False")
        exact_indices = None if cell_indices is None else _require_cell_indices(cell_indices)
        self.send_message(
            {
                "type": "setVisibility",
                "cells": exact_indices,
                "visible": visible,
            }
        )

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

        Safe to call multiple times after a successful stop.
        """
        failures: list[BaseException] = []
        if self._displayed:
            try:
                self.send_message({"type": "freeze"})
            except BaseException as error:
                failures.append(error)

        if self._server:
            try:
                self._server.stop()
            except BaseException as error:
                failures.append(error)
            else:
                self._server = None

        # Clear hooks
        self._hooks.clear()

        # Unregister from message routing
        _unregister_viewer_for_messages(self)

        # Remove from active viewers set to prevent double cleanup
        _active_viewers.discard(self)

        logger.info(f"{self.__class__.__name__} stopped")
        if failures:
            details = "; ".join(f"{type(error).__name__}: {error}" for error in failures)
            raise RuntimeError(f"Cellucid viewer shutdown failed: {details}") from failures[0]

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
            raise RuntimeError("The registered viewer no longer exists")
        resolved._handle_frontend_message(event)

    register_event_callback(
        viewer._viewer_id,
        viewer._viewer_token,
        _deliver,
    )
    try:
        viewer._routing_finalizer = weakref.finalize(
            viewer,
            unregister_event_callback,
            viewer._viewer_id,
            viewer._viewer_token,
        )
    except BaseException:
        unregister_event_callback(
            viewer._viewer_id,
            viewer._viewer_token,
        )
        raise
    viewer._message_routing_registered = True
    logger.debug(f"Registered HTTP event callback for viewer {viewer._viewer_id}")


def _unregister_viewer_for_messages(viewer: BaseViewer):
    """Unregister a viewer from message routing."""
    finalizer = viewer._routing_finalizer
    if finalizer is not None and finalizer.alive:
        finalizer()
    elif viewer._message_routing_registered:
        unregister_event_callback(
            viewer._viewer_id,
            viewer._viewer_token,
        )
    viewer._routing_finalizer = None
    viewer._message_routing_registered = False


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
        *,
        client_server_url: str | None = None,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
    ):
        """
        Initialize the viewer.

        Args:
            data_dir: Path to the cellucid dataset directory
            port: Port for the data server (auto-selected if None)
            height: Height of the embedded viewer in pixels
            auto_open: Automatically display when created
            client_server_url: Exact browser-reachable data server base URL.
            web_source_url: Origin publishing the web asset inventory.
            web_cache_dir: Directory holding the active verified web build.
        """
        super().__init__(
            port=port,
            height=height,
            client_server_url=client_server_url,
            web_source_url=web_source_url,
            web_cache_dir=web_cache_dir,
        )

        try:
            self.data_dir = Path(data_dir).resolve()

            if not self.data_dir.exists():
                raise FileNotFoundError(f"Data directory not found: {self.data_dir}")

            self._start_server()
            self._activate()

            if auto_open and self._context["in_jupyter"]:
                self.display()
        except BaseException:
            self._rollback_failed_construction()
            raise

    def _start_server(self):
        """Start the background data server."""
        self._server = CellucidServer(
            data_dir=self.data_dir,
            port=self.port,
            host="127.0.0.1",
            open_browser=False,
            quiet=True,
            serve_web_ui=True,
            web_source_url=self.web_source_url,
            web_cache_dir=self.web_cache_dir,
        )
        self._server.start_background()
        self.port = self._server.port
        logger.info(f"Started cellucid server at {self._server.url}")

    def __repr__(self) -> str:
        status = "running" if self._server and self._server.is_running() else "stopped"
        return f"CellucidViewer('{self.data_dir}', port={self.port}, status={status})"


def show(
    data_dir: str | Path,
    height: int = 600,
    *,
    client_server_url: str | None = None,
    web_source_url: str = CELLUCID_WEB_URL,
    web_cache_dir: str | Path | None = None,
) -> CellucidViewer:
    """
    Quick function to display a cellucid dataset in a notebook.

    Args:
        data_dir: Path to the dataset directory
        height: Height of the viewer in pixels
        client_server_url: Exact browser-reachable data server base URL.
        web_source_url: Origin publishing the web asset inventory.
        web_cache_dir: Directory holding the active verified web build.

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
        client_server_url=client_server_url,
        web_source_url=web_source_url,
        web_cache_dir=web_cache_dir,
    )


def cleanup_all():
    """Stop all active viewers and their servers."""
    for viewer in list(_active_viewers):
        try:
            viewer.stop()
        except Exception:
            logger.exception("Failed to stop %r during interpreter shutdown", viewer)


# Register cleanup on interpreter exit
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
    - zarr stores materialized by ``anndata.read_zarr``

    Example:
        >>> viewer = AnnDataViewer(
        ...     adata,
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ... )
        >>> viewer.display()
        >>>
        >>> @viewer.on_selection
        ... def analyze_selection(event):
        ...     subset = adata[event['cells']]
        ...     sc.pl.violin(subset, ['gene1', 'gene2'])
        >>>
        >>> # From h5ad file with lazy loading
        >>> viewer = AnnDataViewer(
        ...     "/path/to/data.h5ad",
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ... )
    """

    def __init__(
        self,
        data: str | Path | anndata.AnnData,
        port: int | None = None,
        height: int = 600,
        auto_open: bool = True,
        *,
        latent_key: str | None = None,
        gene_id_column: str | None = None,
        normalize_embeddings: bool = True,
        centroid_outlier_quantile: float = 0.95,
        centroid_min_points: int = 10,
        dataset_name: str,
        dataset_id: str,
        vector_field_default: str | None = None,
        client_server_url: str | None = None,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
    ) -> None:
        """
        Initialize the AnnData viewer.

        Args:
            data: AnnData object or path to h5ad file or zarr directory.
            port: Port for the data server (auto-selected if None).
            height: Height of the embedded viewer in pixels.
            auto_open: Automatically display when created.
            latent_key: Explicit key in ``obsm`` for the latent space.
            gene_id_column: Exact column in ``var`` containing gene identifiers.
                If None, identifiers come from ``var.index``.
            normalize_embeddings: Whether to normalize embeddings into the viewer
                coordinate range.
            centroid_outlier_quantile: Quantile used for categorical centroids.
            centroid_min_points: Minimum category size used for centroid computation.
            dataset_name: Explicit human-readable dataset name.
            dataset_id: Explicit stable dataset identifier.
            vector_field_default: Exact field id required when multiple UMAP
                vector fields exist.
            client_server_url: Exact browser-reachable data server base URL.
            web_source_url: Origin publishing the web asset inventory.
            web_cache_dir: Directory holding the active verified web build.
        """
        super().__init__(
            port=port,
            height=height,
            client_server_url=client_server_url,
            web_source_url=web_source_url,
            web_cache_dir=web_cache_dir,
        )

        try:
            self.data = data

            self._start_server(
                latent_key=latent_key,
                gene_id_column=gene_id_column,
                normalize_embeddings=normalize_embeddings,
                centroid_outlier_quantile=centroid_outlier_quantile,
                centroid_min_points=centroid_min_points,
                dataset_name=dataset_name,
                dataset_id=dataset_id,
                vector_field_default=vector_field_default,
            )
            self._activate()

            if auto_open and self._context["in_jupyter"]:
                self.display()
        except BaseException:
            self._rollback_failed_construction()
            raise

    def _start_server(
        self,
        *,
        latent_key: str | None,
        gene_id_column: str | None,
        normalize_embeddings: bool,
        centroid_outlier_quantile: float,
        centroid_min_points: int,
        dataset_name: str,
        dataset_id: str,
        vector_field_default: str | None,
    ) -> None:
        """Start the AnnData server."""
        from .anndata_server import AnnDataServer

        self._server = AnnDataServer(
            data=self.data,
            port=self.port,
            host="127.0.0.1",
            open_browser=False,
            quiet=True,
            serve_web_ui=True,
            web_source_url=self.web_source_url,
            web_cache_dir=self.web_cache_dir,
            latent_key=latent_key,
            gene_id_column=gene_id_column,
            normalize_embeddings=normalize_embeddings,
            centroid_outlier_quantile=centroid_outlier_quantile,
            centroid_min_points=centroid_min_points,
            dataset_name=dataset_name,
            dataset_id=dataset_id,
            vector_field_default=vector_field_default,
        )
        self._server.start_background()
        self.port = self._server.port
        logger.info(f"Started AnnData server at {self._server.url}")

    @property
    def viewer_url(self) -> str:
        """Get the full viewer URL with anndata flag."""
        from urllib.parse import urlencode

        base = self._get_client_server_url()
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
    data: str | Path | anndata.AnnData,
    height: int = 600,
    *,
    latent_key: str | None = None,
    gene_id_column: str | None = None,
    normalize_embeddings: bool = True,
    centroid_outlier_quantile: float = 0.95,
    centroid_min_points: int = 10,
    dataset_name: str,
    dataset_id: str,
    vector_field_default: str | None = None,
    client_server_url: str | None = None,
    web_source_url: str = CELLUCID_WEB_URL,
    web_cache_dir: str | Path | None = None,
) -> AnnDataViewer:
    """
    Quickly display an AnnData object, h5ad file, or zarr store in a notebook.

    This is the easiest way to visualize AnnData directly without running
    prepare first. However, it's slower than using pre-exported data.

    Args:
        data: AnnData object, path to h5ad file, or path to zarr directory.
        height: Height of the viewer in pixels.
        latent_key: Explicit key in ``obsm`` for the latent space.
        gene_id_column: Exact column in ``var`` containing gene identifiers.
            If None, identifiers come from ``var.index``.
        normalize_embeddings: Whether to normalize embeddings into the viewer
            coordinate range.
        centroid_outlier_quantile: Quantile used for categorical centroids.
        centroid_min_points: Minimum category size used for centroid computation.
        dataset_name: Explicit human-readable dataset name.
        dataset_id: Explicit stable dataset identifier.
        vector_field_default: Exact field id required when multiple UMAP
            vector fields exist.
        client_server_url: Exact browser-reachable data server base URL.
        web_source_url: Origin publishing the web asset inventory.
        web_cache_dir: Directory holding the active verified web build.

    Returns:
        AnnDataViewer instance for interaction via hooks.

    Supported formats:
        - In-memory AnnData objects
        - .h5ad files (HDF5-based, lazy loading via backed mode)
        - .zarr directories materialized by ``anndata.read_zarr``

    Example:
        >>> from cellucid import show_anndata
        >>> viewer = show_anndata(
        ...     adata,
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ... )
        >>>
        >>> @viewer.on_selection
        ... def handle(event):
        ...     subset = adata[event['cells']]
        ...     sc.pl.violin(subset, 'gene1')

        >>> # With custom options
        >>> viewer = show_anndata(
        ...     adata,
        ...     latent_key="X_pca",
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ...     height=800,
        ... )

    Note:
        For production use or sharing, consider using prepare
        to create optimized binary files, then use show() to display them.
    """
    return AnnDataViewer(
        data=data,
        height=height,
        auto_open=True,
        latent_key=latent_key,
        gene_id_column=gene_id_column,
        normalize_embeddings=normalize_embeddings,
        centroid_outlier_quantile=centroid_outlier_quantile,
        centroid_min_points=centroid_min_points,
        dataset_name=dataset_name,
        dataset_id=dataset_id,
        vector_field_default=vector_field_default,
        client_server_url=client_server_url,
        web_source_url=web_source_url,
        web_cache_dir=web_cache_dir,
    )
