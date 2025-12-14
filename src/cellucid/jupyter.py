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
import weakref
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

if TYPE_CHECKING:
    import anndata

from .server import CellucidServer, DEFAULT_PORT
from ._server_base import register_event_callback, unregister_event_callback

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


def _get_notebook_url() -> str | None:
    """Try to get the notebook server URL."""
    try:
        from notebook import notebookapp
        servers = list(notebookapp.list_running_servers())
        if servers:
            return servers[0].get("url")
    except Exception:
        pass

    try:
        from jupyter_server import serverapp
        servers = list(serverapp.list_running_servers())
        if servers:
            return servers[0].get("url")
    except Exception:
        pass

    return None


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
        self._context = _detect_jupyter_context()
        self._displayed = False

        # Initialize hooks system
        self._hooks = HookRegistry()

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
            - source: str - selection method ('lasso', 'click', 'range')

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
            - dimensions: int - embedding dimensionality (2 or 3)

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

        self._hooks.trigger(event, data)

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
        return f"http://127.0.0.1:{self.port}"

    @property
    def viewer_url(self) -> str:
        """Get the full viewer URL. Subclasses can override to add extra params."""
        params = [
            f"remote={self.server_url}",
            f"jupyter=true",
            f"viewerId={self._viewer_id}",
        ]
        return f"https://www.cellucid.com?{'&'.join(params)}"

    def _get_pre_display_html(self) -> str | None:
        """Override to add HTML before the viewer (e.g., warnings)."""
        return None

    def display(self):
        """Display the viewer in a notebook cell."""
        if not self._context["in_jupyter"]:
            print(f"Not in Jupyter environment. Open manually: {self.viewer_url}")
            return

        from IPython.display import display, HTML

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
        return f"""
        <div id="cellucid-viewer-{self._viewer_id}" style="width:100%; height:{self.height}px;">
            <iframe
                id="cellucid-iframe-{self._viewer_id}"
                src="{self.viewer_url}"
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

            // Set up message passing infrastructure
            window.cellucidViewers = window.cellucidViewers || {{}};
            window.cellucidViewers[viewerId] = {{
                sendMessage: function(msg) {{
                    var iframe = document.getElementById('cellucid-iframe-' + viewerId);
                    if (iframe && iframe.contentWindow) {{
                        iframe.contentWindow.postMessage(
                            {{...msg, viewerId: viewerId}},
                            'https://www.cellucid.com'
                        );
                    }}
                }}
            }};

            // Note: Frontend → Python communication uses HTTP POST to /_cellucid/events.
            // The postMessage listener is no longer needed for event routing.
        }})();
        </script>
        """

    def send_message(self, message: dict):
        """
        Send a message to the viewer iframe.

        This is the low-level API for sending commands to the frontend.
        Messages are sent via postMessage to the cellucid.com iframe.

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
# 1. The viewer iframe on cellucid.com POSTs events to the local data server
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

# Global registry of viewers by ID for message routing
_viewer_registry: dict[str, weakref.ref[BaseViewer]] = {}


def _register_viewer_for_messages(viewer: BaseViewer):
    """
    Register a viewer to receive messages from frontend.

    This sets up the HTTP-based event routing which works in ALL environments.
    The frontend POSTs events to /_cellucid/events on the data server,
    which routes them to the appropriate viewer.
    """
    _viewer_registry[viewer._viewer_id] = weakref.ref(viewer)

    # Register HTTP event callback - this works in ALL environments!
    # The data server (which is already running) will route POSTed events here.
    register_event_callback(
        viewer._viewer_id,
        viewer._handle_frontend_message
    )
    logger.debug(f"Registered HTTP event callback for viewer {viewer._viewer_id}")


def _unregister_viewer_for_messages(viewer: BaseViewer):
    """Unregister a viewer from message routing."""
    _viewer_registry.pop(viewer._viewer_id, None)
    unregister_event_callback(viewer._viewer_id)


def _route_message_to_viewer(viewer_id: str, message: dict):
    """Route a message from frontend to the appropriate viewer."""
    viewer_ref = _viewer_registry.get(viewer_id)
    if viewer_ref:
        viewer = viewer_ref()
        if viewer:
            viewer._handle_frontend_message(message)


def _check_hook_support(context: dict) -> bool:
    """
    Check if the current environment supports frontend → Python hooks.

    Returns True if hooks like on_selection will actually fire.
    With HTTP-based event routing, this now returns True for all
    Jupyter environments since the local data server can receive events.
    """
    # HTTP-based routing works in all Jupyter environments
    return context.get('in_jupyter', False)


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
        logger.info(f"Started AnnData server at {self._server.url}")

    @property
    def viewer_url(self) -> str:
        """Get the full viewer URL with anndata flag."""
        params = [
            f"remote={self.server_url}",
            f"anndata=true",
            f"jupyter=true",
            f"viewerId={self._viewer_id}",
        ]
        return f"https://www.cellucid.com?{'&'.join(params)}"

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
