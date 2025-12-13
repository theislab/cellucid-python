"""
Cellucid Jupyter Integration

Provides seamless integration with Jupyter notebooks for visualizing
cellucid datasets directly in notebook cells.

Usage:
    from cellucid.jupyter import CellucidViewer, show

    # Quick visualization of a data directory
    show("/path/to/my_dataset")

    # Or with more control
    viewer = CellucidViewer("/path/to/my_dataset")
    viewer.display()  # Shows in notebook cell

    # For live interaction
    viewer = CellucidViewer("/path/to/my_dataset", interactive=True)
    viewer.display()
    viewer.highlight_cells([1, 2, 3])  # Highlight specific cells
    viewer.set_color_by("cell_type")  # Change coloring
"""

from __future__ import annotations

import json
import logging
import os
import secrets
import threading
import weakref
from pathlib import Path
from typing import Any, Callable

from .server import CellucidServer, CellucidServerAsync, DEFAULT_PORT

logger = logging.getLogger("cellucid.jupyter")

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


class CellucidViewer:
    """
    Interactive cellucid viewer for Jupyter notebooks.

    Embeds the cellucid web viewer in a notebook cell, connected to a
    local data server. Supports bidirectional communication for
    interactive analysis.

    Example:
        >>> viewer = CellucidViewer("/path/to/dataset")
        >>> viewer.display()  # Show in notebook

        # Highlight specific cells
        >>> viewer.highlight_cells([100, 200, 300])

        # Change coloring
        >>> viewer.set_color_by("cluster")

        # Get currently selected cells
        >>> selected = viewer.get_selected_cells()
    """

    def __init__(
        self,
        data_dir: str | Path,
        port: int | None = None,
        height: int = 600,
        interactive: bool = True,
        auto_open: bool = True,
    ):
        """
        Initialize the viewer.

        Args:
            data_dir: Path to the cellucid dataset directory
            port: Port for the data server (auto-selected if None)
            height: Height of the embedded viewer in pixels
            interactive: Enable bidirectional communication
            auto_open: Automatically display when created
        """
        self.data_dir = Path(data_dir).resolve()
        self.port = port or self._find_free_port()
        self.height = height
        self.interactive = interactive

        self._server: CellucidServer | None = None
        self._viewer_id = secrets.token_hex(8)
        self._context = _detect_jupyter_context()
        self._message_handlers: list[Callable] = []
        self._displayed = False

        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {self.data_dir}")

        # Start the server
        self._start_server()

        # Register for cleanup
        _active_viewers.add(self)

        if auto_open and self._context["in_jupyter"]:
            self.display()

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

    @property
    def server_url(self) -> str:
        """Get the data server URL."""
        return f"http://127.0.0.1:{self.port}"

    @property
    def viewer_url(self) -> str:
        """Get the full viewer URL."""
        params = [
            f"remote={self.server_url}",
            f"jupyter=true",
            f"viewerId={self._viewer_id}",
        ]
        return f"https://www.cellucid.com?{'&'.join(params)}"

    def display(self):
        """Display the viewer in a notebook cell."""
        if not self._context["in_jupyter"]:
            print(f"Not in Jupyter environment. Open manually: {self.viewer_url}")
            return

        from IPython.display import display, HTML, IFrame

        if self._context["notebook_type"] == "colab":
            # Google Colab needs special handling
            html = self._generate_colab_html()
            display(HTML(html))
        else:
            # Standard Jupyter/JupyterLab - use IFrame
            display(IFrame(src=self.viewer_url, width="100%", height=self.height))

        self._displayed = True

    def _generate_colab_html(self) -> str:
        """Generate HTML for Google Colab embedding."""
        return f"""
        <div id="cellucid-viewer-{self._viewer_id}" style="width:100%; height:{self.height}px;">
            <iframe
                src="{self.viewer_url}"
                width="100%"
                height="100%"
                frameborder="0"
                allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope"
                allowfullscreen>
            </iframe>
        </div>
        <script>
        // Set up message handler for communication with viewer
        window.addEventListener('message', function(event) {{
            if (event.origin !== 'https://www.cellucid.com') return;
            if (event.data.viewerId !== '{self._viewer_id}') return;

            // Forward to Python via Google Colab's callback mechanism
            if (window.google && window.google.colab) {{
                google.colab.kernel.invokeFunction(
                    'cellucid_message_handler',
                    [event.data],
                    {{}}
                );
            }}
        }});
        </script>
        """

    def _generate_interactive_js(self) -> str:
        """Generate JavaScript for interactive communication."""
        return f"""
        <script>
        (function() {{
            const viewerId = '{self._viewer_id}';
            const serverUrl = '{self.server_url}';

            // Set up message channel with the viewer iframe
            window.cellucidViewers = window.cellucidViewers || {{}};
            window.cellucidViewers[viewerId] = {{
                sendMessage: function(msg) {{
                    const iframe = document.querySelector('#cellucid-viewer-' + viewerId + ' iframe');
                    if (iframe) {{
                        iframe.contentWindow.postMessage(
                            {{...msg, viewerId: viewerId}},
                            'https://www.cellucid.com'
                        );
                    }}
                }}
            }};

            // Listen for messages from the viewer
            window.addEventListener('message', function(event) {{
                if (event.origin !== 'https://www.cellucid.com') return;
                if (event.data.viewerId !== viewerId) return;

                // Handle different message types
                switch (event.data.type) {{
                    case 'selection':
                        console.log('Selected cells:', event.data.cells);
                        break;
                    case 'ready':
                        console.log('Cellucid viewer ready');
                        break;
                }}
            }});
        }})();
        </script>
        """

    def send_message(self, message: dict):
        """
        Send a message to the viewer.

        Args:
            message: Message dictionary to send
        """
        if not self._displayed:
            logger.warning("Viewer not yet displayed")
            return

        from IPython.display import display, Javascript

        js = f"""
        if (window.cellucidViewers && window.cellucidViewers['{self._viewer_id}']) {{
            window.cellucidViewers['{self._viewer_id}'].sendMessage({json.dumps(message)});
        }}
        """
        display(Javascript(js))

    def highlight_cells(self, cell_indices: list[int], color: str = "#ff0000"):
        """
        Highlight specific cells in the viewer.

        Args:
            cell_indices: List of cell indices to highlight
            color: Highlight color (hex or CSS color name)
        """
        self.send_message({
            "type": "highlight",
            "cells": cell_indices,
            "color": color,
        })

    def clear_highlights(self):
        """Clear all cell highlights."""
        self.send_message({"type": "clearHighlights"})

    def set_color_by(self, field: str):
        """
        Set the coloring field.

        Args:
            field: Field name to color by (e.g., 'cluster', 'cell_type')
        """
        self.send_message({
            "type": "setColorBy",
            "field": field,
        })

    def set_visibility(self, cell_indices: list[int] | None = None, visible: bool = True):
        """
        Set visibility of cells.

        Args:
            cell_indices: Cell indices to modify (None = all cells)
            visible: Whether cells should be visible
        """
        self.send_message({
            "type": "setVisibility",
            "cells": cell_indices,
            "visible": visible,
        })

    def reset_view(self):
        """Reset the camera to the default view."""
        self.send_message({"type": "resetCamera"})

    def stop(self):
        """Stop the data server and cleanup."""
        if self._server:
            self._server.stop()
            self._server = None
        logger.info("Cellucid viewer stopped")

    def __del__(self):
        """Cleanup on deletion."""
        self.stop()

    def __repr__(self) -> str:
        status = "running" if self._server and self._server.is_running() else "stopped"
        return f"CellucidViewer('{self.data_dir}', port={self.port}, status={status})"

    def _repr_html_(self) -> str:
        """HTML representation for Jupyter display."""
        if self._context["in_jupyter"]:
            return f'<a href="{self.viewer_url}" target="_blank">Open Cellucid Viewer</a>'
        return repr(self)


def show(
    data_dir: str | Path,
    height: int = 600,
    interactive: bool = True,
) -> CellucidViewer:
    """
    Quick function to display a cellucid dataset in a notebook.

    Args:
        data_dir: Path to the dataset directory
        height: Height of the viewer in pixels
        interactive: Enable bidirectional communication

    Returns:
        CellucidViewer instance for further interaction

    Example:
        >>> from cellucid.jupyter import show
        >>> viewer = show("/path/to/my_dataset")
        >>> viewer.highlight_cells([1, 2, 3])
    """
    return CellucidViewer(
        data_dir=data_dir,
        height=height,
        interactive=interactive,
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
