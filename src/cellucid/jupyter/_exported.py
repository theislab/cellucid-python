"""
The notebook viewer for a directory that :func:`cellucid.prepare` already wrote.

:class:`CellucidViewer` starts a :class:`cellucid.server.CellucidServer` over an
export directory and embeds the viewer against it; :func:`show` is the one-call
form.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from .._server_base import CELLUCID_WEB_URL
from ..server import CellucidServer
from ._base import BaseViewer, logger


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
        allowed_hosts: Sequence[str] | None = None,
    ):
        """
        Initialize the viewer.

        Args:
            data_dir: Path to the cellucid dataset directory. A leading ``~`` is
                expanded to the user's home directory and the path is resolved
                to an absolute path.
            port: Port for the data server (auto-selected if None)
            height: Height of the embedded viewer in pixels
            auto_open: Automatically display when created
            client_server_url: Exact browser-reachable data server base URL.
            web_source_url: Origin publishing the web asset inventory.
            web_cache_dir: Directory holding the active verified web build.
            allowed_hosts: Extra ``Host`` names the data server answers to.
                Required when ``client_server_url`` points at a reverse proxy
                such as ``jupyter-server-proxy``: name that proxy's host here.
        """
        super().__init__(
            port=port,
            height=height,
            client_server_url=client_server_url,
            web_source_url=web_source_url,
            web_cache_dir=web_cache_dir,
            allowed_hosts=allowed_hosts,
        )

        try:
            self.data_dir = Path(data_dir).expanduser().resolve()

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
            allowed_hosts=sorted(self.allowed_hosts),
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
    allowed_hosts: Sequence[str] | None = None,
) -> CellucidViewer:
    """
    Quick function to display a cellucid dataset in a notebook.

    Args:
        data_dir: Path to the dataset directory. A leading ``~`` is expanded to
            the user's home directory and the path is resolved to an absolute
            path.
        height: Height of the viewer in pixels
        client_server_url: Exact browser-reachable data server base URL.
        web_source_url: Origin publishing the web asset inventory.
        web_cache_dir: Directory holding the active verified web build.
        allowed_hosts: Extra ``Host`` names the data server answers to.
            Required when ``client_server_url`` points at a reverse proxy such
            as ``jupyter-server-proxy``: name that proxy's host here.

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
        allowed_hosts=allowed_hosts,
    )
