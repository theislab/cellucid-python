"""
The prepared-data server and its one-call entry point.

:class:`CellucidServer` validates the served root, builds the artifact
inventory for every dataset it found, binds the socket, and runs the handler
either in the foreground or on a background thread. :func:`serve` is the
one-call form used by ``cellucid serve`` and by :func:`cellucid.serve`.
"""

from __future__ import annotations

import threading
import webbrowser
from collections.abc import Sequence
from functools import partial
from http.server import HTTPServer
from pathlib import Path

from .._console import console_print
from .._server_base import (
    CELLUCID_WEB_URL,
    DEFAULT_HOST,
    DEFAULT_PORT,
    _web_cache_dir,
    print_detail,
    print_server_banner,
    print_step,
    print_success,
    require_allowed_hosts,
    require_server_port,
)
from ._artifacts import _build_prepared_artifact_inventory, _read_json_object
from ._datasets import _list_exported_datasets
from ._handler import CORSRequestHandler, logger


class CellucidServer:
    """
    Cellucid data server for serving datasets over HTTP.

    Supports multiple deployment modes:
    - Local: Direct browser access on localhost
    - SSH tunnel: Access via port forwarding from remote server
    - Jupyter: Embedded in notebook environment

    Example:
        server = CellucidServer("/path/to/data")
        server.start()  # Blocking

        # Or non-blocking:
        server.start_background()
        # ... do other things ...
        server.stop()
    """

    def __init__(
        self,
        data_dir: str | Path,
        port: int = DEFAULT_PORT,
        host: str = DEFAULT_HOST,
        open_browser: bool = False,
        quiet: bool = False,
        *,
        serve_web_ui: bool = True,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
        allowed_hosts: Sequence[str] | None = None,
    ):
        """
        Initialize the server.

        Args:
            data_dir: Path to the dataset directory (single dataset or
                multi-dataset). A leading ``~`` is expanded to the user's home
                directory and the path is resolved to an absolute path.
            port: Port to serve on (default: 8765)
            host: Host to bind to (default: 127.0.0.1 for localhost only)
            open_browser: Whether to open the browser on start
            quiet: Suppress info messages
            serve_web_ui: Establish and serve the exact current web build.
            web_source_url: Origin publishing the web asset inventory.
            web_cache_dir: Directory holding the active verified web build.
            allowed_hosts: Extra ``Host`` names this server answers to, needed
                only behind a reverse proxy such as ``jupyter-server-proxy``.
                Each entry is one bare DNS name or IP literal — no port, no
                scheme, no wildcard — and is matched on any port. Declaring
                names also applies the ``Host`` check to a non-loopback bind.
        """
        self.data_dir = Path(data_dir).expanduser().resolve()
        self.port = require_server_port(port)
        self.host = host
        self.allowed_hosts = require_allowed_hosts(allowed_hosts)
        if type(open_browser) is not bool:
            raise TypeError("open_browser must be a boolean")
        if type(quiet) is not bool:
            raise TypeError("quiet must be a boolean")
        self.open_browser = open_browser
        self.quiet = quiet
        if type(serve_web_ui) is not bool:
            raise TypeError("serve_web_ui must be a boolean")
        self.serve_web_ui = serve_web_ui
        self.web_source_url = web_source_url
        self.web_cache_dir = (
            Path(web_cache_dir).expanduser().resolve()
            if web_cache_dir is not None
            else _web_cache_dir()
        )

        # Step 1: Validate dataset
        if not quiet:
            print_step(1, 3, "Validating dataset")
            print_detail("Path", str(self.data_dir))

        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {self.data_dir}")
        if not self.data_dir.is_dir():
            raise NotADirectoryError(f"Prepared-data path must be a directory: {self.data_dir}")

        self._datasets = _list_exported_datasets(self.data_dir)
        if not self._datasets:
            raise ValueError("Prepared-data path does not contain one complete current dataset.")
        self._artifact_inventory = _build_prepared_artifact_inventory(
            self.data_dir,
            self._datasets,
        )

        if not quiet:
            print_success("Dataset valid")

        # Step 2: Load dataset info
        if not quiet:
            print_step(2, 3, "Loading dataset info")
            self._print_dataset_info()
            print_success("Dataset loaded")

        self._server: HTTPServer | None = None
        self._thread: threading.Thread | None = None
        self._running = False
        self._started = False
        self._closed = False
        self._serving = False
        self._background_error: BaseException | None = None

        from .. import __version__

        self.server_info = {
            "version": __version__,
            "host": self.host,
            "port": self.port,
            "mode": "standalone",
        }

    def _print_dataset_info(self):
        """Print information about the dataset."""
        if len(self._datasets) != 1 or self._datasets[0]["path"] != "/":
            print_detail("Datasets", str(len(self._datasets)))
            return
        identity = _read_json_object(
            self.data_dir / "dataset_identity.json",
            label="dataset_identity.json",
        )
        stats = identity.get("stats")
        if not isinstance(stats, dict):
            raise TypeError("dataset_identity.json stats must be a JSON object.")
        n_cells = stats.get("n_cells")
        n_genes = stats.get("n_genes")
        has_connectivity = stats.get("has_connectivity")
        if type(n_cells) is not int or n_cells < 0:
            raise ValueError("dataset_identity.json stats.n_cells must be non-negative.")
        if type(n_genes) is not int or n_genes < 0:
            raise ValueError("dataset_identity.json stats.n_genes must be non-negative.")
        if type(has_connectivity) is not bool:
            raise ValueError("dataset_identity.json stats.has_connectivity must be a boolean.")
        print_detail("Cells", f"{n_cells:,}")
        print_detail("Genes", f"{n_genes:,}")
        print_detail("Connectivity", "yes" if has_connectivity else "no")

    @property
    def url(self) -> str:
        """Get the URL of the currently running server."""
        if not self._running or self._server is None:
            raise RuntimeError(
                "CellucidServer URL is unavailable because the server is not running"
            )
        return f"http://{self.host}:{self.port}"

    @property
    def viewer_url(self) -> str:
        """Open this server's prepared-data catalog in the verified viewer.

        The viewer selects the sole declared dataset when the catalog is
        unique. A multi-dataset catalog requires an exact dataset selection;
        this URL never embeds or guesses an arbitrary catalog entry.
        """
        return f"{self.url}/?source=remote"

    def start(self, blocking: bool = True):
        """Start this single-use server."""
        if type(blocking) is not bool:
            raise TypeError("blocking must be a boolean")
        if self._running:
            raise RuntimeError("Server is already running.")
        if self._started or self._closed:
            raise RuntimeError(
                "CellucidServer is single-use and has been closed. Create a new server instance."
            )
        self._started = True

        try:
            # Step 3: Start server
            if not self.quiet:
                print_step(3, 3, "Starting server")

            if self.serve_web_ui:
                from ..web_cache import ensure_web_ui_cached

                if not self.quiet:
                    print_detail(
                        "Viewer UI generation",
                        "establishing exact configured source",
                    )
                ensure_web_ui_cached(
                    cache_dir=self.web_cache_dir,
                    source_url=self.web_source_url,
                    force=True,
                    show_progress=not self.quiet,
                )
                if not self.quiet:
                    print_success("Viewer UI generation established")

            handler = partial(
                CORSRequestHandler,
                data_dir=self.data_dir,
                server_info=self.server_info,
                datasets=self._datasets,
                artifact_inventory=self._artifact_inventory,
                serve_web_ui=self.serve_web_ui,
                web_cache_dir=self.web_cache_dir,
                allowed_hosts=self.allowed_hosts,
            )

            self._server = HTTPServer((self.host, self.port), handler)
            self.port = require_server_port(self._server.server_address[1])
            self.server_info["port"] = self.port
            self._running = True

            if not self.quiet:
                print_success("Server ready")
                print_server_banner(self.url, self.viewer_url)

            if blocking:
                if self.open_browser and webbrowser.open(self.viewer_url) is not True:
                    raise RuntimeError(f"Could not open the browser for {self.viewer_url}")
                self._serving = True
                try:
                    self._server.serve_forever()
                finally:
                    self._serving = False
                self._finish_serving()
            else:
                self._serve_entered = threading.Event()
                self._thread = threading.Thread(
                    target=self._serve_in_background,
                    daemon=True,
                )
                self._thread.start()
                self._serve_entered.wait()
                if self.open_browser and webbrowser.open(self.viewer_url) is not True:
                    raise RuntimeError(f"Could not open the browser for {self.viewer_url}")
        except BaseException:
            self._rollback_failed_start(shutdown=self._thread is not None)
            raise

    def _serve_in_background(self) -> None:
        """Run the bound server and retain an exact asynchronous failure."""
        self._serving = True
        self._serve_entered.set()
        serving_error: BaseException | None = None
        try:
            server = self._server
            if server is None:
                raise RuntimeError("CellucidServer lost its bound HTTP server before serving")
            server.serve_forever()
        except BaseException as error:
            serving_error = error
        finally:
            self._serving = False
            try:
                self._finish_serving()
            except BaseException as cleanup_error:
                if serving_error is None:
                    serving_error = cleanup_error
                else:
                    logger.exception(
                        "Prepared-data server cleanup also failed after its serving loop failed"
                    )
        self._background_error = serving_error
        if serving_error is not None:
            logger.error(
                "Prepared-data background server failed",
                exc_info=(
                    type(serving_error),
                    serving_error,
                    serving_error.__traceback__,
                ),
            )

    def _finish_serving(self) -> None:
        """Close the socket after the serving loop has ended."""
        self._running = False
        server = self._server
        self._server = None
        self._closed = True
        if server is not None:
            try:
                server.server_close()
            except BaseException as error:
                raise RuntimeError(
                    f"Prepared-data server cleanup failed: {type(error).__name__}: {error}"
                ) from error

    def _rollback_failed_start(self, *, shutdown: bool) -> None:
        """Release acquired resources without replacing the startup exception."""
        self._running = False
        server = self._server
        if server is not None and shutdown and self._serving:
            try:
                server.shutdown()
            except BaseException:
                logger.exception(
                    "Failed to shut down the prepared-data server after startup failed"
                )
        if server is not None:
            try:
                server.server_close()
            except BaseException:
                logger.exception("Failed to close the prepared-data socket after startup failed")
        thread = self._thread
        if thread is not None and thread is not threading.current_thread() and thread.is_alive():
            thread.join()
        self._server = None
        self._thread = None
        self._serving = False
        self._closed = True

    def start_background(self):
        """Start the server in a background thread."""
        self.start(blocking=False)

    def stop(self):
        """Stop this server and release its socket."""
        self._running = False
        failures: list[BaseException] = []
        server = self._server
        thread = self._thread

        if server is not None and self._serving:
            try:
                server.shutdown()
            except BaseException as error:
                failures.append(error)
        if server is not None:
            try:
                server.server_close()
            except BaseException as error:
                failures.append(error)
        self._server = None

        if thread is not None and thread is not threading.current_thread() and thread.is_alive():
            thread.join()
        self._thread = None
        self._serving = False
        self._closed = True

        if not self.quiet:
            console_print("Server stopped")
        if failures:
            details = "; ".join(f"{type(error).__name__}: {error}" for error in failures)
            raise RuntimeError(f"Prepared-data server shutdown failed: {details}") from failures[0]

    def is_running(self) -> bool:
        """Check if the server is running."""
        return self._running

    def wait(self):
        """Wait for the background server to stop."""
        thread = self._thread
        if thread is not None:
            try:
                thread.join()
            except BaseException:
                self.stop()
                raise
        if self._background_error is not None:
            raise self._background_error


def serve(
    data_dir: str | Path,
    port: int = DEFAULT_PORT,
    host: str = DEFAULT_HOST,
    open_browser: bool = True,
    quiet: bool = False,
    *,
    serve_web_ui: bool = True,
    web_source_url: str = CELLUCID_WEB_URL,
    web_cache_dir: str | Path | None = None,
    allowed_hosts: Sequence[str] | None = None,
):
    """
    Serve a cellucid dataset directory.

    This is the main entry point for serving data. It starts an HTTP server
    that serves the dataset files with CORS headers enabled.

    Args:
        data_dir: Path to the dataset directory. A leading ``~`` is expanded to
            the user's home directory and the path is resolved to an absolute
            path.
        port: Port to serve on (default: 8765)
        host: Host to bind to (default: 127.0.0.1)
        open_browser: Whether to open the viewer in browser (default: True)
        quiet: Suppress info messages
        serve_web_ui: Establish and serve the exact current web build.
        web_source_url: Origin publishing the web asset inventory.
        web_cache_dir: Directory holding the active verified web build.
        allowed_hosts: Extra ``Host`` names this server answers to, needed only
            behind a reverse proxy such as ``jupyter-server-proxy``. Each entry
            is one bare DNS name or IP literal — no port, no scheme, no
            wildcard — and is matched on any port.

    Example:
        >>> from cellucid import serve
        >>> serve("/path/to/my_dataset")

        # For remote server access via SSH:
        >>> serve("/path/to/data", host="0.0.0.0")
        # Then on local machine: ssh -L 8765:localhost:8765 remote-server
    """
    server = CellucidServer(
        data_dir=data_dir,
        port=port,
        host=host,
        open_browser=open_browser,
        quiet=quiet,
        serve_web_ui=serve_web_ui,
        web_source_url=web_source_url,
        web_cache_dir=web_cache_dir,
        allowed_hosts=allowed_hosts,
    )
    server.start()
