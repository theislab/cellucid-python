"""
Cellucid Data Server

A lightweight HTTP server for serving pre-exported cellucid datasets.
Supports both local and remote access patterns:

1. Local mode: Run on your machine, open browser locally
2. SSH tunnel mode: Run on remote server, access via SSH port forwarding
3. Jupyter mode: Run alongside Jupyter, embed in notebooks

Usage:
    from cellucid import serve
    serve("/path/to/data", port=8765)

Or via CLI:
    cellucid serve /path/to/data --port 8765

The server provides:
- Static file serving for dataset files
- CORS headers for cross-origin access (needed for web viewer)
- Health check endpoint for connection validation

For serving AnnData directly (without pre-export), use:
    from cellucid import serve_anndata
    serve_anndata("/path/to/data.h5ad")
"""

from __future__ import annotations

import json
import logging
import threading
import webbrowser
from functools import partial
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from urllib.parse import unquote, urlparse

logger = logging.getLogger("cellucid.server")

# Import shared configuration from _server_base to avoid duplication
from ._server_base import (
    CORSMixin,
    DEFAULT_PORT,
    DEFAULT_HOST,
    CELLUCID_WEB_URL,
    ensure_port_available,
)


class CORSRequestHandler(CORSMixin, SimpleHTTPRequestHandler):
    """HTTP handler with CORS support for serving dataset files."""

    allow_caching = True  # Static files can be cached

    def __init__(self, *args, data_dir: Path, server_info: dict, **kwargs):
        self.data_dir = data_dir
        self.server_info = server_info
        # Must call super().__init__ last because it calls do_GET immediately
        super().__init__(*args, directory=str(data_dir), **kwargs)

    def end_headers(self):
        """Add CORS headers to every response."""
        self.add_cors_headers()
        super().end_headers()

    def do_POST(self):
        """Handle POST requests (events from frontend)."""
        if self.handle_event_post():
            return
        # No other POST endpoints - return 404
        self.send_error_response(404, f"POST not supported for path: {self.path}")

    def do_GET(self):
        """Handle GET requests with special endpoints."""
        parsed = urlparse(self.path)
        path = unquote(parsed.path)

        # Health check endpoint
        if path == "/_cellucid/health":
            self.send_json({
                "status": "ok",
                "type": "exported",
                "version": self.server_info.get("version", "unknown"),
                "data_dir": str(self.data_dir),
            })
            return

        # Server info endpoint
        if path == "/_cellucid/info":
            self.send_json(self.server_info)
            return

        # Datasets list endpoint
        if path == "/_cellucid/datasets":
            datasets = self._list_datasets()
            self.send_json({"datasets": datasets})
            return

        # Regular file serving
        super().do_GET()

    def _list_datasets(self) -> list[dict]:
        """List available datasets in the data directory."""
        datasets = []

        # Check if data_dir itself is a dataset
        if self._is_dataset_dir(self.data_dir):
            datasets.append({
                "id": self.data_dir.name,
                "path": "/",
                "name": self._get_dataset_name(self.data_dir),
            })
        else:
            # Look for subdirectories that are datasets
            for subdir in self.data_dir.iterdir():
                if subdir.is_dir() and self._is_dataset_dir(subdir):
                    datasets.append({
                        "id": subdir.name,
                        "path": f"/{subdir.name}/",
                        "name": self._get_dataset_name(subdir),
                    })

        return datasets

    def _is_dataset_dir(self, path: Path) -> bool:
        """Check if a directory is a valid cellucid dataset."""
        # Must have obs_manifest.json
        if not (path / "obs_manifest.json").exists():
            return False
        # Must have at least one points file
        for dim in ["1d", "2d", "3d", "4d"]:
            if (path / f"points_{dim}.bin").exists():
                return True
            if (path / f"points_{dim}.bin.gz").exists():
                return True
        return False

    def _get_dataset_name(self, path: Path) -> str:
        """Get the display name for a dataset."""
        identity_file = path / "dataset_identity.json"
        if identity_file.exists():
            try:
                with open(identity_file) as f:
                    identity = json.load(f)
                    return identity.get("name", path.name)
            except Exception:
                pass
        return path.name

    def log_message(self, format: str, *args):
        """Override to use Python logging instead of stderr."""
        logger.debug("%s - %s", self.address_string(), format % args)


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
    ):
        """
        Initialize the server.

        Args:
            data_dir: Path to the dataset directory (single dataset or multi-dataset)
            port: Port to serve on (default: 8765)
            host: Host to bind to (default: 127.0.0.1 for localhost only)
            open_browser: Whether to open the browser on start
            quiet: Suppress info messages
        """
        self.data_dir = Path(data_dir).resolve()
        self.port = port
        self.host = host
        self.open_browser = open_browser
        self.quiet = quiet

        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {self.data_dir}")

        self._server: HTTPServer | None = None
        self._thread: threading.Thread | None = None
        self._ws_server = None
        self._ws_task = None
        self._running = False

        # Version from package
        try:
            from . import __version__
        except ImportError:
            __version__ = "0.0.0"

        self.server_info = {
            "version": __version__,
            "data_dir": str(self.data_dir),
            "host": self.host,
            "port": self.port,
            "mode": "standalone",
        }

    @property
    def url(self) -> str:
        """Get the server URL."""
        return f"http://{self.host}:{self.port}"

    @property
    def viewer_url(self) -> str:
        """Get the full URL to open the viewer with this server's data."""
        return f"{CELLUCID_WEB_URL}?remote={self.url}"

    def start(self, blocking: bool = True):
        """
        Start the server.

        Args:
            blocking: If True, block until interrupted. If False, start in background.
        """
        if self._running:
            logger.warning("Server is already running")
            return

        # Ensure port is available (finds new one if needed)
        self.port = ensure_port_available(self.host, self.port, self.quiet)
        self.server_info["port"] = self.port

        # Create handler with data directory
        handler = partial(
            CORSRequestHandler,
            data_dir=self.data_dir,
            server_info=self.server_info,
        )

        self._server = HTTPServer((self.host, self.port), handler)
        self._running = True

        if not self.quiet:
            print(f"\n{'=' * 60}")
            print(f"  Cellucid Data Server")
            print(f"{'=' * 60}")
            print(f"  Data directory: {self.data_dir}")
            print(f"  Server URL:     {self.url}")
            print(f"")
            print(f"  Open viewer at:")
            print(f"    {self.viewer_url}")
            print(f"")
            print(f"  For SSH tunnel access, run on your local machine:")
            print(f"    ssh -L {self.port}:localhost:{self.port} <remote-host>")
            print(f"")
            print(f"  Press Ctrl+C to stop the server")
            print(f"{'=' * 60}\n")

        if self.open_browser:
            webbrowser.open(self.viewer_url)

        if blocking:
            try:
                self._server.serve_forever()
            except KeyboardInterrupt:
                if not self.quiet:
                    print("\nShutting down server...")
                self.stop()
        else:
            self._thread = threading.Thread(target=self._server.serve_forever, daemon=True)
            self._thread.start()

    def start_background(self):
        """Start the server in a background thread."""
        self.start(blocking=False)

    def stop(self):
        """Stop the server."""
        if self._server:
            self._server.shutdown()
            self._server = None
        self._running = False
        if not self.quiet:
            logger.info("Server stopped")

    def is_running(self) -> bool:
        """Check if the server is running."""
        return self._running

    def wait(self):
        """Wait for the server to stop (blocks until Ctrl+C or stop())."""
        if self._thread:
            try:
                while self._running:
                    self._thread.join(timeout=1)
            except KeyboardInterrupt:
                self.stop()


def serve(
    data_dir: str | Path,
    port: int = DEFAULT_PORT,
    host: str = DEFAULT_HOST,
    open_browser: bool = True,
    quiet: bool = False,
):
    """
    Serve a cellucid dataset directory.

    This is the main entry point for serving data. It starts an HTTP server
    that serves the dataset files with CORS headers enabled.

    Args:
        data_dir: Path to the dataset directory
        port: Port to serve on (default: 8765)
        host: Host to bind to (default: 127.0.0.1)
        open_browser: Whether to open the viewer in browser (default: True)
        quiet: Suppress info messages

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
    )
    server.start()


def main():
    """CLI entry point for cellucid serve (legacy, kept for direct imports)."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Serve a cellucid dataset directory",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Serve a local dataset
    cellucid serve /path/to/my_dataset

    # Serve on a different port
    cellucid serve /path/to/data --port 9000

    # Serve on all interfaces (for remote access)
    cellucid serve /path/to/data --host 0.0.0.0

    # For SSH tunnel access from remote server:
    # On the server: cellucid serve /path/to/data
    # On local machine: ssh -L 8765:localhost:8765 user@server
    # Then open: https://www.cellucid.com?remote=http://localhost:8765
""",
    )

    parser.add_argument(
        "data_dir",
        type=str,
        help="Path to the dataset directory",
    )
    parser.add_argument(
        "--port", "-p",
        type=int,
        default=DEFAULT_PORT,
        help=f"Port to serve on (default: {DEFAULT_PORT})",
    )
    parser.add_argument(
        "--host", "-H",
        type=str,
        default=DEFAULT_HOST,
        help=f"Host to bind to (default: {DEFAULT_HOST})",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Don't open browser automatically",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress info messages",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    # Configure logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    elif not args.quiet:
        logging.basicConfig(level=logging.INFO)

    serve(
        data_dir=args.data_dir,
        port=args.port,
        host=args.host,
        open_browser=not args.no_browser,
        quiet=args.quiet,
    )


if __name__ == "__main__":
    main()
