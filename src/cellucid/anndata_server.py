"""
AnnData Server for Cellucid

HTTP server that serves AnnData data dynamically in the Cellucid format.
This allows direct visualization of AnnData files without pre-export.

Supports:
- h5ad files: Single HDF5 files (backed mode for lazy loading)
- zarr stores: Directory-based stores (inherently lazy-loaded)
- In-memory AnnData objects

Usage:
    from cellucid.anndata_server import serve_anndata

    # Serve an h5ad file
    serve_anndata("/path/to/data.h5ad")

    # Serve a zarr store (must be a directory)
    serve_anndata("/path/to/data.zarr")

    # Or an in-memory AnnData
    serve_anndata(adata)
"""

from __future__ import annotations

import logging
import re
import socket
import threading
import webbrowser
from functools import partial
from http import HTTPStatus
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union
from urllib.parse import parse_qs, unquote, urlparse

from .anndata_adapter import AnnDataAdapter, _safe_filename_component
from ._server_base import (
    CORSMixin,
    DEFAULT_PORT,
    DEFAULT_HOST,
    CELLUCID_WEB_URL,
    ensure_port_available,
)

if TYPE_CHECKING:
    import anndata

logger = logging.getLogger("cellucid.anndata_server")


class AnnDataRequestHandler(CORSMixin, SimpleHTTPRequestHandler):
    """
    HTTP handler for serving AnnData data in Cellucid format.

    Routes:
        /dataset_identity.json - Dataset metadata
        /obs_manifest.json - Cell metadata manifest
        /var_manifest.json - Gene expression manifest
        /connectivity_manifest.json - Connectivity manifest
        /points_{dim}d.bin - Embedding coordinates
        /obs/{field}.values.f32.bin - Continuous obs field
        /obs/{field}.codes.{ext}.bin - Categorical obs field codes
        /obs/{field}.outliers.f32.bin - Categorical outlier quantiles
        /var/{gene}.values.f32.bin - Gene expression values
        /connectivity/edges.src.bin - Edge sources
        /connectivity/edges.dst.bin - Edge destinations
        /_cellucid/health - Health check
        /_cellucid/info - Server info
        /_cellucid/events - POST endpoint for frontend events (hooks)
    """

    allow_caching = False  # Data is dynamic, don't cache

    def __init__(self, *args, adapter: AnnDataAdapter, server_info: dict, **kwargs):
        self.adapter = adapter
        self.server_info = server_info
        # Don't call super().__init__ with directory since we're serving virtual files
        super(SimpleHTTPRequestHandler, self).__init__(*args, **kwargs)

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

    def do_HEAD(self):
        """Handle HEAD requests (for urlExists checks)."""
        self.do_GET(head_only=True)

    def do_GET(self, head_only: bool = False):
        """Handle GET requests."""
        parsed = urlparse(self.path)
        path = unquote(parsed.path).lstrip("/")

        # Check for gzip support
        accept_encoding = self.headers.get("Accept-Encoding", "")
        supports_gzip = "gzip" in accept_encoding

        try:
            if path == "_cellucid/health":
                self.send_json({
                    "status": "ok",
                    "type": "anndata",
                    "version": self.server_info.get("version", "unknown"),
                    "n_cells": self.adapter.n_cells,
                    "n_genes": self.adapter.n_genes,
                }, head_only)
            elif path == "_cellucid/info":
                self.send_json(self.server_info, head_only)
            elif path == "dataset_identity.json":
                self.send_json(self.adapter.get_dataset_identity(), head_only)
            elif path == "obs_manifest.json":
                self.send_json(self.adapter.get_obs_manifest(), head_only)
            elif path == "var_manifest.json":
                self.send_json(self.adapter.get_var_manifest(), head_only)
            elif path == "connectivity_manifest.json":
                manifest = self.adapter.get_connectivity_manifest()
                if manifest:
                    self.send_json(manifest, head_only)
                else:
                    self.send_error_response(404, "No connectivity data")
            elif path.startswith("points_") and path.endswith(".bin"):
                self._handle_points(path, head_only, supports_gzip)
            elif path.startswith("obs/"):
                self._handle_obs(path, head_only, supports_gzip)
            elif path.startswith("var/"):
                self._handle_var(path, head_only, supports_gzip)
            elif path.startswith("connectivity/"):
                self._handle_connectivity(path, head_only, supports_gzip)
            else:
                self.send_error_response(404, f"Not found: {path}")

        except Exception as e:
            logger.error(f"Error handling {path}: {e}", exc_info=True)
            self.send_error_response(500, str(e))

    def _handle_points(self, path: str, head_only: bool, supports_gzip: bool):
        """Handle points_Xd.bin requests."""
        # Parse dimension from filename: points_1d.bin, points_2d.bin, points_3d.bin
        # Also handle .gz suffix for explicit compressed requests
        clean_path = path[:-3] if path.endswith(".gz") else path
        match = re.match(r"points_(\d)d\.bin", clean_path)
        if not match:
            self.send_error_response(404, f"Invalid points path: {path}")
            return

        dim = int(match.group(1))
        try:
            # Compress if client supports gzip (transparent compression for better network perf)
            compress = supports_gzip
            data = self.adapter.get_points_binary(dim, compress=compress)
            self.send_binary(data, head_only=head_only, compressed=compress)
        except ValueError as e:
            self.send_error_response(404, str(e))

    def _handle_obs(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle obs field requests.

        Supports paths in both formats for compatibility:
        - obs/{field}.values.f32 (prepare format)
        - obs/{field}.values.f32.bin (legacy format)
        - obs/{field}.codes.u8 or obs/{field}.codes.u16
        - obs/{field}.outliers.f32
        """
        # Remove 'obs/' prefix
        filename = path[4:]

        # Strip suffixes for parsing
        if filename.endswith(".gz"):
            filename = filename[:-3]
        if filename.endswith(".bin"):
            filename = filename[:-4]

        # Compress if client supports gzip (transparent compression)
        compress = supports_gzip

        # Parse filename: {field}.{type}.{dtype}
        # Examples: cell_type.codes.u8, n_counts.values.f32, cluster.outliers.f32
        parts = filename.rsplit(".", 2)  # Split from right: [field, type, dtype]
        if len(parts) < 3:
            self.send_error_response(404, f"Invalid obs path format: {path}")
            return

        field_name = parts[0]
        data_type = parts[1]  # 'values', 'codes', or 'outliers'

        # Find the actual field key (may differ from safe filename)
        obs_keys = self.adapter.get_obs_keys()
        actual_key = None
        for key in obs_keys:
            if _safe_filename_component(key) == field_name:
                actual_key = key
                break

        if actual_key is None:
            self.send_error_response(404, f"Obs field not found: {field_name}")
            return

        try:
            if data_type == "values":
                # Continuous field values
                data = self.adapter.get_obs_continuous_values(actual_key, compress=compress)
                self.send_binary(data, head_only=head_only, compressed=compress)
            elif data_type == "codes":
                # Categorical codes
                data, categories, missing = self.adapter.get_obs_categorical_codes(actual_key, compress=compress)
                self.send_binary(data, head_only=head_only, compressed=compress)
            elif data_type == "outliers":
                # Outlier quantiles for categorical field
                data = self.adapter.get_obs_outlier_quantiles(actual_key, compress=compress)
                self.send_binary(data, head_only=head_only, compressed=compress)
            else:
                self.send_error_response(404, f"Unknown obs data type: {data_type} (expected values, codes, or outliers)")
        except Exception as e:
            self.send_error_response(500, str(e))

    def _handle_var(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle var (gene expression) requests.

        Supports paths in both formats for compatibility:
        - var/{gene}.values.f32 (prepare format)
        - var/{gene}.values.f32.bin (legacy format)
        """
        # Remove 'var/' prefix
        filename = path[4:]

        # Strip suffixes for parsing
        if filename.endswith(".gz"):
            filename = filename[:-3]
        if filename.endswith(".bin"):
            filename = filename[:-4]

        # Compress if client supports gzip (transparent compression)
        compress = supports_gzip

        # Check format: {gene}.values.f32
        if not filename.endswith(".values.f32"):
            self.send_error_response(404, f"Invalid var path format: {path} (expected {{gene}}.values.f32)")
            return

        gene_safe = filename[:-11]  # Remove '.values.f32'

        # Find actual gene ID (may differ from safe filename due to special characters)
        gene_ids = self.adapter.get_gene_ids()
        actual_gene = None

        # First try safe filename match
        for gid in gene_ids:
            if _safe_filename_component(gid) == gene_safe:
                actual_gene = gid
                break

        # Fallback: try exact match (in case gene ID is already safe)
        if actual_gene is None and gene_safe in gene_ids:
            actual_gene = gene_safe

        if actual_gene is None:
            self.send_error_response(404, f"Gene not found: {gene_safe}")
            return

        try:
            data = self.adapter.get_gene_expression(actual_gene, compress=compress)
            self.send_binary(data, head_only=head_only, compressed=compress)
        except KeyError as e:
            self.send_error_response(404, str(e))
        except Exception as e:
            self.send_error_response(500, str(e))

    def _handle_connectivity(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle connectivity edge requests.

        Supported paths:
        - connectivity/edges.src.bin - Edge source indices
        - connectivity/edges.dst.bin - Edge destination indices
        """
        filename = path.split("/")[-1]

        # Strip suffixes for parsing
        if filename.endswith(".gz"):
            filename = filename[:-3]

        # Compress if client supports gzip (transparent compression)
        compress = supports_gzip

        try:
            sources_data, dests_data, n_edges, max_neighbors = self.adapter.get_connectivity_edges(compress=compress)

            if filename == "edges.src.bin":
                self.send_binary(sources_data, head_only=head_only, compressed=compress)
            elif filename == "edges.dst.bin":
                self.send_binary(dests_data, head_only=head_only, compressed=compress)
            else:
                self.send_error_response(404, f"Unknown connectivity file: {filename} (expected edges.src.bin or edges.dst.bin)")
        except ValueError as e:
            self.send_error_response(404, str(e))
        except Exception as e:
            self.send_error_response(500, str(e))

    def log_message(self, format: str, *args):
        """Override to use Python logging."""
        logger.debug("%s - %s", self.address_string(), format % args)


class AnnDataServer:
    """
    Server for serving AnnData data in Cellucid format.

    Example:
        server = AnnDataServer(adata)
        server.start()  # Blocking

        # Or non-blocking:
        server.start_background()
        # ... do other things ...
        server.stop()
    """

    def __init__(
        self,
        data: Union[str, Path, "anndata.AnnData"],
        port: int = DEFAULT_PORT,
        host: str = DEFAULT_HOST,
        open_browser: bool = False,
        quiet: bool = False,
        **adapter_kwargs,
    ):
        """
        Initialize the server.

        Parameters
        ----------
        data : str, Path, or AnnData
            Path to h5ad file or AnnData object.
        port : int
            Port to serve on.
        host : str
            Host to bind to.
        open_browser : bool
            Whether to open browser on start.
        quiet : bool
            Suppress info messages.
        **adapter_kwargs
            Additional arguments passed to AnnDataAdapter.
        """
        self.port = port
        self.host = host
        self.open_browser = open_browser
        self.quiet = quiet

        # Create adapter
        if isinstance(data, (str, Path)):
            self.adapter = AnnDataAdapter.from_file(data, **adapter_kwargs)
            self.data_source = str(data)
        else:
            self.adapter = AnnDataAdapter(data, **adapter_kwargs)
            self.data_source = "in-memory AnnData"

        self._server: Optional[HTTPServer] = None
        self._thread: Optional[threading.Thread] = None
        self._running = False

        # Version
        try:
            from . import __version__
        except ImportError:
            __version__ = "0.0.0"

        self.server_info = {
            "version": __version__,
            "type": "anndata",
            "data_source": self.data_source,
            "host": self.host,
            "port": self.port,
            "n_cells": self.adapter.n_cells,
            "n_genes": self.adapter.n_genes,
            "is_backed": self.adapter.is_backed,
        }

    @property
    def url(self) -> str:
        """Get the server URL."""
        return f"http://{self.host}:{self.port}"

    @property
    def viewer_url(self) -> str:
        """Get the full URL to open the viewer."""
        return f"{CELLUCID_WEB_URL}?remote={self.url}&anndata=true"

    def start(self, blocking: bool = True):
        """Start the server."""
        if self._running:
            logger.warning("Server is already running")
            return

        # Ensure port is available (finds new one if needed)
        self.port = ensure_port_available(self.host, self.port, self.quiet)
        self.server_info["port"] = self.port

        # Create handler with adapter
        handler = partial(
            AnnDataRequestHandler,
            adapter=self.adapter,
            server_info=self.server_info,
        )

        self._server = HTTPServer((self.host, self.port), handler)
        self._running = True

        if not self.quiet:
            print(f"\n{'=' * 60}")
            print(f"  Cellucid AnnData Server")
            print(f"{'=' * 60}")
            print(f"  Data source:    {self.data_source}")
            print(f"  Cells:          {self.adapter.n_cells:,}")
            print(f"  Genes:          {self.adapter.n_genes:,}")
            print(f"  Dimensions:     {self.adapter._available_dimensions}")
            print(f"  Backed mode:    {self.adapter.is_backed}")
            print(f"  Server URL:     {self.url}")
            print(f"")
            print(f"  Open viewer at:")
            print(f"    {self.viewer_url}")
            print(f"")
            print(f"  NOTE: Loading data directly from AnnData is slower than")
            print(f"  using prepare. For production use, consider")
            print(f"  exporting your data first.")
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
        """Stop the server and cleanup resources."""
        if self._server:
            self._server.shutdown()
            self._server = None
        self._running = False

        # IMPORTANT: Close the adapter to release memory and file handles
        if hasattr(self, 'adapter') and self.adapter:
            self.adapter.close()

        if not self.quiet:
            logger.info("AnnData server stopped")

    def is_running(self) -> bool:
        """Check if the server is running."""
        return self._running

    def wait(self):
        """Wait for the server to stop."""
        if self._thread:
            try:
                while self._running:
                    self._thread.join(timeout=1)
            except KeyboardInterrupt:
                self.stop()


def serve_anndata(
    data: Union[str, Path, "anndata.AnnData"],
    port: int = DEFAULT_PORT,
    host: str = DEFAULT_HOST,
    open_browser: bool = True,
    quiet: bool = False,
    **kwargs,
) -> AnnDataServer:
    """
    Serve an AnnData object or h5ad file directly.

    This is a convenience function for quickly visualizing AnnData.
    For production use, consider using prepare instead.

    Parameters
    ----------
    data : str, Path, or AnnData
        Path to h5ad file or AnnData object.
    port : int
        Port to serve on (default: 8765).
    host : str
        Host to bind to (default: 127.0.0.1).
    open_browser : bool
        Whether to open browser (default: True).
    quiet : bool
        Suppress info messages.
    **kwargs
        Additional arguments passed to AnnDataAdapter.

    Returns
    -------
    AnnDataServer
        The running server instance.

    Example
    -------
    >>> from cellucid import serve_anndata
    >>> serve_anndata("/path/to/data.h5ad")

    >>> # Or with in-memory AnnData
    >>> import anndata as ad
    >>> adata = ad.read_h5ad("data.h5ad")
    >>> serve_anndata(adata)
    """
    server = AnnDataServer(
        data=data,
        port=port,
        host=host,
        open_browser=open_browser,
        quiet=quiet,
        **kwargs,
    )
    server.start()
    return server


def main():
    """CLI entry point for cellucid serve-anndata (legacy, kept for direct imports)."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Serve AnnData (h5ad or zarr) directly for visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Serve an h5ad file
    cellucid serve-anndata /path/to/data.h5ad

    # Serve a zarr store
    cellucid serve-anndata /path/to/data.zarr

    # Serve on a different port
    cellucid serve-anndata /path/to/data.h5ad --port 9000

Supported formats:
    - .h5ad files: HDF5-based AnnData files
    - .zarr directories: Zarr-based AnnData stores

Note:
    This serves data directly from AnnData, which is slower than using
    pre-exported data. For production use, export your data first,
    then serve the exported directory with:
        cellucid serve /path/to/export/
""",
    )

    parser.add_argument(
        "anndata_path",
        type=str,
        help="Path to h5ad file or zarr directory",
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
        "--no-backed",
        action="store_true",
        help="Load entire file into memory (default: lazy loading for h5ad)",
    )
    parser.add_argument(
        "--latent-key",
        type=str,
        default=None,
        help="Key in obsm for latent space (auto-detected if not specified)",
    )

    args = parser.parse_args()

    # Configure logging
    if not args.quiet:
        logging.basicConfig(level=logging.INFO)

    serve_anndata(
        data=args.anndata_path,
        port=args.port,
        host=args.host,
        open_browser=not args.no_browser,
        quiet=args.quiet,
        backed=not args.no_backed,
        latent_key=args.latent_key,
    )


if __name__ == "__main__":
    main()
