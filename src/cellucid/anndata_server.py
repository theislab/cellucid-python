"""
AnnData Server for Cellucid

HTTP server that serves AnnData data dynamically in the Cellucid format.
This allows direct visualization of AnnData files without pre-export.

Supports:
- h5ad files: Single HDF5 files (backed mode for lazy loading)
- zarr stores: Directory-based stores materialized by anndata.read_zarr
- In-memory AnnData objects

Usage:
    from cellucid.anndata_server import serve_anndata

    # Serve an h5ad file
    serve_anndata(
        "/path/to/data.h5ad",
        dataset_name="Example",
        dataset_id="example",
    )

    # Serve a zarr store (must be a directory)
    serve_anndata(
        "/path/to/data.zarr",
        dataset_name="Example",
        dataset_id="example",
    )

    # Or an in-memory AnnData
    serve_anndata(
        adata,
        dataset_name="Example",
        dataset_id="example",
    )
"""

from __future__ import annotations

import json
import logging
import re
import threading
import webbrowser
from functools import partial
from http import HTTPStatus
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.parse import unquote, urlparse

from ._server_base import (
    CELLUCID_WEB_URL,
    DEFAULT_HOST,
    DEFAULT_PORT,
    WEB_ASSET_INVENTORY_FILENAME,
    CORSMixin,
    _web_cache_dir,
    print_detail,
    print_server_banner,
    print_step,
    print_success,
    require_server_port,
)
from .anndata_adapter import AnnDataAdapter, _classify_anndata_path
from .connectivity_contract import build_connectivity_manifest

if TYPE_CHECKING:
    import anndata

logger = logging.getLogger("cellucid.anndata_server")


_QVALUE_RE = re.compile(r"(?:0(?:\.\d{0,3})?|1(?:\.0{0,3})?)\Z")


def _accepts_gzip(header_value: str) -> bool:
    """Return whether an Accept-Encoding value permits the exact gzip coding."""
    gzip_qualities: list[float] = []
    wildcard_qualities: list[float] = []

    for item in header_value.split(","):
        segments = [segment.strip() for segment in item.split(";")]
        coding = segments[0].lower()
        if not coding:
            continue

        quality = 1.0
        seen_quality = False
        valid = True
        for parameter in segments[1:]:
            if "=" not in parameter:
                valid = False
                break
            name, value = (part.strip() for part in parameter.split("=", 1))
            if name.lower() != "q" or seen_quality or not _QVALUE_RE.fullmatch(value):
                valid = False
                break
            quality = float(value)
            seen_quality = True

        if not valid:
            continue
        if coding == "gzip":
            gzip_qualities.append(quality)
        elif coding == "*":
            wildcard_qualities.append(quality)

    if gzip_qualities:
        return max(gzip_qualities) > 0
    return bool(wildcard_qualities and max(wildcard_qualities) > 0)


class AnnDataRequestHandler(CORSMixin, SimpleHTTPRequestHandler):
    """
    HTTP handler for serving AnnData data in Cellucid format.

    Routes:
        /dataset_identity.json - Dataset metadata
        /obs_manifest.json - Cell metadata manifest
        /var_manifest.json - Gene expression manifest
        /connectivity_manifest.json - Connectivity manifest
        /points_{dim}d.bin - Embedding coordinates
        /vectors/{fieldId}_{dim}d.bin - Vector field displacement vectors
        /obs/{field}.values.f32 - Continuous obs field
        /obs/{field}.codes.{ext} - Categorical obs field codes
        /obs/{field}.outliers.f32 - Categorical outlier quantiles
        /var/{gene}.values.f32 - Gene expression values
        /connectivity/edges.src.bin - Edge sources
        /connectivity/edges.dst.bin - Edge destinations
        /connectivity/edges.weights.f64.bin - Edge weights
        /_cellucid/health - Health check
        /_cellucid/info - Server info
        /_cellucid/events - POST endpoint for frontend events (hooks)
    """

    allow_caching = False  # Data is dynamic, don't cache

    def __init__(
        self,
        *args,
        adapter: AnnDataAdapter,
        server_info: dict,
        metadata_bodies: dict[str, bytes],
        payload_lengths: dict[str, int],
        serve_web_ui: bool,
        web_cache_dir: Path,
        **kwargs,
    ):
        self.adapter = adapter
        self.server_info = server_info
        self.metadata_bodies = metadata_bodies
        self.payload_lengths = payload_lengths
        self.serve_web_ui = serve_web_ui
        self.web_cache_dir = web_cache_dir
        # We serve "virtual" files via explicit route handlers, but still want
        # the normal BaseHTTPRequestHandler lifecycle (parse_request, etc.).
        # `directory` is unused because we never call the default file-serving
        # paths for AnnData mode.
        super().__init__(*args, **kwargs)

    def end_headers(self):
        """Add CORS headers to every response."""
        self.add_cors_headers()
        super().end_headers()

    def do_POST(self):
        """Handle POST requests (events from frontend)."""
        if self.handle_event_post():
            return
        if self.handle_session_bundle_post():
            return
        # No other POST endpoints - return 404
        self.send_error_response(404, f"POST not supported for path: {self.path}")

    def do_HEAD(self):
        """Handle HEAD requests (for urlExists checks)."""
        self.do_GET(head_only=True)

    def do_GET(self, head_only: bool = False):
        """Handle GET requests."""
        parsed = urlparse(self.path)
        raw_path = unquote(parsed.path)

        if self.serve_web_ui and self.handle_web_asset_get(raw_path, head_only=head_only):
            return
        if self.serve_web_ui and (
            raw_path == f"/{WEB_ASSET_INVENTORY_FILENAME}"
            or raw_path == "/assets"
            or raw_path.startswith("/assets/")
        ):
            self.send_error_response(404, "Web asset is not declared by the active build")
            return

        path = raw_path.lstrip("/")

        # HEAD always describes the uncompressed representation so it can be
        # answered from validated metadata without materializing payload bytes.
        accept_encoding = self.headers.get("Accept-Encoding", "")
        supports_gzip = not head_only and _accepts_gzip(accept_encoding)

        try:
            # Root path - redirect to viewer
            if path == "" or path == "index.html":
                self.send_error_response(503, "Cellucid viewer UI unavailable")
                return

            if path in self.metadata_bodies:
                self._send_preencoded_json(self.metadata_bodies[path], head_only)
            elif path == "_cellucid/health":
                self.send_json(
                    {
                        "status": "ok",
                        "type": "anndata",
                        "version": self.server_info["version"],
                        "format": self.server_info["format"],
                        "is_backed": self.server_info["is_backed"],
                        "n_cells": self.adapter.n_cells,
                        "n_genes": self.adapter.n_genes,
                    },
                    head_only,
                )
            elif path == "_cellucid/info":
                self.send_json(self.server_info, head_only)
            elif path == "connectivity_manifest.json":
                self.send_error_response(404, "No connectivity data")
            elif head_only and path in self.payload_lengths:
                self._send_binary_head(self.payload_lengths[path])
            elif path not in self.payload_lengths:
                self.send_error_response(404, f"Not found: {path}")
            elif re.fullmatch(r"points_[123]d\.bin", path):
                self._handle_points(path, head_only, supports_gzip)
            elif path.startswith("obs/"):
                self._handle_obs(path, head_only, supports_gzip)
            elif path.startswith("var/"):
                self._handle_var(path, head_only, supports_gzip)
            elif path.startswith("connectivity/"):
                self._handle_connectivity(path, head_only, supports_gzip)
            elif path.startswith("vectors/"):
                self._handle_vector_fields(path, head_only, supports_gzip)
            else:
                self.send_error_response(404, f"Not found: {path}")

        except Exception:
            logger.exception("Error handling %s", path)
            self.send_error_response(500, "Internal server error")

    def _send_preencoded_json(self, body: bytes, head_only: bool) -> None:
        """Send validated JSON bytes without rebuilding them for HEAD."""
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if not head_only:
            self.wfile.write(body)

    def _send_binary_head(self, content_length: int) -> None:
        """Send headers for an uncompressed scientific payload."""
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(content_length))
        self.end_headers()

    def _handle_points(self, path: str, head_only: bool, supports_gzip: bool):
        """Handle points_Xd.bin requests."""
        match = re.fullmatch(r"points_([123])d\.bin", path)
        if not match:
            self.send_error_response(404, f"Invalid points path: {path}")
            return

        dim = int(match.group(1))
        try:
            # Compress if client supports gzip (transparent compression for better network perf)
            compress = supports_gzip
            data = self.adapter.get_points_binary(dim, compress=compress)
            self.send_binary(
                data,
                head_only=head_only,
                compressed=compress,
                vary_accept_encoding=True,
            )
        except ValueError:
            self.send_error_response(404, "Points data not available")

    def _handle_vector_fields(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle vector field requests under /vectors/.

        Supported paths:
        - vectors/{fieldId}_{dim}d.bin
        """
        match = re.fullmatch(r"vectors/(.+)_([123])d\.bin", path)
        if not match:
            self.send_error_response(404, f"Invalid vectors path: {path}")
            return

        field_id = match.group(1)
        dim = int(match.group(2))

        try:
            data = self.adapter.get_vector_field_binary(
                field_id,
                dim,
                compress=supports_gzip,
            )
            self.send_binary(
                data,
                head_only=head_only,
                compressed=supports_gzip,
                vary_accept_encoding=True,
            )
        except ValueError as e:
            self.send_error_response(404, str(e))

    def _handle_obs(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle obs field requests.

        Supported paths:
        - obs/{field}.values.f32
        - obs/{field}.codes.u8 or obs/{field}.codes.u16
        - obs/{field}.outliers.f32
        """
        # Remove 'obs/' prefix
        filename = path[4:]

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

        actual_key = self.adapter.get_obs_key_for_payload_component(field_name)

        if actual_key is None:
            self.send_error_response(404, f"Obs field not found: {field_name}")
            return

        try:
            if data_type == "values":
                # Continuous field values
                data = self.adapter.get_obs_continuous_values(actual_key, compress=compress)
                self.send_binary(
                    data,
                    head_only=head_only,
                    compressed=compress,
                    vary_accept_encoding=True,
                )
            elif data_type == "codes":
                # Categorical codes
                data, categories, missing = self.adapter.get_obs_categorical_codes(
                    actual_key, compress=compress
                )
                self.send_binary(
                    data,
                    head_only=head_only,
                    compressed=compress,
                    vary_accept_encoding=True,
                )
            elif data_type == "outliers":
                # Outlier quantiles for categorical field
                data = self.adapter.get_obs_outlier_quantiles(actual_key, compress=compress)
                self.send_binary(
                    data,
                    head_only=head_only,
                    compressed=compress,
                    vary_accept_encoding=True,
                )
            else:
                self.send_error_response(
                    404, f"Unknown obs data type: {data_type} (expected values, codes, or outliers)"
                )
        except Exception:
            logger.exception("Error handling obs request: %s", path)
            self.send_error_response(500, "Internal server error")

    def _handle_var(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle var (gene expression) requests.

        Supported paths:
        - var/{gene}.values.f32
        """
        # Remove 'var/' prefix
        filename = path[4:]

        # Compress if client supports gzip (transparent compression)
        compress = supports_gzip

        # Check format: {gene}.values.f32
        if not filename.endswith(".values.f32"):
            self.send_error_response(
                404, f"Invalid var path format: {path} (expected {{gene}}.values.f32)"
            )
            return

        gene_safe = filename[:-11]  # Remove '.values.f32'

        actual_gene = self.adapter.get_gene_id_for_payload_component(gene_safe)

        if actual_gene is None:
            self.send_error_response(404, f"Gene not found: {gene_safe}")
            return

        try:
            data = self.adapter.get_gene_expression(actual_gene, compress=compress)
            self.send_binary(
                data,
                head_only=head_only,
                compressed=compress,
                vary_accept_encoding=True,
            )
        except KeyError:
            self.send_error_response(404, "Gene not found")
        except Exception:
            logger.exception("Error handling var request: %s", path)
            self.send_error_response(500, "Internal server error")

    def _handle_connectivity(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle connectivity edge requests.

        Supported paths:
        - connectivity/edges.src.bin - Edge source indices
        - connectivity/edges.dst.bin - Edge destination indices
        - connectivity/edges.weights.f64.bin - Edge Float64 weights
        """
        filename = path.split("/")[-1]

        # Compress if client supports gzip (transparent compression)
        compress = supports_gzip

        try:
            (
                sources_data,
                dests_data,
                weights_data,
                _n_edges,
                _max_neighbors,
            ) = self.adapter.get_connectivity_edges(compress=compress)

            if filename == "edges.src.bin":
                self.send_binary(
                    sources_data,
                    head_only=head_only,
                    compressed=compress,
                    vary_accept_encoding=True,
                )
            elif filename == "edges.dst.bin":
                self.send_binary(
                    dests_data,
                    head_only=head_only,
                    compressed=compress,
                    vary_accept_encoding=True,
                )
            elif filename == "edges.weights.f64.bin":
                self.send_binary(
                    weights_data,
                    head_only=head_only,
                    compressed=compress,
                    vary_accept_encoding=True,
                )
            else:
                self.send_error_response(
                    404,
                    "Unknown connectivity file: "
                    f"{filename} (expected edges.src.bin, edges.dst.bin, "
                    "or edges.weights.f64.bin)",
                )
        except ValueError:
            self.send_error_response(404, "Connectivity data not available")
        except Exception:
            logger.exception("Error handling connectivity request: %s", path)
            self.send_error_response(500, "Internal server error")

    def log_message(self, format: str, *args):
        """Override to use Python logging."""
        logger.debug("%s - %s", self.address_string(), format % args)


class AnnDataServer:
    """
    Server for serving AnnData data in Cellucid format.

    Examples
    --------
    Start a blocking server::

        AnnDataServer(
            adata,
            dataset_name="Example",
            dataset_id="example",
        ).start()

    Start and stop a background server::

        server = AnnDataServer(
            adata,
            dataset_name="Example",
            dataset_id="example",
        )
        server.start_background()
        server.stop()
    """

    def __init__(
        self,
        data: str | Path | anndata.AnnData,
        port: int = DEFAULT_PORT,
        host: str = DEFAULT_HOST,
        open_browser: bool = False,
        quiet: bool = False,
        *,
        latent_key: str | None = None,
        gene_id_column: str | None = None,
        normalize_embeddings: bool = True,
        centroid_outlier_quantile: float = 0.95,
        centroid_min_points: int = 10,
        dataset_name: str,
        dataset_id: str,
        vector_field_default: str | None = None,
        serve_web_ui: bool = True,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
    ) -> None:
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
        latent_key : str, optional
            Explicit key in ``obsm`` for the latent space.
        gene_id_column : str, optional
            Exact column in ``var`` containing gene identifiers. If None,
            identifiers come from ``var.index``.
        normalize_embeddings : bool
            Whether to normalize embeddings into the viewer coordinate range.
        centroid_outlier_quantile : float
            Quantile used when computing categorical centroids.
        centroid_min_points : int
            Minimum category size used for centroid computation.
        dataset_name : str
            Explicit human-readable dataset name.
        dataset_id : str
            Explicit stable dataset identifier.
        vector_field_default : str, optional
            Exact field id required when multiple UMAP vector fields exist.
        serve_web_ui : bool
            Establish and serve the exact current web build.
        web_source_url : str
            Origin publishing the web asset inventory.
        web_cache_dir : str or Path, optional
            Directory holding the active verified web build.
        """
        self.port = require_server_port(port)
        self.host = host
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
        self._server: HTTPServer | None = None
        self._thread: threading.Thread | None = None
        self._running = False
        self._started = False
        self._closed = False
        self._serving = False
        self._background_error: BaseException | None = None

        # Step 1: Detect format
        format_name: str
        self.data_format: str
        if isinstance(data, str | Path):
            format_name = _classify_anndata_path(data)
            self.data_source = str(data)
            self.data_format = format_name

            if not quiet:
                print_step(1, 4, "Detecting format")
                print_detail("Path", str(data))
                print_detail("Format", format_name)
                print_success("Format detected")
        else:
            # In-memory AnnData
            self.data_source = "in-memory AnnData"
            self.data_format = "in-memory"
            if not quiet:
                print_step(1, 4, "Detecting format")
                print_detail("Source", "in-memory AnnData")
                print_success("Format detected")

        # Step 2: Load file
        if not quiet:
            print_step(2, 4, "Loading AnnData")
            mode = "read-only backed h5ad" if self.data_format == "h5ad" else "in-memory"
            print_detail("Mode", mode)

        if isinstance(data, str | Path):
            self.adapter = AnnDataAdapter.from_file(
                data,
                latent_key=latent_key,
                gene_id_column=gene_id_column,
                normalize_embeddings=normalize_embeddings,
                centroid_outlier_quantile=centroid_outlier_quantile,
                centroid_min_points=centroid_min_points,
                dataset_name=dataset_name,
                dataset_id=dataset_id,
                vector_field_default=vector_field_default,
            )
        else:
            self.adapter = AnnDataAdapter(
                data,
                latent_key=latent_key,
                gene_id_column=gene_id_column,
                normalize_embeddings=normalize_embeddings,
                centroid_outlier_quantile=centroid_outlier_quantile,
                centroid_min_points=centroid_min_points,
                dataset_name=dataset_name,
                dataset_id=dataset_id,
                vector_field_default=vector_field_default,
            )

        try:
            if not quiet:
                print_success("File opened")

            # Step 3: Analyze dataset
            if not quiet:
                print_step(3, 4, "Analyzing dataset")
                self._print_dataset_info()
                print_success("Analysis complete")

            from . import __version__

            self.server_info = {
                "version": __version__,
                "type": "anndata",
                "format": self.data_format,
                "host": self.host,
                "port": self.port,
                "n_cells": self.adapter.n_cells,
                "n_genes": self.adapter.n_genes,
                "is_backed": self.adapter.is_backed,
            }
            (
                self._identity,
                self._metadata_bodies,
                self._payload_lengths,
            ) = self._build_http_contract()
        except BaseException:
            self._closed = True
            try:
                self.adapter.close()
            except BaseException:
                logger.exception(
                    "Failed to close the AnnData adapter after server construction failed"
                )
            raise

    def _build_http_contract(self) -> tuple[dict, dict[str, bytes], dict[str, int]]:
        """Build and cross-check every public direct-AnnData route once."""
        identity = self.adapter.get_dataset_identity()
        obs_manifest = self.adapter.get_obs_manifest()
        var_manifest = self.adapter.get_var_manifest()
        connectivity_manifest = self.adapter.get_connectivity_manifest()

        has_connectivity = identity["stats"]["has_connectivity"]
        n_edges = identity["stats"]["n_edges"]
        if has_connectivity:
            if connectivity_manifest is None:
                raise ValueError("Identity advertises connectivity without a connectivity manifest")
            expected_connectivity_keys = {
                "format",
                "n_cells",
                "n_edges",
                "max_neighbors",
                "index_bytes",
                "index_dtype",
                "sourcesPath",
                "destinationsPath",
                "weightsPath",
                "weight_dtype",
                "weight_bytes",
                "compression",
            }
            if set(connectivity_manifest) != expected_connectivity_keys:
                raise ValueError("Connectivity manifest does not contain the exact current fields")
            if n_edges != connectivity_manifest["n_edges"]:
                raise ValueError("Identity and connectivity manifest disagree on the edge count")
            expected_connectivity_manifest = build_connectivity_manifest(
                n_cells=identity["stats"]["n_cells"],
                n_edges=connectivity_manifest["n_edges"],
                max_neighbors=connectivity_manifest["max_neighbors"],
                index_bytes=connectivity_manifest["index_bytes"],
                index_dtype=connectivity_manifest["index_dtype"],
                compression=None,
            )
            if (
                list(connectivity_manifest) != list(expected_connectivity_manifest)
                or connectivity_manifest != expected_connectivity_manifest
            ):
                raise ValueError(
                    "Connectivity manifest does not match the direct weighted edge contract"
                )
        elif connectivity_manifest is not None or n_edges is not None:
            raise ValueError(
                "Connectivity manifest or edge count exists while connectivity is absent"
            )

        metadata = {
            "dataset_identity.json": identity,
            "obs_manifest.json": obs_manifest,
            "var_manifest.json": var_manifest,
            "_cellucid/datasets": {
                "datasets": [
                    {
                        "id": identity["id"],
                        "path": "/",
                        "name": identity["name"],
                    }
                ]
            },
        }
        if connectivity_manifest is not None:
            metadata["connectivity_manifest.json"] = connectivity_manifest
        metadata_bodies = {
            path: json.dumps(value).encode("utf-8") for path, value in metadata.items()
        }

        n_cells = identity["stats"]["n_cells"]
        payload_lengths: dict[str, int] = {}

        def add_payload(path: str, content_length: int) -> None:
            if not isinstance(path, str) or not path or path.startswith("/"):
                raise ValueError(f"Invalid direct-AnnData payload path: {path!r}")
            if path in payload_lengths:
                raise ValueError(f"Duplicate direct-AnnData payload path: {path!r}")
            if type(content_length) is not int or content_length < 0:
                raise ValueError(f"Invalid payload length for {path!r}: {content_length!r}")
            payload_lengths[path] = content_length

        embeddings = identity["embeddings"]
        for dim in embeddings["available_dimensions"]:
            add_payload(embeddings["files"][f"{dim}d"], n_cells * dim * 4)

        vector_fields = identity.get("vector_fields")
        if vector_fields is not None:
            for field in vector_fields["fields"].values():
                for dim in field["available_dimensions"]:
                    add_payload(field["files"][f"{dim}d"], n_cells * dim * 4)

        obs_schemas = obs_manifest["_obsSchemas"]
        continuous_schema = obs_schemas.get("continuous")
        if continuous_schema is not None:
            for field in obs_manifest["_continuousFields"]:
                key = self.adapter.get_obs_payload_component(field[0])
                add_payload(
                    continuous_schema["pathPattern"].format(key=key),
                    n_cells * 4,
                )

        categorical_schema = obs_schemas.get("categorical")
        if categorical_schema is not None:
            for field in obs_manifest["_categoricalFields"]:
                key = self.adapter.get_obs_payload_component(field[0])
                ext = {"uint8": "u8", "uint16": "u16"}[field[2]]
                item_size = {"uint8": 1, "uint16": 2}[field[2]]
                add_payload(
                    categorical_schema["codesPathPattern"].format(
                        key=key,
                        ext=ext,
                    ),
                    n_cells * item_size,
                )
                outlier_pattern = categorical_schema["outlierPathPattern"]
                if outlier_pattern is not None:
                    add_payload(outlier_pattern.format(key=key), n_cells * 4)

        var_schema = var_manifest["_varSchema"]
        for field in var_manifest["fields"]:
            key = self.adapter.get_gene_payload_component(field[0])
            add_payload(var_schema["pathPattern"].format(key=key), n_cells * 4)

        if connectivity_manifest is not None:
            edge_bytes = connectivity_manifest["n_edges"] * connectivity_manifest["index_bytes"]
            add_payload(connectivity_manifest["sourcesPath"], edge_bytes)
            add_payload(connectivity_manifest["destinationsPath"], edge_bytes)
            weight_bytes = connectivity_manifest["n_edges"] * connectivity_manifest["weight_bytes"]
            add_payload(connectivity_manifest["weightsPath"], weight_bytes)

        return identity, metadata_bodies, payload_lengths

    def _print_dataset_info(self):
        """Print information about the loaded dataset."""
        print_detail("Cells", f"{self.adapter.n_cells:,}")
        print_detail("Genes", f"{self.adapter.n_genes:,}")

        # Get embedding info
        dims = self.adapter._available_dimensions
        dim_str = ", ".join(f"{d}D" for d in sorted(dims)) if dims else "none"
        print_detail("Embeddings", dim_str)

        # Count obs field types
        obs_keys = self.adapter.get_obs_keys()
        n_categorical = sum(1 for k in obs_keys if self.adapter.get_obs_field_kind(k) == "category")
        n_continuous = sum(
            1 for k in obs_keys if self.adapter.get_obs_field_kind(k) == "continuous"
        )
        print_detail("Obs fields", f"{n_categorical} categorical, {n_continuous} continuous")

        # Check connectivity
        has_conn = self.adapter.has_connectivity()
        if has_conn:
            manifest = self.adapter.get_connectivity_manifest()
            if manifest is None:
                raise RuntimeError("Validated connectivity is missing its manifest")
            conn_str = f"yes ({manifest['n_edges']:,} edges)"
        else:
            conn_str = "no"
        print_detail("Connectivity", conn_str)

    @property
    def url(self) -> str:
        """Get the URL of the currently running server."""
        if not self._running or self._server is None:
            raise RuntimeError("AnnDataServer URL is unavailable because the server is not running")
        return f"http://{self.host}:{self.port}"

    @property
    def viewer_url(self) -> str:
        """Get the full URL to open the viewer."""
        return f"{self.url}/?anndata=true"

    def start(self, blocking: bool = True):
        """Start this single-use server."""
        if type(blocking) is not bool:
            raise TypeError("blocking must be a boolean")
        if self._running:
            raise RuntimeError("Server is already running.")
        if self._started or self._closed:
            raise RuntimeError(
                "AnnDataServer is single-use and has been closed. Create a new server instance."
            )
        self._started = True

        try:
            # Step 4: Start server
            if not self.quiet:
                print_step(4, 4, "Starting server")

            if self.serve_web_ui:
                from .web_cache import ensure_web_ui_cached

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
                AnnDataRequestHandler,
                adapter=self.adapter,
                server_info=self.server_info,
                metadata_bodies=self._metadata_bodies,
                payload_lengths=self._payload_lengths,
                serve_web_ui=self.serve_web_ui,
                web_cache_dir=self.web_cache_dir,
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
                raise RuntimeError("AnnDataServer lost its bound HTTP server before serving")
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
                        "AnnData server cleanup also failed after its serving loop failed"
                    )
        self._background_error = serving_error
        if serving_error is not None:
            logger.error(
                "AnnData background server failed",
                exc_info=(
                    type(serving_error),
                    serving_error,
                    serving_error.__traceback__,
                ),
            )

    def _finish_serving(self) -> None:
        """Close the socket and adapter after the serving loop has ended."""
        failures: list[BaseException] = []
        self._running = False
        server = self._server
        self._server = None
        if server is not None:
            try:
                server.server_close()
            except BaseException as error:
                failures.append(error)
        try:
            self.adapter.close()
        except BaseException as error:
            failures.append(error)
        self._closed = True
        if failures:
            details = "; ".join(f"{type(error).__name__}: {error}" for error in failures)
            raise RuntimeError(f"AnnData server cleanup failed: {details}") from failures[0]

    def _rollback_failed_start(self, *, shutdown: bool) -> None:
        """Rollback acquired resources without replacing the startup exception."""
        self._running = False
        server = self._server
        if server is not None and shutdown and self._serving:
            try:
                server.shutdown()
            except BaseException:
                logger.exception("Failed to shut down the AnnData server after startup failed")
        if server is not None:
            try:
                server.server_close()
            except BaseException:
                logger.exception("Failed to close the AnnData socket after startup failed")
        thread = self._thread
        if thread is not None and thread is not threading.current_thread() and thread.is_alive():
            thread.join()
        self._server = None
        self._thread = None
        self._serving = False
        try:
            self.adapter.close()
        except BaseException:
            logger.exception("Failed to close the AnnData adapter after startup failed")
        self._closed = True

    def start_background(self):
        """Start the server in a background thread."""
        self.start(blocking=False)

    def stop(self):
        """Stop the server and cleanup resources."""
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

        try:
            self.adapter.close()
        except BaseException as error:
            failures.append(error)
        self._closed = True

        if not self.quiet:
            print("AnnData server stopped")
        if failures:
            details = "; ".join(f"{type(error).__name__}: {error}" for error in failures)
            raise RuntimeError(f"AnnData server shutdown failed: {details}") from failures[0]

    def is_running(self) -> bool:
        """Check if the server is running."""
        return self._running

    def wait(self):
        """Wait for the server to stop."""
        thread = self._thread
        if thread is not None:
            try:
                thread.join()
            except BaseException:
                self.stop()
                raise
        if self._background_error is not None:
            raise self._background_error


def serve_anndata(
    data: str | Path | anndata.AnnData,
    port: int = DEFAULT_PORT,
    host: str = DEFAULT_HOST,
    open_browser: bool = True,
    quiet: bool = False,
    *,
    latent_key: str | None = None,
    gene_id_column: str | None = None,
    normalize_embeddings: bool = True,
    centroid_outlier_quantile: float = 0.95,
    centroid_min_points: int = 10,
    dataset_name: str,
    dataset_id: str,
    vector_field_default: str | None = None,
    serve_web_ui: bool = True,
    web_source_url: str = CELLUCID_WEB_URL,
    web_cache_dir: str | Path | None = None,
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
    latent_key : str, optional
        Explicit key in ``obsm`` for the latent space.
    gene_id_column : str, optional
        Exact column in ``var`` containing gene identifiers. If None,
        identifiers come from ``var.index``.
    normalize_embeddings : bool
        Whether to normalize embeddings into the viewer coordinate range.
    centroid_outlier_quantile : float
        Quantile used when computing categorical centroids.
    centroid_min_points : int
        Minimum category size used for centroid computation.
    dataset_name : str
        Explicit human-readable dataset name.
    dataset_id : str
        Explicit stable dataset identifier.
    vector_field_default : str, optional
        Exact field id required when multiple UMAP vector fields exist.
    serve_web_ui : bool
        Establish and serve the exact current web build.
    web_source_url : str
        Origin publishing the web asset inventory.
    web_cache_dir : str or Path, optional
        Directory holding the active verified web build.
    Returns
    -------
    AnnDataServer
        The running server instance.

    Example
    -------
    >>> from cellucid import serve_anndata
    >>> serve_anndata(
    ...     "/path/to/data.h5ad",
    ...     dataset_name="Example",
    ...     dataset_id="example",
    ... )

    >>> # Or with in-memory AnnData
    >>> import anndata as ad
    >>> adata = ad.read_h5ad("data.h5ad")
    >>> serve_anndata(
    ...     adata,
    ...     dataset_name="Example",
    ...     dataset_id="example",
    ... )
    """
    server = AnnDataServer(
        data=data,
        port=port,
        host=host,
        open_browser=open_browser,
        quiet=quiet,
        latent_key=latent_key,
        gene_id_column=gene_id_column,
        normalize_embeddings=normalize_embeddings,
        centroid_outlier_quantile=centroid_outlier_quantile,
        centroid_min_points=centroid_min_points,
        dataset_name=dataset_name,
        dataset_id=dataset_id,
        vector_field_default=vector_field_default,
        serve_web_ui=serve_web_ui,
        web_source_url=web_source_url,
        web_cache_dir=web_cache_dir,
    )
    server.start()
    return server
