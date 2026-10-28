"""
AnnData Server for Cellucid

Serves one dataset to the Cellucid viewer, which reads it itself.

The web app reads every format of the input format registry in the browser
(``input_formats.json``): AnnData (``.h5ad``, Zarr), MuData (``.h5mu``, MuData
Zarr), Loom, 10x Cell Ranger ``.h5`` files and ``outs/`` folders, MatrixMarket
folders, and R objects (``.rds``, ``.rda``/``.RData``, gzipped or not). The
formats that allow it are read with HTTP range requests -- only the bytes a
view needs are transferred -- and the single-stream ones (R objects,
MatrixMarket text) are read once from the start. This server therefore
computes nothing per request. It serves:

- the data's bytes under ``/_cellucid/anndata/<name>``: a file with exact
  ``Range``/``206``/``Content-Range``/``HEAD``/``ETag`` answers, sent by the
  kernel (``sendfile``); a Zarr store folder key by key (404 for an absent
  key, consolidated metadata built in memory when the store has none); a
  folder of files (MatrixMarket, ``outs/``) file by file, listed in
  ``/_cellucid/datasets`` because a browser cannot list a folder;
- ``/_cellucid/datasets``, which names the dataset (``dataset_id``,
  ``dataset_name``), where its store is, and the launch choices (embeddings,
  graph, expression, gene names, latent space, vector field, obs columns,
  normalization, centroid rule) the viewer applies when it opens the store;
- the viewer channel (``/_cellucid/ws``, a WebSocket: commands from a
  notebook to its viewer, events back; :mod:`cellucid._viewer_channel`) and
  the viewer's web build.

An AnnData object, and an AnnData file (``.h5ad``, Zarr) whose cells are not
in priority order, are served as a Zarr store in **priority order**
(:mod:`cellucid._served_store`): the layout ``prepare()`` writes, made chunk
by chunk from the source as the browser reads it. So are a MuData file
(``.h5mu``, MuData Zarr: a MuData Zarr store, each modality's rows following
its cells), a Loom file and a 10x Cell Ranger matrix or ``outs/`` folder,
presented as the browser presents the file itself
(:mod:`cellucid._served_formats`). A subsample the viewer opens is then the
first rows of every array, so the bytes sent are in proportion to the cells
opened, whatever the source's layout. A file in priority order already
(``prepare()``'s AnnData or MuData) is served as it is, unless it holds an
array stored with a codec the browser does not decode (an HDF5 filter such as
Blosc or Zstandard, a Zarr codec such as bz2): then it is served through the
same store, which decodes it on the kernel. R objects and MatrixMarket
folders, which the browser reads once from the start, are served as they are.
A MuData object is written once to a temporary Cellucid-ready ``.h5mu`` (by
``prepare()``, in priority order), removed when the server stops or Python
exits.

Usage:
    from cellucid.anndata_server import serve_anndata

    serve_anndata("/path/to/data.h5ad", dataset_name="Example", dataset_id="example")
    serve_anndata("/path/to/data.zarr", dataset_name="Example", dataset_id="example")
    serve_anndata(adata, dataset_name="Example", dataset_id="example")
"""

from __future__ import annotations

import hashlib
import json
import logging
import threading
import weakref
import webbrowser
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace
from functools import partial
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING, Any
from urllib.parse import unquote, urlparse

from ._anndata_files import (
    ANNDATA_KINDS,
    MUDATA_KINDS,
    AnnDataContents,
    AnnDataStore,
    browser_undecodable,
    classify_anndata_path,
    consolidated_zarr_metadata,
    is_mudata,
    read_anndata_contents,
    read_cellucid_identity,
    read_into_memory,
    remove_temporary_file,
    write_temporary_h5mu,
)
from ._console import console_print
from ._contracts import (
    _require_dataset_id,
    _require_dataset_name,
)
from ._file_versions import file_version, path_version
from ._launch_defaults import launch_defaults
from ._plural import plural
from ._server_base import (
    CELLUCID_WEB_URL,
    DEFAULT_HOST,
    DEFAULT_PORT,
    WEB_ASSET_INVENTORY_FILENAME,
    CellucidThreadingHTTPServer,
    CORSMixin,
    _web_cache_dir,
    bind_host_is_wildcard,
    format_http_origin,
    loopback_origin_for_bind,
    machine_origin,
    parse_byte_range,
    print_detail,
    print_server_banner,
    print_step,
    print_success,
    require_allowed_hosts,
    require_server_port,
)
from ._viewer_channel import CHANNEL_ROUTE, accept_viewer_channel
from .data_mapping import EmbeddingKeySpec

if TYPE_CHECKING:
    import anndata
    import mudata

    from ._served_store import PriorityStore

logger = logging.getLogger("cellucid.anndata_server")

#: The one query a direct-AnnData viewer URL carries, in every place such a
#: URL is built, so the loopback form and the network form cannot drift.
VIEWER_QUERY = "/?anndata=true"

#: The route prefix every byte of the served store is read under.
ANNDATA_ROUTE_PREFIX = "_cellucid/anndata/"

#: How long a kept-alive connection may sit idle before its thread closes it.
_IDLE_CONNECTION_SECONDS = 120


def _protocol_capability_document() -> dict[str, Any]:
    """Return the viewer wire capabilities this installation declares.

    The answer is read from :mod:`cellucid.jupyter._wire`, which is the module
    the validator itself reads, so the route cannot publish a capability the
    notebook would then reject. The import is deferred because reaching ``_wire``
    initializes the whole ``cellucid.jupyter`` package; this package keeps
    such imports lazy so the CLI starts fast, and one viewer asks this route
    once.
    """
    from .jupyter import _wire

    return _wire._protocol_capability_document()


# =============================================================================
# HTTP
# =============================================================================


def _route_name(store: AnnDataStore, dataset_id: str) -> str:
    """The last route segment: the dataset id with the ending the format is named by.

    Named by the dataset id, so the local file name is never published; the
    ending is the file's own (``.rds.gz``, ``.loom``), a Zarr store folder's
    ``.zarr``, and none for a folder of files.
    """
    if store.layout == "folder":
        return f"{dataset_id}.zarr" if store.files is None else dataset_id
    return f"{dataset_id}{store.ending}"


@dataclass(frozen=True)
class _ServedStore:
    """What the handler serves: a store on disk as it is, or a store in
    priority order (:mod:`cellucid._served_store`), and its route name."""

    #: The last route segment (:func:`_route_name`; ``<id>.zarr`` for a store in priority order).
    route_name: str
    #: The format id the viewer reads the store as (``input_formats.json``).
    kind: str
    #: ``"file"``, ``"zip"`` or ``"folder"`` (a store in priority order is a Zarr folder).
    layout: str
    #: The store on disk served as it is; None for a store in priority order.
    store: AnnDataStore | None = None
    #: The store in priority order; None for a store on disk served as it is.
    priority: PriorityStore | None = None
    #: (store key, bytes) answered in place of a missing consolidated document.
    generated_metadata: tuple[str, bytes] | None = None
    #: The generated document's validator: a hash of its bytes, computed once.
    generated_etag: str = ""


#: The formats a server serves by path in priority order besides AnnData,
#: read on the kernel (:mod:`cellucid._served_store`, :mod:`cellucid._served_formats`):
#: MuData, Loom, 10x Cell Ranger files and ``outs/`` folders.
PRIORITY_SOURCE_KINDS = MUDATA_KINDS | frozenset({"loom", "tenx-h5", "tenx-outs"})


def _unreadable(error: OSError) -> str:
    """What a dataset list says of a file served as it is that is gone."""
    name = Path(error.filename).name if error.filename else "The served file"
    return f"{name} is no longer readable"


def _canonical_store_key(text: str) -> str | None:
    """One Zarr store key, or None when it is not a plain relative key."""
    if not text or text.startswith("/") or "\\" in text or "\x00" in text:
        return None
    parts = text.split("/")
    if any(part in ("", ".", "..") for part in parts):
        return None
    return text


def _resolve_store_route(served: _ServedStore, route: str) -> Path | tuple[bytes, str] | None:
    """What one route under :data:`ANNDATA_ROUTE_PREFIX` names: the file to
    answer from, bytes and their validator (generated metadata, a key of a
    store in priority order), or None for a route the store does not have.

    A Zarr store folder is answered key by key; an absent key is None, which
    a Zarr reader reads as "absent" (an unwritten chunk is its fill value). A
    key may not leave the store.
    """
    if served.layout != "folder":
        assert served.store is not None
        return served.store.path if route == served.route_name else None
    prefix = f"{served.route_name}/"
    key = _canonical_store_key(route[len(prefix) :]) if route.startswith(prefix) else None
    if key is None:
        return None
    if served.priority is not None:
        body = served.priority.get(key)
        return None if body is None else (body, served.priority.etag)
    generated = served.generated_metadata
    if generated is not None and key == generated[0]:
        return generated[1], served.generated_etag
    assert served.store is not None
    path = served.store.path.joinpath(*PurePosixPath(key).parts)
    try:
        resolved = path.resolve(strict=True)
    except (OSError, RuntimeError):
        return None
    if not resolved.is_file() or not resolved.is_relative_to(served.store.path):
        return None
    return resolved


class AnnDataRequestHandler(CORSMixin, BaseHTTPRequestHandler):
    """
    HTTP handler of an :class:`AnnDataServer`.

    Routes:
        /_cellucid/anndata/<name>[/<key>] - The store's bytes (GET, HEAD, Range)
        /_cellucid/datasets - The dataset, its store, and its launch choices
        /_cellucid/health - Health check
        /_cellucid/info - Server info
        /_cellucid/protocol - The viewer channel, and the commands and events it carries
        /_cellucid/ws - The viewer channel (a WebSocket)
    """

    # Kept-alive connections: a browser reads a file with many range requests.
    protocol_version = "HTTP/1.1"
    timeout = _IDLE_CONNECTION_SECONDS
    # Store bytes are kept by the browser and revalidated (ETag, 304) before
    # reuse: a reload of an unchanged file transfers nothing it has.
    allow_caching = False
    # Headers and body are separate writes; without TCP_NODELAY the body's
    # last segment would wait for the client's delayed ACK of the headers.
    disable_nagle_algorithm = True

    def __init__(
        self,
        *args,
        served: _ServedStore,
        server_info: dict,
        health: dict,
        datasets: Callable[[], dict[str, Any]],
        serve_web_ui: bool,
        web_cache_dir: Path,
        allowed_hosts: frozenset[str],
        **kwargs,
    ):
        self.served = served
        self.server_info = server_info
        self.health = health
        #: The dataset list, made at each request: a file served as it is is
        #: described as it is then (its size and validator).
        self.datasets = datasets
        self.serve_web_ui = serve_web_ui
        self.web_cache_dir = web_cache_dir
        # Read by CORSMixin.parse_request, which runs inside super().__init__().
        self.allowed_hosts = allowed_hosts
        super().__init__(*args, **kwargs)

    def end_headers(self):
        """Add CORS headers to every response."""
        self.add_cors_headers()
        super().end_headers()

    def do_HEAD(self):
        """Answer GET's headers without its body."""
        self._handle_read(head_only=True)

    def do_GET(self):
        """Answer one route."""
        self._handle_read(head_only=False)

    def _handle_read(self, *, head_only: bool) -> None:
        raw_path = unquote(urlparse(self.path).path)
        try:
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
            if path in ("", "index.html"):
                self.send_error_response(503, "Cellucid viewer UI unavailable")
            elif path == "_cellucid/health":
                self.send_json(self.health, head_only)
            elif path == "_cellucid/info":
                self.send_json(self.server_info, head_only)
            elif path == "_cellucid/datasets":
                try:
                    datasets = self.datasets()
                except OSError as error:
                    # The file served as it is is gone (or unreadable).
                    self.send_error_response(404, _unreadable(error))
                    return
                self.send_json(datasets, head_only)
            elif path == "_cellucid/protocol":
                self.send_json(_protocol_capability_document(), head_only)
            elif path == CHANNEL_ROUTE:
                accept_viewer_channel(self)
            elif path.startswith(ANNDATA_ROUTE_PREFIX):
                self._handle_store(path[len(ANNDATA_ROUTE_PREFIX) :], head_only=head_only)
            else:
                self.send_error_response(404, f"Not found: {path}")
        except ConnectionError:
            # The client went away mid-response: a closed tab, a reload, or a
            # cancelled read. There is no socket left to answer on.
            self.close_connection = True
            logger.debug("Client disconnected while receiving %s", raw_path)
        except Exception:
            logger.exception("Error handling %s", raw_path)
            self.close_connection = True
            self.send_error_response(500, "Internal server error")

    def _handle_store(self, route: str, *, head_only: bool) -> None:
        resolved = _resolve_store_route(self.served, route)
        if resolved is None:
            self.send_error_response(404, f"Not found: {ANNDATA_ROUTE_PREFIX}{route}")
        elif isinstance(resolved, Path):
            self._send_file(resolved, head_only=head_only)
        else:
            self._send_bytes(resolved[0], resolved[1], head_only=head_only)

    def _range_interval(self, size: int, etag: str) -> tuple[int, int] | None | bool:
        """The requested interval, None for the whole resource, False when refused.

        A ``Range`` sent with an ``If-Range`` naming another version (a browser
        cache completing bytes it kept of an older file) is answered with the
        whole current resource, as HTTP requires, never with a range of a
        different version.
        """
        headers = self._request_headers()
        values = headers.get_all("Range", [])
        if len(values) > 1:
            return False
        condition = headers.get("If-Range")
        if values and condition is not None and condition.strip() != etag:
            return None
        try:
            return parse_byte_range(values[0] if values else None, size)
        except ValueError:
            return False

    def _not_modified(self, etag: str) -> bool:
        """Answer 304 when the browser already holds this version (``If-None-Match``)."""
        condition = self._request_headers().get("If-None-Match")
        if condition is None:
            return False
        tags = {tag.strip().removeprefix("W/") for tag in condition.split(",")}
        if etag not in tags and "*" not in tags:
            return False
        self.send_response(HTTPStatus.NOT_MODIFIED)
        self.send_header("ETag", etag)
        self.end_headers()
        return True

    def _send_head(
        self,
        *,
        size: int,
        interval: tuple[int, int] | None,
        etag: str,
    ) -> int:
        start, end = interval if interval is not None else (0, size - 1)
        length = 0 if size == 0 else end - start + 1
        self.send_response(HTTPStatus.PARTIAL_CONTENT if interval is not None else HTTPStatus.OK)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(length))
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("ETag", etag)
        if interval is not None:
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.end_headers()
        return length

    def _send_file(self, path: Path, *, head_only: bool) -> None:
        """Answer one file (or one range of it) with the kernel copying its bytes.

        The body goes from the page cache to the socket with ``sendfile``
        (``socket.sendfile``; a plain copy loop where the platform has none):
        no read into Python, no copy per request, nothing held in memory.
        """
        try:
            file = open(path, "rb", buffering=0)  # noqa: SIM115 - closed below
        except OSError:
            self.send_error_response(404, f"{path.name} is no longer readable")
            return
        with file:
            size, etag = file_version(file.fileno())
            if self._not_modified(etag):
                return
            interval = self._range_interval(size, etag)
            if interval is False:
                self.send_range_not_satisfiable(size)
                return
            assert interval is not True
            length = self._send_head(size=size, interval=interval, etag=etag)
            if head_only or length == 0:
                return
            sent = self.connection.sendfile(file, 0 if interval is None else interval[0], length)
            if sent != length:
                # The file shrank after its size was sent: the response
                # cannot be completed, so the connection is closed.
                self.close_connection = True

    def _send_bytes(self, body: bytes, etag: str, *, head_only: bool) -> None:
        if self._not_modified(etag):
            return
        interval = self._range_interval(len(body), etag)
        if interval is False:
            self.send_range_not_satisfiable(len(body))
            return
        assert interval is not True
        length = self._send_head(size=len(body), interval=interval, etag=etag)
        if not head_only and length:
            start = 0 if interval is None else interval[0]
            self._response_writer().write(body[start : start + length])

    def log_message(self, format: str, *args):
        """Override to use Python logging (formatted only when it is shown)."""
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("%s - %s", self.address_string(), format % args)


# =============================================================================
# Server
# =============================================================================


def _describe_layout(store: AnnDataStore) -> str:
    label = str(input_formats_label(store.kind))
    if store.layout == "zip":
        return f"{label}, zipped"
    if store.layout == "folder" and store.zarr_format is not None:
        return f"{label}, Zarr v{store.zarr_format} folder"
    if store.files is not None:
        return f"{label}, a folder of {len(store.files)} {plural(len(store.files), 'file')}"
    return label


def _identity_subject(kind: str) -> str:
    """Why a dataset of format ``kind`` has no identity unless one is given.

    AnnData and MuData declare theirs in ``uns['cellucid']``; every other
    format is read by the browser alone, so the server cannot take one from it.
    """
    if kind in ANNDATA_KINDS or kind in MUDATA_KINDS:
        noun = "MuData" if kind in MUDATA_KINDS else "AnnData"
        return (
            f"this {noun} declares no Cellucid identity of its own "
            "(uns['cellucid']['id'] and ['name'])"
        )
    return (
        f"a {input_formats_label(kind)} is read by the browser alone, so the server "
        "cannot take a Cellucid identity from it"
    )


def input_formats_label(kind: str) -> str:
    """The registry's label of one format id (``"Loom (.loom)"``)."""
    from ._anndata_files import input_formats

    for entry in input_formats()["formats"]:
        if entry["id"] == kind:
            return str(entry["label"])
    raise ValueError(f"Unknown input format {kind!r}.")


class AnnDataServer:
    """
    Server that hands one dataset (a file or folder Cellucid reads, or an
    AnnData or MuData object) to the Cellucid viewer.

    The viewer reads the data itself, with range requests where the format
    allows them; this server names
    the dataset, serves the store's bytes, carries the launch choices, and
    hosts the notebook channel.

    Examples
    --------
    Start a blocking server::

        AnnDataServer(
            "/path/to/data.h5ad",
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
        data: str | Path | anndata.AnnData | mudata.MuData,
        port: int = DEFAULT_PORT,
        host: str = DEFAULT_HOST,
        open_browser: bool = False,
        quiet: bool = False,
        *,
        latent_key: str | None = None,
        gene_id_column: str | None = None,
        normalize_embeddings: bool | None = None,
        centroid_min_points: int | None = None,
        dataset_name: str | None = None,
        dataset_id: str | None = None,
        obs_keys: Sequence[str] | None = None,
        embedding_keys: Mapping[str, EmbeddingKeySpec] | None = None,
        knn_key: str | None = None,
        expression_source: str | None = None,
        vector_field_default: str | None = None,
        modality: str | None = None,
        serve_web_ui: bool = True,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
        allowed_hosts: Sequence[str] | None = None,
    ) -> None:
        """
        Check the store and the launch choices; the server binds on ``start``.

        Parameters
        ----------
        data : str, Path, AnnData or MuData
            A file or folder in any format Cellucid reads (``.h5ad``, a zipped
            Zarr store, a Zarr store folder, ``.h5mu`` or MuData Zarr,
            ``.loom``, a 10x ``.h5`` file or ``outs/`` folder, a MatrixMarket
            folder, ``.rds``/``.rda``/``.RData``, gzipped or not), or an
            AnnData or MuData object. A leading ``~`` in a path is expanded to
            the user's home directory. An AnnData object, and an AnnData,
            MuData, Loom or 10x file whose cells are not in priority order,
            are served as a Zarr store in priority order made from them as
            the viewer reads it (:mod:`cellucid._served_store`); a file
            ``prepare()`` wrote, an R object and a MatrixMarket folder are
            served as they are on disk; a MuData object is written once to a
            temporary Cellucid-ready ``.h5mu``, removed when the server stops.
        port : int
            Port to serve on.
        host : str
            Host to bind to.
        open_browser : bool
            Whether to open browser on start.
        quiet : bool
            Suppress info messages.
        latent_key : str, optional
            ``obsm`` array (a latent space such as scVI or PCA) the viewer's
            outlier filter measures in, instead of the store's own
            ``uns['cellucid']['defaults']['latentKey']``. The viewer reads the
            outlier quantiles the store holds for it, or computes them on first
            use.
        gene_id_column : str, optional
            Exact column in ``var`` (``raw.var`` with ``expression_source="raw"``)
            whose values name the genes. If None, genes are named by the index.
        normalize_embeddings : bool, optional
            Whether the viewer scales embeddings into its coordinate range.
        centroid_min_points : int, optional
            The fewest cells a category needs to get outlier quantiles (10 when
            None).
        dataset_name : str, optional
            Exact non-empty human-readable name without surrounding whitespace
            or control characters. Unicode is preserved. Defaults to
            ``uns['cellucid']['name']``; required when the store has none.
        dataset_id : str, optional
            Portable 1-180 byte ASCII identifier beginning with a letter or
            digit and otherwise using only letters, digits, ``.``, ``_``, or
            ``-``; it cannot end in ``.`` or use a Windows device name.
            Sessions and community annotations are keyed by it. Defaults to
            ``uns['cellucid']['id']``; required when the store has none.
        obs_keys : sequence of str, optional
            Exact ``obs`` columns the viewer shows, in this order. If None,
            every column is shown.
        embedding_keys : mapping, optional
            The embedding of a dimension slot, overriding the default for that
            slot: ``{"2d": "X_tsne"}`` reads the first two columns of
            ``obsm['X_tsne']``, and ``{"2d": ("X_pca", [2, 0])}`` reads the 0-based
            columns named. Slots not named follow the default rule
            (``X_umap_1d``/``X_umap_2d``/``X_umap_3d``, else the first 1- to
            3-column of ``X_umap``, ``X_UMAP``, ``X_wnn.umap`` and
            ``X_ref.umap``); when no slot resolves, the viewer makes a deterministic
            best guess from the other ``obsm`` arrays and reports it.
        knn_key : str, optional
            The ``obsp`` matrix the viewer draws as the neighbour graph, usually
            ``"connectivities"``. None follows the viewer's default rule. The
            graph is read only when a view draws it.
        expression_source : str, optional
            Where gene expression is read: ``"X"``, ``"raw"`` (``raw.X``, with
            genes named from ``raw.var``), or ``"layer:<name>"``; None follows
            the viewer's default rule (``X`` when present, else ``raw`` or a
            layer, reported as guessed).
        vector_field_default : str, optional
            The vector field the overlay opens on: ``<field>`` of an obsm
            array ``<field>_<basis>`` drawn over ``X_<basis>`` (``velocity``
            for ``velocity_umap``). One of its embeddings must be shown when
            the store opens (a chosen or default embedding, or a 1-3 column
            array that is not the latent space).
        modality : str, optional
            The modality a multimodal dataset (a MuData, a Seurat object with
            several assays) opens on, one of its modality names; None opens
            the default one (the first named as gene expression -- ``rna``,
            ``gex``, ``gene expression``, ``gene_expression``, ``sct`` -- else
            the first). Set only
            for multimodal data.
        serve_web_ui : bool
            Establish and serve the exact current web build.
        web_source_url : str
            Origin publishing the web asset inventory.
        web_cache_dir : str or Path, optional
            Directory holding the active verified web build.
        allowed_hosts : sequence of str, optional
            Extra ``Host`` names this server answers to, needed only behind a
            reverse proxy such as ``jupyter-server-proxy``. Each entry is one
            bare DNS name or IP literal — no port, no scheme, no wildcard — and
            is matched on any port. Declaring names also applies the ``Host``
            check to a non-loopback bind.
        """
        self.dataset_name: str = ""
        self.dataset_id: str = ""
        given_name = (
            None if dataset_name is None else _require_dataset_name(dataset_name, label="dataset_name")
        )
        given_id = None if dataset_id is None else _require_dataset_id(dataset_id, label="dataset_id")
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
        self._server: CellucidThreadingHTTPServer | None = None
        self._thread: threading.Thread | None = None
        self._running = False
        self._started = False
        self._closed = False
        self._serving = False
        self._background_error: BaseException | None = None
        #: Releases what serving holds (a store's open file and its files, a
        #: temporary file); idempotent.
        self._release: Callable[[], Any] | None = None
        #: Whether a taken ``port`` is passed over for the next free one
        #: (:class:`AnnDataViewer` sets it when no port was given).
        self.next_free_port = False
        #: The program serving the viewer, as ``/_cellucid/info`` names it and
        #: the viewer shows it (the shared protocol's ``runtimes``): "python",
        #: or "jupyter" for a notebook viewer (:class:`AnnDataViewer` sets it).
        self.runtime = "python"

        try:
            self._prepare(
                data,
                given_id,
                given_name,
                latent_key=latent_key,
                gene_id_column=gene_id_column,
                normalize_embeddings=normalize_embeddings,
                centroid_min_points=centroid_min_points,
                obs_keys=obs_keys,
                embedding_keys=embedding_keys,
                knn_key=knn_key,
                expression_source=expression_source,
                vector_field_default=vector_field_default,
                modality=modality,
            )
        except BaseException:
            self._closed = True
            self._release_store()
            raise

    def _resolve_identity(
        self,
        given_id: str | None,
        given_name: str | None,
        declared: tuple[object, object],
        *,
        subject: str,
    ) -> None:
        """Name the dataset: as given, else as ``uns['cellucid']`` declares it.

        ``subject`` says why no identity was found, in the words of what was
        handed over (:func:`_identity_subject`).
        """
        declared_id, declared_name = declared
        if given_id is None and declared_id is not None:
            given_id = _require_dataset_id(declared_id, label="uns['cellucid']['id']")
        if given_name is None and declared_name is not None:
            given_name = _require_dataset_name(declared_name, label="uns['cellucid']['name']")
        missing = [
            label
            for label, value in (("dataset_name", given_name), ("dataset_id", given_id))
            if value is None
        ]
        if missing:
            raise ValueError(
                f"{' and '.join(missing)} {'are' if len(missing) > 1 else 'is'} required: "
                f"{subject}. Name the dataset, for example dataset_name=\"My study\", "
                "dataset_id=\"my-study-v1\" (--dataset-name/--dataset-id on the command line)."
            )
        assert given_id is not None and given_name is not None
        self.dataset_id = given_id
        self.dataset_name = given_name

    def _prepare(
        self,
        data: str | Path | anndata.AnnData | mudata.MuData,
        given_id: str | None,
        given_name: str | None,
        **choices: Any,
    ) -> None:
        quiet = self.quiet
        if not quiet:
            print_step(1, 4, "Detecting format")
        if isinstance(data, str | Path):
            store = classify_anndata_path(data)
            if not quiet:
                print_detail("Path", str(store.path))
                print_detail("Format", _describe_layout(store))
            self._resolve_identity(
                given_id,
                given_name,
                read_cellucid_identity(store),
                subject=_identity_subject(store.kind),
            )
            # Every HDF5 and Zarr format is served through the store in
            # priority order when the browser cannot decode one of its
            # codecs: the kernel decodes each chunk it is asked for.
            undecodable = browser_undecodable(store)
            if undecodable and not quiet:
                shown = ", ".join(undecodable[:3])
                more = f" and {len(undecodable) - 3} more" if len(undecodable) > 3 else ""
                print_detail("Browser cannot decode", f"{shown}{more}")
            if store.kind in ANNDATA_KINDS:
                self._serve_anndata_file(store, undecodable=bool(undecodable))
            elif store.kind in PRIORITY_SOURCE_KINDS:
                self._serve_foreign_file(store, undecodable=bool(undecodable))
            else:
                # Read once from the start by the browser (an R object, a
                # MatrixMarket folder): served as it is.
                self._served = _ServedStore(_route_name(store, self.dataset_id), store.kind, store.layout, store=store)
        else:
            import anndata as ad

            mudata = is_mudata(data)
            if not mudata and not isinstance(data, ad.AnnData):
                raise TypeError(
                    "data must be a path to a file or folder Cellucid reads, or an AnnData "
                    f"or MuData object; got {type(data).__name__}."
                )
            if not quiet:
                print_detail("Source", f"in-memory {'MuData' if mudata else 'AnnData'}")
                print_detail(
                    "Format",
                    "written once to a temporary .h5mu, in priority order"
                    if mudata
                    else "a Zarr store in priority order, made as the browser reads it",
                )
            uns = data.uns  # type: ignore[union-attr]
            declared = uns.get("cellucid") if hasattr(uns, "get") else None
            self._resolve_identity(
                given_id,
                given_name,
                (declared.get("id"), declared.get("name")) if isinstance(declared, Mapping) else (None, None),
                subject=_identity_subject("h5mu" if mudata else "h5ad"),
            )
            if mudata:
                self._serve_temporary_mudata(data)
            else:
                from ._served_store import PriorityStore

                self._serve_priority(PriorityStore.of_anndata(data))  # type: ignore[arg-type]
        if not quiet:
            print_success("Format detected")
        served = self._served
        self.data_format = served.kind

        if not quiet:
            print_step(2, 4, "Reading metadata")
        # None for a format only the browser reads (an R object, a
        # MatrixMarket folder): what it holds is learned when the browser
        # opens it, and the launch choices are checked then.
        if served.priority is not None:
            contents: AnnDataContents | None = served.priority.contents()
        else:
            assert served.store is not None
            contents = read_anndata_contents(served.store)
        if contents is not None and contents.n_obs <= 0:
            noun = "MuData" if served.kind in MUDATA_KINDS else "AnnData"
            raise ValueError(f"{noun} datasets must contain at least one observation")
        self.contents = contents
        if not quiet:
            self._print_contents()
            print_success("Metadata read")

        if not quiet:
            print_step(3, 4, "Checking launch choices")
        self.defaults = launch_defaults(contents, **choices)
        generated = (
            consolidated_zarr_metadata(served.store)
            if served.store is not None and served.layout == "folder" and served.store.files is None
            else None
        )
        if not quiet:
            for name, value in self.defaults.items():
                print_detail(name, str(value))
            if generated is not None:
                print_detail("Consolidated metadata", f"built in memory ({generated[0]})")
            print_success("Launch choices checked")
        if generated is not None:
            self._served = replace(
                served,
                generated_metadata=generated,
                generated_etag=f'"{hashlib.sha256(generated[1]).hexdigest()[:32]}"',
            )

        from . import __version__

        self.server_info = {
            "version": __version__,
            "type": "anndata",
            "format": served.kind,
            "host": self.host,
            "port": self.port,
            # Null when only the browser can read the data (see above).
            "n_cells": None if contents is None else contents.n_obs,
            "n_genes": None if contents is None else contents.n_vars,
            # The program serving the viewer, which the viewer names: set
            # from ``runtime`` when the server starts.
            "runtime": self.runtime,
        }

    def _serve_anndata_file(self, store: AnnDataStore, *, undecodable: bool) -> None:
        """Serve an AnnData file: as it is when it is in priority order, else in priority order.

        A file in priority order whose codecs the browser cannot decode is
        served through the store in priority order too (in its own order),
        which decodes it on the kernel. A layout before anndata 0.8 is
        served as it is (the browser reads it), or, with such a codec, read
        by anndata into memory first.
        """
        from ._served_store import PriorityStore, file_order

        order = file_order(store)
        if order == "priority" and not undecodable or order == "older" and not undecodable:
            self._served = _ServedStore(_route_name(store, self.dataset_id), store.kind, store.layout, store=store)
            if not self.quiet:
                print_detail(
                    "Serving",
                    "as it is (cells in priority order)" if order == "priority" else "as it is",
                )
            return
        if order == "older":
            priority = PriorityStore.of_anndata(read_into_memory(store))
        else:
            priority = PriorityStore.of_file(store)
        if not self.quiet:
            print_detail(
                "Serving",
                "a Zarr store "
                + ("decoded on the kernel" if order == "priority" else "in priority order")
                + ", made as the browser reads it",
            )
        self._serve_priority(priority)

    def _serve_foreign_file(self, store: AnnDataStore, *, undecodable: bool) -> None:
        """Serve a MuData, Loom or Cell Ranger file in priority order (as the browser presents it).

        A MuData file in priority order already (``prepare()``'s) is served
        as it is, unless the browser cannot decode one of its codecs.
        """
        from ._served_store import PriorityStore, file_order

        if store.kind in MUDATA_KINDS:
            if file_order(store) == "priority" and not undecodable:
                self._served = _ServedStore(_route_name(store, self.dataset_id), store.kind, store.layout, store=store)
                if not self.quiet:
                    print_detail("Serving", "as it is (cells in priority order)")
                return
            priority = PriorityStore.of_mudata_file(store)
        elif store.kind == "loom":
            priority = PriorityStore.of_loom(store)
        else:
            priority = PriorityStore.of_cell_ranger(store)
        if not self.quiet:
            what = "a MuData Zarr store" if priority.layout == "mudata" else "a Zarr store"
            print_detail("Serving", f"{what} in priority order, made as the browser reads it")
        self._serve_priority(priority)

    def _serve_priority(self, priority: PriorityStore) -> None:
        """Serve a store in priority order; it is closed when the server stops or Python exits."""
        kind = "mudata-zarr" if priority.layout == "mudata" else "zarr"
        self._served = _ServedStore(f"{self.dataset_id}.zarr", kind, "folder", priority=priority)
        self._release = priority.close

    def _serve_temporary_mudata(self, mdata: Any) -> None:
        """Write ``mdata`` to the temporary Cellucid-ready ``.h5mu`` this server serves."""
        path = write_temporary_h5mu(
            mdata, dataset_id=self.dataset_id, dataset_name=self.dataset_name
        )
        # Removed on stop, and when Python exits without a stop.
        self._release = weakref.finalize(self, remove_temporary_file, path)
        store = classify_anndata_path(path)
        self._served = _ServedStore(_route_name(store, self.dataset_id), store.kind, store.layout, store=store)

    def _print_contents(self) -> None:
        contents = self.contents
        if contents is None:
            print_detail("Contents", "read by the browser when it opens the data")
            return
        print_detail("Cells", f"{contents.n_obs:,}")
        print_detail("Genes", f"{contents.n_vars:,}")
        print_detail("obsm", ", ".join(sorted(contents.obsm)) or "none")
        print_detail("obs fields", f"{len(contents.obs_columns)}")
        print_detail("obsp", ", ".join(sorted(contents.obsp)) or "none")
        print_detail("layers", ", ".join(sorted(contents.layers)) or "none")

    def answer_request(
        self, method: str, path: str, range_header: str | None
    ) -> tuple[int, dict[str, str], bytes]:
        """Answer one request for the dataset as the HTTP handler does, without a socket.

        This is how a notebook widget's viewer reads the dataset
        (:mod:`cellucid.jupyter._widget`): the routes of
        :class:`AnnDataRequestHandler` but the viewer channel and the web
        build -- the documents (``_cellucid/health``, ``info``, ``datasets``,
        ``protocol``) and the store's bytes, one ``Range`` interval at a time
        (``206`` with ``Content-Range``, ``416`` for an interval the resource
        does not have). ``path`` is relative to the base address; ``method``
        is ``GET`` or ``HEAD``. Returns the status, the headers and the body
        (empty for ``HEAD``).
        """
        if method not in ("GET", "HEAD"):
            return 405, {"Content-Type": "text/plain"}, b"Only GET and HEAD are answered"
        route = unquote(urlparse(path).path).lstrip("/")
        documents = {
            "_cellucid/health": self._health_document,
            "_cellucid/info": lambda: self.server_info,
            "_cellucid/datasets": self._datasets_document,
            "_cellucid/protocol": _protocol_capability_document,
        }
        if route in documents:
            try:
                document = documents[route]()
            except OSError as error:
                body = _unreadable(error).encode()
                return 404, {"Content-Type": "text/plain", "Content-Length": str(len(body))}, body
            body = json.dumps(document).encode("utf-8")
            headers = {"Content-Type": "application/json", "Content-Length": str(len(body))}
            return 200, headers, b"" if method == "HEAD" else body
        resolved = (
            _resolve_store_route(self._served, route[len(ANNDATA_ROUTE_PREFIX) :])
            if route.startswith(ANNDATA_ROUTE_PREFIX)
            else None
        )
        if resolved is None:
            body = f"Not found: {route}".encode()
            return 404, {"Content-Type": "text/plain", "Content-Length": str(len(body))}, body
        try:
            file = None if not isinstance(resolved, Path) else open(resolved, "rb", buffering=0)  # noqa: SIM115 - closed below
        except OSError:
            body = f"{route} is no longer readable".encode()
            return 404, {"Content-Type": "text/plain", "Content-Length": str(len(body))}, body
        try:
            if file is not None:
                size, etag = file_version(file.fileno())
            else:
                assert not isinstance(resolved, Path)
                size, etag = len(resolved[0]), resolved[1]
            try:
                interval = parse_byte_range(range_header, size)
            except ValueError:
                body = b"Requested byte range is not satisfiable"
                return 416, {
                    "Content-Type": "text/plain",
                    "Content-Length": str(len(body)),
                    "Content-Range": f"bytes */{size}",
                    "Accept-Ranges": "bytes",
                }, body
            start, end = interval if interval is not None else (0, size - 1)
            length = 0 if size == 0 else end - start + 1
            headers = {
                "Content-Type": "application/octet-stream",
                "Content-Length": str(length),
                "Accept-Ranges": "bytes",
                "ETag": etag,
            }
            if interval is not None:
                headers["Content-Range"] = f"bytes {start}-{end}/{size}"
            status_code = 206 if interval is not None else 200
            if method == "HEAD" or length == 0:
                return status_code, headers, b""
            if file is None:
                assert not isinstance(resolved, Path)
                return status_code, headers, resolved[0][start : start + length]
            file.seek(start)
            body = file.read(length)
            if len(body) != length:
                # The file shrank after its size was read.
                message = f"{route} changed while it was read".encode()
                return 409, {"Content-Type": "text/plain", "Content-Length": str(len(message))}, message
            return status_code, headers, body
        finally:
            if file is not None:
                file.close()

    @property
    def file_route(self) -> str:
        """The store's route, relative to the server's base URL."""
        suffix = "/" if self._served.layout == "folder" else ""
        return f"{ANNDATA_ROUTE_PREFIX}{self._served.route_name}{suffix}"

    def _datasets_document(self) -> dict[str, Any]:
        served = self._served
        file: dict[str, Any] = {
            "url": self.file_route,
            "kind": served.kind,
            "layout": served.layout,
        }
        store = served.store
        if store is not None and served.layout != "folder":
            # The file as it is now: its size, and the validator its bytes are
            # answered with, which names this version of it to the viewer
            # (what it keeps of the file is filed under it).
            file["size"], file["etag"] = path_version(store.path)
        elif store is not None and store.files is not None:
            # A browser cannot list a folder: its files, their exact sizes and validators.
            file["files"] = [
                {"path": name, "size": size, "etag": etag}
                for name, (size, etag) in ((name, path_version(store.path / name)) for name in store.files)
            ]
        return {
            "datasets": [
                {
                    "id": self.dataset_id,
                    "name": self.dataset_name,
                    "anndata": {"file": file, "defaults": dict(self.defaults)},
                }
            ]
        }

    def _health_document(self) -> dict[str, Any]:
        info = self.server_info
        return {
            "status": "ok",
            "type": "anndata",
            "version": info["version"],
            "format": info["format"],
            "n_cells": info["n_cells"],
            "n_genes": info["n_genes"],
        }

    def _release_store(self) -> None:
        release = self._release
        self._release = None
        if release is not None:
            release()

    @property
    def url(self) -> str:
        """Get the URL of the currently running server.

        A bind to every interface reports the loopback origin, because the
        wildcard is not an address anything can open. The addresses other
        machines use are reported separately by ``network_urls``.
        """
        if not self._running or self._server is None:
            raise RuntimeError("AnnDataServer URL is unavailable because the server is not running")
        if bind_host_is_wildcard(self.host):
            return loopback_origin_for_bind(self.host, self.port)
        return format_http_origin(self.host, self.port)

    @property
    def viewer_url(self) -> str:
        """Get the full URL to open the viewer."""
        return f"{self.url}{VIEWER_QUERY}"

    @property
    def network_urls(self) -> list[tuple[str, str]]:
        """Return the origins another machine can open, when any exist."""
        if not self._running or self._server is None:
            return []
        if not bind_host_is_wildcard(self.host):
            return []
        origin = machine_origin(self.port)
        if origin is None:
            return []
        return [("Machine URL:  ", origin), ("Viewer URL:   ", f"{origin}{VIEWER_QUERY}")]

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

            self._server = self._bind()
            self.port = require_server_port(self._server.server_address[1])
            self.server_info["port"] = self.port
            self.server_info["runtime"] = self.runtime
            # The documents name the bound port, so they are built after binding.
            self._server.RequestHandlerClass = self._handler()
            self._running = True

            if not self.quiet:
                print_success("Server ready")
                print_server_banner(
                    self.url,
                    self.viewer_url,
                    network_urls=self.network_urls,
                )

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

    def _bind(self) -> CellucidThreadingHTTPServer:
        """Bind ``port``; with ``next_free_port``, the next free port when it is taken.

        A notebook viewer asks for the default port, the only local origin
        Google sign-in accepts, and takes the next free one when another
        viewer holds it (as cellucid-r's ``cellucid_view()`` does), saying so.
        """
        import errno

        port = self.port
        while True:
            try:
                return CellucidThreadingHTTPServer((self.host, port), self._handler())
            except OSError as error:
                taken = error.errno in (errno.EADDRINUSE, errno.EACCES) or getattr(error, "winerror", None) == 10013
                if not (self.next_free_port and taken and 0 < port < 65535):
                    raise
                port += 1

    def _handler(self) -> Any:
        return partial(
            AnnDataRequestHandler,
            served=self._served,
            server_info=dict(self.server_info),
            health=self._health_document(),
            datasets=self._datasets_document,
            serve_web_ui=self.serve_web_ui,
            web_cache_dir=self.web_cache_dir,
            allowed_hosts=self.allowed_hosts,
        )

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
        """Close the socket and remove a temporary file after serving ended."""
        self._running = False
        server = self._server
        self._server = None
        self._closed = True
        try:
            if server is not None:
                server.server_close()
        finally:
            self._release_store()

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
        self._release_store()
        self._closed = True

    def start_background(self):
        """Start the server in a background thread."""
        self.start(blocking=False)

    def stop(self):
        """Stop the server and remove a temporary file it wrote."""
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
        self._release_store()
        self._closed = True

        if not self.quiet:
            console_print("AnnData server stopped")
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
    data: str | Path | anndata.AnnData | mudata.MuData,
    port: int = DEFAULT_PORT,
    host: str = DEFAULT_HOST,
    open_browser: bool = True,
    quiet: bool = False,
    *,
    latent_key: str | None = None,
    gene_id_column: str | None = None,
    normalize_embeddings: bool | None = None,
    centroid_min_points: int | None = None,
    dataset_name: str | None = None,
    dataset_id: str | None = None,
    obs_keys: Sequence[str] | None = None,
    embedding_keys: Mapping[str, EmbeddingKeySpec] | None = None,
    knn_key: str | None = None,
    expression_source: str | None = None,
    vector_field_default: str | None = None,
    modality: str | None = None,
    serve_web_ui: bool = True,
    web_source_url: str = CELLUCID_WEB_URL,
    web_cache_dir: str | Path | None = None,
    allowed_hosts: Sequence[str] | None = None,
) -> AnnDataServer:
    """
    Serve a dataset (a file, a folder, an AnnData or MuData object) to the viewer, which reads it itself.

    The viewer reads the data in place with range requests where the format
    allows it (only what a view needs is transferred), and once from the
    start for an R object or a MatrixMarket folder; nothing is converted
    first.

    Parameters
    ----------
    data : str, Path, AnnData or MuData
        A file or folder in any format Cellucid reads (``.h5ad``, AnnData
        Zarr folder or ``.zip``, ``.h5mu``, MuData Zarr, ``.loom``, a 10x
        ``.h5`` file or ``outs/`` folder, a MatrixMarket folder,
        ``.rds``/``.rda``/``.RData``, gzipped or not), or an AnnData or MuData
        object. A leading ``~`` in a path is expanded to the user's home
        directory. An AnnData object, and an AnnData, MuData, Loom or 10x
        file not in priority order, are served as a Zarr store in priority
        order, made as the viewer reads it; a MuData object is written once
        to a temporary ``.h5mu``, removed when the server stops.
    port : int
        Port to serve on (default: 8765).
    host : str
        Host to bind to (default: 127.0.0.1).
    open_browser : bool
        Whether to open browser (default: True).
    quiet : bool
        Suppress info messages.
    latent_key : str, optional
        ``obsm`` array (a latent space such as scVI or PCA) the viewer's outlier
        filter measures in, instead of the store's own
        ``uns['cellucid']['defaults']['latentKey']``. The viewer reads the
        outlier quantiles the store holds for it, or computes them on first use.
    gene_id_column : str, optional
        Exact column in ``var`` (``raw.var`` with ``expression_source="raw"``)
        whose values name the genes. If None, genes are named by the index.
    normalize_embeddings : bool, optional
        Whether the viewer scales embeddings into its coordinate range.
    centroid_min_points : int, optional
        The fewest cells a category needs to get outlier quantiles (10 when
        None).
    dataset_name : str, optional
        Exact non-empty human-readable name without surrounding whitespace or
        control characters. Unicode is preserved. Defaults to
        ``uns['cellucid']['name']``; required when the store has none.
    dataset_id : str, optional
        Portable 1-180 byte ASCII identifier beginning with a letter or digit
        and otherwise using only letters, digits, ``.``, ``_``, or ``-``; it
        cannot end in ``.`` or use a Windows device name. Defaults to
        ``uns['cellucid']['id']``; required when the store has none.
    obs_keys : sequence of str, optional
        Exact ``obs`` columns the viewer shows, in this order. If None, every
        column is shown.
    embedding_keys : mapping, optional
        The embedding of a dimension slot, overriding the default for that
        slot: ``{"2d": "X_tsne"}`` reads the first two columns of
        ``obsm['X_tsne']``, and ``{"2d": ("X_pca", [2, 0])}`` reads the 0-based
        columns named.
    knn_key : str, optional
        The ``obsp`` matrix the viewer draws as the neighbour graph, usually
        ``"connectivities"``; read only when a view draws it.
    expression_source : str, optional
        Where gene expression is read: ``"X"``, ``"raw"`` or
        ``"layer:<name>"``; None follows the viewer's default rule.
    vector_field_default : str, optional
        The vector field the overlay opens on: ``<field>`` of an obsm array
        ``<field>_<basis>`` drawn over ``X_<basis>``, one of whose embeddings
        is shown when the store opens (as for ``AnnDataServer``).
    modality : str, optional
        The modality a multimodal dataset opens on (as for ``AnnDataServer``).
    serve_web_ui : bool
        Establish and serve the exact current web build.
    web_source_url : str
        Origin publishing the web asset inventory.
    web_cache_dir : str or Path, optional
        Directory holding the active verified web build.
    allowed_hosts : sequence of str, optional
        Extra ``Host`` names this server answers to, needed only behind a
        reverse proxy such as ``jupyter-server-proxy``. Each entry is one bare
        DNS name or IP literal — no port, no scheme, no wildcard — and is
        matched on any port.

    Returns
    -------
    AnnDataServer
        The server, after it stops (``start`` blocks until then).

    Example
    -------
    >>> from cellucid import serve_anndata
    >>> serve_anndata(
    ...     "/path/to/data.h5ad",
    ...     dataset_name="Example",
    ...     dataset_id="example",
    ... )

    >>> # Or an in-memory AnnData
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
        centroid_min_points=centroid_min_points,
        dataset_name=dataset_name,
        dataset_id=dataset_id,
        obs_keys=obs_keys,
        embedding_keys=embedding_keys,
        knn_key=knn_key,
        expression_source=expression_source,
        vector_field_default=vector_field_default,
        modality=modality,
        serve_web_ui=serve_web_ui,
        web_source_url=web_source_url,
        web_cache_dir=web_cache_dir,
        allowed_hosts=allowed_hosts,
    )
    server.start()
    return server
