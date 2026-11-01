"""
The viewer's Python side: :class:`AnnDataViewer` and :func:`show_anndata`.

:class:`AnnDataViewer` is the viewer API -- hooks, the event log that
``wait_for_event`` reads, ``viewer.state``, session bundle capture, the
outbound commands, shutdown -- and, on its own, a viewer in a browser tab: it
starts a :class:`cellucid.anndata_server.AnnDataServer` that hands the
AnnData store to the page (read in place with range requests; an in-memory
AnnData, or an AnnData, MuData, Loom or 10x file not in priority order, as a
store in priority order made as the page reads it) and serves the page and its
WebSocket channel. In a notebook, :func:`show_anndata` returns a
:class:`cellucid.jupyter.Widget` instead: the same API, shown in the cell's
output and reached through the notebook's own connection
(:mod:`cellucid.jupyter._widget`).

The module also owns the bridge to the viewer channel in
:mod:`cellucid._viewer_channel`: :func:`_register_viewer_for_messages` installs
a weak-referencing callback so a collected viewer never keeps a route alive,
and :func:`cleanup_all` stops whatever is still running at interpreter exit.
"""

from __future__ import annotations

import atexit
import copy
import html
import json
import logging
import os
import secrets
import sys
import threading
import time
import weakref
from collections import deque
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from .._console import console_print
from .._contracts import _require_dataset_id, _require_dataset_name
from .._plural import plural
from .._server_base import (
    CELLUCID_WEB_URL,
    DEFAULT_PORT,
    _extract_web_build_id,
    _web_cache_dir,
    require_allowed_hosts,
    require_server_port,
)
from .._viewer_channel import (
    ViewerProtocolError,
    cancel_session_bundle_request,
    register_event_callback,
    register_session_bundle_request,
    send_viewer_command,
    unregister_event_callback,
    viewer_is_connected,
)
from ..anndata_session._results import (
    PullSummary,
    Selection,
    ViewerResult,
    _gene_sets_of,
    apply_results,
    field_message,
    read_obs_names,
)
from ..data_mapping import _ANNDATA_MAPPING_KEYS, MODALITY_KEY, EmbeddingKeySpec
from ._hooks import HookCallback, HookRegistry, ViewerState
from ._wire import (
    _require_cell_indices,
    _require_client_server_url,
    _require_data_mapping,
    _require_debug_dataset_list,
    _require_exact_message_text,
    _require_frontend_message,
    _require_inbound_jupyter_event,
)

if TYPE_CHECKING:
    import anndata
    import mudata

    from ..anndata_server import AnnDataServer
    from ._widget import Widget

logger = logging.getLogger("cellucid.jupyter")

# The public method that sends each command, used to report a `command_error`
# in the notebook's own vocabulary: a viewer failure is far easier to act on
# when it names the call that caused it rather than the wire type. Every
# command of `_wire._COMMAND_TYPES` must appear here, and every value must name a
# real `AnnDataViewer` method; both are asserted by the contract tests.
_COMMAND_PYTHON_CALLERS: dict[str, str] = {
    "ping": "viewer.debug_connection()",
    "debug_snapshot": "viewer.debug_connection()",
    "requestSessionBundle": "viewer.request_session_bundle()",
    "highlight": "viewer.highlight_cells(...)",
    "clearHighlights": "viewer.clear_highlights()",
    "setColorBy": "viewer.set_color_by(...)",
    "setVisibility": "viewer.set_visibility(...)",
    "resetCamera": "viewer.reset_view()",
    "freeze": "viewer.stop()",
    "setDataMapping": "viewer.set_data_mapping(...)",
    "setField": "viewer.push(...)",
    "setGeneSet": "viewer.push_gene_set(...)",
}

# How many recent `command_error` payloads `debug_connection()` reports. The
# event log itself is already bounded; this keeps one diagnostic key readable.
_DEBUG_COMMAND_ERROR_LIMIT = 10


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



class _EventWaiter:
    """One ``wait_for_event`` call: what it matches, and the event it was handed."""

    __slots__ = ("event", "matches", "payload")

    def __init__(self, matches: Callable[[str, dict[str, Any]], bool]) -> None:
        self.matches = matches
        self.event: str | None = None
        self.payload: dict[str, Any] | None = None

class SessionBundleRequest:
    """A session bundle the viewer was asked for (:meth:`AnnDataViewer.request_session_bundle`).

    ``done`` says whether the viewer answered; :meth:`result` returns the
    bundle. A request that is never collected is withdrawn by :meth:`cancel`,
    or when the viewer stops.
    """

    def __init__(self, viewer: AnnDataViewer, request_id: str, waiter: _EventWaiter) -> None:
        self._viewer = viewer
        self.request_id = request_id
        self._waiter: _EventWaiter | None = waiter
        self._bundle: Any = None

    @property
    def done(self) -> bool:
        """Whether the viewer answered (with the bundle, or with why it could not)."""
        return self._bundle is not None or (self._waiter is not None and self._waiter.payload is not None)

    def result(self, timeout: float | None = 60.0):
        """
        The bundle, waiting up to ``timeout`` seconds for it (None: forever).

        Returns the :class:`~cellucid.session_bundle.CellucidSessionBundle`
        (the same one on every call). Raises ``TimeoutError`` when it did not
        arrive in time, and ``RuntimeError`` when the viewer could not capture
        it -- or, in a notebook widget, when it has not arrived yet: the
        viewer answers between cells, so collect it in a later cell.
        """
        from ..session_bundle import CellucidSessionBundle

        if self._bundle is not None:
            return self._bundle
        waiter = self._waiter
        if waiter is None:
            raise RuntimeError("This session bundle request was cancelled.")
        viewer = self._viewer
        viewer._await(waiter, timeout, "Timed out waiting for the session bundle")
        self._finish()
        assert waiter.payload is not None
        if waiter.event == "command_error":
            raise RuntimeError(
                "The viewer did not capture the session bundle: " + waiter.payload["reason"]
            )
        path = waiter.payload.get("path")
        if not path:
            raise RuntimeError("Session bundle upload succeeded but no temp path was returned")
        self._bundle = CellucidSessionBundle(Path(path))
        return self._bundle

    def cancel(self) -> None:
        """Withdraw the request; an answer that still arrives is dropped."""
        if self._waiter is None or self._bundle is not None:
            return
        self._finish()

    def _finish(self) -> None:
        waiter = self._waiter
        viewer = self._viewer
        with viewer._event_cv:
            if waiter in viewer._waiters:
                viewer._waiters.remove(waiter)
        if waiter is not None and waiter.payload is None:
            self._waiter = None
        if viewer._message_routing_registered:
            cancel_session_bundle_request(viewer._viewer_id, viewer._viewer_token, self.request_id)

    def __repr__(self) -> str:
        return f"SessionBundleRequest(done={self.done})"


def _probe_served_store(url: str) -> dict[str, Any]:
    """What a served AnnData store answers the probe the viewer opens it with.

    ``HEAD`` for its size and identity, then one byte with ``Range``: a store
    the viewer can read answers ``206`` with ``Content-Range``. A Zarr folder
    is probed at its root metadata.
    """
    import urllib.error
    import urllib.request

    target = url
    if url.endswith("/"):
        # Zarr v2 consolidated metadata, else a v3 root.
        target = f"{url}.zmetadata"
        try:
            with urllib.request.urlopen(urllib.request.Request(target, method="HEAD"), timeout=2):
                pass
        except urllib.error.HTTPError as error:
            if error.code == 404:
                target = f"{url}zarr.json"
        except Exception:
            pass
    probe: dict[str, Any] = {"url": url, "probed": target}
    try:
        head = urllib.request.Request(target, method="HEAD")
        with urllib.request.urlopen(head, timeout=2) as response:
            probe["head"] = {
                "status": response.status,
                "content_length": response.headers.get("Content-Length"),
                "accept_ranges": response.headers.get("Accept-Ranges"),
                "etag": response.headers.get("ETag"),
            }
        ranged = urllib.request.Request(target, headers={"Range": "bytes=0-0"})
        with urllib.request.urlopen(ranged, timeout=2) as response:
            response.read()
            probe["range"] = {
                "status": response.status,
                "content_range": response.headers.get("Content-Range"),
            }
    except Exception as error:
        probe["error"] = str(error)
    return probe


class AnnDataViewer:
    """
    Interactive viewer for AnnData objects, in a browser tab of its own.

    This machine serves the Cellucid page and the data, and the page talks
    back over a WebSocket; open :attr:`viewer_url` in a browser (in a
    notebook, :meth:`display` shows a link to it). The viewer reads the
    AnnData store itself, in place, with range requests: only what a view
    needs is read, and nothing is converted first. Supports bidirectional
    communication via hooks.

    To show the viewer inside a notebook cell, use :func:`show_anndata` or
    :class:`cellucid.jupyter.Widget`: the same API, reached through the
    notebook's own connection.

    Supports every file and folder Cellucid reads (``.h5ad`` and AnnData
    Zarr, ``.h5mu`` and MuData Zarr, ``.loom``, 10x ``.h5`` and ``outs/``,
    MatrixMarket folders, ``.rds``/``.rda``/``.RData``, gzipped or not),
    served as it is on disk -- an AnnData, MuData, Loom or 10x file not in
    priority order, and an in-memory AnnData or MuData, as a store in priority
    order made as the page reads it (an object's expression matrices laid out
    gene-major in the background), removed when the viewer stops.

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
        >>> # From an .h5ad file, read in place
        >>> viewer = AnnDataViewer(
        ...     "/path/to/data.h5ad",
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ... )

    Event Hooks:
        Register callbacks to respond to user interactions in the viewer::

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

            @viewer.on_command_error
            def handle(event):
                print(event['command'], event['reason'])  # a command failed

            @viewer.on_data_mapping
            def handle(event):
                print(event['origin'], event['mapping'])  # what the viewer plots

            @viewer.on_message
            def handle(event):
                print(event)  # raw message dict
    """

    def __init__(
        self,
        data: str | Path | anndata.AnnData | mudata.MuData,
        port: int | None = None,
        auto_open: bool = True,
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
        client_server_url: str | None = None,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
        allowed_hosts: Sequence[str] | None = None,
        cache_bytes: int | None = None,
    ) -> None:
        """
        Initialize the AnnData viewer.

        Args:
            data: AnnData or MuData object, or path to a file or folder
                Cellucid reads (see above). A leading ``~`` in a path is
                expanded to the user's home directory.
            port: Port for the data server. None takes 8765 -- the one local
                port Cellucid registers with Google, so Google Drive files
                open (community annotation's GitHub sign-in works only on
                www.cellucid.com) -- or, when another
                viewer holds it, the next free port, said in a printed note
                (cellucid-r's ``cellucid_view()`` does the same). A port given
                is used or refused.
            auto_open: Show where the viewer is when created: a link in a
                notebook (:meth:`display`), its address on the console.
            latent_key: ``obsm`` array (a latent space such as scVI or PCA)
                the viewer's outlier filter measures in, instead of the data's
                own ``uns['cellucid']['defaults']['latentKey']``. The viewer
                reads the outlier quantiles the data stores for it, or
                computes them on first use.
            gene_id_column: Exact column in ``var`` (``raw.var`` with
                ``expression_source="raw"``) whose values name the genes. If
                None, genes are named by the index.
            normalize_embeddings: Whether the viewer scales embeddings into its
                coordinate range.
            centroid_min_points: The fewest cells a category needs to get
                outlier quantiles (10 when None).
            dataset_name: Exact non-empty human-readable name without
                surrounding whitespace or control characters. Defaults to
                ``uns['cellucid']['name']``, else the file's name, or
                ``AnnData (<n> cells)`` for an object.
            dataset_id: Portable 1-180 byte ASCII identifier beginning with a
                letter or digit and otherwise using letters, digits, ``.``,
                ``_``, or ``-``; no trailing ``.`` or Windows device name.
                Defaults to ``uns['cellucid']['id']``, else made from the file's
                name without its extension, or from the object's label (``AnnData (N cells)``) (said in
                one line: sessions are kept under it).
            obs_keys: Exact ``obs`` columns the viewer shows, in this order. If
                None, every column is shown.
            embedding_keys: The embedding of a dimension slot, such as
                ``{"2d": "X_tsne"}`` or ``{"2d": ("X_pca", [2, 0])}`` with 0-based
                columns. Slots not named follow the default rule, and a
                deterministic best guess applies when no slot resolves.
            knn_key: The ``obsp`` matrix the viewer draws as the neighbour
                graph, such as ``"connectivities"``; read only when a view draws
                it. None follows the viewer's default rule.
            expression_source: ``"X"``, ``"raw"``, or ``"layer:<name>"``;
                None follows the viewer's default rule (``raw`` when it holds
                more genes than ``X`` or ``X`` is scaled, as Scanpy reads
                genes; else ``X``).
            vector_field_default: The vector field the overlay opens on:
                ``<field>`` of an obsm array ``<field>_<basis>`` drawn over
                ``X_<basis>`` (``velocity`` for ``velocity_umap``).
            modality: The modality a multimodal dataset (a MuData) opens on;
                None opens the default one. Set only for multimodal data.
            client_server_url: Exact browser-reachable data server base URL.
            web_source_url: Origin publishing the web asset inventory.
            web_cache_dir: Directory holding the active verified web build.
            allowed_hosts: Extra ``Host`` names the data server answers to.
                Required when ``client_server_url`` points at a reverse proxy
                such as ``jupyter-server-proxy``: name that proxy's host here.
            cache_bytes: What the data server keeps of what it makes for a
                store served in priority order, in bytes (least recently used
                first out); None keeps a 64th of the machine's memory, 0
                nothing (as :func:`cellucid.serve_anndata`'s).
        """
        if dataset_name is not None:
            dataset_name = _require_dataset_name(dataset_name, label="dataset_name")
        if dataset_id is not None:
            dataset_id = _require_dataset_id(dataset_id, label="dataset_id")
        self._init_state(
            port=port,
            client_server_url=client_server_url,
            web_source_url=web_source_url,
            web_cache_dir=web_cache_dir,
            allowed_hosts=allowed_hosts,
        )

        try:
            self.data = data

            self._start_server(
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
                cache_bytes=cache_bytes,
            )
            self._activate()

            if auto_open:
                self.display()
        except BaseException:
            self._rollback_failed_construction()
            raise

    def _init_state(
        self,
        *,
        port: int | None,
        client_server_url: str | None,
        web_source_url: str,
        web_cache_dir: str | Path | None,
        allowed_hosts: Sequence[str] | None,
    ) -> None:
        """Set up everything but the data server: ports, ids, hooks, events."""
        self.port = DEFAULT_PORT if port is None else require_server_port(port)
        self._port_given = port is not None
        self.allowed_hosts = require_allowed_hosts(allowed_hosts)
        self._client_server_url = (
            _require_client_server_url(client_server_url) if client_server_url is not None else None
        )
        self.web_source_url = web_source_url
        self.web_cache_dir = (
            Path(web_cache_dir).expanduser().resolve()
            if web_cache_dir is not None
            else _web_cache_dir()
        )

        self._server: AnnDataServer | None = None
        self._init_viewer_state()

    def _init_viewer_state(self) -> None:
        """Set up what every viewer has, however it is shown: ids, hooks, events."""
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
        # The latest events, for the debug report only. Waiters never read it:
        # each registers before it waits and is handed its event as the event
        # is recorded, so a burst of other events cannot push it out.
        self._recent_events: deque[tuple[int, str, dict[str, Any]]] = deque(maxlen=512)
        self._waiters: list[_EventWaiter] = []
        # Set by stop(): every wait still pending ends at once with this reason.
        self._stopped_reason: str | None = None
        # Every result the viewer holds in its current epoch (one page's
        # account of one dataset it opened), by (kind, id) (`results`); the
        # keys of the epoch's complete set, None until its snapshot ended with
        # a `results` inventory; and the host dataset that inventory named
        # (None: the person opened another dataset in the viewer).
        self._results_epoch: str | None = None
        self._results: dict[tuple[str, str], ViewerResult] = {}
        self._results_inventory: set[tuple[str, str]] | None = None
        self._results_dataset: str | None = None
        # The id of the dataset this viewer serves, and its cell count when
        # this process reads it (None: only the browser reads the format),
        # kept past stop().
        self._served_dataset_id: str | None = None
        self._served_cell_count: int | None = None
        # What waits for the viewer, sent once it can be: columns pushed
        # before its cell count was known (`ready`), and data-mapping fields
        # set before it reported its mapping.
        self._pending_pushes: list[tuple[Any, str]] = []
        self._pending_mapping_fields: dict[str, Any] = {}

    def _start_server(
        self,
        *,
        latent_key: str | None,
        gene_id_column: str | None,
        normalize_embeddings: bool | None,
        centroid_min_points: int | None,
        dataset_name: str | None,
        dataset_id: str | None,
        obs_keys: Sequence[str] | None,
        embedding_keys: Mapping[str, EmbeddingKeySpec] | None,
        knn_key: str | None,
        expression_source: str | None,
        vector_field_default: str | None,
        modality: str | None,
        cache_bytes: int | None,
    ) -> None:
        """Start the AnnData server."""
        from ..anndata_server import AnnDataServer

        self._server = AnnDataServer(
            data=self.data,
            port=self.port,
            host="127.0.0.1",
            open_browser=False,
            quiet=True,
            serve_web_ui=True,
            web_source_url=self.web_source_url,
            web_cache_dir=self.web_cache_dir,
            allowed_hosts=sorted(self.allowed_hosts),
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
            cache_bytes=cache_bytes,
        )
        # The viewer names the program serving it from /_cellucid/info: a
        # notebook kernel's viewer is "Jupyter", any other Python's "Python".
        self._server.runtime = "jupyter" if self._context["in_jupyter"] else "python"
        self._server.next_free_port = not self._port_given
        self._server.start_background()
        if not self._port_given and self._server.port != DEFAULT_PORT:
            console_print(
                f"Port {DEFAULT_PORT} is in use; serving on {self._server.port}. Google Drive "
                f"files open only on port {DEFAULT_PORT}."
            )
        self.port = self._server.port
        self._served_dataset_id = self._server.dataset_id
        cell_count = self._server.server_info["n_cells"]
        self._served_cell_count = cell_count if isinstance(cell_count, int) else None
        logger.info(f"Started AnnData server at {self._server.url}")

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
            - cells: list[int] - indices of selected cells, unique and ascending
            - source: str - selection method ('lasso', 'knn', 'proximity', 'annotation', ...)

        Example::

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

        Example::

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

        Example::

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

        Example::

            @viewer.on_ready
            def handle(event):
                print(f"Loaded {event['n_cells']} cells")
        """
        return self._hooks.on("ready")

    @property
    def on_command_error(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a handler for a failed Python → viewer command.

        Called when a command this viewer sent -- ``set_color_by``,
        ``highlight_cells``, ``set_visibility`` and the rest -- was rejected by
        the viewer instead of applied. The commands are one-way and carry no
        request id, so the command method itself has already returned ``None``
        by the time the failure arrives: this hook is how a notebook learns
        that the viewer is unchanged.

        Registering a handler also silences the notice Cellucid otherwise
        prints to ``sys.stderr``, so reporting happens once, where you chose.

        Event data:
            - command: str - the command type that failed, e.g. ``"setColorBy"``
            - reason: str - what the viewer reported, in one line

        Like every hook, this runs on the viewer channel's thread, not in your
        cell: an exception raised here is logged (logger
        ``cellucid.viewer_channel``) and never reaches the notebook, the
        callbacks after it do not run, and the event is not recorded in
        ``viewer.state`` or handed to ``wait_for_event``. Record the failure
        and raise from your own code.
        There is no acknowledgement for a command that did apply, so an empty
        record is not proof of success -- only a non-empty one is proof of
        failure.

        Example::

            rejected = []

            @viewer.on_command_error
            def handle(event):
                rejected.append(event)

            viewer.set_color_by("cell_type")
            ...
            if rejected:
                raise RuntimeError(rejected[-1]["reason"])
        """
        return self._hooks.on("command_error")

    @property
    def on_data_mapping(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a handler for the viewer's data mapping.

        Called when the viewer reports which stored arrays it plots and reads:
        once after the dataset loads, and again whenever the mapping changes --
        from the viewer's controls, a session, or
        :meth:`set_data_mapping`. The viewer sends it because this package
        lists ``data_mapping`` among the events of ``GET /_cellucid/protocol``.

        Event data:
            - origin: str - why the mapping is what it is: ``"default"``,
              ``"best_guess"``, ``"session"``, ``"user"``, or
              ``"notebook"`` (an answer to :meth:`set_data_mapping`)
            - requestId: str or None - the request a ``"notebook"`` answer
              completes; None for every other origin
            - mapping: dict - the complete data mapping, in canonical form

        Example::

            @viewer.on_data_mapping
            def handle(event):
                print(event["mapping"]["embeddings"])
        """
        return self._hooks.on("data_mapping")

    @property
    def on_message(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a raw message handler.

        Called for every validated current-schema message from the viewer.
        Use it to observe the complete event stream while debugging.

        Event data:
            - event: str - the event type
            - ... other fields depend on event type

        Example::

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
            event: Event name ('selection', 'hover', 'click', 'ready',
                'command_error', 'data_mapping', 'message')
            callback: Function to call when event fires

        Returns:
            The callback function

        Example::

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
        try:
            event, data = _require_inbound_jupyter_event(message)
        except (TypeError, ValueError) as error:
            raise ViewerProtocolError(str(error)) from error
        # The event is recorded before any hook runs: a hook that raises is the
        # notebook's error (the channel logs it), and must not cost
        # viewer.state, wait_for_event or a session request the event.
        self._record_event(event, data)
        try:
            self._hooks.trigger(event, data)
        finally:
            if event == "command_error" and not self._hooks.has_callbacks("command_error"):
                self._report_unhandled_command_error(data)
        # What waited for this event goes now.
        if event == "ready":
            self._send_pending_pushes()
        elif event == "data_mapping":
            self._send_pending_mapping()

    def _report_local_command_error(self, command: str, error: Exception) -> None:
        """A command that waited for the viewer and could not be built once it
        could: reported as the viewer reports a refused command (a
        ``command_error``: ``state.command_error``, the hooks, or printed)."""
        reason = str(error).strip() or type(error).__name__
        self._handle_frontend_message({"type": "command_error", "command": command, "reason": reason[:240]})

    def _report_unhandled_command_error(self, data: dict[str, Any]) -> None:
        """
        Print one failed command to ``sys.stderr`` when nothing is listening.

        Commands are one-way: ``viewer.set_color_by(...)`` returns ``None``
        whether the viewer applied it or rejected it, and the rejection arrives
        on the viewer channel's thread long after the cell finished. Nothing can make the
        call itself raise without blocking every command on a browser round
        trip, so a notebook that registered no ``command_error`` hook would
        otherwise see a successful-looking call and an unchanged viewer.

        This runs after the event has been recorded, so a failure to write does
        not cost the notebook the event.
        """
        command = data["command"]
        console_print(
            f"[cellucid] {_COMMAND_PYTHON_CALLERS[command]} did not take effect "
            f"in the viewer: {data['reason']}\n"
            "[cellucid] The viewer is unchanged. Handle this in code with "
            "@viewer.on_command_error.",
            file=sys.stderr,
        )

    def _record_event(self, event: str, data: dict[str, Any]):
        """Update `viewer.state` and notify any `wait_for_event(...)` callers."""
        now = time.monotonic()
        with self._event_cv:
            self._event_seq += 1
            seq = self._event_seq
            self._recent_events.append((seq, event, data))
            for waiter in self._waiters:
                if waiter.payload is None and waiter.matches(event, data):
                    waiter.event = event
                    waiter.payload = data

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
            elif event == "command_error":
                self.state.command_error = data
            elif event == "data_mapping":
                self.state.data_mapping = data
            elif event in {"result", "result_removed", "results"}:
                self._record_result_event(event, data)

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

        Returns
        -------
        dict
            The payload of the first matching event that arrives after the
            call (an event recorded before it is not returned).

        Raises
        ------
        TimeoutError
            No matching event arrived within ``timeout``.
        """
        self._require_answerable(f"wait_for_event({event!r})")
        waiter = _EventWaiter(
            lambda ev_type, payload: ev_type == event
            and (predicate is None or predicate(payload))
        )
        with self._event_cv:
            self._waiters.append(waiter)
        try:
            self._await(waiter, timeout, f"Timed out waiting for event '{event}'")
        finally:
            with self._event_cv:
                self._waiters.remove(waiter)
        assert waiter.payload is not None
        return waiter.payload

    def _require_answerable(self, what: str) -> None:
        """Refuse, before anything is sent, a call this thread cannot wait for.

        A browser-tab viewer answers on its server's threads, so every thread
        can wait; a notebook widget overrides this (its answers arrive only
        between cells).
        """

    def _await(self, waiter: _EventWaiter, timeout: float | None, message: str) -> None:
        """Block until ``waiter`` is handed its event; raise TimeoutError at the
        deadline, or RuntimeError at once when the viewer is stopped."""
        deadline = None if timeout is None else (time.monotonic() + float(timeout))
        with self._event_cv:
            while waiter.payload is None:
                if self._stopped_reason is not None:
                    what = message.removeprefix("Timed out waiting for ")
                    raise RuntimeError(f"Stopped waiting for {what}: {self._stopped_reason}")
                if deadline is not None and time.monotonic() >= deadline:
                    raise TimeoutError(message)
                remaining = None if deadline is None else max(0.0, deadline - time.monotonic())
                self._event_cv.wait(timeout=remaining)

    def wait_for_ready(self, timeout: float | None = 30.0) -> dict[str, Any]:
        """Return the viewer's ``ready`` payload, waiting for it if it has not arrived.

        Raises ``TimeoutError`` when it does not arrive within ``timeout``
        seconds (None waits forever).
        """
        if self.state.ready is not None:
            return self.state.ready
        return self.wait_for_event("ready", timeout=timeout)

    def request_session_bundle(self) -> SessionBundleRequest:
        """
        Ask the viewer for its current `.cellucid-session` bundle, without waiting.

        The request is sent at once (or as soon as the viewer connects) and
        the viewer answers in the background; collect the bundle with the
        returned request's ``result()``. In a notebook widget the answer
        arrives once the cell that asked has finished, so collect it in a later
        cell::

            request = viewer.request_session_bundle()
            # next cell
            bundle = request.result()
            adata = bundle.apply_to_anndata(adata)

        Returns
        -------
        SessionBundleRequest
            The request; ``request.done`` says whether the answer arrived.
        """
        request_id = secrets.token_hex(16)
        register_session_bundle_request(
            self._viewer_id,
            self._viewer_token,
            request_id,
            ttl_seconds=None,
        )
        # The viewer answers with the bundle, or with a command_error naming
        # the request when it could not make one (or the channel refused
        # what it sent): either completes the request.
        waiter = _EventWaiter(
            lambda ev_type, payload: (
                ev_type == "session_bundle" and payload["requestId"] == request_id
            )
            or (ev_type == "command_error" and payload["command"] == "requestSessionBundle")
        )
        with self._event_cv:
            self._waiters.append(waiter)
        try:
            self.send_message(
                {
                    "type": "requestSessionBundle",
                    "requestId": request_id,
                }
            )
        except BaseException:
            with self._event_cv:
                self._waiters.remove(waiter)
            cancel_session_bundle_request(self._viewer_id, self._viewer_token, request_id)
            raise
        return SessionBundleRequest(self, request_id, waiter)

    @property
    def data_mapping(self) -> dict[str, Any] | None:
        """The data mapping the viewer last reported, or None before it reported one.

        This is a copy of the ``mapping`` of the latest ``data_mapping`` event:
        which embedding fills each dimension slot, the bound vector fields,
        and -- for an AnnData source -- the neighbour graph, the gene names,
        the expression matrix and the latent space. Editing the copy changes
        nothing; use :meth:`set_data_mapping`.
        """
        with self._state_lock:
            event = self.state.data_mapping
            return None if event is None else copy.deepcopy(event["mapping"])

    @property
    def data_mapping_origin(self) -> str | None:
        """Why the viewer's data mapping is what it is, or None before it reported one."""
        with self._state_lock:
            event = self.state.data_mapping
            return None if event is None else str(event["origin"])

    def set_data_mapping(self, **fields: Any) -> None:
        """
        Change what the viewer plots and reads.

        Each keyword replaces one top-level field of the mapping the viewer
        reports (:attr:`data_mapping`): ``embeddings``, ``vectorFields``,
        ``geneNames``, ``knn``, ``expression``, ``log1p``, ``highlyVariable``
        and ``latent`` (and ``modality`` for a MuData). ``kind`` cannot
        change.

        Nothing is waited for -- a notebook widget hears the viewer only
        between cells. The merged mapping is validated and sent at once when
        the viewer has reported its mapping (an invalid value raises here and
        the viewer is untouched), or as soon as it reports one (after its
        dataset opens; a cell run before that, "Run All", is fine). The
        viewer's answer becomes :attr:`data_mapping` (origin ``"notebook"``);
        a mapping it refuses -- an array that fails its load-time check, for
        example -- arrives as a ``command_error`` (``state.command_error``,
        ``@viewer.on_command_error``, printed when no hook handles it).

        Example::

            viewer.set_data_mapping(
                embeddings={
                    "enabled": [{"key": "X_pca", "columns": [0, 1]}],
                    "defaults": {
                        "1d": None,
                        "2d": {"key": "X_pca", "columns": [0, 1]},
                        "3d": None,
                    },
                },
                knn="connectivities",
            )
            # a later cell
            viewer.data_mapping["knn"]
        """
        if "kind" in fields:
            raise TypeError("set_data_mapping cannot change the mapping kind.")
        allowed = (_ANNDATA_MAPPING_KEYS - {"kind"}) | {MODALITY_KEY}
        unknown = sorted(set(fields) - allowed)
        if unknown:
            raise TypeError(
                f"set_data_mapping got unknown field(s) {unknown}; a data "
                f"mapping has {sorted(allowed)}."
            )
        if not fields:
            raise TypeError("Name the fields to change, for example expression='layer:counts'.")
        with self._state_lock:
            self._pending_mapping_fields.update(copy.deepcopy(fields))
            reported = self.state.data_mapping is not None
        if reported:
            self._send_pending_mapping(raise_errors=True)

    def _send_pending_mapping(self, *, raise_errors: bool = False) -> None:
        """Send the data-mapping fields set so far, merged into the mapping the
        viewer last reported. A merged mapping that is invalid raises when the
        call that set the fields sends it (``raise_errors``), and is reported
        as a ``command_error`` when it waited for the viewer."""
        with self._state_lock:
            fields = self._pending_mapping_fields
            event = self.state.data_mapping
            if not fields or event is None:
                return
            current = copy.deepcopy(event["mapping"])
            self._pending_mapping_fields = {}
        try:
            unknown = sorted(set(fields) - set(current))
            if unknown:
                raise TypeError(f"This dataset's data mapping has no field(s) {unknown}")
            mapping = _require_data_mapping({**current, **fields}, label="data mapping")
        except (TypeError, ValueError) as error:
            if raise_errors:
                raise
            self._report_local_command_error("setDataMapping", error)
            return
        self.send_message({"type": "setDataMapping", "requestId": secrets.token_hex(16), "mapping": mapping})

    # =========================================================================
    # ROUND TRIP: results back from the viewer, columns and gene sets pushed
    # =========================================================================
    #
    # Everything here answers from what the viewer already sent, never by
    # asking it: once its dataset opens (and whenever its channel opens
    # again) the viewer sends every result and the `results` inventory, then
    # each change as it happens. A notebook widget hears the viewer only
    # between cells, so a lasso drawn and `viewer.pull(adata)` run in the next
    # cell reads the lasso; a browser-tab viewer is read the same way.

    def _record_result_event(self, event: str, data: dict[str, Any]) -> None:
        """Keep the viewer's results current (called holding the event lock).

        An event of another epoch -- a reloaded page, the viewer shown again
        in another output, another dataset opened in it -- starts afresh: its
        snapshot follows. Within an epoch a ``result`` replaces an older
        revision of the same result, a ``result_removed`` drops it, and a
        ``results`` inventory (a snapshot's end) is the complete set: it
        drops every result it does not list, and names the dataset.
        """
        if data["epoch"] != self._results_epoch:
            self._results_epoch = data["epoch"]
            self._results = {}
            self._results_inventory = None
            self._results_dataset = None
        if event == "result":
            result = ViewerResult.from_event(data)
            key = (result.kind, result.id)
            held = self._results.get(key)
            if held is None or held.revision <= result.revision:
                self._results[key] = result
            if self._results_inventory is not None:
                self._results_inventory.add(key)
        elif event == "result_removed":
            key = (str(data["kind"]), str(data["id"]))
            held = self._results.get(key)
            if held is not None and held.revision < data["revision"]:
                del self._results[key]
                if self._results_inventory is not None:
                    self._results_inventory.discard(key)
        else:
            listed = {(item["kind"], item["id"]) for item in data["items"]}
            for key in [key for key in self._results if key not in listed]:
                del self._results[key]
            self._results_inventory = listed
            self._results_dataset = data["dataset"]

    @property
    def results(self) -> list[ViewerResult]:
        """Every result the viewer holds, as of the last event it sent: highlight
        groups, user-defined fields (labels, clusters), scores, gene sets, the
        visible and the open cells.

        Raises ``RuntimeError`` before the viewer has sent them all once (it
        does when its dataset has opened): run the cell again then; and while
        the viewer shows a dataset this viewer did not serve (one opened in
        it instead), whose cells are not this dataset's.
        """
        with self._state_lock:
            inventory = self._results_inventory
            if inventory is None:
                raise RuntimeError(
                    "The viewer has not sent what it holds yet: it does once its dataset has "
                    "opened. Run this again when the viewer shows the data."
                )
            served = self._served_dataset_id
            if self._results_dataset is None:
                raise RuntimeError(
                    f"The viewer shows a dataset this notebook did not open: another one was opened in it. "
                    f"Open {served!r} in the viewer again (or show the viewer again) to pull from it."
                )
            if served is not None and self._results_dataset != served:
                raise RuntimeError(
                    f"The viewer's results describe the dataset {self._results_dataset!r}, and this viewer "
                    f"serves {served!r}: open {served!r} in the viewer again to pull from it."
                )
            missing = sorted(f"{kind} {result_id}" for kind, result_id in inventory if (kind, result_id) not in self._results)
            if missing:
                raise RuntimeError(
                    f"The viewer's results arrived incomplete ({', '.join(missing[:3])} missing); "
                    "reload the viewer's page to have it send them again."
                )
            return list(self._results.values())

    def pull(
        self,
        adata: Any,
        *,
        highlights: bool = True,
        fields: bool = True,
        scores: bool = True,
        visibility: bool = False,
        gene_sets: bool = True,
        include_pushed: bool = False,
        highlight_prefix: str = "cellucid_highlight__",
        field_prefix: str = "",
        on_conflict: str = "error",
        prune: bool = True,
    ) -> PullSummary:
        """Write everything made in this viewer onto ``adata``, in place, in one call.

        Every highlight group, user-defined field (labels, annotations,
        clusters), score and gene set the viewer holds -- see
        :func:`cellucid.pull` for the columns, their names and dtypes, and
        every option. ``adata`` is the AnnData or MuData the viewer shows,
        cells in the same order. Nothing is waited for: the viewer sent its
        results as they were made, so in a notebook a lasso drawn and this run
        in the next cell finds the lasso.

        Example::

            summary = viewer.pull(adata)
            summary.columns          # what was written, by obs column
            adata.uns["cellucid"]["gene_sets"]
        """
        self._require_viewer_order(adata)
        return apply_results(
            adata,
            self.results,
            highlights=highlights,
            fields=fields,
            scores=scores,
            visibility=visibility,
            gene_sets=gene_sets,
            include_pushed=include_pushed,
            highlight_prefix=highlight_prefix,
            field_prefix=field_prefix,
            on_conflict=on_conflict,
            prune=prune,
        )

    def _require_viewer_order(self, adata: Any) -> None:
        """Refuse an object whose rows are not the viewer's cells in order.

        The viewer names a cell by its original position (:attr:`obs_names`
        are in that order): an object's own row when it was shown in its
        original order; for an object or a file in a verified priority order,
        the original position its record gives the row. An object holding the
        cells in another order -- read from such a file as stored
        (``anndata.read_zarr`` of ``prepare()``'s output), or sorted since --
        would take each result in the wrong rows, so it is refused: by its
        verified record (the order it is in, whatever its names), and by its
        names, which must be the viewer's in order (repeated names
        included). When this process cannot read the viewer's names (a
        format only the browser reads), the record alone tells.
        """
        from ..prepare_data._priority import stored_order_of

        # In a verified stored order (`stored_order`) that moves a row: a
        # record that does not describe the cells leaves each row its own
        # original position.
        positions = stored_order_of(adata).positions if hasattr(adata, "uns") and hasattr(adata, "obs") else None
        stored = positions is not None and not np.array_equal(positions, np.arange(len(positions)))
        names = self.obs_names
        given = pd.Index(adata.obs_names)
        if names is not None:
            if len(given) != len(names):
                return  # refused with the counts
            if given.equals(names) and not stored:
                return
            if names.is_unique and given.is_unique:
                missing = names.get_indexer(given) < 0
                if missing.any():
                    raise ValueError(
                        f"{int(missing.sum())} of this object's cells are not cells of the viewer's dataset "
                        f"(first: {given[missing][0]!r}): pass the object the viewer shows."
                    )
        elif not stored:
            return
        restore = (
            "Put them back first: adata = cellucid.restore_original_order(adata)."
            if stored
            else "Reorder it first: adata = adata[viewer.obs_names].copy()."
        )
        raise ValueError(
            "This object holds the viewer's cells in another order"
            + (" (a Cellucid-ready file's stored order: uns['cellucid']['priority'])" if stored else "")
            + ", and the viewer names cells in their original order. "
            + restore
        )

    @property
    def obs_names(self) -> pd.Index | None:
        """The obs names of the cells this viewer shows, in the order it names
        them (their original order), or None when this process cannot read
        them (a format only the browser reads).

        Taken from what is served (:attr:`AnnDataServer.obs_names
        <cellucid.anndata_server.AnnDataServer.obs_names>`): an object's
        names as they were when it was shown -- its cells put back in their
        original order when it was read from a file in priority order --
        whatever is done to it after; a file's at the version it is served
        as. A viewer that stopped keeps the last ones."""
        server = self._server
        if server is not None:
            names = server.obs_names
            self._obs_names_cache = names
            return names
        cached = getattr(self, "_obs_names_cache", False)
        if cached is False:
            cached = read_obs_names(self.data)
            self._obs_names_cache = cached
        return cached

    @property
    def selection(self) -> Selection | None:
        """The viewer's last selection (lasso, k-NN, proximity, annotation), or
        None before the first: its cells as rows of ``obs`` and, when this
        process has them, their obs names."""
        with self._state_lock:
            event = self.state.selection
        if event is None:
            return None
        cells = np.asarray(event["cells"], dtype=np.int64)
        names = self.obs_names
        return Selection(cells=cells, obs_names=None if names is None else names[cells], source=str(event["source"]))

    def _cell_positions(self, cells: Any) -> np.ndarray | list[int]:
        """Cells given by row (ints) or by obs name, as rows of ``obs``."""
        if isinstance(cells, str | bytes):
            raise TypeError("cell_indices must list the cells (rows or obs names), not be one string")
        if isinstance(cells, pd.Index | pd.Series) or (
            isinstance(cells, list | tuple) and cells and all(isinstance(cell, str) for cell in cells)
        ) or (isinstance(cells, np.ndarray) and cells.dtype.kind in "OUS"):
            names = self.obs_names
            if names is None:
                raise ValueError(
                    "This viewer's dataset has no obs names this process can read; name the cells by row."
                )
            wanted = pd.Index(cells.index if isinstance(cells, pd.Series) and cells.dtype == bool else cells)
            if isinstance(cells, pd.Series) and cells.dtype == bool:
                wanted = cells.index[cells.to_numpy()]
            if not names.is_unique:
                raise ValueError(
                    "This viewer's dataset repeats obs names, so a name does not say which cell it is: "
                    "name the cells by row, or make the names unique (adata.obs_names_make_unique()) and show it again."
                )
            positions = names.get_indexer(wanted)
            if (positions < 0).any():
                missing = wanted[positions < 0]
                raise ValueError(
                    f"{len(missing)} of the cells named are not in this viewer's dataset (first: {missing[0]!r})"
                )
            located: np.ndarray = positions.astype(np.int64)
            return located
        return _require_cell_indices(cells if isinstance(cells, np.ndarray) or type(cells) is list else list(cells))

    def subset(self, adata: Any, which: str | None = None) -> Any:
        """The cells of a highlight group, or of the last selection, as a new AnnData.

        Args:
            adata: The AnnData (or MuData) the viewer shows, cells in its
                order (an object holding them in another order is refused,
                as :meth:`pull` refuses it).
            which: A highlight group, by id (``"highlight_3"``) or by its
                label; ``"visible"`` for the cells the filters leave visible;
                None for the last selection.

        Returns:
            ``adata[cells].copy()`` (``to_memory()`` for a backed AnnData).
        """
        self._require_viewer_order(adata)
        if which is None:
            selection = self.selection
            if selection is None:
                raise RuntimeError("The viewer has sent no selection yet; select cells, or name a highlight group.")
            cells = selection.cells
        else:
            results = self.results
            if which == "visible":
                matches = [result for result in results if result.kind == "visibility"]
            else:
                matches = [
                    result for result in results
                    if result.kind == "highlight" and (result.id == which or result.meta["label"] == which)
                ]
            if len(matches) != 1:
                groups = sorted(
                    f"{result.id} ({result.meta['label']})" for result in results if result.kind == "highlight"
                )
                raise ValueError(
                    f"{which!r} names {len(matches)} highlight groups; name one by id. The viewer holds: "
                    + (", ".join(groups) or "no highlight group")
                )
            if matches[0].n_obs != int(adata.n_obs):
                raise ValueError(
                    f"The viewer shows {matches[0].n_obs} {plural(matches[0].n_obs, 'cell')} and this object has "
                    f"{int(adata.n_obs)}: pass the object the viewer shows, with its cells in the same order"
                )
            found = matches[0].cells
            cells = np.arange(int(adata.n_obs)) if found is None else np.asarray(found, dtype=np.int64)
        view = adata[cells]
        return view.to_memory() if getattr(adata, "isbacked", False) else view.copy()

    def gene_sets(self, *, symbols: bool = False) -> dict[str, list[str]]:
        """Every gene set the viewer holds: the notebook's, the shown marker genes
        result's (one per group) and the shown gene signature's, by name (a name
        used again is ``"name (2)"``, as everywhere in Cellucid).

        Each set's genes are identifiers -- the var index, as
        ``adata.var_names`` holds them -- or, with ``symbols=True``, the
        names the viewer shows for them: those of the var column its data
        mapping names the genes by (such as ``feature_name``), else the
        identifiers themselves. :func:`cellucid.write_gmt` writes them as a
        GMT file."""
        if type(symbols) is not bool:
            raise TypeError("symbols must be exactly True or False")
        return {
            name: list(entry["symbols"] if symbols and "symbols" in entry else entry["genes"])
            for _, name, entry in _gene_sets_of(self.results)
        }

    def push(self, values: Any, key: str | None = None) -> None:
        """Show a column in the viewer as a field: labels as a categorical, numbers
        as a continuous field. Pushing the same key again replaces it.

        Args:
            values: A ``pandas.Series`` -- ``adata.obs["leiden"]``, a score --
                aligned to the viewer's cells by its index (obs names), or an
                array or list of one value per cell in ``obs`` order.
                Categorical, string, object and bool values become a
                categorical (missing values have no category; a bool's
                labels are ``"false"`` and ``"true"``); numbers a
                continuous field (missing: no value).
            key: The field's name in the viewer; the Series' name when None.
                A field of that name the notebook did not push is never
                replaced: the viewer refuses it, and says why as a
                ``command_error`` (printed unless a hook handles it).

        The column is sent at once (or, for a format only the browser reads,
        once the viewer has said how many cells it has: a cell run before
        that, "Run All", is fine); the viewer adds it once its dataset is
        open, and :attr:`results` lists it from then on (in a notebook, from
        the next cell). A column the viewer refuses, or that does not fit
        the dataset once it is known, arrives as a ``command_error``.
        """
        name = key if key is not None else getattr(values, "name", None)
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Name the field: pass key=, or a Series with a name")
        n_obs = self._dataset_cell_count()
        if n_obs is None:
            # Only the browser reads this format: the column goes once the
            # viewer says how many cells it has (its `ready`).
            with self._state_lock:
                self._pending_pushes.append((values, name.strip()))
            return
        message = field_message(values, name.strip(), n_obs=n_obs, obs_names=self.obs_names)
        self.send_message({**message, "requestId": secrets.token_hex(16)})

    def _send_pending_pushes(self) -> None:
        """Send the columns pushed before the viewer's cell count was known; a
        column that does not fit the dataset is reported as a
        ``command_error`` of ``setField``."""
        with self._state_lock:
            pending, self._pending_pushes = self._pending_pushes, []
        n_obs = self._dataset_cell_count()
        if n_obs is None:
            return
        for values, name in pending:
            try:
                message = field_message(values, name, n_obs=n_obs, obs_names=self.obs_names)
            except (TypeError, ValueError) as error:
                self._report_local_command_error("setField", error)
                continue
            self.send_message({**message, "requestId": secrets.token_hex(16)})

    def push_gene_set(self, name: str, genes: Sequence[str], description: str = "") -> None:
        """Give the viewer a gene set: Gene signature and Dot plot offer it
        ("From notebook"), and it is saved with the session. A set of the same
        name the notebook pushed before is replaced. Each gene is kept once, in
        order. Sent at once, as :meth:`push` is."""
        if isinstance(genes, str):
            raise TypeError("genes must list the genes, not be one string")
        clean = list(dict.fromkeys(str(gene).strip() for gene in genes if str(gene).strip()))
        self.send_message(
            {"type": "setGeneSet", "requestId": secrets.token_hex(16), "name": name, "description": description, "genes": clean}
        )

    def _dataset_cell_count(self) -> int | None:
        """The cells of the dataset this viewer serves: read by this process,
        or reported by the viewer (its ``ready``); None before either."""
        if self._served_cell_count is not None:
            return int(self._served_cell_count)
        ready = self.state.ready
        return None if ready is None else int(ready["n_obs"])

    def apply_session_to_anndata(
        self,
        adata: Any,
        request: SessionBundleRequest,
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
        Apply the session a :meth:`request_session_bundle` request brought to an AnnData.

        The one way to get a session, in a notebook widget and a browser-tab
        viewer alike (cellucid-r's ``cellucid_apply_session(object,
        request)``): ask, then apply -- in a notebook, in a later cell, since
        the viewer answers between cells::

            request = viewer.request_session_bundle()
            # next cell
            adata = viewer.apply_session_to_anndata(adata, request)

        ``timeout`` bounds the wait for the bundle (``request.result()``; a
        notebook widget's cell cannot wait, and says so at once when the
        answer has not arrived), ``cleanup_bundle`` deletes the bundle file
        afterwards (``request.result().save(...)`` first keeps a copy), and
        ``expected_dataset_id`` defaults to the id this viewer serves. Every
        other keyword is passed to
        :func:`~cellucid.apply_cellucid_session_to_anndata`.
        """
        if type(cleanup_bundle) is not bool:
            raise TypeError("cleanup_bundle must be exactly True or False")
        if not isinstance(request, SessionBundleRequest) or request._viewer is not self:
            raise TypeError("request must be what this viewer's request_session_bundle() returned")
        bundle = request.result(timeout=timeout)
        try:
            resolved_dataset_id = expected_dataset_id
            if resolved_dataset_id is None:
                resolved_dataset_id = self._served_dataset_id
                if resolved_dataset_id is None:
                    raise RuntimeError("This viewer serves no dataset id; pass expected_dataset_id explicitly.")
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
        - exact id/path-keyed identity probes for every declared dataset
        - whether the viewer's channel (its WebSocket) is connected
        - Python→viewer ping/pong roundtrip over the channel
        - frontend debug snapshot (URL/origin/userAgent as seen by the iframe)
        - ``recent_events["command_errors"]``: the commands the viewer refused
          instead of applying
        """
        report: dict[str, Any] = {
            "viewer_id": self._viewer_id,
            "viewer_url": self.viewer_url,
            "server_url": self.server_url,
            "displayed": self._displayed,
            "notebook_context": dict(self._context),
            "server_running": self._server is not None and self._server.is_running(),
            "channel_connected": (
                self._message_routing_registered
                and viewer_is_connected(self._viewer_id, self._viewer_token)
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
            from ..web_cache import verify_web_ui_cache

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
                datasets = _require_debug_dataset_list(payload)
                report["server_datasets"] = datasets
        except Exception as e:
            report["server_datasets_error"] = str(e)

        # Probe every listed store at its exact route.
        if datasets is not None:
            report["anndata_store_probes"] = {
                dataset["id"]: _probe_served_store(
                    f"{self.server_url.rstrip('/')}/{dataset['anndata']['file']['url']}"
                )
                for dataset in datasets
            }

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
                command_errors: list[dict[str, Any]] = []
                for _seq, ev_type, payload in self._recent_events:
                    by_type[ev_type] = by_type.get(ev_type, 0) + 1
                    if ev_type == "command_error":
                        command_errors.append(payload)
                report["recent_events"] = {
                    "count": len(self._recent_events),
                    "by_type": by_type,
                    # The commands the viewer refused, oldest first. A command
                    # that failed is otherwise indistinguishable from one that
                    # applied and matched nothing.
                    "command_errors": command_errors[-_DEBUG_COMMAND_ERROR_LIMIT:],
                }
        except Exception as e:
            report["recent_events_error"] = str(e)

        # ------------------------------------------------------------------
        # Frontend roundtrip probe (requires the viewer iframe to be alive)
        # ------------------------------------------------------------------
        if not report["channel_connected"]:
            reason = (
                "The viewer's channel is not connected: display it (viewer.display()) "
                "or open viewer.viewer_url, and check that any proxy in between passes "
                "WebSockets"
            )
            report["frontend_roundtrip"] = {"ok": False, "error": reason}
            report["frontend_debug_snapshot"] = {"ok": False, "error": reason}
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
        from ..web_cache import clear_web_cache

        return clear_web_cache(cache_dir=self.web_cache_dir)

    def ensure_web_ui_cached(self, *, force: bool = True, show_progress: bool = True):
        """Establish one complete verified web build in this viewer's cache."""
        from ..web_cache import ensure_web_ui_cached

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
        if self._server is not None:
            return self._server.url
        return f"http://127.0.0.1:{self.port}"

    def _get_client_server_url(self) -> str:
        """Return the one browser URL selected for this viewer."""
        if self._client_server_url is not None:
            return self._client_server_url
        if self._server is None:
            raise RuntimeError("The Cellucid data server has not started")
        return _require_client_server_url(str(self._server.url))

    @property
    def viewer_url(self) -> str:
        """Get the full viewer URL with anndata flag."""
        from urllib.parse import urlencode

        base = self._get_client_server_url()
        query = urlencode(
            {
                "viewer": "true",
                "viewerId": self._viewer_id,
                "viewerToken": self._viewer_token,
                "anndata": "true",
            }
        )
        return f"{base}/?{query}"

    def display(self):
        """Say where the viewer is: a link to its browser tab in a notebook,
        its address on the console otherwise.

        This viewer runs in a browser tab of its own (the page and its channel
        are served by this machine). To show a viewer inside a notebook cell,
        use :func:`show_anndata` or :class:`cellucid.jupyter.Widget`.
        """
        if not self._context["in_jupyter"]:
            console_print(f"Open the Cellucid viewer: {self.viewer_url}")
            return

        from IPython.display import HTML, display

        from ..web_cache import verify_web_ui_cache

        verify_web_ui_cache(
            self.web_cache_dir,
            expected_source_url=self.web_source_url,
        )

        display(HTML(self._repr_html_()))
        self._displayed = True

    def send_message(self, message: dict) -> bool:
        """
        Send one command to the viewer.

        This is the low-level API for sending commands to the frontend. The
        command travels on the viewer's channel (a WebSocket the viewer opens
        to this server). Before the viewer has connected -- the cell output is
        still loading, or the page is reconnecting -- the command waits and is
        sent, in order, the moment it connects.

        A command that names cells (``highlight``, ``setVisibility``) travels
        as one binary frame, its cells little-endian ``uint32``; every other
        one as JSON text.

        Args:
            message: Dict with 'type' key and message-specific data. Cells are
                a list of ints or a one-dimensional NumPy integer array.

        Returns:
            True when the command was sent now, False when it is waiting for
            the viewer to connect.

        Example::

            viewer.send_message({'type': 'highlight', 'cells': [1,2,3], 'color': '#ff0000'})
        """
        serialized_message = _require_frontend_message(message)
        if not self._message_routing_registered:
            raise RuntimeError("The viewer has stopped; it accepts no more commands.")
        return send_viewer_command(self._viewer_id, self._viewer_token, serialized_message)

    # =========================================================================
    # PYTHON → FRONTEND CONVENIENCE METHODS
    # =========================================================================
    # These methods send commands TO the viewer (opposite direction from hooks),
    # on the same viewer channel the events arrive on.

    def highlight_cells(
        self,
        cell_indices: Any,
        color: str = "#ff0000",
        *,
        label: str | None = None,
    ):
        """
        Highlight specific cells in the viewer, as one highlight group.

        Args:
            cell_indices: The cells to highlight, in any order, each once: by
                0-based row of ``obs`` (a list of ints, or a one-dimensional
                NumPy integer array -- the fast form for millions of cells,
                e.g. ``np.flatnonzero(mask)``), or by obs name (a list of
                names, a ``pandas.Index`` such as ``adata.obs_names[mask]``).
                They travel as one binary frame. When the viewer opened a
                subsample, the cells it holds are highlighted.
            color: Hex color string (default: red)
            label: The group's name in the viewer; None names it by its
                cell count.

        Example::

            viewer.highlight_cells([100, 200, 300], color="#00ff00")
            viewer.highlight_cells(adata.obs_names[adata.obs["leiden"] == "3"], label="Leiden 3")
        """
        cells = self._cell_positions(cell_indices)
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
                "cells": cells,
                "color": color,
                "label": label,
            }
        )

    def clear_highlights(self):
        """Clear all cell highlights in the viewer: every group on every highlight page (the pages stay)."""
        self.send_message({"type": "clearHighlights"})

    def set_color_by(self, field: str):
        """
        Set the coloring field in the viewer.

        Args:
            field: Name of the obs column to color by

        Example::

            viewer.set_color_by("cell_type")
        """
        exact_field = _require_exact_message_text(
            field,
            label="field",
            allow_whitespace=True,
        )
        self.send_message({"type": "setColorBy", "field": exact_field})

    def set_visibility(self, cell_indices: list[int] | np.ndarray | None = None, visible: bool = True):
        """
        Set visibility of specific cells.

        Args:
            cell_indices: The cells, in any order, each once and at least
                one: by 0-based row of ``obs`` (a list of ints or a
                one-dimensional NumPy integer array) or by obs name. None for
                every cell.
            visible: Whether cells should be visible

        Example::

            viewer.set_visibility([0, 1, 2], visible=False)  # Hide cells
        """
        if type(visible) is not bool:
            raise TypeError("visible must be exactly True or False")
        self.send_message(
            {
                "type": "setVisibility",
                "cells": None if cell_indices is None else self._cell_positions(cell_indices),
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
        1. Stops the underlying data server (which removes its served store's scratch folder)
        2. Removes the viewer from the active viewers set
        3. Clears all hooks
        4. Logs the cleanup

        Safe to call multiple times after a successful stop.
        """
        failures: list[BaseException] = []
        with self._event_cv:
            self._stopped_reason = "the viewer was stopped"
            self._event_cv.notify_all()
        finalizer = self._routing_finalizer
        # At interpreter exit the routing finalizer may already have run
        # (weakref's own exit hook): the registration is gone, and so is the
        # page's channel, so there is nothing left to tell.
        if self._message_routing_registered and (finalizer is None or finalizer.alive):
            try:
                # Only an open page is told to freeze; a page that is not
                # connected sees its channel close instead.
                send_viewer_command(
                    self._viewer_id,
                    self._viewer_token,
                    _require_frontend_message({"type": "freeze"}),
                    queue=False,
                )
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

    def __repr__(self) -> str:
        status = "running" if self._server and self._server.is_running() else "stopped"
        return f"AnnDataViewer(port={self.port}, status={status})"

    def _repr_html_(self) -> str:
        """A link that opens the viewer in a browser tab."""
        return (
            f'<a href="{html.escape(self.viewer_url, quote=True)}" target="_blank" '
            'rel="noopener">Open the Cellucid viewer</a>'
        )


def show_anndata(
    data: str | Path | anndata.AnnData | mudata.MuData,
    height: int = 600,
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
    web_source_url: str = CELLUCID_WEB_URL,
    cache_bytes: int | None = None,
) -> Widget:
    """
    Show an AnnData or MuData object, or any file Cellucid reads, in a notebook.

    Returns a :class:`cellucid.jupyter.Widget` and shows it in the cell's
    output: in JupyterLab, Jupyter Notebook 7, VS Code, Google Colab and
    Databricks it is displayed at once; in marimo, make it the cell's output
    (the last expression). The widget reaches this kernel through the
    notebook's own connection, so no port, proxy or server address is
    involved. The viewer reads the data itself, without running prepare
    first: in place with range requests where the format allows it, once
    from the start for an R object or a MatrixMarket folder. An in-memory
    AnnData, and an AnnData, MuData, Loom or 10x file not in priority order,
    are served as a Zarr store with their cells in priority order, made as
    the viewer reads it, so a subsample reads only its cells; so is a MuData
    object. An object's expression matrices are laid out gene-major in the
    background (read-only meanwhile), so the viewer opens at once; a cell of
    the notebook run again stops the widget it showed before.

    To open a viewer in a browser tab of its own instead (a script, an SSH
    tunnel), use :class:`AnnDataViewer`, or ``serve_anndata`` /
    ``cellucid serve``.

    Args:
        data: AnnData or MuData object, or path to a file or folder Cellucid
            reads (see "Supported formats"). A leading ``~`` in a path is
            expanded to the user's home directory.
        height: Height of the viewer in pixels.
        latent_key: ``obsm`` array (a latent space such as scVI or PCA) the
            viewer's outlier filter measures in, instead of the data's own
            ``uns['cellucid']['defaults']['latentKey']``. The viewer reads the
            outlier quantiles the data stores for it, or computes them on first
            use.
        gene_id_column: Exact column in ``var`` (``raw.var`` with
            ``expression_source="raw"``) whose values name the genes. If None,
            genes are named by the index.
        normalize_embeddings: Whether the viewer scales embeddings into its
            coordinate range.
        centroid_min_points: The fewest cells a category needs to get outlier
            quantiles (10 when None).
        dataset_name: Exact non-empty human-readable name without surrounding
            whitespace or control characters. Defaults to
            ``uns['cellucid']['name']``, else the file's name, or
            ``AnnData (<n> cells)`` for an object.
        dataset_id: Portable 1-180 byte ASCII identifier beginning with a
            letter or digit and otherwise using letters, digits, ``.``, ``_``,
            or ``-``; no trailing ``.`` or Windows device name. Defaults to
            ``uns['cellucid']['id']``, else made from the file's name without
            its extension, or from the object's label (``AnnData (N cells)``) (said in one line:
            sessions are kept under it).
        obs_keys: Exact ``obs`` columns the viewer shows, in this order. If
            None, every column is shown.
        embedding_keys: The embedding of a dimension slot, such as
            ``{"2d": "X_tsne"}`` or ``{"2d": ("X_pca", [2, 0])}`` with 0-based
            columns. Slots not named follow the default rule, and a
            deterministic best guess applies when no slot resolves.
        knn_key: The ``obsp`` matrix the viewer draws as the neighbour graph,
            such as ``"connectivities"``; read only when a view draws it. None
            follows the viewer's default rule.
        expression_source: ``"X"``, ``"raw"``, or ``"layer:<name>"``; None
            follows the viewer's default rule (``raw`` when it holds more genes
            than ``X`` or ``X`` is scaled, as Scanpy reads genes; else ``X``).
        vector_field_default: The vector field the overlay opens on:
            ``<field>`` of an obsm array ``<field>_<basis>`` drawn over
            ``X_<basis>`` (``velocity`` for ``velocity_umap``).
        modality: The modality a multimodal dataset (a MuData) opens on; None
            opens the default one. Set only for multimodal data.
        web_source_url: Where the browser loads the Cellucid app the widget shows:
            ``https://www.cellucid.com`` (the default), or the address of a
            copy you host. The browser, not the kernel, loads it.
        cache_bytes: What the kernel keeps of what it makes for a store
            served in priority order, in bytes (least recently used first
            out); None keeps a 64th of the machine's memory, 0 nothing (as
            :func:`cellucid.serve_anndata`'s).

    Returns:
        The :class:`cellucid.jupyter.Widget`, for hooks and commands.

    Supported formats:
        - .h5ad, AnnData Zarr (folder or .zip), .h5mu, MuData Zarr, .loom,
          10x .h5 files and Cell Ranger outs/ folders: a file prepare() wrote
          read in place, any other served in priority order as the viewer
          reads it
        - MatrixMarket folders (a .mtx[.gz] file of any name, such as GEO's
          GSM..._matrix.mtx.gz) and R objects (.rds, .rda, .RData, gzipped or
          not), read once
        - In-memory AnnData and MuData objects, served in priority order as
          the viewer reads them, their expression matrices laid out gene-major
          in the background

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

    Note:
        To publish a dataset for others, prepare it (or publish the .h5ad or
        Zarr store itself) and share a link; see the data loading guide.
    """
    from ._widget import Widget, show_widget

    widget = Widget(
        data,
        height=height,
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
        web_source_url=web_source_url,
        cache_bytes=cache_bytes,
    )
    try:
        show_widget(widget)
    except BaseException:
        widget.stop()
        raise
    return widget


# =============================================================================
# THE VIEWER CHANNEL
# =============================================================================
#
# The viewer page opens one WebSocket to the data server (/_cellucid/ws) and
# authenticates with the viewer id and token from its address. Commands from
# Python and events from the viewer both travel on it
# (cellucid._viewer_channel); each event reaches _handle_frontend_message on
# the socket's thread, in order, and hooks fire there. The page reaches the
# server at the same address it was loaded from, so a notebook proxy that
# carries the page carries the channel too.
#
# =============================================================================


def _register_viewer_for_messages(viewer: AnnDataViewer):
    """
    Register a viewer to receive messages from frontend.

    The viewer's channel (its WebSocket) authenticates against this
    registration, and its events are delivered to the viewer through it.
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
    logger.debug(f"Registered the viewer channel of viewer {viewer._viewer_id}")


def _unregister_viewer_for_messages(viewer: AnnDataViewer):
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


def cleanup_all():
    """Stop all active viewers and their servers."""
    for viewer in list(_active_viewers):
        try:
            viewer.stop()
        except Exception:
            logger.exception("Failed to stop %r during interpreter shutdown", viewer)


# Register cleanup on interpreter exit
atexit.register(cleanup_all)
