"""
The notebook viewer: :class:`AnnDataViewer` and :func:`show_anndata`.

:class:`AnnDataViewer` starts a :class:`cellucid.anndata_server.AnnDataServer`,
which hands the AnnData store to the viewer: the viewer reads it in place with
range requests (an in-memory object is first written to a temporary ``.h5ad``).
The class also owns notebook-context detection, the viewer id/token pair, hook
dispatch, the event log that ``wait_for_event`` reads, session bundle capture,
iframe HTML generation, the outbound command surface, and shutdown.

The module also owns the bridge to the HTTP event router in
:mod:`cellucid._server_base`: :func:`_register_viewer_for_messages` installs a
weak-referencing callback so a collected viewer never keeps a route alive, and
:func:`cleanup_all` stops whatever is still running at interpreter exit.
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
from urllib.parse import urlparse

from .._console import console_print
from .._contracts import _require_dataset_id, _require_dataset_name
from .._server_base import (
    CELLUCID_WEB_URL,
    _extract_web_build_id,
    _web_cache_dir,
    cancel_session_bundle_request,
    register_event_callback,
    register_session_bundle_request,
    require_allowed_hosts,
    require_server_port,
    unregister_event_callback,
)
from ..data_mapping import EmbeddingKeySpec
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

    from ..anndata_server import AnnDataServer

logger = logging.getLogger("cellucid.jupyter")

# The public method that sends each command, used to report a `command_error`
# in the notebook's own vocabulary: a viewer failure is far easier to act on
# when it names the call that caused it rather than the wire type. Every key of
# `_wire._COMMAND_MESSAGE_FIELDS` must appear here, and every value must name a
# real `AnnDataViewer` method; both are asserted by the contract tests.
_COMMAND_PYTHON_CALLERS: dict[str, str] = {
    "ping": "viewer.debug_connection()",
    "debug_snapshot": "viewer.debug_connection()",
    "requestSessionBundle": "viewer.get_session_bundle()",
    "highlight": "viewer.highlight_cells(...)",
    "clearHighlights": "viewer.clear_highlights()",
    "setColorBy": "viewer.set_color_by(...)",
    "setVisibility": "viewer.set_visibility(...)",
    "resetCamera": "viewer.reset_view()",
    "freeze": "viewer.stop()",
    "setDataMapping": "viewer.set_data_mapping(...)",
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
    Interactive viewer for AnnData objects in Jupyter notebooks.

    The viewer reads the AnnData store itself, in place, with range requests:
    only what a view needs is read, and nothing is converted first. Supports
    bidirectional communication via hooks.

    Supports .h5ad files, zipped Zarr stores (.zip) and Zarr store folders,
    served as they are on disk, and in-memory AnnData objects, written once to
    a temporary .h5ad that is removed when the viewer stops.

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
        data: str | Path | anndata.AnnData,
        port: int | None = None,
        height: int = 600,
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
        client_server_url: str | None = None,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
        allowed_hosts: Sequence[str] | None = None,
    ) -> None:
        """
        Initialize the AnnData viewer.

        Args:
            data: AnnData object, or path to an .h5ad file, a zipped Zarr
                store or a Zarr store folder. A leading ``~`` in a path is
                expanded to the user's home directory.
            port: Port for the data server (auto-selected if None).
            height: Height of the embedded viewer in pixels.
            auto_open: Automatically display when created.
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
                ``uns['cellucid']['name']``; required when the data has none.
            dataset_id: Portable 1-180 byte ASCII identifier beginning with a
                letter or digit and otherwise using letters, digits, ``.``,
                ``_``, or ``-``; no trailing ``.`` or Windows device name.
                Defaults to ``uns['cellucid']['id']``; required when the data
                has none.
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
                None follows the viewer's default rule (``X`` when present).
            vector_field_default: The vector field the overlay opens on:
                ``<field>`` of an obsm array ``<field>_<basis>`` drawn over
                ``X_<basis>`` (``velocity`` for ``velocity_umap``).
            client_server_url: Exact browser-reachable data server base URL.
            web_source_url: Origin publishing the web asset inventory.
            web_cache_dir: Directory holding the active verified web build.
            allowed_hosts: Extra ``Host`` names the data server answers to.
                Required when ``client_server_url`` points at a reverse proxy
                such as ``jupyter-server-proxy``: name that proxy's host here.
        """
        if dataset_name is not None:
            dataset_name = _require_dataset_name(dataset_name, label="dataset_name")
        if dataset_id is not None:
            dataset_id = _require_dataset_id(dataset_id, label="dataset_id")
        self._init_state(
            port=port,
            height=height,
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
            )
            self._activate()

            if auto_open and self._context["in_jupyter"]:
                self.display()
        except BaseException:
            self._rollback_failed_construction()
            raise

    def _init_state(
        self,
        *,
        port: int | None,
        height: int,
        client_server_url: str | None,
        web_source_url: str,
        web_cache_dir: str | Path | None,
        allowed_hosts: Sequence[str] | None,
    ) -> None:
        """Set up everything but the data server: ports, ids, hooks, events."""
        self.port = 0 if port is None else require_server_port(port)
        self.allowed_hosts = require_allowed_hosts(allowed_hosts)
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

        self._server: AnnDataServer | None = None
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
        )
        self._server.start_background()
        self.port = self._server.port
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

        Like every hook, this runs on the server request thread, so raising
        here answers the viewer with HTTP 500 and loses the event rather than
        failing your cell. Record the failure and raise from your own code.
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
        event, data = _require_inbound_jupyter_event(
            message,
            expected_viewer_id=self._viewer_id,
        )
        self._hooks.trigger(event, data)
        self._record_event(event, data)
        if event == "command_error" and not self._hooks.has_callbacks("command_error"):
            self._report_unhandled_command_error(data)

    def _report_unhandled_command_error(self, data: dict[str, Any]) -> None:
        """
        Print one failed command to ``sys.stderr`` when nothing is listening.

        Commands are one-way: ``viewer.set_color_by(...)`` returns ``None``
        whether the viewer applied it or rejected it, and the rejection arrives
        on the server thread long after the cell finished. Nothing can make the
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

    def _await(self, waiter: _EventWaiter, timeout: float | None, message: str) -> None:
        """Block until ``waiter`` is handed its event, or raise TimeoutError."""
        deadline = None if timeout is None else (time.monotonic() + float(timeout))
        with self._event_cv:
            while waiter.payload is None:
                if deadline is not None and time.monotonic() >= deadline:
                    raise TimeoutError(message)
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
        from ..session_bundle import CellucidSessionBundle

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

    def set_data_mapping(self, *, timeout: float | None = 30.0, **fields: Any) -> dict[str, Any]:
        """
        Change what the viewer plots and reads, and wait until it applied it.

        Each keyword replaces one top-level field of the mapping the viewer
        last reported (:attr:`data_mapping`): ``embeddings``, ``vectorFields``,
        ``geneNames``, ``knn``, ``expression``, ``log1p``, ``highlyVariable``
        and ``latent``. ``kind`` cannot change. The merged mapping is validated
        before anything is sent, so an unknown field or an invalid value raises
        here and the viewer is untouched.

        The viewer answers with a ``data_mapping`` event of origin
        ``"notebook"`` naming this request; that event's mapping is returned and
        becomes :attr:`data_mapping`. If the viewer refuses the mapping -- an
        array that fails its load-time check, for example -- it reports a
        ``command_error`` and this raises ``RuntimeError`` with its reason.

        Example::

            viewer.set_data_mapping(
                embeddings={
                    "1d": None,
                    "2d": {"key": "X_pca", "columns": [0, 1]},
                    "3d": None,
                },
                knn="connectivities",
            )

        Until the viewer has reported a mapping (it does once the dataset
        loads), this raises ``RuntimeError`` before sending anything.
        """
        current = self.data_mapping
        if current is None:
            raise RuntimeError(
                "The viewer has not reported its data mapping yet. It reports one "
                "after the dataset loads; wait for it with "
                "viewer.wait_for_event('data_mapping')."
            )
        if "kind" in fields:
            raise TypeError("set_data_mapping cannot change the mapping kind.")
        allowed = set(current) - {"kind"}
        unknown = sorted(set(fields) - allowed)
        if unknown:
            raise TypeError(
                f"set_data_mapping got unknown field(s) {unknown}; a data "
                f"mapping has {sorted(allowed)}."
            )
        mapping = _require_data_mapping({**current, **fields}, label="data mapping")

        request_id = secrets.token_hex(16)
        waiter = _EventWaiter(
            lambda ev_type, payload: (
                ev_type == "data_mapping" and payload["requestId"] == request_id
            )
            or (ev_type == "command_error" and payload["command"] == "setDataMapping")
        )
        with self._event_cv:
            self._waiters.append(waiter)
        try:
            self.send_message(
                {
                    "type": "setDataMapping",
                    "requestId": request_id,
                    "mapping": mapping,
                }
            )
            self._await(
                waiter, timeout, "Timed out waiting for the viewer to apply the data mapping"
            )
        finally:
            with self._event_cv:
                self._waiters.remove(waiter)
        assert waiter.payload is not None
        if waiter.event == "command_error":
            raise RuntimeError(
                "The viewer did not apply the data mapping: " + waiter.payload["reason"]
            )
        return copy.deepcopy(waiter.payload["mapping"])

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

        ``timeout`` bounds the wait for the bundle (see
        :meth:`get_session_bundle`), ``cleanup_bundle`` deletes the captured
        bundle file afterwards, and ``expected_dataset_id`` defaults to the id
        this viewer serves. Every other keyword is passed to
        :func:`~cellucid.apply_cellucid_session_to_anndata`.

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
                if server is None:
                    raise RuntimeError(
                        "The viewer has stopped; pass expected_dataset_id explicitly."
                    )
                resolved_dataset_id = server.dataset_id
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
        - Python→Frontend ping/pong roundtrip (postMessage + HTTP events)
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
    def viewer_origin(self) -> str:
        """Origin (scheme+host+port) for postMessage targetOrigin."""
        parsed = urlparse(self.viewer_url)
        return f"{parsed.scheme}://{parsed.netloc}"

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

    def display(self):
        """Display the viewer in a notebook cell."""
        if not self._context["in_jupyter"]:
            console_print(f"Not in Jupyter environment. Open manually: {self.viewer_url}")
            return

        from IPython.display import HTML, display

        from ..web_cache import verify_web_ui_cache

        verify_web_ui_cache(
            self.web_cache_dir,
            expected_source_url=self.web_source_url,
        )

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

            // Note: Frontend → Python communication uses HTTP POST to /_cellucid/events
            // (and /_cellucid/selection for selections).
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

        Example::

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

        Example::

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

        Example::

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

        Example::

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
        1. Stops the underlying data server (which removes a temporary AnnData file)
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

    def __repr__(self) -> str:
        status = "running" if self._server and self._server.is_running() else "stopped"
        return f"AnnDataViewer(port={self.port}, status={status})"

    def _repr_html_(self) -> str:
        """HTML representation for Jupyter display."""
        if self._context["in_jupyter"]:
            return f'<a href="{self.viewer_url}" target="_blank">Open Cellucid Viewer</a>'
        return repr(self)


def show_anndata(
    data: str | Path | anndata.AnnData,
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
    client_server_url: str | None = None,
    web_source_url: str = CELLUCID_WEB_URL,
    web_cache_dir: str | Path | None = None,
    allowed_hosts: Sequence[str] | None = None,
) -> AnnDataViewer:
    """
    Display an AnnData object, .h5ad file, or Zarr store in a notebook.

    The viewer reads the store in place with range requests, without running
    prepare first. An in-memory object is written once to a temporary .h5ad.

    Args:
        data: AnnData object, or path to an .h5ad file, a zipped Zarr store or
            a Zarr store folder. A leading ``~`` in a path is expanded to the
            user's home directory.
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
            ``uns['cellucid']['name']``; required when the data has none.
        dataset_id: Portable 1-180 byte ASCII identifier beginning with a
            letter or digit and otherwise using letters, digits, ``.``, ``_``,
            or ``-``; no trailing ``.`` or Windows device name. Defaults to
            ``uns['cellucid']['id']``; required when the data has none.
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
            follows the viewer's default rule (``X`` when present).
        vector_field_default: The vector field the overlay opens on:
            ``<field>`` of an obsm array ``<field>_<basis>`` drawn over
            ``X_<basis>`` (``velocity`` for ``velocity_umap``).
        client_server_url: Exact browser-reachable data server base URL.
        web_source_url: Origin publishing the web asset inventory.
        web_cache_dir: Directory holding the active verified web build.
        allowed_hosts: Extra ``Host`` names the data server answers to.
            Required when ``client_server_url`` points at a reverse proxy such
            as ``jupyter-server-proxy``: name that proxy's host here.

    Returns:
        AnnDataViewer instance for interaction via hooks.

    Supported formats:
        - .h5ad files, zipped Zarr stores (.zip) and Zarr store folders, read
          in place
        - In-memory AnnData objects, written once to a temporary .h5ad

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
        To publish a dataset for others, prepare it (or publish the .h5ad or
        Zarr store itself) and share a link; see the data loading guide.
    """
    return AnnDataViewer(
        data=data,
        height=height,
        auto_open=True,
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
        client_server_url=client_server_url,
        web_source_url=web_source_url,
        web_cache_dir=web_cache_dir,
        allowed_hosts=allowed_hosts,
    )


# =============================================================================
# MESSAGE ROUTING (Frontend → Python)
# =============================================================================
#
# Bidirectional communication works in ALL environments via HTTP POST:
#
# 1. The viewer iframe POSTs events to the local data server
# 2. The data server (running on localhost) receives the POST at /_cellucid/events
#    (JSON events) or /_cellucid/selection (a selection, as binary cell indices)
# 3. The event is routed to the appropriate viewer via _handle_frontend_message
# 4. Hooks fire (on_selection, on_hover, on_click, etc.)
#
# This works because:
# - The data server has CORS headers allowing requests from cellucid.com
# - HTTP POST works universally (no special Jupyter/Colab APIs needed)
# - The viewer already knows the server URL (passed in URL params)
#
# =============================================================================


def _register_viewer_for_messages(viewer: AnnDataViewer):
    """
    Register a viewer to receive messages from frontend.

    This sets up the HTTP-based event routing which works in ALL environments.
    The frontend POSTs events to /_cellucid/events (and selections, as binary
    cell indices, to /_cellucid/selection) on the data server, which routes
    them to the appropriate viewer.
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
