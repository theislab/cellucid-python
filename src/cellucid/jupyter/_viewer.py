"""
The viewer's Python side: :class:`AnnDataViewer` and :func:`show_anndata`.

:class:`AnnDataViewer` is the viewer API -- hooks, the event log that
``wait_for_event`` reads, ``viewer.state``, session capture, the actions (one
method per action of the catalogue, ``_actions_generated.py``, and
:meth:`~AnnDataViewer.act`), shutdown -- and, on its own, a viewer in a browser tab: it
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
import contextlib
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
from collections.abc import Callable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from .._agent_link import ViewerAgentPort, publish_viewer, withdraw_viewer
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
    field_push_params,
    read_obs_names,
)
from ..data_mapping import EmbeddingKeySpec
from ._actions import action_params, cell_positions, tool_name
from ._actions_generated import ViewerActions
from ._history import LiveHistory, LoadSessionRequest
from ._hooks import HookCallback, HookRegistry, ViewerState
from ._wire import (
    _AGENT_EVENT_TYPES,
    _require_client_server_url,
    _require_debug_dataset_list,
    _require_frontend_message,
    _require_inbound_jupyter_event,
)

if TYPE_CHECKING:
    import anndata
    import mudata

    from ..anndata_server import AnnDataServer
    from ..session import CellucidSession
    from ._widget import Widget

logger = logging.getLogger("cellucid.jupyter")

# The public method that sends each command, used to report a `command_error`
# in the notebook's own vocabulary: a viewer failure is far easier to act on
# when it names the call that caused it rather than the wire type. Every
# command of `_wire._COMMAND_TYPES` must appear here, and every value must name a
# real `AnnDataViewer` method; both are asserted by the contract tests. A
# failed `action` is named by its own method (`viewer.color_by(...)`).
_COMMAND_PYTHON_CALLERS: dict[str, str] = {
    "ping": "viewer.debug_connection()",
    "debug_snapshot": "viewer.debug_connection()",
    "requestSessionBundle": "viewer.request_session()",
    "freeze": "viewer.stop()",
    "action": "viewer.act(...)",
    "load_session": "viewer.load_session(...)",
    # An agent app's commands (`cellucid mcp`): their failures reach the
    # agent that sent them (`ViewerAgentPort.command_failed`), never a notice.
    "agent_list": "cellucid mcp",
    "agent_call": "cellucid mcp",
    "agent_cancel": "cellucid mcp",
    "agent_answer": "cellucid mcp",
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

class ActionError(RuntimeError):
    """The viewer did not apply an action: :attr:`reason` says why, in the
    viewer's words (a name that does not resolve names what does)."""

    def __init__(self, action: str, request_id: str, reason: str) -> None:
        super().__init__(f"The viewer did not apply {action}: {reason}")
        self.action = action
        self.request_id = request_id
        self.reason = reason


class ActionRequest:
    """One action sent to the viewer (:meth:`AnnDataViewer.act`, or the
    action's own method).

    The viewer answers it once: with its outcome -- a read's result -- and the
    history step it made, or with why it did not apply it. ``done`` says
    whether the answer arrived; :meth:`result` returns it. Nothing has to
    collect a request: an action applies whether or not anyone waits for it.
    """

    def __init__(self, viewer: AnnDataViewer, request_id: str, action: str, waiter: _EventWaiter) -> None:
        self._viewer = viewer
        self.request_id = request_id
        #: The action's dotted type (``color.by``).
        self.action = action
        self._waiter: _EventWaiter | None = waiter

    @property
    def done(self) -> bool:
        """Whether the viewer answered (with the outcome, or with why not)."""
        return self._waiter is not None and self._waiter.payload is not None

    @property
    def step(self) -> dict[str, Any] | None:
        """The history step the action made, ``{"id", "label"}``, once it was
        answered; None before, for a read, and for a change that changed
        nothing."""
        waiter = self._waiter
        if waiter is None or waiter.payload is None or waiter.event != "action_result":
            return None
        step = waiter.payload["step"]
        return None if step is None else dict(step)

    def result(self, timeout: float | None = 60.0) -> Any:
        """
        The action's outcome, waiting up to ``timeout`` seconds for it (None:
        forever): a read's result (JSON: dicts, lists, text and numbers), a
        change's outcome (counts and names), or None.

        Raises :class:`ActionError` with the viewer's reason when it did not
        apply the action, ``TimeoutError`` when no answer arrived in time,
        and ``RuntimeError`` when the viewer stopped or the request was
        cancelled -- or, in a notebook widget, when the answer has not arrived
        yet: the viewer answers between cells, so collect it in a later cell.
        """
        waiter = self._waiter
        if waiter is None:
            raise RuntimeError(f"This {self.action} request was cancelled.")
        if waiter.payload is None:
            self._viewer._await(waiter, timeout, f"Timed out waiting for the viewer's answer to {self.action}")
        assert waiter.payload is not None
        if waiter.event == "command_error":
            raise ActionError(self.action, self.request_id, waiter.payload["reason"])
        return copy.deepcopy(waiter.payload["outcome"])

    def cancel(self) -> None:
        """Stop waiting for the answer (the action itself is not withdrawn: the
        viewer may have applied it already); an answer that still arrives is
        dropped."""
        if self._waiter is None or self._waiter.payload is not None:
            return
        self._forget()
        self._waiter = None

    def _forget(self) -> None:
        viewer = self._viewer
        with viewer._event_cv:
            viewer._action_waiters.pop(self.request_id, None)
            viewer._action_types.pop(self.request_id, None)

    def __repr__(self) -> str:
        state = "cancelled" if self._waiter is None else "answered" if self.done else "waiting"
        return f"ActionRequest({self.action!r}, {state})"


class SessionRequest:
    """A session the viewer was asked for (:meth:`AnnDataViewer.request_session`).

    ``done`` says whether the viewer answered; :meth:`result` returns the
    session -- its history, every tracker, as a ``.cellucid-session`` file
    of format 4. A request that is never collected is withdrawn by
    :meth:`cancel`, or when the viewer stops.
    """

    def __init__(self, viewer: AnnDataViewer, request_id: str, waiter: _EventWaiter) -> None:
        self._viewer = viewer
        self.request_id = request_id
        self._waiter: _EventWaiter | None = waiter
        self._session: CellucidSession | None = None

    @property
    def done(self) -> bool:
        """Whether the viewer answered (with the session, or with why it could not)."""
        return self._session is not None or (self._waiter is not None and self._waiter.payload is not None)

    def result(self, timeout: float | None = 60.0) -> CellucidSession:
        """
        The session, waiting up to ``timeout`` seconds for it (None: forever).

        Returns the :class:`~cellucid.CellucidSession` (the same one on every
        call), read and checked whole. Raises ``TimeoutError`` when it did not
        arrive in time, and ``RuntimeError`` when the viewer could not capture
        it -- or, in a notebook widget, when it has not arrived yet: the
        viewer answers between cells, so collect it in a later cell.
        """
        from ..session import CellucidSession

        if self._session is not None:
            return self._session
        waiter = self._waiter
        if waiter is None:
            raise RuntimeError("This session request was cancelled.")
        viewer = self._viewer
        viewer._await(waiter, timeout, "Timed out waiting for the session")
        self._finish()
        assert waiter.payload is not None
        if waiter.event == "command_error":
            raise RuntimeError("The viewer did not capture the session: " + waiter.payload["reason"])
        path = waiter.payload.get("path")
        if not path:
            raise RuntimeError("The session arrived, but no file was kept for it")
        self._session = CellucidSession(Path(path))
        return self._session

    def cancel(self) -> None:
        """Withdraw the request; an answer that still arrives is dropped."""
        if self._waiter is None or self._session is not None:
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
        return f"SessionRequest(done={self.done})"


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


class AnnDataViewer(ViewerActions):
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

            @viewer.on_history_change
            def handle(event):
                print(event['kind'], event['head'])  # the history changed

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
        # before its cell count was known (`ready`), each with the request
        # already handed out for it.
        self._pending_pushes: list[tuple[Any, str, str, tuple[str, str] | None]] = []
        # The actions sent and not answered yet, by request id: what answers
        # one (`action_result`, or the `command_error` of its id) is handed to
        # its waiter; and each one's type, to name it in a failure notice.
        self._action_waiters: dict[str, _EventWaiter] = {}
        self._action_types: dict[str, str] = {}
        # The notebook turn open now (`turn()`): its id and label, carried by
        # every action sent while it is open; None outside one.
        self._turn: tuple[str, str] | None = None
        # This viewer's end of local agent connections (`cellucid mcp`): an
        # agent app's tool calls, carried to the page and answered back. It
        # reaches the viewer weakly, so a viewer let go is not kept alive by
        # the agents endpoint that publishes it.
        viewer_ref = weakref.ref(self)

        def agent_send(message: dict[str, Any]) -> object:
            viewer = viewer_ref()
            if viewer is None:
                raise LookupError("The viewer was let go")
            return viewer.send_message(message)

        def agent_connected() -> bool:
            viewer = viewer_ref()
            return viewer is not None and viewer_is_connected(viewer._viewer_id, viewer._viewer_token)

        def agent_describe() -> dict[str, Any]:
            viewer = viewer_ref()
            return {} if viewer is None else viewer._agent_description()

        self._agent_port = ViewerAgentPort(send=agent_send, connected=agent_connected, describe=agent_describe)

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
        """Publish this fully started viewer to the routing registries, and to
        agent apps on this computer (``cellucid mcp``)."""
        _register_viewer_for_messages(self)
        try:
            _active_viewers.add(self)
        except BaseException:
            _unregister_viewer_for_messages(self)
            raise
        publish_viewer(self._viewer_id, self._agent_port, runtime=self._agent_description()["runtime"])
        weakref.finalize(self, withdraw_viewer, self._viewer_id)

    def _agent_description(self) -> dict[str, Any]:
        """The viewer as an agent connection names it (``cellucid_status``):
        its id, the program serving it, its dataset, and the address it is
        served at (None for a notebook widget, shown in a cell)."""
        server = self._server
        # Where it is served, never its page address: that carries the
        # viewer's token, which an agent is not given.
        try:
            url: str | None = self._get_client_server_url()
        except (AttributeError, RuntimeError):
            url = None
        if type(self).viewer_url is not AnnDataViewer.viewer_url:
            url = None  # a notebook widget: shown in its cell, at no address
        return {
            "viewer": self._viewer_id,
            "runtime": "jupyter" if self._context["in_jupyter"] else "python",
            "dataset": None
            if server is None
            else {"id": server.dataset_id, "name": getattr(server, "dataset_name", None) or server.dataset_id},
            "url": url,
        }

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
        withdraw_viewer(self._viewer_id)
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

        Called when a command this viewer sent -- an action
        (``viewer.color_by(...)``, ``viewer.act(...)``), a session request --
        was refused by the viewer instead of carried out. An action returns its
        request at once, so its method has returned by the time the failure
        arrives: this hook (or the request's ``result()``, which raises the
        same reason) is how a notebook learns that the viewer is unchanged.

        Registering a handler also silences the notice Cellucid otherwise
        prints to ``sys.stderr``, so reporting happens once, where you chose.

        Event data:
            - command: str - the command type that failed: ``"action"``,
              ``"requestSessionBundle"``, ...
            - requestId: str or None - the request it fails (an action's
              ``request.request_id``); None for a command that names none
            - reason: str - what the viewer reported, in one line

        Like every hook, this runs on the viewer channel's thread, not in your
        cell: an exception raised here is logged (logger
        ``cellucid.viewer_channel``) and never reaches the notebook, the
        callbacks after it do not run, and the event is not recorded in
        ``viewer.state`` or handed to ``wait_for_event``. Record the failure
        and raise from your own code.
        An action that did apply is answered too: its request's ``result()``
        returns its outcome.

        Example::

            rejected = []

            @viewer.on_command_error
            def handle(event):
                rejected.append(event)

            viewer.color_by(field={"source": "obs", "key": "cell_type"})
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
        from the viewer's controls, a session, or this notebook's
        :meth:`mapping_set`. The viewer sends it because this package
        lists ``data_mapping`` among the events of ``GET /_cellucid/protocol``.

        Event data:
            - origin: str - why the mapping is what it is: ``"default"``,
              ``"best_guess"``, ``"session"``, ``"user"``, or
              ``"notebook"`` (this notebook's :meth:`mapping_set`)
            - mapping: dict - the complete data mapping, in canonical form

        Example::

            @viewer.on_data_mapping
            def handle(event):
                print(event["mapping"]["embeddings"])
        """
        return self._hooks.on("data_mapping")

    @property
    def on_history_change(self) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator to register a handler for every change of the viewer's history.

        Called when a step is recorded (by the person, an agent, or this
        notebook), the step shown moves (undo, redo, a jump, another branch),
        the history's structure is edited (a rename, a note, a bookmark, a
        delete, a restore), a tracker is started, used, renamed or removed, or
        a session is loaded. ``viewer.history`` keeps the latest.

        Event data:
            - kind: str - ``"step"``, ``"moved"``, ``"edited"``, ``"tracker"`` or ``"loaded"``
            - tracker: str - the id of the tracker in use
            - head: str - the step shown
            - step: str or None - the step recorded (``"step"``), else None
            - branch: str or None - the id of the active branch

        Example::

            @viewer.on_history_change
            def handle(event):
                if event["kind"] == "step":
                    print("recorded", event["step"])
        """
        return self._hooks.on("history_changed")

    @property
    def history(self) -> LiveHistory:
        """The viewer's history, live (:class:`~cellucid.jupyter.LiveHistory`):
        the latest change it told, a hook on each change, and its history read
        as a session. The history is changed through the generated
        ``history_*`` methods (``history_undo()``, ``history_jump(step=...)``,
        ``history_note_add(...)``, ...)."""
        return LiveHistory(self)

    @contextlib.contextmanager
    def turn(self, label: str) -> Iterator[str]:
        """
        Make every action sent inside the ``with`` block one turn of the viewer's history.

        The State panel shows the turn's steps as one group, labelled
        ``label`` and by this notebook, and undoes them together (Undo turn),
        as it does an agent's turn. Yields the turn's id. Turns do not nest.

        Example::

            with viewer.turn("Normalise and cluster"):
                viewer.color_by(field={"source": "obs", "key": "leiden"})
                viewer.highlight_cells(cells=adata.obs["leiden"] == "3", label="Leiden 3")
        """
        if not isinstance(label, str) or not label.strip():
            raise ValueError("A turn needs a label: what it does, in a few words.")
        turn_id = f"python-{secrets.token_hex(6)}"
        with self._state_lock:
            if self._turn is not None:
                raise RuntimeError(f'A turn ("{self._turn[1]}") is open already; turns do not nest.')
            self._turn = (turn_id, " ".join(label.split()))
        try:
            yield turn_id
        finally:
            with self._state_lock:
                self._turn = None

    def load_session(self, session: str | Path | CellucidSession) -> LoadSessionRequest:
        """
        Load a session file into the viewer: its trackers are added beside the
        person's own (nothing of theirs is replaced), and the file's tracker
        shown is shown.

        The session must be of the dataset this viewer shows. Returns the
        :class:`~cellucid.jupyter.LoadSessionRequest`; its result is
        ``{"tracker": {"id", "name"}, "head"}`` (in a notebook, collect it in
        a later cell).

        Example::

            viewer.load_session("analysis.cellucid-session")
        """
        from ..session import CellucidSession

        path = session.path if isinstance(session, CellucidSession) else Path(session).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(str(path))
        data = path.read_bytes()
        request_id = secrets.token_hex(16)
        waiter = _EventWaiter(lambda _event, _payload: False)
        with self._event_cv:
            self._action_waiters[request_id] = waiter
            self._action_types[request_id] = "load_session"
        request = LoadSessionRequest(self, request_id, path.name, waiter)
        try:
            self.send_message({"type": "load_session", "requestId": request_id, "name": path.name, "file": data})
        except BaseException:
            request.cancel()
            raise
        return request

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
        # An agent app's answers (`cellucid mcp`) are its own: they reach the
        # agent port, never the notebook's hooks or state.
        if event in _AGENT_EVENT_TYPES:
            self._agent_port.deliver(event, data)
            return
        if event == "command_error" and data["command"].startswith("agent_"):
            self._agent_port.command_failed(data)
            return
        # The event is recorded before any hook runs: a hook that raises is the
        # notebook's error (the channel logs it), and must not cost
        # viewer.state, wait_for_event or a session request the event.
        self._record_event(event, data)
        try:
            self._hooks.trigger(event, data)
        finally:
            if event == "command_error" and not self._hooks.has_callbacks("command_error"):
                self._report_unhandled_command_error(data)
            if event in {"action_result", "command_error", "session_loaded"} and data["requestId"] is not None:
                with self._state_lock:
                    self._action_types.pop(data["requestId"], None)
        # What waited for this event goes now.
        if event == "ready":
            self._send_pending_pushes()

    def _report_local_command_error(self, request_id: str, error: Exception) -> None:
        """An action that waited for the viewer and could not be built once it
        could: reported as the viewer reports a refused action (a
        ``command_error`` of its request: its request fails, and
        ``state.command_error``, the hooks, or printed)."""
        reason = str(error).strip() or type(error).__name__
        self._handle_frontend_message(
            {"type": "command_error", "command": "action", "requestId": request_id, "reason": reason[:4000]}
        )

    def _report_unhandled_command_error(self, data: dict[str, Any]) -> None:
        """
        Print one failed command to ``sys.stderr`` when nothing is listening.

        An action returns its request at once: ``viewer.color_by(...)`` does
        not wait for the viewer (a notebook widget hears it only between cells),
        and its refusal arrives on the viewer channel's thread after the cell
        finished. A notebook that registered no ``command_error`` hook and never
        collects the request would otherwise see a successful-looking call and
        an unchanged viewer; ``request.result()`` raises the same reason.

        This runs after the event has been recorded, so a failure to write does
        not cost the notebook the event.
        """
        command = data["command"]
        with self._state_lock:
            action = None if data["requestId"] is None else self._action_types.get(data["requestId"])
        caller = (
            f"viewer.{tool_name(action)}(...)" if command == "action" and action is not None
            else _COMMAND_PYTHON_CALLERS[command]
        )
        console_print(
            f"[cellucid] {caller} did not take effect "
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
            if event in {"action_result", "command_error", "session_loaded"} and data["requestId"] is not None:
                answered = self._action_waiters.pop(data["requestId"], None)
                if answered is not None:
                    answered.event = event
                    answered.payload = data

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
            elif event == "history_changed":
                self.state.history = data
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

    def request_session(self) -> SessionRequest:
        """
        Ask the viewer for its session -- its whole history -- without waiting.

        The request is sent at once (or as soon as the viewer connects) and
        the viewer answers in the background with a ``.cellucid-session``
        file of format 4: every tracker of the dataset's history, each step,
        and the values each step shows (conversations of the chat bar are
        left out). Collect it with the returned request's ``result()`` (a
        :class:`~cellucid.CellucidSession`). In a notebook widget the answer
        arrives once the cell that asked has finished, so collect it in a
        later cell::

            request = viewer.request_session()
            # next cell
            session = request.result()
            plan = session.apply(adata)  # what would change, one row each
            plan.apply()

        Returns
        -------
        SessionRequest
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
            or (ev_type == "command_error" and payload["requestId"] == request_id)
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
        return SessionRequest(self, request_id, waiter)

    @property
    def data_mapping(self) -> dict[str, Any] | None:
        """The data mapping the viewer last reported, or None before it reported one.

        This is a copy of the ``mapping`` of the latest ``data_mapping`` event:
        which embedding fills each dimension slot, the bound vector fields,
        and -- for an AnnData source -- the neighbour graph, the gene names,
        the expression matrix and the latent space. Editing the copy changes
        nothing; use :meth:`mapping_set` (the ``mapping.set`` action).
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

    # =========================================================================
    # ACTIONS: everything that changes or reads the viewer's state
    # =========================================================================

    def act(self, action: str, /, **params: Any) -> ActionRequest:
        """
        Apply (or read) one action of the viewer, by its dotted type.

        Every action has a method of its own as well (``viewer.color_by(...)``
        is ``viewer.act("color.by", ...)``); ``act`` is the same for an action
        named in data, and every action is listed by
        :func:`cellucid.jupyter.action_catalog`. Parameters are checked by the
        catalogue and converted the same way for every action: cells
        (``x-cells``) as rows, a boolean mask or obs names; one number per cell
        (``x-values``) as an array; everything else as JSON (NumPy values and
        pandas indexes as lists). A parameter written ``in`` in the catalogue
        is ``in_`` as a method's keyword.

        The action is sent at once (or as soon as the viewer connects) and
        the viewer applies it as a person's control would: one step of its
        history, labelled as this notebook's, undone like any other. Nothing
        waits: the returned request's :meth:`~ActionRequest.result` waits for
        the outcome -- a read's result -- and raises :class:`ActionError` with
        the viewer's reason when it did not apply the action. In a notebook
        widget the answer arrives once the cell has finished, so collect it
        in a later cell::

            request = viewer.fields_list()
            # next cell
            fields = request.result()

        Example::

            viewer.act("color.by", field={"source": "obs", "key": "cell_type"})
            viewer.highlight_cells(cells=adata.obs["leiden"] == "3", label="Leiden 3")
        """
        exact = action_params(action, params, cells=self._action_cells)
        request = self._expect_action(action)
        self._dispatch_action(request, exact)
        return request

    def _action_cells(self, value: Any, label: str) -> np.ndarray:
        """Cells an action names, as original positions (``uint32``)."""
        return cell_positions(
            value, label=label, n_obs=self._dataset_cell_count(), obs_names=lambda: self.obs_names
        )

    def _expect_action(self, action: str) -> ActionRequest:
        """A request for ``action``, its answer awaited from now on."""
        request_id = secrets.token_hex(16)
        waiter = _EventWaiter(lambda _event, _payload: False)
        with self._event_cv:
            self._action_waiters[request_id] = waiter
            self._action_types[request_id] = action
        return ActionRequest(self, request_id, action, waiter)

    def _dispatch_action(self, request: ActionRequest, params: dict[str, Any]) -> None:
        """Send ``request``'s action with ``params``; a send that fails withdraws it."""
        with self._state_lock:
            turn = self._turn
        try:
            self.send_message({
                "type": "action", "requestId": request.request_id, "action": request.action, "params": params,
                "turn": None if turn is None else turn[0], "turnLabel": None if turn is None else turn[1],
            })
        except BaseException:
            request._forget()
            raise

    # =========================================================================
    # ROUND TRIP: results back from the viewer, columns and gene sets pushed    # =========================================================================
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

    def push(self, values: Any, key: str | None = None) -> ActionRequest:
        """Show a column in the viewer as a field: labels as a categorical, numbers
        as a continuous field. Pushing the same key again replaces it.

        The ``field.push`` action, made from a Python column (its own method,
        ``viewer.field_push(...)``, takes the codes and values as the viewer
        stores them).

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
                replaced: the viewer refuses it, and says why.

        Returns:
            The :class:`ActionRequest`: its result is ``{"name", "kind",
            "id"}`` once the viewer added the field (:attr:`results` lists it
            from then on; in a notebook, from the next cell). The column is
            sent at once, or, for a format only the browser reads, once the
            viewer has said how many cells it has (a cell run before that,
            "Run All", is fine). A column the viewer refuses, or that does not
            fit the dataset once it is known, fails the request with the
            reason (printed unless a ``command_error`` hook handles it).
        """
        name = key if key is not None else getattr(values, "name", None)
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Name the field: pass key=, or a Series with a name")
        request = self._expect_action("field.push")
        n_obs = self._dataset_cell_count()
        if n_obs is None:
            # Only the browser reads this format: the column goes once the
            # viewer says how many cells it has (its `ready`).
            with self._state_lock:
                self._pending_pushes.append((values, name.strip(), request.request_id, self._turn))
            return request
        try:
            params = field_push_params(values, name.strip(), n_obs=n_obs, obs_names=self.obs_names)
        except BaseException:
            request._forget()
            raise
        self._dispatch_action(request, params)
        return request

    def _send_pending_pushes(self) -> None:
        """Send the columns pushed before the viewer's cell count was known; a
        column that does not fit the dataset fails its request."""
        with self._state_lock:
            pending, self._pending_pushes = self._pending_pushes, []
        n_obs = self._dataset_cell_count()
        if n_obs is None:
            return
        for values, name, request_id, turn in pending:
            with self._state_lock:
                waiter = self._action_waiters.get(request_id)
            if waiter is None:
                continue
            try:
                params = field_push_params(values, name, n_obs=n_obs, obs_names=self.obs_names)
            except (TypeError, ValueError) as error:
                self._report_local_command_error(request_id, error)
                continue
            self.send_message({
                "type": "action", "requestId": request_id, "action": "field.push", "params": params,
                "turn": None if turn is None else turn[0], "turnLabel": None if turn is None else turn[1],
            })

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
        request: SessionRequest,
        *,
        tracker: str | None = None,
        branch: str | None = None,
        step: str | None = None,
        dry_run: bool = True,
        overwrite: bool = False,
        keep_session: bool = False,
        timeout: float | None = 60.0,
        cleanup_file: bool = True,
    ) -> Any:
        """
        Plan or apply the session a :meth:`request_session` request brought to an AnnData.

        The one way to get the viewer's state into an object, in a notebook
        widget and a browser-tab viewer alike (cellucid-r's
        ``cellucid_apply_session(object, request)``): ask, then apply -- in a
        notebook, in a later cell, since the viewer answers between cells::

            request = viewer.request_session()
            # next cell
            plan = viewer.apply_session_to_anndata(adata, request)   # what would change
            report = plan.apply()

        ``tracker``, ``branch`` and ``step`` choose the state applied: by
        default the step the viewer shows (the active tracker's head); any
        step of any tracker is the state as it was after that step. Every
        other keyword but ``timeout`` (the wait for the session; a notebook
        widget's cell cannot wait, and says so at once when the answer has not
        arrived) and ``cleanup_file`` (delete the session file the viewer
        sent once it is read and planned: the plan, and a record kept with
        ``keep_session``, need it no longer; ``request.result().save(...)``
        first keeps a copy) is :func:`~cellucid.apply_cellucid_session_to_anndata`'s.
        """
        if type(cleanup_file) is not bool:
            raise TypeError("cleanup_file must be exactly True or False")
        if not isinstance(request, SessionRequest) or request._viewer is not self:
            raise TypeError("request must be what this viewer's request_session() returned")
        session = request.result(timeout=timeout)
        try:
            return session.apply(
                adata, tracker=tracker, branch=branch, step=step, dry_run=dry_run, overwrite=overwrite,
                keep_session=keep_session,
            )
        finally:
            if cleanup_file:
                session.path.unlink(missing_ok=True)

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

        An ``action`` whose parameters hold NumPy arrays (cells as ``uint32``
        original positions, a column's values) travels as one binary frame,
        the arrays little-endian; every other command as JSON text. Use
        :meth:`act` (or an action's own method) for actions: it checks the
        parameters by the catalogue and returns the request its answer
        completes.

        Args:
            message: Dict with the command's exact keys (the shared protocol's
                ``commands``).

        Returns:
            True when the command was sent now, False when it is waiting for
            the viewer to connect.

        Example::

            viewer.send_message({"type": "freeze"})
        """
        serialized_message = _require_frontend_message(message)
        if not self._message_routing_registered:
            raise RuntimeError("The viewer has stopped; it accepts no more commands.")
        return send_viewer_command(self._viewer_id, self._viewer_token, serialized_message)

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
        # Agent apps waiting on this viewer hear it stopped, and it is no
        # longer published to them.
        self._agent_port.stopped("the viewer was stopped")
        withdraw_viewer(self._viewer_id)
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

    def _socket(attached: bool) -> None:
        # A page went: the agent calls it was running get no answer now. A
        # page took the channel: what its tool list says may differ.
        resolved = viewer_ref()
        if resolved is None:
            return
        if attached:
            resolved._agent_port.page_back()
        else:
            resolved._agent_port.page_gone()

    register_event_callback(
        viewer._viewer_id,
        viewer._viewer_token,
        _deliver,
        on_socket=_socket,
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
