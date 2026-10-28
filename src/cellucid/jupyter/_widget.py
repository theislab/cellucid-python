"""
The notebook widget: :class:`Widget`, what :func:`show_anndata` returns.

A :class:`Widget` is an `anywidget <https://anywidget.dev>`_ whose output
frames the Cellucid app (``<web_source_url>/?viewer=embedded``) and connects it to
this kernel through the notebook's own connection, the widget's comm. No
port is opened and no address of this machine is involved, so the same cell
works in JupyterLab, Jupyter Notebook 7, VS Code, Google Colab, marimo and
Databricks, locally or on a remote server.

Three hops carry everything, each authenticated by what carries it:

- **page ⇄ widget**: the widget's front end (``static/widget.js``) accepts
  the page's ``ready`` only from its own frame and the app's origin, and
  answers with a ``MessagePort`` and the viewer's id and token, posted to that
  exact origin. The token never appears in an address, a widget trait (which
  a notebook may save) or a log.
- **widget ⇄ kernel**: anywidget messages, with binary buffers. Every message
  names the output (``v``) it belongs to, since one widget can be shown in
  several outputs.
- **in the kernel**: the frames of the viewer channel go through
  :func:`cellucid._viewer_channel.serve_viewer_channel` -- the same
  ``hello``/``welcome``, commands, events, replacement rule and session bundle
  handling as the WebSocket of a browser-tab viewer -- over
  :class:`_CommSocket`; the dataset's requests (its documents and the
  store's bytes, by ``Range``) are answered by
  :meth:`cellucid.anndata_server.AnnDataServer.answer_request`, which also
  backs the HTTP server.

The widget's messages are handled where the kernel handles every message
from the notebook: between cells. While a cell runs, the viewer waits for
bytes it has not read yet, and a call in that cell cannot wait for the
viewer's answer (:meth:`Widget._require_answerable` says so at once instead of
hanging); hooks and ``viewer.state`` keep up whenever the kernel is idle.
"""

from __future__ import annotations

import logging
import queue
import sys
import threading
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any
from urllib.parse import urlsplit

import anywidget
import traitlets

from .._contracts import _require_dataset_id, _require_dataset_name
from .._server_base import CELLUCID_WEB_URL
from .._viewer_channel import MAX_EVENT_TEXT_BYTES, serve_viewer_channel
from .._websocket import (
    CLOSE_GOING_AWAY,
    CLOSE_POLICY_VIOLATION,
    WebSocketClosed,
    WebSocketProtocolError,
)
from ..data_mapping import EmbeddingKeySpec
from ._viewer import AnnDataViewer

if TYPE_CHECKING:
    import anndata
    import mudata

logger = logging.getLogger("cellucid.jupyter")

_STATIC = Path(__file__).with_name("static")

#: Bodies sent to the front end in one message, at most (a message's shape,
#: not a limit on what is read: a larger body is sent in a message of its own).
_ANSWER_MESSAGE_BYTES = 8 * 1024 * 1024

# A WebSocket close code for a text frame larger than the channel accepts.
_CLOSE_TOO_BIG = 1009


def _require_web_source_url(value: object) -> str:
    """The exact origin (and optional path) the browser loads the app from."""
    if not isinstance(value, str):
        raise TypeError("web_source_url must be a string")
    parts = urlsplit(value)
    if (
        parts.scheme not in ("https", "http")
        or not parts.hostname
        or parts.username is not None
        or parts.password is not None
        or parts.query
        or parts.fragment
        or value != value.strip()
        or value.endswith("/")
    ):
        raise ValueError(
            "web_source_url must be an exact http(s) address of the Cellucid app without a query, "
            "fragment, credentials or trailing slash, such as https://www.cellucid.com"
        )
    return value


# =============================================================================
# The viewer channel over the comm
# =============================================================================


class _CommMessage:
    """One binary frame from the page, read as its pieces arrive (``ChannelMessage``).

    The widget's front end sends a frame in pieces, so no single message is
    larger than the notebook's connection carries; the session bundle a
    frame may hold is streamed to disk as it arrives, as over a WebSocket.
    """

    __slots__ = ("_last", "_offset", "_piece", "_socket", "done")

    def __init__(self, socket_: _CommSocket, piece: bytes, last: bool) -> None:
        self._socket = socket_
        self._piece = memoryview(piece)
        self._offset = 0
        self._last = last
        self.done = last and not piece

    def readinto(self, buffer: memoryview) -> int:
        written = 0
        while written < len(buffer) and not self.done:
            if self._offset == len(self._piece):
                if self._last:
                    self.done = True
                    break
                kind, piece, last = self._socket._next_piece()
                self._piece = memoryview(piece)
                self._offset = 0
                self._last = last
                continue
            count = min(len(buffer) - written, len(self._piece) - self._offset)
            buffer[written : written + count] = self._piece[self._offset : self._offset + count]
            self._offset += count
            written += count
            if self._offset == len(self._piece) and self._last:
                self.done = True
        return written

    def read_exact(self, size: int) -> bytearray:
        result = bytearray(size)
        filled = 0
        with memoryview(result) as view:
            while filled < size:
                count = self.readinto(view[filled:])
                if count == 0:
                    raise WebSocketProtocolError(CLOSE_POLICY_VIOLATION, "Binary message ended early")
                filled += count
        return result

    def read_rest(self, initial: int = 0) -> bytearray:
        result = bytearray(max(initial, len(self._piece) - self._offset))
        filled = 0
        while True:
            if filled == len(result):
                if self.done:
                    break
                result.extend(bytes(max(len(result), 64 * 1024)))
            with memoryview(result) as whole, whole[filled:] as view:
                count = self.readinto(view)
            if count == 0:
                break
            filled += count
        del result[filled:]
        return result

    def drain(self) -> None:
        scratch = memoryview(bytearray(256 * 1024))
        while self.readinto(scratch):
            pass


class _CommSocket:
    """One viewer channel socket, carried by the widget's comm (``ChannelSocket``).

    The comm handler (the kernel's thread) feeds it the page's frames; the
    channel's own thread reads them with :meth:`receive`, exactly as it reads
    a WebSocket. Sends go out as widget messages from whichever thread sends.
    """

    def __init__(self, send: Callable[[dict[str, Any], list[bytes] | None], None], view: str, socket_id: int) -> None:
        self._send = send
        self._view = view
        self._id = socket_id
        self._frames: queue.SimpleQueue[tuple[Any, ...]] = queue.SimpleQueue()
        self._closed = threading.Event()
        self.close_code: int | None = None

    @property
    def closed(self) -> bool:
        return self._closed.is_set()

    # Fed by the widget's comm handler.

    def feed_text(self, text: str) -> None:
        self._frames.put(("text", text))

    def feed_binary(self, piece: bytes, last: bool) -> None:
        self._frames.put(("binary", piece, last))

    def feed_close(self, code: int, reason: str) -> None:
        self._frames.put(("close", code, reason))

    # Read by the channel's thread.

    def _next(self) -> tuple[Any, ...]:
        if self._closed.is_set():
            raise WebSocketClosed(self.close_code or CLOSE_GOING_AWAY)
        item = self._frames.get()
        if item[0] == "close":
            # The page closed it (or its output went away): nothing is sent back.
            self.close_code = item[1]
            self._closed.set()
            raise WebSocketClosed(item[1], item[2])
        if item[0] == "closed":
            raise WebSocketClosed(self.close_code or CLOSE_GOING_AWAY)
        return item

    def _next_piece(self) -> tuple[str, bytes, bool]:
        item = self._next()
        if item[0] != "binary":
            raise WebSocketProtocolError(CLOSE_POLICY_VIOLATION, "A binary frame was interrupted")
        return item[0], item[1], item[2]

    def receive(self) -> str | _CommMessage:
        item = self._next()
        if item[0] == "text":
            text: str = item[1]
            if len(text.encode("utf-8")) > MAX_EVENT_TEXT_BYTES:
                raise WebSocketProtocolError(_CLOSE_TOO_BIG, "Text frame too large")
            return text
        return _CommMessage(self, item[1], item[2])

    def send_text(self, text: str) -> None:
        if self._closed.is_set():
            raise WebSocketClosed(self.close_code or CLOSE_GOING_AWAY)
        self._send({"t": "text", "v": self._view, "s": self._id, "text": text}, None)

    def send_binary(self, payload: bytes | bytearray) -> None:
        if self._closed.is_set():
            raise WebSocketClosed(self.close_code or CLOSE_GOING_AWAY)
        self._send({"t": "binary", "v": self._view, "s": self._id}, [bytes(payload)])

    def ping(self) -> None:
        """Nothing to keep open: the notebook's own connection carries the frames."""

    def close(self, code: int, reason: str = "") -> None:
        if self._closed.is_set():
            return
        self.close_code = code
        self._closed.set()
        # Wake the channel's thread, then tell the page.
        self._frames.put(("closed",))
        try:
            self._send({"t": "close", "v": self._view, "s": self._id, "code": code, "reason": reason}, None)
        except Exception:
            logger.debug("Could not tell the widget's page that its channel closed", exc_info=True)

    def finish(self) -> None:
        """Nothing to drain: the comm delivers whole messages."""


# =============================================================================
# The widget
# =============================================================================


def show_widget(widget: Widget) -> None:
    """Show ``widget`` in the cell's output where the notebook allows it.

    An IPython kernel (JupyterLab, Notebook 7, VS Code, Colab, Databricks)
    displays it now; marimo shows the widget the cell evaluates to, so there
    nothing is displayed here. Outside a notebook there is nowhere to show it.
    """
    if widget._context["in_jupyter"]:
        widget.display()
        return
    marimo = sys.modules.get("marimo")
    if marimo is not None and marimo.running_in_notebook():
        return
    raise RuntimeError(
        "show_anndata() shows the viewer in a notebook cell, and this is not a notebook. "
        "Open it in a browser tab instead: cellucid.AnnDataViewer(...), "
        "cellucid.serve_anndata(...) or `cellucid serve`."
    )


class Widget(AnnDataViewer, anywidget.AnyWidget):
    """
    The Cellucid viewer as a notebook widget.

    The same viewer API as :class:`AnnDataViewer` -- hooks
    (``@viewer.on_selection``, ``on_hover``, ``on_click``, ``on_ready``, ...),
    ``viewer.state``, commands (``highlight_cells``, ``set_color_by``, ...),
    session bundles -- shown in a notebook cell's output and reached through
    the notebook's own connection: no port, proxy or server address, so it
    works in JupyterLab, Jupyter Notebook 7, VS Code, Google Colab, marimo
    and Databricks alike. :func:`show_anndata` creates one and shows it;
    a ``Widget`` made directly is shown by making it a cell's output.

    The page is the Cellucid app at ``web_source_url`` (https://www.cellucid.com by
    default), loaded by the browser; the dataset never leaves this kernel
    except to that page, through the widget.

    The viewer's messages reach Python between cells, as every notebook
    widget's do. So, in the cell that asks:

    - ``wait_for_ready()`` and ``wait_for_event()`` cannot wait for the viewer
      and raise at once with what to do instead;
      ``request_session_bundle()`` asks now, and its ``result()`` (or
      ``apply_session_to_anndata(adata, request)``) collects the bundle in a
      later cell;
    - hooks run, and ``viewer.state`` updates, whenever the kernel is idle;
    - commands (``highlight_cells``, ``push``, ``set_data_mapping`` and the
      rest) wait for nothing: they go as soon as the viewer can take them.

    Example::

        viewer = cellucid.Widget(adata)
        viewer            # the cell's output shows the viewer

        @viewer.on_selection
        def handle(event):
            print(len(event["cells"]))
    """

    _esm = _STATIC / "widget.js"
    _css = _STATIC / "widget.css"

    #: Height of the viewer in pixels.
    height = traitlets.Int(600).tag(sync=True)
    #: Where the browser loads the Cellucid app from.
    web_source_url = traitlets.Unicode(CELLUCID_WEB_URL).tag(sync=True)

    def __init__(
        self,
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
    ) -> None:
        """
        Prepare the dataset for the viewer; the viewer starts when the widget is shown.

        Args:
            data: AnnData or MuData object, or path to a file or folder
                Cellucid reads (:func:`show_anndata` lists them).
            height: Height of the viewer in pixels.
            web_source_url: Where the browser loads the Cellucid app:
                ``https://www.cellucid.com``, or the address of a copy you
                host (no query, fragment or trailing slash).

        Every other argument is :func:`show_anndata`'s.
        """
        if type(height) is not int or height <= 0:
            raise ValueError("height must be a positive whole number of pixels")
        if dataset_name is not None:
            dataset_name = _require_dataset_name(dataset_name, label="dataset_name")
        if dataset_id is not None:
            dataset_id = _require_dataset_id(dataset_id, label="dataset_id")
        anywidget.AnyWidget.__init__(self, height=height, web_source_url=_require_web_source_url(web_source_url))
        self._server = None
        self._init_viewer_state()
        #: The page's channel sockets, by (output, socket id).
        self._sockets: dict[tuple[str, int], _CommSocket] = {}
        self._sockets_lock = threading.Lock()
        self._send_lock = threading.Lock()
        #: The thread the kernel hands this widget's messages to: a wait on it
        #: could never be answered (:meth:`_require_answerable`).
        self._comm_thread: threading.Thread | None = None
        #: The execution count of the cell that displayed it (:meth:`display`).
        self._displayed_in_cell: int | None = None
        try:
            self.data = data
            self._prepare_dataset(
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
            )
            self._activate()
            self.on_msg(self._on_page_message)
        except BaseException:
            try:
                self._rollback_failed_construction()
            finally:
                self.close()
            raise

    def _prepare_dataset(self, **choices: Any) -> None:
        """Check the dataset and its launch choices; nothing is bound or served."""
        from ..anndata_server import AnnDataServer

        server = AnnDataServer(
            data=self.data,
            quiet=True,
            serve_web_ui=False,
            **choices,
        )
        # The viewer names the program serving it from _cellucid/info; a
        # widget has no host or port of its own.
        server.runtime = "jupyter"
        server.server_info.update(runtime="jupyter", host=None, port=None)
        self._server = server
        self._served_dataset_id = server.dataset_id
        cell_count = server.server_info["n_cells"]
        self._served_cell_count = cell_count if isinstance(cell_count, int) else None

    # ------------------------------------------------------------------
    # Messages from the widget's front end
    # ------------------------------------------------------------------

    def _send_to_page(self, content: dict[str, Any], buffers: list[bytes] | None) -> None:
        with self._send_lock:
            self.send(content, buffers=buffers)

    def _on_page_message(self, _widget: Any, content: Any, buffers: Sequence[Any]) -> None:
        self._comm_thread = threading.current_thread()
        if not isinstance(content, dict) or not isinstance(content.get("v"), str):
            logger.warning("Ignored a widget message without its output id")
            return
        view = content["v"]
        kind = content.get("t")
        try:
            if kind == "attach":
                # An output shows the widget (however it was displayed).
                self._displayed = True
                self._send_to_page(
                    {
                        "t": "credentials",
                        "v": view,
                        "viewerId": self._viewer_id,
                        "viewerToken": self._viewer_token,
                    },
                    None,
                )
            elif kind == "fetch":
                self._answer(view, content.get("requests"))
            elif kind == "open":
                self._open_socket(view, content.get("s"))
            elif kind in ("text", "binary", "close"):
                socket_ = self._socket(view, content.get("s"))
                if socket_ is None:
                    return
                if kind == "text" and isinstance(content.get("text"), str):
                    socket_.feed_text(content["text"])
                elif kind == "binary" and len(buffers) == 1:
                    socket_.feed_binary(bytes(buffers[0]), content.get("last") is True)
                elif kind == "close":
                    code = content.get("code")
                    socket_.feed_close(code if type(code) is int else CLOSE_GOING_AWAY, "")
                    self._forget_socket(view, content.get("s"))
                else:
                    socket_.feed_close(CLOSE_POLICY_VIOLATION, "malformed frame")
            elif kind == "detach":
                # The output was removed (cleared, re-rendered, the page left).
                for key in [key for key in list(self._sockets) if key[0] == view]:
                    socket_ = self._sockets.pop(key, None)
                    if socket_ is not None:
                        socket_.feed_close(CLOSE_GOING_AWAY, "")
            else:
                logger.warning("Ignored a widget message of unknown kind %r", kind)
        except Exception:
            logger.exception("A widget message could not be handled")

    def _socket(self, view: str, socket_id: object) -> _CommSocket | None:
        if type(socket_id) is not int:
            return None
        with self._sockets_lock:
            return self._sockets.get((view, socket_id))

    def _forget_socket(self, view: str, socket_id: object) -> None:
        with self._sockets_lock:
            self._sockets.pop((view, socket_id), None)  # type: ignore[arg-type]

    def _open_socket(self, view: str, socket_id: object) -> None:
        if type(socket_id) is not int:
            return
        socket_ = _CommSocket(self._send_to_page, view, socket_id)
        with self._sockets_lock:
            previous = self._sockets.pop((view, socket_id), None)
            self._sockets[(view, socket_id)] = socket_
        if previous is not None:
            previous.feed_close(CLOSE_GOING_AWAY, "")
        # The channel's own thread, as a WebSocket has: hooks run there.
        threading.Thread(
            target=serve_viewer_channel,
            args=(socket_,),
            name="cellucid-widget-channel",
            daemon=True,
        ).start()

    def _answer(self, view: str, requests: object) -> None:
        """Answer the page's requests, several bodies to a message."""
        if not isinstance(requests, list):
            return
        server = self._server
        answers: list[dict[str, Any]] = []
        bodies: list[bytes] = []
        size = 0

        def flush() -> None:
            nonlocal answers, bodies, size
            if answers:
                self._send_to_page({"t": "answers", "v": view, "answers": answers}, bodies)
            answers, bodies, size = [], [], 0

        for request in requests:
            if not isinstance(request, dict) or type(request.get("id")) is not int:
                continue
            method, path, range_header = request.get("method"), request.get("path"), request.get("range")
            if server is None:
                status, headers, body = 410, {"Content-Type": "text/plain"}, b"The viewer was stopped"
            elif (
                not isinstance(method, str)
                or not isinstance(path, str)
                or not (range_header is None or isinstance(range_header, str))
            ):
                status, headers, body = 400, {"Content-Type": "text/plain"}, b"Malformed request"
            else:
                try:
                    status, headers, body = server.answer_request(method, path, range_header)
                except Exception as error:
                    logger.exception("A widget request could not be answered")
                    status, headers, body = 500, {"Content-Type": "text/plain"}, str(error).encode()
            if body and size and size + len(body) > _ANSWER_MESSAGE_BYTES:
                flush()
            answer: dict[str, Any] = {"id": request["id"], "status": status, "headers": headers, "buffer": -1}
            if body:
                answer["buffer"] = len(bodies)
                bodies.append(body)
                size += len(body)
            answers.append(answer)
        flush()

    # ------------------------------------------------------------------
    # Waiting, shown, stopped
    # ------------------------------------------------------------------

    def _require_answerable(self, what: str) -> None:
        """Refuse a wait on the thread that hands this widget the viewer's answers.

        The kernel handles a widget's messages between cells, on the thread
        that runs them: a cell that waits there would wait until its timeout
        for an answer that cannot arrive before it ends.
        """
        if threading.current_thread() is not (self._comm_thread or threading.main_thread()):
            return
        if "session bundle" in what or "session_bundle" in what:
            hint = (
                "ask with request = viewer.request_session_bundle() and collect it in the next "
                "cell with request.result(), or apply it there with "
                "viewer.apply_session_to_anndata(adata, request)"
            )
        elif "results" in what:
            hint = "read viewer.results, which the viewer keeps current, in a later cell"
        else:
            hint = (
                "register a hook (@viewer.on_ready, @viewer.on_selection, ...) or read "
                "viewer.state in a later cell"
            )
        raise RuntimeError(
            f"{what} waits for the viewer, and a notebook widget hears the viewer only between "
            f"cells: nothing it says can arrive while this cell runs. Instead, {hint}. (A thread "
            "of your own may wait: the viewer answers it as soon as the cell ends.)"
        )

    def _await(self, waiter: Any, timeout: float | None, message: str) -> None:
        if waiter.payload is None:
            self._require_answerable(message.removeprefix("Timed out waiting for "))
        super()._await(waiter, timeout, message)

    def wait_for_ready(self, timeout: float | None = 30.0) -> dict[str, Any]:
        """Return the viewer's ``ready`` payload, waiting for it if it has not arrived.

        In the cell that showed the widget (or any cell before the viewer
        loaded) this cannot wait -- the viewer loads once the cell ends -- and
        raises at once; ``@viewer.on_ready`` runs when it is ready.
        """
        if self.state.ready is not None:
            return self.state.ready
        self._require_answerable("wait_for_ready()")
        return self.wait_for_event("ready", timeout=timeout)

    def display(self) -> None:
        """Show the viewer in the cell's output (an IPython kernel)."""
        from IPython.display import display

        display(self)
        self._displayed = True
        self._displayed_in_cell = _execution_count()

    def _shown_by_this_cell(self) -> bool:
        count = _execution_count()
        return count is not None and count == self._displayed_in_cell

    def _repr_mimebundle_(self, **kwargs: Any) -> Any:
        # The cell that showed it (show_anndata(...) as its last line) does
        # not show it a second time.
        if self._shown_by_this_cell():
            return {"text/plain": repr(self)}
        self._displayed = True
        return anywidget.AnyWidget._repr_mimebundle_(self, **kwargs)

    def _repr_html_(self) -> str | None:  # type: ignore[override]
        """What a notebook that cannot show widgets shows instead."""
        if self._shown_by_this_cell():
            return None
        return (
            "<p>This output is a Cellucid viewer (a notebook widget). This notebook "
            "cannot show widgets here: open it in JupyterLab, Jupyter Notebook 7, VS Code, "
            "Google Colab, marimo or Databricks, with <code>anywidget</code> installed where "
            "the notebook server runs.</p>"
        )

    @property
    def viewer_url(self) -> str:
        """A widget has no address of its own: its page lives in the cell's output."""
        raise AttributeError(
            "A notebook widget has no viewer address: it is shown in the cell's output. "
            "For a viewer in a browser tab, use cellucid.AnnDataViewer(...)."
        )

    @property
    def server_url(self) -> str:
        """A widget runs no server: the page reads the data through the widget."""
        raise AttributeError(
            "A notebook widget runs no server: the page reads the data through the widget."
        )

    def clear_web_cache(self) -> Path:
        """A widget keeps no web build: the browser loads the app from ``web_source_url``."""
        raise AttributeError(
            "A notebook widget keeps no web build: the browser loads the app from web_source_url."
        )

    def ensure_web_ui_cached(self, *, force: bool = True, show_progress: bool = True) -> Any:
        """A widget keeps no web build: the browser loads the app from ``web_source_url``."""
        raise AttributeError(
            "A notebook widget keeps no web build: the browser loads the app from web_source_url."
        )

    def debug_connection(self, timeout: float | None = 5.0) -> dict[str, Any]:
        """
        Return a structured report of this widget's connection.

        Says where the app is loaded from, which outputs are connected and
        whether the channel is open, what the viewer last reported, and the
        commands it refused. It waits for nothing, so it works in any cell.
        """
        with self._sockets_lock:
            outputs = sorted({view for view, _ in self._sockets})
        with self._event_cv:
            by_type: dict[str, int] = {}
            command_errors: list[dict[str, Any]] = []
            for _seq, ev_type, payload in self._recent_events:
                by_type[ev_type] = by_type.get(ev_type, 0) + 1
                if ev_type == "command_error":
                    command_errors.append(payload)
        from .._viewer_channel import viewer_is_connected

        return {
            "viewer_id": self._viewer_id,
            "transport": "notebook widget",
            "web_source_url": self.web_source_url,
            "displayed": self._displayed,
            "notebook_context": dict(self._context),
            "connected_outputs": outputs,
            "channel_connected": (
                self._message_routing_registered
                and viewer_is_connected(self._viewer_id, self._viewer_token)
            ),
            "dataset": None
            if self._server is None
            else {"id": self._server.dataset_id, "format": self._server.data_format},
            "state": {
                "ready": self.state.ready,
                "last_event_type": self.state.last_event_type,
                "last_updated_at": self.state.last_updated_at,
            },
            "recent_events": {
                "count": len(self._recent_events),
                "by_type": by_type,
                "command_errors": command_errors[-10:],
            },
        }

    def stop(self) -> None:
        """Stop the viewer: its page keeps the last view, and the temporary file is removed.

        The output stays on screen, frozen; Python stops answering it.
        """
        try:
            super().stop()
        finally:
            with self._sockets_lock:
                sockets = list(self._sockets.values())
                self._sockets.clear()
            for socket_ in sockets:
                socket_.close(CLOSE_GOING_AWAY, "The viewer was stopped")

    def __repr__(self) -> str:
        status = "running" if self._message_routing_registered else "stopped"
        dataset = None if self._server is None else self._server.dataset_id
        return f"cellucid.Widget(dataset={dataset!r}, status={status})"


def _execution_count() -> int | None:
    try:
        from IPython import get_ipython
    except ImportError:
        return None
    shell = get_ipython()
    return None if shell is None else getattr(shell, "execution_count", None)


__all__ = ["Widget", "show_widget"]
