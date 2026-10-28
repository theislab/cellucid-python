"""The notebook widget (``cellucid.jupyter._widget``) against a mock comm.

The widget's front end is played by this test: it sends the widget the
messages ``static/widget.js`` relays from the page (``attach``, ``fetch``,
``open``, ``text``, ``binary`` in pieces, ``close``, ``detach``) and reads what
the widget sends back. Everything the kernel does with them is real: the
dataset's documents and bytes come from
:meth:`cellucid.anndata_server.AnnDataServer.answer_request`, and the frames
go through :func:`cellucid._viewer_channel.serve_viewer_channel`, the channel
a browser-tab viewer's WebSocket runs.

Events are built from the shared protocol fixture's examples, so this file
follows the protocol as it is.
"""

from __future__ import annotations

import json
import struct
import threading
import time
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid.jupyter import Widget
from cellucid.jupyter._widget import _ANSWER_MESSAGE_BYTES, show_widget
from cellucid.prepare_data._priority import ORIGINAL_POSITION_COLUMN, CellPriority, prefix_levels
from cellucid.session_bundle import SESSION_BUNDLE_MAGIC

PROTOCOL = json.loads((Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text())


def _example_event(kind: str) -> dict[str, Any]:
    return next(event for event in PROTOCOL["examples"]["events"] if event["type"] == kind)


def _binary_frame(header: dict[str, Any], body: bytes) -> bytes:
    encoded = json.dumps(header, separators=(",", ":")).encode()
    return struct.pack("<I", len(encoded)) + encoded + body


def _bundle() -> bytes:
    """The smallest valid session bundle: the magic and a manifest with no chunks."""
    manifest = json.dumps(
        {
            "createdAt": "2026-10-27T00:00:00.000Z",
            "datasetFingerprint": {
                "datasetId": None,
                "name": None,
                "cellCount": 0,
                "varCount": 0,
                "cellOrder": {
                    "embedding": {"key": "X_umap", "columns": [0, 1]},
                    "dimension": 2,
                    "digest": "0123456789abcdef",
                },
                "dataMapping": {
                    "kind": "anndata",
                    "embeddings": {
                        "enabled": [{"key": "X_umap", "columns": [0, 1]}],
                        "defaults": {"1d": None, "2d": {"key": "X_umap", "columns": [0, 1]}, "3d": None},
                    },
                    "knn": None,
                    "geneNames": "index",
                    "expression": "X",
                    "log1p": {},
                    "highlyVariable": None,
                    "latent": None,
                    "vectorFields": [],
                },
            },
            "source": None,
            "requires": {"obs": [], "var": []},
            "chunks": [],
        },
        separators=(",", ":"),
    ).encode()
    return SESSION_BUNDLE_MAGIC + struct.pack("<I", len(manifest)) + manifest


class FrontEnd:
    """The widget's front end: one output (``view``) of the widget."""

    def __init__(self, widget: Widget, view: str = "out-1") -> None:
        self.widget = widget
        self.view = view
        self.sent: list[tuple[dict[str, Any], list[bytes]]] = []
        self._lock = threading.Lock()

    def receive(self, content: dict[str, Any], buffers: list[bytes] | None = None) -> None:
        with self._lock:
            self.sent.append((content, [bytes(buffer) for buffer in buffers or []]))

    def post(self, content: dict[str, Any], buffers: list[bytes] | None = None) -> None:
        self.widget._on_page_message(self.widget, {**content, "v": self.view}, [memoryview(b) for b in buffers or []])

    def messages(self, kind: str | None = None) -> list[tuple[dict[str, Any], list[bytes]]]:
        with self._lock:
            return [
                (content, buffers)
                for content, buffers in self.sent
                if content["v"] == self.view and (kind is None or content["t"] == kind)
            ]

    def wait_for(self, predicate, timeout: float = 5.0):
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            for content, buffers in self.messages():
                if predicate(content):
                    return content, buffers
            time.sleep(0.01)
        raise AssertionError(f"the widget did not send it; it sent {[c for c, _ in self.messages()]}")

    def fetch(self, *requests: dict[str, Any]) -> list[tuple[dict[str, Any], bytes]]:
        before = len(self.messages("answers"))
        self.post({"t": "fetch", "requests": list(requests)})
        answered = []
        for content, buffers in self.messages("answers")[before:]:
            for answer in content["answers"]:
                answered.append((answer, buffers[answer["buffer"]] if answer["buffer"] >= 0 else b""))
        return answered

    def connect(self, socket_id: int = 0) -> None:
        widget = self.widget
        self.post({"t": "open", "s": socket_id})
        self.post(
            {
                "t": "text",
                "s": socket_id,
                "text": json.dumps(
                    {"type": "hello", "viewerId": widget._viewer_id, "viewerToken": widget._viewer_token}
                ),
            }
        )
        self.wait_for(lambda content: content["t"] == "text" and json.loads(content["text"]) == {"type": "welcome"})


@pytest.fixture
def widget(tmp_path: Path):
    adata = ad.AnnData(
        X=np.arange(12, dtype=np.float32).reshape(4, 3),
        obs=pd.DataFrame({"cell_type": pd.Categorical(["a", "b", "a", "b"])}, index=list("wxyz")),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
        obsm={"X_umap": np.arange(8, dtype=np.float32).reshape(4, 2)},
    )
    # In priority order (prepare()'s layout), so the file is served as it is
    # and its bytes are what the page reads.
    priority = CellPriority.of(adata.n_obs, 0)
    assert priority is not None
    adata = adata[priority.order].copy()
    adata.obs[ORIGINAL_POSITION_COLUMN] = priority.original.astype(np.uint32)
    adata.uns["cellucid"] = {"priority": priority.record(prefix_levels(4))}
    path = tmp_path / "demo.h5ad"
    adata.write_h5ad(path)
    made = Widget(path, dataset_name="Demo", dataset_id="demo", height=480)
    try:
        yield made
    finally:
        made.stop()
        made.close()


@pytest.fixture
def front(widget: Widget):
    front_end = FrontEnd(widget)
    widget.send = front_end.receive  # type: ignore[method-assign]
    return front_end


def test_the_credentials_go_only_to_the_output_that_asked_and_never_into_the_widget_state(widget, front):
    assert widget.height == 480
    assert widget.web_source_url == "https://www.cellucid.com"
    state = json.dumps({key: str(value) for key, value in widget.get_state().items()})
    assert widget._viewer_token not in state
    front.post({"t": "attach"})
    (content, buffers), = front.messages("credentials")
    assert content == {
        "t": "credentials",
        "v": "out-1",
        "viewerId": widget._viewer_id,
        "viewerToken": widget._viewer_token,
    }
    assert buffers == []


def test_the_page_reads_the_documents_and_the_store_by_range(widget, front, tmp_path):
    answers = front.fetch(
        {"id": 1, "method": "GET", "path": "_cellucid/health", "range": None},
        {"id": 2, "method": "GET", "path": "_cellucid/info", "range": None},
        {"id": 3, "method": "GET", "path": "_cellucid/datasets", "range": None},
        {"id": 4, "method": "GET", "path": "_cellucid/protocol", "range": None},
    )
    assert [answer["status"] for answer, _ in answers] == [200, 200, 200, 200]
    health, info, datasets, protocol = (json.loads(body) for _, body in answers)
    assert health["status"] == "ok" and health["n_cells"] == 4 and health["n_genes"] == 3
    assert info["runtime"] == "jupyter" and info["host"] is None and info["port"] is None
    assert protocol == PROTOCOL["protocolDocument"]
    file = datasets["datasets"][0]["anndata"]["file"]
    assert datasets["datasets"][0]["id"] == "demo" and file["kind"] == "h5ad"

    stored = (tmp_path / "demo.h5ad").read_bytes()
    assert file["size"] == len(stored)
    (ranged, body), (whole, whole_body), (head, head_body), (suffix, suffix_body) = front.fetch(
        {"id": 5, "method": "GET", "path": file["url"], "range": "bytes=8-23"},
        {"id": 6, "method": "GET", "path": file["url"], "range": None},
        {"id": 7, "method": "HEAD", "path": file["url"], "range": None},
        {"id": 8, "method": "GET", "path": file["url"], "range": "bytes=-4"},
    )
    assert ranged["status"] == 206 and body == stored[8:24]
    assert ranged["headers"]["Content-Range"] == f"bytes 8-23/{len(stored)}"
    assert ranged["headers"]["ETag"] == whole["headers"]["ETag"]
    # The list states the version the file's bytes are answered with.
    assert file["etag"] == whole["headers"]["ETag"]
    assert whole["status"] == 200 and whole_body == stored
    assert head["status"] == 200 and head_body == b"" and head["headers"]["Content-Length"] == str(len(stored))
    assert suffix["status"] == 206 and suffix_body == stored[-4:]

    (missing, _), (refused, _), (outside, _) = front.fetch(
        {"id": 9, "method": "GET", "path": "_cellucid/anndata/other.h5ad", "range": None},
        {"id": 10, "method": "POST", "path": "_cellucid/health", "range": None},
        {"id": 11, "method": "GET", "path": file["url"], "range": f"bytes={len(stored)}-"},
    )
    assert missing["status"] == 404
    assert refused["status"] == 405
    assert outside["status"] == 416 and outside["headers"]["Content-Range"] == f"bytes */{len(stored)}"


def test_answers_that_arrive_together_share_messages_of_bounded_size(widget, front, monkeypatch):
    import cellucid.jupyter._widget as module

    monkeypatch.setattr(module, "_ANSWER_MESSAGE_BYTES", 64)
    url = "_cellucid/anndata/demo.h5ad"
    before = len(front.messages("answers"))
    front.fetch(*({"id": index, "method": "GET", "path": url, "range": f"bytes=0-{39 + index}"} for index in range(4)))
    messages = front.messages("answers")[before:]
    # 40, 41, 42, 43 bytes with 64 bytes to a message: one body per message.
    assert [len(content["answers"]) for content, _ in messages] == [1, 1, 1, 1]
    assert [len(buffers[0]) for _, buffers in messages] == [40, 41, 42, 43]
    assert _ANSWER_MESSAGE_BYTES == 8 * 1024 * 1024


def test_the_channel_runs_over_the_comm_commands_first_then_events_in_order(widget, front):
    widget.set_color_by("cell_type")  # sent before the page connected: it waits
    front.connect()
    front.wait_for(lambda c: c["t"] == "text" and json.loads(c["text"]) == {"type": "setColorBy", "field": "cell_type"})

    ready = _example_event("ready")
    front.post({"t": "text", "s": 0, "text": json.dumps(ready)})
    deadline = time.monotonic() + 5
    while widget.state.ready is None and time.monotonic() < deadline:
        time.sleep(0.01)
    assert widget.state.ready is not None

    # A selection, sent in pieces as the front end splits a binary frame.
    example = next(e for e in PROTOCOL["examples"]["binaryEvents"] if e["header"]["type"] == "selection")
    frame = bytes.fromhex(example["hex"])
    seen = []
    widget.on_selection(seen.append)
    front.post({"t": "binary", "s": 0, "last": False}, [frame[:3]])
    front.post({"t": "binary", "s": 0, "last": False}, [frame[3:11]])
    front.post({"t": "binary", "s": 0, "last": True}, [frame[11:]])
    deadline = time.monotonic() + 5
    while not seen and time.monotonic() < deadline:
        time.sleep(0.01)
    assert [list(event["cells"]) for event in seen] == [example["cells"]]

    # A command that names cells is a binary frame, sent whole as one buffer.
    command = PROTOCOL["examples"]["binaryCommands"][0]
    widget.send_message({**command["header"], "cells": command["cells"]})
    content, buffers = front.wait_for(lambda c: c["t"] == "binary")
    assert content["s"] == 0 and buffers == [bytes.fromhex(command["hex"])]


def test_a_wrong_token_is_refused_and_a_second_output_replaces_the_first(widget, front):
    front.post({"t": "open", "s": 0})
    front.post({"t": "text", "s": 0, "text": json.dumps({"type": "hello", "viewerId": widget._viewer_id, "viewerToken": "x" * 32})})
    content, _ = front.wait_for(lambda c: c["t"] == "close")
    assert content["code"] == 4401

    first = FrontEnd(widget, "out-1")
    second = FrontEnd(widget, "out-2")

    def both(content, buffers=None):
        first.receive(content, buffers)
        second.receive(content, buffers)

    widget.send = both  # type: ignore[method-assign]
    first.connect(socket_id=1)
    second.connect(socket_id=0)
    content, _ = first.wait_for(lambda c: c["t"] == "close")
    assert content["code"] == 4409 and content["s"] == 1
    assert not second.messages("close")


def test_a_cell_cannot_wait_for_the_viewer_but_another_thread_can(widget, front):
    front.post({"t": "attach"})  # the kernel hands the widget's messages to this thread
    with pytest.raises(RuntimeError, match="only between cells"):
        widget.wait_for_event("hover", timeout=30)
    with pytest.raises(RuntimeError, match="@viewer.on_ready"):
        widget.wait_for_ready(timeout=30)
    request = widget.request_session_bundle()
    with pytest.raises(RuntimeError, match="later cell|next cell"):
        request.result(timeout=30)
    request.cancel()

    front.connect()
    caught: list[Any] = []
    waiting = threading.Thread(target=lambda: caught.append(widget.wait_for_event("hover", timeout=10)))
    waiting.start()
    time.sleep(0.1)
    hover = _example_event("hover")
    front.post({"t": "text", "s": 0, "text": json.dumps(hover)})
    waiting.join(10)
    assert caught and caught[0]["cell"] == hover["cell"]


def test_a_session_bundle_is_asked_in_one_cell_and_collected_in_the_next(widget, front):
    front.post({"t": "attach"})
    front.connect()
    request = widget.request_session_bundle()
    content, _ = front.wait_for(lambda c: c["t"] == "text" and json.loads(c["text"])["type"] == "requestSessionBundle")
    request_id = json.loads(content["text"])["requestId"]
    assert request.done is False
    with pytest.raises(RuntimeError, match="only between cells"):
        request.result()

    # Between cells: the page sends the bundle, in pieces.
    frame = _binary_frame({"type": "session_bundle", "requestId": request_id}, _bundle())
    front.post({"t": "binary", "s": 0, "last": False}, [frame[:20]])
    front.post({"t": "binary", "s": 0, "last": True}, [frame[20:]])
    deadline = time.monotonic() + 5
    while not request.done and time.monotonic() < deadline:
        time.sleep(0.01)
    bundle = request.result()
    try:
        assert bundle.path.read_bytes() == _bundle()
        assert request.result() is bundle
    finally:
        bundle.path.unlink()


def test_stopping_freezes_the_page_closes_its_channel_and_answers_no_more_reads(widget, front):
    front.connect()
    widget.stop()
    content, _ = front.wait_for(lambda c: c["t"] == "close")
    assert content["code"] == 1001
    texts = [json.loads(c["text"])["type"] for c, _ in front.messages("text")]
    assert texts[-1] == "freeze"
    (answer, body), = front.fetch({"id": 1, "method": "GET", "path": "_cellucid/health", "range": None})
    assert answer["status"] == 410
    assert repr(widget) == "cellucid.Widget(dataset=None, status=stopped)"


def test_a_removed_output_closes_its_channel(widget, front):
    front.connect()
    registration_socket = widget._sockets[("out-1", 0)]
    front.post({"t": "detach"})
    deadline = time.monotonic() + 5
    while not registration_socket.closed and time.monotonic() < deadline:
        time.sleep(0.01)
    assert registration_socket.closed
    assert ("out-1", 0) not in widget._sockets


@pytest.mark.parametrize(
    "web_source_url",
    ["https://www.cellucid.com/", "ftp://x.org", "https://x.org/?a=1", "https://u:p@x.org", " https://x.org", "x"],
)
def test_the_app_address_is_exact(web_source_url, tmp_path):
    adata = ad.AnnData(X=np.zeros((2, 1), np.float32), obsm={"X_umap": np.zeros((2, 2), np.float32)})
    with pytest.raises(ValueError, match="web_source_url"):
        Widget(adata, dataset_name="Demo", dataset_id="demo", web_source_url=web_source_url)


def test_a_self_hosted_app_address_is_accepted(tmp_path):
    adata = ad.AnnData(X=np.zeros((2, 1), np.float32), obsm={"X_umap": np.zeros((2, 2), np.float32)})
    made = Widget(adata, dataset_name="Demo", dataset_id="demo", web_source_url="https://cellucid.example.org/app")
    try:
        assert made.web_source_url == "https://cellucid.example.org/app"
        with pytest.raises(AttributeError, match="no viewer address"):
            made.viewer_url  # noqa: B018
    finally:
        made.stop()
        made.close()


def test_outside_a_notebook_there_is_nowhere_to_show_it(widget):
    widget._context["in_jupyter"] = False
    with pytest.raises(RuntimeError, match="browser tab"):
        show_widget(widget)


def test_the_cell_that_showed_it_does_not_show_it_twice(widget, monkeypatch):
    import cellucid.jupyter._widget as module

    count = {"value": 7}
    monkeypatch.setattr(module, "_execution_count", lambda: count["value"])
    shown = []
    monkeypatch.setattr("IPython.display.display", lambda obj: shown.append(obj))
    widget.display()
    assert shown == [widget]
    # The same cell's last expression: text only, no second view.
    assert widget._repr_mimebundle_() == {"text/plain": repr(widget)}
    assert widget._repr_html_() is None
    # A later cell shows it again, as any widget.
    count["value"] = 8
    bundle = widget._repr_mimebundle_()
    bundle = bundle[0] if isinstance(bundle, tuple) else bundle
    assert "application/vnd.jupyter.widget-view+json" in bundle
    assert "notebook widget" in widget._repr_html_()


# The widget's front end, byte for byte cellucid-r's inst/widget/: pinned here
# and in cellucid-r's tests/testthat/test-viewer-widget.R.
WIDGET_JS_SHA256 = "e296bc2cfba93560f24d35792fe84934b707f5bf3fca82083fa4c39db9ebbc04"
WIDGET_CSS_SHA256 = "3b9acad3ebc04ece34acf580446a72ed7f73f3497b747539a11bc33fe281d82b"


def test_the_front_end_is_cellucid_rs_byte_for_byte():
    import hashlib

    from cellucid.jupyter._widget import _STATIC

    assert hashlib.sha256((_STATIC / "widget.js").read_bytes()).hexdigest() == WIDGET_JS_SHA256
    assert hashlib.sha256((_STATIC / "widget.css").read_bytes()).hexdigest() == WIDGET_CSS_SHA256


def test_the_front_end_module_ships_with_the_package():
    from cellucid.jupyter._widget import _STATIC

    source = (_STATIC / "widget.js").read_text()
    assert "export default { render }" in source
    # The token is asked of the kernel, kept in memory and posted only to the app's origin.
    assert "localStorage" not in source and "console.log" not in source
    assert "postMessage(\n      { [EMBED_KEY]: 'connect'" in source
    assert "event.source !== frame.contentWindow || event.origin !== appOrigin" in source
