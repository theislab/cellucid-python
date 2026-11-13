"""The round trip in a notebook widget, against a mock comm.

A widget hears the viewer only between cells (the kernel reads comm messages
when it is idle), so nothing in the round trip may ask the viewer and wait.
The page here sends what ``app/viewer-host/results.js`` sends -- a snapshot
(every result, then the ``results`` inventory) and each change after -- over
the widget's comm, between "cells"; the next "cell" pulls, subsets, lists the
gene sets and pushes, all at once, on the thread a cell runs on.
"""

from __future__ import annotations

import json
import struct
import time
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid.jupyter import Widget

NAMES = list("wxyz")


def _frame(header: dict[str, Any], body: bytes) -> bytes:
    encoded = json.dumps(header, separators=(",", ":")).encode()
    return struct.pack("<I", len(encoded)) + encoded + body


def _result(kind: str, result_id: str, meta: dict[str, Any], cells: list[int] | None, values: np.ndarray | None,
            revision: int) -> bytes:
    value_type = None if values is None else {"u1": "uint8", "f4": "float32"}[values.dtype.str[1:]]
    body = (b"" if cells is None else np.asarray(cells, "<u4").tobytes()) + (b"" if values is None else values.tobytes())
    header = {"type": "result", "epoch": "0123456789abcdef", "revision": revision, "kind": kind, "id": result_id, "nObs": 4,
              "meta": meta, "cells": None if cells is None else len(cells), "values": value_type}
    return _frame(header, body)


class Page:
    """The widget's front end: one output, its channel open."""

    def __init__(self, widget: Widget) -> None:
        self.widget = widget
        self.sent: list[dict[str, Any]] = []
        widget.send = lambda content, buffers=None: self.sent.append(  # type: ignore[method-assign]
            {**content, "buffers": [bytes(buffer) for buffer in buffers or []]}
        )

    def post(self, content: dict[str, Any], buffers: list[bytes] | None = None) -> None:
        self.widget._on_page_message(self.widget, {**content, "v": "out-1"}, [memoryview(b) for b in buffers or []])

    def connect(self) -> None:
        self.post({"t": "attach"})  # the kernel hands the widget's messages to this thread
        self.post({"t": "open", "s": 0})
        self.post({"t": "text", "s": 0, "text": json.dumps(
            {"type": "hello", "viewerId": self.widget._viewer_id, "viewerToken": self.widget._viewer_token})})

    def text(self, event: dict[str, Any]) -> None:
        self.post({"t": "text", "s": 0, "text": json.dumps(event)})

    def binary(self, frame: bytes) -> None:
        self.post({"t": "binary", "s": 0, "last": True}, [frame])


def _settle(predicate, timeout: float = 5.0) -> None:
    """Between cells: the kernel has read what the page sent."""
    deadline = time.monotonic() + timeout
    while not predicate():
        if time.monotonic() > deadline:
            raise AssertionError("the widget did not take the page's messages")
        time.sleep(0.01)


@pytest.fixture
def adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.zeros((4, 2), dtype=np.float32),
        obs=pd.DataFrame({"cell_type": pd.Categorical(["a", "b", "a", "b"])}, index=NAMES),
        var=pd.DataFrame(index=["g1", "g2"]),
        obsm={"X_umap": np.zeros((4, 2), dtype=np.float32)},
    )


@pytest.fixture
def widget(tmp_path: Path, adata: ad.AnnData):
    path = tmp_path / "demo.h5ad"
    adata.write_h5ad(path)
    made = Widget(path, dataset_name="Demo", dataset_id="demo")
    try:
        yield made
    finally:
        made.stop()
        made.close()


def test_a_lasso_then_pull_in_the_next_cell_needs_no_waiting(widget: Widget, adata: ad.AnnData) -> None:
    page = Page(widget)
    page.connect()
    # Cell 1 ends; the page opened its dataset and sends its snapshot.
    page.text({"type": "ready", "n_cells": 4, "n_obs": 4, "dimensions": 2})
    page.binary(_result("highlight", "highlight_1",
                        {"label": "Lasso", "tool": "lasso", "pageId": "page_1", "pageName": "Page 1",
                         "color": "#e41a1c", "enabled": True}, [1, 3], None, 1))
    page.binary(_result("categorical", "user_cat_1", {"key": "annotation", "categories": ["B", "T"], "colors": ["#ff9900", "#0077c8"], "origin": "user"},
                        None, np.array([0, 1, 255, 0], dtype=np.uint8), 2))
    page.binary(_result("geneSet", "gene_set_1",
                        {"name": "T", "description": "", "genes": ["CD3E"], "symbols": ["CD3E"], "origin": "markers"}, [], None, 3))
    page.text({"type": "results", "epoch": "0123456789abcdef", "dataset": "demo", "revision": 3, "nObs": 4, "items": [
        {"kind": "highlight", "id": "highlight_1", "revision": 1},
        {"kind": "categorical", "id": "user_cat_1", "revision": 2},
        {"kind": "geneSet", "id": "gene_set_1", "revision": 3},
    ]})
    _settle(lambda: widget._results_inventory is not None and len(widget._results) == 3)

    # Cell 2: everything reads what arrived, at once, on the cell's thread.
    started = time.monotonic()
    summary = widget.pull(adata)
    assert time.monotonic() - started < 1.0
    assert adata.obs["cellucid_highlight__highlight_1"].tolist() == [False, True, False, True]
    assert adata.obs["annotation"].astype(object).where(adata.obs["annotation"].notna(), None).tolist() == [
        "B", "T", None, "B"]
    assert set(summary.columns) == {"cellucid_highlight__highlight_1", "annotation"}
    assert list(widget.subset(adata, "Lasso").obs_names) == ["x", "z"]
    assert widget.gene_sets() == {"T": ["CD3E"]}

    # Cell 2 also pushes: the frame goes out at once, nothing waits.
    request = widget.push(adata.obs["cell_type"], key="cell_type_nb")
    widget.gene_set_save(name="B", genes=["MS4A1"])
    pushed = [message for message in page.sent if message["t"] == "binary"]
    header_length = struct.unpack("<I", pushed[-1]["buffers"][0][:4])[0]
    header = json.loads(pushed[-1]["buffers"][0][4 : 4 + header_length])
    assert header["type"] == "action" and header["action"] == "field.push" and header["requestId"] == request.request_id
    assert header["params"]["key"] == "cell_type_nb" and header["params"]["categories"] == ["a", "b"]
    assert any(
        m["t"] == "text" and json.loads(m["text"]).get("action") == "gene-set.save" for m in page.sent
    )
    # Its answer arrives between cells: this cell cannot wait for it.
    with pytest.raises(RuntimeError, match="later cell"):
        request.result(timeout=30)
    request.cancel()

    # Between cells: the person deletes the lasso; the next pull drops its column.
    page.text({"type": "result_removed", "epoch": "0123456789abcdef", "revision": 4, "kind": "highlight", "id": "highlight_1"})
    _settle(lambda: ("highlight", "highlight_1") not in widget._results)
    widget.pull(adata)
    assert "cellucid_highlight__highlight_1" not in adata.obs


def test_before_the_snapshot_a_pull_says_what_to_do_at_once(widget: Widget, adata: ad.AnnData) -> None:
    Page(widget).connect()
    started = time.monotonic()
    with pytest.raises(RuntimeError, match="has not sent what it holds yet"):
        widget.pull(adata)
    assert time.monotonic() - started < 1.0


def test_a_widget_refuses_results_of_a_dataset_it_did_not_serve(widget: Widget, adata: ad.AnnData) -> None:
    assert widget._served_dataset_id == "demo"
    page = Page(widget)
    page.connect()
    # The person opened another dataset in the viewer: a new epoch, no host dataset.
    page.text({"type": "results", "epoch": "1111111111111111", "dataset": None, "revision": 0, "nObs": 2700, "items": []})
    _settle(lambda: widget._results_epoch == "1111111111111111")
    with pytest.raises(RuntimeError, match="did not open: another one was opened in it.*'demo'"):
        widget.pull(adata)
    page.text({"type": "results", "epoch": "2222222222222222", "dataset": "other", "revision": 0, "nObs": 4, "items": []})
    _settle(lambda: widget._results_epoch == "2222222222222222")
    with pytest.raises(RuntimeError, match="describe the dataset 'other', and this viewer serves 'demo'"):
        widget.gene_sets()


def test_in_a_widget_a_session_is_asked_for_in_one_cell_and_applied_in_the_next(widget: Widget, adata: ad.AnnData) -> None:
    Page(widget).connect()
    request = widget.request_session()
    started = time.monotonic()
    with pytest.raises(RuntimeError, match="later cell|next cell"):
        widget.apply_session_to_anndata(adata, request)
    with pytest.raises(TypeError, match="request_session"):
        widget.apply_session_to_anndata(adata, "a session")
    assert time.monotonic() - started < 1.0
    request.cancel()


def test_a_request_collected_in_a_later_cell_says_it_has_not_arrived_yet(
    widget: Widget, adata: ad.AnnData, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Asked in one cell and collected in the next, as the two-cell pattern
    says, an answer that has not come is one the viewer has not sent yet --
    never "ask in one cell and collect it in the next", which is what was
    done."""
    from cellucid.jupyter import _viewer, _widget

    cell = {"count": 1}
    for module in (_viewer, _widget):
        monkeypatch.setattr(module, "_execution_count", lambda: cell["count"])
    request = widget.request_session()
    cell["count"] = 2
    with pytest.raises(RuntimeError, match="session has not arrived yet: the viewer is not shown"):
        request.result()
    widget._displayed = True
    with pytest.raises(RuntimeError, match="the viewer has not connected"):
        request.result()
    Page(widget).connect()
    _settle(lambda: widget.debug_connection()["channel_connected"])
    with pytest.raises(RuntimeError, match="the viewer is still working on it.*Run this cell again"):
        widget.apply_session_to_anndata(adata, request)
    # Collected in the cell that asked: the two-cell advice.
    cell["count"] = 3
    again = widget.request_session()
    with pytest.raises(RuntimeError, match="next cell"):
        again.result()
    request.cancel()
    again.cancel()


@pytest.mark.parametrize("served", ["file", "object"])
def test_a_viewer_given_another_dataset_id_takes_its_own_session_back(
    tmp_path: Path, adata: ad.AnnData, served: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A viewer given dataset_id= serves the data under that id, and saves
    its session under it; applied back through the viewer, the session goes
    to the object it served, which declares its own id. An object
    declaring another id is refused."""
    import _session_files as files

    from cellucid.anndata_session._apply import cell_order_digest

    adata.uns["cellucid"] = {"id": "pancreas"}
    path = tmp_path / "pancreas.h5ad"
    adata.write_h5ad(path)
    made = Widget(path if served == "file" else adata, dataset_name="Pancreas", dataset_id="pancreas-notebook")
    try:
        assert (made._served_dataset_id, made._served_declared_id) == ("pancreas-notebook", "pancreas")
        record = files.fingerprint(cell_count=4, var_count=2, dataset_id="pancreas-notebook",
                                   digest=cell_order_digest(adata.obsm["X_umap"]))
        made_field = files.user_category("user_cat_1", "stage", ["early", "late"], [0, 1, 0, 1])
        session = files.session_file(
            tmp_path / "viewer.cellucid-session", fingerprint_record=record,
            trackers=[{"id": "t1", "name": "Tracker 1", "root": files.state(),
                       "steps": [{"set": {"field-overlays": files.overlays(user_fields=[made_field])}}]}],
        )
        from cellucid import read_session

        def answered(request: Any) -> None:
            monkeypatch.setattr(request, "result", lambda timeout=None: read_session(session))

        request = made.request_session()
        answered(request)
        plan = made.apply_session_to_anndata(adata, request, cleanup_file=False)
        assert [row.status for row in plan.rows] == ["write"]
        plan.apply()
        assert adata.obs["stage"].tolist() == ["early", "late", "early", "late"]
        other = adata.copy()
        other.uns["cellucid"] = {"id": "garcia"}
        request = made.request_session()
        answered(request)
        with pytest.raises(ValueError, match="saved on the dataset 'pancreas-notebook', and this object is 'garcia'"):
            made.apply_session_to_anndata(other, request, cleanup_file=False)
    finally:
        made.stop()
        made.close()
