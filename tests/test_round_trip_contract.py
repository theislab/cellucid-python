"""The round trip: what a person made in the viewer, pulled into an AnnData in
one call, and columns and gene sets pushed into the viewer.

The viewer's side is played by a fake that sends what ``app/viewer-host-results.js``
sends -- a snapshot (every result, then the ``results`` inventory), then each
change -- through the notebook's own event path (``_handle_frontend_message``),
frames and all. Nothing in the round trip asks the viewer and waits: every
call reads what it already sent.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from _serverless_viewer import serverless_viewer

import cellucid
from cellucid.anndata_session._results import field_message, read_obs_names
from cellucid.jupyter import _wire as wire

N_OBS = 6
#: The epoch of the fake page's results (one page's account of one dataset).
EPOCH = "0123456789abcdef"
OBS_NAMES = [f"cell-{index}" for index in range(N_OBS)]


def _adata() -> ad.AnnData:
    # CELLxGENE-shaped: the var index is Ensembl ids, the symbols in feature_name.
    return ad.AnnData(
        X=np.zeros((N_OBS, 2), dtype=np.float32),
        obs=pd.DataFrame({"leiden": pd.Categorical(["0", "1", "0", "2", "1", "0"])}, index=OBS_NAMES),
        var=pd.DataFrame({"feature_name": ["CD3E", "MS4A1"]}, index=["ENSG00000198851", "ENSG00000156738"]),
    )


def _result(
    kind: str,
    result_id: str,
    meta: dict[str, Any],
    *,
    cells: list[int] | None,
    values: np.ndarray | None = None,
    revision: int = 1,
    request_id: str | None = None,
    n_obs: int = N_OBS,
    epoch: str = EPOCH,
) -> dict[str, Any]:
    """One ``result`` event as the channel hands it over: header fields and body."""
    body = b""
    if cells is not None:
        body += np.asarray(cells, dtype="<u4").tobytes()
    value_type = None
    if values is not None:
        value_type = {"u1": "uint8", "u2": "uint16", "u4": "uint32", "f4": "float32", "f8": "float64"}[
            values.dtype.str[1:]
        ]
        body += values.astype(values.dtype.newbyteorder("<")).tobytes()
    return {
        "type": "result",
        "requestId": request_id,
        "epoch": epoch,
        "revision": revision,
        "kind": kind,
        "id": result_id,
        "nObs": n_obs,
        "meta": meta,
        "cells": None if cells is None else len(cells),
        "values": value_type,
        "body": bytearray(body),
    }


def _highlight(result_id: str, cells: list[int], **meta: Any) -> dict[str, Any]:
    return _result(
        "highlight",
        result_id,
        {"label": meta.get("label", "Lasso selection"), "tool": "lasso", "pageId": "page_1",
         "pageName": "Page 1", "color": "#e41a1c", "enabled": True},
        cells=cells,
    )


def _viewer_results(*, open_cells: list[int] | None = None) -> list[dict[str, Any]]:
    """What a viewer holds: two highlights, a label field, a score, a pushed
    column, a gene set, the visible cells and the open cells."""
    if open_cells is None:
        open_count = N_OBS
        cats = np.array([0, 1, 0, 255, 1, 0], dtype=np.uint8)
    else:
        open_count = len(open_cells)
        cats = np.array([0, 1, 255, 1][:open_count], dtype=np.uint8)
    return [
        _result("subset", "open", {"digest": "0123456789abcdef"}, cells=open_cells),
        _highlight("highlight_1", [0, 2]),
        _highlight("highlight_2", [3], label="Picked"),
        _result("categorical", "user_cat_1", {"key": "annotation", "categories": ["B", "T"], "origin": "user"},
                cells=open_cells, values=cats),
        _result("continuous", "signature:T score", {"key": "T score", "origin": "analysis"},
                cells=[1, 3], values=np.array([0.5, -1.0], dtype=np.float32)),
        _result("continuous", "user_obs_value_1", {"key": "leiden_copy", "origin": "notebook"},
                cells=open_cells, values=np.zeros(open_count, dtype=np.float32)),
        _result("geneSet", "gene_set_1", {"name": "T", "description": "", "genes": ["ENSG00000198851"],
                                          "symbols": ["CD3E"], "origin": "notebook"},
                cells=[]),
        _result("visibility", "visible", {}, cells=[0, 1, 2]),
    ]


def send_snapshot(viewer: Any, held: list[dict[str, Any]], *, epoch: str = EPOCH, dataset: str = "demo") -> None:
    """What the viewer sends once its dataset opened: every result, then the inventory."""
    items = []
    for revision, event in enumerate(held, start=1):
        viewer._handle_frontend_message({**event, "epoch": epoch, "revision": revision})
        items.append({"kind": event["kind"], "id": event["id"], "revision": revision})
    viewer._handle_frontend_message(
        {"type": "results", "epoch": epoch, "dataset": dataset, "revision": len(held), "nObs": N_OBS, "items": items}
    )


class FakeViewer:
    """The viewer: its snapshot sent, and its answers to pushed columns and gene sets."""

    def __init__(self, viewer: Any, held: list[dict[str, Any]]) -> None:
        self.viewer = viewer
        self.sent: list[Any] = []
        viewer.state.ready = {"n_cells": N_OBS, "n_obs": N_OBS, "dimensions": 2}
        viewer._obs_names_cache = pd.Index(OBS_NAMES)
        send_snapshot(viewer, held)

    def send(self, message: dict[str, Any]) -> bool:
        frame = wire._require_frontend_message(message)
        self.sent.append((message, frame))
        if message["type"] == "setField":
            meta = (
                {"key": message["key"], "categories": message["categories"], "origin": "notebook"}
                if message["kind"] == "categorical"
                else {"key": message["key"], "origin": "notebook"}
            )
            cells = None if message["cells"] is None else np.asarray(message["cells"]).tolist()
            self.viewer._handle_frontend_message(
                _result(message["kind"], f"pushed:{message['key']}", meta, cells=cells,
                        values=np.asarray(message["values"]), request_id=message["requestId"], revision=50)
            )
        elif message["type"] == "setGeneSet":
            self.viewer._handle_frontend_message(
                _result("geneSet", "gene_set_1", {"name": message["name"], "description": message["description"],
                                                  "genes": message["genes"], "symbols": None,
                                                  "origin": "notebook"},
                        cells=[], request_id=message["requestId"], revision=51)
            )
        return True


@pytest.fixture
def connected(monkeypatch: pytest.MonkeyPatch):
    def make(held: list[dict[str, Any]] | None = None):
        viewer = serverless_viewer()
        fake = FakeViewer(viewer, _viewer_results() if held is None else held)
        monkeypatch.setattr(viewer, "send_message", fake.send)
        return viewer, fake

    return make


def test_one_pull_writes_every_highlight_label_and_score_with_one_naming_scheme(connected) -> None:
    viewer, _fake = connected()
    adata = _adata()
    summary = viewer.pull(adata)
    obs = adata.obs
    assert obs["cellucid_highlight__highlight_1"].tolist() == [True, False, True, False, False, False]
    assert obs["cellucid_highlight__highlight_1"].dtype == bool
    assert obs["cellucid_highlight__highlight_2"].tolist() == [False, False, False, True, False, False]
    assert list(obs["annotation"].cat.categories) == ["B", "T"]
    assert obs["annotation"].astype(object).where(obs["annotation"].notna(), None).tolist() == [
        "B", "T", "B", None, "T", "B",
    ]
    assert obs["T score"].dtype == np.float32
    np.testing.assert_array_equal(obs["T score"].to_numpy(), [np.nan, 0.5, np.nan, -1.0, np.nan, np.nan])
    # A column the notebook pushed comes back only when asked for, and the
    # visible cells only when asked for.
    assert "leiden_copy" not in obs and "cellucid_visible" not in obs
    assert summary.skipped == {"continuous:user_obs_value_1": "pushed from a notebook (include_pushed=True brings it back)"}
    assert set(summary.columns) == {
        "cellucid_highlight__highlight_1", "cellucid_highlight__highlight_2", "annotation", "T score",
    }
    assert summary.columns["cellucid_highlight__highlight_2"]["label"] == "Picked"
    assert adata.uns["cellucid"]["gene_sets"] == {"T": {
        "genes": ["ENSG00000198851"], "symbols": ["CD3E"], "description": "", "origin": "notebook",
        "pulled_from": "viewer", "id": "gene_set_1"}}
    # The genes name rows of var as this object names them (CELLxGENE: Ensembl ids).
    assert set(adata.uns["cellucid"]["gene_sets"]["T"]["genes"]) <= set(adata.var_names)
    assert list(adata.var.loc[adata.uns["cellucid"]["gene_sets"]["T"]["genes"], "feature_name"]) == ["CD3E"]
    assert set(adata.uns["cellucid"]["pulled"]) == set(summary.columns)
    assert adata.uns["cellucid"]["pulled"]["annotation"] == {
        "kind": "categorical", "id": "user_cat_1", "pulled_from": "viewer", "key": "annotation", "origin": "user"}


def test_a_second_pull_replaces_its_own_columns_and_removes_the_gone_ones(connected) -> None:
    viewer, _fake = connected()
    adata = _adata()
    viewer.pull(adata, visibility=True)
    assert adata.obs["cellucid_visible"].tolist() == [True, True, True, False, False, False]
    # The person deletes a highlight: the viewer says so, and the next pull removes its column.
    viewer._handle_frontend_message({"type": "result_removed", "epoch": EPOCH, "revision": 20, "kind": "highlight", "id": "highlight_2"})
    summary = viewer.pull(adata)
    assert "cellucid_highlight__highlight_2" not in adata.obs
    # The visible cells were left out of this pull, not removed in the viewer.
    assert "cellucid_visible" in adata.obs
    assert summary.removed == ["cellucid_highlight__highlight_2"]
    assert "cellucid_highlight__highlight_1" in adata.obs


def test_prune_removes_only_what_the_viewer_no_longer_holds(connected) -> None:
    viewer, _fake = connected()
    adata = _adata()
    viewer.pull(adata, include_pushed=True)
    assert "leiden_copy" in adata.obs
    # Results left out of a pull keep their columns: highlights, pushed fields.
    summary = viewer.pull(adata, highlights=False)
    assert summary.removed == [] and "cellucid_highlight__highlight_1" in adata.obs and "leiden_copy" in adata.obs
    # The viewer drops its gene set: the next pull drops it from uns too.
    viewer._handle_frontend_message({"type": "result_removed", "epoch": EPOCH, "revision": 30, "kind": "geneSet",
                                     "id": "gene_set_1"})
    summary = viewer.pull(adata)
    assert summary.removed_gene_sets == ["T"] and adata.uns["cellucid"]["gene_sets"] == {}


def test_a_file_pull_and_a_viewer_pull_keep_each_others_columns(connected, tmp_path: Path) -> None:
    viewer, _fake = connected()
    adata = _adata()
    (tmp_path / "groups.csv").write_text("obs_name,groups\nkind,categorical\ncell-0,a\ncell-1,b\n")
    cellucid.pull(adata, tmp_path / "groups.csv")
    summary = viewer.pull(adata)
    assert "groups" in adata.obs and summary.removed == []
    assert adata.uns["cellucid"]["pulled"]["groups"]["pulled_from"] == "groups.csv"
    # A file column of a name the viewer wants is another source's: a conflict.
    (tmp_path / "annotation.csv").write_text("obs_name,annotation\nkind,categorical\ncell-0,a\n")
    with pytest.raises(ValueError, match="'annotation'"):
        cellucid.pull(adata, tmp_path / "annotation.csv")


def test_a_suffixed_column_is_written_again_where_it_was(connected) -> None:
    viewer, _fake = connected()
    adata = _adata()
    adata.obs["annotation"] = "mine"
    for _ in range(3):
        summary = viewer.pull(adata, on_conflict="suffix")
        assert "annotation_2" in summary.columns and summary.removed == []
    assert [column for column in adata.obs if column.startswith("annotation")] == ["annotation", "annotation_2"]
    for _ in range(2):
        viewer.pull(adata, on_conflict="suffix", prune=False)
    assert [column for column in adata.obs if column.startswith("annotation")] == ["annotation", "annotation_2"]


def test_a_column_cellucid_did_not_write_is_a_conflict_settled_by_on_conflict(connected) -> None:
    viewer, _fake = connected()
    adata = _adata()
    adata.obs["annotation"] = "mine"
    with pytest.raises(ValueError, match="'annotation'.*on_conflict='replace'"):
        viewer.pull(adata)
    assert list(adata.obs.columns) == ["leiden", "annotation"]
    viewer.pull(adata, on_conflict="suffix")
    assert (adata.obs["annotation"] == "mine").all()
    assert "annotation_2" in adata.obs
    summary = viewer.pull(adata, on_conflict="replace")
    assert list(adata.obs["annotation"].cat.categories) == ["B", "T"]
    # The same result, now under the name asked for: its suffixed column goes.
    assert summary.removed == ["annotation_2"] and "annotation_2" not in adata.obs


def test_a_subsampled_viewer_leaves_its_unopened_cells_unknown(connected) -> None:
    viewer, _fake = connected(_viewer_results(open_cells=[0, 2, 3, 4]))
    adata = _adata()
    viewer.pull(adata)
    highlight = adata.obs["cellucid_highlight__highlight_1"]
    assert str(highlight.dtype) == "boolean"
    assert highlight.tolist() == [True, pd.NA, True, False, False, pd.NA]
    annotation = adata.obs["annotation"].astype(object)
    assert annotation.where(annotation.notna(), None).tolist() == ["B", None, "T", None, "T", None]


def test_a_pull_into_another_dataset_is_refused(connected) -> None:
    viewer, _fake = connected()
    with pytest.raises(ValueError, match="dataset of 6 cells and this object has 5"):
        viewer.pull(_adata()[:5].copy())


def test_before_its_snapshot_or_with_one_short_the_viewer_holds_no_results() -> None:
    viewer = serverless_viewer()
    with pytest.raises(RuntimeError, match="has not sent what it holds yet"):
        viewer.pull(_adata())
    viewer._handle_frontend_message(
        {"type": "results", "epoch": EPOCH, "dataset": "demo", "revision": 1, "nObs": N_OBS,
         "items": [{"kind": "highlight", "id": "highlight_1", "revision": 1}]}
    )
    with pytest.raises(RuntimeError, match="incomplete"):
        viewer.gene_sets()


def test_live_results_keep_the_newest_revision_and_a_resync_drops_the_rest() -> None:
    viewer = serverless_viewer()
    send_snapshot(viewer, [])
    viewer._handle_frontend_message({**_highlight("highlight_1", [0]), "revision": 3})
    viewer._handle_frontend_message({**_highlight("highlight_1", [1]), "revision": 2})
    assert [result.cells.tolist() for result in viewer.results] == [[0]]
    viewer._handle_frontend_message({**_highlight("highlight_2", [2]), "revision": 4})
    viewer._handle_frontend_message({"type": "result_removed", "epoch": EPOCH, "revision": 5, "kind": "highlight", "id": "highlight_1"})
    assert [result.id for result in viewer.results] == ["highlight_2"]
    viewer._handle_frontend_message(
        {"type": "results", "epoch": EPOCH, "dataset": "demo", "revision": 6, "nObs": N_OBS, "items": []}
    )
    assert viewer.results == []


def test_a_new_page_replaces_the_old_ones_results_whatever_their_revisions() -> None:
    # The review's case: page 1 sent revisions 15 and 16; the page reloads (or
    # the viewer is shown again in another output) and counts from 1 again.
    viewer = serverless_viewer()
    send_snapshot(viewer, [_highlight("highlight_1", [1, 2, 3]), _result("visibility", "visible", {}, cells=[0, 1])])
    viewer._handle_frontend_message({**_highlight("highlight_1", [1, 2, 3]), "revision": 15})
    viewer._handle_frontend_message({**_result("visibility", "visible", {}, cells=[0, 1]), "revision": 16})
    second = "fedcba9876543210"
    send_snapshot(viewer, [_highlight("highlight_1", [4, 5]), _result("visibility", "visible", {}, cells=[0, 1, 2, 3, 4])],
                  epoch=second)
    assert {result.kind: result.cells.tolist() for result in viewer.results} == {
        "highlight": [4, 5], "visibility": [0, 1, 2, 3, 4]}
    # A late frame of the old page is not mixed in: it starts its own epoch,
    # incomplete until its snapshot ends, and is refused until then.
    viewer._handle_frontend_message({**_highlight("highlight_9", [0]), "revision": 40})
    with pytest.raises(RuntimeError, match="has not sent what it holds yet"):
        _ = viewer.results


def test_results_of_another_dataset_are_refused_until_the_hosts_is_open_again() -> None:
    viewer = serverless_viewer()
    viewer._served_dataset_id = "demo"
    send_snapshot(viewer, [_highlight("highlight_1", [0])])
    assert len(viewer.results) == 1
    viewer._handle_frontend_message(
        {"type": "results", "epoch": "1111111111111111", "dataset": None, "revision": 0, "nObs": 2700, "items": []}
    )
    with pytest.raises(RuntimeError, match="did not open: another one was opened in it.*'demo'"):
        viewer.pull(_adata())
    with pytest.raises(RuntimeError, match="did not open"):
        viewer.gene_sets()
    send_snapshot(viewer, [_highlight("highlight_1", [0])], epoch="2222222222222222", dataset="other")
    with pytest.raises(RuntimeError, match="describe the dataset 'other', and this viewer serves 'demo'"):
        viewer.subset(_adata(), "highlight_1")
    send_snapshot(viewer, [_highlight("highlight_2", [1])], epoch="3333333333333333")
    assert [result.id for result in viewer.results] == ["highlight_2"]


def test_selection_subset_and_gene_sets(connected) -> None:
    viewer, _fake = connected()
    viewer._handle_frontend_message({"type": "selection", "cells": bytearray(np.array([1, 4], "<u4").tobytes()),
                                     "source": "lasso"})
    selection = viewer.selection
    assert selection is not None and len(selection) == 2
    assert list(selection.obs_names) == ["cell-1", "cell-4"]
    adata = _adata()
    assert list(viewer.subset(adata).obs_names) == ["cell-1", "cell-4"]
    assert list(viewer.subset(adata, "Picked").obs_names) == ["cell-3"]
    assert list(viewer.subset(adata, "highlight_1").obs_names) == ["cell-0", "cell-2"]
    assert list(viewer.subset(adata, "visible").obs_names) == ["cell-0", "cell-1", "cell-2"]
    with pytest.raises(ValueError, match="names 0 highlight groups.*highlight_1 \\(Lasso selection\\)"):
        viewer.subset(adata, "nothing")
    # By identifier (the var index), or by the names shown: here those of a
    # var column (feature_name) the viewer's data mapping names the genes by.
    assert viewer.gene_sets() == {"T": ["ENSG00000198851"]}
    assert viewer.gene_sets(symbols=True) == {"T": ["CD3E"]}
    # Named by the var index, the viewer sends no symbols: the names shown
    # are the identifiers.
    viewer._handle_frontend_message(_result("geneSet", "markers:B", {
        "name": "B", "description": "Marker genes", "genes": ["ENSG00000156738"], "symbols": None, "origin": "markers",
    }, cells=[], revision=40))
    assert viewer.gene_sets()["B"] == ["ENSG00000156738"]
    assert viewer.gene_sets(symbols=True)["B"] == ["ENSG00000156738"]


def _pushed(viewer: Any, result_id: str) -> Any:
    return next(result for result in viewer.results if result.id == result_id)


def test_push_aligns_a_column_by_obs_names_and_sends_it_at_once(connected) -> None:
    viewer, fake = connected()
    adata = _adata()
    assert viewer.push(adata.obs["leiden"]) is None
    result = _pushed(viewer, "pushed:leiden")
    assert result.kind == "categorical" and result.meta["categories"] == ["0", "1", "2"]
    assert result.values.tolist() == [0, 1, 0, 2, 1, 0]
    # Some cells, in another order: the cells named carry their values.
    viewer.push(pd.Series([2.5, -1.0], index=["cell-4", "cell-1"], name="score"))
    result = _pushed(viewer, "pushed:score")
    assert result.kind == "continuous"
    assert result.cells.tolist() == [1, 4] and result.values.tolist() == [-1.0, 2.5]
    with pytest.raises(ValueError, match="not cells of the viewer's dataset"):
        viewer.push(pd.Series([1.0], index=["nope"], name="x"))
    with pytest.raises(ValueError, match="Name the field"):
        viewer.push(np.zeros(N_OBS))
    viewer.push_gene_set("T cells", ["CD3E", "CD3E", " CD3D "], "markers")
    assert viewer.gene_sets()["T cells"] == ["CD3E", "CD3D"]
    frames = [frame for message, frame in fake.sent if message["type"] == "setField"]
    assert all(isinstance(frame, bytearray) for frame in frames)


def test_a_refused_push_is_reported_where_the_notebook_looks(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    viewer = serverless_viewer()
    viewer.state.ready = {"n_cells": N_OBS, "n_obs": N_OBS, "dimensions": 2}
    viewer._obs_names_cache = pd.Index(OBS_NAMES)

    def refuse(message: dict[str, Any]) -> bool:
        viewer._handle_frontend_message(
            {"type": "command_error", "command": "setField", "reason": 'A field named "leiden" exists already'}
        )
        return True

    monkeypatch.setattr(viewer, "send_message", refuse)
    viewer.push(_adata().obs["leiden"])
    assert "viewer.push(...) did not take effect in the viewer" in capsys.readouterr().err



def test_run_all_a_push_before_the_viewer_opened_goes_when_it_says_its_cells(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # A format only the browser reads: the cell count comes with `ready`.
    viewer = serverless_viewer()
    viewer._obs_names_cache = pd.Index(OBS_NAMES)
    sent: list[dict[str, Any]] = []
    monkeypatch.setattr(viewer, "send_message", sent.append)
    errors: list[dict[str, Any]] = []
    viewer.on_command_error(errors.append)
    assert viewer.push(_adata().obs["leiden"]) is None
    viewer.push(np.arange(2), key="short")
    assert sent == []
    viewer._handle_frontend_message({"type": "ready", "n_cells": N_OBS, "n_obs": N_OBS, "dimensions": 2})
    assert [message["key"] for message in sent] == ["leiden"]
    assert [error["command"] for error in errors] == ["setField"]
    assert "2 values" in errors[0]["reason"]
    # A dataset this process reads: its cell count is known, nothing waits.
    served = serverless_viewer()
    served._served_cell_count = N_OBS
    served._obs_names_cache = pd.Index(OBS_NAMES)
    out: list[dict[str, Any]] = []
    monkeypatch.setattr(served, "send_message", out.append)
    served.push(_adata().obs["leiden"])
    assert len(out) == 1


@pytest.mark.parametrize(
    ("values", "kind", "expected"),
    [
        (pd.Series([True, False, True, True, False, False]), "categorical", [1, 0, 1, 1, 0, 0]),
        (pd.Series(["a", None, "b", "a", "b", "a"], dtype="string"), "categorical", [0, 255, 1, 0, 1, 0]),
        (np.arange(N_OBS, dtype=np.int64), "continuous", [0, 1, 2, 3, 4, 5]),
        (pd.Series([0.5, np.nan, 1, 1, 1, 1]), "continuous", [0.5, np.nan, 1, 1, 1, 1]),
    ],
)
def test_push_dtype_rules(values: Any, kind: str, expected: list[float]) -> None:
    message = field_message(values, "k", n_obs=N_OBS, obs_names=None)
    assert message["kind"] == kind
    np.testing.assert_array_equal(message["values"], np.asarray(expected, dtype=message["values"].dtype))
    wire._require_frontend_message(message | {"requestId": "r"})


def test_push_refuses_what_it_cannot_carry() -> None:
    with pytest.raises(ValueError, match="infinity"):
        field_message(np.array([1.0, np.inf, 0, 0, 0, 0]), "k", n_obs=N_OBS, obs_names=None)
    with pytest.raises(ValueError, match="5 values"):
        field_message(np.zeros(5), "k", n_obs=N_OBS, obs_names=None)


def test_highlight_cells_takes_obs_names(connected) -> None:
    viewer, fake = connected()
    viewer.highlight_cells(pd.Index(["cell-5", "cell-1"]), label="Two")
    message, frame = fake.sent[-1]
    assert message["label"] == "Two" and sorted(message["cells"].tolist()) == [1, 5]
    with pytest.raises(ValueError, match="not in this viewer's dataset"):
        viewer.highlight_cells(["cell-9"])


def test_gmt_round_trip(tmp_path: Path) -> None:
    path = cellucid.write_gmt(
        {"T cells": ["CD3E", "CD3D", "CD3E"], "B\tcells": ["MS4A1"], "B\r\ncells": ["CD79A"], "a  b": ["X"]},
        tmp_path / "sets.gmt", descriptions={"T cells": "markers"},
    )
    # Each run of tabs and line breaks one space, never other whitespace; a
    # name that then repeats one is written as read_gmt names it.
    assert path.read_text() == "T cells\tmarkers\tCD3E\tCD3D\nB cells\t\tMS4A1\nB cells (2)\t\tCD79A\na  b\t\tX\n"
    assert cellucid.read_gmt(path) == {
        "T cells": ["CD3E", "CD3D"], "B cells": ["MS4A1"], "B cells (2)": ["CD79A"], "a  b": ["X"]}


def _shared_fixture(name: str, digest: str) -> Path:
    path = Path(__file__).parent / "fixtures" / name
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, f"{name} changed: update it in all three repositories"
    return path


def test_the_shared_gmt_fixture_reads_as_the_web_app_and_r_read_it() -> None:
    path = _shared_fixture("gene_sets.gmt", "071d3e007164ff4d9c311820ba3a1be3061b7fd4b80e59142df17690c37b8bfd")
    assert cellucid.read_gmt(path) == {
        "T cells": ["CD3E", "CD3D"], "Set 2": ["X"], "T cells (2)": ["CD2"], "B": ["MS4A1"]}


def test_the_shared_obs_column_csv_fixture_keeps_every_declared_kind() -> None:
    path = _shared_fixture("obs_columns.csv", "b0da8ee441532a3ac2199e9e22f827c82002c4e63ce1b7bec5fc3d4800ec41e9")
    adata = ad.AnnData(obs=pd.DataFrame(index=["c0", "c,1", "c2", "c3", "c4", "c5"]))
    summary = cellucid.pull(adata, path)
    obs = adata.obs
    assert str(obs["cellucid_highlight__highlight_1"].dtype) == "boolean"
    assert obs["cellucid_highlight__highlight_1"].tolist() == [True, False, True, pd.NA, False, pd.NA]
    # Labels that read as numbers, or as True/False, stay labels.
    assert list(obs["leiden"].cat.categories) == ["0", "1", "10"]
    leiden = obs["leiden"].astype(object)
    assert leiden.where(leiden.notna(), None).tolist() == ["0", "1", "10", "1", None, None]
    assert list(obs["flag"].cat.categories) == ["True", "False"]
    flag = obs["flag"].astype(object)
    assert flag.where(flag.notna(), None).tolist() == ["True", "False", "True", None, "False", None]
    np.testing.assert_array_equal(obs["score"].to_numpy(), [0.5, np.nan, -1.25, 3.0, 1000.0, np.nan])
    assert {info["kind"] for info in summary.columns.values()} == {"highlight", "categorical", "continuous"}


def test_the_web_apps_exports_pull_in_one_call(tmp_path: Path) -> None:
    (tmp_path / "highlight.csv").write_text(
        "obs_name,cellucid_highlight__highlight_1\nkind,highlight\n" + "".join(
            f"cell-{index},{'true' if index in {0, 2} else 'false'}\n" for index in range(N_OBS)))
    (tmp_path / "annotation.csv").write_text("obs_name,annotation\nkind,categorical\ncell-0,T\ncell-1,B\ncell-5,T\n")
    (tmp_path / "score.csv").write_text("obs_name,T score\nkind,continuous\ncell-3,1.5\ncell-1,\n")
    cellucid.write_gmt({"T": ["CD3E"]}, tmp_path / "markers.gmt")
    adata = _adata()
    summary = cellucid.pull(adata, tmp_path)
    obs = adata.obs
    assert obs["cellucid_highlight__highlight_1"].tolist() == [True, False, True, False, False, False]
    annotation = obs["annotation"].astype(object)
    assert annotation.where(annotation.notna(), None).tolist() == ["T", "B", None, None, None, "T"]
    np.testing.assert_array_equal(obs["T score"].to_numpy(), [np.nan, np.nan, np.nan, 1.5, np.nan, np.nan])
    assert summary.gene_sets == {"T": ["CD3E"]}
    # Again: the columns and gene sets it wrote are replaced, not conflicts.
    again = cellucid.pull(adata, tmp_path, on_conflict="suffix")
    assert set(again.columns) == set(summary.columns)
    assert list(adata.uns["cellucid"]["gene_sets"]) == ["T"]
    (tmp_path / "bad.csv").write_text("obs_name,x\nkind,continuous\ncell-99,1\n")
    with pytest.raises(ValueError, match="not in this object"):
        cellucid.pull(adata, tmp_path / "bad.csv")
    (tmp_path / "bad.csv").write_text("obs_name,x\ncell-1,1\n")
    with pytest.raises(ValueError, match="kind row"):
        cellucid.pull(adata, tmp_path / "bad.csv")
    (tmp_path / "bad.csv").write_text("obs_name,x\nkind,continuous\ncell-1,one\n")
    with pytest.raises(ValueError, match="continuous column 'x' holds 'one'"):
        cellucid.pull(adata, tmp_path / "bad.csv")
    (tmp_path / "bad.csv").write_text("obs_name,x\nkind,label\ncell-1,one\n")
    with pytest.raises(ValueError, match="of kind 'label'"):
        cellucid.pull(adata, tmp_path / "bad.csv")


def test_gene_sets_of_one_name_are_named_by_one_rule() -> None:
    from cellucid.anndata_session._results import ViewerResult, _gene_sets_of

    def gene_set(result_id: str, origin: str) -> ViewerResult:
        return ViewerResult(kind="geneSet", id=result_id, revision=1, n_obs=1,
                            meta={"name": "T", "description": "", "genes": [result_id], "origin": origin},
                            cells=np.empty(0, dtype=np.uint32), values=None)

    named = _gene_sets_of([gene_set("a", "markers"), gene_set("b", "markers"), gene_set("c", "notebook")])
    assert [name for _, name, _ in named] == ["T", "T (2)", "T (3)"]


def test_a_bool_column_is_pushed_with_the_labels_r_pushes() -> None:
    message = field_message(pd.Series([True, False, None], dtype="boolean"), "flag", n_obs=3, obs_names=None)
    assert message["categories"] == ["false", "true"]
    assert message["values"].tolist() == [1, 0, 255]
    message = field_message(np.array([True, True]), "flag", n_obs=2, obs_names=None)
    assert message["categories"] == ["true"]


def test_pull_reads_obs_names_of_a_file_alone(tmp_path: Path) -> None:
    path = tmp_path / "data.h5ad"
    _adata().write_h5ad(path)
    assert list(read_obs_names(path)) == OBS_NAMES
    assert list(read_obs_names(_adata())) == OBS_NAMES


def test_mudata_takes_the_pull_on_its_global_obs(connected) -> None:
    md = pytest.importorskip("mudata")
    viewer, _fake = connected()
    mdata = md.MuData({"rna": _adata(), "atac": _adata()}, axis=-1)
    if int(mdata.n_obs) != N_OBS:
        pytest.skip("the two modalities share their cells")
    viewer.pull(mdata)
    assert "cellucid_highlight__highlight_1" in mdata.obs


def test_the_protocol_document_lists_every_result_kind() -> None:
    fixture = json.loads((Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text())
    assert set(fixture["results"]["kinds"]) == set(wire._RESULT_META_FIELDS)


def _prepared(tmp_path: Path, file_format: str) -> tuple[Path, ad.AnnData]:
    """A Cellucid-ready file of 40 cells (priority order) and the object it was written from."""
    n_obs = 40
    rng = np.random.default_rng(3)
    adata = ad.AnnData(
        X=rng.poisson(1.0, size=(n_obs, 3)).astype(np.float32),
        obs=pd.DataFrame({"group": pd.Categorical(rng.choice(["a", "b"], n_obs))},
                         index=[f"cell-{index}" for index in range(n_obs)]),
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    adata.obsm["X_umap"] = rng.normal(size=(n_obs, 2)).astype(np.float32)
    path = tmp_path / f"prepared.{file_format}"
    cellucid.prepare(path, format=file_format, data=adata, dataset_id="prepared", dataset_name="Prepared")
    return path, adata


@pytest.mark.parametrize("file_format", ["h5ad", "zarr", "zarr.zip"])
def test_a_priority_ordered_files_obs_names_are_read_in_original_order(tmp_path: Path, file_format: str) -> None:
    # The viewer names a cell by its original position, whatever order the
    # file stores the cells in: names read from the file follow that order.
    path, adata = _prepared(tmp_path, file_format)
    if file_format != "zarr.zip":
        stored = ad.read_h5ad(path) if file_format == "h5ad" else ad.read_zarr(path)
        assert list(stored.obs_names) != list(adata.obs_names), "prepare() stores cells in priority order"
    assert list(read_obs_names(path)) == list(adata.obs_names)


def test_a_viewer_of_a_priority_ordered_file_names_cells_by_original_position(
    connected, tmp_path: Path
) -> None:
    path, adata = _prepared(tmp_path, "zarr")
    viewer, fake = connected([])
    viewer.data = str(path)
    del viewer._obs_names_cache  # read from the file, as a real viewer of it does
    viewer.state.ready = {"n_cells": 40, "n_obs": 40, "dimensions": 2}
    viewer._handle_frontend_message({"type": "selection", "cells": bytearray(np.array([3, 17], "<u4").tobytes()),
                                     "source": "lasso"})
    selection = viewer.selection
    assert selection is not None and list(selection.obs_names) == ["cell-3", "cell-17"]
    assert viewer._cell_positions(pd.Index(["cell-17", "cell-3"])).tolist() == [17, 3]
    # A column pushed by obs names lands on the cells of those names.
    message = field_message(pd.Series([1.0, 2.0], index=["cell-17", "cell-3"]), "score", n_obs=40,
                            obs_names=viewer.obs_names)
    assert message["cells"].tolist() == [3, 17] and message["values"].tolist() == [2.0, 1.0]
    # The object restored to its original order is what pull and subset take;
    # the file read as stored is refused, never written in the wrong rows.
    stored = ad.read_zarr(path)
    with pytest.raises(ValueError, match="restore_original_order"):
        viewer.subset(stored)
    with pytest.raises(ValueError, match="restore_original_order"):
        viewer.pull(stored)
    with pytest.raises(ValueError, match="restore_original_order"):
        cellucid.pull(stored, viewer)
    restored = cellucid.restore_original_order(stored)
    assert list(viewer.subset(restored).obs_names) == ["cell-3", "cell-17"]


def test_the_object_a_viewer_shows_is_pulled_into_in_its_own_order(connected) -> None:
    # An object is served in its own row order, an earlier file's position
    # column or not: pulling into that very object is in order.
    viewer, _fake = connected()
    adata = _adata()
    adata.obs["cellucid_original_position"] = np.arange(N_OBS, dtype=np.uint32)[::-1]
    viewer.data = adata
    viewer.pull(adata)
    assert adata.obs["cellucid_highlight__highlight_1"].tolist()[:3] == [True, False, True]
    # The same cells in another order (sorted since, say) are refused, never
    # written in the wrong rows; another dataset's cells too.
    with pytest.raises(ValueError, match="another order.*adata\\[viewer.obs_names\\]"):
        viewer.pull(_adata()[::-1].copy())
    with pytest.raises(ValueError, match="another order"):
        viewer.subset(_adata()[::-1].copy(), "highlight_1")
    renamed = _adata()
    renamed.obs_names = [f"other-{index}" for index in range(N_OBS)]
    with pytest.raises(ValueError, match="6 of this object's cells are not cells of the viewer's dataset"):
        viewer.pull(renamed)
