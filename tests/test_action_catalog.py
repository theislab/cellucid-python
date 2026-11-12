"""The action catalog and the viewer methods generated from it.

``src/cellucid/action_catalog.json`` is the web app's action catalog
document byte for byte (written by the web repository's
``node scripts/generate-action-catalog.mjs``, which also writes cellucid-r's
copy); the web app (``tests/action-catalog-document-contract.test.mjs``) and
cellucid-r (``tests/testthat/test-actions.R``) pin the same SHA-256.
``src/cellucid/jupyter/_actions_generated.py`` is one viewer method per
action, written from it by ``python scripts/generate_actions.py``; every
method sends the action through ``AnnDataViewer.act``, which checks and
converts parameters by the catalog, the same way for every action.
"""

from __future__ import annotations

import hashlib
import importlib.util
import inspect
import json
import keyword
import struct
import threading
import warnings
from importlib import resources
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest
from _serverless_viewer import serverless_viewer
from _web_root import rendered_by_web

import cellucid.jupyter as jupyter
import cellucid.jupyter._actions_generated as generated
import cellucid.jupyter._wire as wire
from cellucid.jupyter import (
    OMIT,
    ActionError,
    ActionNotRecordedWarning,
    AnnDataViewer,
    action_catalog,
)

ROOT = Path(__file__).resolve().parents[1]

#: The digest of ``action_catalog.json``; the web app and cellucid-r pin the same one.
ACTION_CATALOG_SHA256 = "725d32e5357cebbcff332be5a599667cf42ebb63ce55bd73e39ce4a1469556c4"

REGENERATE = (
    "regenerate it: node scripts/generate-action-catalog.mjs in the web repository (it writes the three "
    "copies and prints the digest), then python scripts/generate_actions.py here, and pin the digest"
)


def _generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("generate_actions", ROOT / "scripts" / "generate_actions.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _catalog_text() -> str:
    return resources.files("cellucid").joinpath("action_catalog.json").read_text(encoding="utf-8")


def test_the_packaged_catalogue_is_the_shared_document() -> None:
    data = resources.files("cellucid").joinpath("action_catalog.json").read_bytes()
    assert hashlib.sha256(data).hexdigest() == ACTION_CATALOG_SHA256, REGENERATE
    document = json.loads(data)
    assert document["schema"] == 1
    assert set(action_catalog()) == {entry["type"] for entry in document["actions"]}


def test_the_generated_methods_are_what_the_catalogue_makes_now() -> None:
    expected = _generator().render(_catalog_text())
    current = (ROOT / "src" / "cellucid" / "jupyter" / "_actions_generated.py").read_text(encoding="utf-8")
    assert current == expected, "src/cellucid/jupyter/_actions_generated.py is stale: run python scripts/generate_actions.py"
    assert _generator().main(["--check"]) == 0


def _keyword(name: str) -> str:
    return f"{name}_" if keyword.iskeyword(name) else name


def test_every_action_is_one_method_with_the_actions_parameters() -> None:
    methods = {name for name, value in vars(generated.ViewerActions).items() if inspect.isfunction(value)}
    tools = {entry["tool"]: entry for entry in action_catalog().values()}
    assert methods == set(tools)
    for tool, entry in tools.items():
        signature = inspect.signature(getattr(AnnDataViewer, tool))
        parameters = list(signature.parameters.values())[1:]
        assert all(parameter.kind is inspect.Parameter.KEYWORD_ONLY for parameter in parameters), tool
        schema = entry["params"]
        required = set(schema.get("required", []))
        assert [parameter.name for parameter in parameters] == [
            _keyword(name) for name in [*schema.get("required", []), *(n for n in schema.get("properties", {}) if n not in required)]
        ], tool
        for parameter in parameters:
            name = parameter.name.removesuffix("_") if parameter.name.endswith("_") else parameter.name
            if name in required:
                assert parameter.default is inspect.Parameter.empty, f"{tool}({parameter.name})"
            else:
                default = schema["properties"][name].get("default", OMIT)
                if not isinstance(default, bool | int | float | str | type(None)):
                    default = OMIT
                assert parameter.default == default or parameter.default is default, f"{tool}({parameter.name})"
        doc = getattr(AnnDataViewer, tool).__doc__
        assert doc is not None and f"``{entry['type']}``" in doc


def test_no_generated_method_hides_or_is_hidden_by_another_attribute() -> None:
    import anywidget

    from cellucid.jupyter._widget import Widget

    tools = {entry["tool"] for entry in action_catalog().values()}
    others: set[str] = set()
    for owner in (AnnDataViewer, Widget):
        for base in owner.__mro__:
            if base is generated.ViewerActions:
                continue
            others |= set(vars(base))
    others |= set(dir(anywidget.AnyWidget))
    assert tools & others == set()
    # The widget's methods are the generated ones.
    assert Widget.color_by is generated.ViewerActions.color_by


class _Sent:
    """What a viewer sent: each command, and the frame it became."""

    def __init__(self, viewer: AnnDataViewer, monkeypatch: pytest.MonkeyPatch) -> None:
        self.messages: list[dict[str, Any]] = []
        self.frames: list[str | bytearray] = []

        def send(message: dict[str, Any]) -> bool:
            self.messages.append(message)
            self.frames.append(wire._require_frontend_message(message))
            return True

        monkeypatch.setattr(viewer, "send_message", send)


@pytest.fixture
def viewer(monkeypatch: pytest.MonkeyPatch) -> tuple[AnnDataViewer, _Sent]:
    shown = serverless_viewer()
    shown._obs_names_cache = pd.Index([f"cell-{index}" for index in range(6)])
    shown.state.ready = {"n_cells": 6, "n_obs": 6, "dimensions": 2}
    return shown, _Sent(shown, monkeypatch)


def _decode(frame: bytes | bytearray) -> tuple[dict[str, Any], list[list[float]]]:
    (length,) = struct.unpack("<I", frame[:4])
    header = json.loads(bytes(frame[4 : 4 + length]))
    body = bytes(frame[4 + length :])
    arrays = []
    offset = 0
    for entry in header["arrays"]:
        dtype = wire._VALUE_DTYPES[entry["type"]]
        arrays.append(np.frombuffer(body, dtype=dtype, count=entry["length"], offset=offset).tolist())
        offset += entry["length"] * dtype.itemsize
    assert offset == len(body)
    return header, arrays


def test_cells_are_rows_a_mask_or_obs_names_in_the_order_given(viewer: tuple[AnnDataViewer, _Sent]) -> None:
    shown, sent = viewer
    shown.highlight_cells(cells=[4, 1])
    shown.highlight_cells(cells=np.array([5, 0, 3], dtype=np.int64))
    shown.highlight_cells(cells=[False, True, False, True, False, False])
    shown.highlight_cells(cells=np.array([True, False, False, False, False, True]))
    shown.highlight_cells(cells=pd.Index(["cell-2", "cell-0"]))
    shown.highlight_cells(cells=pd.Series([True, False, True], index=["cell-5", "cell-1", "cell-3"]))
    assert [message["params"]["cells"].tolist() for message in sent.messages] == [
        [4, 1], [5, 0, 3], [1, 3], [0, 5], [2, 0], [5, 3],
    ]
    assert {message["params"]["cells"].dtype.str for message in sent.messages} == {"<u4"}
    # A binary frame, its cells uint32 in the order given.
    header, arrays = _decode(sent.frames[1])
    assert header["action"] == "highlight.cells" and header["params"] == {"tool": "annotation"}
    assert header["arrays"] == [{"path": ["cells"], "type": "uint32", "length": 3}] and arrays == [[5, 0, 3]]
    for cells, match in [
        ([1, 1], "more than once"),
        ([-1], "non-negative"),
        ([6], "beyond the dataset's 6 cells"),
        ([True, False], "mask of 2 cells"),
        (["cell-9"], "not in the viewer's dataset"),
        ([1, "cell-2"], "not a mix"),
        ([], "at least 1 cell"),
    ]:
        with pytest.raises((TypeError, ValueError), match=match):
            shown.highlight_cells(cells=cells)
    assert len(sent.messages) == 6


def test_values_travel_as_their_value_type(viewer: tuple[AnnDataViewer, _Sent]) -> None:
    shown, sent = viewer
    for values, dtype in [
        (np.array([0, 1, 255, 2, 0, 1], dtype=np.uint8), "uint8"),
        (np.array([True, False, True, True, False, False]), "uint8"),
        (np.arange(6, dtype=np.int64), "float64"),
        (np.ones(6, dtype=np.float16), "float32"),
        ([0.5, 1.0, 2.0, 3.0, 4.0, 5.0], "float64"),
    ]:
        shown.field_push(key="k", kind="continuous", categories=None, n_obs=6, cells=None, values=values)
        assert sent.messages[-1]["params"]["values"].dtype.name == dtype
    header, arrays = _decode(sent.frames[0])
    assert header["params"] == {"key": "k", "kind": "continuous", "categories": None, "n_obs": 6, "cells": None}
    assert header["arrays"] == [{"path": ["values"], "type": "uint8", "length": 6}]
    assert arrays == [[0, 1, 255, 2, 0, 1]]


def test_parameters_are_checked_by_the_catalogue_before_anything_is_sent(viewer: tuple[AnnDataViewer, _Sent]) -> None:
    shown, sent = viewer
    request = shown.act("history.get")
    assert sent.messages[-1]["params"] == {}
    assert isinstance(sent.frames[-1], str) and json.loads(sent.frames[-1]) == {
        "type": "action", "requestId": request.request_id, "action": "history.get", "params": {},
        "turn": None, "turnLabel": None,
    }
    # OMIT leaves a parameter out; the generated default is sent as given.
    shown.history_get()
    assert sent.messages[-1]["params"] == {"all_branches": False, "limit": 200}
    shown.act("highlight.cells", cells=[1], label=OMIT)
    assert set(sent.messages[-1]["params"]) == {"cells"}
    # A keyword in Python is written with an underscore and sent as written in the catalogue.
    shown.analysis_markers_get(in_=2)
    assert sent.messages[-1]["params"] == {"in": 2}
    # NumPy values and pandas indexes become JSON.
    shown.gene_set_save(name="T", genes=pd.Index(["CD3E", "CD3D"]), description="markers")
    assert sent.messages[-1]["params"] == {"name": "T", "genes": ["CD3E", "CD3D"], "description": "markers"}
    shown.camera_reset(view=np.str_("live"))
    assert sent.messages[-1]["params"] == {"view": "live"}
    count = len(sent.messages)
    with pytest.raises(TypeError, match=r"has no parameter 'colour'; its parameters are: cells, color"):
        shown.act("highlight.cells", cells=[1], colour="#ff0000")
    with pytest.raises(TypeError, match="needs 'field'"):
        shown.act("color.by")
    with pytest.raises(ValueError, match="There is no action 'colour.by'; did you mean 'color.by'"):
        shown.act("colour.by", field=None)
    with pytest.raises(TypeError, match="dotted type"):
        shown.act(7)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="JSON-like"):
        shown.camera_reset(view=object())  # type: ignore[arg-type]
    assert len(sent.messages) == count


def _answer(
    shown: AnnDataViewer, request_id: str, outcome: Any, step: dict[str, str] | None, unrecorded: str | None = None
) -> None:
    shown._handle_frontend_message({
        "type": "action_result", "requestId": request_id, "step": step, "unrecorded": unrecorded,
        "body": bytearray(json.dumps(outcome).encode("utf-8")),
    })


def test_a_request_is_answered_once_with_its_outcome_or_its_refusal(viewer: tuple[AnnDataViewer, _Sent]) -> None:
    shown, _sent = viewer
    request = shown.color_by(field={"source": "obs", "key": "cell_type"})
    assert not request.done and request.step is None
    assert repr(request) == "ActionRequest('color.by', waiting)"
    _answer(shown, request.request_id, {"shown": "cell_type"}, {"id": "s3", "label": "Colored by cell_type"})
    assert request.done and request.result(timeout=0) == {"shown": "cell_type"}
    assert request.step == {"id": "s3", "label": "Colored by cell_type"}
    # The outcome handed out is a copy.
    request.result()["shown"] = "changed"
    assert request.result() == {"shown": "cell_type"}

    refused = shown.color_by(field={"source": "obs", "key": "nope"})
    shown._handle_frontend_message({
        "type": "command_error", "command": "action", "requestId": refused.request_id,
        "reason": 'There is no field "nope" in obs',
    })
    with pytest.raises(ActionError, match='did not apply color.by: There is no field "nope" in obs') as error:
        refused.result()
    assert error.value.reason == 'There is no field "nope" in obs' and error.value.request_id == refused.request_id
    assert refused.step is None
    # Answered requests are forgotten: fire-and-forget actions do not pile up.
    assert shown._action_waiters == {} and shown._action_types == {}

    waiting = shown.history_get()
    with pytest.raises(TimeoutError, match="viewer's answer to history.get"):
        waiting.result(timeout=0.01)
    waiting.cancel()
    assert shown._action_waiters == {}
    with pytest.raises(RuntimeError, match="canceled"):
        waiting.result()
    # An answer to a request nobody awaits any more is dropped.
    _answer(shown, waiting.request_id, None, None)
    assert repr(waiting) == "ActionRequest('history.get', canceled)"


def test_a_change_that_made_no_step_says_why_and_warns_when_collected(viewer: tuple[AnnDataViewer, _Sent]) -> None:
    """The viewer says why a change made no step (``unrecorded``): the request
    keeps it, and collecting the outcome warns with the sentence the shared
    protocol fixture gives for it; a step, or a read, warns nothing."""
    shown, _sent = viewer
    sentences = json.loads((ROOT / "tests" / "fixtures" / "viewer_protocol.json").read_text(encoding="utf-8"))["action"][
        "unrecorded"
    ]
    assert sentences == wire._ACTION_UNRECORDED
    for reason, sentence in sentences.items():
        request = shown.color_by(field={"source": "obs", "key": "cell_type"})
        assert request.unrecorded is None, "not answered yet"
        _answer(shown, request.request_id, {"field": "cell_type"}, None, reason)
        assert request.step is None and request.unrecorded == reason
        with pytest.warns(ActionNotRecordedWarning) as caught:
            assert request.result(timeout=0) == {"field": "cell_type"}
        assert [str(warning.message) for warning in caught] == [f"viewer.color_by(...): {sentence}"]
        notice = caught[0].message
        assert isinstance(notice, ActionNotRecordedWarning)
        assert (notice.action, notice.request_id, notice.unrecorded) == ("color.by", request.request_id, reason)
        assert caught[0].filename == __file__, "the warning names the line that collected it"
    assert issubclass(ActionNotRecordedWarning, UserWarning)

    made = shown.color_by(field={"source": "obs", "key": "cell_type"})
    _answer(shown, made.request_id, None, {"id": "s3", "label": "Colored by cell_type"})
    read = shown.fields_list()
    _answer(shown, read.request_id, {"fields": []}, None)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert made.result(timeout=0) is None and made.unrecorded is None
        assert read.result(timeout=0) == {"fields": []} and read.unrecorded is None


def test_why_no_step_is_one_of_the_closed_set_and_never_beside_a_step() -> None:
    def event(step: dict[str, str] | None, unrecorded: object) -> dict[str, Any]:
        return {"type": "action_result", "requestId": "r1", "step": step, "unrecorded": unrecorded, "body": bytearray(b"null")}

    for reason in wire._ACTION_UNRECORDED:
        _kind, payload = wire._require_inbound_jupyter_event(event(None, reason))
        assert payload["unrecorded"] == reason
    with pytest.raises(ValueError, match="unrecorded must be one of unchanged, retracted"):
        wire._require_inbound_jupyter_event(event(None, "forgotten"))
    with pytest.raises(ValueError, match="unrecorded must be None when the action made a step"):
        wire._require_inbound_jupyter_event(event({"id": "s3", "label": "Colored by leiden"}, "unchanged"))
    with pytest.raises(ValueError, match="unrecorded must be one of"):
        wire._require_inbound_jupyter_event(event(None, 3))


def test_another_thread_waits_for_the_answer_and_a_stopped_viewer_ends_the_wait(
    viewer: tuple[AnnDataViewer, _Sent],
) -> None:
    shown, _sent = viewer
    request = shown.fields_list()
    caught: list[Any] = []
    worker = threading.Thread(target=lambda: caught.append(request.result(timeout=10)))
    worker.start()
    _answer(shown, request.request_id, {"fields": []}, None)
    worker.join(10)
    assert caught == [{"fields": []}]

    pending = shown.fields_list()
    failures: list[BaseException] = []

    def wait() -> None:
        try:
            pending.result(timeout=10)
        except BaseException as error:
            failures.append(error)

    worker = threading.Thread(target=wait)
    worker.start()
    with shown._event_cv:
        shown._stopped_reason = "the viewer was stopped"
        shown._event_cv.notify_all()
    worker.join(10)
    assert len(failures) == 1 and isinstance(failures[0], RuntimeError)
    assert "Stopped waiting for the viewer's answer to fields.list" in str(failures[0])


def test_the_catalogue_is_listed_for_a_notebook() -> None:
    catalog = action_catalog()
    assert catalog["color.by"]["tool"] == "color_by"
    assert callable(getattr(jupyter.AnnDataViewer, catalog["color.by"]["tool"]))


def test_cells_by_row_never_read_the_obs_names() -> None:
    """Rows and masks name cells without the viewer's obs names: those are
    read (from the served file) only for cells named by name."""
    from cellucid.jupyter._actions import cell_positions

    def unreadable() -> None:
        raise AssertionError("the obs names were read for rows")

    rows = cell_positions(np.array([5, 0, 2]), label="cells", n_obs=6, obs_names=unreadable)
    assert rows.tolist() == [5, 0, 2] and rows.dtype == np.uint32
    mask = cell_positions([True, False, True], label="cells", n_obs=3, obs_names=unreadable)
    assert mask.tolist() == [0, 2]
    named = cell_positions(["c2", "c0"], label="cells", n_obs=3, obs_names=lambda: pd.Index(["c0", "c1", "c2"]))
    assert named.tolist() == [2, 0]


def test_the_catalog_is_the_one_the_web_apps_definitions_make_now() -> None:
    """Against the action definitions themselves (``renderActionCatalog``),
    not only a pinned digest: a definition changed without regenerating the
    copies fails here, as it fails the web suite."""
    shipped = resources.files("cellucid").joinpath("action_catalog.json").read_text(encoding="utf-8")
    assert rendered_by_web("scripts/generate-action-catalog.mjs", "renderActionCatalog") == shipped, (
        "src/cellucid/action_catalog.json is not what the web app's definitions make: run "
        "node scripts/generate-action-catalog.mjs in the web app, then python scripts/generate_actions.py here."
    )
