"""The data-mapping half of the Jupyter wire: ``setDataMapping`` and ``data_mapping``.

The viewer reports its data mapping with a ``data_mapping`` event, and a
notebook changes it with the ``setDataMapping`` command, which -- unlike the
one-way commands -- is answered: by a ``data_mapping`` event of origin
``notebook`` naming the request, or by a ``command_error`` naming the command.
Both directions carry a complete DataMapping, validated by the same rule
sessions use, so a mapping the viewer would refuse never crosses the wire.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import pytest
from _serverless_viewer import serverless_viewer

import cellucid.jupyter as jupyter
import cellucid.jupyter._viewer as jupyter_viewer
import cellucid.jupyter._wire as wire

CONFORMANCE = Path(__file__).parent / "fixtures" / "data_mapping_conformance.json"

ANNDATA_MAPPING: dict[str, Any] = {
    "kind": "anndata",
    "embeddings": {"enabled": [{"key": "X_umap", "columns": [0, 1]}], "defaults": {"1d": None, "2d": {"key": "X_umap", "columns": [0, 1]}, "3d": None}},
    "knn": None,
    "geneNames": "index",
    "expression": "X",
    "log1p": {},
    "highlyVariable": None,
    "latent": None,
    "vectorFields": [],
}


def _invalid_mapping_cases() -> list[dict[str, Any]]:
    cases = json.loads(CONFORMANCE.read_text(encoding="utf-8"))["dataMapping"]
    return [case for case in cases if "error" in case["expected"]]


def _event(viewer_id: str, **overrides: Any) -> dict[str, Any]:
    event: dict[str, Any] = {
        "type": "data_mapping",
        "requestId": None,
        "origin": "default",
        "mapping": copy.deepcopy(ANNDATA_MAPPING),
    }
    event.update(overrides)
    return event


# ---------------------------------------------------------------------------
# Declarations
# ---------------------------------------------------------------------------


def test_the_wire_declares_exactly_these_field_sets() -> None:
    assert wire._COMMAND_MESSAGE_FIELDS["setDataMapping"] == frozenset(
        {"type", "requestId", "mapping"}
    )
    assert wire._INBOUND_EVENT_FIELDS["data_mapping"] == frozenset(
        {"type", "requestId", "origin", "mapping"}
    )
    origins = wire._DATA_MAPPING_ORIGINS
    assert origins == {"default", "best_guess", "session", "user", "notebook"}


def test_the_capability_document_lists_both_halves() -> None:
    document = wire._protocol_capability_document()
    assert "setDataMapping" in document["commands"]
    assert "data_mapping" in document["events"]
    assert jupyter_viewer._COMMAND_PYTHON_CALLERS["setDataMapping"] == "viewer.set_data_mapping(...)"


# ---------------------------------------------------------------------------
# Outbound: setDataMapping
# ---------------------------------------------------------------------------


def test_set_data_mapping_serializes_the_canonical_mapping() -> None:
    reordered = dict(reversed(list(ANNDATA_MAPPING.items())))
    serialized = wire._require_frontend_message(
        {"type": "setDataMapping", "requestId": "request-1", "mapping": reordered}
    )
    decoded = json.loads(serialized)
    assert decoded == {"type": "setDataMapping", "requestId": "request-1", "mapping": ANNDATA_MAPPING}
    assert list(decoded["mapping"]) == list(ANNDATA_MAPPING)


@pytest.mark.parametrize(
    ("message", "error"),
    [
        ({"type": "setDataMapping", "mapping": ANNDATA_MAPPING}, TypeError),
        ({"type": "setDataMapping", "requestId": "r", "mapping": ANNDATA_MAPPING, "x": 1}, TypeError),
        ({"type": "setDataMapping", "requestId": "", "mapping": ANNDATA_MAPPING}, TypeError),
        ({"type": "setDataMapping", "requestId": 7, "mapping": ANNDATA_MAPPING}, TypeError),
        ({"type": "setDataMapping", "requestId": "r", "mapping": None}, ValueError),
    ],
)
def test_set_data_mapping_rejects_malformed_commands(message: dict[str, Any], error: type) -> None:
    with pytest.raises(error):
        wire._require_frontend_message(message)


@pytest.mark.parametrize("case", _invalid_mapping_cases(), ids=lambda case: case["name"])
def test_set_data_mapping_refuses_every_invalid_conformance_mapping(case: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        wire._require_frontend_message(
            {"type": "setDataMapping", "requestId": "r", "mapping": case["value"]}
        )


# ---------------------------------------------------------------------------
# Inbound: data_mapping
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("origin", ["default", "best_guess", "session", "user"])
def test_a_reported_mapping_carries_no_request_id(origin: str) -> None:
    event_type, payload = wire._require_inbound_jupyter_event(
        _event("v", origin=origin))
    assert event_type == "data_mapping"
    assert payload == {"requestId": None, "origin": origin, "mapping": ANNDATA_MAPPING}


def test_a_notebook_answer_names_its_request() -> None:
    _event_type, payload = wire._require_inbound_jupyter_event(
        _event("v", origin="notebook", requestId="request-1"))
    assert payload["requestId"] == "request-1"


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"origin": "notebook", "requestId": None}, "requestId"),
        ({"origin": "user", "requestId": "request-1"}, "requestId must be null"),
        ({"origin": "unknown"}, "origin must be one of"),
        ({"origin": None}, "origin must be one of"),
        ({"origin": "notebook", "requestId": ""}, "requestId"),
        ({"mapping": {"kind": "h5ad"}}, "kind"),
        ({"mapping": {**ANNDATA_MAPPING, "extra": True}}, "unknown extra"),
    ],
)
def test_a_malformed_report_is_refused(overrides: dict[str, Any], match: str) -> None:
    with pytest.raises((TypeError, ValueError), match=match):
        wire._require_inbound_jupyter_event(_event("v", **overrides))


@pytest.mark.parametrize("case", _invalid_mapping_cases(), ids=lambda case: case["name"])
def test_a_report_with_an_invalid_conformance_mapping_is_refused(case: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        wire._require_inbound_jupyter_event(
            _event("v", mapping=case["value"]))


def test_a_report_missing_a_field_is_refused() -> None:
    event = _event("v")
    del event["origin"]
    with pytest.raises(ValueError, match="missing origin"):
        wire._require_inbound_jupyter_event(event)


def test_command_error_may_name_set_data_mapping() -> None:
    _event_type, payload = wire._require_inbound_jupyter_event(
        {
            "type": "command_error",
            "command": "setDataMapping",
            "reason": "obsp['distances'] is not symmetric",
        },
    )
    assert payload["command"] == "setDataMapping"


# ---------------------------------------------------------------------------
# The viewer: data_mapping, data_mapping_origin, on_data_mapping
# ---------------------------------------------------------------------------


def test_the_viewer_holds_the_latest_reported_mapping() -> None:
    viewer = serverless_viewer()
    assert viewer.data_mapping is None
    assert viewer.data_mapping_origin is None
    hooked: list[dict[str, Any]] = []
    viewer.on_data_mapping(hooked.append)

    viewer._handle_frontend_message(_event(viewer._viewer_id, origin="best_guess"))

    assert viewer.data_mapping == ANNDATA_MAPPING
    assert viewer.data_mapping_origin == "best_guess"
    assert hooked == [{"requestId": None, "origin": "best_guess", "mapping": ANNDATA_MAPPING}]
    assert viewer.state.data_mapping == hooked[0]

    # The property is a copy: editing it does not change what the viewer holds.
    copied = viewer.data_mapping
    assert copied is not None
    copied["embeddings"]["defaults"]["2d"] = None
    assert viewer.data_mapping == ANNDATA_MAPPING


# ---------------------------------------------------------------------------
# The viewer: set_data_mapping
# ---------------------------------------------------------------------------


def _viewer_with_mapping(mapping: dict[str, Any]) -> jupyter.AnnDataViewer:
    viewer = serverless_viewer()
    viewer._handle_frontend_message(_event(viewer._viewer_id, mapping=copy.deepcopy(mapping)))
    return viewer


def test_set_data_mapping_before_the_viewer_reported_one_goes_when_it_does(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # "Run All": the cell runs before the viewer opened its dataset.
    viewer = serverless_viewer()
    sent: list[dict[str, Any]] = []
    monkeypatch.setattr(viewer, "send_message", sent.append)
    assert viewer.set_data_mapping(knn="connectivities") is None
    assert sent == []
    viewer._handle_frontend_message(_event(viewer._viewer_id, mapping=copy.deepcopy(ANNDATA_MAPPING)))
    (command,) = sent
    assert command["type"] == "setDataMapping"
    assert command["mapping"] == {**ANNDATA_MAPPING, "knn": "connectivities"}


def test_a_waiting_mapping_the_viewer_cannot_take_is_reported_as_a_command_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = serverless_viewer()
    errors: list[dict[str, Any]] = []
    viewer.on_command_error(errors.append)
    sent: list[dict[str, Any]] = []
    monkeypatch.setattr(viewer, "send_message", sent.append)
    viewer.set_data_mapping(expression="layers:x")
    viewer._handle_frontend_message(_event(viewer._viewer_id, mapping=copy.deepcopy(ANNDATA_MAPPING)))
    assert sent == []
    assert [error["command"] for error in errors] == ["setDataMapping"]
    assert "expression" in errors[0]["reason"]


@pytest.mark.parametrize(
    ("mapping", "fields", "error", "match"),
    [
        (ANNDATA_MAPPING, {"kind": "export"}, TypeError, "cannot change the mapping kind"),
        (ANNDATA_MAPPING, {"embedding": {}}, TypeError, "unknown field"),
        (ANNDATA_MAPPING, {"expression": "layers:x"}, ValueError, "expression"),
        (
            ANNDATA_MAPPING,
            {"embeddings": {"enabled": [{"key": "X_pca", "columns": [1, 1]}], "defaults": {"1d": None, "2d": {"key": "X_pca", "columns": [1, 1]}, "3d": None}}},
            ValueError,
            "different column",
        ),
        (
            ANNDATA_MAPPING,
            {"embeddings": {"enabled": [], "defaults": {"1d": None, "2d": None, "3d": None}}},
            ValueError,
            "at least one",
        ),
    ],
)
def test_set_data_mapping_validates_before_sending(
    monkeypatch: pytest.MonkeyPatch,
    mapping: dict[str, Any],
    fields: dict[str, Any],
    error: type,
    match: str,
) -> None:
    viewer = _viewer_with_mapping(mapping)
    sent: list[dict[str, Any]] = []
    monkeypatch.setattr(viewer, "send_message", sent.append)
    with pytest.raises(error, match=match):
        viewer.set_data_mapping(**fields)
    assert sent == []


def test_set_data_mapping_sends_the_merged_mapping_and_its_answer_is_the_live_mapping(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = _viewer_with_mapping(ANNDATA_MAPPING)
    sent: list[str] = []
    monkeypatch.setattr(viewer, "send_message", lambda message: sent.append(wire._require_frontend_message(message)))
    returned = viewer.set_data_mapping(
        embeddings={"enabled": [{"key": "X_pca", "columns": [3, 1]}], "defaults": {"1d": None, "2d": {"key": "X_pca", "columns": [3, 1]}, "3d": None}},
        knn="connectivities",
    )
    assert returned is None
    expected = {
        **ANNDATA_MAPPING,
        "embeddings": {"enabled": [{"key": "X_pca", "columns": [3, 1]}], "defaults": {"1d": None, "2d": {"key": "X_pca", "columns": [3, 1]}, "3d": None}},
        "knn": "connectivities",
    }
    (serialized,) = sent
    command = json.loads(serialized)
    assert set(command) == {"type", "requestId", "mapping"}
    assert command["type"] == "setDataMapping"
    assert len(command["requestId"]) == 32
    assert command["mapping"] == expected
    # Nothing waited: the answer arrives later (between cells) and is the live mapping.
    assert viewer.data_mapping == ANNDATA_MAPPING
    viewer._handle_frontend_message(
        _event(viewer._viewer_id, origin="notebook", requestId=command["requestId"], mapping=command["mapping"])
    )
    assert viewer.data_mapping == expected
    assert viewer.data_mapping_origin == "notebook"


def test_a_refused_mapping_arrives_as_a_command_error() -> None:
    viewer = _viewer_with_mapping(ANNDATA_MAPPING)
    errors: list[dict[str, Any]] = []
    viewer.on_command_error(errors.append)
    viewer._handle_frontend_message(
        {"type": "command_error", "command": "setDataMapping", "reason": "obsp['distances'] is not symmetric"}
    )
    assert errors[0]["reason"] == "obsp['distances'] is not symmetric"
    assert viewer.state.command_error["command"] == "setDataMapping"
    assert viewer.data_mapping == ANNDATA_MAPPING
