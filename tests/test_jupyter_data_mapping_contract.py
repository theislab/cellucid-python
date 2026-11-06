"""The data mapping on the Jupyter wire: the ``data_mapping`` event and the
``mapping.set`` action.

The viewer reports its data mapping with a ``data_mapping`` event: once its
dataset opens, and whenever the mapping changes (``origin`` says why;
``notebook`` for this notebook's own ``mapping.set``). A notebook changes it
with the ``mapping.set`` action (``viewer.mapping_set(...)``), which names
only the groups it changes and is answered like every action: with its
outcome, or with the ``command_error`` of its request. A reported mapping is
validated by the same rule sessions use.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import pytest
from _serverless_viewer import serverless_viewer

import cellucid.jupyter as jupyter
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


def _event(**overrides: Any) -> dict[str, Any]:
    event: dict[str, Any] = {
        "type": "data_mapping",
        "origin": "default",
        "mapping": copy.deepcopy(ANNDATA_MAPPING),
    }
    event.update(overrides)
    return event


# ---------------------------------------------------------------------------
# Declarations
# ---------------------------------------------------------------------------


def test_the_wire_declares_exactly_these_field_sets() -> None:
    assert wire._INBOUND_EVENT_FIELDS["data_mapping"] == frozenset({"type", "origin", "mapping"})
    # Changing the mapping is an action like every other: no command of its own.
    assert "setDataMapping" not in wire._COMMAND_TYPES
    assert "mapping.set" in jupyter.action_catalog()


def test_the_capability_document_lists_the_event() -> None:
    document = wire._protocol_capability_document()
    assert "data_mapping" in document["events"]
    assert "action" in document["commands"]


# ---------------------------------------------------------------------------
# Inbound: data_mapping
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("origin", ["default", "best_guess", "session", "user", "notebook"])
def test_a_reported_mapping_names_why_it_is_what_it_is(origin: str) -> None:
    event_type, payload = wire._require_inbound_jupyter_event(_event(origin=origin))
    assert event_type == "data_mapping"
    assert payload == {"origin": origin, "mapping": ANNDATA_MAPPING}


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"origin": "unknown"}, "origin must be one of"),
        ({"origin": None}, "origin must be one of"),
        # Schema 2's answer to setDataMapping named its request; none does now.
        ({"requestId": None}, "unknown requestId"),
        ({"mapping": {"kind": "h5ad"}}, "kind"),
        ({"mapping": {**ANNDATA_MAPPING, "extra": True}}, "unknown extra"),
    ],
)
def test_a_malformed_report_is_refused(overrides: dict[str, Any], match: str) -> None:
    with pytest.raises((TypeError, ValueError), match=match):
        wire._require_inbound_jupyter_event(_event(**overrides))


@pytest.mark.parametrize("case", _invalid_mapping_cases(), ids=lambda case: case["name"])
def test_a_report_with_an_invalid_conformance_mapping_is_refused(case: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        wire._require_inbound_jupyter_event(_event(mapping=case["value"]))


def test_a_report_missing_a_field_is_refused() -> None:
    event = _event()
    del event["origin"]
    with pytest.raises(ValueError, match="missing origin"):
        wire._require_inbound_jupyter_event(event)


# ---------------------------------------------------------------------------
# The viewer: data_mapping, data_mapping_origin, on_data_mapping
# ---------------------------------------------------------------------------


def test_the_viewer_holds_the_latest_reported_mapping() -> None:
    viewer = serverless_viewer()
    assert viewer.data_mapping is None
    assert viewer.data_mapping_origin is None
    hooked: list[dict[str, Any]] = []
    viewer.on_data_mapping(hooked.append)

    viewer._handle_frontend_message(_event(origin="best_guess"))

    assert viewer.data_mapping == ANNDATA_MAPPING
    assert viewer.data_mapping_origin == "best_guess"
    assert hooked == [{"origin": "best_guess", "mapping": ANNDATA_MAPPING}]
    assert viewer.state.data_mapping == hooked[0]

    # The property is a copy: editing it does not change what the viewer holds.
    copied = viewer.data_mapping
    assert copied is not None
    copied["embeddings"]["defaults"]["2d"] = None
    assert viewer.data_mapping == ANNDATA_MAPPING


# ---------------------------------------------------------------------------
# The viewer: mapping_set (the mapping.set action)
# ---------------------------------------------------------------------------


def test_mapping_set_names_only_the_groups_it_changes_and_its_answer_is_the_live_mapping(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = serverless_viewer()
    sent: list[dict[str, Any]] = []
    monkeypatch.setattr(viewer, "send_message", sent.append)
    request = viewer.mapping_set(knn="connectivities", expression="layer:counts")
    assert sent == [{
        "type": "action", "requestId": request.request_id, "action": "mapping.set",
        "params": {"expression": "layer:counts", "knn": "connectivities"}, "turn": None, "turnLabel": None,
    }]
    # Null is a value of its own ("no graph"), never confused with leaving
    # the group out.
    viewer.mapping_set(knn=None)
    assert sent[-1]["params"] == {"knn": None}
    wire._require_frontend_message(sent[-1])

    # The viewer reports the mapping it now shows (origin notebook), then
    # answers the action.
    live = {**ANNDATA_MAPPING, "knn": "connectivities", "expression": "layer:counts"}
    viewer._handle_frontend_message(_event(origin="notebook", mapping=live))
    viewer._handle_frontend_message({
        "type": "action_result", "requestId": request.request_id,
        "step": {"id": "s2", "label": "Changed the data mapping: neighbour graph, expression"},
        "body": bytearray(b'{"groups":2}'),
    })
    assert request.result(timeout=0) == {"groups": 2}
    assert request.step == {"id": "s2", "label": "Changed the data mapping: neighbour graph, expression"}
    assert viewer.data_mapping == live
    assert viewer.data_mapping_origin == "notebook"


def test_a_refused_mapping_fails_its_request_and_reaches_the_hooks() -> None:
    viewer = serverless_viewer()
    errors: list[dict[str, Any]] = []
    viewer.on_command_error(errors.append)
    request = viewer._expect_action("mapping.set")
    viewer._handle_frontend_message({
        "type": "command_error", "command": "action", "requestId": request.request_id,
        "reason": "obsp['distances'] is not symmetric",
    })
    assert errors == [{"command": "action", "requestId": request.request_id, "reason": "obsp['distances'] is not symmetric"}]
    with pytest.raises(jupyter.ActionError, match="not symmetric"):
        request.result(timeout=0)


def test_mapping_set_refuses_a_group_the_mapping_does_not_have(monkeypatch: pytest.MonkeyPatch) -> None:
    viewer = serverless_viewer()
    monkeypatch.setattr(viewer, "send_message", lambda _message: True)
    with pytest.raises(TypeError, match="has no parameter 'geneNames'"):
        viewer.act("mapping.set", geneNames="index")
    with pytest.raises(TypeError, match="unexpected keyword argument 'kind'"):
        viewer.mapping_set(kind="anndata")  # type: ignore[call-arg]
