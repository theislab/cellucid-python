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


def test_set_data_mapping_needs_a_reported_mapping() -> None:
    viewer = serverless_viewer()
    with pytest.raises(RuntimeError, match="has not reported its data mapping yet"):
        viewer.set_data_mapping(knn="connectivities")


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


def test_set_data_mapping_sends_the_merged_mapping_and_returns_the_answer(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = _viewer_with_mapping(ANNDATA_MAPPING)
    sent: list[str] = []

    def answer(message: dict[str, Any]) -> None:
        sent.append(wire._require_frontend_message(message))
        viewer._handle_frontend_message(
            _event(
                viewer._viewer_id,
                origin="notebook",
                requestId=message["requestId"],
                mapping=message["mapping"],
            )
        )

    monkeypatch.setattr(viewer, "send_message", answer)
    applied = viewer.set_data_mapping(
        embeddings={"enabled": [{"key": "X_pca", "columns": [3, 1]}], "defaults": {"1d": None, "2d": {"key": "X_pca", "columns": [3, 1]}, "3d": None}},
        knn="connectivities",
    )

    expected = {
        **ANNDATA_MAPPING,
        "embeddings": {"enabled": [{"key": "X_pca", "columns": [3, 1]}], "defaults": {"1d": None, "2d": {"key": "X_pca", "columns": [3, 1]}, "3d": None}},
        "knn": "connectivities",
    }
    assert applied == expected
    assert viewer.data_mapping == expected
    assert viewer.data_mapping_origin == "notebook"
    (serialized,) = sent
    command = json.loads(serialized)
    assert set(command) == {"type", "requestId", "mapping"}
    assert command["type"] == "setDataMapping"
    assert len(command["requestId"]) == 32
    assert command["mapping"] == expected


def test_set_data_mapping_ignores_answers_to_other_requests(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = _viewer_with_mapping(ANNDATA_MAPPING)

    def answer_someone_else(message: dict[str, Any]) -> None:
        viewer._handle_frontend_message(
            _event(viewer._viewer_id, origin="notebook", requestId="another-request")
        )
        viewer._handle_frontend_message(_event(viewer._viewer_id, origin="user"))

    monkeypatch.setattr(viewer, "send_message", answer_someone_else)
    with pytest.raises(TimeoutError, match="apply the data mapping"):
        viewer.set_data_mapping(knn="connectivities", timeout=0.05)


def test_set_data_mapping_raises_the_viewer_refusal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = _viewer_with_mapping(ANNDATA_MAPPING)
    viewer.on_command_error(lambda _event: None)

    def refuse(message: dict[str, Any]) -> None:
        viewer._handle_frontend_message(
            {
                "type": "command_error",
                "command": "setDataMapping",
                "reason": "obsp['distances'] is not symmetric",
            }
        )

    monkeypatch.setattr(viewer, "send_message", refuse)
    with pytest.raises(RuntimeError, match="did not apply the data mapping: obsp"):
        viewer.set_data_mapping(knn="distances")
    assert viewer.data_mapping == ANNDATA_MAPPING


def test_an_earlier_refusal_does_not_fail_a_later_request(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    viewer = _viewer_with_mapping(ANNDATA_MAPPING)
    viewer.on_command_error(lambda _event: None)
    viewer._handle_frontend_message(
        {
            "type": "command_error",
            "command": "setDataMapping",
            "reason": "an older refusal",
        }
    )

    def answer(message: dict[str, Any]) -> None:
        viewer._handle_frontend_message(
            _event(
                viewer._viewer_id,
                origin="notebook",
                requestId=message["requestId"],
                mapping=message["mapping"],
            )
        )

    monkeypatch.setattr(viewer, "send_message", answer)
    assert viewer.set_data_mapping(latent="X_pca")["latent"] == "X_pca"
