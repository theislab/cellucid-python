"""
Exact value checks for the Jupyter viewer wire protocol.

Every message that crosses the boundary between Python and the viewer iframe is
checked here before it is trusted: outbound commands are validated and
serialized by :func:`_require_frontend_message`, inbound events are validated
and destructured by :func:`_require_inbound_jupyter_event`, and the diagnostics
path validates the dataset list it probes. The two declarations those checks
read from are also what :func:`_protocol_capability_document` publishes on
``GET /_cellucid/protocol``, so a viewer can discover what this installation
accepts instead of assuming it.

The checks are exact rather than permissive -- ``type(value) is str`` instead of
``isinstance``, exact key sets instead of subsets -- because a viewer running in
a browser is not a trusted caller. This module holds the standard library only,
so the protocol can be reasoned about without loading the viewer machinery.
"""

from __future__ import annotations

import json
import math
from typing import Any
from urllib.parse import urlparse

# Every command Python may send into the viewer iframe, mapped to the complete
# key set that command must carry. This mapping is the single declaration of
# what a notebook command is: :func:`_require_frontend_message` validates
# outbound commands against it, and an inbound ``command_error`` event may only
# name a command that appears here.
_COMMAND_MESSAGE_FIELDS: dict[str, frozenset[str]] = {
    "ping": frozenset({"type", "requestId"}),
    "debug_snapshot": frozenset({"type", "requestId"}),
    "requestSessionBundle": frozenset({"type", "requestId"}),
    "highlight": frozenset({"type", "cells", "color"}),
    "clearHighlights": frozenset({"type"}),
    "setColorBy": frozenset({"type", "field"}),
    "setVisibility": frozenset({"type", "cells", "visible"}),
    "resetCamera": frozenset({"type"}),
    "freeze": frozenset({"type"}),
}

# Every event the viewer may send back, mapped to the complete key set that
# event must carry. This mapping is the single declaration of what an inbound
# event is: :func:`_require_inbound_jupyter_event` validates against it and
# :func:`_protocol_capability_document` publishes its keys, so a viewer can ask
# what this installation accepts instead of guessing from a version number.
_INBOUND_EVENT_FIELDS: dict[str, frozenset[str]] = {
    "selection": frozenset({"type", "viewerId", "cells", "source"}),
    "hover": frozenset({"type", "viewerId", "cell", "position"}),
    "click": frozenset({"type", "viewerId", "cell", "button", "shift", "ctrl"}),
    "ready": frozenset({"type", "viewerId", "n_cells", "dimensions"}),
    "pong": frozenset({"type", "viewerId", "requestId", "t"}),
    "debug_snapshot": frozenset(
        {
            "type",
            "viewerId",
            "requestId",
            "ts",
            "locationHref",
            "origin",
            "serverUrl",
            "connected",
            "parentOrigin",
            "userAgent",
        }
    ),
    "session_bundle": frozenset(
        {
            "type",
            "viewerId",
            "requestId",
            "status",
            "bytes",
            "path",
        }
    ),
    "command_error": frozenset({"type", "viewerId", "command", "reason"}),
}

# The optional pair a prepared export adds to its ``/_cellucid/datasets`` entry
# when it ships a published default session. ``cellucid.server._state`` writes
# it, ``cellucid.server._artifacts`` reads it back to decide which sidecars the
# dataset may serve, and the viewer's local-demo source requires exactly this
# manifest name -- so an entry carrying the pair is correct, not malformed, and
# the diagnostic validator has to accept it or it loses the identity probes.
#
# The names are restated here rather than imported because this module holds
# the standard library only; ``tests/test_diagnostic_dataset_list_contract.py``
# asserts the two declarations agree, which is what was missing when the
# published-state pair was added and the validator was not.
_DATASET_ENTRY_FIELDS = frozenset({"id", "path", "name"})
_PUBLISHED_STATE_FIELDS = frozenset({"state_manifest", "state_sha256"})
_PUBLISHED_STATE_MANIFEST_NAME = "state-snapshots.json"
_SHA256_HEX_DIGITS = frozenset("0123456789abcdef")

# Longest ``command_error`` reason the viewer may report. The viewer truncates
# its own text well below this, so the bound refuses an unbounded string
# without turning a change to the viewer's presentation limit into a protocol
# break.
_COMMAND_ERROR_REASON_LIMIT = 500


def _protocol_capability_document() -> dict[str, list[str]]:
    """Return what this installation accepts on the viewer wire.

    This is the body of ``GET /_cellucid/protocol``, and it exists because the
    two halves of the protocol ship in different releases: the Python package
    loads the viewer from ``https://www.cellucid.com`` without pinning a build,
    so a viewer newer than the notebook's ``cellucid`` is the normal state the
    day a build deploys. A viewer that emits an event this installation does not
    declare does not degrade to silence -- the validator raises, the event
    endpoint logs a traceback into the kernel and answers ``500``, and the
    viewer throws. So the viewer has to ask first, and no field already in
    flight can carry the answer: the config, the health record and every command
    acknowledgement are validated against exact key sets on the viewer side, so
    adding a field to any of them breaks builds that are already deployed. A
    route an older build never requests is the only channel that stays quiet.

    Both keys are lists that grow, never a version to compare against. Growth
    happens inside the arrays so this document's own key set can stay fixed, and
    a viewer therefore tests for the capability it needs rather than inferring
    it from a release number.

    ``events``
        Every inbound event type :func:`_require_inbound_jupyter_event` accepts.
        A viewer sends only the ones named here.
    ``commands``
        Every command this installation may send, which is also exactly the set
        of names a ``command_error`` may report as having failed.
    """
    return {
        "events": sorted(_INBOUND_EVENT_FIELDS),
        "commands": sorted(_COMMAND_MESSAGE_FIELDS),
    }


def _require_client_server_url(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("client_server_url must be a string")
    if not value or value != value.strip() or any(character.isspace() for character in value):
        raise ValueError("client_server_url must be a non-empty URL without surrounding whitespace")
    if value.endswith("/"):
        raise ValueError("client_server_url must not end with '/'")
    parsed = urlparse(value)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        raise ValueError("client_server_url must be an absolute HTTP or HTTPS URL")
    try:
        if parsed.hostname is None:
            raise ValueError("client_server_url must contain a hostname")
        _ = parsed.port
    except ValueError as error:
        raise ValueError("client_server_url contains an invalid host or port") from error
    if parsed.username is not None or parsed.password is not None:
        raise ValueError("client_server_url must not contain credentials")
    if parsed.query or parsed.fragment:
        raise ValueError("client_server_url must not contain a query or fragment")
    return value


def _require_exact_message_text(
    value: object,
    *,
    label: str,
    allow_whitespace: bool = False,
) -> str:
    if (
        type(value) is not str
        or not value
        or value != value.strip()
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
        or (not allow_whitespace and any(character.isspace() for character in value))
    ):
        raise TypeError(f"{label} must be exact non-empty text")
    return value


def _require_exact_json_value(
    value: object,
    *,
    label: str,
    ancestors: set[int] | None = None,
) -> None:
    if value is None or type(value) in {str, bool, int}:
        return
    if type(value) is float:
        if not math.isfinite(value):
            raise ValueError(f"{label} numbers must be finite JSON numbers")
        return
    if ancestors is None:
        ancestors = set()
    identity = id(value)
    if identity in ancestors:
        raise ValueError(f"{label} must not contain JSON cycles")
    ancestors.add(identity)
    try:
        if type(value) is list:
            for index, entry in enumerate(value):
                _require_exact_json_value(
                    entry,
                    label=f"{label}[{index}]",
                    ancestors=ancestors,
                )
            return
        if type(value) is dict:
            for key, entry in value.items():
                if type(key) is not str:
                    raise TypeError(f"{label} object keys must be native strings")
                _require_exact_json_value(
                    entry,
                    label=f"{label}.{key}",
                    ancestors=ancestors,
                )
            return
        raise TypeError(f"{label} must contain only exact JSON values")
    finally:
        ancestors.remove(identity)


def _require_debug_dataset_list(payload: object) -> list[dict[str, str]]:
    """Validate the exact dataset-list response used by diagnostics."""
    if type(payload) is not dict or set(payload) != {"datasets"}:
        raise TypeError("Dataset-list response must contain exactly 'datasets'")
    raw_datasets = payload["datasets"]
    if type(raw_datasets) is not list:
        raise TypeError("Dataset-list response 'datasets' must be a native list")

    datasets: list[dict[str, str]] = []
    dataset_ids: set[str] = set()
    dataset_paths: set[str] = set()
    for index, entry in enumerate(raw_datasets):
        if type(entry) is not dict or set(entry) not in (
            _DATASET_ENTRY_FIELDS,
            _DATASET_ENTRY_FIELDS | _PUBLISHED_STATE_FIELDS,
        ):
            raise TypeError(
                f"Dataset-list entry {index} must contain exactly id, name, and path, "
                "and may add state_manifest and state_sha256 only as a complete pair"
            )
        dataset_id = _require_exact_message_text(
            entry["id"],
            label=f"Dataset-list entry {index} id",
            allow_whitespace=True,
        )
        dataset_name = _require_exact_message_text(
            entry["name"],
            label=f"Dataset-list entry {index} name",
            allow_whitespace=True,
        )
        path = entry["path"]
        if (
            type(path) is not str
            or not path.startswith("/")
            or not path.endswith("/")
            or any(
                not (character.isascii() and (character.isalnum() or character in "/._-"))
                for character in path
            )
            or (
                path != "/" and any(segment in {"", ".", ".."} for segment in path[1:-1].split("/"))
            )
        ):
            raise TypeError(
                f"Dataset-list entry {index} path must be one exact absolute "
                "portable directory route"
            )
        if path != "/" and path.count("/") != 2:
            raise TypeError(
                f"Dataset-list entry {index} path must identify one served dataset directory"
            )
        if dataset_id in dataset_ids:
            raise ValueError(f"Dataset-list response duplicates id {dataset_id!r}")
        if path in dataset_paths:
            raise ValueError(f"Dataset-list response duplicates path {path!r}")
        dataset_ids.add(dataset_id)
        dataset_paths.add(path)
        validated = {
            "id": dataset_id,
            "path": path,
            "name": dataset_name,
        }
        if set(entry) >= _PUBLISHED_STATE_FIELDS:
            manifest = entry["state_manifest"]
            if manifest != _PUBLISHED_STATE_MANIFEST_NAME:
                raise ValueError(
                    f"Dataset-list entry {index} state_manifest must be exactly "
                    f"{_PUBLISHED_STATE_MANIFEST_NAME}"
                )
            digest = entry["state_sha256"]
            if (
                type(digest) is not str
                or len(digest) != 64
                or any(character not in _SHA256_HEX_DIGITS for character in digest)
            ):
                raise ValueError(
                    f"Dataset-list entry {index} state_sha256 must be exactly 64 "
                    "lowercase hexadecimal characters"
                )
            validated["state_manifest"] = manifest
            validated["state_sha256"] = digest
        datasets.append(validated)
    return datasets


def _require_cell_indices(value: object, *, label: str = "cell_indices") -> list[int]:
    if type(value) is not list:
        raise TypeError(f"{label} must be a native list")
    seen: set[int] = set()
    for index, cell_index in enumerate(value):
        if type(cell_index) is not int or cell_index < 0:
            raise TypeError(f"{label}[{index}] must be a non-negative native integer")
        if cell_index in seen:
            raise ValueError(f"{label} must not contain duplicate cell indices")
        seen.add(cell_index)
    return value


def _require_frontend_message(message: object) -> str:
    if type(message) is not dict:
        raise TypeError("message must be a native dictionary")
    _require_exact_json_value(message, label="message")
    message_type = _require_exact_message_text(
        message.get("type"),
        label="message type",
    )
    if "viewerId" in message or "viewerToken" in message:
        raise TypeError("message must not supply viewer routing credentials")

    expected = _COMMAND_MESSAGE_FIELDS.get(message_type)
    if expected is None:
        raise ValueError(f"Unknown current Jupyter command type {message_type!r}")
    if set(message) != expected:
        fields = ", ".join(sorted(expected))
        raise TypeError(f"{message_type} message must contain exactly {fields}")

    if message_type in {"ping", "debug_snapshot", "requestSessionBundle"}:
        _require_exact_message_text(
            message["requestId"],
            label=f"{message_type} requestId",
        )
    elif message_type == "highlight":
        highlight_cells = _require_cell_indices(
            message["cells"],
            label="highlight cells",
        )
        if not highlight_cells:
            raise ValueError("highlight cells must be a non-empty list")
        color = message["color"]
        if (
            type(color) is not str
            or len(color) != 7
            or color[0] != "#"
            or any(character not in "0123456789abcdefABCDEF" for character in color[1:])
        ):
            raise ValueError("highlight color must be an exact six-digit hex color")
    elif message_type == "setColorBy":
        _require_exact_message_text(
            message["field"],
            label="setColorBy field",
            allow_whitespace=True,
        )
    elif message_type == "setVisibility":
        cells = message["cells"]
        if cells is not None:
            _require_cell_indices(cells, label="setVisibility cells")
        if type(message["visible"]) is not bool:
            raise TypeError("setVisibility visible must be exactly True or False")

    return json.dumps(
        message,
        ensure_ascii=True,
        allow_nan=False,
        separators=(",", ":"),
    )


def _require_nonnegative_event_integer(value: object, *, label: str) -> int:
    if type(value) is not int:
        raise TypeError(f"{label} must be a native integer")
    if value < 0:
        raise ValueError(f"{label} must be non-negative")
    return value


def _require_exact_event_fields(
    message: dict[str, Any],
    expected: frozenset[str],
    *,
    event_type: str,
) -> None:
    missing = sorted(expected - set(message))
    unknown = sorted(set(message) - expected)
    if missing or unknown:
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if unknown:
            details.append("unknown " + ", ".join(unknown))
        raise ValueError(
            f"Inbound Jupyter {event_type} event has invalid fields ({'; '.join(details)})"
        )


def _require_inbound_jupyter_event(
    message: object,
    *,
    expected_viewer_id: str,
) -> tuple[str, dict[str, Any]]:
    if type(message) is not dict:
        raise TypeError("Inbound Jupyter event must be a native dictionary")
    _require_exact_json_value(message, label="Inbound Jupyter event")
    event_type = _require_exact_message_text(
        message.get("type"),
        label="Inbound Jupyter event type",
    )
    viewer_id = _require_exact_message_text(
        message.get("viewerId"),
        label="Inbound Jupyter event viewerId",
    )
    if viewer_id != expected_viewer_id:
        raise ValueError("Inbound Jupyter event viewerId does not match this viewer")

    expected_fields = _INBOUND_EVENT_FIELDS.get(event_type)
    if expected_fields is None:
        raise ValueError(f"Unknown inbound Jupyter event type {event_type!r}")
    _require_exact_event_fields(
        message,
        expected_fields,
        event_type=event_type,
    )

    if event_type == "selection":
        _require_cell_indices(message["cells"], label="selection cells")
        _require_exact_message_text(
            message["source"],
            label="selection source",
        )
    elif event_type == "hover":
        cell = message["cell"]
        if cell is not None:
            _require_nonnegative_event_integer(cell, label="hover cell")
        position = message["position"]
        if position is not None:
            if type(position) is not dict:
                raise TypeError("hover position must be a native dictionary or None")
            _require_exact_event_fields(
                position,
                frozenset({"x", "y", "z"}),
                event_type="hover position",
            )
            for axis in ("x", "y", "z"):
                coordinate = position[axis]
                if type(coordinate) not in {int, float} or not math.isfinite(coordinate):
                    raise TypeError(f"hover position {axis} must be a finite number")
    elif event_type == "click":
        _require_nonnegative_event_integer(message["cell"], label="click cell")
        button = message["button"]
        if type(button) is not int or button not in {0, 1, 2}:
            raise TypeError("click button must be the native integer 0, 1, or 2")
        if type(message["shift"]) is not bool or type(message["ctrl"]) is not bool:
            raise TypeError("click shift and ctrl must be native booleans")
    elif event_type == "ready":
        _require_nonnegative_event_integer(message["n_cells"], label="ready n_cells")
        if type(message["dimensions"]) is not int or message["dimensions"] not in {
            1,
            2,
            3,
        }:
            raise TypeError("ready dimensions must be the native integer 1, 2, or 3")
    elif event_type == "pong":
        _require_exact_message_text(
            message["requestId"],
            label="pong requestId",
        )
        _require_nonnegative_event_integer(message["t"], label="pong t")
    elif event_type == "debug_snapshot":
        for key in (
            "requestId",
            "ts",
            "locationHref",
            "origin",
            "serverUrl",
            "parentOrigin",
        ):
            _require_exact_message_text(
                message[key],
                label=f"debug_snapshot {key}",
                allow_whitespace=True,
            )
        if type(message["connected"]) is not bool:
            raise TypeError("debug_snapshot connected must be a native boolean")
        user_agent = message["userAgent"]
        if user_agent is not None:
            _require_exact_message_text(
                user_agent,
                label="debug_snapshot userAgent",
                allow_whitespace=True,
            )
    elif event_type == "command_error":
        # The viewer names the failed command from its own frozen list of the
        # command types this protocol defines. Accepting only that list is what
        # keeps browser-supplied text out of the notebook's failure notice.
        command = _require_exact_message_text(
            message["command"],
            label="command_error command",
        )
        if command not in _COMMAND_MESSAGE_FIELDS:
            raise ValueError(
                "command_error command must name one current Jupyter command, "
                f"not {command!r}"
            )
        reason = _require_exact_message_text(
            message["reason"],
            label="command_error reason",
            allow_whitespace=True,
        )
        if len(reason) > _COMMAND_ERROR_REASON_LIMIT:
            raise ValueError(
                "command_error reason must be at most "
                f"{_COMMAND_ERROR_REASON_LIMIT} characters"
            )
    else:
        if message["status"] != "ok":
            raise ValueError("session_bundle status must equal 'ok'")
        _require_exact_message_text(
            message["requestId"],
            label="session_bundle requestId",
        )
        _require_nonnegative_event_integer(
            message["bytes"],
            label="session_bundle bytes",
        )
        _require_exact_message_text(
            message["path"],
            label="session_bundle path",
            allow_whitespace=True,
        )

    return (
        event_type,
        {key: value for key, value in message.items() if key not in {"type", "viewerId"}},
    )
