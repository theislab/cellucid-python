"""
Exact value checks for the viewer channel.

Every message that crosses the viewer channel (:mod:`cellucid._viewer_channel`,
one WebSocket per viewer) is checked here before it is trusted: outbound
commands are validated and serialized by :func:`_require_frontend_message` (a
command that names cells as one binary frame, every other one as JSON text),
inbound events are validated and destructured by
:func:`_require_inbound_jupyter_event`, and the diagnostics path validates the
dataset list it probes. The declarations those checks read from are also what
:func:`_protocol_capability_document` publishes on ``GET /_cellucid/protocol``,
so a viewer can discover what this installation accepts instead of assuming
it. They are the cross-language rules of ``viewer-protocol.json`` (the web
app's ``data/viewer-protocol.js``), and a test holds them equal.

The checks are exact rather than permissive -- ``type(value) is str`` instead of
``isinstance``, exact key sets instead of subsets -- because a viewer running in
a browser is not a trusted caller. This module imports only the standard
library, NumPy (to read a binary ``selection`` body and write a command's
cells in one vectorized pass),
and :mod:`cellucid.data_mapping`, which is standard-library only itself, so the
protocol can be reasoned about without loading the viewer machinery.
"""

from __future__ import annotations

import json
import math
from typing import Any
from urllib.parse import urlparse

import numpy as np

from .._plural import plural
from ..data_mapping import DataMappingError, validate_data_mapping

# Every command Python sends to the viewer as a JSON text frame, mapped to the
# complete key set that command must carry (the shared fixture's ``commands``).
_COMMAND_MESSAGE_FIELDS: dict[str, frozenset[str]] = {
    "ping": frozenset({"type", "requestId"}),
    "debug_snapshot": frozenset({"type", "requestId"}),
    "requestSessionBundle": frozenset({"type", "requestId"}),
    "clearHighlights": frozenset({"type"}),
    "setColorBy": frozenset({"type", "field"}),
    "resetCamera": frozenset({"type"}),
    "freeze": frozenset({"type"}),
    "setDataMapping": frozenset({"type", "requestId", "mapping"}),
}

# Every command that names cells, which Python sends as one binary frame: a
# ``uint32`` little-endian header length, the UTF-8 JSON header (exactly these
# keys, the shared fixture's ``binaryCommands``), padded with spaces so the
# body starts on a 4-byte boundary, then the cells as little-endian ``uint32``,
# strictly ascending. ``setVisibility`` of every cell is ``all: true`` and no
# cells. Two million cells are 8 MB written in one pass, where the same cells
# as JSON text were 14 MB to build and parse.
_BINARY_COMMAND_HEADER_FIELDS: dict[str, frozenset[str]] = {
    "highlight": frozenset({"type", "color"}),
    "setVisibility": frozenset({"type", "visible", "all"}),
}

# What a caller of ``send_message`` gives for a binary command: its cells under
# ``cells`` (None for every cell of a ``setVisibility``), beside the header's
# fields.
_BINARY_COMMAND_MESSAGE_FIELDS: dict[str, frozenset[str]] = {
    "highlight": frozenset({"type", "cells", "color"}),
    "setVisibility": frozenset({"type", "cells", "visible"}),
}

# Every command Python may send. This is the single declaration of what a
# notebook command is: :func:`_require_frontend_message` validates outbound
# commands against it, and an inbound ``command_error`` event may only name a
# command that appears here.
_COMMAND_TYPES: frozenset[str] = frozenset(_COMMAND_MESSAGE_FIELDS) | frozenset(
    _BINARY_COMMAND_MESSAGE_FIELDS
)

# Every event the viewer may send back, mapped to the complete key set that
# event must carry. This mapping is the single declaration of what an inbound
# event is: :func:`_require_inbound_jupyter_event` validates against it and
# :func:`_protocol_capability_document` publishes its keys, so a viewer can ask
# what this installation accepts instead of guessing from a version number.
# Events carry no viewer id: the channel's socket is one viewer's.
#
# Two events do not arrive as JSON text frames: the channel builds
# ``selection`` from a binary frame (its ``cells`` field is the frame's body,
# little-endian uint32, until :func:`_decode_selection_cells` reads it) and
# ``session_bundle`` from a binary frame it streams to a temporary file.
_INBOUND_EVENT_FIELDS: dict[str, frozenset[str]] = {
    "selection": frozenset({"type", "cells", "source"}),
    "hover": frozenset({"type", "cell", "position"}),
    "click": frozenset({"type", "cell", "button", "shift", "ctrl"}),
    "ready": frozenset({"type", "n_cells", "dimensions"}),
    "pong": frozenset({"type", "requestId", "t"}),
    "debug_snapshot": frozenset(
        {
            "type",
            "requestId",
            "ts",
            "locationHref",
            "origin",
            "serverUrl",
            "connected",
            "userAgent",
        }
    ),
    "session_bundle": frozenset({"type", "requestId", "status", "bytes", "path"}),
    "command_error": frozenset({"type", "command", "reason"}),
    "data_mapping": frozenset({"type", "requestId", "origin", "mapping"}),
}

# Why the viewer's data mapping is what it is, as a ``data_mapping`` event
# reports it. ``notebook`` is the one origin that answers a ``setDataMapping``
# command, so it is the one origin that carries a ``requestId``.
_DATA_MAPPING_ORIGINS = frozenset(
    {"default", "best_guess", "session", "user", "notebook"}
)

#: The one entry an AnnData server lists: the store the viewer reads in place.
_ANNDATA_ENTRY_FIELDS = frozenset({"id", "name", "anndata"})

# Longest ``command_error`` reason the viewer may report. The viewer truncates
# its own text well below this, so the bound refuses an unbounded string
# without turning a change to the viewer's presentation limit into a protocol
# break.
_COMMAND_ERROR_REASON_LIMIT = 500


#: Where the viewer channel is, relative to the server's base address.
_CHANNEL_ROUTE = "_cellucid/ws"


def _protocol_capability_document() -> dict[str, Any]:
    """Return the viewer channel and what this installation accepts on it.

    This is the body of ``GET /_cellucid/protocol``, and it exists because the
    two halves of the protocol ship in different releases: the Python package
    loads the viewer from ``https://www.cellucid.com`` without pinning a build,
    so a viewer newer than the notebook's ``cellucid`` is the normal state the
    day a build deploys. A viewer that sends an event this installation does
    not declare would have its channel closed as a protocol violation, so the
    viewer asks first.

    The lists grow, never a version to compare against, so a viewer tests for
    the capability it needs rather than inferring it from a release number.

    ``channel``
        Where the viewer opens its WebSocket, relative to the server's base
        address.
    ``events``
        Every inbound event type :func:`_require_inbound_jupyter_event` accepts.
        A viewer sends only the ones named here.
    ``commands``
        Every command this installation may send, which is also exactly the set
        of names a ``command_error`` may report as having failed.
    ``binaryCommands``
        The commands it sends as binary frames: those that name cells.
    """
    return {
        "channel": _CHANNEL_ROUTE,
        "commands": sorted(_COMMAND_TYPES),
        "binaryCommands": sorted(_BINARY_COMMAND_HEADER_FIELDS),
        "events": sorted(_INBOUND_EVENT_FIELDS),
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


_EXACT_JSON_SCALAR_TYPES = frozenset({str, bool, int})


def _require_exact_json_value(
    value: object,
    *,
    label: str,
    ancestors: set[int] | None = None,
) -> None:
    if value is None or type(value) in _EXACT_JSON_SCALAR_TYPES:
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
                # Scalars are checked inline: a label is built only for an
                # entry that fails or that nests further, never per element.
                if entry is None or type(entry) in _EXACT_JSON_SCALAR_TYPES:
                    continue
                if type(entry) is float and math.isfinite(entry):
                    continue
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


def _require_debug_dataset_list(payload: object) -> list[dict[str, Any]]:
    """Validate the exact dataset-list response used by diagnostics.

    An AnnData server lists one dataset: its id, its name, and the store the
    viewer reads (``anndata``: the store's ``file`` under
    ``_cellucid/anndata/`` and the launch ``defaults``).
    """
    if type(payload) is not dict or set(payload) != {"datasets"}:
        raise TypeError("Dataset-list response must contain exactly 'datasets'")
    raw_datasets = payload["datasets"]
    if type(raw_datasets) is not list:
        raise TypeError("Dataset-list response 'datasets' must be a native list")

    datasets: list[dict[str, Any]] = []
    dataset_ids: set[str] = set()
    for index, entry in enumerate(raw_datasets):
        if type(entry) is not dict or set(entry) != _ANNDATA_ENTRY_FIELDS:
            raise TypeError(f"Dataset-list entry {index} must contain exactly id, name, and anndata")
        dataset_id = _require_exact_message_text(
            entry["id"], label=f"Dataset-list entry {index} id", allow_whitespace=True
        )
        if dataset_id in dataset_ids:
            raise ValueError(f"Dataset-list response duplicates id {dataset_id!r}")
        dataset_ids.add(dataset_id)
        anndata = entry["anndata"]
        file = anndata.get("file") if type(anndata) is dict else None
        if (
            type(anndata) is not dict
            or set(anndata) != {"file", "defaults"}
            or type(anndata["defaults"]) is not dict
            or type(file) is not dict
            or type(file.get("url")) is not str
            or not file["url"].startswith("_cellucid/anndata/")
        ):
            raise TypeError(
                f"Dataset-list entry {index} anndata must name one served store "
                "(file with a _cellucid/anndata/ url, and defaults)"
            )
        _require_exact_json_value(anndata, label=f"Dataset-list entry {index} anndata")
        datasets.append(
            {
                "id": dataset_id,
                "name": _require_exact_message_text(
                    entry["name"],
                    label=f"Dataset-list entry {index} name",
                    allow_whitespace=True,
                ),
                "anndata": anndata,
            }
        )
    return datasets


def _require_cell_indices(value: object, *, label: str = "cell_indices") -> np.ndarray:
    """Validate cells a command names; return them as ascending ``uint32``.

    ``value`` is a list of native ints or a one-dimensional NumPy integer
    array, holding at least one cell, each a non-negative index below 2**32,
    none twice. Cells given in any order are sorted once (the wire carries
    them strictly ascending); an array already ascending is not copied when it
    is ``uint32`` already. A list is checked element by element (a bool or a
    float is refused); an array by its dtype, in vectorized passes.
    """
    if type(value) is list:
        for index, cell_index in enumerate(value):
            if type(cell_index) is not int or cell_index < 0:
                raise TypeError(f"{label}[{index}] must be a non-negative native integer")
        try:
            cells = np.array(value, dtype=np.int64) if value else np.empty(0, dtype=np.int64)
        except OverflowError as error:
            raise ValueError(f"{label} must be cell indices below 2**32") from error
    elif isinstance(value, np.ndarray):
        if value.ndim != 1 or value.dtype.kind not in "iu":
            raise TypeError(f"{label} must be a one-dimensional NumPy integer array, not {value.dtype} of {value.ndim} {plural(value.ndim, 'dimension')}")
        cells = value
    else:
        raise TypeError(f"{label} must be a native list of ints or a one-dimensional NumPy integer array")
    if cells.size == 0:
        raise ValueError(f"{label} must be non-empty: name at least one cell")
    if cells.dtype.kind == "i" and int(cells.min()) < 0:
        raise TypeError(f"{label} must be non-negative cell indices")
    if int(cells.max()) > _MAX_CELL_INDEX:
        raise ValueError(f"{label} must be cell indices below 2**32")
    if cells.size > 1 and not bool(np.all(cells[1:] > cells[:-1])):
        cells = np.sort(cells)
        if bool(np.any(cells[1:] == cells[:-1])):
            raise ValueError(f"{label} must not contain duplicate cell indices")
    return cells.astype(_SELECTION_CELL_DTYPE, copy=False)


#: The largest cell index the wire carries: a ``uint32``.
_MAX_CELL_INDEX = 2**32 - 1


# One selected cell in a ``selection`` frame's body: a little-endian unsigned
# 32-bit cell index. The body is the whole selection, strictly ascending, so
# it is unique by construction and checked in one pass.
_SELECTION_CELL_DTYPE = np.dtype("<u4")


def _require_selection_cell_array(body: object, *, label: str = "selection cells") -> np.ndarray:
    """Validate one binary selection body and return a view of its cells.

    The body is a ``selection`` frame's: one or more little-endian ``uint32``
    cell indices, strictly ascending. Strict order is what makes the indices
    unique, so the check is one vectorized comparison rather than a set of
    every index. The returned array is a view of ``body`` (no copy).
    """
    if not isinstance(body, bytes | bytearray):
        raise TypeError(f"{label} must be the exact bytes of one selection body")
    itemsize = _SELECTION_CELL_DTYPE.itemsize
    if not body or len(body) % itemsize:
        raise ValueError(
            f"{label} must be one or more little-endian uint32 cell indices "
            f"({itemsize} {plural(itemsize, 'byte')} each), not {len(body)} {plural(len(body), 'byte')}"
        )
    cells: np.ndarray = np.frombuffer(body, dtype=_SELECTION_CELL_DTYPE)
    if cells.size > 1 and not bool(np.all(cells[1:] > cells[:-1])):
        raise ValueError(f"{label} must be strictly ascending, so each cell appears once")
    return cells


def _decode_selection_cells(body: object, *, label: str = "selection cells") -> list[int]:
    """Return the cell indices of one binary selection body as native ints."""
    cells: list[int] = _require_selection_cell_array(body, label=label).tolist()
    return cells


def _require_data_mapping(value: object, *, label: str) -> dict[str, Any]:
    """Validate one data mapping and return its canonical copy.

    The mapping is checked with :func:`cellucid.data_mapping.validate_data_mapping`,
    the same rule sessions and the web app apply, so a mapping the viewer would
    refuse is refused here before it crosses the wire. A mapping that passes is
    replaced by its canonical copy -- the one that is serialized outbound and
    stored inbound -- so the order of keys is never a property of the caller.
    """
    try:
        return validate_data_mapping(value, label=label)
    except DataMappingError as error:
        raise ValueError(str(error)) from error


def _binary_command_frame(header: dict[str, Any], cells: np.ndarray | None) -> bytearray:
    """One binary command frame: header length, padded JSON header, cells."""
    text = json.dumps(header, ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    header_bytes = (text + " " * (-len(text) % 4)).encode("ascii")
    body_offset = 4 + len(header_bytes)
    count = 0 if cells is None else cells.size
    frame = bytearray(body_offset + 4 * count)
    frame[:4] = len(header_bytes).to_bytes(4, "little")
    frame[4:body_offset] = header_bytes
    if cells is not None:
        # The cells are written straight into the frame: one pass, no copy
        # of them beside it.
        np.frombuffer(frame, dtype=_SELECTION_CELL_DTYPE, count=count, offset=body_offset)[:] = cells
    return frame


def _require_frontend_message(message: object) -> str | bytearray:
    """Validate one command and serialize it as the frame the viewer reads.

    A command that names cells (``highlight``, ``setVisibility``) becomes one
    binary frame (bytes); every other one a JSON text frame (str).
    """
    if type(message) is not dict:
        raise TypeError("message must be a native dictionary")
    message_type = _require_exact_message_text(
        message.get("type"),
        label="message type",
    )
    if "viewerId" in message or "viewerToken" in message:
        raise TypeError("message must not supply viewer routing credentials")

    binary = _BINARY_COMMAND_MESSAGE_FIELDS.get(message_type)
    expected = binary if binary is not None else _COMMAND_MESSAGE_FIELDS.get(message_type)
    if expected is None:
        raise ValueError(f"Unknown viewer command type {message_type!r}")
    if set(message) != expected:
        fields = ", ".join(sorted(expected))
        raise TypeError(f"{message_type} message must contain exactly {fields}")

    if message_type == "highlight":
        cells = _require_cell_indices(message["cells"], label="highlight cells")
        color = message["color"]
        if (
            type(color) is not str
            or len(color) != 7
            or color[0] != "#"
            or any(character not in "0123456789abcdefABCDEF" for character in color[1:])
        ):
            raise ValueError("highlight color must be an exact six-digit hex color")
        return _binary_command_frame({"type": "highlight", "color": color}, cells)
    if message_type == "setVisibility":
        visible = message["visible"]
        if type(visible) is not bool:
            raise TypeError("setVisibility visible must be exactly True or False")
        if message["cells"] is None:
            return _binary_command_frame(
                {"type": "setVisibility", "visible": visible, "all": True}, None
            )
        cells = _require_cell_indices(message["cells"], label="setVisibility cells")
        return _binary_command_frame(
            {"type": "setVisibility", "visible": visible, "all": False}, cells
        )

    _require_exact_json_value(message, label="message")
    if message_type in {"ping", "debug_snapshot", "requestSessionBundle"}:
        _require_exact_message_text(
            message["requestId"],
            label=f"{message_type} requestId",
        )
    elif message_type == "setColorBy":
        _require_exact_message_text(
            message["field"],
            label="setColorBy field",
            allow_whitespace=True,
        )
    elif message_type == "setDataMapping":
        _require_exact_message_text(
            message["requestId"],
            label="setDataMapping requestId",
        )
        message = {
            **message,
            "mapping": _require_data_mapping(
                message["mapping"],
                label="setDataMapping mapping",
            ),
        }

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
            f"Inbound viewer {event_type} event has invalid fields ({'; '.join(details)})"
        )


def _require_inbound_jupyter_event(message: object) -> tuple[str, dict[str, Any]]:
    """Validate one event from the viewer channel; return its type and payload.

    The payload is the event without ``type``: what the notebook's hooks
    receive.
    """
    if type(message) is not dict:
        raise TypeError("Inbound viewer event must be a native dictionary")
    # A selection's cells are the binary body the channel received, not JSON;
    # they are read by `_decode_selection_cells` below.
    _require_exact_json_value(
        (
            {key: value for key, value in message.items() if key != "cells"}
            if message.get("type") == "selection"
            else message
        ),
        label="Inbound viewer event",
    )
    event_type = _require_exact_message_text(
        message.get("type"),
        label="Inbound viewer event type",
    )

    expected_fields = _INBOUND_EVENT_FIELDS.get(event_type)
    if expected_fields is None:
        raise ValueError(f"Unknown inbound viewer event type {event_type!r}")
    _require_exact_event_fields(
        message,
        expected_fields,
        event_type=event_type,
    )

    if event_type == "selection":
        _require_exact_message_text(
            message["source"],
            label="selection source",
        )
        message = {
            **message,
            "cells": _decode_selection_cells(message["cells"]),
        }
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
        if command not in _COMMAND_TYPES:
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
    elif event_type == "data_mapping":
        origin = message["origin"]
        if type(origin) is not str or origin not in _DATA_MAPPING_ORIGINS:
            raise ValueError(
                "data_mapping origin must be one of "
                f"{sorted(_DATA_MAPPING_ORIGINS)}, not {origin!r}"
            )
        request_id = message["requestId"]
        if origin == "notebook":
            # Only an answer to setDataMapping says "notebook", and an answer
            # names the request it answers.
            _require_exact_message_text(
                request_id,
                label="data_mapping requestId",
            )
        elif request_id is not None:
            raise ValueError(
                "data_mapping requestId must be null unless the event answers "
                "setDataMapping (origin 'notebook')"
            )
        message = {
            **message,
            "mapping": _require_data_mapping(
                message["mapping"],
                label="data_mapping mapping",
            ),
        }
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
        {key: value for key, value in message.items() if key != "type"},
    )
