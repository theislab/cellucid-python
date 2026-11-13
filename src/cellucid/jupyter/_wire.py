"""
Exact value checks for the viewer channel.

Every message that crosses the viewer channel (:mod:`cellucid._viewer_channel`,
one WebSocket per viewer) is checked here before it is trusted: outbound
commands are validated and serialized by :func:`_require_frontend_message` (an
``action`` whose parameters hold arrays as one binary frame, every other
command as JSON text),
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
library, NumPy (to read a binary ``selection`` body and write an action's
arrays in one vectorized pass),
and :mod:`cellucid.data_mapping`, which is standard-library only itself, so the
protocol can be reasoned about without loading the viewer machinery.
"""

from __future__ import annotations

import json
import math
import re
from typing import Any
from urllib.parse import urlparse

import numpy as np

from .._plural import plural
from ..data_mapping import DataMappingError, validate_data_mapping

# Every command Python sends to the viewer as a JSON text frame, mapped to the
# complete key set that command must carry (the shared fixture's ``commands``).
# Everything that changes or reads the viewer's state is an ``action`` of the
# catalog (``action_catalog.json``): ``action`` names it by its dotted type,
# ``params`` are its parameters, and the viewer answers with one
# ``action_result`` or ``command_error`` carrying the same ``requestId``. An
# action names the notebook turn it belongs to (``turn``, its id, and
# ``turnLabel``, what it was asked to do; both null outside a turn): the
# viewer undoes a turn's steps together.
_COMMAND_MESSAGE_FIELDS: dict[str, frozenset[str]] = {
    "ping": frozenset({"type", "requestId"}),
    "debug_snapshot": frozenset({"type", "requestId"}),
    "requestSessionBundle": frozenset({"type", "requestId"}),
    "freeze": frozenset({"type"}),
    "action": frozenset({"type", "requestId", "action", "params", "turn", "turnLabel"}),
    # An agent app on this computer driving the viewer through ``cellucid mcp``
    # (:mod:`cellucid.mcp`): its tool list, a call, a stop, and the person's
    # answer given in the app. The viewer runs each call as the agent through a
    # local connection, under the person's permissions, and answers it with
    # the agent events below.
    "agent_list": frozenset({"type", "requestId", "connection"}),
    "agent_call": frozenset({"type", "requestId", "tool", "args", "agent", "connection"}),
    "agent_cancel": frozenset({"type", "requestId"}),
    "agent_answer": frozenset({"type", "requestId", "allowed"}),
}

#: The keys of an agent command's ``connection``: the local connection's id
#: (one MCP session of one app) and the app's own name, or None.
_AGENT_CONNECTION_FIELDS = frozenset({"id", "client"})

#: A local connection's id (the shared fixture's ``agent.connectionIdPattern``).
_AGENT_CONNECTION_ID_PATTERN = re.compile(r"^[A-Za-z0-9_-]{1,64}$")

#: An agent's tool name, as the tool list names it (``color_by``).
_AGENT_TOOL_PATTERN = re.compile(r"^[a-z][a-z0-9_]*$")

# An action whose parameters hold arrays (cells as original positions, a
# column's values) is one binary frame: a ``uint32`` little-endian header
# length, the UTF-8 JSON header (exactly these keys, the shared fixture's
# ``binaryCommands``), padded with spaces so the body starts on a 4-byte
# boundary, then each array's elements, little-endian, in the order
# ``arrays`` lists them (``path``: the keys and list indices leading to it in
# ``params``, ending in a key ``params`` does not hold; ``type``; ``length``).
# Two million cells are 8 MB written in one pass, where the same cells as
# JSON text are a longer string to build and parse. ``load_session`` pushes a
# session file into the viewer: its header names the request and the file
# (``name``, without a folder), padded the same way, and its body is the
# file; the viewer answers with ``session_loaded`` or a ``command_error``.
_BINARY_COMMAND_HEADER_FIELDS: dict[str, frozenset[str]] = {
    "action": frozenset({"type", "requestId", "action", "params", "turn", "turnLabel", "arrays"}),
    "load_session": frozenset({"type", "requestId", "name"}),
}

#: The keys of one entry of a binary action's ``arrays``, in their order.
_ACTION_ARRAY_FIELDS = ("path", "type", "length")

#: The keys of an ``action_result``'s ``step``: the history step it made.
_ACTION_STEP_FIELDS = frozenset({"id", "label"})

#: Why a change made no step of the history (an ``action_result``'s
#: ``unrecorded``; None beside a step and for an action that is no change, a
#: read), each with the sentence said after the call that asked
#: (``viewer.color_by(...): ...``). The shared fixture's ``action.unrecorded``,
#: word for word (the web app's ``VIEWER_ACTION_UNRECORDED``, cellucid-r's
#: ``.viewer_action_unrecorded``); an agent is told the same reasons.
_ACTION_UNRECORDED: dict[str, str] = {
    "unchanged": "it changed nothing (it was so already), so no step was recorded.",
    "retracted": "it put back what the step just before it changed, so that step was taken out and no step was recorded.",
    "off": (
        "done, but the History settings do not record this kind of change, so no step was recorded "
        "and Undo does not take it back."
    ),
    "paused": "done, but the history's recording is paused, so no step was recorded and Undo does not take it back.",
    "elsewhere": "done, but another tracker of the history was shown meanwhile, so no step was recorded.",
    "no-history": "done, but the dataset's history has not started yet, so no step was recorded.",
}

#: A dotted action type, as the catalog names one (``color.by``).
_ACTION_TYPE_PATTERN = re.compile(r"^[a-z][a-z0-9]*(?:-[a-z0-9]+)*(?:\.[a-z][a-z0-9]*(?:-[a-z0-9]+)*)+$")

#: The result kinds the viewer sends, with the exact keys of each ``meta``
#: (the shared fixture's ``results.kinds``).
_RESULT_META_FIELDS: dict[str, frozenset[str]] = {
    "highlight": frozenset({"label", "tool", "pageId", "pageName", "color", "enabled"}),
    "categorical": frozenset({"key", "categories", "colors", "origin"}),
    "continuous": frozenset({"key", "origin"}),
    "geneSet": frozenset({"name", "description", "genes", "symbols", "origin"}),
    "visibility": frozenset(),
    "subset": frozenset({"digest"}),
}

#: A color the viewer shows a category in: ``#rrggbb`` in lowercase.
_SHOWN_COLOR = re.compile(r"#[0-9a-f]{6}")

#: The tools a highlight group is made with.
_HIGHLIGHT_TOOLS = frozenset(
    {"annotation", "category", "combined", "knn", "lasso", "proximity", "range"}
)

#: Who made a field or gene set.
_RESULT_ORIGINS = frozenset({"analysis", "markers", "notebook", "user"})

#: Each value type of a result body or an action's array, as a little-endian NumPy dtype.
_VALUE_DTYPES: dict[str, np.dtype] = {
    "uint8": np.dtype("<u1"),
    "uint16": np.dtype("<u2"),
    "uint32": np.dtype("<u4"),
    "float32": np.dtype("<f4"),
    "float64": np.dtype("<f8"),
}
_CODE_TYPES = frozenset({"uint8", "uint16", "uint32"})
_FLOAT_TYPES = frozenset({"float32", "float64"})

# Every command Python may send. This is the single declaration of what a
# notebook command is: :func:`_require_frontend_message` validates outbound
# commands against it, and an inbound ``command_error`` event may only name a
# command that appears here.
_COMMAND_TYPES: frozenset[str] = frozenset(_COMMAND_MESSAGE_FIELDS) | frozenset(
    _BINARY_COMMAND_HEADER_FIELDS
)

# Every event the viewer may send back, mapped to the complete key set that
# event must carry. This mapping is the single declaration of what an inbound
# event is: :func:`_require_inbound_jupyter_event` validates against it and
# :func:`_protocol_capability_document` publishes its keys, so a viewer can ask
# what this installation accepts instead of guessing from a version number.
# Events carry no viewer id: the channel's socket is one viewer's.
#
# Six events do not arrive as JSON text frames: the channel builds
# ``selection`` from a binary frame (its ``cells`` field is the frame's body,
# little-endian uint32, until :func:`_decode_selection_cells` reads it),
# ``result`` from a binary frame (its header's fields and the frame's
# ``body``, read by :func:`_decode_result`), ``action_result`` from a binary
# frame (its header's fields and the frame's ``body``, the action's outcome as
# UTF-8 JSON, read by :func:`_decode_action_result`: a read's result is as
# large as it asks for), ``agent_tools`` and ``agent_result`` from binary
# frames whose body is UTF-8 JSON (:func:`_decode_agent_json`), and
# ``session_bundle`` from a binary frame it streams to a temporary file.
_INBOUND_EVENT_FIELDS: dict[str, frozenset[str]] = {
    "selection": frozenset({"type", "cells", "source"}),
    "result": frozenset(
        {"type", "epoch", "revision", "kind", "id", "nObs", "meta", "cells", "values", "body"}
    ),
    "action_result": frozenset({"type", "requestId", "step", "unrecorded", "body"}),
    "hover": frozenset({"type", "cell", "position"}),
    "click": frozenset({"type", "cell", "button", "shift", "ctrl"}),
    "ready": frozenset({"type", "n_cells", "n_obs", "dimensions"}),
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
    "command_error": frozenset({"type", "command", "requestId", "reason"}),
    "data_mapping": frozenset({"type", "origin", "mapping"}),
    "color_by": frozenset({"type", "field", "source", "symbol", "cells_read"}),
    "result_removed": frozenset({"type", "epoch", "revision", "kind", "id"}),
    "results": frozenset({"type", "epoch", "dataset", "revision", "nObs", "items"}),
    "history_changed": frozenset({"type", "kind", "tracker", "head", "step", "branch"}),
    "session_loaded": frozenset({"type", "requestId", "tracker", "head"}),
    # The answers to an agent's commands (``cellucid mcp``): its tools and a
    # call's MCP result arrive as binary frames whose body is UTF-8 JSON
    # (``tools``, ``result`` once read); the person being asked about a call,
    # and their answer.
    "agent_tools": frozenset({"type", "requestId", "body"}),
    "agent_result": frozenset({"type", "requestId", "body"}),
    "agent_asking": frozenset({"type", "requestId", "message", "detail"}),
    "agent_answered": frozenset({"type", "requestId", "allowed"}),
}

#: Every agent event: delivered to the viewer's agent port
#: (:class:`cellucid._agent_link.ViewerAgentPort`), never to the notebook's hooks.
_AGENT_EVENT_TYPES = frozenset({"agent_tools", "agent_result", "agent_asking", "agent_answered"})

#: What a ``history_changed`` says changed (the shared fixture's
#: ``history.changeKinds``): a step recorded, the step shown moved, the
#: history's structure edited, a tracker started, used, renamed or removed, a
#: session file loaded.
_HISTORY_CHANGE_KINDS = ("step", "moved", "edited", "tracker", "loaded")

# One page's account of one dataset it opened (viewer-protocol.json
# ``results.epochPattern``): every result event carries it, and a host keeps
# the results of one epoch.
_RESULT_EPOCH_PATTERN = re.compile(r"^[0-9a-f]{16}$")

# Why the viewer's data mapping is what it is, as a ``data_mapping`` event
# reports it (``notebook``: this host's own ``mapping.set`` action).
_DATA_MAPPING_ORIGINS = frozenset(
    {"default", "best_guess", "session", "user", "notebook"}
)

#: The one entry an AnnData server lists: the store the viewer reads in place.
_ANNDATA_ENTRY_FIELDS = frozenset({"id", "name", "anndata"})

# Longest ``command_error`` reason the viewer may report. The viewer truncates
# its own text well below this (1000 characters: an action's refusal names
# what it expected), so the bound refuses an unbounded string without turning
# a change to the viewer's presentation limit into a protocol break.
_COMMAND_ERROR_REASON_LIMIT = 4000


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
        The commands it may send as binary frames: an ``action`` whose
        parameters hold arrays, and ``load_session`` (a session file).
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


def _binary_command_frame(header: dict[str, Any], parts: list[np.ndarray]) -> bytearray:
    """One binary command frame: header length, the JSON header padded to a
    multiple of 4 bytes, then each array's little-endian bytes in order."""
    text = json.dumps(header, ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    header_bytes = (text + " " * (-len(text) % 4)).encode("ascii")
    offset = 4 + len(header_bytes)
    frame = bytearray(offset + sum(part.nbytes for part in parts))
    frame[:4] = len(header_bytes).to_bytes(4, "little")
    frame[4:offset] = header_bytes
    for part in parts:
        # Each array is written straight into the frame: one pass, no copy
        # of it beside the frame.
        np.frombuffer(frame, dtype=part.dtype, count=part.size, offset=offset)[:] = part
        offset += part.nbytes
    return frame


def _action_array_type(array: np.ndarray, *, label: str) -> str:
    """The wire type of an action's array: its value type, or refused."""
    if array.ndim != 1:
        raise TypeError(f"{label} must be one-dimensional, not of {array.ndim} {plural(array.ndim, 'dimension')}")
    for name, dtype in _VALUE_DTYPES.items():
        if array.dtype.kind == dtype.kind and array.dtype.itemsize == dtype.itemsize:
            return name
    raise TypeError(f"{label} must be uint8, uint16, uint32, float32 or float64, not {array.dtype}")


def _split_action_params(
    params: dict[str, Any],
) -> tuple[dict[str, Any], list[tuple[list[str | int], str, np.ndarray]]]:
    """An action's parameters with every NumPy array taken out: the JSON rest,
    and each array with its path, in the order met (keys in their order). An
    array goes under an object's key, never directly in a list."""
    arrays: list[tuple[list[str | int], str, np.ndarray]] = []

    def visit(node: Any, path: list[str | int], in_list: bool) -> Any:
        if isinstance(node, np.ndarray):
            if in_list:
                where = ".".join(str(key) for key in path)
                raise TypeError(f"action params hold an array inside a list at {where}; put it under a key")
            name = _action_array_type(node, label=f"action params {'.'.join(str(key) for key in path)}")
            arrays.append((path, name, node.astype(_VALUE_DTYPES[name], copy=False)))
            return _ARRAY_TAKEN
        if type(node) is list:
            return [visit(item, [*path, index], True) for index, item in enumerate(node)]
        if type(node) is dict:
            out: dict[str, Any] = {}
            for key, child in node.items():
                if type(key) is not str:
                    raise TypeError("action params object keys must be native strings")
                kept = visit(child, [*path, key], False)
                if kept is not _ARRAY_TAKEN:
                    out[key] = kept
            return out
        return node

    rest = visit(params, [], False)
    _require_exact_json_value(rest, label="action params")
    return rest, arrays


#: The place of an array taken out of an action's parameters.
_ARRAY_TAKEN = object()


def _action_frame(message: dict[str, Any]) -> str | bytearray:
    """One ``action`` command: JSON text, or a binary frame when its
    parameters hold arrays (the viewer reads both the same way)."""
    request_id = _require_exact_message_text(message["requestId"], label="action requestId")
    action = message["action"]
    if type(action) is not str or _ACTION_TYPE_PATTERN.fullmatch(action) is None:
        raise ValueError(f"action must be a dotted action type such as 'color.by', not {action!r}")
    params = message["params"]
    if type(params) is not dict:
        raise TypeError("action params must be a native dictionary")
    turn, turn_label = message["turn"], message["turnLabel"]
    if (turn is None) != (turn_label is None):
        raise ValueError("an action names its turn and the turn's label together, or neither")
    if turn is not None:
        _require_exact_message_text(turn, label="action turn")
        _require_exact_message_text(turn_label, label="action turnLabel", allow_whitespace=True)
    rest, arrays = _split_action_params(params)
    if not arrays:
        return json.dumps(
            {"type": "action", "requestId": request_id, "action": action, "params": rest, "turn": turn, "turnLabel": turn_label},
            ensure_ascii=True,
            allow_nan=False,
            separators=(",", ":"),
        )
    header = {
        "type": "action",
        "requestId": request_id,
        "action": action,
        "params": rest,
        "turn": turn,
        "turnLabel": turn_label,
        "arrays": [{"path": path, "type": name, "length": int(array.size)} for path, name, array in arrays],
    }
    return _binary_command_frame(header, [array for _, _, array in arrays])


def _load_session_frame(message: dict[str, Any]) -> bytearray:
    """One ``load_session`` command: the header (the request and the file's
    name) padded to 4 bytes, then the session file."""
    request_id = _require_exact_message_text(message["requestId"], label="load_session requestId")
    name = _require_exact_message_text(message["name"], label="load_session name", allow_whitespace=True)
    if "/" in name or "\\" in name:
        raise ValueError("load_session names the file by its name, without a folder")
    data = message["file"]
    if type(data) not in (bytes, bytearray) or not data:
        raise TypeError("load_session carries the session file's bytes")
    return _binary_command_frame(
        {"type": "load_session", "requestId": request_id, "name": name}, [np.frombuffer(bytes(data), dtype=np.uint8)],
    )


def _require_text_list(value: object, *, label: str, nonempty: bool, whitespace: bool) -> list[str]:
    """A list of text, each once (no whitespace inside an entry unless allowed)."""
    if type(value) is not list or (nonempty and not value):
        raise TypeError(f"{label} must be a {'non-empty ' if nonempty else ''}list of text")
    seen: set[str] = set()
    for index, entry in enumerate(value):
        if type(entry) is not str or (not whitespace and (not entry or entry != entry.strip() or any(c.isspace() for c in entry))):
            raise TypeError(f"{label}[{index}] must be text")
        if entry in seen:
            raise ValueError(f"{label} must not repeat {entry!r}")
        seen.add(entry)
    return value


def _require_frontend_message(message: object) -> str | bytearray:
    """Validate one command and serialize it as the frame the viewer reads.

    An ``action`` whose parameters hold NumPy arrays becomes one binary frame
    (bytes); every other command a JSON text frame (str).
    """
    if type(message) is not dict:
        raise TypeError("message must be a native dictionary")
    message_type = _require_exact_message_text(
        message.get("type"),
        label="message type",
    )
    if "viewerId" in message or "viewerToken" in message:
        raise TypeError("message must not supply viewer routing credentials")

    if message_type == "load_session":
        if set(message) != {"type", "requestId", "name", "file"}:
            raise TypeError("load_session message must contain exactly file, name, requestId, type")
        return _load_session_frame(message)
    expected = _COMMAND_MESSAGE_FIELDS.get(message_type)
    if expected is None:
        raise ValueError(f"Unknown viewer command type {message_type!r}")
    if set(message) != expected:
        fields = ", ".join(sorted(expected))
        raise TypeError(f"{message_type} message must contain exactly {fields}")
    if message_type == "action":
        return _action_frame(message)

    _require_exact_json_value(message, label="message")
    if message_type.startswith("agent_"):
        _require_agent_command(message)
    if message_type in {"ping", "debug_snapshot", "requestSessionBundle"}:
        _require_exact_message_text(
            message["requestId"],
            label=f"{message_type} requestId",
        )
    return json.dumps(
        message,
        ensure_ascii=True,
        allow_nan=False,
        separators=(",", ":"),
    )


def _require_agent_connection(value: object, *, label: str) -> dict[str, Any]:
    if type(value) is not dict or set(value) != _AGENT_CONNECTION_FIELDS:
        raise TypeError(f"{label} connection must be {{id, client}}")
    if type(value["id"]) is not str or not _AGENT_CONNECTION_ID_PATTERN.match(value["id"]):
        raise ValueError(f'{label} connection id must be 1 to 64 letters, digits, "_" or "-"')
    client = value["client"]
    if client is not None:
        _require_exact_message_text(client, label=f"{label} connection client", allow_whitespace=True)
    return value


def _require_agent_command(message: dict[str, Any]) -> None:
    """Check an agent command's fields beyond its keys (the viewer checks them
    the same way: ``requireViewerAgentCommand``)."""
    message_type = message["type"]
    _require_exact_message_text(message["requestId"], label=f"{message_type} requestId")
    if message_type in {"agent_list", "agent_call"}:
        _require_agent_connection(message["connection"], label=message_type)
    if message_type == "agent_call":
        if type(message["tool"]) is not str or not _AGENT_TOOL_PATTERN.match(message["tool"]):
            raise ValueError("agent_call names its tool as the tool list does, such as 'color_by'")
        if type(message["args"]) is not dict:
            raise TypeError("agent_call args must be one object")
        _require_exact_message_text(message["agent"], label="agent_call agent", allow_whitespace=True)
    if message_type == "agent_answer" and type(message["allowed"]) is not bool:
        raise TypeError("agent_answer allowed must be True or False")


def _decode_agent_json(message: dict[str, Any], key: str) -> dict[str, Any]:
    """Read one binary agent event (``agent_tools``, ``agent_result``): its
    request and its body, UTF-8 JSON, as ``key``."""
    event_type = message["type"]
    request_id = _require_exact_message_text(message["requestId"], label=f"{event_type} requestId")
    body = message["body"]
    if not isinstance(body, bytes | bytearray):
        raise TypeError(f"An {event_type} body must be bytes")
    try:
        value = json.loads(bytes(body).decode("utf-8"), parse_constant=_reject_constant)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"The body of {event_type} {request_id} is not UTF-8 JSON") from error
    _require_exact_json_value(value, label=f"{event_type} body")
    if event_type == "agent_tools":
        if type(value) is not list or any(type(tool) is not dict or type(tool.get("name")) is not str for tool in value):
            raise TypeError("agent_tools body must be a list of tools, each with its name")
    elif type(value) is not dict or type(value.get("content")) is not list or type(value.get("isError")) is not bool:
        raise TypeError("agent_result body must be an MCP tool result: {content, isError, ...}")
    return {"type": event_type, "requestId": request_id, key: value}


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


def _require_epoch(value: object, *, label: str) -> str:
    if type(value) is not str or _RESULT_EPOCH_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} must be 16 lowercase hexadecimal characters")
    return value


def _require_positive_revision(value: object, *, label: str) -> int:
    if type(value) is not int or value < 1:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _require_result_meta(kind: str, meta: object) -> dict[str, Any]:
    """Check one result's ``meta`` against the rules of its kind."""
    if type(meta) is not dict:
        raise TypeError(f"A {kind} result's meta must be an object")
    _require_exact_event_fields(meta, _RESULT_META_FIELDS[kind], event_type=f"{kind} result meta")
    if "origin" in meta and meta["origin"] not in _RESULT_ORIGINS:
        raise ValueError(f"A {kind} result's origin must be one of {sorted(_RESULT_ORIGINS)}")
    if kind == "highlight":
        _require_exact_message_text(meta["label"], label="highlight result label", allow_whitespace=True)
        if meta["tool"] not in _HIGHLIGHT_TOOLS:
            raise ValueError("highlight result tool must be one of the highlight tools")
        _require_exact_message_text(meta["pageId"], label="highlight result pageId")
        _require_exact_message_text(meta["pageName"], label="highlight result pageName", allow_whitespace=True)
        color = meta["color"]
        if type(color) is not str or len(color) != 7 or color[0] != "#" or any(
            character not in "0123456789abcdefABCDEF" for character in color[1:]
        ):
            raise ValueError("highlight result color must be a six-digit hex color")
        if type(meta["enabled"]) is not bool:
            raise TypeError("highlight result enabled must be true or false")
    elif kind in {"categorical", "continuous"}:
        _require_exact_message_text(meta["key"], label=f"{kind} result key", allow_whitespace=True)
        if kind == "categorical":
            _require_text_list(meta["categories"], label="categorical result categories", nonempty=True, whitespace=True)
            # The colors it is shown in, in the view in focus: one #rrggbb a category, in their order.
            colors = meta["colors"]
            if type(colors) is not list or len(colors) != len(meta["categories"]) or any(
                type(color) is not str or not _SHOWN_COLOR.fullmatch(color) for color in colors
            ):
                raise ValueError("categorical result colors must be one #rrggbb color (lowercase) a category, in order")
    elif kind == "geneSet":
        _require_exact_message_text(meta["name"], label="gene set result name", allow_whitespace=True)
        if type(meta["description"]) is not str:
            raise TypeError("gene set result description must be text")
        _require_text_list(meta["genes"], label="gene set result genes", nonempty=False, whitespace=False)
        symbols = meta["symbols"]
        if symbols is not None and (
            type(symbols) is not list or len(symbols) != len(meta["genes"]) or any(type(symbol) is not str for symbol in symbols)
        ):
            raise ValueError("gene set result symbols must name each of its genes, in order, or be null")
    elif kind == "subset":
        digest = meta["digest"]
        if type(digest) is not str or len(digest) != 16 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("subset result digest must be 16 lowercase hex digits")
    return meta


def _decode_result(message: dict[str, Any]) -> dict[str, Any]:
    """Read one ``result`` event: its header checked against its kind, and its
    body read as ``cells`` (original positions, uint32 ascending, or None for
    every cell) and ``values`` (a NumPy array, or None)."""
    _require_epoch(message["epoch"], label="result epoch")
    _require_positive_revision(message["revision"], label="result revision")
    kind = message["kind"]
    if kind not in _RESULT_META_FIELDS:
        raise ValueError(f"Unknown result kind {kind!r}")
    result_id = _require_exact_message_text(message["id"], label="result id", allow_whitespace=True)
    n_obs = _require_nonnegative_event_integer(message["nObs"], label="result nObs")
    meta = _require_result_meta(kind, message["meta"])
    cell_count = message["cells"]
    if cell_count is not None:
        _require_nonnegative_event_integer(cell_count, label="result cells")
    value_type = message["values"]
    if value_type is not None and value_type not in _VALUE_DTYPES:
        raise ValueError("result values must be a value type or null")
    expected_values = {
        "highlight": None,
        "categorical": _CODE_TYPES,
        "continuous": _FLOAT_TYPES,
        "geneSet": None,
        "visibility": None,
        "subset": None,
    }[kind]
    if (expected_values is None) != (value_type is None) or (
        expected_values is not None and value_type not in expected_values
    ):
        raise ValueError(f"A {kind} result carries the wrong values")
    if kind == "highlight" and (cell_count is None or cell_count < 1):
        raise ValueError("A highlight result names its cells, at least one")
    if kind == "geneSet" and cell_count != 0:
        raise ValueError("A gene set result names no cells")
    if kind in {"categorical", "continuous", "subset"} and cell_count == 0:
        raise ValueError(f"A {kind} result names its cells, or every cell")
    if kind in {"visibility", "subset"} and result_id != {"visibility": "visible", "subset": "open"}[kind]:
        raise ValueError(f"A {kind} result has the id {'visible' if kind == 'visibility' else 'open'}")
    body = message["body"]
    if not isinstance(body, bytes | bytearray):
        raise TypeError("A result body must be bytes")
    cell_bytes = 0 if cell_count is None else 4 * cell_count
    value_count = 0 if value_type is None else (n_obs if cell_count is None else cell_count)
    value_bytes = 0 if value_type is None else value_count * _VALUE_DTYPES[value_type].itemsize
    if len(body) != cell_bytes + value_bytes:
        raise ValueError(
            f"The {kind} result {result_id!r} body must be {cell_bytes + value_bytes} "
            f"{plural(cell_bytes + value_bytes, 'byte')}, not {len(body)}"
        )
    cells = None
    if cell_count is not None:
        cells = np.frombuffer(body, dtype=_SELECTION_CELL_DTYPE, count=cell_count).astype(np.uint32)
        if cells.size > 1 and not bool(np.all(cells[1:] > cells[:-1])):
            raise ValueError(f"The {kind} result {result_id!r} must list its cells strictly ascending")
        if cells.size and int(cells[-1]) >= n_obs:
            raise ValueError(f"The {kind} result {result_id!r} names a cell beyond the dataset's {n_obs}")
    values = None
    if value_type is not None:
        dtype = _VALUE_DTYPES[value_type]
        values = np.frombuffer(body, dtype=dtype, count=value_count, offset=cell_bytes).astype(dtype.newbyteorder("="))
        if kind == "categorical":
            missing = int(np.iinfo(values.dtype).max)
            if len(meta["categories"]) > missing or (
                values.size and bool(np.any((values >= len(meta["categories"])) & (values != missing)))
            ):
                raise ValueError(f"The categorical result {result_id!r} has a code that names no category")
        elif values.size and bool(np.any(np.isinf(values))):
            raise ValueError(f"The continuous result {result_id!r} must hold finite values or NaN")
    return {
        "type": "result",
        "epoch": message["epoch"],
        "revision": message["revision"],
        "kind": kind,
        "id": result_id,
        "nObs": n_obs,
        "meta": meta,
        "cells": cells,
        "values": values,
    }


def _reject_constant(name: str) -> Any:
    raise ValueError(f"an action's outcome holds {name}, which JSON cannot carry")


def _decode_action_result(message: dict[str, Any]) -> dict[str, Any]:
    """Read one ``action_result`` event: the request it answers, the history
    step it made (``{id, label}`` or None), why a change made none
    (``unrecorded``: a key of :data:`_ACTION_UNRECORDED`, or None) and its
    outcome, the UTF-8 JSON of the frame's body."""
    request_id = _require_exact_message_text(message["requestId"], label="action_result requestId")
    step = message["step"]
    if step is not None:
        if type(step) is not dict:
            raise TypeError("action_result step must be an object or null")
        _require_exact_event_fields(step, _ACTION_STEP_FIELDS, event_type="action_result step")
        _require_exact_message_text(step["id"], label="action_result step id")
        _require_exact_message_text(step["label"], label="action_result step label", allow_whitespace=True)
    unrecorded = message["unrecorded"]
    if unrecorded is not None and (type(unrecorded) is not str or unrecorded not in _ACTION_UNRECORDED):
        raise ValueError(f"action_result unrecorded must be one of {', '.join(_ACTION_UNRECORDED)}, or None")
    if step is not None and unrecorded is not None:
        raise ValueError("action_result unrecorded must be None when the action made a step")
    body = message["body"]
    if not isinstance(body, bytes | bytearray):
        raise TypeError("An action_result body must be bytes")
    try:
        outcome = json.loads(bytes(body).decode("utf-8"), parse_constant=_reject_constant)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"The outcome of action request {request_id} is not UTF-8 JSON") from error
    _require_exact_json_value(outcome, label="action_result outcome")
    return {"type": "action_result", "requestId": request_id, "step": step, "unrecorded": unrecorded, "outcome": outcome}


def _require_inbound_jupyter_event(message: object) -> tuple[str, dict[str, Any]]:
    """Validate one event from the viewer channel; return its type and payload.

    The payload is the event without ``type``: what the notebook's hooks
    receive. An event nested more deeply than the interpreter follows (its
    JSON body, or its fields) is refused as malformed (``ValueError``), as
    every other event that breaks the rules is.
    """
    try:
        return _checked_inbound_event(message)
    except RecursionError as error:
        raise ValueError("Inbound viewer event is nested more deeply than this package reads") from error


def _checked_inbound_event(message: object) -> tuple[str, dict[str, Any]]:
    if type(message) is not dict:
        raise TypeError("Inbound viewer event must be a native dictionary")
    # A selection's cells and a result's body are the binary body the
    # channel received, not JSON; they are read by `_decode_selection_cells`
    # and `_decode_result` below.
    binary_key = {
        "selection": "cells",
        "result": "body",
        "action_result": "body",
        "agent_tools": "body",
        "agent_result": "body",
    }.get(str(message.get("type")))
    _require_exact_json_value(
        (
            {key: value for key, value in message.items() if key != binary_key}
            if binary_key is not None
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
    elif event_type == "result":
        message = _decode_result(message)
    elif event_type == "action_result":
        message = _decode_action_result(message)
    elif event_type == "agent_tools":
        message = _decode_agent_json(message, "tools")
    elif event_type == "agent_result":
        message = _decode_agent_json(message, "result")
    elif event_type == "agent_asking":
        _require_exact_message_text(message["requestId"], label="agent_asking requestId")
        _require_exact_message_text(message["message"], label="agent_asking message", allow_whitespace=True)
        if message["detail"] is not None:
            _require_exact_message_text(message["detail"], label="agent_asking detail", allow_whitespace=True)
    elif event_type == "agent_answered":
        _require_exact_message_text(message["requestId"], label="agent_answered requestId")
        if type(message["allowed"]) is not bool:
            raise TypeError("agent_answered allowed must be a native boolean")
    elif event_type == "result_removed":
        _require_epoch(message["epoch"], label="result_removed epoch")
        _require_positive_revision(message["revision"], label="result_removed revision")
        if message["kind"] not in _RESULT_META_FIELDS:
            raise ValueError("result_removed kind must be a result kind")
        _require_exact_message_text(message["id"], label="result_removed id", allow_whitespace=True)
    elif event_type == "results":
        _require_epoch(message["epoch"], label="results epoch")
        if message["dataset"] is not None:
            _require_exact_message_text(message["dataset"], label="results dataset")
        _require_nonnegative_event_integer(message["revision"], label="results revision")
        _require_nonnegative_event_integer(message["nObs"], label="results nObs")
        items = message["items"]
        if type(items) is not list:
            raise TypeError("results items must be a list")
        if message["dataset"] is None and items:
            raise ValueError("results of a dataset the host did not open list nothing")
        for index, item in enumerate(items):
            if type(item) is not dict:
                raise TypeError(f"results items[{index}] must be an object")
            _require_exact_event_fields(item, frozenset({"kind", "id", "revision"}), event_type="results item")
            if item["kind"] not in _RESULT_META_FIELDS:
                raise ValueError(f"results items[{index}] kind must be a result kind")
            _require_exact_message_text(item["id"], label=f"results items[{index}] id", allow_whitespace=True)
            _require_positive_revision(item["revision"], label=f"results items[{index}] revision")
    elif event_type == "color_by":
        field, source, symbol = message["field"], message["source"], message["symbol"]
        if (field is None) != (source is None):
            raise ValueError("color_by names a field and its source, or neither")
        if field is not None:
            _require_exact_message_text(field, label="color_by field", allow_whitespace=True)
            if source not in {"obs", "var"}:
                raise ValueError("color_by source must be 'obs' or 'var'")
        if symbol is not None:
            if source != "var":
                raise ValueError("color_by names a symbol only for a gene")
            _require_exact_message_text(symbol, label="color_by symbol", allow_whitespace=True)
        # A gene read for a sample of the open cells: the first ``cells_read`` of them.
        if message["cells_read"] is not None:
            if source != "var":
                raise ValueError("color_by names the cells read only for a gene")
            _require_nonnegative_event_integer(message["cells_read"], label="color_by cells_read")
    elif event_type == "ready":
        _require_nonnegative_event_integer(message["n_cells"], label="ready n_cells")
        _require_nonnegative_event_integer(message["n_obs"], label="ready n_obs")
        if message["n_cells"] > message["n_obs"]:
            raise ValueError("ready n_cells must not exceed n_obs")
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
        # The request the failed command named, or None (``freeze`` names
        # none, and a frame read no further than its type).
        if message["requestId"] is not None:
            _require_exact_message_text(message["requestId"], label="command_error requestId")
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
        message = {
            **message,
            "mapping": _require_data_mapping(
                message["mapping"],
                label="data_mapping mapping",
            ),
        }
    elif event_type == "history_changed":
        if message["kind"] not in _HISTORY_CHANGE_KINDS:
            raise ValueError(f"history_changed kind must be one of {list(_HISTORY_CHANGE_KINDS)}, not {message['kind']!r}")
        for key in ("tracker", "head"):
            _require_exact_message_text(message[key], label=f"history_changed {key}")
        for key in ("step", "branch"):
            if message[key] is not None:
                _require_exact_message_text(message[key], label=f"history_changed {key}")
    elif event_type == "session_loaded":
        _require_exact_message_text(message["requestId"], label="session_loaded requestId")
        tracker = message["tracker"]
        if type(tracker) is not dict or set(tracker) != {"id", "name"}:
            raise TypeError("session_loaded tracker must be {id, name}")
        _require_exact_message_text(tracker["id"], label="session_loaded tracker id")
        _require_exact_message_text(tracker["name"], label="session_loaded tracker name", allow_whitespace=True)
        _require_exact_message_text(message["head"], label="session_loaded head")
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
