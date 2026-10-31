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
import re
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
    "setGeneSet": frozenset({"type", "requestId", "name", "description", "genes"}),
}

# Every command that names cells or carries a column, which Python sends as
# one binary frame: a ``uint32`` little-endian header length, the UTF-8 JSON
# header (exactly these keys, the shared fixture's ``binaryCommands``), padded
# with spaces so the body starts on a 4-byte boundary, then the cells as
# little-endian ``uint32``, strictly ascending (a ``setField`` then its
# values). ``setVisibility`` of every cell is ``all: true`` and no cells. Two
# million cells are 8 MB written in one pass, where the same cells as JSON
# text were 14 MB to build and parse. Cells are original positions: rows of
# the dataset's ``obs``, whatever subsample the viewer opened.
_BINARY_COMMAND_HEADER_FIELDS: dict[str, frozenset[str]] = {
    "highlight": frozenset({"type", "color", "label"}),
    "setVisibility": frozenset({"type", "visible", "all"}),
    "setField": frozenset(
        {"type", "requestId", "key", "kind", "categories", "nObs", "cells", "values"}
    ),
}

# What a caller of ``send_message`` gives for a binary command: its cells under
# ``cells`` (None for every cell of a ``setVisibility``, and for a ``setField``
# whose values cover every cell), beside the header's fields; a ``setField``
# gives its values under ``values`` as a NumPy array, whose dtype the header
# names.
_BINARY_COMMAND_MESSAGE_FIELDS: dict[str, frozenset[str]] = {
    "highlight": frozenset({"type", "cells", "color", "label"}),
    "setVisibility": frozenset({"type", "cells", "visible"}),
    "setField": frozenset(
        {"type", "requestId", "key", "kind", "categories", "nObs", "cells", "values"}
    ),
}

#: The result kinds the viewer sends, with the exact keys of each ``meta``
#: (the shared fixture's ``results.kinds``).
_RESULT_META_FIELDS: dict[str, frozenset[str]] = {
    "highlight": frozenset({"label", "tool", "pageId", "pageName", "color", "enabled"}),
    "categorical": frozenset({"key", "categories", "origin"}),
    "continuous": frozenset({"key", "origin"}),
    "geneSet": frozenset({"name", "description", "genes", "symbols", "origin"}),
    "visibility": frozenset(),
    "subset": frozenset({"digest"}),
}

#: The tools a highlight group is made with.
_HIGHLIGHT_TOOLS = frozenset(
    {"annotation", "category", "combined", "knn", "lasso", "proximity", "range"}
)

#: Who made a field or gene set.
_RESULT_ORIGINS = frozenset({"analysis", "markers", "notebook", "user"})

#: The kinds of field ``setField`` pushes.
_FIELD_KINDS = frozenset({"categorical", "continuous"})

#: Each value type of a result or ``setField`` body, as a little-endian NumPy dtype.
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
    _BINARY_COMMAND_MESSAGE_FIELDS
)

# Every event the viewer may send back, mapped to the complete key set that
# event must carry. This mapping is the single declaration of what an inbound
# event is: :func:`_require_inbound_jupyter_event` validates against it and
# :func:`_protocol_capability_document` publishes its keys, so a viewer can ask
# what this installation accepts instead of guessing from a version number.
# Events carry no viewer id: the channel's socket is one viewer's.
#
# Three events do not arrive as JSON text frames: the channel builds
# ``selection`` from a binary frame (its ``cells`` field is the frame's body,
# little-endian uint32, until :func:`_decode_selection_cells` reads it),
# ``result`` from a binary frame (its header's fields and the frame's
# ``body``, read by :func:`_decode_result`) and ``session_bundle`` from a
# binary frame it streams to a temporary file.
_INBOUND_EVENT_FIELDS: dict[str, frozenset[str]] = {
    "selection": frozenset({"type", "cells", "source"}),
    "result": frozenset(
        {"type", "requestId", "epoch", "revision", "kind", "id", "nObs", "meta", "cells", "values", "body"}
    ),
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
    "command_error": frozenset({"type", "command", "reason"}),
    "data_mapping": frozenset({"type", "requestId", "origin", "mapping"}),
    "color_by": frozenset({"type", "field", "source", "symbol", "cells_read"}),
    "result_removed": frozenset({"type", "epoch", "revision", "kind", "id"}),
    "results": frozenset({"type", "epoch", "dataset", "revision", "nObs", "items"}),
}

# One page's account of one dataset it opened (viewer-protocol.json
# ``results.epochPattern``): every result event carries it, and a host keeps
# the results of one epoch.
_RESULT_EPOCH_PATTERN = re.compile(r"^[0-9a-f]{16}$")

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


def _binary_command_frame(
    header: dict[str, Any], cells: np.ndarray | None, values: np.ndarray | None = None
) -> bytearray:
    """One binary command frame: header length, padded JSON header, cells, values."""
    text = json.dumps(header, ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    header_bytes = (text + " " * (-len(text) % 4)).encode("ascii")
    body_offset = 4 + len(header_bytes)
    count = 0 if cells is None else cells.size
    values_offset = body_offset + 4 * count
    frame = bytearray(values_offset + (0 if values is None else values.nbytes))
    frame[:4] = len(header_bytes).to_bytes(4, "little")
    frame[4:body_offset] = header_bytes
    if cells is not None:
        # The cells are written straight into the frame: one pass, no copy
        # of them beside it.
        np.frombuffer(frame, dtype=_SELECTION_CELL_DTYPE, count=count, offset=body_offset)[:] = cells
    if values is not None:
        frame[values_offset:] = values.tobytes()
    return frame


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


def _require_field_values(
    values: object, *, kind: str, categories: list[str] | None, count: int, label: str
) -> tuple[str, np.ndarray]:
    """A pushed column's values: codes (a code type's largest value is no
    category) or floats (NaN is no value, never an infinity), one per cell."""
    if not isinstance(values, np.ndarray) or values.ndim != 1:
        raise TypeError(f"{label} values must be a one-dimensional NumPy array")
    dtype_name = next(
        (
            name
            for name, dtype in _VALUE_DTYPES.items()
            if values.dtype.kind == dtype.kind and values.dtype.itemsize == dtype.itemsize
        ),
        None,
    )
    if dtype_name is None:
        raise TypeError(f"{label} values must be uint8, uint16, uint32, float32 or float64, not {values.dtype}")
    if values.size != count:
        raise ValueError(f"{label} holds {values.size} {plural(values.size, 'value')} for {count} {plural(count, 'cell')}")
    if kind == "categorical":
        if dtype_name not in _CODE_TYPES:
            raise TypeError(f"{label} codes must be uint8, uint16 or uint32")
        assert categories is not None
        missing = int(np.iinfo(values.dtype).max)
        if len(categories) > missing:
            raise ValueError(f"{label} has more categories than its {dtype_name} code type holds")
        if values.size and bool(np.any((values >= len(categories)) & (values != missing))):
            raise ValueError(f"{label} has a code that names no category")
    else:
        if dtype_name not in _FLOAT_TYPES:
            raise TypeError(f"{label} values must be float32 or float64")
        if bool(np.any(np.isinf(values))):
            raise ValueError(f"{label} values must be finite numbers or NaN (no value)")
    return dtype_name, values.astype(_VALUE_DTYPES[dtype_name], copy=False)


def _set_field_frame(message: dict[str, Any]) -> bytearray:
    """One ``setField`` frame: a column for the viewer to add as a field."""
    request_id = _require_exact_message_text(message["requestId"], label="setField requestId")
    key = _require_exact_message_text(message["key"], label="setField key", allow_whitespace=True)
    kind = message["kind"]
    if kind not in _FIELD_KINDS:
        raise ValueError("setField kind must be 'categorical' or 'continuous'")
    categories = message["categories"]
    if kind == "categorical":
        categories = _require_text_list(categories, label="setField categories", nonempty=True, whitespace=True)
    elif categories is not None:
        raise TypeError("A continuous setField has no categories")
    n_obs = message["nObs"]
    if type(n_obs) is not int or not 1 <= n_obs <= _MAX_CELL_INDEX + 1:
        raise ValueError("setField nObs must be the dataset's cell count")
    raw_cells = message["cells"]
    cells = None if raw_cells is None else _require_cell_indices(raw_cells, label="setField cells")
    if cells is not None and int(cells[-1]) >= n_obs:
        raise ValueError(f"setField cells must be below the dataset's {n_obs} {plural(n_obs, 'cell')}")
    count = n_obs if cells is None else int(cells.size)
    dtype_name, values = _require_field_values(
        message["values"], kind=kind, categories=categories, count=count, label=f"setField {key!r}"
    )
    if cells is not None:
        # The cells go out ascending; each value goes with its cell.
        order = np.argsort(np.asarray(raw_cells, dtype=np.int64), kind="stable")
        values = values[order]
    header = {
        "type": "setField",
        "requestId": request_id,
        "key": key,
        "kind": kind,
        "categories": categories,
        "nObs": n_obs,
        "cells": None if cells is None else int(cells.size),
        "values": dtype_name,
    }
    return _binary_command_frame(header, cells, values)


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
        label = message["label"]
        if label is not None:
            _require_exact_message_text(label, label="highlight label", allow_whitespace=True)
        return _binary_command_frame({"type": "highlight", "color": color, "label": label}, cells)
    if message_type == "setField":
        return _set_field_frame(message)
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
    elif message_type == "setGeneSet":
        _require_exact_message_text(message["requestId"], label="setGeneSet requestId")
        _require_exact_message_text(message["name"], label="setGeneSet name", allow_whitespace=True)
        if type(message["description"]) is not str:
            raise TypeError("setGeneSet description must be text")
        _require_text_list(message["genes"], label="setGeneSet genes", nonempty=False, whitespace=False)
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
            raise ValueError("highlight result color must be a six-digit hex colour")
        if type(meta["enabled"]) is not bool:
            raise TypeError("highlight result enabled must be true or false")
    elif kind in {"categorical", "continuous"}:
        _require_exact_message_text(meta["key"], label=f"{kind} result key", allow_whitespace=True)
        if kind == "categorical":
            _require_text_list(meta["categories"], label="categorical result categories", nonempty=True, whitespace=True)
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
    request_id = message["requestId"]
    if request_id is not None:
        _require_exact_message_text(request_id, label="result requestId")
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
        "requestId": request_id,
        "epoch": message["epoch"],
        "revision": message["revision"],
        "kind": kind,
        "id": result_id,
        "nObs": n_obs,
        "meta": meta,
        "cells": cells,
        "values": values,
    }


def _require_inbound_jupyter_event(message: object) -> tuple[str, dict[str, Any]]:
    """Validate one event from the viewer channel; return its type and payload.

    The payload is the event without ``type``: what the notebook's hooks
    receive.
    """
    if type(message) is not dict:
        raise TypeError("Inbound viewer event must be a native dictionary")
    # A selection's cells and a result's body are the binary body the
    # channel received, not JSON; they are read by `_decode_selection_cells`
    # and `_decode_result` below.
    binary_key = {"selection": "cells", "result": "body"}.get(str(message.get("type")))
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
