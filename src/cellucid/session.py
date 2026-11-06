"""
Cellucid session files (``.cellucid-session``, session format 4), read in Python.

A session file holds a dataset's history: its trackers, each a tree of steps
(every action a person, an agent or a notebook took, by name) on named
branches, with groups, notes, bookmarks, turns and the journal of the
history's own edits; the value of every state domain each step left behind;
and, when the person chose to include them, the conversations kept with a
tracker. No picture of a step is ever stored. Every record is the history model's
(``history_model.schema.json``, the web app's ``app/history/model.js``; read
as :mod:`cellucid.history`). It is the web app's
``app/session/history-file.js`` written down; the framing is
``bundle/format.js``:

1) the signature ``CELLUCID_SESSION/4\n`` (the format version is part of it);
2) the manifest's length (u32 little-endian) and the manifest (UTF-8 JSON);
3) the chunks in manifest order, each its length (u32 little-endian) and its
   stored bytes.

**Manifest**: ``createdAt``, ``app`` (the Cellucid version that wrote it and
the SHA-256 of its action catalogue), ``actions`` (the action types its steps
use: one stamp per type and parameters version, with the digest of that
version's parameter schema), ``scope`` (what the share chose to include),
``datasetFingerprint`` (the dataset it was saved on), ``source`` (where that
dataset opens from, or null), ``requires`` (the obs columns and genes the
shown state uses), ``chrome`` (page chrome kept outside the history: the
pressed button of each pressed-button group, ``{pressed: {<group id>:
<button id>}}``, such as the selection tool; ``{pressed: {}}`` from a file
written without a page), ``cells`` when it was saved on a
subsample, ``genes`` when genes were read for a sample of the open cells,
``activeTracker``, ``trackers`` (each the model's ``Tracker``, its steps,
journal and conversations as chunk ids) and ``chunks``.

**Chunks**, every eager chunk before any lazy one and, within a priority,
the kinds in this order:

- ``steps`` (``steps/<n>/<k>``, eager): a run of one tracker's ``Step``
  records as JSON, in id order; a tracker lists its steps chunks, read in
  order as one list. ``before``, ``after`` and ``reads`` name domain values
  by their hash; ``params`` and ``outcome`` are skeletons.
- ``parts`` (``parts/<priority>/<k>``): ``{<name>: <skeleton>}``, sub-values
  stored once (below); eager for what the shown values name, else lazy.
- ``array`` (``array/<digest>``): one typed array a value or a step holds,
  stored once by its digest ``<type>:<length>:<sha256 of its little-endian
  bytes>``: ``raw`` elements or a ``cell-set``; eager for what the shown
  values name, else lazy.
- ``values`` (``values/<priority>/<k>``): ``{<hash>: <skeleton>}``; eager:
  the values the active head shows; lazy: every other value a step or a
  journal entry refers to. At least one eager one, except in a file shared
  as its steps only (``scope.content`` ``steps``), which holds none.
- ``journal`` (``journal/<n>``, lazy): one tracker's ``JournalEntry``
  records, their steps' params as skeletons.
- ``chat`` (``chat/<n>/<m>``, lazy): one conversation's transcript.

``steps``, ``parts`` and ``values`` are JSON cut into pieces of about 16 MiB
of text (one list or map in several chunks). Each chunk is stored gzip'd
(``codec: "gzip"``) or as it is (``"none"``, its ``storedBytes`` equal to its
``uncompressedBytes``).

A **skeleton** is the value as JSON with each typed array as ``{"$array":
<digest>}``, each number JSON cannot hold as ``{"$number": "NaN" |
"Infinity" | "-Infinity" | "-0"}``, and each object or list stored once as
a part as ``{"$part": <name>}`` (the SHA-256 of its skeleton text; a part
names the parts inside it). A value's hash is the SHA-256 of the value
joined, as canonical text (keys sorted by UTF-16 code units, numbers and
strings written as JavaScript's ``JSON.stringify`` writes them), as
``app/history/value-pool.js`` hashes it: never of its parts, so a part's name
is not trusted. A part is joined once and shared by every value that names
it; a part that names itself (through others), and a value whose text would
pass :data:`MAX_VALUE_TEXT` characters, are refused before anything is
expanded.

**Untrusted input.** :class:`CellucidSession` reads and checks the whole file
when it is made, as the web app does before it applies anything: every chunk
is inflated to exactly its declared size, every array is checked against its
digest, every value against its hash once joined, every record against the
history model and every tree whole (one root, branches that are lines from
their fork step to their tip, unique names, the head on the active branch's
path), what the share left out is absent, the values the active head shows
are in the eager chunks with every part and array they name, and a value,
part or array nothing refers to, one a step needs that the file lacks, and
trailing bytes are refused, each with a sentence that says what is wrong. A file of another
session format is refused with the web app's own sentence, and a file written
by a newer Cellucid than this package, which this package does not read,
with one that says so (:func:`newer_session_message`), never as damaged.
"""

from __future__ import annotations

import base64
import binascii
import copy
import hashlib
import json
import math
import re
import shutil
import struct
import zlib
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, BinaryIO
from urllib.parse import urlsplit

import numpy as np

from . import _history_model as model
from . import _history_rules as _rules
from ._history_model import ActionStamp, AppStamp, Requires, Scope
from ._plural import plural
from .data_mapping import DataMappingError, validate_data_mapping, validate_embedding_ref
from .history import Tracker, build_tracker
from .session_codecs import decode_cell_set, decode_uvarint

#: The session format this release reads, and its file signature.
SESSION_FORMAT_VERSION = 4
_SIGNATURE_PREFIX = b"CELLUCID_SESSION/"
#: The bytes every session file of this format starts with.
SESSION_SIGNATURE = _SIGNATURE_PREFIX + f"{SESSION_FORMAT_VERSION}\n".encode("ascii")
#: The web app's sentence for a file without a session signature
#: (``NOT_A_SESSION_FILE_MESSAGE`` in ``bundle/format.js``).
NOT_A_SESSION_FILE_MESSAGE = (
    f"This file is not a Cellucid session file of format {SESSION_FORMAT_VERSION} "
    f"(it does not start with CELLUCID_SESSION/{SESSION_FORMAT_VERSION})."
)


def other_session_format_message(version: int) -> str:
    """The web app's sentence for a session file of another format version.

    One saved by an earlier Cellucid (a saved state, where this version keeps
    the history of every step), or by a newer one (``otherSessionFormatMessage``
    in ``bundle/format.js``).
    """
    if version < SESSION_FORMAT_VERSION:
        return (
            f"This session was saved by an earlier version of Cellucid (session format {version}), "
            "which this version no longer reads: a session file now holds the history of every step "
            f"(session format {SESSION_FORMAT_VERSION}). It was not opened. Open it in the version "
            "of Cellucid that saved it, or re-create the session in this version."
        )
    return (
        f"This session was saved by a newer version of Cellucid (session format {version}), "
        "which this version cannot read. Open it in that version of Cellucid."
    )


def session_signature_refusal(magic: bytes | bytearray) -> str:
    """Why a file whose first bytes are ``magic`` is not read.

    ``magic`` is as many bytes as this format's signature holds (fewer for a
    shorter file): another session format's signature, or none
    (``sessionSignatureRefusal`` in ``bundle/format.js``).
    """
    signature = bytes(magic)
    match = (
        re.match(rb"([0-9]+)(?:\n|\Z)", signature[len(_SIGNATURE_PREFIX) :])
        if signature.startswith(_SIGNATURE_PREFIX)
        else None
    )
    version = None if match is None else int(match.group(1))
    if version is None or version in (0, SESSION_FORMAT_VERSION):
        return NOT_A_SESSION_FILE_MESSAGE
    return other_session_format_message(version)


U32_BYTES = 4
#: The most bytes one frame (the manifest, or one stored chunk) can hold: its
#: length is a u32, so this is the format's limit, not a policy. Every inflate
#: is bounded by the chunk's own declared ``uncompressedBytes``.
MAX_FRAME_BYTES = 0xFFFF_FFFF
#: The largest integer the web app reads exactly from JSON.
MAX_SAFE_INTEGER = (1 << 53) - 1
_HASH_RE = re.compile(r"[0-9a-f]{64}\Z")
_DIGEST_RE = re.compile(r"([A-Za-z0-9]+):(0|[1-9][0-9]*):([0-9a-f]{64})\Z")
_HEX16_RE = re.compile(r"[0-9a-f]{16}\Z")
_BASE64_RE = re.compile(r"[A-Za-z0-9+/]*={0,2}\Z")
_UTC_TIMESTAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z\Z")
_CONTROL_CHARACTER_RE = re.compile(r"[\x00-\x1f\x7f]")
_ARRAY_INDEX_KEY_RE = re.compile(r"(?:0|[1-9][0-9]*)\Z")
_LONE_SURROGATE_RE = re.compile("[\ud800-\udfff]")

_CHUNK_BASE_KEYS = frozenset({"id", "kind", "priority", "codec", "storedBytes", "uncompressedBytes"})
_CHUNK_KEYS_BY_KIND = {
    "steps": _CHUNK_BASE_KEYS | {"tracker"},
    "parts": _CHUNK_BASE_KEYS,
    "array": _CHUNK_BASE_KEYS | {"digest", "encoding"},
    "values": _CHUNK_BASE_KEYS,
    "journal": _CHUNK_BASE_KEYS | {"tracker"},
    "chat": _CHUNK_BASE_KEYS | {"tracker", "conversation"},
}
#: The kinds in the order a file holds them within a priority (``KIND_ORDER``).
_KIND_ORDER = ("steps", "parts", "array", "values", "journal", "chat")
_ARRAY_ENCODINGS = ("raw", "cell-set")
_FINGERPRINT_KEYS = ("datasetId", "name", "cellCount", "varCount", "cellOrder", "dataMapping")
_CELL_ORDER_KEYS = frozenset({"embedding", "dimension", "digest"})
_SOURCE_KEYS_BY_KIND = {
    "sample": frozenset({"kind", "dataset"}),
    "link": frozenset({"kind", "locator", "dataset"}),
    "remote": frozenset({"kind", "url", "dataset"}),
    "file": frozenset({"kind", "name", "dataset"}),
}
_SUBSAMPLE_METHODS = ("uniform", "blocks", "stratified", "spatial")
_SUBSAMPLE_SPEC_KEYS = {
    "uniform": (),
    "blocks": ("blockSize",),
    "stratified": ("field", "allocation", "minPerCategory"),
    "spatial": ("embedding", "columns", "gridSize"),
}
_STRATIFIED_ALLOCATIONS = ("proportional", "balanced")
_MAX_UINT32 = 0xFFFF_FFFF

#: The typed arrays a value may hold, by the name their digest starts with
#: (``ARRAY_TYPES`` in ``app/history/value-parts.js``), and the little-endian
#: numpy type of their elements.
_ARRAY_TYPES: dict[str, np.dtype[Any]] = {
    "Int8Array": np.dtype("i1"),
    "Uint8Array": np.dtype("u1"),
    "Uint8ClampedArray": np.dtype("u1"),
    "Int16Array": np.dtype("<i2"),
    "Uint16Array": np.dtype("<u2"),
    "Int32Array": np.dtype("<i4"),
    "Uint32Array": np.dtype("<u4"),
    "Float32Array": np.dtype("<f4"),
    "Float64Array": np.dtype("<f8"),
    "BigInt64Array": np.dtype("<i8"),
    "BigUint64Array": np.dtype("<u8"),
    "Float16Array": np.dtype("<f2"),
}
_ARRAY_KEY = "$array"
_NUMBER_KEY = "$number"
#: The key that names a part (a sub-value stored once) in a skeleton (``PART_KEY``).
_PART_KEY = "$part"
#: The longest skeleton text, in UTF-16 code units, one value may expand to
#: (``MAX_VALUE_TEXT`` in ``value-parts.js``: V8's longest string). A file
#: whose parts would expand a value past it is refused before it is expanded.
MAX_VALUE_TEXT = 2**29 - 24
_SPECIAL_NUMBERS = {"NaN": math.nan, "Infinity": math.inf, "-Infinity": -math.inf, "-0": -0.0}


def _damaged(what: str) -> ValueError:
    """A file that breaks the format: says what, for the person who chose it."""
    return ValueError(f"This session file is damaged: {what}")


# ---------------------------------------------------------------------------
# JSON as JavaScript reads and writes it
# ---------------------------------------------------------------------------


def _refuse_constant(name: str) -> Any:
    raise ValueError(f"{name} is not JSON")


def _parse_int(text: str) -> int | float:
    # JSON numbers are JavaScript numbers: "-0" is negative zero, and digits
    # past what an integer text may hold read as the number they name.
    if text == "-0":
        return -0.0
    try:
        return int(text)
    except ValueError:
        return float(text)


def _parse_json(data: bytes, label: str) -> Any:
    """``JSON.parse`` of UTF-8 bytes, as ``TextDecoder('utf-8', {fatal: true})``
    reads them (one leading byte-order mark dropped); NaN and Infinity are not
    JSON."""
    try:
        text = data.decode("utf-8")
        if text.startswith("﻿"):
            text = text[1:]
        return json.loads(text, parse_constant=_refuse_constant, parse_int=_parse_int)
    except (UnicodeDecodeError, ValueError, RecursionError) as error:
        raise _damaged(f"{label} is not JSON.") from error


def _utf16_key(value: str) -> bytes:
    """Sort key matching JavaScript string order (UTF-16 code units)."""
    return value.encode("utf-16-be", "surrogatepass")


def _js_string(value: str) -> str:
    """``JSON.stringify`` of a string: a lone surrogate as ``\\uXXXX``."""
    text = json.dumps(value, ensure_ascii=False)
    return _LONE_SURROGATE_RE.sub(lambda match: f"\\u{ord(match.group(0)):04x}", text)


def _js_number(value: int | float) -> str:
    """``JSON.stringify`` of a finite JavaScript number (``Number::toString``)."""
    if isinstance(value, int):
        if -(1 << 53) <= value <= 1 << 53:
            return str(value)
        value = float(value)
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    text = repr(abs(value))
    mantissa, _, exponent_text = text.partition("e")
    whole, _, fraction = mantissa.partition(".")
    all_digits = whole + fraction
    leading = len(all_digits) - len(all_digits.lstrip("0"))
    digits = all_digits.lstrip("0").rstrip("0")
    k = len(digits)
    n = len(whole) + (int(exponent_text) if exponent_text else 0) - leading
    if k <= n <= 21:
        return sign + digits + "0" * (n - k)
    if 0 < n <= 21:
        return sign + digits[:n] + "." + digits[n:]
    if -6 < n <= 0:
        return sign + "0." + "0" * (-n) + digits
    exponent = n - 1
    exponent_part = ("+" if exponent >= 0 else "-") + str(abs(exponent))
    if k == 1:
        return sign + digits + "e" + exponent_part
    return sign + digits[0] + "." + digits[1:] + "e" + exponent_part


def _is_number(value: Any) -> bool:
    return isinstance(value, int | float) and not isinstance(value, bool)


def _is_finite_number(value: Any) -> bool:
    if not _is_number(value):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


def _is_safe_integer(value: Any) -> bool:
    """``Number.isSafeInteger``: an integral JSON number of at most 2**53 - 1."""
    if isinstance(value, bool):
        return False
    if isinstance(value, int):
        return abs(value) <= MAX_SAFE_INTEGER
    return isinstance(value, float) and value.is_integer() and abs(value) <= MAX_SAFE_INTEGER


def _special_number_name(value: float) -> str | None:
    """The name of a number JSON cannot hold, or None (``specialNumberName``)."""
    if value == 0 and math.copysign(1.0, value) < 0:
        return "-0"
    if math.isnan(value):
        return "NaN"
    if math.isinf(value):
        return "Infinity" if value > 0 else "-Infinity"
    return None


def _js_key_order(keys: list[str]) -> list[str]:
    """The order JavaScript keeps an object's keys in: array indices first,
    ascending, then the rest as they were made."""
    indices = [key for key in keys if _ARRAY_INDEX_KEY_RE.match(key) and int(key) < _MAX_UINT32]
    taken = set(indices)
    return sorted(indices, key=int) + [key for key in keys if key not in taken]


def _js_stringify(value: Any) -> str:
    """``JSON.stringify`` of parsed JSON (finite numbers), keys in JavaScript's order."""
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return _js_string(value)
    if _is_number(value):
        return _js_number(value) if _is_finite_number(value) else "null"
    if isinstance(value, list):
        return "[" + ",".join(_js_stringify(item) for item in value) + "]"
    if isinstance(value, dict):
        return "{" + ",".join(
            f"{_js_string(key)}:{_js_stringify(value[key])}" for key in _js_key_order(list(value))
        ) + "}"
    raise TypeError(f"{type(value).__name__} is not JSON")


# ---------------------------------------------------------------------------
# Values: skeletons, arrays and hashes
# ---------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class _Array:
    """One typed array a value holds: its JavaScript type, its digest, its
    elements (read-only, little-endian)."""

    type: str
    digest: str
    data: np.ndarray


def _parse_array_digest(digest: Any) -> tuple[str, int]:
    """What an array digest says: its type and length (``parseArrayDigest``)."""
    match = _DIGEST_RE.match(digest) if isinstance(digest, str) else None
    if match is None or match.group(1) not in _ARRAY_TYPES or int(match.group(2)) > MAX_SAFE_INTEGER:
        raise ValueError(f'"{str(digest)[:80]}" is not the digest of an array a state value may hold')
    return match.group(1), int(match.group(2))


def _array_digest(type_name: str, data: np.ndarray) -> str:
    raw = np.ascontiguousarray(data, dtype=_ARRAY_TYPES[type_name]).tobytes()
    return f"{type_name}:{data.size}:{hashlib.sha256(raw).hexdigest()}"


def _canonical(value: Any, out: list[str]) -> None:
    """The canonical text of a value (``canonical`` in ``value-pool.js``)."""
    if value is None:
        out.append("null")
    elif isinstance(value, bool):
        out.append("true" if value else "false")
    elif isinstance(value, str):
        out.append(_js_string(value))
    elif _is_number(value):
        name = None if isinstance(value, int) else _special_number_name(value)
        out.append(_js_number(value) if name is None else f'{{"{_NUMBER_KEY}":"{name}"}}')
    elif isinstance(value, _Array):
        out.append(f'{{"{_ARRAY_KEY}":{_js_string(value.digest)}}}')
    elif isinstance(value, list):
        out.append("[")
        for index, item in enumerate(value):
            if index:
                out.append(",")
            _canonical(item, out)
        out.append("]")
    elif isinstance(value, dict):
        out.append("{")
        for index, key in enumerate(sorted(value, key=_utf16_key)):
            if index:
                out.append(",")
            out.append(_js_string(key))
            out.append(":")
            _canonical(value[key], out)
        out.append("}")
    else:
        raise TypeError(f"A state value holds a {type(value).__name__}")


def _hash_value(value: Any) -> str:
    """The content hash of one domain value (``hashValue`` in ``value-pool.js``)."""
    parts: list[str] = []
    _canonical(value, parts)
    return hashlib.sha256("".join(parts).encode("utf-8")).hexdigest()


def _placeholder_key(node: dict[str, Any], path: str) -> str | None:
    """The placeholder an object of a skeleton is (``$array``, ``$number``,
    ``$part``), or None; a malformed one is refused (``placeholderKey``)."""
    key = next((candidate for candidate in (_ARRAY_KEY, _NUMBER_KEY, _PART_KEY) if candidate in node), None)
    if key is None:
        return None
    if len(node) != 1:
        raise TypeError(f"A stored value has a malformed placeholder at {path}")
    if key == _ARRAY_KEY and not isinstance(node[_ARRAY_KEY], str):
        raise TypeError(f"A stored value names an array wrongly at {path}")
    if key == _PART_KEY and (not isinstance(node[_PART_KEY], str) or _HASH_RE.match(node[_PART_KEY]) is None):
        raise TypeError(f"A stored value names a part wrongly at {path}")
    return key


def _join(skeleton: Any, array_of: Any, part_of: Any = None, path: str = "value") -> Any:
    """The value a skeleton describes (``joinValue`` in ``value-parts.js``):
    its arrays from ``array_of(digest)``, its parts from ``part_of(name)``
    (None where no part may be named)."""
    if skeleton is None or isinstance(skeleton, bool | str):
        return skeleton
    if _is_number(skeleton):
        if not _is_finite_number(skeleton):
            raise TypeError(f"A stored value has a non-finite number at {path}")
        return skeleton
    if isinstance(skeleton, list):
        return [_join(item, array_of, part_of, f"{path}[{index}]") for index, item in enumerate(skeleton)]
    if not isinstance(skeleton, dict):
        raise TypeError(f"A stored value is not JSON at {path}")
    key = _placeholder_key(skeleton, path)
    if key == _ARRAY_KEY:
        return array_of(skeleton[_ARRAY_KEY])
    if key == _PART_KEY:
        if part_of is None:
            raise TypeError(f"A stored value names a part where none may be at {path}")
        return part_of(skeleton[_PART_KEY])
    if key == _NUMBER_KEY:
        name = skeleton[_NUMBER_KEY]
        if not isinstance(name, str) or name not in _SPECIAL_NUMBERS:
            raise TypeError(f"A stored value has an unknown number at {path}")
        return _SPECIAL_NUMBERS[name]
    return {child_key: _join(child, array_of, part_of, f"{path}.{child_key}") for child_key, child in skeleton.items()}


def _utf16_length(text: str) -> int:
    """A string's length as JavaScript counts it (UTF-16 code units)."""
    return len(text.encode("utf-16-le", "surrogatepass")) // 2


#: What ``JSON.stringify`` writes a string with escapes for.
_ESCAPED_RE = re.compile(r'[\x00-\x1f"\\\ud800-\udfff]')


def _json_string_length(text: str) -> int:
    """The length of ``JSON.stringify(text)`` in UTF-16 code units."""
    if text.isascii() and _ESCAPED_RE.search(text) is None:
        return len(text) + 2
    return _utf16_length(_js_string(text))


class _PartJoiner:
    """Joining over stored parts, untrusted input (``partJoiner`` in ``value-parts.js``).

    ``part_of(name)`` joins a part once and gives the same value to every
    value that names it. A part that names itself (through others) and a
    value whose text would pass :data:`MAX_VALUE_TEXT` are refused before
    anything is expanded: :meth:`expanded_length` counts a skeleton's text
    with each part counted once from its own text. A part's name is never
    trusted: a value is checked against its own hash once joined.
    """

    def __init__(self, skeleton_of: Any, array_of: Any) -> None:
        self._skeleton_of = skeleton_of
        self._array_of = array_of
        self._joined: dict[str, Any] = {}
        self._lengths: dict[str, int] = {}
        self._visiting: set[str] = set()
        #: The parts joined so far (a reader refuses parts nothing named).
        self.used: set[str] = set()

    def _part_length(self, name: str) -> int:
        if name in self._lengths:
            return self._lengths[name]
        if name in self._visiting:
            raise TypeError("A stored value names a part inside itself")
        self._visiting.add(name)
        length = self._text_length(self._skeleton_of(name))
        self._visiting.discard(name)
        self._lengths[name] = length
        return length

    def _text_length(self, node: Any) -> int:
        """The skeleton's ``JSON.stringify`` length, each part expanded."""
        if node is None:
            return 4
        if isinstance(node, bool):
            return 4 if node else 5
        if isinstance(node, str):
            return _json_string_length(node)
        if _is_number(node):
            return len(_js_number(node)) if _is_finite_number(node) else 4
        if isinstance(node, list):
            return 2 + max(len(node) - 1, 0) + sum(self._text_length(item) for item in node)
        if isinstance(node, dict):
            # A placeholder's shape is checked before what it names is looked for.
            reserved = _placeholder_key(node, "value")
            if reserved == _PART_KEY:
                return self._part_length(node[_PART_KEY])
            return 2 + max(len(node) - 1, 0) + sum(
                _json_string_length(key) + 1 + self._text_length(child) for key, child in node.items()
            )
        return 0

    def expanded_length(self, skeleton: Any) -> int:
        """The text length ``skeleton`` expands to; refused past :data:`MAX_VALUE_TEXT`."""
        length = self._text_length(skeleton)
        if length > MAX_VALUE_TEXT:
            raise TypeError(
                f"A stored value expands to {length} {plural(length, 'character')} of text, "
                f"more than a page can hold ({MAX_VALUE_TEXT})"
            )
        return length

    def part_of(self, name: str) -> Any:
        if name in self._joined:
            return self._joined[name]
        value = _join(self._skeleton_of(name), self._array_of, self.part_of)
        self._joined[name] = value
        self.used.add(name)
        return value

    def join(self, skeleton: Any) -> Any:
        """``skeleton`` joined, after its expanded size is checked."""
        self.expanded_length(skeleton)
        return _join(skeleton, self._array_of, self.part_of)


def _plain(value: Any) -> Any:
    """A value as this package hands it out: typed arrays as read-only numpy arrays."""
    if isinstance(value, _Array):
        return value.data
    if isinstance(value, list):
        return [_plain(item) for item in value]
    if isinstance(value, dict):
        return {key: _plain(child) for key, child in value.items()}
    return value


def _described(value: Any) -> Any:
    """A value as JSON text can hold it: each typed array as ``{"array": type,
    "length": n}``, each number JSON cannot hold as ``{"$number": name}``."""
    if isinstance(value, _Array):
        return {"array": value.type, "length": int(value.data.size)}
    if isinstance(value, float):
        name = _special_number_name(value)
        return value if name is None else {_NUMBER_KEY: name}
    if isinstance(value, list):
        return [_described(item) for item in value]
    if isinstance(value, dict):
        return {key: _described(child) for key, child in value.items()}
    return value


# ---------------------------------------------------------------------------
# The manifest
# ---------------------------------------------------------------------------


def _require_exact_keys(value: Any, keys: frozenset[str] | set[str] | tuple[str, ...], label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise _damaged(f"{label} must be an object.")
    expected = set(keys)
    missing = sorted(expected - set(value), key=_utf16_key)
    unexpected = sorted(set(value) - expected, key=_utf16_key)
    if missing or unexpected:
        details = [*([f"missing: {', '.join(missing)}"] if missing else []),
                   *([f"unexpected: {', '.join(unexpected)}"] if unexpected else [])]
        raise _damaged(
            f"{label} must contain exactly: {', '.join(sorted(expected, key=_utf16_key))} ({'; '.join(details)})."
        )
    return value


def _require_text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise _damaged(f"{label} must be a non-empty text.")
    return value


def _require_nullable_text(value: Any, label: str) -> str | None:
    return None if value is None else _require_text(value, label)


def _require_trimmed(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise _damaged(f"{label} must be a nonempty trimmed string.")
    return value


def _require_finite(value: Any, label: str) -> int | float:
    if not _is_finite_number(value):
        raise _damaged(f"{label} must be a finite number.")
    number: int | float = value
    return number


def _require_safe_integer(value: Any, label: str, *, minimum: int = 0, maximum: int = MAX_SAFE_INTEGER) -> int:
    if not _is_safe_integer(value) or not minimum <= value <= maximum:
        raise _damaged(f"{label} must be a safe integer from {minimum} through {maximum}.")
    return int(value)


def _require_json(value: Any, label: str) -> Any:
    """A JSON value with finite numbers only (what a manifest field and a transcript hold)."""
    if value is None or isinstance(value, bool | str):
        return value
    if _is_number(value):
        return _require_finite(value, label)
    if isinstance(value, list):
        for index, item in enumerate(value):
            _require_json(item, f"{label}[{index}]")
        return value
    if isinstance(value, dict):
        for key, child in value.items():
            _require_json(child, f"{label}.{key}")
        return value
    raise _damaged(f"{label} must be JSON.")


def _check_created_at(value: Any) -> None:
    valid = isinstance(value, str) and _UTC_TIMESTAMP_RE.match(value) is not None
    if valid:
        try:
            datetime.strptime(value, "%Y-%m-%dT%H:%M:%S.%fZ")
        except ValueError:
            valid = False
    if not valid:
        raise _damaged("its manifest createdAt must be a canonical ISO time.")


def _validate_fingerprint(value: Any) -> dict[str, Any]:
    """Check one dataset fingerprint (``assertDatasetFingerprint``) and return a copy.

    ``dataMapping`` names the embeddings, gene names, graph, expression matrix
    and latent space the session was saved with, in canonical form;
    ``cellOrder`` names the one embedding whose stored values the viewer
    digested, each cell by its original position (the web app's
    ``data/cell-order-digest.js``). Its shape is checked here; the viewer
    checks the digest when it restores a session, and an apply checks it when
    the object has that embedding (``anndata_session/_apply.py``
    ``cell_order_digest``).
    """
    label = "Session manifest datasetFingerprint"
    fingerprint = _require_exact_keys(value, _FINGERPRINT_KEYS, label)
    for key in ("datasetId", "name"):
        if fingerprint[key] is not None:
            _require_trimmed(fingerprint[key], f"{label} {key}")
    _require_safe_integer(fingerprint["cellCount"], f"{label} cellCount")
    _require_safe_integer(fingerprint["varCount"], f"{label} varCount")
    try:
        mapping = validate_data_mapping(fingerprint["dataMapping"], label=f"{label} dataMapping")
    except DataMappingError as error:
        raise _damaged(str(error)) from error
    if _js_stringify(mapping) != _js_stringify(fingerprint["dataMapping"]):
        raise _damaged(f"{label} dataMapping must be in canonical form.")
    order_label = f"{label} cellOrder"
    cell_order = _require_exact_keys(fingerprint["cellOrder"], _CELL_ORDER_KEYS, order_label)
    dimension = _require_safe_integer(cell_order["dimension"], f"{order_label} dimension", minimum=1, maximum=3)
    digest = _require_trimmed(cell_order["digest"], f"{order_label} digest")
    if not _HEX16_RE.match(digest):
        raise _damaged(f"{order_label} digest must be 16 lowercase hexadecimal characters.")
    try:
        embedding = validate_embedding_ref(cell_order["embedding"], dimension=dimension, label=f"{order_label} embedding")
    except DataMappingError as error:
        raise _damaged(str(error)) from error
    out = copy.deepcopy(fingerprint)
    out["cellCount"] = int(fingerprint["cellCount"])
    out["varCount"] = int(fingerprint["varCount"])
    out["dataMapping"] = mapping
    out["cellOrder"] = {"embedding": embedding, "dimension": dimension, "digest": digest}
    return out


def _validate_source(value: Any) -> dict[str, str] | None:
    """Check where the session's dataset opens from (``assertDatasetSourceRecord``)."""
    label = "Session manifest source"
    if value is None:
        return None
    kind = value.get("kind") if isinstance(value, dict) else None
    keys = _SOURCE_KEYS_BY_KIND.get(kind) if isinstance(kind, str) else None
    if keys is None:
        raise _damaged(f'{label} kind must be "sample", "link", "remote" or "file".')
    _require_exact_keys(value, keys, label)
    for key in ("dataset", "locator", "name", "url"):
        if key not in keys:
            continue
        text = value[key]
        if (
            not isinstance(text, str)
            or not text
            or text != text.strip()
            or _CONTROL_CHARACTER_RE.search(text) is not None
        ):
            raise _damaged(f"{label} {key} must be a non-empty, unpadded string.")
    if kind == "remote":
        url = urlsplit(value["url"])
        if url.scheme.lower() not in ("http", "https") or not url.netloc:
            raise _damaged(f"{label} url must be an http(s) address.")
    return dict(value)


def _validate_requires(value: Any) -> dict[str, list[str]]:
    """Check ``requires`` (``assertSessionRequires``): each list sorted as
    JavaScript sorts strings, unique, of trimmed keys."""
    label = "Session manifest requires"
    _require_exact_keys(value, ("obs", "var"), label)
    out: dict[str, list[str]] = {}
    for source in ("obs", "var"):
        keys = value[source]
        if not isinstance(keys, list):
            raise _damaged(f"{label} {source} must be an array of field keys.")
        for index, key in enumerate(keys):
            _require_trimmed(key, f"{label} {source}[{index}]")
            if index > 0 and not _utf16_key(keys[index - 1]) < _utf16_key(key):
                raise _damaged(f"{label} {source} must be sorted and without duplicates.")
        out[source] = list(keys)
    return out


def _require_uint32(value: Any, label: str, *, minimum: int = 0) -> int:
    if not _is_safe_integer(value) or not minimum <= value <= _MAX_UINT32:
        raise _damaged(f"{label} must be an integer in [{minimum}, 2^32 - 1].")
    return int(value)


def _validate_subsample_spec(spec: Any, label: str) -> dict[str, Any]:
    """Check a subsample specification (``requireSubsampleSpec`` in ``data/subsample.js``)."""
    if not isinstance(spec, dict):
        raise _damaged(f"{label} must be an object.")
    method = spec.get("method")
    if method not in _SUBSAMPLE_METHODS:
        raise _damaged(f"{label} names an unknown subsample method.")
    allowed = ("method", "nObs", "k", "seed", *_SUBSAMPLE_SPEC_KEYS[method])
    _require_exact_keys(spec, allowed, label)
    out: dict[str, Any] = {
        "method": method,
        "nObs": _require_uint32(spec["nObs"], f"{label} nObs"),
        "k": _require_uint32(spec["k"], f"{label} k", minimum=1),
        "seed": _require_uint32(spec["seed"], f"{label} seed"),
    }
    if out["k"] >= out["nObs"]:
        raise _damaged(f"{label} keeps every cell; it is no subsample.")
    if method == "blocks":
        out["blockSize"] = _require_uint32(spec["blockSize"], f"{label} blockSize", minimum=1)
    if method == "stratified":
        if not isinstance(spec["field"], str) or not spec["field"]:
            raise _damaged(f"{label} must name its obs field.")
        if spec["allocation"] not in _STRATIFIED_ALLOCATIONS:
            raise _damaged(f"{label} allocation must be proportional or balanced.")
        out["field"] = spec["field"]
        out["allocation"] = spec["allocation"]
        out["minPerCategory"] = _require_uint32(spec["minPerCategory"], f"{label} minPerCategory")
    if method == "spatial":
        if not isinstance(spec["embedding"], str) or not spec["embedding"]:
            raise _damaged(f"{label} must name its embedding.")
        columns = spec["columns"]
        if (
            not isinstance(columns, list)
            or not 1 <= len(columns) <= 3
            or not all(_is_safe_integer(column) and column >= 0 for column in columns)
        ):
            raise _damaged(f"{label} must name one to three embedding columns.")
        out["embedding"] = spec["embedding"]
        out["columns"] = [int(column) for column in columns]
        out["gridSize"] = _require_uint32(spec["gridSize"], f"{label} gridSize", minimum=1)
    return out


def _forgiving_base64(text: str) -> bytes | None:
    """``atob``: padding optional, a length of one more than a multiple of four refused."""
    if len(text) % 4 == 0 and text.endswith("="):
        text = text[:-2] if text.endswith("==") else text[:-1]
    if len(text) % 4 == 1 or "=" in text:
        return None
    try:
        return base64.b64decode(text + "=" * (-len(text) % 4), validate=True)
    except (binascii.Error, ValueError):
        return None


def _decode_cells(value: Any, cell_count: int) -> dict[str, Any]:
    """The ``cells`` record of a session saved on a subsample, checked and
    decoded (``decodeCellSubsetRecord``): its specification, the kept cells'
    original positions (ascending, int64) and their cell-set digest."""
    from .subsample import cell_set_digest
    from .subsample import decode_cell_set as decode_subsample_cells

    label = "Session manifest cells"
    record = _require_exact_keys(value, ("spec", "kept", "digest"), label)
    spec = _validate_subsample_spec(record["spec"], f"{label} spec")
    kept_text = record["kept"]
    if not isinstance(kept_text, str) or _BASE64_RE.match(kept_text) is None:
        raise _damaged(f"{label} kept must be base64 text.")
    digest = record["digest"]
    if not isinstance(digest, str) or _HEX16_RE.match(digest) is None:
        raise _damaged(f"{label} digest must be 16 lowercase hex digits.")
    raw = _forgiving_base64(kept_text)
    if raw is None:
        raise _damaged(f"{label} kept must be base64 text.")
    try:
        kept = decode_subsample_cells(raw, spec["nObs"])
    except ValueError as error:
        raise _damaged(f"{label} kept: {error}.") from error
    count = spec["nObs"] if kept is None else int(kept.size)
    if count != spec["k"]:
        raise _damaged(f"{label} lists {count} {plural(count, 'cell')}; its specification keeps {spec['k']}.")
    assert kept is not None
    if cell_set_digest(spec["nObs"], kept) != digest:
        raise _damaged(f"{label} digest does not match its cells.")
    if spec["k"] != cell_count:
        raise _damaged("its manifest cells must keep exactly the fingerprint cell count.")
    return {"spec": spec, "kept": kept, "digest": digest}


def _decode_genes(value: Any, cell_count: int) -> dict[str, Any]:
    """The ``genes`` record (``decodeGeneSampleRecord``): genes were read for
    the first ``cells`` open cells in their order (``order`` names a
    ``blocks`` cell order, or null for a subsample's priority order), whose
    original positions have the cell-set digest ``digest``."""
    label = "Session manifest genes"
    record = _require_exact_keys(value, ("cells", "order", "digest"), label)
    cells = record["cells"]
    if not _is_safe_integer(cells) or not 1 <= cells < cell_count:
        raise _damaged(f"{label} cells must be a count of open cells from 1 to {cell_count - 1}.")
    order = record["order"]
    if order is not None:
        _require_exact_keys(order, ("method", "seed", "blockSize"), f"{label} order")
        if order["method"] != "blocks":
            raise _damaged(f"{label} order must be the priority order of a blocks subsample.")
        seed, block_size = order["seed"], order["blockSize"]
        if not _is_safe_integer(seed) or not 0 <= seed <= _MAX_UINT32:
            raise _damaged(f"{label} order seed must be a uint32.")
        if not _is_safe_integer(block_size) or block_size < 1:
            raise _damaged(f"{label} order blockSize must be a positive count.")
        order = {"method": "blocks", "seed": int(seed), "blockSize": int(block_size)}
    digest = record["digest"]
    if not isinstance(digest, str) or _HEX16_RE.match(digest) is None:
        raise _damaged(f"{label} digest must be 16 lowercase hex digits.")
    return {"cells": int(cells), "order": order, "digest": digest}


def _inflate_exact(stored: bytes, expected_bytes: int, chunk_id: str) -> bytes:
    """One gzip member inflated to exactly ``expected_bytes``, never more."""
    decompressor = zlib.decompressobj(wbits=16 + zlib.MAX_WBITS)
    output = bytearray()
    remaining_input = stored
    try:
        while True:
            output.extend(decompressor.decompress(remaining_input, expected_bytes - len(output) + 1))
            if len(output) > expected_bytes:
                break
            if decompressor.unconsumed_tail:
                remaining_input = decompressor.unconsumed_tail
                continue
            break
        if len(output) <= expected_bytes:
            output.extend(decompressor.flush(expected_bytes - len(output) + 1))
    except zlib.error as error:
        raise _damaged(f'the chunk "{chunk_id}" is not a gzip stream ({error}).') from error
    if len(output) > expected_bytes or (decompressor.eof and len(output) != expected_bytes):
        raise _damaged(f'the chunk "{chunk_id}" holds another size than its manifest says.')
    if not decompressor.eof:
        raise _damaged(f'the chunk "{chunk_id}" is a gzip stream that ends early.')
    if decompressor.unused_data:
        raise _damaged(f'the chunk "{chunk_id}" holds bytes after its gzip stream.')
    return bytes(output)


class _Frames:
    """The file's bytes, read in order, each read checked against what is left."""

    def __init__(self, handle: BinaryIO, size: int, position: int, digest: Any) -> None:
        self._handle = handle
        self._size = size
        self.position = position
        #: The SHA-256 of every byte read so far (the signature included by the caller).
        self._digest = digest

    @property
    def remaining(self) -> int:
        return self._size - self.position

    def read(self, count: int, what: str) -> bytes:
        if count > self.remaining:
            raise _damaged(
                f"{what} needs {count} {plural(count, 'byte')} and only {self.remaining} remain "
                "(is the file truncated?)."
            )
        data = self._handle.read(count)
        if len(data) != count:
            raise _damaged(f"{what} needs {count} {plural(count, 'byte')} and the file ended (is it truncated?).")
        self.position += count
        self._digest.update(data)
        return data

    def u32(self, what: str) -> int:
        return int(struct.unpack("<I", self.read(U32_BYTES, what))[0])


# ---------------------------------------------------------------------------
# Newer files
# ---------------------------------------------------------------------------

_VERSION_RE = re.compile(r"(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)")


def _version_tuple(text: str) -> tuple[int, int, int] | None:
    """A Cellucid version's numbers (``0.9.7``, ``0.9.7.dev1`` -> ``(0, 9, 7)``), or None."""
    match = _VERSION_RE.match(text)
    return None if match is None else (int(match.group(1)), int(match.group(2)), int(match.group(3)))


def package_version() -> str:
    """This cellucid package's version."""
    from importlib import metadata

    return metadata.version("cellucid")


def newer_session_message(version: str, package: str) -> str:
    """Why a session file written by a newer Cellucid than this package is not read.

    The same sentence in cellucid-python and cellucid-r: a file whose records
    this package does not know, written by a newer web app, is newer, never
    "damaged".
    """
    return (
        f"This session was written by a newer Cellucid ({version}) than this package reads ({package}); "
        "update the cellucid package to read it."
    )


def _newer_than_package(manifest: Any) -> str | None:
    """The writing version when the manifest names one newer than this package, else None."""
    app = manifest.get("app") if isinstance(manifest, dict) else None
    version = app.get("version") if isinstance(app, dict) else None
    if not isinstance(version, str):
        return None
    written = _version_tuple(version)
    current = _version_tuple(package_version())
    return version if written is not None and current is not None and written > current else None


# ---------------------------------------------------------------------------
# The manifest
# ---------------------------------------------------------------------------

#: The manifest's header keys (``cells`` and ``genes`` too, for a subsample).
_HEADER_KEYS = ("createdAt", "app", "actions", "scope", "datasetFingerprint", "source", "requires", "chrome")
#: A manifest tracker: the model's ``Tracker``, its steps, journal and conversations as chunk ids.
_TRACKER_KEYS = tuple(key for key, _, _, _ in model.FIELDS["Tracker"])


def _check_tracker_record(record: Any, index: int) -> None:
    """A manifest tracker's own keys and chunk ids (its records are checked with its steps)."""
    label = f"tracker {index}"
    _require_exact_keys(record, _TRACKER_KEYS, f"Session manifest {label}")
    _require_text(record["id"], f"{label} id")
    if not isinstance(record["steps"], list) or not record["steps"]:
        raise _damaged(f"{label} steps must list its steps chunks.")
    for position, chunk_id in enumerate(record["steps"]):
        _require_text(chunk_id, f"{label} steps[{position}]")
    if record["journal"] is not None:
        _require_text(record["journal"], f"{label} journal")
    if not isinstance(record["chats"], list):
        raise _damaged(f"{label} chats must be a list of chunk ids.")
    for position, chunk_id in enumerate(record["chats"]):
        _require_text(chunk_id, f"{label} chats[{position}]")


def _check_chunk_record(entry: Any, index: int) -> None:
    label = f"chunk {index}"
    kind = entry.get("kind") if isinstance(entry, dict) else None
    keys = _CHUNK_KEYS_BY_KIND.get(kind) if isinstance(kind, str) else None
    if keys is None:
        raise _damaged(f"{label} kind must be one of {', '.join(_KIND_ORDER)}.")
    _require_exact_keys(entry, keys, f"Session manifest {label}")
    _require_text(entry["id"], f"{label} id")
    if entry["priority"] not in ("eager", "lazy"):
        raise _damaged(f"{label} priority must be eager or lazy.")
    if entry["codec"] not in ("gzip", "none"):
        raise _damaged(f"{label} codec must be gzip or none.")
    if not _is_safe_integer(entry["storedBytes"]) or not 0 <= entry["storedBytes"] <= MAX_FRAME_BYTES:
        raise _damaged(f"{label} storedBytes must be a frame length.")
    if not _is_safe_integer(entry["uncompressedBytes"]) or entry["uncompressedBytes"] < 0:
        raise _damaged(f"{label} uncompressedBytes must be a byte count.")
    if entry["codec"] == "none" and entry["storedBytes"] != entry["uncompressedBytes"]:
        raise _damaged(f"{label} is stored as it is, at another size.")
    if kind in ("steps", "chat", "journal"):
        _require_text(entry["tracker"], f"{label} tracker")
    if kind == "steps" and entry["priority"] != "eager":
        raise _damaged(f"{label}: steps are read first.")
    if kind in ("chat", "journal") and entry["priority"] != "lazy":
        raise _damaged(f"{label}: a {kind} is read last.")
    if kind == "chat":
        _require_text(entry["conversation"], f"{label} conversation")
    if kind == "array":
        try:
            type_name, _ = _parse_array_digest(entry["digest"])
        except ValueError as error:
            raise _damaged(f"{label}: {error}.") from error
        if entry["encoding"] not in _ARRAY_ENCODINGS or (
            entry["encoding"] == "cell-set" and type_name != "Uint32Array"
        ):
            raise _damaged(f"{label} encoding must be raw, or cell-set for a Uint32Array.")


def _stable_id(value: Any) -> bool:
    return isinstance(value, str) and value != "" and value == value.strip()


def _check_chrome(chrome: Any) -> dict[str, Any]:
    """The page chrome a session carries outside its history (``assertPageChrome``):
    ``{pressed: {<group id>: <button id>}}``, the pressed button of each of the
    web app's pressed-button groups (the selection tool), which no step
    records; ``{pressed: {}}`` from a file written without a page."""
    if not isinstance(chrome, dict) or list(chrome) != ["pressed"]:
        raise _damaged("its page chrome must be exactly {pressed}.")
    if not isinstance(chrome["pressed"], dict):
        raise _damaged("its page chrome must map pressed-button groups to their pressed button.")
    for group, button in chrome["pressed"].items():
        # The web app's rule for an id it puts back (``utils/stable-id.js``): nonempty, no white space around it.
        if not _stable_id(group) or not _stable_id(button):
            raise _damaged(f'its page chrome must name the pressed button of the group "{group}" by its id.')
    return chrome


def _check_model_record(name: str, value: Any, label: str) -> Any:
    """One header record against the history model (``app``, ``scope``, an action stamp)."""
    try:
        _rules.check_record(name, value, label)
    except _rules.HistoryModelError as error:
        raise _damaged(f"{error}.") from error
    return value


@dataclass(frozen=True)
class _Header:
    fingerprint: dict[str, Any]
    source: dict[str, str] | None
    requires: dict[str, list[str]]
    cells: dict[str, Any] | None
    genes: dict[str, Any] | None


def _check_manifest(manifest: Any) -> _Header:
    """Check a whole manifest before any chunk is read (``assertHistoryManifest``):
    its dataset part, what wrote it and what it shares, its trackers' chunk
    ids, and that its chunks are exactly what the trackers and the share name,
    eager before lazy."""
    if not isinstance(manifest, dict):
        raise _damaged("its manifest must be a JSON object.")
    header_keys = [*_HEADER_KEYS, *(["cells"] if "cells" in manifest else []), *(["genes"] if "genes" in manifest else [])]
    _require_exact_keys(manifest, (*header_keys, "activeTracker", "trackers", "chunks"), "Session manifest")
    _check_created_at(manifest["createdAt"])
    _check_model_record("AppStamp", manifest["app"], "session.app")
    if not isinstance(manifest["actions"], list):
        raise _damaged("session.actions must be a list.")
    for index, stamp in enumerate(manifest["actions"]):
        _check_model_record("ActionStamp", stamp, f"session.actions[{index}]")
    scope = _check_model_record("Scope", manifest["scope"], "session.scope")
    fingerprint = _validate_fingerprint(manifest["datasetFingerprint"])
    source = _validate_source(manifest["source"])
    if scope["source"] == "name" and source is not None:
        raise _damaged("it was shared with the dataset's name only (scope.source), and names where it opens from.")
    requires = _validate_requires(manifest["requires"])
    _check_chrome(manifest["chrome"])
    cells = _decode_cells(manifest["cells"], fingerprint["cellCount"]) if "cells" in manifest else None
    genes = _decode_genes(manifest["genes"], fingerprint["cellCount"]) if "genes" in manifest else None
    trackers, chunks = manifest["trackers"], manifest["chunks"]
    if not isinstance(trackers, list) or not trackers:
        raise _damaged("it holds no tracker.")
    if not isinstance(chunks, list):
        raise _damaged("its manifest has no chunk list.")
    for index, record in enumerate(trackers):
        _check_tracker_record(record, index)
    for index, entry in enumerate(chunks):
        _check_chunk_record(entry, index)
    tracker_ids: set[str] = set()
    for record in trackers:
        if record["id"] in tracker_ids:
            raise _damaged(f'it holds the tracker "{record["id"]}" twice.')
        tracker_ids.add(record["id"])
        if record["journal"] is not None and not scope["journal"]:
            raise _damaged(f'tracker "{record["id"]}" keeps its journal, which the share left out (scope.journal).')
        if record["chats"] and not scope["chats"]:
            raise _damaged(f'tracker "{record["id"]}" keeps conversations, which the share left out (scope.chats).')
    if not isinstance(manifest["activeTracker"], str) or manifest["activeTracker"] not in tracker_ids:
        raise _damaged("its active tracker is not one it holds.")
    chunk_by_id: dict[str, dict[str, Any]] = {}
    last = (0, 0)
    eager_values = 0
    digests: set[str] = set()
    for entry in chunks:
        if entry["id"] in chunk_by_id:
            raise _damaged(f'it lists the chunk "{entry["id"]}" twice.')
        chunk_by_id[entry["id"]] = entry
        at = (0 if entry["priority"] == "eager" else 1, _KIND_ORDER.index(entry["kind"]))
        if at < last:
            raise _damaged(f'its chunk "{entry["id"]}" is out of order (eager before lazy; {", ".join(_KIND_ORDER)}).')
        last = at
        if entry["kind"] == "values" and entry["priority"] == "eager":
            eager_values += 1
        if entry["kind"] == "array":
            if entry["digest"] in digests:
                raise _damaged("it stores one array twice.")
            digests.add(entry["digest"])
    if scope["content"] == "steps":
        if any(entry["kind"] == "values" for entry in chunks):
            raise _damaged("it was shared as its steps only (scope.content), and holds values.")
    elif eager_values == 0:
        raise _damaged("it holds no values of the step it shows.")
    named: set[str] = set()
    position = {entry["id"]: index for index, entry in enumerate(chunks)}
    for record in trackers:
        for chunk_id in record["steps"]:
            steps = chunk_by_id.get(chunk_id)
            if steps is None or steps["kind"] != "steps" or steps["tracker"] != record["id"] or chunk_id in named:
                raise _damaged(f'tracker "{record["id"]}" names a steps chunk that is not its own.')
            named.add(chunk_id)
        # A tracker's steps are read as one list, in the order it names its chunks: the file holds them in that order.
        if any(position[record["steps"][index]] < position[record["steps"][index - 1]] for index in range(1, len(record["steps"]))):
            raise _damaged(f'tracker "{record["id"]}" names its steps chunks in another order than the file holds them.')
        if record["journal"] is not None:
            journal = chunk_by_id.get(record["journal"])
            if journal is None or journal["kind"] != "journal" or journal["tracker"] != record["id"] or record["journal"] in named:
                raise _damaged(f'tracker "{record["id"]}" names a journal chunk that is not its own.')
            named.add(record["journal"])
        conversations: set[str] = set()
        for chunk_id in record["chats"]:
            chat = chunk_by_id.get(chunk_id)
            if (
                chat is None
                or chat["kind"] != "chat"
                or chat["tracker"] != record["id"]
                or chat["conversation"] in conversations
                or chunk_id in named
            ):
                raise _damaged(f'tracker "{record["id"]}" names a conversation chunk that is not its own.')
            conversations.add(chat["conversation"])
            named.add(chunk_id)
    for entry in chunks:
        if entry["kind"] in ("steps", "chat", "journal") and entry["id"] not in named:
            raise _damaged(f'its chunk "{entry["id"]}" belongs to no tracker.')
    return _Header(fingerprint=fingerprint, source=source, requires=requires, cells=cells, genes=genes)


# ---------------------------------------------------------------------------
# Reading a file
# ---------------------------------------------------------------------------


def _walk_times(name: str, value: Any, visit: Any, path: str) -> None:
    """Call ``visit(time, path)`` for every time field of a record and the records in it."""
    if value is None:
        return
    for key, _, spec, _ in model.FIELDS[name]:
        if key not in value:
            continue
        item = value[key]
        kind = spec["type"]
        if kind == "time":
            visit(item, f"{path}.{key}")
        elif kind == "record" and item is not None and spec["name"] not in model.MAP_RECORDS:
            _walk_times(spec["name"], item, visit, f"{path}.{key}")
        elif kind == "array" and spec["item"]["type"] == "record" and item is not None:
            for index, entry in enumerate(item):
                _walk_times(spec["item"]["name"], entry, visit, f"{path}.{key}[{index}]")


_DAY_MS = 86_400_000


@dataclass(frozen=True)
class _Contents:
    manifest: dict[str, Any]
    header: _Header
    #: Every tracker as the model's JSON records (params and outcome joined).
    trackers: list[dict[str, Any]]
    values: dict[str, Any]
    sha256: str


class _Reading:
    """One file being read: what its chunks gave so far, checked as the web app checks it.

    Chunks are taken in the file's order; once they all are, :meth:`finish`
    joins the values the active head shows from the eager chunks alone (as
    the web app shows the head before it reads its lazy chunks), then every
    other value and every step's ``params`` and ``outcome`` (which may name
    arrays and parts only the lazy chunks hold).
    """

    def __init__(self, manifest: dict[str, Any]) -> None:
        self.manifest = manifest
        self.arrays: dict[str, _Array] = {}
        #: The digests and part names the lazy chunks hold.
        self.lazy_arrays: set[str] = set()
        self.lazy_parts: set[str] = set()
        self.referenced_arrays: set[str] = set()
        self.parts: dict[str, Any] = {}
        self.values: dict[str, Any] = {}
        #: Every stored value's skeleton by its hash, and those of the eager chunks.
        self.value_skeletons: dict[str, Any] = {}
        self.eager_hashes: list[str] = []
        self.step_chunks: dict[str, list[Any]] = {record["id"]: [] for record in manifest["trackers"]}
        self.journal_chunks: dict[str, Any] = {}
        self.chats: dict[str, dict[str, Any]] = {}
        #: While the eager values are joined, a part or an array only the lazy chunks hold is refused.
        self.eager_only = False
        self.joiner = _PartJoiner(self.skeleton_of, self.array_of)

    def array_of(self, digest: str) -> _Array:
        array = self.arrays.get(digest)
        if array is None:
            raise _damaged(f"a value names an array the file does not hold ({digest[:40]}…).")
        if self.eager_only and digest in self.lazy_arrays:
            raise _damaged(f"a value of the step it shows names an array that is not among its first chunks ({digest[:40]}…).")
        self.referenced_arrays.add(digest)
        return array

    def skeleton_of(self, name: str) -> Any:
        if name not in self.parts:
            raise _damaged(f"a value names a part the file does not hold ({name[:12]}…).")
        if self.eager_only and name in self.lazy_parts:
            raise _damaged(f"a value of the step it shows names a part that is not among its first chunks ({name[:12]}…).")
        return self.parts[name]

    def join(self, skeleton: Any, label: str) -> Any:
        try:
            return self.joiner.join(skeleton)
        except TypeError as error:
            raise _damaged(f"{label}: {error}.") from error

    def take_array(self, meta: dict[str, Any], payload: bytes) -> None:
        type_name, length = _parse_array_digest(meta["digest"])
        dtype = _ARRAY_TYPES[type_name]
        if meta["encoding"] == "cell-set":
            try:
                bound, offset = decode_uvarint(payload, 0)
                data = decode_cell_set(payload[offset:], cell_count=bound)
            except (TypeError, ValueError) as error:
                raise _damaged(f'the array chunk "{meta["id"]}" is not a cell set ({error}).') from error
            data = data.astype(dtype)
        else:
            if len(payload) != length * dtype.itemsize:
                raise _damaged(f'the array chunk "{meta["id"]}" has the wrong length.')
            data = np.frombuffer(payload, dtype=dtype)
        if data.size != length or _array_digest(type_name, data) != meta["digest"]:
            raise _damaged(f'the array chunk "{meta["id"]}" does not match its digest.')
        data.flags.writeable = False
        self.arrays[meta["digest"]] = _Array(type=type_name, digest=meta["digest"], data=data)
        if meta["priority"] == "lazy":
            self.lazy_arrays.add(meta["digest"])

    @staticmethod
    def take_map(target: dict[str, Any], record: Any, label: str) -> list[str]:
        """A ``parts`` or ``values`` chunk's ``{<sha256>: <skeleton>}``, each name once in the file."""
        if not isinstance(record, dict):
            raise _damaged(f"{label} must map names to skeletons.")
        for name, skeleton in record.items():
            if _HASH_RE.match(name) is None or name in target:
                raise _damaged(f"{label} holds a name twice or malformed.")
            target[name] = skeleton
        return list(record)

    def take(self, meta: dict[str, Any], payload: bytes) -> None:
        kind = meta["kind"]
        if kind == "steps":
            label = f'the steps of tracker "{meta["tracker"]}"'
            steps = _parse_json(payload, label)
            if not isinstance(steps, list):
                raise _damaged(f"{label} must be a list.")
            # A tracker's steps chunks are pieces of one list, in the file's order (the order it names them).
            self.step_chunks[meta["tracker"]].extend(steps)
        elif kind == "parts":
            names = self.take_map(self.parts, _parse_json(payload, "its parts"), "its parts")
            if meta["priority"] == "lazy":
                self.lazy_parts.update(names)
        elif kind == "journal":
            self.journal_chunks[meta["tracker"]] = _parse_json(payload, f'the journal of tracker "{meta["tracker"]}"')
        elif kind == "array":
            self.take_array(meta, payload)
        elif kind == "values":
            names = self.take_map(self.value_skeletons, _parse_json(payload, "its values"), "its values")
            if meta["priority"] == "eager":
                self.eager_hashes.extend(names)
        else:
            label = f'the conversation "{meta["conversation"]}"'
            self.chats[meta["id"]] = {"conversation": meta["conversation"], "transcript": _require_json(_parse_json(payload, label), label)}

    def join_values(self, hashes: Any, label: str) -> None:
        """Join each value and check it against its hash (a part's name is never trusted)."""
        for value_hash in hashes:
            value = self.join(self.value_skeletons[value_hash], label)
            if _hash_value(value) != value_hash:
                raise _damaged(f"{label} holds a value that does not match its hash.")
            self.values[value_hash] = value

    def joined_step(self, step: Any, label: str) -> Any:
        """A stored step with its params and outcome joined to their arrays and parts."""
        if not isinstance(step, dict):
            raise _damaged(f"{label} must be a Step record.")
        out = dict(step)
        for key in ("params", "outcome"):
            if key in step:
                out[key] = self.join(step[key], label)
        return out

    def joined_patch(self, patch: Any, label: str) -> Any:
        if not isinstance(patch, dict) or not isinstance(patch.get("steps"), list):
            return patch
        return {**patch, "steps": [self.joined_step(step, f"{label}.steps[{index}]") for index, step in enumerate(patch["steps"])]}

    def trackers(self) -> list[dict[str, Any]]:
        """Every tracker as the model's records: its steps, journal and conversations read in."""
        out = []
        for record in self.manifest["trackers"]:
            stored = self.step_chunks[record["id"]]
            steps = [self.joined_step(step, f'tracker "{record["id"]}" step {index}') for index, step in enumerate(stored)]
            journal: list[Any] = []
            if record["journal"] is not None:
                entries = self.journal_chunks.get(record["id"])
                if not isinstance(entries, list):
                    raise _damaged(f'the journal of tracker "{record["id"]}" must be a list.')
                for index, entry in enumerate(entries):
                    label = f'tracker "{record["id"]}" journal entry {index}'
                    if not isinstance(entry, dict):
                        raise _damaged(f"{label} must be a JournalEntry record.")
                    journal.append({**entry, "before": self.joined_patch(entry.get("before"), f"{label}.before"),
                                    "after": self.joined_patch(entry.get("after"), f"{label}.after")})
            chats = [[self.chats[chunk_id]["conversation"], self.chats[chunk_id]["transcript"]] for chunk_id in record["chats"]]
            out.append({**record, "steps": steps, "journal": journal, "chats": chats})
        return out

    def finish(self, header: _Header) -> list[dict[str, Any]]:
        """Every check that needs the whole file: the model's rules over every
        tracker, the share's promises, and every value, part and array
        accounted for."""
        scope = self.manifest["scope"]
        values_shared = scope["content"] != "steps"
        eager = set(self.eager_hashes)
        self.eager_only = True
        self.join_values(self.eager_hashes, "its eager values")
        self.eager_only = False
        self.join_values([value_hash for value_hash in self.value_skeletons if value_hash not in eager], "its other values")
        trackers = self.trackers()
        session = {
            "createdAt": self.manifest["createdAt"], "app": self.manifest["app"], "actions": self.manifest["actions"],
            "scope": scope, "datasetFingerprint": self.manifest["datasetFingerprint"], "source": self.manifest["source"],
            "requires": self.manifest["requires"], "cells": self.manifest.get("cells"), "genes": self.manifest.get("genes"),
            "activeTracker": self.manifest["activeTracker"], "trackers": trackers,
        }
        try:
            _rules.validate_session(session)
        except _rules.HistoryModelError as error:
            raise _damaged(f"{error}.") from error
        self._check_share(trackers, scope)
        referenced: set[str] = set()
        for tracker in trackers:
            for step in tracker["steps"]:
                referenced.update(step["before"].values())
                referenced.update(step["after"].values())
            for entry in tracker["journal"]:
                if entry["purged"]:
                    continue
                for patch in (entry["before"], entry["after"]):
                    for step in patch["steps"]:
                        referenced.update(step["before"].values())
                        referenced.update(step["after"].values())
        if values_shared:
            active = next(tracker for tracker in trackers if tracker["id"] == self.manifest["activeTracker"])
            for value_hash in _values_at({step["id"]: step for step in active["steps"]}, active["headId"]).values():
                if value_hash not in eager:
                    raise _damaged("the values of the step it shows are not among its first chunks.")
            if any(value_hash not in self.values for value_hash in referenced):
                raise _damaged("a step refers to a value the file does not hold.")
            if any(value_hash not in referenced for value_hash in self.values):
                raise _damaged("it holds a value no step refers to.")
        if any(name not in self.joiner.used for name in self.parts):
            raise _damaged("it holds a part nothing refers to.")
        if any(digest not in self.referenced_arrays for digest in self.arrays):
            raise _damaged("it holds an array nothing refers to.")
        return trackers

    @staticmethod
    def _check_share(trackers: list[dict[str, Any]], scope: dict[str, Any]) -> None:
        """What the share left out is absent; its times are as it chose them."""
        for tracker in trackers:
            label = f'tracker "{tracker["id"]}"'

            def visit(time: Any, path: str, label: str = label) -> None:
                if scope["times"] == "none" and time is not None:
                    raise _damaged(f"{label} holds a time at {path}, and the share left times out (scope.times).")
                # A note never edited has no edit time, whatever the share kept.
                if scope["times"] != "none" and time is None and not path.endswith(".editedAt"):
                    raise _damaged(f"{label} has no time at {path}, and the share kept times (scope.times).")
                if scope["times"] == "dates" and time is not None and time % _DAY_MS != 0:
                    raise _damaged(f"{label} holds a time of day at {path}, and the share kept dates only (scope.times).")

            _walk_times("Tracker", tracker, visit, "tracker")


def _values_at(by_id: dict[str, dict[str, Any]], step_id: str) -> dict[str, str]:
    """The value hash of every domain at ``step_id`` (the after-values on the path from the root)."""
    path: list[dict[str, Any]] = []
    step: dict[str, Any] | None = by_id[step_id]
    while step is not None:
        path.append(step)
        step = None if step["parent"] is None else by_id[step["parent"]]
    hashes: dict[str, str] = {}
    for item in reversed(path):
        hashes.update(item["after"])
    return hashes


def _read_contents(path: Path) -> _Contents:
    size = path.stat().st_size
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        magic = handle.read(len(SESSION_SIGNATURE))
        digest.update(magic)
        if len(magic) < len(SESSION_SIGNATURE):
            raise ValueError(NOT_A_SESSION_FILE_MESSAGE)
        if magic != SESSION_SIGNATURE:
            raise ValueError(session_signature_refusal(magic))
        frames = _Frames(handle, size, len(SESSION_SIGNATURE), digest)
        manifest_length = frames.u32("its manifest length")
        if manifest_length == 0:
            raise _damaged("its manifest length is 0 bytes.")
        if manifest_length > frames.remaining:
            raise _damaged(
                f"its manifest length of {manifest_length} {plural(manifest_length, 'byte')} is more than the "
                "file holds (is it truncated?)."
            )
        manifest = _parse_json(frames.read(manifest_length, "its manifest"), "its manifest")
        newer = _newer_than_package(manifest)
        try:
            if not isinstance(manifest, dict):
                raise _damaged("its manifest must be a JSON object.")
            if not isinstance(manifest.get("chunks"), list):
                raise _damaged("its manifest chunks must be an array.")
            header = _check_manifest(manifest)
            reading = _Reading(manifest)
            for meta in manifest["chunks"]:
                chunk_id = meta["id"]
                stored_length = frames.u32(f'the chunk "{chunk_id}"')
                if stored_length != meta["storedBytes"]:
                    raise _damaged(
                        f'the chunk "{chunk_id}" is framed as {stored_length} {plural(stored_length, "byte")} and its '
                        f"manifest says {meta['storedBytes']}."
                    )
                stored = frames.read(stored_length, f'the chunk "{chunk_id}"')
                if meta["codec"] == "gzip":
                    payload = _inflate_exact(stored, int(meta["uncompressedBytes"]), chunk_id)
                elif len(stored) != meta["uncompressedBytes"]:
                    raise _damaged(f'the chunk "{chunk_id}" holds another size than its manifest says.')
                else:
                    payload = stored
                reading.take(meta, payload)
            if frames.remaining != 0 or handle.read(1):
                raise _damaged("it holds bytes after its last chunk.")
            trackers = reading.finish(header)
        except ValueError as error:
            if newer is not None:
                raise ValueError(newer_session_message(newer, package_version())) from error
            raise
    return _Contents(
        manifest=manifest,
        header=header,
        trackers=trackers,
        values=reading.values,
        sha256=digest.hexdigest(),
    )


# ---------------------------------------------------------------------------
# The public handle
# ---------------------------------------------------------------------------


class CellucidSession:
    """A ``.cellucid-session`` file (session format 4): its trackers, branches, steps and states.

    Making one reads the whole file and checks all of it (see the module
    notes); a file of another session format, or one that breaks the format,
    raises ``ValueError`` with a sentence that says why (a file written by a
    newer Cellucid than this package says so, never that it is damaged), and
    nothing about it is read lazily later. :func:`read_session` makes one.

    Every record is the history model's (:mod:`cellucid.history`): the
    session's header (``created_at``, ``app``, ``actions``, ``scope``,
    ``dataset_fingerprint``, ``source``, ``requires``, ``chrome``, ``cells``, ``genes``,
    ``active_tracker``) and its ``trackers``, each a
    :class:`~cellucid.history.Tracker` with its branches, its steps (a table,
    and each a :class:`~cellucid.history.Step`), its turns, notes and journal.

    Parameters
    ----------
    path : str or Path
        The ``.cellucid-session`` file; ``~`` is expanded. A missing file
        raises ``FileNotFoundError``.

    Examples
    --------
    >>> session = cellucid.read_session("annotation.cellucid-session")
    >>> tracker = session.tracker()          # the tracker the viewer showed
    >>> tracker.steps[["id", "branch", "label", "origin"]]
    >>> tracker.branch("Main").steps         # the path to its tip
    >>> tracker.step("s5").state()["gene-sets"]
    >>> plan = session.apply(adata)          # what applying the head would change
    >>> report = plan.apply()
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path).expanduser().resolve()
        if not self.path.is_file():
            raise FileNotFoundError(str(self.path))
        self._contents = _read_contents(self.path)
        manifest = self._contents.manifest
        self.created_at: str = manifest["createdAt"]
        """Canonical ISO time of the save (UTC)."""
        self.app: AppStamp = model.build("AppStamp", manifest["app"])
        """What wrote the file: the Cellucid version and the SHA-256 of its action catalogue."""
        self.actions: tuple[ActionStamp, ...] = tuple(model.build("ActionStamp", stamp) for stamp in manifest["actions"])
        """The action types its steps use, each with the digest of its parameter schema, sorted by type."""
        self.scope: Scope = model.build("Scope", manifest["scope"])
        """What the share or save chose to include (trackers, branches, content, results, notes, chats, journal, times, source)."""
        self.active_tracker: str = manifest["activeTracker"]
        """The id of the tracker the session shows."""
        header = self._contents.header
        subsample = header.cells
        facts = {
            "name": header.fingerprint["name"],
            "cells": header.fingerprint["cellCount"] if subsample is None else subsample["spec"]["nObs"],
            "genes": header.fingerprint["varCount"],
            "citation": None,
            "open_cells": None if subsample is None else {
                "open": subsample["spec"]["k"], "of": subsample["spec"]["nObs"],
                "method": subsample["spec"]["method"], "seed": subsample["spec"]["seed"],
            },
            "version": self.app.version,
        }
        self.trackers: tuple[Tracker, ...] = tuple(
            build_tracker(
                record, value=self._value, values_included=self.scope.content != "steps", dataset=facts,
            )
            for record in self._contents.trackers
        )
        """Every tracker, in the file's order (one in Recently deleted has ``removed`` set)."""

    # -- the header ---------------------------------------------------------

    @property
    def sha256(self) -> str:
        """The SHA-256 of the file's bytes (what an apply records as ``session_sha256``)."""
        return self._contents.sha256

    @property
    def manifest(self) -> dict[str, Any]:
        """The manifest as the file holds it (a copy)."""
        return copy.deepcopy(self._contents.manifest)

    @property
    def dataset_fingerprint(self) -> dict[str, Any]:
        """The dataset the session was saved on: ``datasetId``, ``name``,
        ``cellCount`` (the cells open), ``varCount``, ``cellOrder`` and
        ``dataMapping`` (a copy)."""
        return copy.deepcopy(self._contents.header.fingerprint)

    @property
    def source(self) -> dict[str, str] | None:
        """Where the session's dataset was opened from, or ``None``.

        ``{"kind": "sample", "dataset"}``, ``{"kind": "link", "locator",
        "dataset"}`` (``locator`` is the ``?data=`` value), ``{"kind":
        "remote", "url", "dataset"}`` or ``{"kind": "file", "name",
        "dataset"}`` (the local file or folder name); ``None`` for a session
        saved in a notebook, or shared with the dataset's name only.
        """
        source = self._contents.header.source
        return None if source is None else dict(source)

    @property
    def requires(self) -> Requires:
        """The dataset fields the shown state uses: ``obs`` columns and ``var``
        genes (named under the session's data mapping) a view shows or filters
        on, or a highlight group is bound to, by the key the dataset gives it.
        Applying another step checks that step's own (:meth:`apply`)."""
        requires = self._contents.header.requires
        return Requires(obs=tuple(requires["obs"]), var=tuple(requires["var"]))

    @property
    def cells(self) -> dict[str, Any] | None:
        """The subsample the session was saved on, or ``None`` for every cell.

        ``{"spec": {...}, "kept": numpy.ndarray, "digest": str}``: the
        subsample's specification (``method``, ``nObs``, ``k``, ``seed`` and
        the method's own keys), the kept cells' original positions (rows of
        ``obs``, ascending, int64) and their cell-set digest. The same cells
        are ``cellucid.subsample_cells(...)`` of the specification.
        """
        cells = self._contents.header.cells
        if cells is None:
            return None
        return {"spec": copy.deepcopy(cells["spec"]), "kept": cells["kept"].copy(), "digest": cells["digest"]}

    @property
    def genes(self) -> dict[str, Any] | None:
        """The open cells genes were read for, or ``None`` when every open cell's were.

        ``{"cells": int, "order": dict | None, "digest": str}``: genes were
        read for the first ``cells`` open cells in their order (a subsample's
        priority order, or, every cell being open, the order of the ``blocks``
        subsample ``order`` names); ``digest`` is the cell-set digest of their
        original positions.
        """
        genes = self._contents.header.genes
        return None if genes is None else copy.deepcopy(genes)

    @property
    def chrome(self) -> dict[str, Any]:
        """The page chrome the session carries outside its history (a copy).

        ``{"pressed": {<group id>: <button id>}}``: the pressed button of
        each of the web app's pressed-button groups, which it puts back after
        a load and no step records, such as the selection tool
        (``{"highlight-mode": "highlight-mode-lasso"}``). ``{"pressed": {}}``
        for a file written without a page. It is not part of the data:
        :meth:`apply` ignores it.
        """
        return copy.deepcopy(self._contents.manifest["chrome"])

    # -- the history --------------------------------------------------------

    def tracker(self, tracker: str | None = None) -> Tracker:
        """A tracker by its name or its id (``None``: the one the session shows).

        Names are compared by their key (case and extra white space aside)
        among the trackers not removed; a removed tracker is named by its id.
        An unknown name is refused, naming the trackers there are.
        """
        if tracker is None:
            return next(item for item in self.trackers if item.id == self.active_tracker)
        if not isinstance(tracker, str):
            raise TypeError("tracker must be a tracker's name or id (a string) or None")
        found = next((item for item in self.trackers if item.id == tracker), None)
        if found is None:
            key = _rules.name_key(tracker)
            found = next((item for item in self.trackers if item.removed is None and _rules.name_key(item.name) == key), None)
        if found is None:
            raise ValueError(_rules.no_tracker(tracker, [item.name for item in self.trackers if item.removed is None]))
        return found

    def conversations(self, tracker: str | Tracker | None = None) -> list[dict[str, Any]]:
        """The conversations kept with a tracker (its name, id, the tracker,
        or ``None``: the one shown), in its order: ``{"conversation": id,
        "transcript": JSON}`` each. A file holds them only when the person
        chose to include them."""
        found = tracker if isinstance(tracker, Tracker) else self.tracker(tracker)
        return [{"conversation": conversation, "transcript": copy.deepcopy(transcript)} for conversation, transcript in found.chats]

    def _value(self, value_hash: str) -> Any:
        """One value as the file holds it (typed arrays as :class:`_Array`)."""
        return self._contents.values[value_hash]

    # -- files and objects --------------------------------------------------

    def save(self, dest: str | Path) -> Path:
        """Copy the session file to ``dest`` and return the destination path.

        The copy is written beside ``dest`` and moved into place once
        complete, so an interrupted save never leaves a partial file there;
        saving a session onto its own file leaves it as it is.
        """
        dest_path = Path(dest).expanduser().resolve()
        if dest_path == self.path:
            return dest_path
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        partial = dest_path.with_name(f".{dest_path.name}.partial")
        try:
            shutil.copyfile(self.path, partial)
            partial.replace(dest_path)
        except BaseException:
            partial.unlink(missing_ok=True)
            raise
        return dest_path

    def apply(
        self,
        adata: Any,
        *,
        tracker: str | None = None,
        branch: str | None = None,
        step: str | None = None,
        dry_run: bool = True,
        overwrite: bool = False,
        keep_session: bool = False,
    ) -> Any:
        """What applying one step of this session to ``adata`` would change, or the change made.

        Same as :func:`~cellucid.apply_cellucid_session_to_anndata` with this
        session; see it for the parameters, the plan's rules and the errors.
        """
        from .anndata_session import apply_cellucid_session_to_anndata

        return apply_cellucid_session_to_anndata(
            adata, self, tracker=tracker, branch=branch, step=step, dry_run=dry_run, overwrite=overwrite,
            keep_session=keep_session,
        )

    def __repr__(self) -> str:
        count = len(self.trackers)
        return f"CellucidSession({str(self.path)!r}: {count} {plural(count, 'tracker')}, active {self.active_tracker!r})"


def read_session(path: str | Path) -> CellucidSession:
    """Read a ``.cellucid-session`` file whole and check it (:class:`CellucidSession`).

    Examples
    --------
    >>> session = cellucid.read_session("annotation.cellucid-session")
    >>> session.tracker().steps.head()
    """
    return CellucidSession(path)
