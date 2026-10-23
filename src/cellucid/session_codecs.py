"""
Internal codecs used by `.cellucid-session` bundles (session format 2).

These mirror the web app's codecs exactly:

- `cellucid/assets/js/app/session/codecs/varint.js` (unsigned LEB128)
- `cellucid/assets/js/app/session/codecs/cell-set.js` (highlight memberships)
- `cellucid/assets/js/app/session/contributors/user-defined-codes.js`
  (user-defined categorical codes)
- `cellucid/assets/js/app/session/codecs/typed-json.js` (analysis settings
  and results)

Every decoder is vectorized with numpy: a membership of half a million cells is
decoded without a Python loop per cell.
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np

MAX_SAFE_INTEGER = (1 << 53) - 1
MAX_UINT32 = (1 << 32) - 1
MAX_UINT16 = (1 << 16) - 1

CELL_SET_SORTED_DELTA = 0
CELL_SET_BITSET = 1
CELL_SET_RUNS = 2
CELL_SET_REFERENCE_MASK = 3
CELL_SET_REFERENCE_ALL = 4

_CODES_RAW = 0
_CODES_RUNS = 1
_CODES_PACKED = 2

#: Each user-defined code width: the little-endian dtype and its largest value,
#: keyed by the JavaScript typed-array name the bundle declares.
_CODE_WIDTHS: dict[str, tuple[str, int]] = {
    "Uint8Array": ("u1", 0xFF),
    "Uint16Array": ("<u2", MAX_UINT16),
    "Uint32Array": ("<u4", MAX_UINT32),
}
_NATIVE_CODE_DTYPES: dict[str, type[np.unsignedinteger]] = {
    "Uint8Array": np.uint8,
    "Uint16Array": np.uint16,
    "Uint32Array": np.uint32,
}

#: Typed-JSON dtype codes (`codecs/typed-json.js`).
_TYPED_JSON_DTYPES: dict[int, str] = {
    1: "u1",
    2: "<u2",
    3: "<u4",
    4: "<i4",
    5: "<f4",
    6: "<f8",
}
_TYPED_KEY = "$typed"


def _as_bytes(data: bytes | bytearray | memoryview, label: str) -> bytes:
    if not isinstance(data, bytes | bytearray | memoryview):
        raise TypeError(f"{label}: data must be bytes-like")
    return bytes(data)


def decode_uvarint(data: bytes | bytearray | memoryview, offset: int = 0) -> tuple[int, int]:
    """
    Decode one canonical unsigned LEB128 varint.

    Returns
    -------
    value, next_offset
    """
    if not isinstance(data, bytes | bytearray | memoryview):
        raise TypeError("decode_uvarint: data must be bytes-like")
    if isinstance(offset, bool) or not isinstance(offset, int):
        raise TypeError("decode_uvarint: offset must be an integer")
    if offset < 0 or offset >= len(data):
        raise ValueError("decode_uvarint: offset is outside the payload")

    value = 0
    shift = 0
    idx = offset
    length = len(data)
    encoded_bytes = 0

    while True:
        if idx >= length:
            raise ValueError("decode_uvarint: truncated")
        b = data[idx]
        idx += 1
        encoded_bytes += 1

        value |= (b & 0x7F) << shift
        if (b & 0x80) == 0:
            if value > MAX_SAFE_INTEGER:
                raise ValueError("decode_uvarint: value exceeds the exact integer range")
            minimum_bytes = max(1, (value.bit_length() + 6) // 7)
            if encoded_bytes != minimum_bytes:
                raise ValueError("decode_uvarint: noncanonical encoding")
            return value, idx

        shift += 7
        if encoded_bytes >= 8 or shift > 49:
            raise ValueError("decode_uvarint: varint too long")


def _decode_uvarint_array(region: np.ndarray, count: int, *, label: str) -> np.ndarray:
    """Decode exactly ``count`` canonical uvarints filling ``region`` (uint8)."""
    if count == 0:
        if region.size != 0:
            raise ValueError(f"{label}: trailing bytes")
        return np.empty(0, dtype=np.uint64)
    ends = np.flatnonzero((region & 0x80) == 0)
    if ends.size > count or (ends.size == count and ends[-1] != region.size - 1):
        raise ValueError(f"{label}: trailing bytes")
    if ends.size < count:
        raise ValueError(f"{label}: truncated (expected {count} varints, found {ends.size})")
    starts = np.empty(count, dtype=np.int64)
    starts[0] = 0
    starts[1:] = ends[:-1] + 1
    lengths = ends - starts + 1
    if int(lengths.max()) > 7:
        raise ValueError(f"{label}: varint too long")
    if np.any((lengths > 1) & (region[ends] == 0)):
        raise ValueError(f"{label}: noncanonical varint encoding")
    position = np.arange(region.size, dtype=np.int64) - np.repeat(starts, lengths)
    contributions = (region.astype(np.uint64) & np.uint64(0x7F)) << (
        position.astype(np.uint64) * np.uint64(7)
    )
    values = np.add.reduceat(contributions, starts)
    if values.size and int(values.max()) > MAX_SAFE_INTEGER:
        raise ValueError(f"{label}: value exceeds the exact integer range")
    return values.astype(np.uint64)


def digest_cell_set(indices: np.ndarray) -> int:
    """
    Order-independent 32-bit digest of a set of cell indices.

    The sum, modulo 2**32, of one avalanche-mixed word per index -- the web
    app's ``digestCellSet`` (`codecs/cell-set.js`), computed vectorized with
    32-bit wrap-around.
    """
    x = np.asarray(indices, dtype=np.uint32)
    with np.errstate(over="ignore"):
        x = (x ^ (x >> np.uint32(16))) * np.uint32(0x45D9F3B)
        x = (x ^ (x >> np.uint32(16))) * np.uint32(0x45D9F3B)
        x = x ^ (x >> np.uint32(16))
    return int(x.astype(np.uint64).sum() % (1 << 32))


def cell_set_needs_reference(data: bytes | bytearray | memoryview) -> bool:
    """Whether a cell-set payload is stored relative to a category's cells."""
    payload = _as_bytes(data, "cell_set_needs_reference")
    if not payload:
        raise ValueError("cell_set_needs_reference: empty payload")
    return payload[0] in (CELL_SET_REFERENCE_MASK, CELL_SET_REFERENCE_ALL)


def _bitset_members(bits: np.ndarray, count: int, *, label: str) -> np.ndarray:
    flags = np.unpackbits(bits, bitorder="little")
    if np.any(flags[count:]):
        raise ValueError(f"{label}: bits are set past its length")
    return np.flatnonzero(flags[:count])


def decode_cell_set(
    data: bytes | bytearray | memoryview,
    *,
    cell_count: int,
    reference: np.ndarray | None = None,
) -> np.ndarray:
    """
    Decode one highlight membership (``highlights/cells/*``) into a strictly
    increasing ``uint32`` array.

    Byte 0 is the encoding: 0 sorted-delta, 1 bitset, 2 runs, 3 reference-mask,
    4 reference-all. Encodings 3 and 4 need ``reference``: the cells of the
    group's category, ascending; its count and digest are checked against the
    ones the writer recorded.
    """
    payload = _as_bytes(data, "decode_cell_set")
    if isinstance(cell_count, bool) or not isinstance(cell_count, int) or cell_count < 0:
        raise TypeError("decode_cell_set: cell_count must be a non-negative integer")
    if not payload:
        raise ValueError("decode_cell_set: empty payload")
    raw = np.frombuffer(payload, dtype=np.uint8)
    encoding = payload[0]

    if encoding in (CELL_SET_REFERENCE_MASK, CELL_SET_REFERENCE_ALL):
        if reference is None:
            raise ValueError(
                "decode_cell_set: this cell set is stored relative to a category; "
                "its cells are required"
            )
        ref = np.asarray(reference, dtype=np.uint32)
        reference_count, offset = decode_uvarint(payload, 1)
        if offset + 4 > len(payload):
            raise ValueError("decode_cell_set: truncated")
        digest = int.from_bytes(payload[offset : offset + 4], "little")
        if reference_count != ref.size or digest != digest_cell_set(ref):
            raise ValueError(
                "The category this selection was saved relative to no longer holds "
                "the same cells."
            )
        body = offset + 4
        if encoding == CELL_SET_REFERENCE_ALL:
            if body != len(payload):
                raise ValueError("decode_cell_set: trailing bytes")
            return ref.copy()
        if len(payload) - body != (reference_count + 7) // 8:
            raise ValueError("decode_cell_set: reference mask has the wrong length")
        members = _bitset_members(raw[body:], reference_count, label="decode_cell_set mask")
        selected: np.ndarray = ref[members]
        return selected.astype(np.uint32)

    if encoding == CELL_SET_BITSET:
        if len(payload) != 1 + (cell_count + 7) // 8:
            raise ValueError(
                "decode_cell_set: bitset has the wrong length for the dataset cell count"
            )
        return _bitset_members(raw[1:], cell_count, label="decode_cell_set bitset").astype(
            np.uint32
        )

    if encoding == CELL_SET_SORTED_DELTA:
        count, offset = decode_uvarint(payload, 1)
        if count > cell_count:
            raise ValueError(
                f"decode_cell_set: {count} cells exceed the dataset's {cell_count}"
            )
        deltas = _decode_uvarint_array(raw[offset:], count, label="decode_cell_set")
        if count > 1 and np.any(deltas[1:] == 0):
            raise ValueError("decode_cell_set: indices must be strictly increasing")
        indices = np.cumsum(deltas, dtype=np.uint64)
        if count and int(indices[-1]) >= cell_count:
            raise ValueError("decode_cell_set: index exceeds the dataset cell count")
        return indices.astype(np.uint32)

    if encoding == CELL_SET_RUNS:
        run_count, offset = decode_uvarint(payload, 1)
        pairs = _decode_uvarint_array(raw[offset:], 2 * run_count, label="decode_cell_set runs")
        gaps = pairs[0::2].astype(np.int64)
        lengths = pairs[1::2].astype(np.int64)
        if run_count and (np.any(lengths < 1) or np.any(gaps[1:] < 1)):
            raise ValueError("decode_cell_set: runs must be non-empty and separated")
        ends = np.cumsum(gaps + lengths)
        starts = ends - lengths
        if run_count and int(ends[-1]) > cell_count:
            raise ValueError("decode_cell_set: run exceeds the dataset cell count")
        total = int(lengths.sum())
        offsets = np.arange(total, dtype=np.int64) - np.repeat(np.cumsum(lengths) - lengths, lengths)
        return (np.repeat(starts, lengths) + offsets).astype(np.uint32)

    raise ValueError(f"decode_cell_set: encoding {encoding} is not one this version reads")


def decode_user_defined_codes(
    data: bytes | bytearray | memoryview,
    *,
    expected_length: int,
    expected_codes_type: str,
) -> np.ndarray:
    """
    Decode a user-defined categorical codes chunk (pre-gzip payload).

    Layout: ``u8 encoding``, ``uvarint length``, then

    - 0 raw: ``length`` little-endian codes at the declared width;
    - 1 runs: ``uvarint runCount``, then ``runCount`` × (uvarint value,
      uvarint run); runs are positive and adjacent values differ;
    - 2 packed: ``u8 bits`` (1, 2 or 4), then ``ceil(length × bits / 8)``
      bytes; code ``i`` occupies bits ``[i·bits, (i+1)·bits)`` counted from the
      least significant bit of byte 0; unused high bits are zero.
    """
    payload = _as_bytes(data, "decode_user_defined_codes")
    if isinstance(expected_length, bool) or not isinstance(expected_length, int):
        raise TypeError("decode_user_defined_codes: expected_length must be an integer")
    if expected_length < 0 or expected_length > MAX_UINT32:
        raise ValueError(
            "decode_user_defined_codes: expected_length exceeds Uint32 capacity"
        )
    if expected_codes_type not in _CODE_WIDTHS:
        raise ValueError(
            "decode_user_defined_codes: expected_codes_type must be "
            "Uint8Array, Uint16Array or Uint32Array"
        )
    if len(payload) < 2:
        raise ValueError("decode_user_defined_codes: payload too short")

    encoding = payload[0]
    length, offset = decode_uvarint(payload, 1)
    if length != expected_length:
        raise ValueError(
            "decode_user_defined_codes: declared length "
            f"{length} does not match expected length {expected_length}"
        )
    dtype, maximum_value = _CODE_WIDTHS[expected_codes_type]
    native = _NATIVE_CODE_DTYPES[expected_codes_type]

    if encoding == _CODES_RAW:
        nbytes = length * np.dtype(dtype).itemsize
        if len(payload) - offset != nbytes:
            raise ValueError(
                "decode_user_defined_codes: raw payload has the wrong length or trailing bytes"
            )
        return np.frombuffer(payload, dtype=dtype, count=length, offset=offset).astype(native)

    if encoding == _CODES_PACKED:
        if offset >= len(payload):
            raise ValueError("decode_user_defined_codes: packed payload is truncated")
        bits = payload[offset]
        if bits not in (1, 2, 4):
            raise ValueError("decode_user_defined_codes: packed codes must use 1, 2 or 4 bits")
        body = np.frombuffer(payload, dtype=np.uint8, offset=offset + 1)
        if body.size != (length * bits + 7) // 8:
            raise ValueError(
                "decode_user_defined_codes: packed payload has the wrong length or trailing bytes"
            )
        flags = np.unpackbits(body, bitorder="little")
        if np.any(flags[length * bits :]):
            raise ValueError("decode_user_defined_codes: packed codes set bits past their length")
        groups = flags[: length * bits].reshape(length, bits).astype(np.uint32)
        weights = (np.uint32(1) << np.arange(bits, dtype=np.uint32)).astype(np.uint32)
        return (groups @ weights).astype(native)

    if encoding != _CODES_RUNS:
        raise ValueError(f"decode_user_defined_codes: unknown encoding type {encoding}")
    run_count, offset = decode_uvarint(payload, offset)
    if length == 0:
        if run_count != 0:
            raise ValueError(
                "decode_user_defined_codes: empty runs payload must contain zero runs"
            )
    elif run_count == 0 or run_count > length:
        raise ValueError(
            "decode_user_defined_codes: run count cannot fill the declared length"
        )
    pairs = _decode_uvarint_array(
        np.frombuffer(payload, dtype=np.uint8, offset=offset),
        2 * run_count,
        label="decode_user_defined_codes runs",
    )
    values = pairs[0::2]
    runs = pairs[1::2]
    if run_count and int(values.max()) > maximum_value:
        raise ValueError("decode_user_defined_codes: run value exceeds its declared code type")
    if run_count and int(runs.min()) == 0:
        raise ValueError("decode_user_defined_codes: run length must be positive")
    if run_count > 1 and np.any(values[1:] == values[:-1]):
        raise ValueError("decode_user_defined_codes: adjacent runs must have distinct values")
    if int(runs.sum()) != length:
        raise ValueError("decode_user_defined_codes: runs do not fill the declared length")
    return np.repeat(values, runs.astype(np.int64)).astype(native)


def decode_typed_json(data: bytes | bytearray | memoryview) -> Any:
    """
    Decode typed JSON (``analysis/windows``, ``analysis/docked``).

    Layout: ``u32 LE jsonByteLength``, UTF-8 JSON in which the k-th typed array
    is ``{"$typed": k}``, then per array ``u8 dtype`` (1 u8, 2 u16, 3 u32,
    4 i32, 5 f32, 6 f64), ``uvarint count`` and its little-endian elements.
    Typed arrays decode to numpy arrays; every array must be named exactly
    once.
    """
    payload = _as_bytes(data, "decode_typed_json")
    if len(payload) < 4:
        raise ValueError("decode_typed_json: payload is shorter than its header")
    json_length = int.from_bytes(payload[:4], "little")
    if 4 + json_length > len(payload):
        raise ValueError("decode_typed_json: truncated JSON")
    try:
        value = json.loads(payload[4 : 4 + json_length].decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"decode_typed_json: invalid JSON: {error}") from error

    arrays: list[np.ndarray] = []
    offset = 4 + json_length
    while offset < len(payload):
        code = payload[offset]
        offset += 1
        dtype = _TYPED_JSON_DTYPES.get(code)
        if dtype is None:
            raise ValueError(f"decode_typed_json: dtype {code} is not one this version reads")
        count, offset = decode_uvarint(payload, offset)
        nbytes = count * np.dtype(dtype).itemsize
        if offset + nbytes > len(payload):
            raise ValueError("decode_typed_json: truncated array")
        arrays.append(np.frombuffer(payload, dtype=dtype, count=count, offset=offset).copy())
        offset += nbytes

    used: set[int] = set()

    def restore(node: Any, path: str) -> Any:
        if isinstance(node, float) and not np.isfinite(node):
            raise ValueError(f"decode_typed_json: {path} must be finite")
        if isinstance(node, list):
            return [restore(item, f"{path}[{index}]") for index, item in enumerate(node)]
        if isinstance(node, dict):
            if _TYPED_KEY in node:
                k = node[_TYPED_KEY]
                if (
                    len(node) != 1
                    or isinstance(k, bool)
                    or not isinstance(k, int)
                    or not 0 <= k < len(arrays)
                    or k in used
                ):
                    raise ValueError(
                        f"decode_typed_json: {path} names a typed array the payload does not hold"
                    )
                used.add(k)
                return arrays[k]
            return {key: restore(child, f"{path}.{key}") for key, child in node.items()}
        return node

    restored = restore(value, "value")
    if len(used) != len(arrays):
        raise ValueError("decode_typed_json: payload holds arrays its value does not name")
    return restored
