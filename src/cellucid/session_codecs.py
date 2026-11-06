"""
The codecs of a ``.cellucid-session`` file's array chunks (session format 4).

These mirror the web app's exactly:

- ``cellucid/assets/js/app/session/codecs/varint.js`` (unsigned LEB128)
- ``cellucid/assets/js/app/session/codecs/cell-set.js`` (a strictly increasing
  set of cells, stored in whichever of its encodings is smallest)

An ``array`` chunk with encoding ``cell-set`` is a uvarint bound (one past the
largest cell) followed by the cell-set bytes (:mod:`cellucid.session`). Every
decoder is vectorized with numpy: a set of half a million cells is decoded
without a Python loop per cell.
"""

from __future__ import annotations

import numpy as np

from ._plural import plural

MAX_SAFE_INTEGER = (1 << 53) - 1

CELL_SET_SORTED_DELTA = 0
CELL_SET_BITSET = 1
CELL_SET_RUNS = 2


def decode_uvarint(data: bytes | bytearray | memoryview, offset: int = 0) -> tuple[int, int]:
    """Decode one canonical unsigned LEB128 varint at ``offset``.

    Returns
    -------
    value, next_offset
    """
    if not isinstance(data, bytes | bytearray | memoryview):
        raise TypeError("decode_uvarint: data must be bytes-like")
    if isinstance(offset, bool) or not isinstance(offset, int):
        raise TypeError("decode_uvarint: offset must be an integer")
    if offset < 0 or offset >= len(data):
        raise ValueError("decode_uvarint: truncated varint")

    value = 0
    shift = 0
    idx = offset
    length = len(data)
    encoded_bytes = 0

    while True:
        if idx >= length:
            raise ValueError("decode_uvarint: truncated varint")
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


def _bitset_members(bits: np.ndarray, count: int, *, label: str) -> np.ndarray:
    flags = np.unpackbits(bits, bitorder="little")
    if np.any(flags[count:]):
        raise ValueError(f"{label}: bits are set past its length")
    return np.flatnonzero(flags[:count])


def decode_cell_set(data: bytes | bytearray | memoryview, *, cell_count: int) -> np.ndarray:
    """Decode one cell set into a strictly increasing ``uint32`` array.

    Byte 0 is the encoding: 0 sorted-delta (a uvarint count, then each cell's
    gap from the previous one), 1 bitset (``ceil(cell_count / 8)`` bytes, cell
    ``i`` is bit ``i & 7`` of byte ``i >> 3``), 2 runs (a uvarint run count,
    then each run's gap and length). Every cell is below ``cell_count``. Any
    other encoding is refused.
    """
    if not isinstance(data, bytes | bytearray | memoryview):
        raise TypeError("decode_cell_set: data must be bytes-like")
    payload = bytes(data)
    if isinstance(cell_count, bool) or not isinstance(cell_count, int) or cell_count < 0:
        raise TypeError("decode_cell_set: cell_count must be a non-negative integer")
    if not payload:
        raise ValueError("decode_cell_set: empty payload")
    raw = np.frombuffer(payload, dtype=np.uint8)
    encoding = payload[0]

    if encoding == CELL_SET_BITSET:
        if len(payload) != 1 + (cell_count + 7) // 8:
            raise ValueError("decode_cell_set: bitset has the wrong length for its cell count")
        return _bitset_members(raw[1:], cell_count, label="decode_cell_set bitset").astype(np.uint32)

    if encoding == CELL_SET_SORTED_DELTA:
        count, offset = decode_uvarint(payload, 1)
        if count > cell_count:
            raise ValueError(
                f"decode_cell_set: {count} {plural(count, 'cell exceeds', 'cells exceed')} its {cell_count}"
            )
        deltas = _decode_uvarint_array(raw[offset:], count, label="decode_cell_set")
        if count > 1 and np.any(deltas[1:] == 0):
            raise ValueError("decode_cell_set: indices must be strictly increasing")
        indices = np.cumsum(deltas, dtype=np.uint64)
        if count and int(indices[-1]) >= cell_count:
            raise ValueError("decode_cell_set: index exceeds its cell count")
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
            raise ValueError("decode_cell_set: run exceeds its cell count")
        total = int(lengths.sum())
        offsets = np.arange(total, dtype=np.int64) - np.repeat(np.cumsum(lengths) - lengths, lengths)
        return (np.repeat(starts, lengths) + offsets).astype(np.uint32)

    raise ValueError(f"decode_cell_set: encoding {encoding} is not one this version reads")
