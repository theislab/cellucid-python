"""
Internal codecs used by `.cellucid-session` bundles.

These mirror the lightweight JS codecs in:
- `cellucid/assets/js/app/session/codecs/varint.js`
- `cellucid/assets/js/app/session/codecs/delta-varint.js`
- `cellucid/assets/js/app/session/contributors/user-defined-codes.js`
"""

from __future__ import annotations

import numpy as np


def decode_uvarint(data: bytes | bytearray | memoryview, offset: int = 0) -> tuple[int, int]:
    """
    Decode an unsigned LEB128-style varint.

    Returns
    -------
    value, next_offset
    """
    if not isinstance(data, (bytes, bytearray, memoryview)):
        raise TypeError("decode_uvarint: data must be bytes-like")

    value = 0
    shift = 0
    idx = offset
    length = len(data)

    while True:
        if idx >= length:
            raise ValueError("decode_uvarint: truncated")
        b = data[idx]
        idx += 1

        value |= (b & 0x7F) << shift
        if (b & 0x80) == 0:
            return value, idx

        shift += 7
        if shift > 63:
            raise ValueError("decode_uvarint: varint too long")


def decode_delta_uvarint(
    data: bytes | bytearray | memoryview,
    *,
    max_count: int | None = None,
    max_index: int | None = None,
) -> np.ndarray:
    """
    Decode delta+uvarint encoded sorted indices (pre-gzip payload).

    Format:
    - count (uvarint)
    - count deltas (uvarint), where:
        - idx0 = delta0
        - idxi = idx(i-1) + deltai
    """
    count, offset = decode_uvarint(data, 0)
    if max_count is not None and count > max_count:
        raise ValueError(f"decode_delta_uvarint: count {count} exceeds max_count {max_count}")

    out = np.empty(count, dtype=np.uint32)
    acc = 0
    for i in range(count):
        delta, offset = decode_uvarint(data, offset)
        acc += delta
        if max_index is not None and acc > max_index:
            raise ValueError(f"decode_delta_uvarint: index {acc} exceeds max_index {max_index}")
        out[i] = acc
    return out


def decode_user_defined_codes(data: bytes | bytearray | memoryview) -> np.ndarray:
    """
    Decode a user-defined categorical codes chunk (pre-gzip payload).

    Encoding (byte layout):
    - 1 byte: encodingType
      - 0 = raw Uint8
      - 1 = raw Uint16 (little-endian)
      - 2 = RLE pairs encoded as uvarints (value, runLength)
    - codesLength (uvarint)
    - payload:
      - raw: `codesLength * bytesPerElement` bytes
      - RLE: pairCount (uvarint) followed by (value uvarint, runLength uvarint) pairs
    """
    if not isinstance(data, (bytes, bytearray, memoryview)):
        raise TypeError("decode_user_defined_codes: data must be bytes-like")
    if len(data) < 2:
        raise ValueError("decode_user_defined_codes: payload too short")

    enc = data[0]
    length, offset = decode_uvarint(data, 1)

    if enc == 0:
        end = offset + length
        if end > len(data):
            raise ValueError("decode_user_defined_codes: truncated raw u8 payload")
        return np.frombuffer(bytes(data[offset:end]), dtype=np.uint8)

    if enc == 1:
        nbytes = length * 2
        end = offset + nbytes
        if end > len(data):
            raise ValueError("decode_user_defined_codes: truncated raw u16 payload")
        return np.frombuffer(bytes(data[offset:end]), dtype="<u2")

    if enc == 2:
        pair_count, offset = decode_uvarint(data, offset)
        pairs: list[tuple[int, int]] = []
        needs_u16 = False
        for _ in range(pair_count):
            value, offset = decode_uvarint(data, offset)
            run, offset = decode_uvarint(data, offset)
            if value >= 256:
                needs_u16 = True
            pairs.append((value, run))

        out = np.empty(length, dtype=np.uint16 if needs_u16 else np.uint8)
        i = 0
        for value, run in pairs:
            if run <= 0:
                continue
            end = min(length, i + run)
            out[i:end] = value
            i = end
            if i >= length:
                break

        if i != length:
            raise ValueError("decode_user_defined_codes: RLE did not fill expected length")

        return out

    raise ValueError(f"decode_user_defined_codes: unknown encoding type {enc}")

