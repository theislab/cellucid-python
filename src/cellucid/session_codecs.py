"""
Internal codecs used by `.cellucid-session` bundles.

These mirror the lightweight JS codecs in:
- `cellucid/assets/js/app/session/codecs/varint.js`
- `cellucid/assets/js/app/session/codecs/delta-varint.js`
- `cellucid/assets/js/app/session/contributors/user-defined-codes.js`
"""

from __future__ import annotations

import numpy as np

MAX_SAFE_INTEGER = (1 << 53) - 1
MAX_UINT32 = (1 << 32) - 1
MAX_UINT16 = (1 << 16) - 1


def decode_uvarint(data: bytes | bytearray | memoryview, offset: int = 0) -> tuple[int, int]:
    """
    Decode an unsigned LEB128-style varint.

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
    if max_count is not None and (
        isinstance(max_count, bool) or not isinstance(max_count, int) or max_count < 0
    ):
        raise TypeError("decode_delta_uvarint: max_count must be a non-negative integer")
    if max_index is not None and (
        isinstance(max_index, bool) or not isinstance(max_index, int) or max_index < -1
    ):
        raise TypeError("decode_delta_uvarint: max_index must be an integer of at least -1")

    count, offset = decode_uvarint(data, 0)
    if count > MAX_UINT32:
        raise ValueError("decode_delta_uvarint: count exceeds Uint32 capacity")
    if max_count is not None and count > max_count:
        raise ValueError(f"decode_delta_uvarint: count {count} exceeds max_count {max_count}")

    out = np.empty(count, dtype=np.uint32)
    acc = 0
    for i in range(count):
        delta, offset = decode_uvarint(data, offset)
        if i > 0 and delta == 0:
            raise ValueError("decode_delta_uvarint: indices must be strictly increasing")
        acc += delta
        if acc > MAX_UINT32:
            raise ValueError("decode_delta_uvarint: index exceeds Uint32 capacity")
        if max_index is not None and acc > max_index:
            raise ValueError(f"decode_delta_uvarint: index {acc} exceeds max_index {max_index}")
        out[i] = acc
    if offset != len(data):
        raise ValueError("decode_delta_uvarint: trailing bytes")
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
    if not isinstance(data, bytes | bytearray | memoryview):
        raise TypeError("decode_user_defined_codes: data must be bytes-like")
    if len(data) < 2:
        raise ValueError("decode_user_defined_codes: payload too short")

    enc = data[0]
    length, offset = decode_uvarint(data, 1)

    if enc == 0:
        end = offset + length
        if end > len(data):
            raise ValueError("decode_user_defined_codes: truncated raw u8 payload")
        if end != len(data):
            raise ValueError("decode_user_defined_codes: trailing bytes")
        return np.frombuffer(bytes(data[offset:end]), dtype=np.uint8)

    if enc == 1:
        nbytes = length * 2
        end = offset + nbytes
        if end > len(data):
            raise ValueError("decode_user_defined_codes: truncated raw u16 payload")
        if end != len(data):
            raise ValueError("decode_user_defined_codes: trailing bytes")
        return np.frombuffer(bytes(data[offset:end]), dtype="<u2")

    if enc == 2:
        pair_count, offset = decode_uvarint(data, offset)
        if length == 0:
            if pair_count != 0:
                raise ValueError(
                    "decode_user_defined_codes: empty RLE payload must contain zero pairs"
                )
        elif pair_count == 0 or pair_count > length:
            raise ValueError(
                "decode_user_defined_codes: RLE pair count cannot fill the declared length"
            )
        pairs: list[tuple[int, int]] = []
        needs_u16 = False
        decoded_length = 0
        for _ in range(pair_count):
            value, offset = decode_uvarint(data, offset)
            run, offset = decode_uvarint(data, offset)
            if value > MAX_UINT16:
                raise ValueError(
                    "decode_user_defined_codes: RLE value exceeds Uint16 capacity"
                )
            if run == 0:
                raise ValueError(
                    "decode_user_defined_codes: RLE run length must be positive"
                )
            decoded_length += run
            if decoded_length > length:
                raise ValueError(
                    "decode_user_defined_codes: RLE run total exceeds declared length"
                )
            if value >= 256:
                needs_u16 = True
            pairs.append((value, run))
        if offset != len(data):
            raise ValueError("decode_user_defined_codes: trailing bytes")
        if decoded_length != length:
            raise ValueError("decode_user_defined_codes: RLE did not fill expected length")

        out = np.empty(length, dtype=np.uint16 if needs_u16 else np.uint8)
        i = 0
        for value, run in pairs:
            end = i + run
            out[i:end] = value
            i = end

        return out

    raise ValueError(f"decode_user_defined_codes: unknown encoding type {enc}")
