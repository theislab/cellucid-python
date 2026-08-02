"""The one byte order every binary payload is published in.

The format's dtype names -- ``float32``, ``uint16``, ``float64`` -- carry no
byte-order component, and the web app builds its typed arrays straight from the
bytes it receives, which is host order by definition in JavaScript. There is
therefore no field a reader could check and no error it could raise: a payload
written in the wrong order decodes into plausible coordinates and plausible
expression values that are simply not the ones the producer had.

So the order cannot be inherited from whichever machine happens to run the
export. numpy writes native order unless it is told otherwise -- little-endian
on x86 and ARM, big-endian on s390x -- which makes the rule invisible until it
is wrong. Every point where this package turns an array into payload bytes
converts through here first, exactly as ``cellucid-r`` names ``endian =
"little"`` at every ``writeBin``.

The conversion is skipped whenever the array already holds the published order,
so on a little-endian host nothing is copied and nothing moves.
"""

from typing import Any, Literal

import numpy as np
from numpy.typing import DTypeLike

PAYLOAD_BYTE_ORDER: Literal["<"] = "<"


def little_endian_payload_dtype(dtype: DTypeLike) -> np.dtype[Any]:
    """Return the published storage dtype for one payload scalar type."""
    return np.dtype(dtype).newbyteorder(PAYLOAD_BYTE_ORDER)


def little_endian_payload(data: np.ndarray) -> np.ndarray:
    """Return one C-contiguous little-endian array holding the payload values.

    Both properties are what the reader assumes and neither is checkable from
    the bytes: a row-major buffer read in the published order.
    """
    array = np.asarray(data)
    storage_dtype = little_endian_payload_dtype(array.dtype)
    # ``dtype.str`` always spells the byte order out -- '<', '>', or '|' for the
    # single-byte types that have none -- so this comparison is exact, and an
    # array already in the published order is passed through unconverted.
    if array.dtype.str == storage_dtype.str:
        return np.ascontiguousarray(array)
    return np.ascontiguousarray(array, dtype=storage_dtype)


def little_endian_payload_view(data: np.ndarray) -> memoryview:
    """Expose one payload's published bytes without copying an already-exact array."""
    return little_endian_payload(data).data.cast("B")


def little_endian_payload_bytes(data: np.ndarray) -> bytes:
    """Return one payload's exact published bytes."""
    return little_endian_payload_view(data).tobytes()
