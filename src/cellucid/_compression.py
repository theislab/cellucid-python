"""Exact gzip encoding shared by prepared exports and direct AnnData responses."""

from __future__ import annotations

import gzip
import io
from typing import BinaryIO


def open_deterministic_gzip_writer(
    fileobj: BinaryIO,
    *,
    compresslevel: int,
) -> gzip.GzipFile:
    """Return a writer with no filename field and a fixed Unix-epoch timestamp."""
    return gzip.GzipFile(
        filename="",
        mode="wb",
        compresslevel=compresslevel,
        fileobj=fileobj,
        mtime=0,
    )


def deterministic_gzip_compress(
    data: bytes,
    *,
    compresslevel: int,
) -> bytes:
    """Compress bytes with one cross-platform canonical gzip header."""
    buffer = io.BytesIO()
    with open_deterministic_gzip_writer(
        buffer,
        compresslevel=compresslevel,
    ) as compressed:
        compressed.write(data)
    return buffer.getvalue()
