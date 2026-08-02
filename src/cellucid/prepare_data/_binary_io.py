"""
Writing one payload file.

Every binary artifact the export publishes is written through here, so the
bytes carry the format's little-endian order whatever the exporting host is,
the gzip member is deterministic, the bytes land through a temporary file that
is fsynced and atomically renamed, and no write can silently land on a staged
path that already exists.
"""

import os
import tempfile
from pathlib import Path
from typing import BinaryIO, cast

import numpy as np

from .._byte_order import little_endian_payload_view
from .._compression import open_deterministic_gzip_writer


def _output_path_is_writable(
    path: Path,
    description: str,
) -> bool:
    """Admit one new staged artifact path."""
    if path.exists():
        raise FileExistsError(f"Cannot write {description}: staged path already exists: {path}")
    return True


def _gzip_suffix(compression: int | None) -> str:
    """The exact filename suffix a payload written at this setting carries.

    The compression setting decides a payload's name in two independent
    places -- where the manifest declares the path a reader will fetch, and
    where :func:`_write_binary` chooses the path it writes -- and the two names
    have to be one name. Both sides spell the predicate exactly once, here.
    """
    return ".gz" if compression is not None else ""


def _final_binary_path(path: Path, compression: int | None) -> Path:
    """The exact path :func:`_write_binary` publishes one payload at."""
    return Path(f"{path}{_gzip_suffix(compression)}")


def _write_binary(
    path: Path,
    data: np.ndarray,
    compression: int | None = None,
) -> Path:
    """
    Write binary data, optionally with gzip compression.

    The payload is published row-major and little-endian whatever the host's
    own byte order is, because the format's dtype names carry no byte-order
    component and the reader therefore cannot detect a mismatch.

    Parameters
    ----------
    path : Path
        Output path. If compression is enabled, '.gz' will be appended.
    data : np.ndarray
        Data to write.
    compression : int or None
        Gzip compression level (1-9). None means no compression.

    Returns
    -------
    Path
        Actual path written (may have .gz suffix).
    """
    if compression is not None and (type(compression) is not int or not 1 <= compression <= 9):
        raise ValueError("compression must be None or one integer from 1 to 9.")

    payload = little_endian_payload_view(data)
    final_path = _final_binary_path(path, compression)
    if compression is not None:
        _atomic_write_gzip(
            final_path,
            payload,
            compresslevel=compression,
        )
    else:
        final_path.write_bytes(payload)
    return final_path


def _atomic_write_gzip(
    path: Path,
    payload: memoryview,
    *,
    compresslevel: int,
) -> None:
    """Write one deterministic gzip payload and publish it atomically."""
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_file:
            temporary_path = Path(temporary_file.name)
            with open_deterministic_gzip_writer(
                cast(BinaryIO, temporary_file),
                compresslevel=compresslevel,
            ) as compressed:
                compressed.write(payload)
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
        os.chmod(temporary_path, 0o644)
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
