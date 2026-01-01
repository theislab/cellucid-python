"""
`.cellucid-session` bundle reader (Python).

This implements the same framing format as the web app:

1) MAGIC bytes: b"CELLUCID_SESSION\\n"
2) manifestByteLength (u32 LE)
3) manifest JSON bytes (UTF-8)
4) repeated chunks: [chunkByteLength (u32 LE), chunkBytes...]

The manifest contains chunk metadata in order; chunk payloads are either JSON
or binary and may be gzip-compressed.
"""

from __future__ import annotations

import gzip
import json
import shutil
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO


SESSION_BUNDLE_MAGIC = b"CELLUCID_SESSION\n"
U32_BYTES = 4
MAX_MANIFEST_BYTES = 16 * 1024 * 1024
DEFAULT_MAX_UNCOMPRESSED_CHUNK_BYTES = 512 * 1024 * 1024


@dataclass(frozen=True)
class SessionChunkRef:
    id: str
    meta: dict[str, Any]
    offset: int
    stored_bytes: int


def _read_exact(fp: BinaryIO, n: int) -> bytes:
    data = fp.read(n)
    if data is None or len(data) != n:
        raise ValueError("Unexpected EOF while reading session bundle")
    return data


def _read_u32_le(fp: BinaryIO) -> int:
    return struct.unpack("<I", _read_exact(fp, U32_BYTES))[0]


class CellucidSessionBundle:
    """
    Handle for a `.cellucid-session` file on disk.

    The reader is streaming-friendly: it indexes chunk offsets once, and only
    reads/decompresses chunk payloads on demand.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path).expanduser().resolve()
        if not self.path.is_file():
            raise FileNotFoundError(str(self.path))

        self._manifest: dict[str, Any] | None = None
        self._chunks_by_id: dict[str, SessionChunkRef] | None = None

    def save(self, dest: str | Path) -> Path:
        """Copy the session bundle to `dest` and return the destination path."""
        dest_path = Path(dest).expanduser().resolve()
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(self.path, dest_path)
        return dest_path

    def apply_to_anndata(self, adata: Any, *, inplace: bool = False, **kwargs):
        """
        Apply this bundle to an AnnData.

        See `cellucid.anndata_session.apply_cellucid_session_to_anndata`.
        """
        from .anndata_session import apply_cellucid_session_to_anndata

        return apply_cellucid_session_to_anndata(self, adata, inplace=inplace, **kwargs)

    @property
    def manifest(self) -> dict[str, Any]:
        """Parsed manifest JSON (loaded lazily)."""
        self._ensure_indexed()
        assert self._manifest is not None
        return self._manifest

    @property
    def dataset_fingerprint(self) -> dict[str, Any] | None:
        fp = self.manifest.get("datasetFingerprint")
        return fp if isinstance(fp, dict) else None

    def list_chunk_ids(self) -> list[str]:
        self._ensure_indexed()
        assert self._chunks_by_id is not None
        return list(self._chunks_by_id.keys())

    def get_chunk_meta(self, chunk_id: str) -> dict[str, Any]:
        self._ensure_indexed()
        assert self._chunks_by_id is not None
        ref = self._chunks_by_id.get(chunk_id)
        if ref is None:
            raise KeyError(chunk_id)
        return ref.meta

    def read_chunk_bytes(self, chunk_id: str) -> bytes:
        """Read the stored bytes for a chunk (codec not applied)."""
        self._ensure_indexed()
        assert self._chunks_by_id is not None
        ref = self._chunks_by_id.get(chunk_id)
        if ref is None:
            raise KeyError(chunk_id)
        with self.path.open("rb") as f:
            f.seek(ref.offset)
            return _read_exact(f, ref.stored_bytes)

    def decode_chunk(self, chunk_id: str) -> Any:
        """Decode a chunk payload (apply codec + parse JSON for json-kind chunks)."""
        meta = self.get_chunk_meta(chunk_id)
        stored = self.read_chunk_bytes(chunk_id)
        return self._decode_payload(meta, stored)

    def iter_chunks(self, *, decode: bool = False):
        """Iterate chunks in manifest order."""
        self._ensure_indexed()
        assert self._chunks_by_id is not None
        chunks = self.manifest.get("chunks") or []
        for entry in chunks:
            chunk_id = entry.get("id")
            if not isinstance(chunk_id, str) or not chunk_id:
                continue
            ref = self._chunks_by_id.get(chunk_id)
            if ref is None:
                continue
            if not decode:
                yield ref
            else:
                yield (ref, self._decode_payload(ref.meta, self.read_chunk_bytes(chunk_id)))

    # ---------------------------------------------------------------------
    # Internals
    # ---------------------------------------------------------------------

    def _ensure_indexed(self) -> None:
        if self._manifest is not None and self._chunks_by_id is not None:
            return

        with self.path.open("rb") as f:
            magic = _read_exact(f, len(SESSION_BUNDLE_MAGIC))
            if magic != SESSION_BUNDLE_MAGIC:
                raise ValueError("Not a .cellucid-session file (invalid MAGIC header)")

            manifest_len = _read_u32_le(f)
            if manifest_len <= 0 or manifest_len > MAX_MANIFEST_BYTES:
                raise ValueError(f"Invalid session manifest length: {manifest_len}")

            manifest_bytes = _read_exact(f, manifest_len)
            try:
                manifest = json.loads(manifest_bytes.decode("utf-8"))
            except Exception as e:
                raise ValueError(f"Invalid session manifest JSON: {e}") from e

            chunks_meta = manifest.get("chunks")
            if not isinstance(chunks_meta, list):
                raise ValueError("Invalid session manifest (missing chunks list)")

            chunks_by_id: dict[str, SessionChunkRef] = {}
            for meta in chunks_meta:
                if not isinstance(meta, dict):
                    raise ValueError("Invalid session manifest (chunk meta must be object)")
                chunk_id = meta.get("id")
                if not isinstance(chunk_id, str) or not chunk_id:
                    raise ValueError("Invalid session manifest (chunk id missing)")

                stored_len = _read_u32_le(f)
                offset = f.tell()
                # Bounds check without trusting meta.
                if stored_len < 0:
                    raise ValueError(f"Invalid chunk length for {chunk_id}: {stored_len}")
                f.seek(stored_len, 1)
                chunks_by_id[chunk_id] = SessionChunkRef(
                    id=chunk_id,
                    meta=meta,
                    offset=offset,
                    stored_bytes=stored_len,
                )

        self._manifest = manifest
        self._chunks_by_id = chunks_by_id

    def _decode_payload(self, meta: dict[str, Any], stored: bytes) -> Any:
        codec = meta.get("codec")
        kind = meta.get("kind")

        if codec not in ("none", "gzip"):
            raise ValueError(f"Unsupported session chunk codec: {codec!r}")
        if kind not in ("json", "binary"):
            raise ValueError(f"Unsupported session chunk kind: {kind!r}")

        if codec == "gzip":
            uncompressed = gzip.decompress(stored)
            max_bytes = meta.get("uncompressedBytes")
            if isinstance(max_bytes, int) and max_bytes > 0:
                if len(uncompressed) > max_bytes:
                    raise ValueError(
                        f"Chunk decompressed larger than expected ({len(uncompressed)} > {max_bytes})"
                    )
            elif len(uncompressed) > DEFAULT_MAX_UNCOMPRESSED_CHUNK_BYTES:
                raise ValueError("Chunk decompressed too large (missing uncompressedBytes guard)")
        else:
            uncompressed = stored

        if kind == "binary":
            return uncompressed

        # kind == "json"
        try:
            return json.loads(uncompressed.decode("utf-8"))
        except Exception as e:
            raise ValueError(f"Invalid JSON chunk: {e}") from e
