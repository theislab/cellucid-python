"""
`.cellucid-session` bundle reader (Python).

This implements the same framing format as the web app:

1) MAGIC bytes: b"CELLUCID_SESSION\\n"
2) manifestByteLength (u32 LE)
3) manifest JSON bytes (UTF-8)
4) repeated chunks: [chunkByteLength (u32 LE), chunkBytes...]

The manifest contains exactly ``createdAt``, ``datasetFingerprint``, and ordered
``chunks`` metadata. Chunk payloads are either JSON or binary and may be
gzip-compressed.
"""

from __future__ import annotations

import json
import re
import shutil
import struct
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

SESSION_BUNDLE_MAGIC = b"CELLUCID_SESSION\n"
U32_BYTES = 4
MAX_MANIFEST_BYTES = 16 * 1024 * 1024
DEFAULT_MAX_UNCOMPRESSED_CHUNK_BYTES = 512 * 1024 * 1024
MAX_STORED_CHUNK_BYTES = 512 * 1024 * 1024

_MANIFEST_KEYS = {
    "createdAt",
    "datasetFingerprint",
    "chunks",
}
_FINGERPRINT_KEYS = {"sourceType", "datasetId", "cellCount", "varCount"}
_CHUNK_REQUIRED_KEYS = {
    "id",
    "contributorId",
    "priority",
    "kind",
    "codec",
    "label",
    "datasetDependent",
    "storedBytes",
    "uncompressedBytes",
}
_UTC_TIMESTAMP_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z$")


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
    return int(struct.unpack("<I", _read_exact(fp, U32_BYTES))[0])


def _require_exact_nonnegative_int(value: Any, *, label: str, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"Invalid session manifest ({label} must be an integer)")
    if value < 0:
        raise ValueError(f"Invalid session manifest ({label} must be non-negative)")
    if value > maximum:
        raise ValueError(f"Invalid session manifest ({label} exceeds the {maximum}-byte limit)")
    return int(value)


def _validate_fingerprint(value: Any) -> None:
    if not isinstance(value, dict):
        raise ValueError("Invalid session manifest (datasetFingerprint must be an object)")
    missing = _FINGERPRINT_KEYS - set(value)
    unknown = set(value) - _FINGERPRINT_KEYS
    if missing or unknown:
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(sorted(missing)))
        if unknown:
            details.append("unknown " + ", ".join(sorted(unknown)))
        raise ValueError(
            "Invalid session manifest (datasetFingerprint fields: " + "; ".join(details) + ")"
        )
    for key in ("sourceType", "datasetId"):
        item = value[key]
        if item is not None and (not isinstance(item, str) or not item):
            raise ValueError(
                f"Invalid session manifest (datasetFingerprint.{key} must be "
                "a non-empty string or null)"
            )
    for key in ("cellCount", "varCount"):
        _require_exact_nonnegative_int(
            value[key],
            label=f"datasetFingerprint.{key}",
            maximum=(1 << 32) - 1,
        )


def _validate_manifest(manifest: Any) -> list[dict[str, Any]]:
    if not isinstance(manifest, dict):
        raise ValueError("Invalid session manifest (expected an object)")
    if set(manifest) != _MANIFEST_KEYS:
        root_missing = sorted(_MANIFEST_KEYS - set(manifest))
        root_unknown = sorted(set(manifest) - _MANIFEST_KEYS)
        root_details: list[str] = []
        if root_missing:
            root_details.append("missing " + ", ".join(root_missing))
        if root_unknown:
            root_details.append("unknown " + ", ".join(root_unknown))
        raise ValueError("Invalid session manifest fields (" + "; ".join(root_details) + ")")
    created_at = manifest["createdAt"]
    if not isinstance(created_at, str) or _UTC_TIMESTAMP_RE.fullmatch(created_at) is None:
        raise ValueError(
            "Invalid session manifest (createdAt must be an exact UTC millisecond timestamp)"
        )
    _validate_fingerprint(manifest["datasetFingerprint"])

    chunks = manifest["chunks"]
    if not isinstance(chunks, list):
        raise ValueError("Invalid session manifest (chunks must be an array)")

    validated: list[dict[str, Any]] = []
    chunk_ids: set[str] = set()
    lazy_seen = False
    for position, raw in enumerate(chunks):
        if not isinstance(raw, dict):
            raise ValueError(f"Invalid session manifest (chunk #{position} must be an object)")
        chunk_missing = _CHUNK_REQUIRED_KEYS - set(raw)
        chunk_unknown = set(raw) - _CHUNK_REQUIRED_KEYS
        if chunk_missing or chunk_unknown:
            chunk_details: list[str] = []
            if chunk_missing:
                chunk_details.append("missing " + ", ".join(sorted(chunk_missing)))
            if chunk_unknown:
                chunk_details.append("unknown " + ", ".join(sorted(chunk_unknown)))
            raise ValueError(
                f"Invalid session manifest (chunk #{position} fields: "
                + "; ".join(chunk_details)
                + ")"
            )

        chunk_id = raw["id"]
        if (
            not isinstance(chunk_id, str)
            or not chunk_id
            or chunk_id != chunk_id.strip()
        ):
            raise ValueError(f"Invalid session manifest (chunk #{position} has an invalid id)")
        if chunk_id in chunk_ids:
            raise ValueError(f"Invalid session manifest (duplicate chunk id {chunk_id!r})")
        chunk_ids.add(chunk_id)

        contributor_id = raw["contributorId"]
        if (
            not isinstance(contributor_id, str)
            or not contributor_id
            or contributor_id != contributor_id.strip()
        ):
            raise ValueError(f"Invalid session manifest (chunk {chunk_id!r} contributorId)")
        priority = raw["priority"]
        if priority not in {"eager", "lazy"}:
            raise ValueError(f"Invalid session manifest (chunk {chunk_id!r} priority)")
        if priority == "lazy":
            lazy_seen = True
        elif lazy_seen:
            raise ValueError("Invalid session manifest (all eager chunks must precede lazy chunks)")
        if raw["kind"] not in {"json", "binary"}:
            raise ValueError(f"Invalid session manifest (chunk {chunk_id!r} kind)")
        if raw["codec"] not in {"none", "gzip"}:
            raise ValueError(f"Invalid session manifest (chunk {chunk_id!r} codec)")
        if not isinstance(raw["label"], str) or not raw["label"]:
            raise ValueError(f"Invalid session manifest (chunk {chunk_id!r} label)")
        if type(raw["datasetDependent"]) is not bool:
            raise ValueError(f"Invalid session manifest (chunk {chunk_id!r} datasetDependent)")
        stored_bytes = _require_exact_nonnegative_int(
            raw["storedBytes"],
            label=f"chunk {chunk_id!r} storedBytes",
            maximum=MAX_STORED_CHUNK_BYTES,
        )
        uncompressed_bytes = _require_exact_nonnegative_int(
            raw["uncompressedBytes"],
            label=f"chunk {chunk_id!r} uncompressedBytes",
            maximum=DEFAULT_MAX_UNCOMPRESSED_CHUNK_BYTES,
        )
        if raw["codec"] == "none" and stored_bytes != uncompressed_bytes:
            raise ValueError(
                f"Invalid session manifest (chunk {chunk_id!r} storedBytes and "
                "uncompressedBytes must match for codec 'none')"
            )

        validated.append(raw)
    return validated


def _decompress_gzip_exact(stored: bytes, expected_bytes: int) -> bytes:
    decompressor = zlib.decompressobj(wbits=16 + zlib.MAX_WBITS)
    output = bytearray()
    remaining_input = stored
    while True:
        remaining_capacity = expected_bytes - len(output)
        decoded = decompressor.decompress(remaining_input, remaining_capacity + 1)
        output.extend(decoded)
        if len(output) > expected_bytes:
            raise ValueError("Session chunk uncompressedBytes is smaller than the gzip payload")
        if decompressor.unconsumed_tail:
            remaining_input = decompressor.unconsumed_tail
            if remaining_capacity == 0:
                continue
            continue
        break

    remaining_capacity = expected_bytes - len(output)
    output.extend(decompressor.flush(remaining_capacity + 1))
    if len(output) > expected_bytes:
        raise ValueError("Session chunk uncompressedBytes is smaller than the gzip payload")
    if not decompressor.eof:
        raise ValueError("Invalid gzip session chunk (truncated stream)")
    if decompressor.unused_data:
        raise ValueError("Invalid gzip session chunk (trailing or concatenated stream)")
    if len(output) != expected_bytes:
        raise ValueError(
            "Session chunk uncompressedBytes does not match the decoded payload "
            f"({expected_bytes} declared, {len(output)} decoded)"
        )
    return bytes(output)


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

    def apply_to_anndata(
        self,
        adata: Any,
        *,
        expected_dataset_id: str,
        inplace: bool = False,
        add_highlights: bool = True,
        highlights_prefix: str = "cellucid_highlight__",
        add_user_defined_fields: bool = True,
        user_defined_prefix: str = "",
        include_deleted_user_defined_fields: bool = False,
        store_uns: bool = True,
        return_summary: bool = False,
    ) -> Any:
        """Apply this exact bundle atomically to its matching AnnData dataset."""
        from .anndata_session import apply_cellucid_session_to_anndata

        return apply_cellucid_session_to_anndata(
            self,
            adata,
            expected_dataset_id=expected_dataset_id,
            inplace=inplace,
            add_highlights=add_highlights,
            highlights_prefix=highlights_prefix,
            add_user_defined_fields=add_user_defined_fields,
            user_defined_prefix=user_defined_prefix,
            include_deleted_user_defined_fields=include_deleted_user_defined_fields,
            store_uns=store_uns,
            return_summary=return_summary,
        )

    @property
    def manifest(self) -> dict[str, Any]:
        """Parsed manifest JSON (loaded lazily)."""
        self._ensure_indexed()
        assert self._manifest is not None
        return self._manifest

    @property
    def dataset_fingerprint(self) -> dict[str, Any]:
        fp = self.manifest.get("datasetFingerprint")
        assert isinstance(fp, dict)
        return fp

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
        chunks = self.manifest["chunks"]
        for entry in chunks:
            chunk_id = entry["id"]
            ref = self._chunks_by_id[chunk_id]
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
            except (UnicodeDecodeError, json.JSONDecodeError) as e:
                raise ValueError(f"Invalid session manifest JSON: {e}") from e
            chunks_meta = _validate_manifest(manifest)

            chunks_by_id: dict[str, SessionChunkRef] = {}
            file_size = self.path.stat().st_size
            for meta in chunks_meta:
                chunk_id = meta["id"]
                stored_len = _read_u32_le(f)
                if stored_len != meta["storedBytes"]:
                    raise ValueError(
                        f"Invalid session chunk {chunk_id!r}: storedBytes declares "
                        f"{meta['storedBytes']}, frame contains {stored_len}"
                    )
                if stored_len > MAX_STORED_CHUNK_BYTES:
                    raise ValueError(
                        f"Invalid session chunk {chunk_id!r}: storedBytes exceeds limit"
                    )
                offset = f.tell()
                if offset + stored_len > file_size:
                    raise ValueError(f"Invalid session chunk {chunk_id!r}: frame exceeds file size")
                f.seek(stored_len, 1)
                chunks_by_id[chunk_id] = SessionChunkRef(
                    id=chunk_id,
                    meta=meta,
                    offset=offset,
                    stored_bytes=stored_len,
                )
            if f.tell() != file_size:
                raise ValueError("Invalid session bundle (trailing bytes after final chunk)")

        self._manifest = manifest
        self._chunks_by_id = chunks_by_id

    def _decode_payload(self, meta: dict[str, Any], stored: bytes) -> Any:
        codec = meta["codec"]
        kind = meta["kind"]
        expected_bytes = meta["uncompressedBytes"]

        if codec == "gzip":
            if expected_bytes > DEFAULT_MAX_UNCOMPRESSED_CHUNK_BYTES:
                raise ValueError("Session chunk uncompressedBytes exceeds the limit")
            uncompressed = _decompress_gzip_exact(stored, expected_bytes)
        else:
            if len(stored) != expected_bytes:
                raise ValueError(
                    "Session chunk uncompressedBytes does not match the stored payload"
                )
            uncompressed = stored

        if kind == "binary":
            return uncompressed

        # kind == "json"
        try:
            return json.loads(uncompressed.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as e:
            raise ValueError(f"Invalid JSON chunk: {e}") from e
