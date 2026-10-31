"""
`.cellucid-session` bundle reader (Python).

This implements the same framing format as the web app:

1) MAGIC bytes: b"CELLUCID_SESSION/3\\n" (the session format version)
2) manifestByteLength (u32 LE)
3) manifest JSON bytes (UTF-8)
4) repeated chunks: [chunkByteLength (u32 LE), chunkBytes...]

The manifest contains exactly ``createdAt``, ``datasetFingerprint``, ``source``
(where the dataset can be opened from, or null), ``requires`` (the obs columns
and genes the session uses) and ordered ``chunks`` metadata, ``cells``
when the session was saved on a subsample (which cells were open; see
:attr:`CellucidSessionBundle.cells`), and ``genes`` when genes were read for a
sample of the open cells (see :attr:`CellucidSessionBundle.genes`). Chunk payloads are
either JSON or binary and may be gzip-compressed.
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

from ._plural import plural
from .data_mapping import (
    DataMappingError,
    validate_data_mapping,
    validate_embedding_ref,
)

#: The session format this release reads, and its file signature.
SESSION_FORMAT_VERSION = 3
_SESSION_SIGNATURE_PREFIX = b"CELLUCID_SESSION/"
SESSION_BUNDLE_MAGIC = _SESSION_SIGNATURE_PREFIX + f"{SESSION_FORMAT_VERSION}\n".encode("ascii")
#: The web app's sentence for a file without a session signature
#: (`NOT_A_SESSION_FILE_MESSAGE` in `bundle/format.js`).
NOT_A_SESSION_FILE_MESSAGE = (
    f"This file is not a Cellucid session file of format {SESSION_FORMAT_VERSION} "
    f"(it does not start with CELLUCID_SESSION/{SESSION_FORMAT_VERSION})."
)


def other_session_format_message(version: int) -> str:
    """The web app's sentence for a session file of another format version.

    One saved by an earlier Cellucid (whose cell check this version cannot
    verify), or by a newer one (`otherSessionFormatMessage` in
    `bundle/format.js`).
    """
    if version < SESSION_FORMAT_VERSION:
        return (
            f"This session was saved by an earlier version of Cellucid (session format {version}): "
            "the check that it marks the same cells cannot be verified, so it was not opened. "
            "Re-create it in this version of Cellucid."
        )
    return (
        f"This session was saved by a newer version of Cellucid (session format {version}), "
        "which this version cannot read. Open it in that version of Cellucid."
    )


def session_signature_refusal(magic: bytes | bytearray) -> str:
    """Why a file whose first bytes are ``magic`` is not read.

    ``magic`` is as many bytes as this format's signature holds (fewer for a
    shorter file): another session format's signature, or none
    (`sessionSignatureRefusal` in `bundle/format.js`).
    """
    signature = bytes(magic)
    match = (
        re.match(rb"(\d+)(?:\n|$)", signature[len(_SESSION_SIGNATURE_PREFIX) :])
        if signature.startswith(_SESSION_SIGNATURE_PREFIX)
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

_MANIFEST_KEYS = {
    "createdAt",
    "datasetFingerprint",
    "source",
    "requires",
    "chunks",
}
# Where the session's dataset can be opened from, so loading the session in the
# viewer can open it first: a sample of the page's catalog, a data link (its
# ``?data=`` value) or a Cellucid server; or the name of the local file or
# folder it was opened from, which the viewer can only ask for. Null for a
# notebook.
_SOURCE_KEYS_BY_KIND = {
    "sample": {"kind", "dataset"},
    "link": {"kind", "locator", "dataset"},
    "remote": {"kind", "url", "dataset"},
    "file": {"kind", "name", "dataset"},
}
_REQUIRES_KEYS = {"obs", "var"}
_CONTROL_CHARACTER_RE = re.compile(r"[\x00-\x1f\x7f]")
# The fingerprint ties a session to one dataset and to one data mapping, not to
# the route the dataset was opened through. Every key is required: the viewer
# refuses a session without its cell-order record, and a reader that accepted
# one would apply selections the viewer would refuse. ``name`` is the dataset's
# display name, recorded for messages only.
_FINGERPRINT_KEYS = {
    "datasetId",
    "name",
    "cellCount",
    "varCount",
    "cellOrder",
    "dataMapping",
}
_CELL_ORDER_KEYS = {"embedding", "dimension", "digest"}
_CELL_ORDER_DIGEST_RE = re.compile(r"\A[0-9a-f]{16}\Z")
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


def _validate_fingerprint(value: Any) -> dict[str, Any]:
    """Validate one dataset fingerprint and return it with a canonical mapping.

    ``dataMapping`` names the embeddings, gene names and (for an AnnData
    source) the graph, expression matrix and latent space the session was
    saved with; ``cellOrder`` names the one embedding whose stored values the
    viewer digested, each cell by its original position (the web app's
    ``data/cell-order-digest.js``). Its shape is checked here; the digest is
    the viewer's check that a session is restored on the cells and
    coordinates it was saved on, and applying a session does not repeat it.
    """
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
    for key in ("datasetId", "name"):
        item = value[key]
        if item is not None and (not isinstance(item, str) or not item or item != item.strip()):
            raise ValueError(
                f"Invalid session manifest (datasetFingerprint.{key} must be "
                "a non-empty trimmed string or null)"
            )
    for key in ("cellCount", "varCount"):
        _require_exact_nonnegative_int(
            value[key],
            label=f"datasetFingerprint.{key}",
            maximum=(1 << 32) - 1,
        )
    try:
        mapping = validate_data_mapping(value["dataMapping"], label="datasetFingerprint.dataMapping")
    except DataMappingError as error:
        raise ValueError(f"Invalid session manifest ({error})") from error
    cell_order = _validate_cell_order(value["cellOrder"])
    canonical = dict(value)
    canonical["cellOrder"] = cell_order
    canonical["dataMapping"] = mapping
    return canonical


def _validate_cell_order(value: Any) -> dict[str, Any]:
    """Validate the cell-order record that ties a session to a cell ordering.

    The viewer records a digest of the coordinates a session was saved against,
    so a dataset republished with the same name and counts but a different row
    order cannot silently re-point every saved selection at different cells.
    The digest is over one embedding's stored values (rounded to float32, each
    cell named by its original position, summed so that how the cells were
    opened does not matter), so the record names that embedding -- an AnnData
    ``obsm`` key with the 0-based columns read -- and its dimension.
    """
    if not isinstance(value, dict):
        raise ValueError(
            "Invalid session manifest (datasetFingerprint.cellOrder must be an object)"
        )
    missing = _CELL_ORDER_KEYS - set(value)
    unknown = set(value) - _CELL_ORDER_KEYS
    if missing or unknown:
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(sorted(missing)))
        if unknown:
            details.append("unknown " + ", ".join(sorted(unknown)))
        raise ValueError(
            "Invalid session manifest (datasetFingerprint.cellOrder fields: "
            + "; ".join(details)
            + ")"
        )
    dimension = value["dimension"]
    if type(dimension) is not int or dimension not in (1, 2, 3):
        raise ValueError(
            "Invalid session manifest (datasetFingerprint.cellOrder.dimension "
            "must be exactly 1, 2, or 3)"
        )
    digest = value["digest"]
    if not isinstance(digest, str) or not _CELL_ORDER_DIGEST_RE.match(digest):
        raise ValueError(
            "Invalid session manifest (datasetFingerprint.cellOrder.digest must be "
            "16 lowercase hexadecimal characters)"
        )
    try:
        embedding = validate_embedding_ref(
            value["embedding"],
            dimension=dimension,
            label="datasetFingerprint.cellOrder.embedding",
        )
    except DataMappingError as error:
        raise ValueError(f"Invalid session manifest ({error})") from error
    return {"embedding": embedding, "dimension": dimension, "digest": digest}


def _validate_text(value: Any, label: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value != value.strip()
        or _CONTROL_CHARACTER_RE.search(value) is not None
    ):
        raise ValueError(f"Invalid session manifest ({label} must be a non-empty, unpadded string)")
    return value


def _validate_source(value: Any) -> dict[str, str] | None:
    """Validate the recorded dataset source (the viewer's dataset-source-record)."""
    if value is None:
        return None
    kind = value.get("kind") if isinstance(value, dict) else None
    keys = _SOURCE_KEYS_BY_KIND.get(kind) if isinstance(kind, str) else None
    if keys is None:
        raise ValueError(
            'Invalid session manifest (source kind must be "sample", "link", "remote" or "file")'
        )
    if set(value) != keys:
        raise ValueError(
            "Invalid session manifest (source fields must be exactly " + ", ".join(sorted(keys)) + ")"
        )
    for key in sorted(keys - {"kind"}):
        _validate_text(value[key], f"source.{key}")
    if kind == "remote" and not value["url"].startswith(("https://", "http://")):
        raise ValueError("Invalid session manifest (source.url must be an http(s) address)")
    return dict(value)


def _utf16_key(value: str) -> bytes:
    """Sort key matching JavaScript string order (UTF-16 code units)."""
    return value.encode("utf-16-be", "surrogatepass")


def _validate_requires(value: Any) -> dict[str, list[str]]:
    """Validate the manifest's ``requires``: the dataset fields the session uses.

    Each list holds the dataset's own keys (obs columns; gene names under the
    session's data mapping), sorted as JavaScript sorts strings and unique.
    """
    if not isinstance(value, dict) or set(value) != _REQUIRES_KEYS:
        raise ValueError("Invalid session manifest (requires must have exactly obs and var)")
    out: dict[str, list[str]] = {}
    for source in ("obs", "var"):
        keys = value[source]
        if not isinstance(keys, list):
            raise ValueError(f"Invalid session manifest (requires.{source} must be an array)")
        for index, key in enumerate(keys):
            if not isinstance(key, str) or not key or key != key.strip():
                raise ValueError(
                    f"Invalid session manifest (requires.{source}[{index}] must be a "
                    "non-empty trimmed string)"
                )
            if index > 0 and not _utf16_key(keys[index - 1]) < _utf16_key(key):
                raise ValueError(
                    f"Invalid session manifest (requires.{source} must be sorted and "
                    "without duplicates)"
                )
        out[source] = list(keys)
    return out


def _decode_cells(cells: Any, cell_count: Any) -> dict[str, Any]:
    """The ``cells`` record of a session saved on a subsample, checked and decoded.

    ``{"spec": {...}, "kept": <base64 of the cell-set bytes>, "digest": <16 hex>}``
    (docs: *Subsampling specification*, "Cell set encoding"). Decoded:
    ``kept``, the ascending original positions: the session's per-cell state
    of cell ``i`` is obs row ``kept[i]`` (per-cell chunks name cells by
    original position, ascending, whatever order the viewer opened them in).
    """
    import base64
    import binascii

    from .subsample import SUBSAMPLE_METHODS, cell_set_digest, decode_cell_set

    if not isinstance(cells, dict) or set(cells) != {"spec", "kept", "digest"}:
        raise ValueError("Invalid session manifest (cells must have exactly spec, kept and digest)")
    spec = cells["spec"]
    if not isinstance(spec, dict) or spec.get("method") not in SUBSAMPLE_METHODS:
        raise ValueError("Invalid session manifest (cells.spec must name a subsample method)")
    # Exactly the keys of its method, as the web app's requireSubsampleSpec.
    method_keys = {"uniform": set(), "blocks": {"blockSize"}, "stratified": {"field", "allocation", "minPerCategory"},
                   "spatial": {"embedding", "columns", "gridSize"}}[spec["method"]]
    if set(spec) != {"method", "nObs", "k", "seed"} | method_keys:
        raise ValueError(f"Invalid session manifest (cells.spec of {spec['method']} must have exactly method, nObs, k, seed"
                         + "".join(f", {key}" for key in sorted(method_keys)) + ")")
    n_obs, k = spec.get("nObs"), spec.get("k")
    if not (isinstance(n_obs, int) and isinstance(k, int) and 1 <= k < n_obs <= 0xFFFFFFFF):
        raise ValueError("Invalid session manifest (cells.spec must keep 1 <= k < nObs cells)")
    if k != cell_count:
        raise ValueError("Invalid session manifest (cells must keep exactly the fingerprint cell count)")
    if not isinstance(cells["kept"], str) or not isinstance(cells["digest"], str):
        raise ValueError("Invalid session manifest (cells.kept and cells.digest must be text)")
    try:
        raw = base64.b64decode(cells["kept"], validate=True)
    except (binascii.Error, ValueError) as error:
        raise ValueError("Invalid session manifest (cells.kept must be base64)") from error
    kept = decode_cell_set(raw, n_obs)
    if kept is None or kept.size != k:
        raise ValueError("Invalid session manifest (cells.kept must list exactly k cells)")
    if cell_set_digest(n_obs, kept) != cells["digest"]:
        raise ValueError("Invalid session manifest (cells.digest does not match its cells)")
    if spec.get("method") == "blocks" and not (isinstance(spec.get("blockSize"), int) and spec["blockSize"] >= 1):
        raise ValueError("Invalid session manifest (cells.spec must give a blocks subsample its blockSize)")
    return {"spec": dict(spec), "kept": kept, "digest": cells["digest"]}


def _decode_genes(genes: Any, cell_count: Any) -> dict[str, Any]:
    """The ``genes`` record of a session saved while genes were read for a sample of the open cells.

    ``{"cells": <count>, "order": null | {"method": "blocks", "seed", "blockSize"}, "digest": <16 hex>}``:
    genes were read for the first ``cells`` open cells in their order (the
    order of the ``blocks`` subsample of every cell ``order`` names -- every
    cell's, or a share's spread through a file in its own order, its cells
    filtered -- or, with ``order`` null, a subsample's priority order),
    whose original positions have the cell-set digest ``digest``. Checked as
    the web app checks it.
    """
    if not isinstance(genes, dict) or set(genes) != {"cells", "order", "digest"}:
        raise ValueError("Invalid session manifest (genes must have exactly cells, order and digest)")
    cells = genes["cells"]
    if not (isinstance(cells, int) and not isinstance(cells, bool) and isinstance(cell_count, int) and 1 <= cells < cell_count):
        raise ValueError("Invalid session manifest (genes.cells must be a count of open cells below the fingerprint cell count)")
    order = genes["order"]
    if order is not None:
        if not isinstance(order, dict) or set(order) != {"method", "seed", "blockSize"} or order["method"] != "blocks":
            raise ValueError("Invalid session manifest (genes.order must be {method: 'blocks', seed, blockSize})")
        seed, block_size = order["seed"], order["blockSize"]
        if not (isinstance(seed, int) and not isinstance(seed, bool) and 0 <= seed <= 0xFFFFFFFF):
            raise ValueError("Invalid session manifest (genes.order.seed must be a uint32)")
        if not (isinstance(block_size, int) and not isinstance(block_size, bool) and block_size >= 1):
            raise ValueError("Invalid session manifest (genes.order.blockSize must be a positive count)")
    digest = genes["digest"]
    if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{16}", digest) is None:
        raise ValueError("Invalid session manifest (genes.digest must be 16 lowercase hex digits)")
    return {"cells": cells, "order": None if order is None else dict(order), "digest": digest}


def _validate_manifest(manifest: Any) -> list[dict[str, Any]]:
    if not isinstance(manifest, dict):
        raise ValueError("Invalid session manifest (expected an object)")
    expected = _MANIFEST_KEYS | ({"cells"} if "cells" in manifest else set()) | ({"genes"} if "genes" in manifest else set())
    if set(manifest) != expected:
        root_missing = sorted(expected - set(manifest))
        root_unknown = sorted(set(manifest) - expected)
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
    _validate_source(manifest["source"])
    _validate_requires(manifest["requires"])
    if "cells" in manifest:
        _decode_cells(manifest["cells"], manifest["datasetFingerprint"].get("cellCount"))
    if "genes" in manifest:
        _decode_genes(manifest["genes"], manifest["datasetFingerprint"].get("cellCount"))

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
            maximum=MAX_FRAME_BYTES,
        )
        uncompressed_bytes = _require_exact_nonnegative_int(
            raw["uncompressedBytes"],
            label=f"chunk {chunk_id!r} uncompressedBytes",
            maximum=MAX_SAFE_INTEGER,
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

    Parameters
    ----------
    path : str or Path
        The ``.cellucid-session`` file; ``~`` is expanded. A missing file
        raises ``FileNotFoundError``.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path).expanduser().resolve()
        if not self.path.is_file():
            raise FileNotFoundError(str(self.path))

        self._manifest: dict[str, Any] | None = None
        self._chunks_by_id: dict[str, SessionChunkRef] | None = None

    def save(self, dest: str | Path) -> Path:
        """Copy the session bundle to `dest` and return the destination path.

        The copy is written beside ``dest`` and moved into place once
        complete, so an interrupted save never leaves a partial file there;
        saving a bundle onto its own file leaves it as it is.
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
        """Apply this exact bundle atomically to its matching AnnData dataset.

        Same as :func:`~cellucid.apply_cellucid_session_to_anndata` with this
        bundle; see it for the parameters, the return value, and the errors.
        """
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
        """The manifest's ``datasetFingerprint``: the dataset the session was saved on."""
        fp = self.manifest.get("datasetFingerprint")
        assert isinstance(fp, dict)
        return fp

    @property
    def dataset_source(self) -> dict[str, str] | None:
        """Where the session's dataset was opened from, or ``None``.

        ``{"kind": "sample", "dataset"}``, ``{"kind": "link", "locator",
        "dataset"}`` (``locator`` is the ``?data=`` value), ``{"kind":
        "remote", "url", "dataset"}`` or ``{"kind": "file", "name",
        "dataset"}`` (the local file or folder name); ``None`` for a session
        saved in a notebook.
        """
        source = self.manifest["source"]
        return None if source is None else dict(source)

    @property
    def cells(self) -> dict[str, Any] | None:
        """The subsample the session was saved on, or ``None`` for every cell.

        ``{"spec": {...}, "kept": numpy.ndarray, "digest": str}``: the
        subsample's specification (``method``, ``nObs``, ``k``, ``seed`` and the
        method's own keys), the kept cells' original positions (rows of
        ``obs``, ascending, int64) and their cell-set digest. Every per-cell
        chunk of such a session is in the positions of these cells: cell ``i``
        of a chunk is ``obs`` row ``kept[i]``. The same cells are
        ``cellucid.subsample_cells(...)`` of the specification.
        """
        if "cells" not in self.manifest:
            return None
        return _decode_cells(self.manifest["cells"], self.dataset_fingerprint.get("cellCount"))

    @property
    def genes(self) -> dict[str, Any] | None:
        """The open cells genes were read for, or ``None`` when every open cell's were.

        ``{"cells": int, "order": dict | None, "digest": str}``: a large dataset
        stored cell by cell opens every cell on the map and reads a gene for a
        sample of them, the first ``cells`` open cells in their order (a
        subsample's priority order, or, every cell being open, the order of
        the ``blocks`` subsample ``order`` names); ``digest`` is the cell-set
        digest of their original positions. Gene values the session holds
        (none of its per-cell chunks do) were read for those cells only;
        opening the session in the viewer reads genes for them again.
        """
        if "genes" not in self.manifest:
            return None
        return _decode_genes(self.manifest["genes"], self.dataset_fingerprint.get("cellCount"))

    @property
    def required_fields(self) -> dict[str, list[str]]:
        """The dataset fields the session uses: ``{"obs": [...], "var": [...]}``.

        Every obs column and gene (named under the session's data mapping) that
        a view shows or filters on, or that a highlight group is bound to.
        """
        requires = self.manifest["requires"]
        return {"obs": list(requires["obs"]), "var": list(requires["var"])}

    def list_chunk_ids(self) -> list[str]:
        """The ids of the bundle's chunks, in manifest order."""
        self._ensure_indexed()
        assert self._chunks_by_id is not None
        return list(self._chunks_by_id.keys())

    def get_chunk_meta(self, chunk_id: str) -> dict[str, Any]:
        """The manifest entry of one chunk; ``KeyError`` for an unknown id."""
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
            magic = f.read(len(SESSION_BUNDLE_MAGIC))
            if magic != SESSION_BUNDLE_MAGIC:
                raise ValueError(session_signature_refusal(magic))

            manifest_len = _read_u32_le(f)
            if manifest_len <= 0:
                raise ValueError(f"Invalid session manifest length: {manifest_len}")
            if f.tell() + manifest_len > self.path.stat().st_size:
                raise ValueError(
                    f"Invalid session manifest length: {manifest_len} {plural(manifest_len, 'byte')} is more "
                    "than the session file holds"
                )

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
