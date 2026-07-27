"""
Cellucid Data Server

A lightweight HTTP server for serving pre-exported cellucid datasets.
Supports both local and remote access patterns:

1. Local mode: Run on your machine, open browser locally
2. SSH tunnel mode: Run on remote server, access via SSH port forwarding
3. Jupyter mode: Run alongside Jupyter, embed in notebooks

Usage:
    from cellucid import serve
    serve("/path/to/data", port=8765)

Or via CLI:
    cellucid serve /path/to/data --port 8765

The server provides:
- Static file serving for dataset files
- CORS headers for cross-origin access (needed for web viewer)
- Health check endpoint for connection validation

For serving AnnData directly (without pre-export), use:
    from cellucid import serve_anndata
    serve_anndata(
        "/path/to/data.h5ad",
        dataset_name="Example",
        dataset_id="example",
    )
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import stat
import threading
import webbrowser
from dataclasses import dataclass
from functools import partial
from http import HTTPStatus
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from urllib.parse import urlparse

from ._server_base import (
    CELLUCID_WEB_URL,
    DEFAULT_HOST,
    DEFAULT_PORT,
    WEB_ASSET_INVENTORY_FILENAME,
    CORSMixin,
    _web_cache_dir,
    print_detail,
    print_server_banner,
    print_step,
    print_success,
    require_server_port,
)
from .prepare_data import (
    _require_nonempty_string,
    _require_portable_filename_component,
)

logger = logging.getLogger("cellucid.server")

_PUBLISHED_STATE_MANIFEST = "state-snapshots.json"
_PUBLISHED_STATE_BUNDLE = "default.cellucid-session"
_PUBLISHED_STATE_MANIFEST_BYTES = 43
_MAX_PUBLISHED_STATE_BUNDLE_BYTES = 32 * 1024


def _same_regular_file(left: os.stat_result, right: os.stat_result) -> bool:
    return (
        stat.S_ISREG(left.st_mode)
        and stat.S_ISREG(right.st_mode)
        and left.st_size == right.st_size
        and left.st_mtime_ns == right.st_mtime_ns
        and left.st_dev == right.st_dev
        and left.st_ino == right.st_ino
    )


def _read_published_state_file(
    dataset_root: Path,
    filename: str,
    *,
    minimum_bytes: int,
    maximum_bytes: int,
) -> bytes:
    """Read one unchanged in-root regular state sidecar without following links."""
    candidate = dataset_root / filename
    try:
        before = candidate.lstat()
    except OSError as error:
        raise ValueError(
            f"Published sample state sidecar is unreadable: {filename}"
        ) from error
    if stat.S_ISLNK(before.st_mode):
        raise ValueError(
            f"Published sample state sidecar must not be a symbolic link: {filename}"
        )
    if not stat.S_ISREG(before.st_mode):
        raise ValueError(
            f"Published sample state sidecar must be a regular file: {filename}"
        )
    if before.st_size < minimum_bytes or before.st_size > maximum_bytes:
        if minimum_bytes == maximum_bytes:
            raise ValueError(
                f"Published sample state {filename} must be exactly "
                f"{minimum_bytes} bytes."
            )
        raise ValueError(
            f"Published sample state {filename} size must be from "
            f"{minimum_bytes} through {maximum_bytes} bytes."
        )

    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(dataset_root)
    except (OSError, ValueError) as error:
        raise ValueError(
            f"Published sample state sidecar must remain inside its dataset root: {filename}"
        ) from error
    if resolved != candidate:
        raise ValueError(
            f"Published sample state sidecar must not traverse a symbolic link: {filename}"
        )

    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0)
    nofollow = getattr(os, "O_NOFOLLOW", 0)
    if nofollow:
        flags |= nofollow
    try:
        descriptor = os.open(candidate, flags)
    except OSError as error:
        raise ValueError(
            f"Published sample state sidecar could not be opened safely: {filename}"
        ) from error

    try:
        opened = os.fstat(descriptor)
        if not _same_regular_file(before, opened):
            raise ValueError(
                f"Published sample state sidecar changed during validation: {filename}"
            )
        remaining = opened.st_size
        chunks: list[bytes] = []
        while remaining:
            chunk = os.read(descriptor, min(64 * 1024, remaining))
            if not chunk:
                raise ValueError(
                    f"Published sample state sidecar ended during validation: {filename}"
                )
            chunks.append(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            raise ValueError(
                f"Published sample state sidecar grew during validation: {filename}"
            )
        after = candidate.lstat()
        if not _same_regular_file(opened, after):
            raise ValueError(
                f"Published sample state sidecar changed during validation: {filename}"
            )
        return b"".join(chunks)
    finally:
        os.close(descriptor)


def _read_optional_published_state(dataset_dir: Path) -> dict[str, str]:
    """Validate and describe the exact optional published sample state pair."""
    dataset_root = dataset_dir.resolve(strict=True)
    manifest_path = dataset_root / _PUBLISHED_STATE_MANIFEST
    bundle_path = dataset_root / _PUBLISHED_STATE_BUNDLE
    manifest_present = manifest_path.is_symlink() or manifest_path.exists()
    bundle_present = bundle_path.is_symlink() or bundle_path.exists()
    if not manifest_present and not bundle_present:
        return {}
    if not manifest_present or not bundle_present:
        raise ValueError(
            "Published sample state requires the exact pair "
            f"{_PUBLISHED_STATE_MANIFEST} and {_PUBLISHED_STATE_BUNDLE}."
        )

    extra_bundles = sorted(
        child.name
        for child in dataset_root.iterdir()
        if child.name.endswith(".cellucid-session")
        and child.name != _PUBLISHED_STATE_BUNDLE
    )
    if extra_bundles:
        raise ValueError(
            "Published sample state has extra session bundles; exactly one "
            f"{_PUBLISHED_STATE_BUNDLE} is permitted: {', '.join(extra_bundles)}"
        )

    manifest_bytes = _read_published_state_file(
        dataset_root,
        _PUBLISHED_STATE_MANIFEST,
        minimum_bytes=_PUBLISHED_STATE_MANIFEST_BYTES,
        maximum_bytes=_PUBLISHED_STATE_MANIFEST_BYTES,
    )
    try:
        manifest = json.loads(manifest_bytes.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(
            f"Published sample state {_PUBLISHED_STATE_MANIFEST} "
            "must contain readable UTF-8 JSON."
        ) from error
    if not isinstance(manifest, dict):
        raise TypeError(
            f"Published sample state {_PUBLISHED_STATE_MANIFEST} "
            "must contain a JSON object."
        )
    if set(manifest) != {"states"}:
        raise ValueError(
            f"Published sample state {_PUBLISHED_STATE_MANIFEST} "
            "must contain the exact keys: states."
        )
    if type(manifest["states"]) is not list or manifest["states"] != [
        _PUBLISHED_STATE_BUNDLE
    ]:
        raise ValueError(
            f"Published sample state {_PUBLISHED_STATE_MANIFEST} states "
            f"must contain exactly {_PUBLISHED_STATE_BUNDLE}."
        )

    bundle_bytes = _read_published_state_file(
        dataset_root,
        _PUBLISHED_STATE_BUNDLE,
        minimum_bytes=1,
        maximum_bytes=_MAX_PUBLISHED_STATE_BUNDLE_BYTES,
    )
    return {
        "state_manifest": _PUBLISHED_STATE_MANIFEST,
        "state_sha256": hashlib.sha256(bundle_bytes).hexdigest(),
    }


def _read_exported_dataset_entry(
    dataset_dir: Path,
    *,
    public_path: str,
) -> dict[str, str]:
    """Validate one complete prepared-dataset root and return its catalog entry."""
    required_files = (
        dataset_dir / "dataset_identity.json",
        dataset_dir / "obs_manifest.json",
    )
    point_files: list[Path] = []
    for dimension in (1, 2, 3):
        candidates = [
            dataset_dir / f"points_{dimension}d.bin",
            dataset_dir / f"points_{dimension}d.bin.gz",
        ]
        existing = [candidate for candidate in candidates if candidate.exists()]
        if len(existing) > 1:
            raise ValueError(
                f"{dataset_dir.name!r} is not a complete exported dataset: "
                f"both compressed and uncompressed {dimension}D points exist."
            )
        point_files.extend(existing)

    if (
        not dataset_dir.is_dir()
        or any(not path.is_file() for path in required_files)
        or not point_files
        or any(path.stat().st_size == 0 for path in point_files)
    ):
        raise ValueError(f"{dataset_dir.name!r} is not a complete exported dataset.")

    try:
        identity = json.loads(required_files[0].read_text(encoding="utf-8"))
        obs_manifest = json.loads(required_files[1].read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(
            f"{dataset_dir.name!r} is not a complete exported dataset with "
            "readable UTF-8 JSON metadata."
        ) from error
    if not isinstance(identity, dict):
        raise TypeError("dataset_identity.json must contain a JSON object.")
    if not isinstance(obs_manifest, dict):
        raise TypeError("obs_manifest.json must contain a JSON object.")
    if type(identity.get("version")) is not int or identity["version"] != 2:
        raise ValueError("dataset_identity.json version must be exactly 2.")

    dataset_id = _require_nonempty_string(
        identity.get("id"),
        label="dataset_id",
    )
    dataset_name = _require_nonempty_string(
        identity.get("name"),
        label="dataset_name",
    )
    for label, value in (("dataset_id", dataset_id), ("dataset_name", dataset_name)):
        if value != value.strip() or any(
            ord(character) < 32 or ord(character) == 127 for character in value
        ):
            raise ValueError(
                f"{label} must be exact text without surrounding whitespace or control characters."
            )

    entry = {
        "id": dataset_id,
        "path": public_path,
        "name": dataset_name,
    }
    entry.update(_read_optional_published_state(dataset_dir))
    return entry


def _list_exported_datasets(data_dir: Path) -> list[dict[str, str]]:
    """Classify one exact prepared dataset or a root containing only datasets."""
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        raise NotADirectoryError(f"Prepared-data path must be a directory: {data_dir}")

    direct_markers = (
        data_dir / "dataset_identity.json",
        data_dir / "obs_manifest.json",
        *data_dir.glob("points_*d.bin"),
        *data_dir.glob("points_*d.bin.gz"),
    )
    if any(path.exists() for path in direct_markers):
        return [_read_exported_dataset_entry(data_dir, public_path="/")]

    subdirectories = sorted(path for path in data_dir.iterdir() if path.is_dir())
    if not subdirectories:
        return []

    candidate_flags = [
        any(
            marker.exists()
            for marker in (
                subdir / "dataset_identity.json",
                subdir / "obs_manifest.json",
                *subdir.glob("points_*d.bin"),
                *subdir.glob("points_*d.bin.gz"),
            )
        )
        for subdir in subdirectories
    ]
    if not any(candidate_flags):
        return []

    entries: list[dict[str, str]] = []
    seen_ids: set[str] = set()
    for subdir in subdirectories:
        _require_portable_filename_component(
            subdir.name,
            label="Dataset directory",
        )
        entry = _read_exported_dataset_entry(
            subdir,
            public_path=f"/{subdir.name}/",
        )
        if entry["id"] in seen_ids:
            raise ValueError(f"duplicate dataset id {entry['id']!r}.")
        seen_ids.add(entry["id"])
        entries.append(entry)
    return entries


@dataclass(frozen=True)
class _PreparedArtifact:
    path: Path
    size: int
    mtime_ns: int
    device: int
    inode: int
    content_type: str


def _require_artifact_path(value: object, *, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise TypeError(f"{label} must be a non-empty string.")
    if value != value.strip() or value.startswith("/") or "\\" in value:
        raise ValueError(f"{label} must be one exact relative POSIX path.")
    parts = value.split("/")
    if any(
        part in {"", ".", ".."}
        or any(ord(character) < 32 or ord(character) == 127 for character in part)
        for part in parts
    ):
        raise ValueError(f"{label} must be one exact relative POSIX path.")
    return value


def _read_json_object(path: Path, *, label: str) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} must contain readable UTF-8 JSON.") from error
    if not isinstance(value, dict):
        raise TypeError(f"{label} must contain a JSON object.")
    return value


def _expand_artifact_pattern(
    pattern: object,
    *,
    key: str,
    label: str,
    extension: str | None = None,
) -> str:
    pattern_value = _require_artifact_path(pattern, label=label)
    value = pattern_value.replace("{key}", key)
    if extension is not None:
        value = value.replace("{ext}", extension)
    if "{" in value or "}" in value:
        raise ValueError(f"{label} contains an unsupported placeholder.")
    return _require_artifact_path(value, label=label)


def _declared_dataset_artifacts(
    dataset_dir: Path,
    *,
    include_published_state: bool,
) -> set[str]:
    """Return every artifact declared by one current prepared generation."""
    identity_path = dataset_dir / "dataset_identity.json"
    obs_manifest_path = dataset_dir / "obs_manifest.json"
    identity = _read_json_object(identity_path, label="dataset_identity.json")
    obs_manifest = _read_json_object(obs_manifest_path, label="obs_manifest.json")
    paths = {"dataset_identity.json", "obs_manifest.json"}

    embeddings = identity.get("embeddings")
    if not isinstance(embeddings, dict) or not isinstance(embeddings.get("files"), dict):
        raise ValueError("dataset_identity.json embeddings.files must be a JSON object.")
    for dimension, artifact_path in embeddings["files"].items():
        paths.add(
            _require_artifact_path(
                artifact_path,
                label=f"dataset_identity.json embeddings.files.{dimension}",
            )
        )

    schemas = obs_manifest.get("_obsSchemas")
    continuous_fields = obs_manifest.get("_continuousFields")
    categorical_fields = obs_manifest.get("_categoricalFields")
    if (
        not isinstance(schemas, dict)
        or not isinstance(continuous_fields, list)
        or not isinstance(categorical_fields, list)
    ):
        raise ValueError("obs_manifest.json must declare exact compact field schemas.")

    continuous_schema = schemas.get("continuous")
    if continuous_fields and not isinstance(continuous_schema, dict):
        raise ValueError("Continuous observation fields require their exact schema.")
    if not continuous_fields and continuous_schema is not None:
        raise ValueError("Continuous observation schema exists without any fields.")
    continuous_schema_values = continuous_schema if isinstance(continuous_schema, dict) else {}
    for index, field in enumerate(continuous_fields):
        if not isinstance(field, list) or not field:
            raise ValueError(f"obs_manifest.json continuous field {index} is invalid.")
        key = _require_portable_filename_component(
            field[0],
            label=f"Observation field {index}",
        )
        paths.add(
            _expand_artifact_pattern(
                continuous_schema_values.get("pathPattern"),
                key=key,
                label="obs continuous pathPattern",
            )
        )

    categorical_schema = schemas.get("categorical")
    if categorical_fields and not isinstance(categorical_schema, dict):
        raise ValueError("Categorical observation fields require their exact schema.")
    if not categorical_fields and categorical_schema is not None:
        raise ValueError("Categorical observation schema exists without any fields.")
    categorical_schema_values = categorical_schema if isinstance(categorical_schema, dict) else {}
    dtype_extensions = {"uint8": "u8", "uint16": "u16"}
    for index, field in enumerate(categorical_fields):
        if not isinstance(field, list) or len(field) < 3:
            raise ValueError(f"obs_manifest.json categorical field {index} is invalid.")
        key = _require_portable_filename_component(
            field[0],
            label=f"Observation field {index}",
        )
        dtype = field[2]
        if dtype not in dtype_extensions:
            raise ValueError(f"obs_manifest.json categorical field {index} has an invalid dtype.")
        paths.add(
            _expand_artifact_pattern(
                categorical_schema_values.get("codesPathPattern"),
                key=key,
                extension=dtype_extensions[dtype],
                label="obs categorical codesPathPattern",
            )
        )
        outlier_pattern = categorical_schema_values.get("outlierPathPattern")
        if outlier_pattern is not None:
            paths.add(
                _expand_artifact_pattern(
                    outlier_pattern,
                    key=key,
                    label="obs categorical outlierPathPattern",
                )
            )

    stats = identity.get("stats")
    if not isinstance(stats, dict):
        raise ValueError("dataset_identity.json stats must be a JSON object.")
    n_genes = stats.get("n_genes")
    if type(n_genes) is not int or n_genes < 0:
        raise ValueError("dataset_identity.json stats.n_genes must be non-negative.")
    var_manifest_path = dataset_dir / "var_manifest.json"
    if n_genes > 0:
        paths.add("var_manifest.json")
        var_manifest = _read_json_object(var_manifest_path, label="var_manifest.json")
        var_schema = var_manifest.get("_varSchema")
        fields = var_manifest.get("fields")
        if not isinstance(var_schema, dict) or not isinstance(fields, list):
            raise ValueError("var_manifest.json must declare one exact field schema.")
        for index, field in enumerate(fields):
            if not isinstance(field, list) or not field:
                raise ValueError(f"var_manifest.json field {index} is invalid.")
            key = _require_portable_filename_component(
                field[0],
                label=f"Gene field {index}",
            )
            paths.add(
                _expand_artifact_pattern(
                    var_schema.get("pathPattern"),
                    key=key,
                    label="var pathPattern",
                )
            )
        if len(fields) != n_genes:
            raise ValueError("var_manifest.json field count must match identity stats.n_genes.")
    elif var_manifest_path.exists():
        raise ValueError("var_manifest.json is present while identity stats.n_genes is zero.")

    has_connectivity = stats.get("has_connectivity")
    if type(has_connectivity) is not bool:
        raise ValueError("dataset_identity.json stats.has_connectivity must be a boolean.")
    connectivity_manifest_path = dataset_dir / "connectivity_manifest.json"
    if has_connectivity:
        paths.add("connectivity_manifest.json")
        connectivity = _read_json_object(
            connectivity_manifest_path,
            label="connectivity_manifest.json",
        )
        for key in ("sourcesPath", "destinationsPath", "weightsPath"):
            paths.add(
                _require_artifact_path(
                    connectivity.get(key),
                    label=f"connectivity_manifest.json {key}",
                )
            )
    elif connectivity_manifest_path.exists():
        raise ValueError("connectivity_manifest.json is present while connectivity is disabled.")

    vector_fields = identity.get("vector_fields")
    if vector_fields is not None:
        if not isinstance(vector_fields, dict) or not isinstance(
            vector_fields.get("fields"),
            dict,
        ):
            raise ValueError("dataset_identity.json vector_fields.fields must be an object.")
        for field_id, field in vector_fields["fields"].items():
            if not isinstance(field, dict) or not isinstance(field.get("files"), dict):
                raise ValueError(f"dataset_identity.json vector field {field_id!r} is invalid.")
            for dimension, artifact_path in field["files"].items():
                paths.add(
                    _require_artifact_path(
                        artifact_path,
                        label=f"vector field {field_id!r} file {dimension}",
                    )
                )
    if include_published_state:
        paths.add(_PUBLISHED_STATE_MANIFEST)
        paths.add(_PUBLISHED_STATE_BUNDLE)
    return paths


def _build_prepared_artifact_inventory(
    data_dir: Path,
    datasets: list[dict[str, str]],
) -> dict[str, _PreparedArtifact]:
    root = data_dir.resolve(strict=True)
    inventory: dict[str, _PreparedArtifact] = {}
    for dataset in datasets:
        public_path = dataset["path"]
        prefix = "" if public_path == "/" else public_path.strip("/")
        dataset_dir = root if not prefix else root / prefix
        has_state_manifest = "state_manifest" in dataset
        has_state_sha256 = "state_sha256" in dataset
        if has_state_manifest != has_state_sha256:
            raise ValueError(
                "Prepared dataset state_manifest and state_sha256 must be declared together."
            )
        if has_state_manifest:
            current_state = _read_optional_published_state(dataset_dir)
            expected_state = {
                "state_manifest": dataset["state_manifest"],
                "state_sha256": dataset["state_sha256"],
            }
            if current_state != expected_state:
                raise ValueError(
                    "Published sample state changed while the artifact inventory was built."
                )
        for relative_path in sorted(
            _declared_dataset_artifacts(
                dataset_dir,
                include_published_state=has_state_manifest,
            )
        ):
            request_path = relative_path if not prefix else f"{prefix}/{relative_path}"
            candidate = dataset_dir / relative_path
            current = dataset_dir
            for part in Path(relative_path).parts:
                current = current / part
                if current.is_symlink():
                    raise ValueError(
                        f"Declared artifact must not traverse a symbolic link: {request_path}"
                    )
            try:
                resolved = candidate.resolve(strict=True)
                resolved.relative_to(root)
            except (FileNotFoundError, ValueError) as error:
                raise ValueError(
                    f"Declared artifact is missing or outside the export root: {request_path}"
                ) from error
            metadata = resolved.stat()
            if not stat.S_ISREG(metadata.st_mode):
                raise ValueError(f"Declared artifact must be a regular file: {request_path}")
            if request_path in inventory:
                raise ValueError(f"Prepared artifact path is duplicated: {request_path}")
            inventory[request_path] = _PreparedArtifact(
                path=resolved,
                size=metadata.st_size,
                mtime_ns=metadata.st_mtime_ns,
                device=metadata.st_dev,
                inode=metadata.st_ino,
                content_type=(
                    "application/json"
                    if relative_path.endswith(".json")
                    else "application/octet-stream"
                ),
            )
    return inventory


def _parse_byte_range(value: str | None, size: int) -> tuple[int, int] | None:
    if value is None:
        return None
    if type(size) is not int or size < 0:
        raise ValueError("Artifact size must be a non-negative integer.")
    match = re.fullmatch(r"bytes=(\d{0,20})-(\d{0,20})", value)
    if match is None or (not match.group(1) and not match.group(2)):
        raise ValueError("Range must contain one exact byte interval.")
    if size == 0:
        raise ValueError("An empty artifact has no satisfiable byte range.")
    if match.group(1):
        start = int(match.group(1))
        end = int(match.group(2)) if match.group(2) else size - 1
        if start >= size or end < start:
            raise ValueError("Range is outside the artifact.")
        return start, min(end, size - 1)
    suffix_length = int(match.group(2))
    if suffix_length <= 0:
        raise ValueError("Range suffix length must be positive.")
    return max(0, size - suffix_length), size - 1


class CORSRequestHandler(CORSMixin, SimpleHTTPRequestHandler):
    """Serve only endpoints and immutable prepared artifacts in the active contract."""

    allow_caching = True  # Static files can be cached

    def __init__(
        self,
        *args,
        data_dir: Path,
        server_info: dict,
        datasets: list[dict[str, str]],
        artifact_inventory: dict[str, _PreparedArtifact],
        serve_web_ui: bool,
        web_cache_dir: Path,
        **kwargs,
    ):
        self.data_dir = data_dir
        self.server_info = server_info
        self.datasets = datasets
        self.artifact_inventory = artifact_inventory
        self.serve_web_ui = serve_web_ui
        self.web_cache_dir = web_cache_dir
        # Must call super().__init__ last because it calls do_GET immediately
        super().__init__(*args, directory=str(data_dir), **kwargs)

    def end_headers(self):
        """Add CORS headers to every response."""
        self.add_cors_headers()
        super().end_headers()

    def do_POST(self):
        """Handle POST requests (events from frontend)."""
        if self.handle_event_post():
            return
        if self.handle_session_bundle_post():
            return
        # No other POST endpoints - return 404
        self.send_error_response(404, f"POST not supported for path: {self.path}")

    def do_GET(self):
        """Serve one exact GET endpoint or prepared artifact."""
        self._handle_read_request(head_only=False)

    def do_HEAD(self):
        """Serve GET metadata without a response body."""
        self._handle_read_request(head_only=True)

    def _handle_read_request(self, *, head_only: bool) -> None:
        path = self._canonical_request_path()
        if path is None:
            self.send_error_response(404, "Request path is not in the active contract")
            return

        if path == "/_cellucid/health":
            self.send_json(
                {
                    "status": "ok",
                    "type": "exported",
                    "version": self.server_info["version"],
                },
                head_only=head_only,
            )
            return

        if path == "/_cellucid/info":
            self.send_json(self.server_info, head_only=head_only)
            return

        if path == "/_cellucid/datasets":
            self.send_json({"datasets": self.datasets}, head_only=head_only)
            return

        if self.serve_web_ui and self.handle_web_asset_get(path, head_only=head_only):
            return
        if self.serve_web_ui and (
            path == f"/{WEB_ASSET_INVENTORY_FILENAME}"
            or path == "/assets"
            or path.startswith("/assets/")
        ):
            self.send_error_response(404, "Web asset is not declared by the active build")
            return

        if path in {"/", "/index.html"}:
            self.send_error_response(503, "Cellucid viewer UI unavailable")
            return

        self._serve_prepared_artifact(path[1:], head_only=head_only)

    def _canonical_request_path(self) -> str | None:
        """Return one exact ASCII origin-form path without decoding aliases."""
        parsed = urlparse(self.path)
        exact_target = parsed.path
        if parsed.query:
            exact_target += f"?{parsed.query}"
        if (
            parsed.scheme
            or parsed.netloc
            or parsed.params
            or parsed.fragment
            or exact_target != self.path
            or (parsed.query and parsed.path != "/")
            or not parsed.path.startswith("/")
            or "\\" in parsed.path
            or "%" in parsed.path
        ):
            return None
        try:
            parsed.path.encode("ascii")
        except UnicodeEncodeError:
            return None
        if any(
            ord(character) < 33 or ord(character) == 127
            for character in parsed.path
            if character != "/"
        ):
            return None
        if parsed.path != "/":
            parts = parsed.path[1:].split("/")
            if any(part in {"", ".", ".."} for part in parts):
                return None
        return parsed.path

    @staticmethod
    def _artifact_metadata_matches(
        metadata: os.stat_result,
        artifact: _PreparedArtifact,
    ) -> bool:
        return (
            stat.S_ISREG(metadata.st_mode)
            and metadata.st_size == artifact.size
            and metadata.st_mtime_ns == artifact.mtime_ns
            and metadata.st_dev == artifact.device
            and metadata.st_ino == artifact.inode
        )

    def _open_prepared_artifact(self, artifact: _PreparedArtifact) -> int:
        """Open one unchanged regular artifact without following its final symlink."""
        try:
            if artifact.path.is_symlink() or artifact.path.resolve(strict=True) != artifact.path:
                raise OSError("Prepared artifact path changed.")
            before = artifact.path.stat()
            if not self._artifact_metadata_matches(before, artifact):
                raise OSError("Prepared artifact metadata changed.")

            flags = os.O_RDONLY | getattr(os, "O_BINARY", 0)
            nofollow = getattr(os, "O_NOFOLLOW", 0)
            if nofollow:
                flags |= nofollow
            descriptor = os.open(artifact.path, flags)
        except (OSError, RuntimeError):
            raise

        try:
            opened = os.fstat(descriptor)
            after = artifact.path.stat()
            if not self._artifact_metadata_matches(
                opened, artifact
            ) or not self._artifact_metadata_matches(after, artifact):
                raise OSError("Prepared artifact changed while it was opened.")
        except BaseException:
            os.close(descriptor)
            raise
        return descriptor

    def _serve_prepared_artifact(self, request_path: str, *, head_only: bool) -> None:
        artifact = self.artifact_inventory.get(request_path)
        if artifact is None:
            self.send_error_response(404, "Prepared artifact is not declared")
            return

        range_values = self._request_headers().get_all("Range", [])
        if len(range_values) > 1:
            self._send_range_error(artifact.size)
            return
        try:
            interval = _parse_byte_range(
                range_values[0] if range_values else None,
                artifact.size,
            )
        except ValueError:
            self._send_range_error(artifact.size)
            return

        try:
            descriptor = self._open_prepared_artifact(artifact)
        except OSError:
            self.send_error_response(
                HTTPStatus.CONFLICT,
                "Prepared artifact changed after server validation",
            )
            return

        start, end = interval if interval is not None else (0, artifact.size - 1)
        content_length = 0 if artifact.size == 0 else end - start + 1
        self.send_response(HTTPStatus.PARTIAL_CONTENT if interval is not None else HTTPStatus.OK)
        self.send_header("Content-Type", artifact.content_type)
        self.send_header("Content-Length", str(content_length))
        self.send_header("Accept-Ranges", "bytes")
        if interval is not None:
            self.send_header(
                "Content-Range",
                f"bytes {start}-{end}/{artifact.size}",
            )
        self.end_headers()

        try:
            if not head_only and content_length:
                os.lseek(descriptor, start, os.SEEK_SET)
                remaining = content_length
                while remaining:
                    chunk = os.read(descriptor, min(1024 * 1024, remaining))
                    if not chunk:
                        raise OSError("Prepared artifact ended during response.")
                    self._response_writer().write(chunk)
                    remaining -= len(chunk)
        finally:
            os.close(descriptor)

    def _send_range_error(self, artifact_size: int) -> None:
        body = b"Requested byte range is not satisfiable"
        self.send_response(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Content-Range", f"bytes */{artifact_size}")
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        if self.command != "HEAD":
            self._response_writer().write(body)

    def _list_datasets(self) -> list[dict]:
        """Validate and return the prepared datasets under ``data_dir``."""
        return _list_exported_datasets(self.data_dir)

    def _is_dataset_dir(self, path: Path) -> bool:
        """Check if a directory is a valid cellucid dataset."""
        try:
            _read_exported_dataset_entry(path, public_path="/")
        except (OSError, TypeError, ValueError):
            return False
        return True

    def _get_dataset_identity_fields(self, path: Path) -> tuple[str, str]:
        """Return (dataset_id, dataset_name) for a dataset directory."""
        entry = _read_exported_dataset_entry(path, public_path="/")
        return entry["id"], entry["name"]

    def log_message(self, format: str, *args):
        """Override to use Python logging instead of stderr."""
        logger.debug("%s - %s", self.address_string(), format % args)


class CellucidServer:
    """
    Cellucid data server for serving datasets over HTTP.

    Supports multiple deployment modes:
    - Local: Direct browser access on localhost
    - SSH tunnel: Access via port forwarding from remote server
    - Jupyter: Embedded in notebook environment

    Example:
        server = CellucidServer("/path/to/data")
        server.start()  # Blocking

        # Or non-blocking:
        server.start_background()
        # ... do other things ...
        server.stop()
    """

    def __init__(
        self,
        data_dir: str | Path,
        port: int = DEFAULT_PORT,
        host: str = DEFAULT_HOST,
        open_browser: bool = False,
        quiet: bool = False,
        *,
        serve_web_ui: bool = True,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
    ):
        """
        Initialize the server.

        Args:
            data_dir: Path to the dataset directory (single dataset or multi-dataset)
            port: Port to serve on (default: 8765)
            host: Host to bind to (default: 127.0.0.1 for localhost only)
            open_browser: Whether to open the browser on start
            quiet: Suppress info messages
            serve_web_ui: Establish and serve the exact current web build.
            web_source_url: Origin publishing the web asset inventory.
            web_cache_dir: Directory holding the active verified web build.
        """
        self.data_dir = Path(data_dir).resolve()
        self.port = require_server_port(port)
        self.host = host
        if type(open_browser) is not bool:
            raise TypeError("open_browser must be a boolean")
        if type(quiet) is not bool:
            raise TypeError("quiet must be a boolean")
        self.open_browser = open_browser
        self.quiet = quiet
        if type(serve_web_ui) is not bool:
            raise TypeError("serve_web_ui must be a boolean")
        self.serve_web_ui = serve_web_ui
        self.web_source_url = web_source_url
        self.web_cache_dir = (
            Path(web_cache_dir).expanduser().resolve()
            if web_cache_dir is not None
            else _web_cache_dir()
        )

        # Step 1: Validate dataset
        if not quiet:
            print_step(1, 3, "Validating dataset")
            print_detail("Path", str(self.data_dir))

        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {self.data_dir}")
        if not self.data_dir.is_dir():
            raise NotADirectoryError(f"Prepared-data path must be a directory: {self.data_dir}")

        self._datasets = _list_exported_datasets(self.data_dir)
        if not self._datasets:
            raise ValueError("Prepared-data path does not contain one complete current dataset.")
        self._artifact_inventory = _build_prepared_artifact_inventory(
            self.data_dir,
            self._datasets,
        )

        if not quiet:
            print_success("Dataset valid")

        # Step 2: Load dataset info
        if not quiet:
            print_step(2, 3, "Loading dataset info")
            self._print_dataset_info()
            print_success("Dataset loaded")

        self._server: HTTPServer | None = None
        self._thread: threading.Thread | None = None
        self._running = False
        self._started = False
        self._closed = False
        self._serving = False
        self._background_error: BaseException | None = None

        from . import __version__

        self.server_info = {
            "version": __version__,
            "host": self.host,
            "port": self.port,
            "mode": "standalone",
        }

    def _print_dataset_info(self):
        """Print information about the dataset."""
        if len(self._datasets) != 1 or self._datasets[0]["path"] != "/":
            print_detail("Datasets", str(len(self._datasets)))
            return
        identity = _read_json_object(
            self.data_dir / "dataset_identity.json",
            label="dataset_identity.json",
        )
        stats = identity.get("stats")
        if not isinstance(stats, dict):
            raise TypeError("dataset_identity.json stats must be a JSON object.")
        n_cells = stats.get("n_cells")
        n_genes = stats.get("n_genes")
        has_connectivity = stats.get("has_connectivity")
        if type(n_cells) is not int or n_cells < 0:
            raise ValueError("dataset_identity.json stats.n_cells must be non-negative.")
        if type(n_genes) is not int or n_genes < 0:
            raise ValueError("dataset_identity.json stats.n_genes must be non-negative.")
        if type(has_connectivity) is not bool:
            raise ValueError("dataset_identity.json stats.has_connectivity must be a boolean.")
        print_detail("Cells", f"{n_cells:,}")
        print_detail("Genes", f"{n_genes:,}")
        print_detail("Connectivity", "yes" if has_connectivity else "no")

    @property
    def url(self) -> str:
        """Get the URL of the currently running server."""
        if not self._running or self._server is None:
            raise RuntimeError(
                "CellucidServer URL is unavailable because the server is not running"
            )
        return f"http://{self.host}:{self.port}"

    @property
    def viewer_url(self) -> str:
        """Open this server's prepared-data catalog in the verified viewer.

        The viewer selects the sole declared dataset when the catalog is
        unique. A multi-dataset catalog requires an exact dataset selection;
        this URL never embeds or guesses an arbitrary catalog entry.
        """
        return f"{self.url}/?source=remote"

    def start(self, blocking: bool = True):
        """Start this single-use server."""
        if type(blocking) is not bool:
            raise TypeError("blocking must be a boolean")
        if self._running:
            raise RuntimeError("Server is already running.")
        if self._started or self._closed:
            raise RuntimeError(
                "CellucidServer is single-use and has been closed. Create a new server instance."
            )
        self._started = True

        try:
            # Step 3: Start server
            if not self.quiet:
                print_step(3, 3, "Starting server")

            if self.serve_web_ui:
                from .web_cache import ensure_web_ui_cached

                if not self.quiet:
                    print_detail(
                        "Viewer UI generation",
                        "establishing exact configured source",
                    )
                ensure_web_ui_cached(
                    cache_dir=self.web_cache_dir,
                    source_url=self.web_source_url,
                    force=True,
                    show_progress=not self.quiet,
                )
                if not self.quiet:
                    print_success("Viewer UI generation established")

            handler = partial(
                CORSRequestHandler,
                data_dir=self.data_dir,
                server_info=self.server_info,
                datasets=self._datasets,
                artifact_inventory=self._artifact_inventory,
                serve_web_ui=self.serve_web_ui,
                web_cache_dir=self.web_cache_dir,
            )

            self._server = HTTPServer((self.host, self.port), handler)
            self.port = require_server_port(self._server.server_address[1])
            self.server_info["port"] = self.port
            self._running = True

            if not self.quiet:
                print_success("Server ready")
                print_server_banner(self.url, self.viewer_url)

            if blocking:
                if self.open_browser and webbrowser.open(self.viewer_url) is not True:
                    raise RuntimeError(f"Could not open the browser for {self.viewer_url}")
                self._serving = True
                try:
                    self._server.serve_forever()
                finally:
                    self._serving = False
                self._finish_serving()
            else:
                self._serve_entered = threading.Event()
                self._thread = threading.Thread(
                    target=self._serve_in_background,
                    daemon=True,
                )
                self._thread.start()
                self._serve_entered.wait()
                if self.open_browser and webbrowser.open(self.viewer_url) is not True:
                    raise RuntimeError(f"Could not open the browser for {self.viewer_url}")
        except BaseException:
            self._rollback_failed_start(shutdown=self._thread is not None)
            raise

    def _serve_in_background(self) -> None:
        """Run the bound server and retain an exact asynchronous failure."""
        self._serving = True
        self._serve_entered.set()
        serving_error: BaseException | None = None
        try:
            server = self._server
            if server is None:
                raise RuntimeError("CellucidServer lost its bound HTTP server before serving")
            server.serve_forever()
        except BaseException as error:
            serving_error = error
        finally:
            self._serving = False
            try:
                self._finish_serving()
            except BaseException as cleanup_error:
                if serving_error is None:
                    serving_error = cleanup_error
                else:
                    logger.exception(
                        "Prepared-data server cleanup also failed after its serving loop failed"
                    )
        self._background_error = serving_error
        if serving_error is not None:
            logger.error(
                "Prepared-data background server failed",
                exc_info=(
                    type(serving_error),
                    serving_error,
                    serving_error.__traceback__,
                ),
            )

    def _finish_serving(self) -> None:
        """Close the socket after the serving loop has ended."""
        self._running = False
        server = self._server
        self._server = None
        self._closed = True
        if server is not None:
            try:
                server.server_close()
            except BaseException as error:
                raise RuntimeError(
                    f"Prepared-data server cleanup failed: {type(error).__name__}: {error}"
                ) from error

    def _rollback_failed_start(self, *, shutdown: bool) -> None:
        """Release acquired resources without replacing the startup exception."""
        self._running = False
        server = self._server
        if server is not None and shutdown and self._serving:
            try:
                server.shutdown()
            except BaseException:
                logger.exception(
                    "Failed to shut down the prepared-data server after startup failed"
                )
        if server is not None:
            try:
                server.server_close()
            except BaseException:
                logger.exception("Failed to close the prepared-data socket after startup failed")
        thread = self._thread
        if thread is not None and thread is not threading.current_thread() and thread.is_alive():
            thread.join()
        self._server = None
        self._thread = None
        self._serving = False
        self._closed = True

    def start_background(self):
        """Start the server in a background thread."""
        self.start(blocking=False)

    def stop(self):
        """Stop this server and release its socket."""
        self._running = False
        failures: list[BaseException] = []
        server = self._server
        thread = self._thread

        if server is not None and self._serving:
            try:
                server.shutdown()
            except BaseException as error:
                failures.append(error)
        if server is not None:
            try:
                server.server_close()
            except BaseException as error:
                failures.append(error)
        self._server = None

        if thread is not None and thread is not threading.current_thread() and thread.is_alive():
            thread.join()
        self._thread = None
        self._serving = False
        self._closed = True

        if not self.quiet:
            print("Server stopped")
        if failures:
            details = "; ".join(f"{type(error).__name__}: {error}" for error in failures)
            raise RuntimeError(f"Prepared-data server shutdown failed: {details}") from failures[0]

    def is_running(self) -> bool:
        """Check if the server is running."""
        return self._running

    def wait(self):
        """Wait for the background server to stop."""
        thread = self._thread
        if thread is not None:
            try:
                thread.join()
            except BaseException:
                self.stop()
                raise
        if self._background_error is not None:
            raise self._background_error


def serve(
    data_dir: str | Path,
    port: int = DEFAULT_PORT,
    host: str = DEFAULT_HOST,
    open_browser: bool = True,
    quiet: bool = False,
    *,
    serve_web_ui: bool = True,
    web_source_url: str = CELLUCID_WEB_URL,
    web_cache_dir: str | Path | None = None,
):
    """
    Serve a cellucid dataset directory.

    This is the main entry point for serving data. It starts an HTTP server
    that serves the dataset files with CORS headers enabled.

    Args:
        data_dir: Path to the dataset directory
        port: Port to serve on (default: 8765)
        host: Host to bind to (default: 127.0.0.1)
        open_browser: Whether to open the viewer in browser (default: True)
        quiet: Suppress info messages
        serve_web_ui: Establish and serve the exact current web build.
        web_source_url: Origin publishing the web asset inventory.
        web_cache_dir: Directory holding the active verified web build.

    Example:
        >>> from cellucid import serve
        >>> serve("/path/to/my_dataset")

        # For remote server access via SSH:
        >>> serve("/path/to/data", host="0.0.0.0")
        # Then on local machine: ssh -L 8765:localhost:8765 remote-server
    """
    server = CellucidServer(
        data_dir=data_dir,
        port=port,
        host=host,
        open_browser=open_browser,
        quiet=quiet,
        serve_web_ui=serve_web_ui,
        web_source_url=web_source_url,
        web_cache_dir=web_cache_dir,
    )
    server.start()
