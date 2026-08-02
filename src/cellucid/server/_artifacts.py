"""
The exact inventory of files one prepared dataset is allowed to serve.

Nothing under a served root is offered because it happens to exist: the
manifests are read, every payload they declare is expanded into an exact
relative path, and the inventory is the set of those paths. A request that does
not name one of them is a 404 even if the file is there, which is what keeps an
unrelated file dropped into an export directory unreachable.
"""

from __future__ import annotations

import json
import stat
from dataclasses import dataclass
from pathlib import Path

from ._state import (
    _PUBLISHED_STATE_BUNDLE,
    _PUBLISHED_STATE_MANIFEST,
    _read_optional_published_state,
)


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
    index: int,
    label: str,
    extension: str | None = None,
) -> str:
    pattern_value = _require_artifact_path(pattern, label=label)
    value = pattern_value.replace("{index}", str(index))
    if extension is not None:
        value = value.replace("{ext}", extension)
    if "{" in value or "}" in value:
        raise ValueError(f"{label} contains an unsupported placeholder.")
    return _require_artifact_path(value, label=label)


def _require_payload_index(value: object, *, label: str) -> int:
    """Require one declared payload index, which is also the payload filename."""
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} must declare a non-negative integer payload index.")
    return value


def _require_dense_payload_indices(indices: list[int], *, label: str) -> None:
    """Require one axis to declare exactly the indices 0..N-1, each once.

    Two entries sharing an index name one file, so one field's payload is read
    under the other field's name. The set of expanded paths hides that -- it
    simply loses one member -- so the indices are checked directly.
    """
    if sorted(indices) != list(range(len(indices))):
        raise ValueError(
            f"{label} payload indices must be exactly 0..{len(indices) - 1}, each used once."
        )


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
    # obs/ carries both arrays, so their payload indices share one space.
    obs_payload_indices: list[int] = []
    for position, field in enumerate(continuous_fields):
        if not isinstance(field, list) or len(field) < 2:
            raise ValueError(f"obs_manifest.json continuous field {position} is invalid.")
        payload_index = _require_payload_index(
            field[0],
            label=f"Observation field {position}",
        )
        obs_payload_indices.append(payload_index)
        paths.add(
            _expand_artifact_pattern(
                continuous_schema_values.get("pathPattern"),
                index=payload_index,
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
    for position, field in enumerate(categorical_fields):
        if not isinstance(field, list) or len(field) < 4:
            raise ValueError(f"obs_manifest.json categorical field {position} is invalid.")
        payload_index = _require_payload_index(
            field[0],
            label=f"Observation field {position}",
        )
        obs_payload_indices.append(payload_index)
        dtype = field[3]
        if dtype not in dtype_extensions:
            raise ValueError(f"obs_manifest.json categorical field {position} has an invalid dtype.")
        paths.add(
            _expand_artifact_pattern(
                categorical_schema_values.get("codesPathPattern"),
                index=payload_index,
                extension=dtype_extensions[dtype],
                label="obs categorical codesPathPattern",
            )
        )
        outlier_pattern = categorical_schema_values.get("outlierPathPattern")
        if outlier_pattern is not None:
            paths.add(
                _expand_artifact_pattern(
                    outlier_pattern,
                    index=payload_index,
                    label="obs categorical outlierPathPattern",
                )
            )
    _require_dense_payload_indices(obs_payload_indices, label="obs_manifest.json")

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
        var_payload_indices: list[int] = []
        for position, field in enumerate(fields):
            if not isinstance(field, list) or len(field) < 2:
                raise ValueError(f"var_manifest.json field {position} is invalid.")
            payload_index = _require_payload_index(
                field[0],
                label=f"Gene field {position}",
            )
            var_payload_indices.append(payload_index)
            paths.add(
                _expand_artifact_pattern(
                    var_schema.get("pathPattern"),
                    index=payload_index,
                    label="var pathPattern",
                )
            )
        _require_dense_payload_indices(var_payload_indices, label="var_manifest.json")
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
