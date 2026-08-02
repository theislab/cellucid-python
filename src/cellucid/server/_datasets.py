"""
Discovering the exported datasets a served directory declares.

A served root is either one export directory or a directory of them. Each
candidate is admitted only with a readable ``dataset_identity.json`` whose id
and name pass the same contracts the writer applied, and whose directory name
is a portable filename component, because that name becomes the URL path the
viewer requests.
"""

from __future__ import annotations

import json
from pathlib import Path

from .._contracts import (
    _require_dataset_id,
    _require_dataset_name,
    _require_portable_filename_component,
)
from ._state import _read_optional_published_state


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

    dataset_id = _require_dataset_id(
        identity.get("id"),
        label="dataset_id",
    )
    dataset_name = _require_dataset_name(
        identity.get("name"),
        label="dataset_name",
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
    data_dir = Path(data_dir).expanduser()
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
