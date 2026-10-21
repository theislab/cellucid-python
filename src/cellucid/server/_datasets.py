"""
Discovering the exported datasets a served directory declares.

A served root is either one export directory or a directory of them. An export
directory is recognised by its ``dataset_identity.json``, ``obs_manifest.json``
and ``embeddings/`` directory. Each candidate is admitted only with a readable
version-3 ``dataset_identity.json`` whose id and name pass the same contracts
the writer applied, and whose directory name is a portable filename component,
because that name becomes the URL path the viewer requests. Which embedding
files it must hold is decided by its declaration, which the artifact inventory
reads.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import NotRequired, TypedDict

from .._contracts import (
    EMBEDDINGS_DIRNAME,
    IDENTITY_VERSION,
    _require_dataset_id,
    _require_dataset_name,
    _require_portable_filename_component,
)
from ._state import _read_optional_published_state


class _DatasetEntry(TypedDict):
    """One ``/_cellucid/datasets`` listing entry.

    The three state fields are present together or not at all: they advertise
    a published default session by its manifest name, digest and exact length.
    """

    id: str
    path: str
    name: str
    state_manifest: NotRequired[str]
    state_sha256: NotRequired[str]
    state_bytes: NotRequired[int]


def _read_exported_dataset_entry(
    dataset_dir: Path,
    *,
    public_path: str,
) -> _DatasetEntry:
    """Validate one complete prepared-dataset root and return its catalog entry."""
    required_files = (
        dataset_dir / "dataset_identity.json",
        dataset_dir / "obs_manifest.json",
    )
    if not dataset_dir.is_dir() or any(not path.is_file() for path in required_files):
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
    if type(identity.get("version")) is not int or identity["version"] != IDENTITY_VERSION:
        raise ValueError(
            f"dataset_identity.json version must be exactly {IDENTITY_VERSION}. "
            "Re-export the dataset with the current prepare() or cellucid_prepare()."
        )
    # Checked after the version, so an export written in an earlier format is
    # told to re-export rather than that a directory is missing.
    embeddings_dir = dataset_dir / EMBEDDINGS_DIRNAME
    if (
        embeddings_dir.is_symlink()
        or not embeddings_dir.is_dir()
        or not any(embeddings_dir.iterdir())
    ):
        raise ValueError(
            f"{dataset_dir.name!r} is not a complete exported dataset: it has no "
            f"{EMBEDDINGS_DIRNAME}/ payloads."
        )

    dataset_id = _require_dataset_id(
        identity.get("id"),
        label="dataset_id",
    )
    dataset_name = _require_dataset_name(
        identity.get("name"),
        label="dataset_name",
    )

    entry: _DatasetEntry = {
        "id": dataset_id,
        "path": public_path,
        "name": dataset_name,
    }
    state = _read_optional_published_state(dataset_dir)
    if state is not None:
        entry["state_manifest"] = state["state_manifest"]
        entry["state_sha256"] = state["state_sha256"]
        entry["state_bytes"] = state["state_bytes"]
    return entry


def _declares_exported_dataset(directory: Path) -> bool:
    """Whether a directory holds any marker of one export generation."""
    return any(
        marker.exists()
        for marker in (
            directory / "dataset_identity.json",
            directory / "obs_manifest.json",
            directory / EMBEDDINGS_DIRNAME,
        )
    )


def _list_exported_datasets(data_dir: Path) -> list[_DatasetEntry]:
    """Classify one exact prepared dataset or a root containing only datasets."""
    data_dir = Path(data_dir).expanduser()
    if not data_dir.is_dir():
        raise NotADirectoryError(f"Prepared-data path must be a directory: {data_dir}")

    if _declares_exported_dataset(data_dir):
        return [_read_exported_dataset_entry(data_dir, public_path="/")]

    subdirectories = sorted(path for path in data_dir.iterdir() if path.is_dir())
    if not subdirectories:
        return []

    candidate_flags = [_declares_exported_dataset(subdir) for subdir in subdirectories]
    if not any(candidate_flags):
        return []

    entries: list[_DatasetEntry] = []
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
