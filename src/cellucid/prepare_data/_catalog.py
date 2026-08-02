"""
The datasets.json catalog over a directory of exports.

A catalog is not an export: it validates the identity files that separate
export runs already published, refuses two datasets that a URL or a
case-insensitive filesystem could confuse, and republishes the manifest the
web viewer reads to discover them.
"""

import json
import os
import re
import tempfile
from pathlib import Path

from .._console import console_print
from .._contracts import (
    _WINDOWS_RESERVED_COMPONENTS,
    _require_dataset_id,
    _require_dataset_name,
    _require_display_text,
)
from ._manifest import _require_manifest_reads_back


def generate_datasets_manifest(
    exports_dir: str | Path = "./exports",
    *,
    default_dataset: str,
) -> Path:
    """
    Validate an exports directory and atomically publish its datasets.json manifest.

    This utility helps maintain the datasets.json manifest file that the frontend uses
    to discover available demo datasets. Run this after adding or removing datasets.

    Parameters
    ----------
    exports_dir : Path or str
        Directory containing dataset subdirectories. Defaults to ``./exports``,
        resolved against the current working directory at call time. A leading
        ``~`` is expanded to the user's home directory.
    default_dataset : str
        Exact ID of the dataset selected by default.

    Returns
    -------
    Path
        Path to the generated datasets.json file.

    Example
    -------
    >>> from cellucid.prepare_data import generate_datasets_manifest
    >>> generate_datasets_manifest("./exports", default_dataset="my_dataset")
    """
    exports_dir = Path(exports_dir).expanduser()
    if not exports_dir.exists():
        raise FileNotFoundError(f"Exports directory not found: {exports_dir}")
    if not exports_dir.is_dir():
        raise NotADirectoryError(f"Exports path must be a directory: {exports_dir}")
    default_dataset = _require_dataset_id(
        default_dataset,
        label="default_dataset",
    )

    console_print(f"Scanning {exports_dir} for datasets...")

    datasets: list[dict[str, object]] = []
    dataset_ids: set[str] = set()
    directory_names: dict[str, str] = {}

    for subdir in sorted(exports_dir.iterdir()):
        if not subdir.is_dir():
            continue

        directory_name = subdir.name
        if (
            len(directory_name.encode("utf-8")) > 180
            or directory_name in {".", ".."}
            or re.fullmatch(r"[A-Za-z0-9_.-]+", directory_name, flags=re.ASCII) is None
            or directory_name.endswith(".")
            or directory_name.split(".", 1)[0].upper() in _WINDOWS_RESERVED_COMPONENTS
        ):
            raise ValueError(
                f"Dataset directory {directory_name!r} is not a portable URL/path component."
            )
        directory_collision_key = directory_name.casefold()
        if directory_collision_key in directory_names:
            raise ValueError(
                f"Dataset directories {directory_names[directory_collision_key]!r} and "
                f"{directory_name!r} collide case-insensitively."
            )
        directory_names[directory_collision_key] = directory_name

        identity_file = subdir / "dataset_identity.json"
        if not identity_file.is_file():
            raise ValueError(
                f"Dataset directory {directory_name!r} has no dataset_identity.json file."
            )

        try:
            identity = json.loads(identity_file.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise ValueError(f"{identity_file} must contain readable UTF-8 JSON.") from error
        if not isinstance(identity, dict):
            raise TypeError(f"{identity_file} must contain a JSON object.")
        if type(identity.get("version")) is not int or identity["version"] != 2:
            raise ValueError(f"{identity_file} version must be exactly 2.")

        dataset_id = _require_dataset_id(
            identity.get("id"),
            label=f"{identity_file} dataset_id",
        )
        if dataset_id in dataset_ids:
            raise ValueError(f"duplicate dataset id {dataset_id!r}.")
        dataset_ids.add(dataset_id)

        dataset_name = _require_dataset_name(
            identity.get("name"),
            label=f"{identity_file} dataset_name",
        )
        dataset_entry: dict[str, object] = {
            "id": dataset_id,
            "path": f"{directory_name}/",
            "name": dataset_name,
        }
        if "description" in identity:
            description = identity["description"]
            if not isinstance(description, str):
                raise TypeError(f"{identity_file} description must be a string.")
            dataset_entry["description"] = _require_display_text(
                description,
                label=f"{identity_file} description",
                allow_empty=True,
            )

        stats = identity.get("stats")
        if not isinstance(stats, dict):
            raise TypeError(f"{identity_file} stats must be a JSON object.")
        for count_key in ("n_cells", "n_genes"):
            count = stats.get(count_key)
            if type(count) is not int or count < 0 or count > (1 << 53) - 1:
                raise ValueError(
                    f"{identity_file} stats.{count_key} must be a non-negative safe integer."
                )
            dataset_entry[count_key] = count

        datasets.append(dataset_entry)
        console_print(f"  ✓ Found dataset: {dataset_entry['name']} ({dataset_entry['id']})")

    if not datasets:
        raise ValueError(f"Exports directory contains no datasets: {exports_dir}")
    if default_dataset not in dataset_ids:
        raise ValueError(
            f"default_dataset {default_dataset!r} is missing from the validated datasets."
        )

    manifest = {"version": 1, "default": default_dataset, "datasets": datasets}
    manifest_bytes = json.dumps(manifest, indent=2).encode("utf-8")
    manifest_path = exports_dir / "datasets.json"

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=".datasets.json.",
        suffix=".tmp",
        dir=exports_dir,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(manifest_bytes)
            output.flush()
            os.fsync(output.fileno())
        # The catalog is the first file the viewer reads, and it reaches the
        # browser through JSON.parse exactly as the four export manifests do.
        # Verified here, on the temporary file, so a catalog whose bytes lost a
        # value never replaces the one that is already published.
        _require_manifest_reads_back(temporary_path, manifest)
        os.replace(temporary_path, manifest_path)
    except BaseException:
        if temporary_path.exists():
            temporary_path.unlink()
        raise

    console_print(f"✓ Wrote datasets manifest with {len(datasets)} datasets to {manifest_path}")
    console_print(f"  Default dataset: {default_dataset}")

    return manifest_path
