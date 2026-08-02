"""
Reading the published-state sidecars an export directory may carry.

A prepared dataset can ship a default session: a tiny manifest naming one
bundle, and the bundle itself. Both are read here under exact limits and
without following links -- the file is stat-ed before and after the read and
rejected unless it is the same unchanged regular file, so a sidecar cannot be
swapped for a symlink out of the served root between the check and the read.
"""

from __future__ import annotations

import hashlib
import json
import os
import stat
from pathlib import Path

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
        raise ValueError(f"Published sample state sidecar is unreadable: {filename}") from error
    if stat.S_ISLNK(before.st_mode):
        raise ValueError(f"Published sample state sidecar must not be a symbolic link: {filename}")
    if not stat.S_ISREG(before.st_mode):
        raise ValueError(f"Published sample state sidecar must be a regular file: {filename}")
    if before.st_size < minimum_bytes or before.st_size > maximum_bytes:
        if minimum_bytes == maximum_bytes:
            raise ValueError(
                f"Published sample state {filename} must be exactly {minimum_bytes} bytes."
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
            raise ValueError(f"Published sample state sidecar grew during validation: {filename}")
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
        if child.name.endswith(".cellucid-session") and child.name != _PUBLISHED_STATE_BUNDLE
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
            f"Published sample state {_PUBLISHED_STATE_MANIFEST} must contain readable UTF-8 JSON."
        ) from error
    if not isinstance(manifest, dict):
        raise TypeError(
            f"Published sample state {_PUBLISHED_STATE_MANIFEST} must contain a JSON object."
        )
    if set(manifest) != {"states"}:
        raise ValueError(
            f"Published sample state {_PUBLISHED_STATE_MANIFEST} "
            "must contain the exact keys: states."
        )
    if type(manifest["states"]) is not list or manifest["states"] != [_PUBLISHED_STATE_BUNDLE]:
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
