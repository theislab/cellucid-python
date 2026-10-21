"""
Reading the published-state sidecars an export directory may carry.

A prepared dataset can ship a default session: a tiny manifest naming one
bundle, and the bundle itself. The manifest has exactly one valid form; the
bundle has no size limit, and its exact length and digest are what the listing
advertises. Both are read here without following links -- the file is stat-ed before and after the read and
rejected unless it is the same unchanged regular file, so a sidecar cannot be
swapped for a symlink out of the served root between the check and the read.
"""

from __future__ import annotations

import hashlib
import json
import os
import stat
from collections.abc import Callable
from pathlib import Path
from typing import TypedDict

_PUBLISHED_STATE_MANIFEST = "state-snapshots.json"
_PUBLISHED_STATE_BUNDLE = "default.cellucid-session"
# The one valid manifest, ``{ "states": ["default.cellucid-session"] }``, is
# exactly this long as written: a format constant, not a size policy.
_PUBLISHED_STATE_MANIFEST_BYTES = 43


class _PublishedState(TypedDict):
    """The listing fields that advertise one published default session."""

    state_manifest: str
    state_sha256: str
    state_bytes: int


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
    maximum_bytes: int | None,
    sink: Callable[[bytes], None],
) -> int:
    """Stream one unchanged in-root regular state sidecar into ``sink``.

    Links are not followed, and the file must be the same unchanged regular
    file before, during and after the read. Returns the byte length read.
    """
    candidate = dataset_root / filename
    try:
        before = candidate.lstat()
    except OSError as error:
        raise ValueError(f"Published sample state sidecar is unreadable: {filename}") from error
    if stat.S_ISLNK(before.st_mode):
        raise ValueError(f"Published sample state sidecar must not be a symbolic link: {filename}")
    if not stat.S_ISREG(before.st_mode):
        raise ValueError(f"Published sample state sidecar must be a regular file: {filename}")
    if before.st_size < minimum_bytes or (
        maximum_bytes is not None and before.st_size > maximum_bytes
    ):
        if minimum_bytes == maximum_bytes:
            raise ValueError(
                f"Published sample state {filename} must be exactly {minimum_bytes} bytes."
            )
        if maximum_bytes is None:
            raise ValueError(
                f"Published sample state {filename} must be at least {minimum_bytes} bytes."
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
        while remaining:
            chunk = os.read(descriptor, min(64 * 1024, remaining))
            if not chunk:
                raise ValueError(
                    f"Published sample state sidecar ended during validation: {filename}"
                )
            sink(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            raise ValueError(f"Published sample state sidecar grew during validation: {filename}")
        after = candidate.lstat()
        if not _same_regular_file(opened, after):
            raise ValueError(
                f"Published sample state sidecar changed during validation: {filename}"
            )
        return opened.st_size
    finally:
        os.close(descriptor)


def _read_optional_published_state(dataset_dir: Path) -> _PublishedState | None:
    """Validate and describe the exact optional published sample state pair."""
    dataset_root = dataset_dir.resolve(strict=True)
    manifest_path = dataset_root / _PUBLISHED_STATE_MANIFEST
    bundle_path = dataset_root / _PUBLISHED_STATE_BUNDLE
    manifest_present = manifest_path.is_symlink() or manifest_path.exists()
    bundle_present = bundle_path.is_symlink() or bundle_path.exists()
    if not manifest_present and not bundle_present:
        return None
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

    manifest_chunks: list[bytes] = []
    _read_published_state_file(
        dataset_root,
        _PUBLISHED_STATE_MANIFEST,
        minimum_bytes=_PUBLISHED_STATE_MANIFEST_BYTES,
        maximum_bytes=_PUBLISHED_STATE_MANIFEST_BYTES,
        sink=manifest_chunks.append,
    )
    manifest_bytes = b"".join(manifest_chunks)
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

    # Hashed as it streams, so a bundle of any size is described without
    # holding it in memory.
    bundle_digest = hashlib.sha256()
    bundle_bytes = _read_published_state_file(
        dataset_root,
        _PUBLISHED_STATE_BUNDLE,
        minimum_bytes=1,
        maximum_bytes=None,
        sink=bundle_digest.update,
    )
    return {
        "state_manifest": _PUBLISHED_STATE_MANIFEST,
        "state_sha256": bundle_digest.hexdigest(),
        "state_bytes": bundle_bytes,
    }
