"""
Durable filesystem primitives shared by the export transaction.

Every operation here either refuses a path that is not what it claims to be --
a symbolic link, a Windows reparse point, a hard-linked or irregular file -- or
orders a creation, rename, and removal durably enough that an interrupted
publication can still be recovered by inspection alone.
"""

import os
import shutil
import stat
import sys
from pathlib import Path


def _path_lstat(path: Path) -> os.stat_result | None:
    """Inspect one path without following a symbolic link."""
    try:
        return os.lstat(path)
    except FileNotFoundError:
        return None


def _is_windows_reparse_point(path_stat: os.stat_result) -> bool:
    """Return whether a Windows path is any reparse-point alias."""
    return sys.platform == "win32" and bool(
        getattr(path_stat, "st_file_attributes", 0)
        & getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0)
    )


def _require_export_directory_or_absent(path: Path, *, label: str) -> bool:
    """Require one reserved generation path to be absent or a real directory."""
    path_stat = _path_lstat(path)
    if path_stat is None:
        return False
    if (
        stat.S_ISLNK(path_stat.st_mode)
        or _is_windows_reparse_point(path_stat)
        or not stat.S_ISDIR(path_stat.st_mode)
    ):
        raise RuntimeError(f"{label} must be an ordinary non-symbolic directory or absent: {path}")
    return True


def _require_export_regular_file(path: Path, *, label: str) -> os.stat_result:
    """Require one reserved control path to be one unlinked regular file."""
    path_stat = _path_lstat(path)
    if path_stat is None:
        raise RuntimeError(f"{label} is missing: {path}")
    if (
        stat.S_ISLNK(path_stat.st_mode)
        or _is_windows_reparse_point(path_stat)
        or not stat.S_ISREG(path_stat.st_mode)
        or path_stat.st_nlink != 1
    ):
        raise RuntimeError(f"{label} must be one non-linked, non-symbolic regular file: {path}")
    return path_stat


def _fsync_export_directory(path: Path) -> None:
    """Durably order sibling creation, rename, and removal where supported."""
    if sys.platform == "win32":
        return
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_DIRECTORY", 0)
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _rename_export_path(source: Path, destination: Path) -> None:
    """Rename one transaction path and durably publish the directory entry."""
    source.rename(destination)
    _fsync_export_directory(destination.parent)


def _remove_export_tree(path: Path) -> None:
    """Remove one validated transaction-owned directory durably."""
    if not _require_export_directory_or_absent(
        path,
        label="Export transaction directory",
    ):
        return
    shutil.rmtree(path)
    _fsync_export_directory(path.parent)


def _remove_export_control_file(path: Path) -> None:
    """Remove one validated transaction-owned control file durably."""
    _require_export_regular_file(path, label="Export transaction control file")
    path.unlink()
    _fsync_export_directory(path.parent)
