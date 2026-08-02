"""
The exclusive, non-waiting writer lock for one export target.

Two writers publishing into one directory would interleave renames, so a
generation is held under an advisory byte-range lock on a sidecar file, taken
with the same primitive cellucid-r uses so the two exporters exclude each
other. The lock file identity is captured and re-validated after acquisition,
because an unlink/create race can hand out the same inode twice, and a
process-local registry plus fork handlers stop one process from mistaking its
own inherited descriptor for a free lock.
"""

import errno
import os
import stat
import sys
import threading
from collections.abc import Iterator
from contextlib import contextmanager, suppress
from pathlib import Path
from typing import NamedTuple

from ._filesystem import _is_windows_reparse_point

if sys.platform == "win32":
    import msvcrt
else:
    import fcntl

_EXPORT_LOCK_REGISTRY_GUARD = threading.RLock()
_EXPORT_LOCK_REGISTRY: dict[str, tuple[Path, int]] = {}
_EXPORT_LOCK_REGISTRY_PID = os.getpid()


class _ExportLockGeneration(NamedTuple):
    """Stable metadata used to distinguish reused lock-file identities."""

    device: int
    inode: int
    mode: int
    link_count: int
    owner: int
    group: int
    size: int
    modified_ns: int
    changed_ns: int
    birth_ns: int | None
    filesystem_generation: int
    flags: int
    file_attributes: int
    reparse_tag: int


def _acquire_export_lock_registry_before_fork() -> None:
    """Prevent a fork from inheriting an unregistered coordination descriptor."""
    _EXPORT_LOCK_REGISTRY_GUARD.acquire()


def _release_export_lock_registry_after_fork_in_parent() -> None:
    """Release the registry snapshot held across a parent-side fork."""
    _EXPORT_LOCK_REGISTRY_GUARD.release()


def _reset_export_lock_registry_after_fork() -> None:
    """Close inherited descriptors and drop record-lock claims after a fork."""
    global _EXPORT_LOCK_REGISTRY_GUARD
    global _EXPORT_LOCK_REGISTRY
    global _EXPORT_LOCK_REGISTRY_PID

    inherited_descriptors = tuple(
        descriptor for _lock_path, descriptor in _EXPORT_LOCK_REGISTRY.values()
    )
    _EXPORT_LOCK_REGISTRY_GUARD = threading.RLock()
    _EXPORT_LOCK_REGISTRY = {}
    _EXPORT_LOCK_REGISTRY_PID = os.getpid()
    for descriptor in inherited_descriptors:
        with suppress(OSError):
            os.close(descriptor)


_register_at_fork = getattr(os, "register_at_fork", None)
if _register_at_fork is not None:
    _register_at_fork(
        before=_acquire_export_lock_registry_before_fork,
        after_in_parent=_release_export_lock_registry_after_fork_in_parent,
        after_in_child=_reset_export_lock_registry_after_fork,
    )


def _export_lock_generation(path_stat: os.stat_result) -> _ExportLockGeneration:
    """Capture identity plus mutation metadata for one inspected lock file.

    Device and inode alone are insufficient: filesystems may immediately reuse
    an inode after an unlink/create race. Change, birth, and filesystem
    generation metadata distinguish that replacement where the platform
    exposes them, while size and modification metadata also detect same-inode
    mutation between inspection and acquisition.
    """
    birth_ns_value = getattr(path_stat, "st_birthtime_ns", None)
    if birth_ns_value is None:
        birth_seconds = getattr(path_stat, "st_birthtime", None)
        if birth_seconds is not None:
            birth_ns_value = round(float(birth_seconds) * 1_000_000_000)
    return _ExportLockGeneration(
        device=int(path_stat.st_dev),
        inode=int(path_stat.st_ino),
        mode=int(path_stat.st_mode),
        link_count=int(path_stat.st_nlink),
        owner=int(getattr(path_stat, "st_uid", 0)),
        group=int(getattr(path_stat, "st_gid", 0)),
        size=int(path_stat.st_size),
        modified_ns=int(path_stat.st_mtime_ns),
        changed_ns=int(path_stat.st_ctime_ns),
        birth_ns=int(birth_ns_value) if birth_ns_value is not None else None,
        filesystem_generation=int(getattr(path_stat, "st_gen", 0)),
        flags=int(getattr(path_stat, "st_flags", 0)),
        file_attributes=int(getattr(path_stat, "st_file_attributes", 0)),
        reparse_tag=int(getattr(path_stat, "st_reparse_tag", 0)),
    )


def _same_export_lock_identity(
    descriptor_stat: os.stat_result,
    path_stat: os.stat_result,
) -> bool:
    """Compare the stable filesystem identity exposed by both stat APIs.

    Windows does not guarantee that every metadata field returned for an open
    descriptor is synchronized with a fresh path lookup.  In particular,
    timestamp fields have changed semantics across supported Python releases.
    Device and file ID are the cross-interface identity; generation metadata
    is compared only between path lookups made through the same API.
    """
    return (
        int(descriptor_stat.st_dev),
        int(descriptor_stat.st_ino),
    ) == (
        int(path_stat.st_dev),
        int(path_stat.st_ino),
    )


def _canonical_export_lock_key(lock_path: Path) -> str:
    """Return one process-local identity for path aliases to a lock inode."""
    return os.path.normcase(os.path.realpath(os.fspath(lock_path)))


def _open_and_reserve_process_export_lock(
    lock_path: Path,
) -> tuple[str, int, _ExportLockGeneration] | None:
    """Open a lock path only after excluding same-process inode aliases."""
    global _EXPORT_LOCK_REGISTRY_PID

    process_id = os.getpid()
    with _EXPORT_LOCK_REGISTRY_GUARD:
        if process_id != _EXPORT_LOCK_REGISTRY_PID:
            _EXPORT_LOCK_REGISTRY.clear()
            _EXPORT_LOCK_REGISTRY_PID = process_id

        lock_key = _canonical_export_lock_key(lock_path)
        if lock_path.is_symlink():
            raise RuntimeError(f"Export lock path must not be a symbolic link: {lock_path}")
        try:
            path_stat = os.lstat(lock_path)
        except FileNotFoundError:
            path_stat = None
        if path_stat is not None and (
            not stat.S_ISREG(path_stat.st_mode)
            or path_stat.st_nlink != 1
            or _is_windows_reparse_point(path_stat)
        ):
            raise RuntimeError(
                "Export lock path must identify one non-linked regular non-symbolic "
                f"file: {lock_path}"
            )
        if lock_key in _EXPORT_LOCK_REGISTRY:
            return None
        for claimed_path, _descriptor in _EXPORT_LOCK_REGISTRY.values():
            try:
                if os.path.samefile(lock_path, claimed_path):
                    return None
            except FileNotFoundError:
                pass

        create_lock_file = path_stat is None
        try:
            try:
                descriptor, generation = _open_export_lock_descriptor(
                    lock_path,
                    create=create_lock_file,
                    expected_stat=path_stat,
                )
            except FileExistsError:
                if not create_lock_file:
                    raise
                raced_path_stat = os.lstat(lock_path)
                descriptor, generation = _open_export_lock_descriptor(
                    lock_path,
                    create=False,
                    expected_stat=raced_path_stat,
                )
        except FileNotFoundError as error:
            raise RuntimeError(
                f"Export lock path changed while establishing ownership: {lock_path}"
            ) from error
        _EXPORT_LOCK_REGISTRY[lock_key] = (lock_path, descriptor)
        return lock_key, descriptor, generation


def _release_process_export_lock(
    lock_key: str,
    descriptor: int | None,
) -> None:
    """Release one matching process-local target claim."""
    with _EXPORT_LOCK_REGISTRY_GUARD:
        active_claim = _EXPORT_LOCK_REGISTRY.get(lock_key)
        if active_claim is not None and active_claim[1] == descriptor:
            _EXPORT_LOCK_REGISTRY.pop(lock_key, None)


def _validate_export_lock_descriptor(
    lock_path: Path,
    descriptor: int,
    *,
    expected_generation: _ExportLockGeneration | None,
) -> _ExportLockGeneration:
    """Require one descriptor to retain the inspected lock-file generation."""
    descriptor_stat = os.fstat(descriptor)
    try:
        path_stat = os.lstat(lock_path)
    except FileNotFoundError as error:
        raise RuntimeError(
            f"Export lock path changed while establishing ownership: {lock_path}"
        ) from error
    is_windows_reparse_point = _is_windows_reparse_point(path_stat)
    path_generation = _export_lock_generation(path_stat)
    if (
        (expected_generation is not None and path_generation != expected_generation)
        or not _same_export_lock_identity(descriptor_stat, path_stat)
    ):
        raise RuntimeError(f"Export lock path changed while establishing ownership: {lock_path}")
    if (
        not stat.S_ISREG(descriptor_stat.st_mode)
        or descriptor_stat.st_nlink != 1
        or path_stat.st_nlink != 1
        or stat.S_ISLNK(path_stat.st_mode)
        or is_windows_reparse_point
    ):
        raise RuntimeError(
            f"Export lock path must identify one non-linked regular non-symbolic file: {lock_path}"
        )
    return path_generation


def _open_export_lock_descriptor(
    lock_path: Path,
    *,
    create: bool,
    expected_stat: os.stat_result | None,
) -> tuple[int, _ExportLockGeneration]:
    """Open one persistent, empty, non-symbolic regular coordination file."""
    if lock_path.is_symlink():
        raise RuntimeError(f"Export lock path must not be a symbolic link: {lock_path}")

    flags = os.O_RDWR
    if create:
        flags |= os.O_CREAT | os.O_EXCL
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(lock_path, flags, 0o600)
    try:
        generation = _validate_export_lock_descriptor(
            lock_path,
            descriptor,
            expected_generation=(
                _export_lock_generation(expected_stat) if expected_stat is not None else None
            ),
        )
    except BaseException as error:
        try:
            os.close(descriptor)
        except OSError as close_error:
            error.add_note(
                "The invalid export lock descriptor also failed to close: "
                f"{type(close_error).__name__}: {close_error}"
            )
        raise
    return descriptor, generation


if sys.platform == "win32":

    def _acquire_export_lock_descriptor(descriptor: int) -> None:
        """Acquire the R-interoperable Windows byte range without waiting."""
        os.lseek(descriptor, 0, os.SEEK_SET)
        msvcrt.locking(descriptor, msvcrt.LK_NBLCK, 1)

    def _release_export_lock_descriptor(descriptor: int) -> None:
        """Release the exact Windows byte range acquired above."""
        os.lseek(descriptor, 0, os.SEEK_SET)
        msvcrt.locking(descriptor, msvcrt.LK_UNLCK, 1)

else:

    def _acquire_export_lock_descriptor(descriptor: int) -> None:
        """Acquire the R-interoperable POSIX byte range without waiting."""
        fcntl.lockf(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB, 0, 0, os.SEEK_SET)


    def _release_export_lock_descriptor(descriptor: int) -> None:
        """Release the exact POSIX byte range acquired above."""
        fcntl.lockf(descriptor, fcntl.LOCK_UN, 0, 0, os.SEEK_SET)


def _is_export_lock_contention(error: OSError) -> bool:
    """Distinguish a live owner from filesystem or permission failures."""
    contention_errors = {errno.EACCES, errno.EAGAIN}
    if sys.platform == "win32":
        contention_errors.add(getattr(errno, "EDEADLOCK", errno.EDEADLK))
    return error.errno in contention_errors


@contextmanager
def _exclusive_export_generation(target_dir: Path) -> Iterator[None]:
    """Hold one non-waiting, process-owned writer lock for an export target."""
    lock_path = target_dir.parent / f".{target_dir.name}.cellucid.lock"
    reservation = _open_and_reserve_process_export_lock(lock_path)
    if reservation is None:
        raise RuntimeError(f"An export generation is already active for {target_dir}.")

    lock_key, descriptor, opened_generation = reservation
    owner_pid = os.getpid()
    descriptor_locked = False
    try:
        try:
            _acquire_export_lock_descriptor(descriptor)
        except OSError as error:
            if _is_export_lock_contention(error):
                raise RuntimeError(
                    f"An export generation is already active for {target_dir}."
                ) from error
            raise
        descriptor_locked = True
        _validate_export_lock_descriptor(
            lock_path,
            descriptor,
            expected_generation=opened_generation,
        )
        yield
    finally:
        if os.getpid() == owner_pid:
            active_error = sys.exc_info()[1]
            cleanup_errors: list[OSError] = []
            cleanup_cancellation: BaseException | None = None
            os_lock_released = not descriptor_locked
            try:
                if descriptor is not None and descriptor_locked:
                    try:
                        _release_export_lock_descriptor(descriptor)
                        os_lock_released = True
                    except OSError as unlock_error:
                        cleanup_errors.append(unlock_error)
            except BaseException as unlock_cancellation:
                cleanup_cancellation = unlock_cancellation
            finally:
                if descriptor is not None:
                    try:
                        os.close(descriptor)
                        os_lock_released = True
                    except OSError as close_error:
                        cleanup_errors.append(close_error)
                    except BaseException as close_cancellation:
                        if cleanup_cancellation is None:
                            cleanup_cancellation = close_cancellation
                        else:
                            cleanup_cancellation.add_note(
                                "Closing the export lock descriptor also failed: "
                                f"{type(close_cancellation).__name__}: "
                                f"{close_cancellation}"
                            )

            if os_lock_released:
                _release_process_export_lock(
                    lock_key,
                    descriptor,
                )

            if cleanup_cancellation is not None:
                for cleanup_error in cleanup_errors:
                    cleanup_cancellation.add_note(
                        "Export lock cleanup also failed: "
                        f"{type(cleanup_error).__name__}: {cleanup_error}"
                    )
                raise cleanup_cancellation

            if cleanup_errors:
                cleanup_message = (
                    f"Failed to release the export generation lock for {target_dir}: "
                    + "; ".join(f"{type(error).__name__}: {error}" for error in cleanup_errors)
                )
                if active_error is not None:
                    active_error.add_note(cleanup_message)
                else:
                    raise RuntimeError(cleanup_message) from cleanup_errors[0]
