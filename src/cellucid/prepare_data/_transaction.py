"""
The write-ahead export transaction.

Publishing a generation is a rename, not a copy: the staged directory takes the
target's place and the previous generation is renamed aside and removed. A
journal recording the transaction identity, and whether a target existed, is
written and fsynced before any of that starts, so an interrupted run leaves a
state that can be read back and resolved without guessing whether to commit or
roll back. The journal bytes are exact and canonical because cellucid-r writes
and reads the same file.
"""

import json
import os
import re
import secrets
from pathlib import Path

from ._filesystem import (
    _fsync_export_directory,
    _path_lstat,
    _remove_export_control_file,
    _remove_export_tree,
    _rename_export_path,
    _require_export_directory_or_absent,
    _require_export_regular_file,
)

_EXPORT_TRANSACTION_FORMAT = "cellucid-export-transaction"
_EXPORT_TRANSACTION_VERSION = 1
_EXPORT_TRANSACTION_ID_PATTERN = re.compile(r"^[0-9a-f]{32}$", flags=re.ASCII)


def _export_transaction_paths(
    target_dir: Path,
    transaction_id: str,
) -> tuple[Path, Path, Path, Path]:
    """Derive every reserved path from one validated transaction identity."""
    if not _EXPORT_TRANSACTION_ID_PATTERN.fullmatch(transaction_id):
        raise ValueError("Export transaction identity is not canonical.")
    parent = target_dir.parent
    stem = f".{target_dir.name}"
    return (
        parent / f"{stem}.cellucid-transaction.json",
        parent / f"{stem}.cellucid-transaction.json.tmp",
        parent / f"{stem}.cellucid-stage-{transaction_id}",
        parent / f"{stem}.cellucid-backup-{transaction_id}",
    )


def _export_transaction_control_paths(target_dir: Path) -> tuple[Path, Path]:
    """Return the fixed journal and its fixed atomic-write temporary path."""
    journal, journal_temp, _stage, _backup = _export_transaction_paths(
        target_dir,
        "0" * 32,
    )
    return journal, journal_temp


def _serialize_export_transaction(
    transaction_id: str,
    *,
    had_target: bool,
) -> bytes:
    """Serialize the exact Python/R interoperable transaction descriptor."""
    if not _EXPORT_TRANSACTION_ID_PATTERN.fullmatch(transaction_id):
        raise ValueError("Export transaction identity is not canonical.")
    if type(had_target) is not bool:
        raise TypeError("Export transaction had_target must be exactly boolean.")
    had_target_json = "true" if had_target else "false"
    return (
        f'{{"format":"{_EXPORT_TRANSACTION_FORMAT}",'
        f'"version":{_EXPORT_TRANSACTION_VERSION},'
        f'"transaction_id":"{transaction_id}",'
        f'"had_target":{had_target_json}}}\n'
    ).encode("ascii")


def _reject_duplicate_json_members(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    """Reject ambiguous transaction journals instead of keeping a last value."""
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate export transaction journal member: {key!r}.")
        result[key] = value
    return result


def _read_export_transaction(journal_path: Path) -> tuple[str, bool]:
    """Read and strictly validate one cross-language transaction journal."""
    path_stat = _require_export_regular_file(
        journal_path,
        label="Export transaction journal",
    )
    if path_stat.st_size > 512:
        raise RuntimeError(f"Export transaction journal is unexpectedly large: {journal_path}")
    try:
        journal_bytes = journal_path.read_bytes()
        payload = json.loads(
            journal_bytes,
            object_pairs_hook=_reject_duplicate_json_members,
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise RuntimeError(f"Export transaction journal is malformed: {journal_path}") from error
    if type(payload) is not dict or set(payload) != {
        "format",
        "version",
        "transaction_id",
        "had_target",
    }:
        raise RuntimeError(f"Export transaction journal has an invalid schema: {journal_path}")
    transaction_id = payload["transaction_id"]
    had_target = payload["had_target"]
    if (
        payload["format"] != _EXPORT_TRANSACTION_FORMAT
        or type(payload["version"]) is not int
        or payload["version"] != _EXPORT_TRANSACTION_VERSION
        or type(transaction_id) is not str
        or not _EXPORT_TRANSACTION_ID_PATTERN.fullmatch(transaction_id)
        or type(had_target) is not bool
        or journal_bytes
        != _serialize_export_transaction(
            transaction_id,
            had_target=had_target,
        )
    ):
        raise RuntimeError(f"Export transaction journal is not canonical: {journal_path}")
    return transaction_id, had_target


def _require_active_export_transaction(
    journal_path: Path,
    transaction_id: str,
    *,
    had_target: bool,
) -> None:
    """Require the journal to retain the exact active transaction owner."""
    journal_transaction_id, journal_had_target = _read_export_transaction(
        journal_path,
    )
    if journal_transaction_id != transaction_id or journal_had_target is not had_target:
        raise RuntimeError("Export transaction journal does not describe the active transaction.")


def _discard_export_transaction_temp(journal_temp: Path) -> None:
    """Recover a process death during the journal's atomic file write."""
    if _path_lstat(journal_temp) is None:
        return
    _remove_export_control_file(journal_temp)


def _recover_export_transaction(target_dir: Path) -> None:
    """Resolve every valid interrupted publication state before a new write."""
    journal_path, journal_temp = _export_transaction_control_paths(target_dir)
    _discard_export_transaction_temp(journal_temp)
    if _path_lstat(journal_path) is None:
        return

    transaction_id, had_target = _read_export_transaction(journal_path)
    _journal, _journal_temp, staging_dir, backup_dir = _export_transaction_paths(
        target_dir,
        transaction_id,
    )
    target_exists = _require_export_directory_or_absent(
        target_dir,
        label="Export target",
    )
    stage_exists = _require_export_directory_or_absent(
        staging_dir,
        label="Staged export generation",
    )
    backup_exists = _require_export_directory_or_absent(
        backup_dir,
        label="Prior export generation",
    )
    state = (target_exists, stage_exists, backup_exists)

    if had_target:
        if state == (True, False, False):
            pass
        elif state == (True, True, False):
            _remove_export_tree(staging_dir)
        elif state == (False, True, True):
            _rename_export_path(backup_dir, target_dir)
            _remove_export_tree(staging_dir)
        elif state == (True, False, True):
            _remove_export_tree(backup_dir)
        else:
            raise RuntimeError(
                "Export transaction cannot be recovered without guessing whether "
                f"to commit or roll back: target/stage/backup state is {state}."
            )
    else:
        if state == (False, False, False):
            pass
        elif state == (False, True, False):
            _remove_export_tree(staging_dir)
        elif state == (True, False, False):
            pass
        else:
            raise RuntimeError(
                "Initial export transaction cannot be recovered without guessing "
                f"whether to commit or roll back: target/stage/backup state is {state}."
            )

    _require_active_export_transaction(
        journal_path,
        transaction_id,
        had_target=had_target,
    )
    _remove_export_control_file(journal_path)


def _new_export_transaction_id(target_dir: Path) -> str:
    """Choose one unpredictable identity whose derived artifact paths are free."""
    for _attempt in range(128):
        transaction_id = secrets.token_hex(16)
        _journal, _journal_temp, stage, backup = _export_transaction_paths(
            target_dir,
            transaction_id,
        )
        if _path_lstat(stage) is None and _path_lstat(backup) is None:
            return transaction_id
    raise RuntimeError("Could not allocate a unique export transaction identity.")


def _write_export_transaction(
    target_dir: Path,
    transaction_id: str,
    *,
    had_target: bool,
) -> None:
    """Atomically and durably publish one write-ahead transaction journal."""
    journal_path, journal_temp, _stage, _backup = _export_transaction_paths(
        target_dir,
        transaction_id,
    )
    if _path_lstat(journal_path) is not None or _path_lstat(journal_temp) is not None:
        raise RuntimeError(f"Export transaction control path is already occupied for {target_dir}.")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(journal_temp, flags, 0o600)
    journal_bytes = _serialize_export_transaction(
        transaction_id,
        had_target=had_target,
    )
    try:
        with os.fdopen(descriptor, "wb", closefd=True) as stream:
            descriptor = -1
            stream.write(journal_bytes)
            stream.flush()
            os.fsync(stream.fileno())
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    _rename_export_path(journal_temp, journal_path)


def _begin_export_transaction(
    target_dir: Path,
) -> tuple[str, bool, Path, Path]:
    """Journal ownership before creating any staged generation directory."""
    transaction_id = _new_export_transaction_id(target_dir)
    had_target = _require_export_directory_or_absent(
        target_dir,
        label="Export target",
    )
    _write_export_transaction(
        target_dir,
        transaction_id,
        had_target=had_target,
    )
    _journal, _journal_temp, staging_dir, backup_dir = _export_transaction_paths(
        target_dir,
        transaction_id,
    )
    staging_dir.mkdir()
    _fsync_export_directory(target_dir.parent)
    return transaction_id, had_target, staging_dir, backup_dir


def _publish_export_generation(
    staging_dir: Path,
    target_dir: Path,
    *,
    transaction_id: str,
    had_target: bool,
) -> None:
    """Publish and settle one journal-owned complete export generation."""
    journal_path, journal_temp, expected_stage, backup_dir = _export_transaction_paths(
        target_dir,
        transaction_id,
    )
    if staging_dir != expected_stage:
        raise RuntimeError("Staged export path does not belong to the active transaction.")
    if _path_lstat(journal_temp) is not None:
        raise RuntimeError(
            f"Export transaction temporary control path reappeared before publication: "
            f"{journal_temp}"
        )
    _require_active_export_transaction(
        journal_path,
        transaction_id,
        had_target=had_target,
    )
    if not _require_export_directory_or_absent(
        staging_dir,
        label="Staged export generation",
    ):
        raise RuntimeError(f"Staged export directory is missing: {staging_dir}")
    if had_target:
        if not _require_export_directory_or_absent(
            target_dir,
            label="Export target",
        ):
            raise RuntimeError(f"Prior export generation is missing: {target_dir}")
        if _path_lstat(backup_dir) is not None:
            raise RuntimeError(f"Export backup path is already occupied: {backup_dir}")
        _rename_export_path(target_dir, backup_dir)
    elif _path_lstat(target_dir) is not None:
        raise RuntimeError(f"Initial export target appeared during publication: {target_dir}")

    _rename_export_path(staging_dir, target_dir)
    if had_target:
        _remove_export_tree(backup_dir)
    _require_active_export_transaction(
        journal_path,
        transaction_id,
        had_target=had_target,
    )
    _remove_export_control_file(journal_path)
