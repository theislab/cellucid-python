#!/usr/bin/env python3
"""Normalize one source distribution into a reproducible tar.gz archive."""

from __future__ import annotations

import argparse
import binascii
import copy
import io
import os
import struct
import tarfile
import tempfile
import unicodedata
from pathlib import Path, PurePosixPath

CANONICAL_FILE_MODE = 0o644
CANONICAL_EXECUTABLE_MODE = 0o755
CANONICAL_DIRECTORY_MODE = 0o755
DEFLATE_STORED_BLOCK_BYTES = 65_535
MAX_NORMALIZED_SDIST_BYTES = 1_048_576
MAX_SDIST_MEMBERS = 4_096
WINDOWS_FORBIDDEN_CHARACTERS = frozenset('<>:"\\|?*')
WINDOWS_RESERVED_NAMES = frozenset(
    {"CON", "PRN", "AUX", "NUL"}
    | {f"COM{index}" for index in range(1, 10)}
    | {f"LPT{index}" for index in range(1, 10)}
)


def _resolve_sdist(path: Path) -> Path:
    if path.is_dir():
        candidates = sorted(path.glob("cellucid-*.tar.gz"))
        if len(candidates) != 1:
            raise ValueError(f"{path} must contain exactly one cellucid source distribution")
        path = candidates[0]
    if not path.is_file() or not path.name.startswith("cellucid-"):
        raise ValueError(f"{path} is not a Cellucid source distribution")
    if not path.name.endswith(".tar.gz"):
        raise ValueError(f"{path} must end in .tar.gz")
    return path


def _validated_members(archive: tarfile.TarFile) -> list[tarfile.TarInfo]:
    members = archive.getmembers()
    if not members:
        raise ValueError("source distribution is empty")
    if len(members) > MAX_SDIST_MEMBERS:
        raise ValueError("source distribution has too many members")

    names: set[str] = set()
    logical_names: set[str] = set()
    roots: set[str] = set()
    payload_bytes = 0
    for member in members:
        name = PurePosixPath(member.name)
        parts = member.name.split("/")
        if (
            name.is_absolute()
            or not name.parts
            or name.as_posix() != member.name
            or any(not part or part in {".", ".."} for part in parts)
        ):
            raise ValueError(f"unsafe source-distribution member: {member.name!r}")
        for part in parts:
            if (
                unicodedata.normalize("NFC", part) != part
                or part.endswith((" ", "."))
                or any(character in WINDOWS_FORBIDDEN_CHARACTERS for character in part)
                or any(ord(character) < 32 or ord(character) == 127 for character in part)
                or len(part.encode("utf-8")) > 255
                or part.split(".", 1)[0].upper() in WINDOWS_RESERVED_NAMES
            ):
                raise ValueError(f"nonportable source-distribution member: {member.name!r}")
        if member.name in names:
            raise ValueError(f"duplicate source-distribution member: {member.name!r}")
        logical_name = member.name.casefold()
        if logical_name in logical_names:
            raise ValueError(f"aliased source-distribution member: {member.name!r}")
        if not (member.isdir() or member.isfile()):
            raise ValueError(f"unsupported source-distribution member: {member.name!r}")
        names.add(member.name)
        logical_names.add(logical_name)
        roots.add(name.parts[0])
        payload_bytes += member.size
        if payload_bytes > MAX_NORMALIZED_SDIST_BYTES:
            raise ValueError("source-distribution payload exceeds the release limit")

    if len(roots) != 1:
        raise ValueError("source distribution must contain exactly one root directory")
    root = next(iter(roots))
    root_members = [member for member in members if member.name == root]
    if len(root_members) != 1 or not root_members[0].isdir():
        raise ValueError("source distribution root must be one declared directory")
    return sorted(members, key=lambda member: member.name)


def _canonical_member(member: tarfile.TarInfo) -> tarfile.TarInfo:
    normalized = copy.copy(member)
    normalized.mtime = 0
    normalized.uid = 0
    normalized.gid = 0
    normalized.uname = ""
    normalized.gname = ""
    normalized.pax_headers = {}
    normalized.devmajor = 0
    normalized.devminor = 0
    if normalized.isdir():
        normalized.mode = CANONICAL_DIRECTORY_MODE
    elif member.mode & 0o111:
        normalized.mode = CANONICAL_EXECUTABLE_MODE
    else:
        normalized.mode = CANONICAL_FILE_MODE
    return normalized


def _write_deterministic_gzip(target, payload: bytes) -> None:
    target.write(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff")
    if not payload:
        target.write(b"\x01\x00\x00\xff\xff")
    else:
        offset = 0
        while offset < len(payload):
            block = payload[offset : offset + DEFLATE_STORED_BLOCK_BYTES]
            offset += len(block)
            target.write(b"\x01" if offset == len(payload) else b"\x00")
            target.write(struct.pack("<HH", len(block), len(block) ^ 0xFFFF))
            target.write(block)
    target.write(
        struct.pack(
            "<II",
            binascii.crc32(payload) & 0xFFFFFFFF,
            len(payload) & 0xFFFFFFFF,
        )
    )


def normalize_sdist(path: Path) -> Path:
    """Rewrite *path* atomically with canonical member order and metadata."""
    source_path = _resolve_sdist(path.resolve())
    source_mode = source_path.stat().st_mode & 0o777
    temporary_path: Path | None = None
    try:
        with tarfile.open(source_path, mode="r:gz") as source_archive:
            members = _validated_members(source_archive)
            tar_bytes = io.BytesIO()
            with tarfile.open(
                fileobj=tar_bytes,
                mode="w",
                format=tarfile.PAX_FORMAT,
            ) as target_archive:
                for member in members:
                    canonical_member = _canonical_member(member)
                    payload = source_archive.extractfile(member) if member.isfile() else None
                    target_archive.addfile(canonical_member, payload)
            with tempfile.NamedTemporaryFile(
                dir=source_path.parent,
                prefix=f".{source_path.name}.",
                suffix=".tmp",
                delete=False,
            ) as temporary_file:
                temporary_path = Path(temporary_file.name)
                _write_deterministic_gzip(temporary_file, tar_bytes.getvalue())
                temporary_file.flush()
                os.fsync(temporary_file.fileno())
                if temporary_file.tell() > MAX_NORMALIZED_SDIST_BYTES:
                    raise ValueError("normalized source distribution exceeds the release limit")
        os.chmod(temporary_path, source_mode)
        os.replace(temporary_path, source_path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return source_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Normalize exactly one Cellucid source distribution."
    )
    parser.add_argument(
        "path",
        type=Path,
        help="A cellucid-*.tar.gz file or a directory containing exactly one.",
    )
    normalized = normalize_sdist(parser.parse_args().path)
    print(f"Normalized source distribution: {normalized}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
