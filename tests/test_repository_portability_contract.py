from __future__ import annotations

import re
import unicodedata
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PORTABLE_ROOTS = (
    REPOSITORY_ROOT / ".github",
    REPOSITORY_ROOT / "docs",
    REPOSITORY_ROOT / "examples",
    REPOSITORY_ROOT / "scripts",
    REPOSITORY_ROOT / "src",
    REPOSITORY_ROOT / "tests",
)
WINDOWS_RESERVED_NAME = re.compile(
    r"^(?:CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\.|$)",
    re.IGNORECASE,
)
WINDOWS_INVALID_CHARACTERS = set('<>:"|?*')


def _repository_entries() -> list[Path]:
    entries: list[Path] = []
    for root in PORTABLE_ROOTS:
        entries.append(root)
        entries.extend(root.rglob("*"))
    entries.extend(
        REPOSITORY_ROOT / filename
        for filename in (
            "CHANGELOG.md",
            "README.md",
            "pyproject.toml",
        )
    )
    return entries


def test_python_repository_paths_are_one_portable_regular_file_graph() -> None:
    entries = _repository_entries()
    symlinks = [
        str(path.relative_to(REPOSITORY_ROOT))
        for path in entries
        if path.is_symlink()
    ]
    assert symlinks == []

    folded_paths: dict[str, str] = {}
    collisions: list[str] = []
    invalid_paths: list[str] = []
    for path in entries:
        relative = path.relative_to(REPOSITORY_ROOT).as_posix()
        folded = unicodedata.normalize("NFC", relative).casefold()
        previous = folded_paths.get(folded)
        if previous is not None and previous != relative:
            collisions.append(f"{previous} <> {relative}")
        folded_paths[folded] = relative

        for segment in path.relative_to(REPOSITORY_ROOT).parts:
            if (
                segment.endswith((".", " "))
                or WINDOWS_RESERVED_NAME.match(segment)
                or any(character in WINDOWS_INVALID_CHARACTERS for character in segment)
            ):
                invalid_paths.append(relative)
                break

    assert collisions == []
    assert invalid_paths == []
