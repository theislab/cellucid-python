from __future__ import annotations

import hashlib
import io
import tarfile
import tomllib
from pathlib import Path

import pytest
from scripts.normalize_sdist import MAX_NORMALIZED_SDIST_BYTES, normalize_sdist
from scripts.validate_release import (
    CURRENT_RUNTIME_REQUIREMENTS,
    EXCLUDED_SDIST_DIRECTORIES,
    REPRODUCIBLE_BUILD_EPOCH,
    _validate_recipe_text,
    validate_recipe,
    validate_release,
    validate_sdist,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def _write_source_distribution(
    path: Path,
    *,
    timestamp: int,
    extra_members: tuple[str, ...] = (),
    payload: bytes = b"version = '0.9.1'\n",
) -> None:
    with tarfile.open(path, mode="w:gz") as archive:
        directory = tarfile.TarInfo("cellucid-0.9.1")
        directory.type = tarfile.DIRTYPE
        directory.mode = 0o775
        directory.mtime = timestamp
        directory.uid = 501
        directory.gid = 20
        directory.uname = "builder"
        directory.gname = "staff"
        archive.addfile(directory)

        project = tarfile.TarInfo("cellucid-0.9.1/pyproject.toml")
        project.size = len(payload)
        project.mode = 0o664
        project.mtime = timestamp
        project.uid = 501
        project.gid = 20
        project.uname = "builder"
        project.gname = "staff"
        archive.addfile(project, io.BytesIO(payload))

        for extra_member in extra_members:
            unsafe = tarfile.TarInfo(extra_member)
            unsafe.size = 1
            archive.addfile(unsafe, io.BytesIO(b"x"))


def _digest_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_manifest_prunes_the_tests_directory() -> None:
    """The one line that holds the tests out of the sdist must stay present.

    ``scripts`` is excluded for an unrelated reason -- it is not a package --
    so only ``tests`` has a rule here to lose.
    """
    manifest = (REPOSITORY_ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "prune tests" in manifest.split("\n")


@pytest.mark.parametrize("excluded_directory", EXCLUDED_SDIST_DIRECTORIES)
def test_source_distribution_may_not_ship_excluded_directories(
    tmp_path: Path, excluded_directory: str
) -> None:
    """A matching digest must not be enough to pass an archive that ships them.

    Deleting ``MANIFEST.in`` republishes ``tests/`` to PyPI, and before this
    check the release gate still passed once the digest was re-recorded.
    """
    archive_path = tmp_path / "cellucid-0.9.1.tar.gz"
    _write_source_distribution(
        archive_path,
        timestamp=1_700_000_000,
        extra_members=(f"cellucid-0.9.1/{excluded_directory}/test_smuggled.py",),
    )

    with pytest.raises(SystemExit) as excluded_failure:
        validate_sdist(archive_path, "0.9.1", _digest_of(archive_path))

    assert excluded_directory in str(excluded_failure.value)
    assert "test_smuggled.py" in str(excluded_failure.value)


def test_source_distribution_without_excluded_directories_passes(tmp_path: Path) -> None:
    archive_path = tmp_path / "cellucid-0.9.1.tar.gz"
    _write_source_distribution(archive_path, timestamp=1_700_000_000)

    assert validate_sdist(archive_path, "0.9.1", _digest_of(archive_path)) == archive_path


def test_source_distribution_excludes_only_whole_directory_names(tmp_path: Path) -> None:
    """A file whose name merely starts with an excluded word is not excluded."""
    archive_path = tmp_path / "cellucid-0.9.1.tar.gz"
    _write_source_distribution(
        archive_path,
        timestamp=1_700_000_000,
        extra_members=("cellucid-0.9.1/src/cellucid/testing_helpers.py",),
    )

    assert validate_sdist(archive_path, "0.9.1", _digest_of(archive_path)) == archive_path


def test_source_distribution_normalization_is_byte_reproducible(tmp_path: Path) -> None:
    first = tmp_path / "cellucid-first.tar.gz"
    second = tmp_path / "cellucid-second.tar.gz"
    _write_source_distribution(first, timestamp=1_700_000_000)
    _write_source_distribution(second, timestamp=1_800_000_000)

    normalize_sdist(first)
    normalize_sdist(second)

    first_bytes = first.read_bytes()
    assert first_bytes == second.read_bytes()
    assert len(first_bytes) <= MAX_NORMALIZED_SDIST_BYTES
    assert int.from_bytes(first_bytes[4:8], byteorder="little") == 0
    assert (
        hashlib.sha256(first_bytes).hexdigest() == hashlib.sha256(second.read_bytes()).hexdigest()
    )

    with tarfile.open(first, mode="r:gz") as archive:
        members = archive.getmembers()
    assert [member.name for member in members] == sorted(member.name for member in members)
    assert all(member.mtime == 0 for member in members)
    assert all(member.uid == member.gid == 0 for member in members)
    assert all(member.uname == member.gname == "" for member in members)
    assert members[0].mode == 0o755
    assert members[1].mode == 0o644


@pytest.mark.parametrize(
    "member_name",
    [
        "../outside.txt",
        "cellucid-0.9.1/./alias.txt",
        "cellucid-0.9.1/back\\slash.txt",
        "C:/drive.txt",
        "cellucid-0.9.1/control\u0001.txt",
        "cellucid-0.9.1/CON.txt",
        "cellucid-0.9.1/trailing.",
        "cellucid-0.9.1/e\u0301.txt",
    ],
)
def test_source_distribution_normalization_rejects_noncanonical_members_atomically(
    tmp_path: Path,
    member_name: str,
) -> None:
    path = tmp_path / "cellucid-unsafe.tar.gz"
    _write_source_distribution(
        path,
        timestamp=1_700_000_000,
        extra_members=(member_name,),
    )
    original = path.read_bytes()

    with pytest.raises(ValueError, match="source-distribution member"):
        normalize_sdist(path)

    assert path.read_bytes() == original
    assert not (tmp_path / "outside.txt").exists()


def test_source_distribution_normalization_rejects_case_aliases_atomically(
    tmp_path: Path,
) -> None:
    path = tmp_path / "cellucid-alias.tar.gz"
    _write_source_distribution(
        path,
        timestamp=1_700_000_000,
        extra_members=(
            "cellucid-0.9.1/Case.txt",
            "cellucid-0.9.1/case.txt",
        ),
    )
    original = path.read_bytes()

    with pytest.raises(ValueError, match="aliased source-distribution member"):
        normalize_sdist(path)

    assert path.read_bytes() == original


def test_source_distribution_normalization_enforces_the_release_size_limit(
    tmp_path: Path,
) -> None:
    path = tmp_path / "cellucid-large.tar.gz"
    _write_source_distribution(
        path,
        timestamp=1_700_000_000,
        payload=b"x" * (MAX_NORMALIZED_SDIST_BYTES + 1),
    )
    original = path.read_bytes()

    with pytest.raises(ValueError, match="exceeds the release limit"):
        normalize_sdist(path)

    assert path.read_bytes() == original


def test_checked_in_downstream_recipe_is_one_exact_current_contract() -> None:
    version = validate_release(tag=None)
    digest = validate_recipe(version)

    assert version == "0.9.1"
    assert len(digest) == 64
    assert digest == digest.lower()

    project = tomllib.loads((REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8"))[
        "project"
    ]
    assert tuple(project["dependencies"]) == CURRENT_RUNTIME_REQUIREMENTS
    assert len(project["dependencies"]) == len(set(project["dependencies"]))


def test_ci_builds_use_one_canonical_zip_compatible_epoch() -> None:
    assert REPRODUCIBLE_BUILD_EPOCH == "315532800"
    for relative_path in (
        ".github/workflows/test.yml",
        ".github/workflows/pypi-publish.yml",
    ):
        workflow = (REPOSITORY_ROOT / relative_path).read_text(encoding="utf-8")
        assert workflow.count(f'SOURCE_DATE_EPOCH: "{REPRODUCIBLE_BUILD_EPOCH}"') == 1


@pytest.mark.parametrize(
    "mutation",
    [
        lambda recipe: recipe + '\n{% set version = "0.9.1" %}\n',
        lambda recipe: recipe.replace(
            "  sha256: ",
            "  sha256: "
            "0000000000000000000000000000000000000000000000000000000000000000\n"
            "  sha256: ",
            1,
        ),
        lambda recipe: recipe.replace(
            "  url: ",
            "  url: https://example.invalid/cellucid.tar.gz\n  url: ",
            1,
        ),
        lambda recipe: recipe.replace(
            "  run:\n",
            "  run:\n    - python >=3.11,<3.15\n",
            1,
        ),
        lambda recipe: recipe.replace(
            "    - numpy >=1.26,<3\n    - pandas",
            "    - pandas >=2.1.0,!=2.1.2,<3\n    - numpy >=1.26,<3\n    - pandas",
            1,
        ),
    ],
)
def test_downstream_recipe_rejects_duplicate_or_reordered_identity(
    mutation,
) -> None:
    recipe = (REPOSITORY_ROOT / "scripts/publishing/meta.yaml").read_text(encoding="utf-8")

    with pytest.raises(SystemExit):
        _validate_recipe_text(mutation(recipe), "0.9.1")
