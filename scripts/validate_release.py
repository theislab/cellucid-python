#!/usr/bin/env python3
"""Validate the exact Cellucid Python release contract."""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import tarfile
import tomllib
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CURRENT_PYTHON_RANGE = ">=3.11,<3.15"
CURRENT_BUILD_REQUIREMENTS = ("setuptools==83.0.0",)
REPRODUCIBLE_BUILD_EPOCH = "315532800"
CURRENT_RUNTIME_REQUIREMENTS = (
    "numpy>=1.26,<3",
    "pandas>=2.1.0,!=2.1.2,<3",
    "scipy>=1.12,!=1.17.0,<2",
    "tqdm>=4.45",
    "anndata>=0.12.19,<0.13",
    "zarr>=3.1.4,<4",
    "numcodecs>=0.16.3,<0.17",
    "ipython>=7.23",
    "jupyter-server-proxy>=4.1",
)
CURRENT_CONDA_HOST_REQUIREMENTS = (
    "python >=3.11,<3.15",
    "setuptools ==83.0.0",
    "pip",
)
CURRENT_CONDA_RUNTIME_REQUIREMENTS = (
    "python >=3.11,<3.15",
    "numpy >=1.26,<3",
    "pandas >=2.1.0,!=2.1.2,<3",
    "scipy >=1.12,!=1.17.0,<2",
    "tqdm >=4.45",
    "anndata >=0.12.19,<0.13",
    "zarr >=3.1.4,<4",
    "numcodecs >=0.16.3,<0.17",
    "ipython >=7.23",
    "jupyter-server-proxy >=4.1",
)
HISTORICAL_RELEASE_TAGS = {
    "0.0.9": "v0.9.0",
    "0.0.1a0": "v0.0.1a0",
}


def _read(relative_path: str) -> str:
    return (REPOSITORY_ROOT / relative_path).read_text(encoding="utf-8")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def _recipe_requirements(recipe: str, section: str, next_section: str) -> tuple[str, ...]:
    matches = list(
        re.finditer(
            (
                rf"(?ms)^  {re.escape(section)}:\n"
                rf"(?P<body>.*?)(?=^(?:  )?{re.escape(next_section)}:)"
            ),
            recipe,
        )
    )
    _require(
        len(matches) == 1,
        f"meta.yaml must have exactly one {section} requirement section",
    )
    body_lines = matches[0].group("body").splitlines()
    _require(
        all(not line or line.startswith("    - ") for line in body_lines),
        f"meta.yaml {section} requirement section has an invalid entry",
    )
    return tuple(line.removeprefix("    - ") for line in body_lines if line.startswith("    - "))


def _validate_recipe_text(recipe: str, version: str) -> str:
    version_matches = re.findall(r'(?m)^\{% set version = "([^"]+)" %\}$', recipe)
    _require(
        version_matches == [version],
        f"meta.yaml must declare release {version} exactly once",
    )
    expected_url = (
        "https://pypi.io/packages/source/{{ name[0] }}/{{ name }}/cellucid-{{ version }}.tar.gz"
    )
    url_matches = re.findall(r"(?m)^  url: ([^\r\n]+)$", recipe)
    _require(
        url_matches == [expected_url],
        "meta.yaml must declare the exact versioned PyPI source URL once",
    )
    digest_matches = re.findall(r"(?m)^  sha256: ([0-9a-f]{64})$", recipe)
    _require(
        len(digest_matches) == 1,
        "meta.yaml must declare exactly one lowercase SHA-256 source digest",
    )
    _require(
        _recipe_requirements(recipe, "host", "run") == CURRENT_CONDA_HOST_REQUIREMENTS,
        "meta.yaml host requirements do not match the validated current contract",
    )
    _require(
        _recipe_requirements(recipe, "run", "test") == CURRENT_CONDA_RUNTIME_REQUIREMENTS,
        "meta.yaml runtime requirements do not match the validated current contract",
    )
    return digest_matches[0]


def validate_recipe(version: str) -> str:
    """Validate the downstream recipe and return its source-distribution digest."""
    return _validate_recipe_text(_read("scripts/publishing/meta.yaml"), version)


#: Top-level directories that must never reach the public source distribution.
#: ``tests`` is held out by the single ``prune tests`` line in ``MANIFEST.in``;
#: ``scripts`` is absent only because it is not a package and no rule adds it.
#: Neither exclusion is self-announcing, so both are asserted on the artifact.
EXCLUDED_SDIST_DIRECTORIES = ("tests", "scripts")


def _sdist_excluded_members(archive: tarfile.TarFile, version: str) -> tuple[str, ...]:
    prefixes = tuple(
        f"cellucid-{version}/{directory}/" for directory in EXCLUDED_SDIST_DIRECTORIES
    )
    return tuple(
        sorted(name for name in archive.getnames() if name.startswith(prefixes))
    )


def validate_sdist(path: Path, version: str, expected_digest: str) -> Path:
    """Require one normalized source distribution to match the recipe digest."""
    if path.is_dir():
        candidates = sorted(path.glob("cellucid-*.tar.gz"))
        _require(
            len(candidates) == 1,
            f"{path} must contain exactly one Cellucid source distribution",
        )
        path = candidates[0]
    _require(path.is_file(), f"source distribution does not exist: {path}")
    _require(
        path.name == f"cellucid-{version}.tar.gz",
        f"source distribution must be named cellucid-{version}.tar.gz",
    )
    actual_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    _require(
        actual_digest == expected_digest,
        (f"source-distribution SHA-256 {actual_digest} does not match meta.yaml {expected_digest}"),
    )
    with tarfile.open(path, mode="r:gz") as archive:
        excluded = _sdist_excluded_members(archive, version)
    _require(
        not excluded,
        (
            "source distribution must not ship "
            + " or ".join(EXCLUDED_SDIST_DIRECTORIES)
            + f": found {', '.join(excluded)}"
        ),
    )
    return path


def _validate_release_history(changelog: str, verify_local_tags: bool) -> None:
    for distribution_version, tag in HISTORICAL_RELEASE_TAGS.items():
        expected_link = (
            f"[{distribution_version}]: "
            f"https://github.com/theislab/cellucid-python/releases/tag/{tag}"
        )
        _require(
            changelog.count(expected_link) == 1,
            f"CHANGELOG.md must map {distribution_version} to {tag} exactly once",
        )
    if verify_local_tags:
        result = subprocess.run(
            ["git", "tag", "--list"],
            cwd=REPOSITORY_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        local_tags = set(result.stdout.splitlines())
        missing = sorted(set(HISTORICAL_RELEASE_TAGS.values()) - local_tags)
        _require(not missing, f"release-history tags are missing locally: {missing}")


def validate_release(tag: str | None, *, verify_local_tags: bool = False) -> str:
    """Validate package metadata and return the declared release version."""
    metadata = tomllib.loads(_read("pyproject.toml"))
    build_system = metadata["build-system"]
    _require(
        tuple(build_system["requires"]) == CURRENT_BUILD_REQUIREMENTS,
        "build-system requirements do not match the validated exact backend",
    )
    _require(
        build_system["build-backend"] == "setuptools.build_meta",
        "build backend must be exactly setuptools.build_meta",
    )
    project = metadata["project"]
    version = project["version"]
    _require(
        isinstance(version, str) and re.fullmatch(r"\d+\.\d+\.\d+", version) is not None,
        "pyproject.toml must declare one stable three-component version",
    )
    _require(
        project["requires-python"] == CURRENT_PYTHON_RANGE,
        f"requires-python must be exactly {CURRENT_PYTHON_RANGE}",
    )
    _require(
        tuple(project["dependencies"]) == CURRENT_RUNTIME_REQUIREMENTS,
        "runtime requirements do not match the validated current contract",
    )

    if tag is not None:
        _require(tag == f"v{version}", f"release tag {tag!r} must equal v{version}")

    changelog = _read("CHANGELOG.md")
    _require(f"## [{version}]" in changelog, "CHANGELOG.md has no current release heading")
    _require(
        "## [Unreleased]" not in changelog and "[Unreleased]:" not in changelog,
        "CHANGELOG.md must contain release entries only",
    )
    _require(
        f"[{version}]: https://github.com/theislab/cellucid-python/releases/tag/v{version}"
        in changelog,
        "CHANGELOG.md has no current release link",
    )
    _validate_release_history(changelog, verify_local_tags)

    for workflow_path in (
        ".github/workflows/test.yml",
        ".github/workflows/pypi-publish.yml",
    ):
        workflow = _read(workflow_path)
        _require(
            workflow.count(f'SOURCE_DATE_EPOCH: "{REPRODUCIBLE_BUILD_EPOCH}"') == 1,
            (
                f"{workflow_path} must set the one canonical wheel-build "
                f"SOURCE_DATE_EPOCH to {REPRODUCIBLE_BUILD_EPOCH}"
            ),
        )

    exact_text: dict[str, tuple[str, ...]] = {
        "CITATION.cff": (f"version: {version}",),
        "docs/conf.py": (f'release = version = "{version}"  # CELLUCID_VERSION',),
        "docs/user_guide/python_package/installation.md": (
            f'pip install "cellucid=={version}"  # CELLUCID_VERSION',
        ),
        "docs/user_guide/python_package/h_developer_docs/04_build_install_and_packaging.md": (
            f'`version = "{version}"` (beta) <!-- CELLUCID_VERSION -->',
        ),
        "docs/user_guide/python_package/h_developer_docs/14_release_process.md": (
            f"`cellucid` is currently in beta (`{version}`). <!-- CELLUCID_VERSION -->",
            # The commands a maintainer copies. A bump that misses these tags and
            # publishes the previous version.
            f'git tag -a v{version} -m "Cellucid Python {version}"',
            f"git push origin v{version}",
            f"python scripts/validate_release.py --tag v{version}",
            f"python -c \"import cellucid; assert cellucid.__version__ == '{version}'\"",
            f'python -m pip install "cellucid=={version}"',
        ),
        "docs/user_guide/r_package/installation.md": (
            f"Version `{version}` from `packageVersion` <!-- CELLUCID_VERSION -->",
        ),
        (
            "docs/user_guide/python_package/c_data_preparation_api/"
            "09_output_format_specification_exports_directory.md"
        ): (f"`cellucid_data_version` is `{version}`. <!-- CELLUCID_VERSION -->",),
    }
    for relative_path, required_texts in exact_text.items():
        contents = _read(relative_path)
        for required_text in required_texts:
            _require(
                required_text in contents,
                f"{relative_path} does not declare release {version} "
                f"at {required_text!r}",
            )
    citation = _read("CITATION.cff")
    _require(
        "date-released:" not in citation,
        "CITATION.cff must not claim a release date before publication",
    )
    _require(
        f"## [{version}]\n" in changelog,
        "CHANGELOG.md current version heading must not claim a publication date",
    )

    validate_recipe(version)
    return version


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate Cellucid Python release metadata and an optional Git tag."
    )
    parser.add_argument("--tag", help="Exact pushed Git tag, including the v prefix.")
    parser.add_argument(
        "--verify-local-tags",
        action="store_true",
        help="Require every historical GitHub release link to name a fetched tag.",
    )
    parser.add_argument(
        "--sdist",
        type=Path,
        help="Normalized source distribution, or directory containing exactly one.",
    )
    arguments = parser.parse_args()
    version = validate_release(
        arguments.tag,
        verify_local_tags=arguments.verify_local_tags,
    )
    if arguments.sdist is not None:
        expected_digest = validate_recipe(version)
        validated_sdist = validate_sdist(
            arguments.sdist,
            version,
            expected_digest,
        )
        print(f"Source distribution matches the recipe: {validated_sdist}")
    print(f"Cellucid Python release contract is valid for v{version}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
