from __future__ import annotations

import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPOSITORY_ROOT / "docs"
PYTHON_GUIDE = DOCS_ROOT / "user_guide" / "python_package"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _package_version() -> str:
    pyproject = _read(REPOSITORY_ROOT / "pyproject.toml")
    match = re.search(
        r'(?m)^version = "([^"]+)"\s+# CELLUCID_VERSION$',
        pyproject,
    )
    assert match is not None
    return match.group(1)


def test_python_installation_is_a_first_class_guide_page() -> None:
    installation_path = PYTHON_GUIDE / "installation.md"
    installation = _read(installation_path)
    home = _read(PYTHON_GUIDE / "index.md")
    landing = _read(PYTHON_GUIDE / "a_landing_pages" / "index.md")

    assert installation_path.is_file()
    assert not (PYTHON_GUIDE / "a_landing_pages" / "02_installation.md").exists()
    for exact in (
        "# Installation",
        "macOS/Linux:",
        "Windows PowerShell:",
        "Windows Command Prompt:",
        "pip install cellucid",
        "python -c \"import cellucid; print(cellucid.__version__)\"",
        "cellucid --version",
        "cellucid serve --help",
    ):
        assert exact in installation

    assert ":link: installation" in home
    assert "\ninstallation\na_landing_pages/index\n" in home
    assert ":link: ../installation" in landing


def test_python_home_and_changelog_publish_exact_active_release_status() -> None:
    version = _package_version()
    home = _read(PYTHON_GUIDE / "index.md")
    installation = _read(PYTHON_GUIDE / "installation.md")
    changelog = _read(REPOSITORY_ROOT / "CHANGELOG.md")

    assert f"Active package version — {version}" in home
    assert (
        f"Version **{version}** is the current Cellucid Python source and documentation"
        in home
    )
    assert f"**Active package version: {version}.**" in installation
    assert "version and the PyPI submission release" in home
    assert "is authoritative for" in home
    assert f"## [{version}]\n" in changelog
    assert f"Version {version} is the PyPI submission release." in changelog
    assert "## [Unreleased]" not in changelog
    assert "[Unreleased]:" not in changelog
