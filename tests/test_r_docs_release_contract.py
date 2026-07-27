from __future__ import annotations

import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPOSITORY_ROOT / "docs"
R_GUIDE = DOCS_ROOT / "user_guide" / "r_package"


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


def test_r_installation_is_a_first_class_guide_page() -> None:
    installation_path = R_GUIDE / "installation.md"
    installation = _read(installation_path)
    home = _read(R_GUIDE / "index.md")
    landing = _read(R_GUIDE / "a_landing_pages" / "index.md")
    docs_and_validator = "\n".join(
        path.read_text(encoding="utf-8")
        for root in (DOCS_ROOT, REPOSITORY_ROOT / "scripts")
        for path in root.rglob("*")
        if path.is_file() and path.suffix in {".md", ".py"} and "_build" not in path.parts
    )

    assert installation_path.is_file()
    assert not (R_GUIDE / "a_landing_pages" / "02_installation.md").exists()
    for exact in (
        "# Installation",
        "**Active package version: 0.9.1.**",
        "install.packages(\"cellucid\")",
        "remotes::install_github(\"theislab/cellucid-r\")",
        "packageVersion(\"cellucid\")",
    ):
        assert exact in installation

    assert ":link: installation" in home
    assert "\ninstallation\na_landing_pages/index\n" in home
    assert ":link: ../installation" in landing
    assert "r_package/a_landing_pages/02_installation" not in docs_and_validator
    assert "../a_landing_pages/02_installation" not in docs_and_validator


def test_r_home_describes_the_active_submission_without_claiming_publication() -> None:
    version = _package_version()
    home = _read(R_GUIDE / "index.md")
    installation = _read(R_GUIDE / "installation.md")
    normalized_installation = " ".join(installation.split())

    assert f"Active package version — {version}" in home
    assert (
        f"Version **{version}** is the active Cellucid for R source and documentation"
        in home
    )
    assert "and the CRAN submission release" in home
    assert "is authoritative for registry availability" in home
    assert f"**Active package version: {version}.**" in installation
    assert "when the index lists `cellucid` 0.9.1" in normalized_installation
    assert "published on CRAN" not in home
    assert "available on CRAN" not in home
