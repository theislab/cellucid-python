from __future__ import annotations

from importlib.metadata import version
from pathlib import Path

from packaging.requirements import Requirement
from packaging.version import Version

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_core_metadata_declares_one_current_bounded_direct_data_runtime() -> None:
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    anndata_declaration = Requirement("anndata>=0.12.19,<0.13")
    zarr_declaration = Requirement("zarr>=3.1.4,<4")
    numcodecs_declaration = Requirement("numcodecs>=0.16.3,<0.17")
    assert pyproject.count('"anndata>=0.12.19,<0.13",') == 1
    assert pyproject.count('"zarr>=3.1.4,<4",') == 1
    assert pyproject.count('"numcodecs>=0.16.3,<0.17",') == 1
    assert Version(version("anndata")) in anndata_declaration.specifier
    assert Version(version("zarr")) in zarr_declaration.specifier
    assert Version(version("numcodecs")) in numcodecs_declaration.specifier


