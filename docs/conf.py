"""Sphinx configuration for the cellucid Python package."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))

project = "cellucid"
author = "cellucid contributors"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "myst_parser",
]

autosummary_generate = True
templates_path = ["_templates"]
exclude_patterns: list[str] = ["_build"]

html_theme = "sphinx_rtd_theme"

try:
    from cellucid import __version__
except Exception:  # pragma: no cover - import may fail during initial setup
    release = version = "0.0.0"
else:
    release = version = __version__
