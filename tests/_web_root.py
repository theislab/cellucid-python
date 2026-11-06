"""The web app's and cellucid-r's checkouts, for the tests that compare with them or run their readers.

``CELLUCID_WEB_ROOT`` names it (a worktree of the ``cellucid`` repository);
without it, the ``cellucid`` checkout beside this repository is used. A root
that is named, or present, must hold what a test reads from it: a file it
lacks fails the test (a shared file renamed or deleted on one side is a
divergence, never a skip). Only when no web checkout exists at all (this
repository checked out alone, as in its own CI) is the test skipped, saying so.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

REPOSITORY = Path(__file__).resolve().parents[1]


def web_root() -> Path | None:
    """The web app's checkout, or None when there is none."""
    named = os.environ.get("CELLUCID_WEB_ROOT")
    if named:
        root = Path(named)
        if not (root / "assets" / "js").is_dir():
            pytest.fail(f"CELLUCID_WEB_ROOT ({root}) is not a checkout of the web app (it has no assets/js).")
        return root
    root = REPOSITORY.parent / "cellucid"
    return root if (root / "assets" / "js").is_dir() else None


def web_file(relative: str) -> Path:
    """A file of the web app's checkout: skipped without a checkout, failed when the checkout lacks it."""
    root = web_root()
    if root is None:
        pytest.skip("The web app is not checked out beside this repository (set CELLUCID_WEB_ROOT to compare with it).")
    path = root / relative
    if not path.is_file():
        pytest.fail(f"The web app's checkout ({root}) has no {relative}.")
    return path


def sibling_file(repository: str, relative: str) -> Path:
    """A file of a sibling repository beside this one (``cellucid-r``): skipped
    without that checkout, failed when the checkout lacks it."""
    root = REPOSITORY.parent / repository
    if not root.is_dir():
        pytest.skip(f"{repository} is not checked out beside this repository.")
    path = root / relative
    if not path.is_file():
        pytest.fail(f"The {repository} checkout ({root}) has no {relative}.")
    return path
