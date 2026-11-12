"""The sibling repositories' checkouts, for the tests that compare with them or run their readers.

``CELLUCID_WEB_ROOT`` names the web app's checkout (a worktree of the
``cellucid`` repository) and ``CELLUCID_R_ROOT`` cellucid-r's; without them,
the checkout beside this repository is used. A root that is named, or present,
must hold what a test reads from it: a file it lacks fails the test (a shared
file renamed or deleted on one side is a divergence, never a skip).

When no checkout of a sibling exists at all (this repository checked out
alone) the test is skipped, saying so, unless ``CELLUCID_SIBLINGS_REQUIRED``
is set: this repository's CI checks out both siblings and sets it, so there a
missing checkout fails the test instead of turning the comparison into a skip.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

REPOSITORY = Path(__file__).resolve().parents[1]

#: Each sibling repository: the variable that names its checkout, and the
#: entry a checkout of it must hold.
_SIBLINGS = {
    "cellucid": ("CELLUCID_WEB_ROOT", Path("assets") / "js"),
    "cellucid-r": ("CELLUCID_R_ROOT", Path("DESCRIPTION")),
}


def sibling_root(repository: str) -> Path | None:
    """The checkout of ``repository`` (``cellucid`` for the web app, or
    ``cellucid-r``), or None when there is none and none is required."""
    variable, marker = _SIBLINGS[repository]
    named = os.environ.get(variable)
    if named:
        root = Path(named)
        if not (root / marker).exists():
            pytest.fail(f"{variable} ({root}) is not a checkout of {repository} (it has no {marker.as_posix()}).")
        return root
    root = REPOSITORY.parent / repository
    if (root / marker).exists():
        return root
    if os.environ.get("CELLUCID_SIBLINGS_REQUIRED"):
        pytest.fail(
            f"CELLUCID_SIBLINGS_REQUIRED is set but {repository} is not checked out beside this repository "
            f"(or named by {variable})."
        )
    return None


def web_root() -> Path | None:
    """The web app's checkout, or None when there is none."""
    return sibling_root("cellucid")


def sibling_file(repository: str, relative: str) -> Path:
    """A file of a sibling repository's checkout: skipped without that
    checkout (unless one is required), failed when the checkout lacks it."""
    root = sibling_root(repository)
    if root is None:
        pytest.skip(f"{repository} is not checked out beside this repository (set {_SIBLINGS[repository][0]} to compare with it).")
    path = root / relative
    if not path.is_file():
        pytest.fail(f"The {repository} checkout ({root}) has no {relative}.")
    return path


def sibling_dir(repository: str, relative: str) -> Path:
    """A directory of a sibling repository's checkout, as :func:`sibling_file`."""
    root = sibling_root(repository)
    if root is None:
        pytest.skip(f"{repository} is not checked out beside this repository (set {_SIBLINGS[repository][0]} to compare with it).")
    path = root / relative
    if not path.is_dir():
        pytest.fail(f"The {repository} checkout ({root}) has no {relative}/.")
    return path


def node_executable(purpose: str) -> str:
    """The ``node`` that runs the web app's code for ``purpose``: skipped
    without Node, unless the siblings are required (as in this repository's
    CI), where its absence fails the test."""
    import shutil

    node = shutil.which("node")
    if node is not None:
        return node
    if os.environ.get("CELLUCID_SIBLINGS_REQUIRED"):
        pytest.fail(f"CELLUCID_SIBLINGS_REQUIRED is set but Node is not installed: {purpose} cannot run.")
    pytest.skip(f"Node is not installed: {purpose} cannot run.")


def web_file(relative: str) -> Path:
    """A file of the web app's checkout, as :func:`sibling_file`."""
    return sibling_file("cellucid", relative)


#: Runs one export of a web generator script and prints what it returns.
_RENDER = (
    "import { pathToFileURL } from 'node:url';"
    # argv[1] is a word, never the script: the script runs its own main()
    # when argv[1] names it.
    "const [script, name] = process.argv.slice(2);"
    "const module = await import(pathToFileURL(script).href);"
    "process.stdout.write(module[name]());"
)


def rendered_by_web(script: str, function: str) -> str:
    """What the web app's generator ``script`` renders now from its source
    (``renderHistoryModel`` from ``app/history/model.js``,
    ``renderActionCatalog`` from the action definitions): a shared file this
    package holds is compared with the source itself, not only with the web
    app's committed copy, so a source changed without regenerating the copies
    fails here too. Skipped without a web checkout or Node, unless the
    siblings are required."""
    import subprocess

    path = web_file(script)
    node = node_executable("the web app's generator")
    completed = subprocess.run(
        [node, "--input-type=module", "-e", _RENDER, "render", str(path), function],
        capture_output=True, text=True, encoding="utf-8", check=False, timeout=120,
    )
    if completed.returncode != 0:
        pytest.fail(f"{script} {function}() failed: {completed.stderr}")
    return completed.stdout
