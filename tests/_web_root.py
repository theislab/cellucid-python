"""The sibling repositories' checkouts, for the tests that compare with them or run their readers.

``CELLUCID_WEB_ROOT`` names the web app's checkout (a worktree of the
``cellucid`` repository) and ``CELLUCID_R_ROOT`` cellucid-r's; without them,
the checkout beside this repository is used. A checkout that is there must
hold what a test reads from it: a file it lacks fails the test (a shared file
renamed or deleted on one side is a divergence, never a skip), and a named
folder that is not a checkout of that repository fails it too.

Where a sibling is not there (this repository checked out alone, an sdist, an
installed package, a variable naming a folder that does not exist) the test is
skipped, saying where it looked, unless ``CELLUCID_REQUIRE_SIBLINGS`` is
``1``: this repository's CI checks out both siblings at ``main`` and sets it,
so there a missing checkout fails the test instead of turning the comparison
into a skip.

A shared file is held to its other copy byte for byte by
:func:`assert_shared_file`. The web repository is the home of every file the
three repositories share; this repository is the home of the files only the
two packages share, which cellucid-r copies. No digest is written anywhere:
when two copies differ, one command in the web repository makes every copy of
every shared file agree again (:data:`SYNC_SHARED_FILES`).
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from _repository_source import CHILD_PROCESS_FENCE

REPOSITORY = Path(__file__).resolve().parents[1]

#: Each sibling repository: the variable that names its checkout, and the
#: entry a checkout of it must hold.
_SIBLINGS = {
    "cellucid": ("CELLUCID_WEB_ROOT", Path("assets") / "js"),
    "cellucid-r": ("CELLUCID_R_ROOT", Path("DESCRIPTION")),
}


def siblings_required() -> bool:
    """Whether a missing sibling checkout (or Node) fails a test instead of skipping it."""
    return os.environ.get("CELLUCID_REQUIRE_SIBLINGS") == "1"


def _not_there(repository: str) -> str:
    """Where ``repository`` was looked for and not found, as one sentence."""
    variable, _ = _SIBLINGS[repository]
    named = os.environ.get(variable)
    if named:
        return f"{repository} is not there: {variable} names {named}, which does not exist."
    return f"{repository} is not checked out beside this repository (set {variable} to name its checkout)."


def sibling_root(repository: str) -> Path | None:
    """The checkout of ``repository`` (``cellucid`` for the web app, or
    ``cellucid-r``), or None when it is not there and none is required."""
    variable, marker = _SIBLINGS[repository]
    named = os.environ.get(variable)
    root = Path(named) if named else REPOSITORY.parent / repository
    if (root / marker).exists():
        return root
    if named and root.exists():
        pytest.fail(f"{variable} ({root}) is not a checkout of {repository} (it has no {marker.as_posix()}).")
    if siblings_required():
        pytest.fail(f"CELLUCID_REQUIRE_SIBLINGS is 1, and {_not_there(repository)}")
    return None


def sibling_file(repository: str, relative: str) -> Path:
    """A file of a sibling repository's checkout: skipped without that
    checkout (unless one is required), failed when the checkout lacks it."""
    root = sibling_root(repository)
    if root is None:
        pytest.skip(_not_there(repository))
    path = root / relative
    if not path.is_file():
        pytest.fail(f"The {repository} checkout ({root}) has no {relative}.")
    return path


def sibling_dir(repository: str, relative: str) -> Path:
    """A directory of a sibling repository's checkout, as :func:`sibling_file`."""
    root = sibling_root(repository)
    if root is None:
        pytest.skip(_not_there(repository))
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
    if siblings_required():
        pytest.fail(f"CELLUCID_REQUIRE_SIBLINGS is 1 but Node is not installed: {purpose} cannot run.")
    pytest.skip(f"Node is not installed: {purpose} cannot run.")


def web_file(relative: str) -> Path:
    """A file of the web app's checkout, as :func:`sibling_file`."""
    return sibling_file("cellucid", relative)


#: The one command that makes every copy of every shared file agree, run in
#: the web repository: it writes the generated files again and copies every
#: home into the packages.
SYNC_SHARED_FILES = "node scripts/sync-shared-files.mjs"


def assert_shared_file(ours: str, theirs: str, *, repository: str = "cellucid", then: str = "") -> None:
    """This repository's ``ours`` is byte for byte ``repository``'s ``theirs``
    (both relative to their repository's root): the web app's file for a file
    whose home is the web repository, cellucid-r's copy for one whose home is
    here. Skipped without that checkout, unless the siblings are required
    (:func:`sibling_file`). ``then`` names what to run after the copies agree
    again, for a file generated code is written from."""
    other = sibling_file(repository, theirs)
    if (REPOSITORY / ours).read_bytes() != other.read_bytes():
        pytest.fail(
            f"{ours} is not byte for byte {repository}'s {theirs}. Run `{SYNC_SHARED_FILES}` in the web "
            f"repository (cellucid): it makes every copy of every shared file agree{then}."
        )


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
        capture_output=True, text=True, encoding="utf-8", check=False, timeout=CHILD_PROCESS_FENCE,
    )
    if completed.returncode != 0:
        pytest.fail(f"{script} {function}() failed: {completed.stderr}")
    return completed.stdout
