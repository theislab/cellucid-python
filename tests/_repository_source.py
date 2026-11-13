"""The environment a test subprocess needs to run this checkout's ``cellucid``.

A child interpreter or ``mypy`` run resolves ``cellucid`` through its own
import path, not the test process's ``sys.path``. Left alone it finds whatever
``cellucid`` the interpreter has installed, which may be another checkout, and
the test then exercises code that is not under test. Every subprocess that
imports or type-checks ``cellucid`` takes its environment from here, and its
wait from ``CHILD_PROCESS_FENCE``.
"""

from __future__ import annotations

import os
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[1] / "src"

#: Seconds a test waits for a child process whose work only takes longer on a loaded machine:
#: an interpreter starting, importing numpy, pandas and ``cellucid``, writing a dataset with
#: ``prepare``, or Node reading one. A wait is as long as the machine needs; this is the fence
#: against a hang (the suite runs without ``pytest-timeout``), so it is generous. A test about
#: time itself keeps its own short bound where it stands.
CHILD_PROCESS_FENCE = 600


def repository_source_environment(**overrides: str) -> dict[str, str]:
    """Return ``os.environ`` with this checkout's ``src`` first on both paths."""
    environment = os.environ.copy()
    for variable in ("PYTHONPATH", "MYPYPATH"):
        inherited = environment.get(variable)
        environment[variable] = (
            f"{SOURCE_ROOT}{os.pathsep}{inherited}" if inherited else str(SOURCE_ROOT)
        )
    environment.update(overrides)
    return environment
