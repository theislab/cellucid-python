"""The environment a test subprocess needs to run this checkout's ``cellucid``.

A child interpreter or ``mypy`` run resolves ``cellucid`` through its own
import path, not the test process's ``sys.path``. Left alone it finds whatever
``cellucid`` the interpreter has installed, which may be another checkout, and
the test then exercises code that is not under test. Every subprocess that
imports or type-checks ``cellucid`` takes its environment from here.
"""

from __future__ import annotations

import os
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[1] / "src"


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
