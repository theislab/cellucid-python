"""The suite tests this checkout's ``cellucid``, in process and in subprocesses.

An interpreter with another checkout installed in editable mode would
otherwise import that copy, and the suite would pass or fail on code that is
not under test.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from _repository_source import SOURCE_ROOT, repository_source_environment

import cellucid


def test_the_test_process_imports_this_checkout() -> None:
    assert Path(cellucid.__file__).resolve().is_relative_to(SOURCE_ROOT)


def test_a_test_subprocess_imports_this_checkout() -> None:
    completed = subprocess.run(
        [sys.executable, "-c", "import cellucid; print(cellucid.__file__)"],
        env=repository_source_environment(),
        capture_output=True,
        text=True,
        check=True,
    )
    assert Path(completed.stdout.strip()).resolve().is_relative_to(SOURCE_ROOT)
