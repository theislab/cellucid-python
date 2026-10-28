"""The files Python still maps, found on any operating system.

Windows removes no file that is memory-mapped, and lets no other program
replace one: a served store must let every mapping go before it removes its
scratch folder, and let its source file go when it closes. Other systems
remove a mapped file without complaint, so a test asserts there is no
mapping left, which fails on every system where Windows would fail.
"""

from __future__ import annotations

import gc
from pathlib import Path

import numpy as np


def mapped_files(under: Path) -> list[Path]:
    """The files under ``under`` (a folder, or one file) that a live ``np.memmap`` maps.

    No collection is run first: an array left in a reference cycle is still
    mapped until a collection happens to run, which is not a release.
    """
    root = under.resolve()
    found = []
    for value in gc.get_objects():
        if isinstance(value, np.memmap) and value.filename is not None:
            path = Path(value.filename).resolve()
            if path == root or path.is_relative_to(root):
                found.append(path)
    return sorted(set(found))
