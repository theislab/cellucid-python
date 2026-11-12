"""
Scanpy's colors of an obs column: ``uns['<column>_colors']``, one color a
category in the column's category order (``sc.pl.*`` writes them; the viewer
opens a column in them).

One implementation for every writer that carries them: ``prepare()``
(:mod:`cellucid.prepare_data`, a column written keeps its colors) and the
session/pull writers (:mod:`cellucid.anndata_session._columns`, a column's
colors follow it). cellucid-r keeps the same pair in ``R/round-trip.R``
(``.colors_key``, ``.text_colors``).
"""

from __future__ import annotations

from typing import Any

import numpy as np


def colors_key(column: str) -> str:
    """Where Scanpy keeps a column's colors: ``uns['<column>_colors']``."""
    return f"{column}_colors"


def text_colors(value: Any) -> list[str] | None:
    """Colors in Scanpy's form, a one-dimensional list of text (a list, a
    tuple or an array), as a list; None for anything else."""
    if isinstance(value, np.ndarray):
        if value.ndim != 1:
            return None
        items = value.tolist()
    elif isinstance(value, list | tuple):
        items = list(value)
    else:
        return None
    return items if all(isinstance(item, str) for item in items) else None
