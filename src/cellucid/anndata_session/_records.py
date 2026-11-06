"""
The records the planners hand back.

:class:`_ColumnPlan` is one prepared ``obs`` column a pull writes;
:class:`_Candidate` is one column a session's step holds (a highlight group or
a user-defined field), with its values when the session applies it and why
not when it does not; the apply plan settles each against the object
(:mod:`._apply`). Nothing is written until every candidate is built, which is
what makes an application atomic.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class _ColumnPlan:
    name: str
    values: pd.Series
    metadata: dict[str, Any]


@dataclass(frozen=True)
class _Candidate:
    """One column a step's state holds.

    ``kind`` is the plan row's (``user-field`` or ``highlight``), ``name`` the
    column it is written to, ``values`` the column (None when the session
    itself does not apply it: a field deleted in it, ``note`` says why),
    ``cells`` the cells it marks or labels, ``metadata`` what a pull records
    beside it (its result's kind and id).
    """

    kind: str
    domain: str
    name: str
    values: pd.Series | None
    cells: int | None
    note: str | None
    metadata: dict[str, Any]
