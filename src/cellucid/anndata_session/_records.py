"""
The records the planners hand back.

:class:`ApplySummary` is the result of one successful application.
:class:`_ColumnPlan` is one prepared ``obs`` column, held un-applied until every
plan has been built, which is what makes application atomic.
:class:`_CurrentChunkInventory` is the validated chunk listing.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class ApplySummary:
    """Columns materialized by one successful, atomic session application."""

    added_obs_columns: list[str]


@dataclass(frozen=True)
class _ColumnPlan:
    name: str
    values: pd.Series
    metadata: dict[str, Any]


@dataclass(frozen=True)
class _CurrentChunkInventory:
    ids: set[str]
    metadata_by_id: dict[str, dict[str, Any]]
    analysis_artifact_ids: tuple[str, ...]
