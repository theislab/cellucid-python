"""
The ``obs`` columns Cellucid writes onto an AnnData (or a MuData's global
``obs``), and the one writer that puts them there.

Both ways results come back from the viewer use these: a saved session
(:func:`~cellucid.apply_cellucid_session_to_anndata`) and a pull from a live
viewer or the web app's exports (:func:`~cellucid.pull`). One rule each:

- a highlight group, the filter's visible cells: ``bool``, ``True`` for its
  cells; when the viewer opened a subsample, the nullable ``boolean`` dtype
  with ``<NA>`` for every cell it did not open (unknown, not "not
  selected");
- a categorical field: ``category`` with the viewer's categories in its order,
  missing (``NaN``) where a cell has no category or was not open;
- a continuous field or score: ``float32`` (``float64`` when the values came
  as such), ``NaN`` where a cell has no value or was not open.

Cells are original positions: rows of the object's ``obs`` in order.
:func:`write_columns` puts every planned column on at once, after every one
was built, so a failure leaves the object as it was.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .._plural import plural
from ._records import _ColumnPlan


def membership_column(
    n_obs: int,
    cells: np.ndarray | None,
    *,
    open_cells: np.ndarray | None,
    index: pd.Index,
) -> pd.Series:
    """``True`` for ``cells`` (every cell when None); ``<NA>`` outside ``open_cells``."""
    mask = np.ones(n_obs, dtype=bool) if cells is None else np.zeros(n_obs, dtype=bool)
    if cells is not None:
        mask[np.asarray(cells, dtype=np.int64)] = True
    if open_cells is None:
        return pd.Series(mask, index=index, copy=False)
    known = np.zeros(n_obs, dtype=bool)
    known[np.asarray(open_cells, dtype=np.int64)] = True
    if cells is None:
        # Every open cell: the visible set of a subsample is its open cells.
        mask = known.copy()
    values = pd.arrays.BooleanArray(mask, ~known)
    return pd.Series(values, index=index, copy=False)


def categorical_column(
    codes: np.ndarray,
    categories: Sequence[Any],
    *,
    cells: np.ndarray | None,
    n_obs: int,
    missing_code: int | None,
    index: pd.Index,
) -> pd.Series:
    """A categorical from codes over ``cells`` (every cell when None).

    ``missing_code`` is the code meaning "no category" (a code type's largest
    value), or None when every code names one. A cell outside ``cells`` has
    no category.
    """
    category_codes = np.asarray(codes).astype(np.int64, copy=True)
    if missing_code is not None:
        category_codes[category_codes == missing_code] = -1
    invalid = np.flatnonzero((category_codes < -1) | (category_codes >= len(categories)))
    if invalid.size:
        raise ValueError(
            f"Code {int(category_codes[invalid[0]])} at position {int(invalid[0])} names no category "
            f"of the {len(categories)} declared"
        )
    if cells is not None:
        placed = np.full(n_obs, -1, dtype=np.int64)
        placed[np.asarray(cells, dtype=np.int64)] = category_codes
        category_codes = placed
    elif category_codes.size != n_obs:
        raise ValueError(f"{category_codes.size} {plural(category_codes.size, 'code')} do not cover the {n_obs} {plural(n_obs, 'cell')}")
    categorical = pd.Categorical.from_codes(
        category_codes.tolist() if category_codes.size == 0 else category_codes,
        categories=pd.Index(list(categories)),
        ordered=False,
    )
    return pd.Series(categorical, index=index)


def continuous_column(
    values: np.ndarray,
    *,
    cells: np.ndarray | None,
    n_obs: int,
    index: pd.Index,
) -> pd.Series:
    """Numbers over ``cells`` (every cell when None), ``NaN`` elsewhere."""
    array = np.asarray(values)
    dtype = np.float64 if array.dtype == np.float64 else np.float32
    if cells is None:
        if array.size != n_obs:
            raise ValueError(f"{array.size} {plural(array.size, 'value')} do not cover the {n_obs} {plural(n_obs, 'cell')}")
        return pd.Series(array.astype(dtype, copy=True), index=index)
    placed = np.full(n_obs, np.nan, dtype=dtype)
    placed[np.asarray(cells, dtype=np.int64)] = array
    return pd.Series(placed, index=index)


def without_unknown(values: pd.Series, known: np.ndarray | None) -> pd.Series:
    """A column with every cell outside ``known`` (ascending original
    positions; None: every cell) missing: no category for a categorical,
    ``NaN`` for numbers. A session that knows the state of some cells only
    (the cells the viewer had open) says nothing about the others."""
    if known is None:
        return values
    unknown = np.ones(len(values), dtype=bool)
    unknown[np.asarray(known, dtype=np.int64)] = False
    if isinstance(values.dtype, pd.CategoricalDtype):
        codes = values.cat.codes.to_numpy().astype(np.int64, copy=True)
        codes[unknown] = -1
        return pd.Series(pd.Categorical.from_codes(codes, dtype=values.dtype), index=values.index)
    numbers = values.to_numpy(copy=True)
    if numbers.dtype.kind != "f":
        raise TypeError("Only a categorical or a floating-point column has a missing value for an unknown cell")
    numbers[unknown] = np.nan
    return pd.Series(numbers, index=values.index)


def write_columns(
    target: Any,
    plans: Sequence[_ColumnPlan],
    *,
    remove: Sequence[str] = (),
    replace: Sequence[_ColumnPlan] = (),
    rename: Mapping[str, str] | None = None,
) -> None:
    """Change ``target.obs`` in one assignment: drop the columns in ``remove``,
    put each of ``replace`` in place of the column of its name, rename the
    columns ``rename`` names (old name to new, all at once, each keeping its
    place), then put every planned column on (replacing a column of the same
    name). Every column is built before this runs, so ``obs`` changes once."""
    next_obs = target.obs.copy(deep=False)
    for name in remove:
        if name in next_obs.columns:
            next_obs = next_obs.drop(columns=[name])
    for plan in replace:
        next_obs[plan.name] = plan.values
    if rename:
        next_obs = next_obs.rename(columns=dict(rename))
    for plan in plans:
        next_obs[plan.name] = plan.values
    target.obs = next_obs
