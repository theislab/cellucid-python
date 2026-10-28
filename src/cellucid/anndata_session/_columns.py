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

from collections.abc import Sequence
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


def place_rows(values: pd.Series, rows: np.ndarray, index: pd.Index) -> pd.Series:
    """A column of some rows (``values``, one per row of ``rows``) placed in a
    column of every row: the others unknown -- ``<NA>`` for a membership
    (nullable ``boolean``), missing for a category, ``NaN`` for numbers."""
    rows = np.asarray(rows, dtype=np.int64)
    n_obs = len(index)
    if pd.api.types.is_bool_dtype(values.dtype):
        flags = np.zeros(n_obs, dtype=bool)
        known = np.zeros(n_obs, dtype=bool)
        flags[rows] = values.to_numpy(dtype=bool, na_value=False)
        known[rows] = values.notna().to_numpy()
        return pd.Series(pd.arrays.BooleanArray(flags, ~known), index=index, copy=False)
    if isinstance(values.dtype, pd.CategoricalDtype):
        codes = np.full(n_obs, -1, dtype=np.int64)
        codes[rows] = values.cat.codes.to_numpy()
        return pd.Series(pd.Categorical.from_codes(codes, dtype=values.dtype), index=index)
    numbers = np.full(n_obs, np.nan, dtype=values.to_numpy().dtype if values.dtype.kind == "f" else np.float64)
    numbers[rows] = values.to_numpy()
    return pd.Series(numbers, index=index)


def write_columns(target: Any, plans: Sequence[_ColumnPlan], *, remove: Sequence[str] = ()) -> None:
    """Put every planned column on ``target.obs`` (replacing a column of the
    same name) and drop the columns in ``remove``, in one assignment: the
    columns are all built before this runs, and ``obs`` changes once."""
    next_obs = target.obs.copy(deep=False)
    for name in remove:
        if name in next_obs.columns:
            next_obs = next_obs.drop(columns=[name])
    for plan in plans:
        next_obs[plan.name] = plan.values
    target.obs = next_obs
