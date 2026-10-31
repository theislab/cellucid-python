"""
Latent-space outlier quantiles, which :func:`cellucid.prepare` stores with
the latent space it is given.

For each categorical obs column, each cell's quantile says how far it lies
from the centre of its category in the latent space (an ``obsm`` array such
as scVI or PCA, never the 2D/3D layout): the share of its category at that
distance or closer. The web app's outlier filter hides the cells above a
threshold. The rule, for each category with at least ``min_points`` cells,
over the latent values as stored (float32, widened exactly to float64):

1. the centroid is the column sums, accumulated cell by cell in cell order,
   divided by the cell count;
2. a cell's distance is the square root of its squared differences from the
   centroid, accumulated column by column in column order;
3. its quantile is the number of the category's cells at that distance or
   closer, divided by the cell count, both as float32.

A cell whose latent row is not finite (NaN or an infinity in any column: a
latent space computed on some of the cells, a modality that did not measure
the cell) has no place in the latent space: it is left out as a cell with a
missing value is, adding nothing to its category's centroid or size, and
``min_points`` counts a category's placed cells. Every other cell (a missing
value, an unplaced cell, or a category under ``min_points``) is NaN. The web
app computes the same bytes in the browser for a file that stores
none (``assets/js/data/latent-outliers.js``), cellucid-r's
``cellucid_prepare()`` stores the same bytes, and both are checked against an
independent numpy reference (``scripts/generate-latent-outlier-conformance.py``
in the web repository).

Which columns: every obs column the viewer shows as categorical -- a
categorical, a text column, a boolean column, a nullable boolean column --
partitioned by its values (a missing value belongs to no category). A column
none of whose categories reaches ``min_points`` cells gets no entry.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

#: The fewest cells a category needs for its cells to get outlier quantiles,
#: unless ``centroid_min_points`` says otherwise; the web app's default too.
DEFAULT_MIN_POINTS = 10


def latent_outlier_quantiles(
    latent: np.ndarray, codes: np.ndarray, n_categories: int, min_points: int
) -> np.ndarray:
    """One float32 quantile per cell by the module rule; NaN where none.

    ``latent`` is ``(n_cells, width)`` float32; ``codes`` one integer per
    cell, a code outside ``[0, n_categories)`` being a missing value, as is a
    cell whose latent row is not finite.
    """
    values = np.asarray(latent)
    if values.dtype != np.float32 or values.ndim != 2:
        raise TypeError("The latent space must be a two-dimensional float32 array.")
    n_cells, width = values.shape
    codes = np.asarray(codes, dtype=np.int64)
    if codes.shape != (n_cells,):
        raise ValueError(f"Category codes must be one integer per cell ({n_cells}).")
    quantiles = np.full(n_cells, np.nan, dtype=np.float32)
    valid = (codes >= 0) & (codes < n_categories) & np.isfinite(values).all(axis=1)
    counts = np.bincount(codes[valid], minlength=n_categories)
    scored = counts >= min_points
    # The scored cells, in cell order, and their categories.
    cells = np.flatnonzero(valid & scored[np.where(valid, codes, 0)])
    if cells.size == 0:
        return quantiles
    group = codes[cells]
    n_of_cell = counts[group]

    # Step 1, then step 2, column by column. ``bincount`` with weights adds
    # them one by one in the order given, so each category's sum is
    # accumulated cell by cell in cell order (``np.sum`` adds pairwise).
    squared = np.zeros(cells.size, dtype=np.float64)
    for column in range(width):
        column_values = values[cells, column].astype(np.float64)
        centroids = np.bincount(group, weights=column_values, minlength=n_categories) / np.maximum(counts, 1)
        difference = column_values - centroids[group]
        squared += difference * difference
    distances = np.sqrt(squared)

    # Step 3: within each category, the cells at a distance or closer are the
    # cells up to the end of that distance's run in (category, distance) order.
    order = np.lexsort((distances, group))
    sorted_group = group[order]
    sorted_distance = distances[order]
    run_starts = np.empty(cells.size, dtype=bool)
    run_starts[0] = True
    run_starts[1:] = (sorted_group[1:] != sorted_group[:-1]) | (sorted_distance[1:] != sorted_distance[:-1])
    run_of = np.cumsum(run_starts) - 1
    run_ends = np.append(np.flatnonzero(run_starts)[1:], cells.size)
    # A category's first position in that order: the scored cells before it.
    scored_before = np.concatenate(([0], np.cumsum(np.where(scored, counts, 0))))
    at_most = np.empty(cells.size, dtype=np.int64)
    at_most[order] = run_ends[run_of] - scored_before[sorted_group]
    quantiles[cells] = at_most.astype(np.float32) / n_of_cell.astype(np.float32)
    return quantiles


def _partition(value: Any) -> tuple[np.ndarray, int] | None:
    """A stored obs column as category codes (-1 missing) and their count; None if not categorical."""
    if isinstance(value, pd.Categorical):
        return np.asarray(value.codes, dtype=np.int64), len(value.categories)
    if isinstance(value, pd.arrays.BooleanArray):
        codes = value.to_numpy(dtype=np.int64, na_value=-1)
        return codes, 2
    if isinstance(value, np.ndarray) and value.dtype.kind == "b":
        return value.astype(np.int64), 2
    if isinstance(value, np.ndarray) and value.dtype.kind == "O":
        _, codes = np.unique(value.astype(str), return_inverse=True)
        return codes.astype(np.int64).reshape(-1), int(codes.max(initial=-1)) + 1
    return None


def stored_outlier_quantiles(
    obs: dict[str, Any], latent_key: str, latent: np.ndarray, min_points: int
) -> dict[str, Any]:
    """``uns['cellucid']['outlierQuantiles']``: every categorical obs column's quantiles in ``latent``."""
    fields: dict[str, np.ndarray] = {}
    placed = np.isfinite(latent).all(axis=1)
    for key, value in obs.items():
        partition = _partition(value)
        if partition is None:
            continue
        codes, n_categories = partition
        valid = codes[(codes >= 0) & (codes < n_categories) & placed]
        if not (np.bincount(valid, minlength=n_categories) >= min_points).any():
            continue
        fields[key] = latent_outlier_quantiles(latent, codes, n_categories, min_points)
    return {"latentKey": latent_key, "minPoints": min_points, "fields": fields}
