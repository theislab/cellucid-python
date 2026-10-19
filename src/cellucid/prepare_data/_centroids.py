"""
Category centroids and latent-space outlier quantiles.

Centroids are measured on the normalized embedding, because they are label
positions the viewer draws; the outlier quantile of a cell is measured in the
latent space, because that is where distance to a category's centre means
something biological.
"""

import numpy as np

from .._contracts import JsonScalar


def _compute_centroids_for_field(
    coords: np.ndarray,
    codes: np.ndarray,
    categories: list[JsonScalar],
    outlier_quantile: float = 0.95,
    min_points: int = 10,
) -> list[dict]:
    """
    Compute centroids per category with outlier removal (based on embedding coords for display).

    Works with any dimensionality (1D, 2D, 3D, etc.) - the position will have the same
    number of dimensions as the input coords.
    """
    if coords.shape[0] != codes.shape[0]:
        raise ValueError("coords and codes must have the same length.")

    centroids: list[dict[str, object]] = []

    if not (0.5 < outlier_quantile < 1.0):
        raise ValueError("outlier_quantile must be strictly between 0.5 and 1.0")

    for code, label in enumerate(categories):
        mask = codes == code
        idx = np.nonzero(mask)[0]
        n = idx.size
        if n < min_points:
            continue

        pts = coords[idx, :]  # (n, ndim)
        center = pts.mean(axis=0)

        if n > min_points:
            dists = np.linalg.norm(pts - center, axis=1)
            thr = float(np.quantile(dists, outlier_quantile))
            inlier_mask = dists <= thr
            n_in = int(inlier_mask.sum())
            if n_in >= min_points:
                pts_in = pts[inlier_mask, :]
                center = pts_in.mean(axis=0)
                used_count = n_in
            else:
                used_count = n
        else:
            used_count = n

        centroids.append(
            {
                "category": label,
                "position": center.astype(float).tolist(),
                "n_points": int(used_count),
            }
        )

    return centroids


def _compute_centroids_for_all_embeddings(
    coordinates_by_id: dict[str, np.ndarray],
    codes: np.ndarray,
    categories: list[JsonScalar],
    outlier_quantile: float = 0.95,
    min_points: int = 10,
) -> dict[str, list[dict]]:
    """
    Compute centroids on every exported embedding.

    Returns a dictionary keyed by embedding id, in the order the ids were given.
    """
    return {
        embedding_id: _compute_centroids_for_field(
            coords, codes, categories, outlier_quantile, min_points
        )
        for embedding_id, coords in coordinates_by_id.items()
    }


def _compute_latent_space_quantiles(
    latent: np.ndarray,
    codes: np.ndarray,
    categories: list[JsonScalar],
    min_points: int = 10,
) -> np.ndarray:
    """Per-cell outlier quantiles: how far a cell lies from its category's centre.

    The one rule ``prepare()``, the AnnData server and the web app's
    ``data/latent-outliers.js`` share, bit for bit. For each category with at
    least ``min_points`` cells, over the latent values as stored (float32,
    widened exactly to float64):

    1. the centroid is the column sums, accumulated cell by cell in cell
       order, divided by the cell count;
    2. a cell's distance is the square root of its squared differences from
       the centroid, accumulated column by column in column order;
    3. its quantile is the number of the category's cells at that distance or
       closer, divided by the cell count, both as float32.

    Every other cell (a missing code, or a category under ``min_points``) is
    NaN.
    """
    n_cells = latent.shape[0]
    quantiles = np.full(n_cells, np.nan, dtype=np.float32)
    values = np.asarray(latent, dtype=np.float32)

    for code, _label in enumerate(categories):
        idx = np.nonzero(codes == code)[0]
        n = idx.size
        if n < min_points:
            continue

        # One column at a time, so a category never holds more than a few
        # float64 vectors. A cumulative sum adds in cell order (numpy's sum
        # would add pairwise), so every implementation can reproduce it.
        squared = np.zeros(n, dtype=np.float64)
        for column in range(values.shape[1]):
            column_values = values[idx, column].astype(np.float64)
            centroid = np.cumsum(column_values)[-1] / n
            difference = column_values - centroid
            squared += difference * difference
        dists = np.sqrt(squared)
        sorted_dists = np.sort(dists)
        ranks = np.searchsorted(sorted_dists, dists, side="right")
        quantiles[idx] = ranks.astype(np.float32) / np.float32(n)

    return quantiles
