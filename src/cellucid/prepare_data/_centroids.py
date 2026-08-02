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


def _compute_centroids_for_all_dimensions(
    embeddings: dict[int, np.ndarray],
    codes: np.ndarray,
    categories: list[JsonScalar],
    outlier_quantile: float = 0.95,
    min_points: int = 10,
) -> dict[int, list[dict]]:
    """
    Compute centroids for each available dimension.

    Returns a dictionary keyed by dimension (1, 2, 3) with centroid lists for each.
    """
    centroids_by_dim = {}
    for dim, coords in embeddings.items():
        centroids_by_dim[dim] = _compute_centroids_for_field(
            coords, codes, categories, outlier_quantile, min_points
        )
    return centroids_by_dim


def _compute_latent_space_quantiles(
    latent: np.ndarray,
    codes: np.ndarray,
    categories: list[JsonScalar],
    min_points: int = 10,
) -> np.ndarray:
    """
    Compute per-cell outlier quantiles based on latent space distances to category centroids.
    """
    n_cells = latent.shape[0]
    quantiles = np.full(n_cells, np.nan, dtype=np.float32)

    for code, _label in enumerate(categories):
        mask = codes == code
        idx = np.nonzero(mask)[0]
        n = idx.size

        if n < min_points:
            continue

        pts = latent[idx, :]
        centroid = pts.mean(axis=0)
        dists = np.linalg.norm(pts - centroid, axis=1)
        sorted_dists = np.sort(dists)
        ranks = np.searchsorted(sorted_dists, dists, side="right")
        cell_quantiles = ranks.astype(np.float32) / n
        quantiles[idx] = cell_quantiles

    return quantiles
