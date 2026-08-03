"""
Numeric input contracts and the one coordinate normalization.

Every array the exporter accepts passes through here first, so the rule that
values must stay finite once they reach the viewer's float32 domain is stated
once. Embeddings are the exception that proves it: their extents are measured
at the input's own precision, because casting first collapses any axis whose
spread is finer than float32 resolution at its magnitude.
"""

from typing import cast

import numpy as np
from scipy import sparse

from ..continuous_payload_diagnosis import describe_non_finite


def _require_finite_float32_array(
    values: object,
    *,
    label: str,
) -> np.ndarray:
    """Validate real finite values that remain finite after float32 conversion."""
    dense = cast(sparse.spmatrix, values).toarray() if sparse.issparse(values) else values
    array = np.asarray(dense)
    if array.dtype.kind not in {"i", "u", "f"}:
        raise TypeError(f"{label} must contain real numeric values.")
    if not np.isfinite(array).all():
        raise ValueError(
            f"{label} must contain only finite values. "
            + describe_non_finite(array, label)
        )
    with np.errstate(over="ignore", invalid="ignore"):
        float32_array = array.astype(np.float32, copy=False)
    if not np.isfinite(float32_array).all():
        raise ValueError(
            f"{label} contains values outside the finite float32 range. "
            + describe_non_finite(array, label)
        )
    # Underflow is the other end of the same range, and it is the quiet one: a
    # nonzero value below the smallest float32 subnormal becomes 0.0, which is
    # finite, so an overflow-only check publishes a field whose every value has
    # been replaced by zero. Measured on eight distinct inputs near 1e-320: one
    # distinct value written. The R writer refuses exactly this, and one export
    # format cannot mean two things.
    if np.any((array != 0) & (float32_array == 0)):
        raise ValueError(
            f"{label} contains values outside the finite float32 range. "
            + describe_non_finite(array, label)
        )
    return float32_array


def _require_finite_embedding_source(
    values: object,
    *,
    label: str,
) -> np.ndarray:
    """Validate one embedding without discarding precision before normalization.

    Exported coordinates are float32, but the normalization extents have to be
    computed at the input's own precision. Casting first collapses any axis
    whose spread is finer than float32 resolution at its magnitude -- a UMAP
    around 1e8 loses unit spacing entirely and the axis is exported as a
    constant, silently flattening the embedding to a line. The R exporter
    normalizes from the source doubles for this reason.

    The float32 range check is kept: the centre is exported as float32, so a
    magnitude that cannot survive the conversion is still rejected here.
    """
    dense = cast(sparse.spmatrix, values).toarray() if sparse.issparse(values) else values
    array = np.asarray(dense)
    if array.dtype.kind not in {"i", "u", "f"}:
        raise TypeError(f"{label} must contain real numeric values.")
    if not np.isfinite(array).all():
        raise ValueError(
            f"{label} must contain only finite values. "
            + describe_non_finite(array, label)
        )
    with np.errstate(over="ignore", invalid="ignore"):
        as_float32 = array.astype(np.float32, copy=False)
        if not np.isfinite(as_float32).all():
            raise ValueError(
                f"{label} contains values outside the finite float32 range. "
                + describe_non_finite(array, label)
            )
        # Underflow is the other end of the same range, and it is the quiet one:
        # a nonzero value below the smallest float32 subnormal becomes 0.0,
        # which is finite, so an overflow-only check publishes a column whose
        # every value has been replaced by zero. The R writer refuses exactly
        # this, and one export format cannot mean two things.
        if np.any((array != 0) & (as_float32 == 0)):
            raise ValueError(
                f"{label} contains values outside the finite float32 range. "
                + describe_non_finite(array, label)
            )
    return array.astype(np.float64, copy=False)


def _require_continuous_obs_values(
    values: object,
    *,
    key: str,
    n_cells: int,
) -> np.ndarray:
    """Require one exact finite float32 observation vector."""
    array = _require_finite_float32_array(
        values,
        label=f"Continuous obs field {key!r}",
    )
    if array.ndim != 1 or array.shape[0] != n_cells:
        raise ValueError(
            f"Continuous obs field {key!r} must have shape ({n_cells},), got {array.shape}."
        )
    return array


def _normalize_finite_float32_embedding(
    embedding: np.ndarray,
    *,
    label: str,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Normalize one validated embedding with a required nonzero range.

    The normalized coordinates stay in float64. They are published as float32,
    but the rounding belongs to the single place that writes the bytes:
    rounding here would also round the centroids that are measured from these
    coordinates and published as JSON doubles. cellucid-r normalizes in double
    and rounds at the write too, so this is what keeps the two writers equal.
    """
    if embedding.ndim != 2 or embedding.shape[0] == 0:
        raise ValueError(f"{label} must be a non-empty 2D array.")
    working = embedding.astype(np.float64)
    axis_mins = working.min(axis=0)
    axis_maxs = working.max(axis=0)
    max_range = float((axis_maxs - axis_mins).max())
    if max_range <= 0:
        raise ValueError(f"{label} has no coordinate variation and cannot be normalized.")
    center = (axis_mins + axis_maxs) / 2
    scale_factor = 2.0 / max_range
    normalized = (working - center) * scale_factor
    with np.errstate(over="ignore", invalid="ignore"):
        if not np.isfinite(normalized.astype(np.float32)).all():
            raise ValueError(f"{label} normalization produced non-finite coordinates.")
    return normalized, center, scale_factor, max_range
