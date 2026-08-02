"""Embedding normalization must not discard precision before measuring extents.

Exported coordinates are float32, but the extents have to be computed at the
input's own precision. Casting first collapses any axis whose spread is finer
than float32 resolution at its magnitude, exporting that axis as a constant and
silently flattening the embedding to a line.
"""

from __future__ import annotations

import numpy as np
import pytest

from cellucid.prepare_data._arrays import (
    _normalize_finite_float32_embedding,
    _require_finite_embedding_source,
)


def _normalize(embedding: np.ndarray) -> np.ndarray:
    source = _require_finite_embedding_source(embedding, label="X_umap_2d")
    normalized, _center, _scale, _range = _normalize_finite_float32_embedding(
        source,
        label="X_umap_2d",
    )
    return normalized


def test_unit_spacing_survives_a_large_coordinate_magnitude() -> None:
    # float32 spacing at 1e8 is 8, so every one of these rounds to the same
    # value if the cast happens before the extents are measured.
    embedding = np.array(
        [[1e8, 0.0], [1e8 + 1.0, 1.0], [1e8 + 2.0, 2.0], [1e8 + 3.0, 3.0]],
        dtype=np.float64,
    )

    normalized = _normalize(embedding)

    assert len(set(normalized[:, 0].tolist())) == 4, (
        "the first axis collapsed to a constant: four distinct coordinates were "
        "exported as one, flattening the embedding to a line"
    )
    # Both axes carry the same spacing, so a correct normalization maps them
    # identically -- the points lie on the diagonal, not on a vertical line.
    assert normalized[:, 0].tolist() == normalized[:, 1].tolist()
    assert normalized[:, 0].min() == pytest.approx(-1.0)
    assert normalized[:, 0].max() == pytest.approx(1.0)


def test_normalized_coordinates_keep_the_input_precision() -> None:
    # The payload is float32, but the rounding belongs to the one place that
    # writes the bytes: rounding here would round the centroids measured from
    # these coordinates too, and cellucid-r rounds only at the write.
    embedding = np.array([[0.0, 0.0], [1.0, 2.0]], dtype=np.float64)
    source = _require_finite_embedding_source(embedding, label="X_umap_2d")
    normalized, center, _scale, _range = _normalize_finite_float32_embedding(
        source,
        label="X_umap_2d",
    )
    assert normalized.dtype == np.float64
    assert center.dtype == np.float64


@pytest.mark.parametrize(
    ("embedding", "message"),
    [
        (np.array([[1e300, 0.0], [1e301, 1.0]]), "outside the finite float32 range"),
        # Underflow is the other end of the same range, and it is the quiet
        # one: every nonzero value becomes 0.0, which is finite, so an
        # overflow-only check publishes a column that has lost its content.
        # The R writer refuses this; Python must agree, or one export format
        # has two meanings.
        (
            np.array([[1e-320, 0.0], [2e-320, 1.0]]),
            "outside the finite float32 range",
        ),
        (np.array([[np.nan, 0.0], [1.0, 1.0]]), "must contain only finite values"),
        (np.array([[np.inf, 0.0], [1.0, 1.0]]), "must contain only finite values"),
    ],
)
def test_validation_guards_still_reject(embedding: np.ndarray, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _require_finite_embedding_source(embedding, label="X_umap_2d")


def test_a_wholly_constant_embedding_is_still_rejected() -> None:
    source = _require_finite_embedding_source(np.ones((4, 2)), label="X_umap_2d")
    with pytest.raises(ValueError, match="no coordinate variation"):
        _normalize_finite_float32_embedding(source, label="X_umap_2d")


def test_a_float32_input_is_unchanged_by_the_wider_validation() -> None:
    # Callers who already hold float32 must get exactly what they had before.
    embedding = np.array([[0.0, 0.0], [1.5, 3.25], [2.5, 1.75]], dtype=np.float32)
    normalized = _normalize(embedding)
    direct, _c, _s, _r = _normalize_finite_float32_embedding(
        embedding.astype(np.float64),
        label="X_umap_2d",
    )
    assert normalized.tolist() == direct.tolist()


@pytest.mark.parametrize(
    ("values", "published"),
    [
        # Eight distinct inputs; float32 writes one distinct value. Every
        # nonzero value is replaced by zero, and zero is finite, so an
        # overflow-only check let this through.
        (np.array([1e-320 * (index + 1) for index in range(8)]), False),
        # A genuine float32 subnormal survives the conversion and must still be
        # accepted - the rule is about values that are lost, not small ones.
        (np.array([0.0] + [1e-40] * 7), True),
        # No nonzero value, so nothing is lost.
        (np.zeros(8), True),
        (np.arange(1.0, 9.0), True),
    ],
)
def test_continuous_obs_agrees_with_the_r_writer_on_float32_underflow(
    values: np.ndarray,
    published: bool,
) -> None:
    """The two writers must accept and refuse the same obs columns.

    `binary-payloads.R` refuses any nonzero value below 2^-149. Python checked
    only the overflow end, so one export format had two meanings.
    """
    from cellucid.prepare_data._arrays import _require_continuous_obs_values

    if published:
        result = _require_continuous_obs_values(values, key="tiny", n_cells=8)
        assert result.shape == (8,)
        return

    with pytest.raises(ValueError, match="outside the finite float32 range"):
        _require_continuous_obs_values(values, key="tiny", n_cells=8)
