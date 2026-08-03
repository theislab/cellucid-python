"""One rule governs every dimensional quantity Cellucid accepts.

The dimension is always *declared* -- by the argument name ``X_umap_<n>d`` or by
the ``<field>_umap_<n>d`` key -- and never inferred from an array's width. A
declared 1D quantity, however, accepts the plain ``(n_cells,)`` vector that is
its natural spelling, because AnnData itself performs that reshape when such an
array is put in ``obsm``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import prepare
from cellucid.vector_fields import validate_vector_fields

N_CELLS = 8


def _obs() -> pd.DataFrame:
    return pd.DataFrame(
        {"leiden": pd.Categorical(["0", "1"] * (N_CELLS // 2))},
        index=[f"cell-{index}" for index in range(N_CELLS)],
    )


def _points_1d() -> np.ndarray:
    return np.linspace(-1.0, 1.0, N_CELLS, dtype=np.float32)


@pytest.mark.parametrize("key", ["velocity_umap", "velocity", "velocity_umap_4d", "velocity_2d"])
def test_a_vector_field_key_without_an_exact_dimension_suffix_is_rejected(key: str) -> None:
    with pytest.raises(ValueError, match=r"<field>_umap_<1\|2\|3>d"):
        validate_vector_fields(
            {key: np.zeros((N_CELLS, 2), dtype=np.float32)},
            n_cells=N_CELLS,
            available_dimensions=[2],
        )


def test_a_declared_1d_vector_field_accepts_a_plain_vector_as_its_single_column() -> None:
    validated = validate_vector_fields(
        {"velocity_umap_1d": _points_1d()},
        n_cells=N_CELLS,
        available_dimensions=[1],
    )
    vectors = validated.fields["velocity_umap"][1]
    assert vectors.shape == (N_CELLS, 1)
    # A validated vector is an input to the export, not a payload: it keeps the
    # caller's precision so the export rounds to float32 exactly once.
    assert vectors.dtype == np.float64
    assert np.array_equal(vectors[:, 0], _points_1d())


@pytest.mark.parametrize("declared", [2, 3])
def test_a_plain_vector_never_satisfies_a_higher_declared_dimension(declared: int) -> None:
    with pytest.raises(ValueError, match="must be a 2D array"):
        validate_vector_fields(
            {f"velocity_umap_{declared}d": _points_1d()},
            n_cells=N_CELLS,
            available_dimensions=[declared],
        )


def test_a_declared_1d_vector_field_still_requires_one_value_per_cell() -> None:
    with pytest.raises(ValueError, match="must have exactly"):
        validate_vector_fields(
            {"velocity_umap_1d": np.zeros(N_CELLS + 1, dtype=np.float32)},
            n_cells=N_CELLS,
            available_dimensions=[1],
        )


def test_prepare_accepts_a_plain_vector_for_the_declared_1d_embedding(tmp_path: Path) -> None:
    out_dir = tmp_path / "exports"
    prepare(
        latent_space=np.zeros((N_CELLS, 2), dtype=np.float32),
        obs=_obs(),
        out_dir=out_dir,
        obs_categorical_dtype="uint8",
        dataset_name="Declared dimensions",
        dataset_id="declared-dimensions",
        created_at="2026-01-01T00:00:00Z",
        centroid_min_points=1,
        X_umap_1d=_points_1d(),
        vector_fields={"velocity_umap_1d": _points_1d()},
    )

    points = np.frombuffer((out_dir / "points_1d.bin").read_bytes(), dtype="<f4")
    assert points.size == N_CELLS
    vectors = np.frombuffer(
        (out_dir / "vectors" / "0_1d.bin").read_bytes(),
        dtype="<f4",
    )
    assert vectors.size == N_CELLS


def test_prepare_still_refuses_a_plain_vector_for_a_higher_embedding(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="X_umap_2d must be a 2D array"):
        prepare(
            latent_space=np.zeros((N_CELLS, 2), dtype=np.float32),
            obs=_obs(),
            out_dir=tmp_path / "exports",
            obs_categorical_dtype="uint8",
            dataset_name="Declared dimensions",
            dataset_id="declared-dimensions",
            X_umap_2d=_points_1d(),
        )


def test_the_anndata_path_and_prepare_agree_on_a_plain_1d_declaration() -> None:
    adata = AnnData(
        X=np.zeros((N_CELLS, 1), dtype=np.float32),
        obs=_obs(),
        var=pd.DataFrame(index=["GeneA"]),
    )
    # AnnData stores the plain vector as the (n_cells, 1) column it declares,
    # which is exactly what prepare now accepts for the same declaration.
    adata.obsm["X_umap_1d"] = _points_1d()
    adata.obsm["velocity_umap_1d"] = _points_1d()
    assert adata.obsm["X_umap_1d"].shape == (N_CELLS, 1)

    adapter = AnnDataAdapter(
        adata,
        dataset_name="Declared dimensions",
        dataset_id="declared-dimensions",
        centroid_min_points=1,
        serve_vector_fields=True,
    )
    try:
        assert adapter.get_embedding(1).shape == (N_CELLS, 1)
        metadata = adapter.get_dataset_identity()["vector_fields"]
        assert list(metadata["fields"]) == ["velocity_umap"]
    finally:
        adapter.close()


def test_an_unsuffixed_obsm_vector_key_is_never_discovered_as_a_vector_field() -> None:
    adata = AnnData(
        X=np.zeros((N_CELLS, 1), dtype=np.float32),
        obs=_obs(),
        var=pd.DataFrame(index=["GeneA"]),
        obsm={
            "X_umap_2d": np.tile(_points_1d()[:, None], (1, 2)),
            "velocity_umap": np.zeros((N_CELLS, 2), dtype=np.float32),
        },
    )
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Declared dimensions",
        dataset_id="declared-dimensions",
        centroid_min_points=1,
        serve_vector_fields=True,
    )
    try:
        assert "vector_fields" not in adapter.get_dataset_identity()
    finally:
        adapter.close()
