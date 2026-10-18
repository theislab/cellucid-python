"""One rule governs every dimensional quantity Cellucid accepts.

An embedding's dimension is the number of columns it is read from -- every
column of a 1- to 3-column array, or the ``columns`` an ``Embedding`` names --
and a vector field's dimension is the dimension of the embedding it is bound
to. Nothing is inferred from a name. A one-column quantity, however, accepts the
plain ``(n_cells,)`` vector that is its natural spelling, because AnnData itself
performs that reshape when such an array is put in ``obsm``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import Embedding, prepare
from cellucid.vector_fields import validate_vector_fields

N_CELLS = 8


def _obs() -> pd.DataFrame:
    return pd.DataFrame(
        {"leiden": pd.Categorical(["0", "1"] * (N_CELLS // 2))},
        index=[f"cell-{index}" for index in range(N_CELLS)],
    )


def _points_1d() -> np.ndarray:
    return np.linspace(-1.0, 1.0, N_CELLS, dtype=np.float32)


def test_a_vector_field_bound_to_an_embedding_that_does_not_exist_is_rejected() -> None:
    with pytest.raises(ValueError, match="not one of the exported embeddings"):
        validate_vector_fields(
            {"velocity": {"X_tsne": np.zeros((N_CELLS, 2), dtype=np.float32)}},
            n_cells=N_CELLS,
            embedding_dimensions={"X_umap": 2},
        )


def test_a_1d_bound_vector_field_accepts_a_plain_vector_as_its_single_column() -> None:
    validated = validate_vector_fields(
        {"velocity": {"X_line": _points_1d()}},
        n_cells=N_CELLS,
        embedding_dimensions={"X_line": 1},
    )
    vectors = validated.fields["velocity"]["X_line"]
    assert vectors.shape == (N_CELLS, 1)
    # A validated vector is an input to the export, not a payload: it keeps the
    # caller's precision so the export rounds to float32 exactly once.
    assert vectors.dtype == np.float64
    assert np.array_equal(vectors[:, 0], _points_1d())


@pytest.mark.parametrize("dimension", [2, 3])
def test_a_plain_vector_never_satisfies_a_higher_bound_dimension(dimension: int) -> None:
    with pytest.raises(ValueError, match="must be a 2D array"):
        validate_vector_fields(
            {"velocity": {"X_umap": _points_1d()}},
            n_cells=N_CELLS,
            embedding_dimensions={"X_umap": dimension},
        )


def test_a_1d_bound_vector_field_still_requires_one_value_per_cell() -> None:
    with pytest.raises(ValueError, match="must have exactly"):
        validate_vector_fields(
            {"velocity": {"X_line": np.zeros(N_CELLS + 1, dtype=np.float32)}},
            n_cells=N_CELLS,
            embedding_dimensions={"X_line": 1},
        )


def test_prepare_accepts_a_plain_vector_as_a_1d_embedding(tmp_path: Path) -> None:
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
        embeddings={"X_line": _points_1d()},
        vector_fields={"velocity": {"X_line": _points_1d()}},
    )

    points = np.frombuffer((out_dir / "embeddings" / "0_1d.bin").read_bytes(), dtype="<f4")
    assert points.size == N_CELLS
    vectors = np.frombuffer(
        (out_dir / "vectors" / "0_0.bin").read_bytes(),
        dtype="<f4",
    )
    assert vectors.size == N_CELLS


def test_prepare_refuses_a_plain_vector_read_as_a_higher_embedding(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="'X_umap'.*needs at least 2 columns; this array has 1"):
        prepare(
            latent_space=np.zeros((N_CELLS, 2), dtype=np.float32),
            obs=_obs(),
            out_dir=tmp_path / "exports",
            obs_categorical_dtype="uint8",
            dataset_name="Declared dimensions",
            dataset_id="declared-dimensions",
            embeddings={"X_umap": Embedding(_points_1d(), columns=[0, 1])},
        )
    assert not (tmp_path / "exports").exists()


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
        assert adapter.get_embedding("X_umap_1d").shape == (N_CELLS, 1)
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
