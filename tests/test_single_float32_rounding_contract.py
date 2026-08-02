"""Every float32 payload is rounded exactly once, at the write.

A coordinate and a vector component are both float32 on disk, but the value
that is rounded must be the scaled one, computed from the caller's own numbers.
Rounding the input first and the product second rounds twice, and every
component whose two roundings disagree lands one ULP away from the correctly
rounded value.

That is not only an accuracy question. cellucid-r scales the caller's doubles
and rounds once, so a second rounding here is exactly what made the two writers
publish different bytes from one input: on a 600-cell fixture, 25-28% of every
``vectors/*`` component differed by one ULP while ``points_*`` -- rounded once
on both sides -- was byte-identical. These tests pin the single rounding on
every payload the two writers share.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import prepare

N_CELLS = 96


def _source_arrays() -> tuple[np.ndarray, np.ndarray]:
    """Return one embedding and one vector field that expose a double rounding.

    The magnitudes are irrational multiples that keep every float64 mantissa
    bit occupied, so rounding to float32 before scaling changes the published
    value for a large fraction of the components rather than for none.
    """
    rng = np.random.default_rng(20260730)
    embedding = rng.normal(loc=3.7182818284590452, scale=4.1231056256176605, size=(N_CELLS, 2))
    embedding += 1e-9 * rng.normal(size=(N_CELLS, 2))
    vectors = rng.normal(scale=0.3183098861837907, size=(N_CELLS, 2))
    vectors += 1e-9 * rng.normal(size=(N_CELLS, 2))
    return embedding, vectors


def _normalization(embedding: np.ndarray) -> tuple[np.ndarray, float]:
    axis_mins = embedding.min(axis=0)
    axis_maxs = embedding.max(axis=0)
    center = (axis_mins + axis_maxs) / 2
    scale_factor = 2.0 / float((axis_maxs - axis_mins).max())
    return center, scale_factor


def _obs(categories: int = 4) -> pd.DataFrame:
    codes = np.arange(N_CELLS) % categories
    labels = [f"Group {index}" for index in range(categories)]
    return pd.DataFrame(
        {"cell_type": pd.Categorical.from_codes(codes, categories=labels)},
        index=[f"cell-{index}" for index in range(N_CELLS)],
    )


def _prepare(out_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    embedding, vectors = _source_arrays()
    prepare(
        latent_space=embedding.copy(),
        obs=_obs(),
        X_umap_2d=embedding,
        out_dir=out_dir,
        dataset_name="Single rounding",
        dataset_id="single-rounding",
        obs_categorical_dtype="uint16",
        centroid_min_points=4,
        vector_fields={"velocity_umap_2d": vectors},
    )
    return embedding, vectors


def test_a_vector_payload_is_the_scaled_value_rounded_once(tmp_path: Path) -> None:
    embedding, vectors = _prepare(tmp_path)
    _center, scale_factor = _normalization(embedding)

    published = np.fromfile(tmp_path / "vectors" / "0_2d.bin", dtype="<f4")
    rounded_once = (vectors * scale_factor).astype(np.float32).reshape(-1)
    rounded_twice = (
        (vectors.astype(np.float32).astype(np.float64) * scale_factor)
        .astype(np.float32)
        .reshape(-1)
    )

    assert np.array_equal(published, rounded_once)
    # The fixture is only evidence if the two models actually disagree on it.
    differing = int(np.count_nonzero(rounded_once != rounded_twice))
    assert differing > 0
    assert not np.array_equal(published, rounded_twice)


def test_a_points_payload_is_the_normalized_value_rounded_once(tmp_path: Path) -> None:
    embedding, _vectors = _prepare(tmp_path)
    center, scale_factor = _normalization(embedding)

    published = np.fromfile(tmp_path / "points_2d.bin", dtype="<f4")
    rounded_once = ((embedding - center) * scale_factor).astype(np.float32).reshape(-1)

    assert np.array_equal(published, rounded_once)


def test_a_centroid_is_measured_from_unrounded_coordinates(tmp_path: Path) -> None:
    embedding, _vectors = _prepare(tmp_path)
    center, scale_factor = _normalization(embedding)
    coordinates = (embedding - center) * scale_factor

    manifest = json.loads((tmp_path / "obs_manifest.json").read_text(encoding="utf-8"))
    centroids = manifest["_categoricalFields"][0][5]["2"]
    assert len(centroids) == 4

    codes = np.arange(N_CELLS) % 4
    for code, centroid in enumerate(centroids):
        points = coordinates[codes == code, :]
        exact = points.mean(axis=0)
        distances = np.linalg.norm(points - exact, axis=1)
        threshold = float(np.quantile(distances, 0.95))
        inliers = points[distances <= threshold, :]
        expected = inliers.mean(axis=0) if inliers.shape[0] >= 4 else exact
        assert centroid["position"] == expected.tolist()

        # The same measurement taken from float32-rounded coordinates is a
        # different number, so the assertion above is not vacuous.
        rounded = coordinates.astype(np.float32)[codes == code, :].astype(np.float64)
        rounded_distances = np.linalg.norm(rounded - rounded.mean(axis=0), axis=1)
        rounded_threshold = float(np.quantile(rounded_distances, 0.95))
        rounded_inliers = rounded[rounded_distances <= rounded_threshold, :]
        assert centroid["position"] != rounded_inliers.mean(axis=0).tolist()


def test_the_served_payloads_carry_the_prepared_bytes(tmp_path: Path) -> None:
    # One AnnData served live and the same AnnData exported must publish the
    # same bytes, or the viewer shows different values for one dataset.
    embedding, vectors = _prepare(tmp_path)
    adata = ad.AnnData(
        X=np.empty((N_CELLS, 0), dtype=np.float32),
        obs=_obs(),
        var=pd.DataFrame(index=pd.Index([], dtype=str)),
    )
    adata.obsm["X_umap_2d"] = embedding
    adata.obsm["velocity_umap_2d"] = vectors
    adapter = AnnDataAdapter(
        adata,
        latent_key=None,
        centroid_min_points=4,
        dataset_name="Single rounding",
        dataset_id="single-rounding",
    )
    try:
        assert adapter.get_points_binary(2) == (tmp_path / "points_2d.bin").read_bytes()
        assert (
            adapter.get_vector_field_binary("velocity_umap", 2)
            == (tmp_path / "vectors" / "0_2d.bin").read_bytes()
        )
    finally:
        adapter.close()


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_a_vector_field_of_any_dimension_is_rounded_once(
    tmp_path: Path,
    dimension: int,
) -> None:
    rng = np.random.default_rng(11 + dimension)
    embedding = rng.normal(scale=2.7182818284590452, size=(N_CELLS, dimension))
    vectors = rng.normal(scale=0.5772156649015329, size=(N_CELLS, dimension))
    out_dir = tmp_path / f"{dimension}d"
    prepare(
        latent_space=embedding.copy(),
        obs=_obs(),
        out_dir=out_dir,
        dataset_name="Single rounding",
        dataset_id="single-rounding",
        obs_categorical_dtype="uint16",
        centroid_min_points=4,
        vector_fields={f"velocity_umap_{dimension}d": vectors},
        **{f"X_umap_{dimension}d": embedding},
    )
    _center, scale_factor = _normalization(embedding)
    published = np.fromfile(out_dir / "vectors" / f"0_{dimension}d.bin", dtype="<f4")
    assert np.array_equal(published, (vectors * scale_factor).astype(np.float32).reshape(-1))
