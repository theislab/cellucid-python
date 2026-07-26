from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import prepare

INVALID_VECTORS = [
    pytest.param(
        np.array([0.0, np.nan, 2.0], dtype=np.float64),
        id="nan",
    ),
    pytest.param(
        np.array([0.0, np.inf, 2.0], dtype=np.float64),
        id="infinity",
    ),
    pytest.param(
        np.array([0.0, np.finfo(np.float64).max, 2.0], dtype=np.float64),
        id="outside-float32",
    ),
    pytest.param(
        np.array([0.0 + 0.0j, 1.0 + 2.0j, 2.0 + 0.0j], dtype=np.complex128),
        id="complex",
    ),
]


def _embedding() -> np.ndarray:
    return np.array(
        [[0.0, 0.0], [1.0, 0.25], [0.25, 1.0]],
        dtype=np.float32,
    )


def _obs() -> pd.DataFrame:
    return pd.DataFrame(
        {"group": pd.Categorical(["A", "A", "B"])},
        index=["cell-1", "cell-2", "cell-3"],
    )


def _adata(
    *,
    expression: np.ndarray | None = None,
    latent: np.ndarray | None = None,
) -> ad.AnnData:
    values = (
        np.array([[1.0], [2.0], [3.0]], dtype=np.float32)
        if expression is None
        else expression.reshape(3, 1)
    )
    adata = ad.AnnData(
        X=values,
        obs=_obs(),
        var=pd.DataFrame(index=["GeneA"]),
    )
    adata.obsm["X_umap_2d"] = _embedding()
    adata.obsm["X_latent"] = (
        np.array([[0.0, 0.0], [1.0, 0.5], [2.0, 1.0]], dtype=np.float32)
        if latent is None
        else latent
    )
    return adata


@pytest.mark.parametrize("invalid", INVALID_VECTORS)
def test_prepare_rejects_invalid_latent_values_atomically(
    tmp_path: Path,
    invalid: np.ndarray,
) -> None:
    output = tmp_path / "must-not-exist"
    latent = np.column_stack([invalid, np.array([0.0, 1.0, 2.0])])

    with pytest.raises((TypeError, ValueError), match=r"latent_space.*real|finite|float32"):
        prepare(
            latent_space=latent,
            obs=_obs(),
            X_umap_2d=_embedding(),
            out_dir=output,
            dataset_name="Invalid latent",
            dataset_id="invalid-latent",
            obs_categorical_dtype="uint16",
            centroid_min_points=1,
        )

    assert not output.exists()


@pytest.mark.parametrize("invalid", INVALID_VECTORS)
def test_prepare_rejects_invalid_gene_values_atomically(
    tmp_path: Path,
    invalid: np.ndarray,
) -> None:
    output = tmp_path / "must-not-exist"

    with pytest.raises((TypeError, ValueError), match=r"GeneA.*real|finite|float32"):
        prepare(
            latent_space=_embedding(),
            obs=_obs(),
            var=pd.DataFrame(index=["GeneA"]),
            gene_expression=invalid.reshape(3, 1),
            X_umap_2d=_embedding(),
            out_dir=output,
            dataset_name="Invalid expression",
            dataset_id="invalid-expression",
            obs_categorical_dtype="uint16",
            centroid_min_points=1,
        )

    assert not output.exists()


@pytest.mark.parametrize("invalid", INVALID_VECTORS)
def test_direct_anndata_rejects_invalid_latent_before_quantile_publication(
    invalid: np.ndarray,
) -> None:
    latent = np.column_stack([invalid, np.array([0.0, 1.0, 2.0])])
    adapter = AnnDataAdapter(
        _adata(latent=latent),
        latent_key="X_latent",
        dataset_name="Invalid latent",
        dataset_id="invalid-latent",
        centroid_min_points=1,
    )

    with pytest.raises((TypeError, ValueError), match=r"X_latent.*real|finite|float32"):
        adapter.get_obs_outlier_quantiles("group")


@pytest.mark.parametrize("invalid", INVALID_VECTORS)
def test_direct_anndata_rejects_invalid_gene_before_binary_publication(
    invalid: np.ndarray,
) -> None:
    adapter = AnnDataAdapter(
        _adata(expression=invalid),
        dataset_name="Invalid expression",
        dataset_id="invalid-expression",
    )

    with pytest.raises((TypeError, ValueError), match=r"GeneA.*real|finite|float32"):
        adapter.get_gene_expression("GeneA")


def test_valid_latent_and_gene_values_preserve_exact_finite_float32_semantics() -> None:
    adapter = AnnDataAdapter(
        _adata(),
        latent_key="X_latent",
        dataset_name="Valid scientific values",
        dataset_id="valid-scientific-values",
        centroid_min_points=1,
    )

    expression = np.frombuffer(
        adapter.get_gene_expression("GeneA"),
        dtype="<f4",
    )
    quantiles = np.frombuffer(
        adapter.get_obs_outlier_quantiles("group"),
        dtype="<f4",
    )
    assert expression.tolist() == [1.0, 2.0, 3.0]
    assert np.isfinite(quantiles).all()
    assert ((quantiles >= 0.0) & (quantiles <= 1.0)).all()
