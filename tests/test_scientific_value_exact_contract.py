from __future__ import annotations

from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest

from cellucid import prepare

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


#: A NaN or an infinity in a row of a per-cell array is a cell the array does
#: not place (web, Python and R alike): it is written as given, not refused.
UNPLACED_VECTORS = [param for param in INVALID_VECTORS if param.id in {"nan", "infinity"}]
REFUSED_LATENT_VECTORS = [param for param in INVALID_VECTORS if param.id not in {"nan", "infinity"}]


@pytest.mark.parametrize("invalid", REFUSED_LATENT_VECTORS)
def test_prepare_rejects_invalid_latent_values_atomically(
    tmp_path: Path,
    invalid: np.ndarray,
) -> None:
    output = tmp_path / "must-not-exist.h5ad"
    latent = np.column_stack([invalid, np.array([0.0, 1.0, 2.0])])

    with pytest.raises((TypeError, ValueError), match=r"Embedding 'X_latent'.*(real|finite)"):
        prepare(
            output,
            format="h5ad",
            obs=_obs(),
            embeddings={"X_umap": _embedding(), "X_latent": latent},
            latent_key="X_latent",
            dataset_name="Invalid latent",
            dataset_id="invalid-latent",
        )

    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("invalid", UNPLACED_VECTORS)
def test_a_latent_row_with_nan_or_infinity_is_an_unplaced_cell(
    tmp_path: Path,
    invalid: np.ndarray,
) -> None:
    latent = np.column_stack([invalid, np.array([0.0, 1.0, 2.0])])

    written = anndata.read_h5ad(
        prepare(
            tmp_path / "unplaced.h5ad",
            format="h5ad",
            obs=_obs(),
            embeddings={"X_umap": _embedding(), "X_latent": latent},
            latent_key="X_latent",
            centroid_min_points=1,
            dataset_name="Unplaced latent",
            dataset_id="unplaced-latent",
        )
    )

    # prepare() may store the cells in another order: compare them by name.
    rows = written.obs_names.get_indexer(_obs().index)
    stored = np.asarray(written.obsm["X_latent"])[rows]
    assert np.array_equal(stored, latent.astype(np.float32), equal_nan=True)
    quantiles = np.asarray(written.uns["cellucid"]["outlierQuantiles"]["fields"]["group"])[rows]
    assert np.isnan(quantiles[1]) and not np.isnan(quantiles[[0, 2]]).any()


@pytest.mark.parametrize("invalid", INVALID_VECTORS)
def test_prepare_rejects_invalid_gene_values_atomically(
    tmp_path: Path,
    invalid: np.ndarray,
) -> None:
    output = tmp_path / "must-not-exist.h5ad"

    with pytest.raises((TypeError, ValueError), match=r"gene_expression.*(real|finite)"):
        prepare(
            output,
            format="h5ad",
            obs=_obs(),
            var=pd.DataFrame(index=["GeneA"]),
            gene_expression=invalid.reshape(3, 1),
            embeddings={"X_umap": _embedding()},
            dataset_name="Invalid expression",
            dataset_id="invalid-expression",
        )

    assert list(tmp_path.iterdir()) == []
