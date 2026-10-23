from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid import prepare


@pytest.mark.parametrize(
    ("value", "error_type"),
    [
        (True, TypeError),
        (False, TypeError),
        (1.0, TypeError),
        (np.int64(1), TypeError),
        ("1", TypeError),
        (0, ValueError),
        (-1, ValueError),
    ],
)
def test_prepare_rejects_invalid_centroid_minimum_before_target_mutation(
    tmp_path: Path,
    value: object,
    error_type: type[Exception],
) -> None:
    output = tmp_path / "generation.zarr"
    embedding = np.array(
        [[0.0, 0.0], [1.0, 0.25], [0.25, 1.0]],
        dtype=np.float32,
    )
    with pytest.raises(error_type, match="centroid_min_points"):
        prepare(
            output,
            format="zarr",
            obs=pd.DataFrame({"group": pd.Categorical(["A", "A", "B"])}),
            embeddings={"X_umap": embedding},
            dataset_name="Centroid contract",
            dataset_id="centroid-contract",
            centroid_min_points=value,  # type: ignore[arg-type]
        )
    assert not output.exists()
    assert list(tmp_path.iterdir()) == []


def test_prepare_accepts_one_positive_native_centroid_minimum(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation.zarr"
    embedding = np.array(
        [[0.0, 0.0], [1.0, 0.25], [0.25, 1.0]],
        dtype=np.float32,
    )
    prepare(
        output,
        format="zarr",
        obs=pd.DataFrame({"group": pd.Categorical(["A", "A", "B"])}),
        embeddings={"X_umap": embedding},
        dataset_name="Centroid contract",
        dataset_id="centroid-contract",
        centroid_min_points=1,
    )
    defaults = ad.read_zarr(output).uns["cellucid"]["defaults"]
    assert defaults["centroidMinPoints"] == 1
