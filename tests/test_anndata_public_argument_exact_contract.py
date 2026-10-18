from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import prepare


def _adapter() -> AnnDataAdapter:
    embedding = np.array(
        [[0.0, 0.0], [1.0, 0.25], [0.25, 1.0]],
        dtype=np.float32,
    )
    return AnnDataAdapter(
        AnnData(
            X=np.array(
                [[1.0, 3.0], [2.0, 4.0], [3.0, 5.0]],
                dtype=np.float32,
            ),
            obs=pd.DataFrame(
                {
                    "score": [0.25, 0.5, 0.75],
                    "group": pd.Categorical(["A", "A", "B"]),
                },
                index=["cell-1", "cell-2", "cell-3"],
            ),
            var=pd.DataFrame(index=["GeneA", "GeneB"]),
            obsm={
                "X_umap_2d": embedding,
                "X_latent": embedding.copy(),
                "velocity_umap_2d": np.array(
                    [[0.1, 0.0], [0.0, 0.1], [-0.1, 0.0]],
                    dtype=np.float32,
                ),
            },
            obsp={
                "connectivities": np.array(
                    [
                        [0.0, 0.5, 0.0],
                        [0.5, 0.0, 0.25],
                        [0.0, 0.25, 0.0],
                    ],
                    dtype=np.float64,
                )
            },
        ),
        latent_key="X_latent",
        dataset_name="Exact public arguments",
        dataset_id="exact-public-arguments",
        centroid_min_points=1,
        serve_vector_fields=True,
    )


@pytest.mark.parametrize("embedding_id", [True, 2, 2.0, np.str_("X_umap_2d"), None, b"X_umap_2d"])
@pytest.mark.parametrize("operation", ["embedding", "binary", "vector"])
def test_every_public_embedding_argument_requires_one_native_string(
    embedding_id: object,
    operation: str,
) -> None:
    adapter = _adapter()
    if isinstance(embedding_id, np.str_):
        # A numpy string is a str subclass; it names the served embedding.
        if operation == "embedding":
            assert adapter.get_embedding(embedding_id).shape == (3, 2)
        return
    with pytest.raises(TypeError, match="embedding_id.*string"):
        if operation == "embedding":
            adapter.get_embedding(embedding_id)  # type: ignore[arg-type]
        elif operation == "binary":
            adapter.get_embedding_binary(embedding_id)  # type: ignore[arg-type]
        else:
            adapter.get_vector_field_binary(
                "velocity_umap",
                embedding_id,  # type: ignore[arg-type]
            )


@pytest.mark.parametrize("embedding_id", ["", "X_umap", "X_latent", "X_umap_2d[1,0]"])
@pytest.mark.parametrize("operation", ["embedding", "binary", "vector"])
def test_every_public_embedding_argument_rejects_an_unserved_id(
    embedding_id: str,
    operation: str,
) -> None:
    adapter = _adapter()
    with pytest.raises(ValueError, match="embedding_id|not served|not bound"):
        if operation == "embedding":
            adapter.get_embedding(embedding_id)
        elif operation == "binary":
            adapter.get_embedding_binary(embedding_id)
        else:
            adapter.get_vector_field_binary("velocity_umap", embedding_id)


def _compressed_call(
    adapter: AnnDataAdapter,
    operation: str,
    compress: object,
) -> object:
    if operation == "points":
        return adapter.get_embedding_binary("X_umap_2d", compress=compress)  # type: ignore[arg-type]
    if operation == "vector":
        return adapter.get_vector_field_binary(
            "velocity_umap",
            "X_umap_2d",
            compress=compress,  # type: ignore[arg-type]
        )
    if operation == "continuous":
        return adapter.get_obs_continuous_values(
            "score",
            compress=compress,  # type: ignore[arg-type]
        )
    if operation == "categorical":
        return adapter.get_obs_categorical_codes(
            "group",
            compress=compress,  # type: ignore[arg-type]
        )
    if operation == "outliers":
        return adapter.get_obs_outlier_quantiles(
            "group",
            compress=compress,  # type: ignore[arg-type]
        )
    if operation == "gene":
        return adapter.get_gene_expression(
            "GeneA",
            compress=compress,  # type: ignore[arg-type]
        )
    if operation == "connectivity":
        return adapter.get_connectivity_edges(
            compress=compress,  # type: ignore[arg-type]
        )
    raise AssertionError(f"unknown test operation {operation}")


@pytest.mark.parametrize(
    "compress",
    [0, 1, "false", "true", None, np.bool_(False)],
)
@pytest.mark.parametrize(
    "operation",
    [
        "points",
        "vector",
        "continuous",
        "categorical",
        "outliers",
        "gene",
        "connectivity",
    ],
)
def test_every_public_binary_getter_requires_one_native_boolean(
    compress: object,
    operation: str,
) -> None:
    adapter = _adapter()
    with pytest.raises(TypeError, match="compress.*True or False"):
        _compressed_call(adapter, operation, compress)


@pytest.mark.parametrize(
    "operation",
    ["points", "vector", "continuous", "outliers", "gene"],
)
def test_exact_boolean_compression_preserves_raw_and_gzip_bytes(
    operation: str,
) -> None:
    adapter = _adapter()
    raw = _compressed_call(adapter, operation, False)
    compressed = _compressed_call(adapter, operation, True)
    assert isinstance(raw, bytes)
    assert isinstance(compressed, bytes)
    assert gzip.decompress(compressed) == raw


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
    output = tmp_path / "generation"
    embedding = np.array(
        [[0.0, 0.0], [1.0, 0.25], [0.25, 1.0]],
        dtype=np.float32,
    )
    with pytest.raises(error_type, match="centroid_min_points"):
        prepare(
            latent_space=embedding.copy(),
            obs=pd.DataFrame({"group": pd.Categorical(["A", "A", "B"])}),
            embeddings={"X_umap": embedding},
            out_dir=output,
            dataset_name="Centroid contract",
            dataset_id="centroid-contract",
            obs_categorical_dtype="uint16",
            centroid_min_points=value,  # type: ignore[arg-type]
        )
    assert not output.exists()
    assert list(tmp_path.iterdir()) == []


def test_prepare_accepts_one_positive_native_centroid_minimum(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation"
    embedding = np.array(
        [[0.0, 0.0], [1.0, 0.25], [0.25, 1.0]],
        dtype=np.float32,
    )
    prepare(
        latent_space=embedding.copy(),
        obs=pd.DataFrame({"group": pd.Categorical(["A", "A", "B"])}),
        embeddings={"X_umap": embedding},
        out_dir=output,
        dataset_name="Centroid contract",
        dataset_id="centroid-contract",
        obs_categorical_dtype="uint16",
        centroid_min_points=1,
    )
    assert (output / "dataset_identity.json").is_file()
