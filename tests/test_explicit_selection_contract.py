from __future__ import annotations

import inspect
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import prepare
from cellucid.vector_fields import (
    add_transition_drift_to_obsm,
    compute_transition_drift,
)


def _embedding() -> np.ndarray:
    return np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        dtype=np.float32,
    )


def _adata(*, explicit: bool, generic: bool) -> ad.AnnData:
    adata = ad.AnnData(
        X=np.empty((3, 0), dtype=np.float32),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "B", "A"])},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(index=pd.Index([], dtype=str)),
    )
    if explicit:
        adata.obsm["X_umap_2d"] = _embedding()
    if generic:
        adata.obsm["X_umap"] = _embedding()
    return adata


def test_generic_umap_is_not_an_embedding_declaration() -> None:
    with pytest.raises(ValueError, match="X_umap_1d.*X_umap_2d.*X_umap_3d"):
        AnnDataAdapter(
            _adata(explicit=False, generic=True),
            dataset_name="Generic UMAP",
            dataset_id="generic-umap",
        )

    adapter = AnnDataAdapter(
        _adata(explicit=True, generic=True),
        dataset_name="Explicit UMAP",
        dataset_id="explicit-umap",
    )
    identity = adapter.get_dataset_identity()
    assert identity["embeddings"]["available_dimensions"] == [2]
    assert "umap_resolution" not in identity["embeddings"]


def test_transition_helper_requires_an_explicit_umap_key() -> None:
    with pytest.raises(ValueError, match="X_umap_1d.*X_umap_2d.*X_umap_3d"):
        add_transition_drift_to_obsm(
            _adata(explicit=False, generic=True),
            np.eye(3, dtype=np.float32),
            normalize_rows=False,
        )


def test_transition_row_normalization_is_an_explicit_required_choice() -> None:
    for public_function in (compute_transition_drift, add_transition_drift_to_obsm):
        parameter = inspect.signature(public_function).parameters["normalize_rows"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is inspect.Parameter.empty

    with pytest.raises(TypeError, match="normalize_rows"):
        compute_transition_drift(np.eye(3, dtype=np.float32), _embedding())

    with pytest.raises(TypeError, match="normalize_rows"):
        add_transition_drift_to_obsm(
            _adata(explicit=True, generic=False),
            np.eye(3, dtype=np.float32),
        )


def _prepare(
    out_dir: Path,
    *,
    obs_categorical_dtype: object = inspect.Parameter.empty,
) -> None:
    kwargs: dict[str, object] = {
        "latent_space": _embedding(),
        "obs": pd.DataFrame(
            {"cluster": pd.Categorical(["A", "B", "A"])},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        "X_umap_2d": _embedding(),
        "out_dir": out_dir,
        "dataset_name": "Categorical dtype",
        "dataset_id": "categorical-dtype",
        "centroid_min_points": 1,
    }
    if obs_categorical_dtype is not inspect.Parameter.empty:
        kwargs["obs_categorical_dtype"] = obs_categorical_dtype
    prepare(**kwargs)


def test_prepare_requires_one_exact_categorical_storage_dtype(tmp_path: Path) -> None:
    parameter = inspect.signature(prepare).parameters["obs_categorical_dtype"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is inspect.Parameter.empty

    with pytest.raises(TypeError, match="obs_categorical_dtype"):
        _prepare(tmp_path / "missing")
    assert not (tmp_path / "missing").exists()

    with pytest.raises(ValueError, match="obs_categorical_dtype"):
        _prepare(tmp_path / "automatic", obs_categorical_dtype="auto")
    assert not (tmp_path / "automatic").exists()

    _prepare(tmp_path / "uint16", obs_categorical_dtype="uint16")
    manifest = json.loads(
        (tmp_path / "uint16" / "obs_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["_categoricalFields"][0][2] == "uint16"
    identity = json.loads(
        (tmp_path / "uint16" / "dataset_identity.json").read_text(encoding="utf-8")
    )
    assert identity["export_settings"]["obs_categorical_dtype"] == "uint16"
