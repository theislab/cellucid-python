from __future__ import annotations

import inspect

import anndata as ad
import numpy as np
import pandas as pd
import pytest

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


def test_transition_helper_reads_the_basis_it_is_given() -> None:
    # The drift is computed over the basis named, and named for it.
    written = add_transition_drift_to_obsm(
        _adata(explicit=False, generic=True),
        np.eye(3, dtype=np.float32),
        normalize_rows=False,
    )
    assert written == "T_fwd_umap"
    written = add_transition_drift_to_obsm(
        _adata(explicit=True, generic=False),
        np.eye(3, dtype=np.float32),
        basis="umap_2d",
        normalize_rows=False,
    )
    assert written == "T_fwd_umap_2d"

    with pytest.raises(ValueError, match=r"obsm\['X_umap'\]"):
        add_transition_drift_to_obsm(
            _adata(explicit=False, generic=False),
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
