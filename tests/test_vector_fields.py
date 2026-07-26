import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from cellucid.vector_fields import (
    add_transition_drift_to_obsm,
    compute_transition_drift,
)


def test_compute_transition_drift_identity_dense():
    rng = np.random.default_rng(0)
    X = rng.standard_normal((5, 2), dtype=np.float32)
    T = np.eye(5, dtype=np.float32)
    drift = compute_transition_drift(T, X, normalize_rows=False)
    assert drift.dtype == np.float32
    assert drift.shape == X.shape
    assert np.allclose(drift, 0)


def test_compute_transition_drift_row_normalization_dense():
    rng = np.random.default_rng(1)
    X = rng.standard_normal((5, 3), dtype=np.float32)
    T = 2.0 * np.eye(5, dtype=np.float32)  # not row-stochastic, but proportional to identity
    drift = compute_transition_drift(T, X, normalize_rows=True)
    assert np.allclose(drift, 0)


def test_compute_transition_drift_sparse():
    rng = np.random.default_rng(2)
    X = rng.standard_normal((5, 2), dtype=np.float32)
    T = sparse.eye(5, format="csr", dtype=np.float32) * 3.0
    drift = compute_transition_drift(T, X, normalize_rows=True)
    assert np.allclose(drift, 0, atol=1e-6)


def test_add_transition_drift_to_obsm_infers_dim_from_explicit_key():
    adata = ad.AnnData(X=np.zeros((3, 0), dtype=np.float32))
    adata.obsm["X_umap_3d"] = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float32
    )
    T = np.eye(3, dtype=np.float32)
    key = add_transition_drift_to_obsm(
        adata,
        T,
        basis="umap",
        field_prefix="T_fwd",
        normalize_rows=False,
    )
    assert key == "T_fwd_umap_3d"
    assert key in adata.obsm
    assert adata.obsm[key].shape == (3, 3)
    assert np.allclose(adata.obsm[key], 0)


def test_add_transition_drift_to_obsm_rejects_generic_umap_key():
    adata = ad.AnnData(X=np.zeros((3, 0), dtype=np.float32))
    adata.obsm["X_umap"] = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    T = np.eye(3, dtype=np.float32)
    with pytest.raises(ValueError, match="X_umap_1d.*X_umap_2d.*X_umap_3d"):
        add_transition_drift_to_obsm(
            adata,
            T,
            basis="umap",
            field_prefix="T_fwd",
            normalize_rows=False,
        )


def test_add_transition_drift_to_obsm_uses_explicit_key_and_ignores_generic_key():
    adata = ad.AnnData(X=np.zeros((3, 0), dtype=np.float32))
    adata.obsm["X_umap_2d"] = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    adata.obsm["X_umap"] = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    T = np.eye(3, dtype=np.float32)
    key = add_transition_drift_to_obsm(
        adata,
        T,
        basis="umap",
        field_prefix="T_fwd",
        normalize_rows=False,
    )
    assert key == "T_fwd_umap_2d"


def test_add_transition_drift_requires_dim_when_multiple_embeddings_exist():
    adata = ad.AnnData(X=np.zeros((3, 0), dtype=np.float32))
    adata.obsm["X_umap_2d"] = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        dtype=np.float32,
    )
    adata.obsm["X_umap_3d"] = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        dtype=np.float32,
    )

    with pytest.raises(ValueError, match="dim must be explicit"):
        add_transition_drift_to_obsm(
            adata,
            np.eye(3, dtype=np.float32),
            basis="umap",
            field_prefix="T_fwd",
            normalize_rows=False,
        )


def test_compute_transition_drift_rejects_zero_row_sum_instead_of_substitution():
    embedding = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        dtype=np.float32,
    )
    transition = np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float32,
    )

    with pytest.raises(ValueError, match="nonzero finite row sum"):
        compute_transition_drift(
            transition,
            embedding,
            normalize_rows=True,
        )
