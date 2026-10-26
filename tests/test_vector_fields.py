import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from cellucid.data_mapping import parse_vector_field_key
from cellucid.vector_fields import (
    add_transition_drift_to_obsm,
    compute_transition_drift,
)


def test_compute_transition_drift_identity_dense():
    rng = np.random.default_rng(0)
    X = rng.standard_normal((5, 2), dtype=np.float32)
    T = np.eye(5, dtype=np.float32)
    drift = compute_transition_drift(T, X, normalize_rows=False)
    # Drift keeps the precision it was computed in; it is rounded to float32
    # once, when it is written.
    assert drift.dtype == np.float64
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


def _adata_with(**obsm: np.ndarray) -> ad.AnnData:
    adata = ad.AnnData(X=np.zeros((3, 0), dtype=np.float32))
    for key, values in obsm.items():
        adata.obsm[key] = values
    return adata


TRIANGLE_2D = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
TRIANGLE_3D = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float32)


def test_add_transition_drift_to_obsm_names_the_drift_for_its_basis():
    adata = _adata_with(X_umap=TRIANGLE_2D)
    key = add_transition_drift_to_obsm(adata, np.eye(3, dtype=np.float32), normalize_rows=False)
    assert key == "T_fwd_umap"
    assert adata.obsm[key].shape == (3, 2)
    assert np.allclose(adata.obsm[key], 0)


@pytest.mark.parametrize(
    ("basis", "embedding", "expected"),
    [
        ("umap_3d", TRIANGLE_3D, "T_fwd_umap_3d"),
        ("branches", TRIANGLE_3D, "T_fwd_branches"),
        ("pca", np.arange(30, dtype=np.float32).reshape(3, 10), "T_fwd_pca"),
    ],
)
def test_add_transition_drift_to_obsm_works_on_any_basis(basis, embedding, expected):
    adata = _adata_with(X_umap=TRIANGLE_2D, **{f"X_{basis}": embedding})
    T = np.array([[0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=np.float32)
    key = add_transition_drift_to_obsm(adata, T, basis=basis, normalize_rows=True)
    assert key == expected
    # The drift is as wide as the embedding, in its coordinates.
    np.testing.assert_allclose(adata.obsm[key], T @ embedding - embedding)
    parsed = parse_vector_field_key(key, set(adata.obsm.keys()))
    assert parsed is not None
    assert (parsed.field, parsed.embedding_key) == ("T_fwd", f"X_{basis}")


def test_add_transition_drift_to_obsm_names_the_missing_embedding():
    adata = _adata_with(X_umap_2d=TRIANGLE_2D, X_umap_3d=TRIANGLE_3D)
    with pytest.raises(ValueError, match=r"obsm\['X_umap'\].*X_umap_2d.*X_umap_3d"):
        add_transition_drift_to_obsm(adata, np.eye(3, dtype=np.float32), normalize_rows=False)


def test_add_transition_drift_to_obsm_refuses_a_key_that_reads_as_another_field():
    adata = _adata_with(X_umap=TRIANGLE_2D, X_fwd_umap=TRIANGLE_2D)
    with pytest.raises(ValueError, match=r"reads as the field 'T' on obsm\['X_fwd_umap'\]"):
        add_transition_drift_to_obsm(adata, np.eye(3, dtype=np.float32), normalize_rows=False)
    with pytest.raises(ValueError, match="starts with 'X_'"):
        add_transition_drift_to_obsm(
            adata, np.eye(3, dtype=np.float32), field_prefix="X", normalize_rows=False
        )
    assert set(adata.obsm.keys()) == {"X_umap", "X_fwd_umap"}


def test_add_transition_drift_to_obsm_keeps_an_existing_key_unless_told():
    adata = _adata_with(X_umap=TRIANGLE_2D, T_fwd_umap=np.ones((3, 2), dtype=np.float32))
    with pytest.raises(KeyError, match="overwrite=True"):
        add_transition_drift_to_obsm(adata, np.eye(3, dtype=np.float32), normalize_rows=False)
    key = add_transition_drift_to_obsm(
        adata, np.eye(3, dtype=np.float32), normalize_rows=False, overwrite=True
    )
    assert np.allclose(adata.obsm[key], 0)


def test_compute_transition_drift_takes_any_embedding_width():
    X = np.arange(15, dtype=np.float32).reshape(3, 5)
    drift = compute_transition_drift(np.eye(3, dtype=np.float32), X, normalize_rows=False)
    assert drift.shape == (3, 5)


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
