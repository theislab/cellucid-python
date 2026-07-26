from __future__ import annotations

import inspect
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import prepare
from cellucid.vector_fields import add_transition_drift_to_obsm


def _embedding_2d() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        dtype=np.float32,
    )


def _embedding_3d() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=np.float32,
    )


def _vectors_2d(offset: float = 0.0) -> np.ndarray:
    return np.array(
        [
            [0.1 + offset, 0.2],
            [0.3, 0.4 + offset],
            [0.5, 0.6],
        ],
        dtype=np.float32,
    )


def _prepare_kwargs(
    out_dir: Path,
    vector_fields: dict[object, object],
) -> dict[str, object]:
    embedding = _embedding_2d()
    return {
        "latent_space": embedding.copy(),
        "obs": pd.DataFrame(index=["cell-1", "cell-2", "cell-3"]),
        "X_umap_2d": embedding,
        "out_dir": out_dir,
        "dataset_name": "Exact vector fields",
        "dataset_id": "exact-vector-fields",
        "obs_categorical_dtype": "uint16",
        "centroid_min_points": 1,
        "vector_fields": vector_fields,
    }


def _adata(vector_fields: dict[str, object]) -> ad.AnnData:
    adata = ad.AnnData(
        X=np.empty((3, 0), dtype=np.float32),
        obs=pd.DataFrame(index=["cell-1", "cell-2", "cell-3"]),
        var=pd.DataFrame(index=pd.Index([], dtype=str)),
    )
    adata.obsm["X_umap_2d"] = _embedding_2d()
    for key, value in vector_fields.items():
        adata.obsm[key] = value
    return adata


def _adapter(
    vector_fields: dict[str, object],
    *,
    vector_field_default: object = None,
) -> AnnDataAdapter:
    return AnnDataAdapter(
        _adata(vector_fields),
        dataset_name="Exact vector fields",
        dataset_id="exact-vector-fields",
        vector_field_default=vector_field_default,
    )


def test_vector_declarations_require_an_explicit_dimension_suffix(
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, {"velocity_umap": _vectors_2d()})

    with pytest.raises(ValueError, match=r"<field>_umap_<1\|2\|3>d"):
        prepare(**kwargs)
    assert not out_dir.exists()

    adapter = _adapter({"velocity_umap": _vectors_2d()})
    assert "vector_fields" not in adapter.get_dataset_identity()


def test_transition_drift_exposes_no_unsuffixed_output_switch() -> None:
    assert "explicit_dim_suffix" not in inspect.signature(add_transition_drift_to_obsm).parameters


def test_prepare_single_vector_field_is_unambiguous_and_metadata_is_exact(
    tmp_path: Path,
) -> None:
    vectors = _vectors_2d()
    prepare(**_prepare_kwargs(tmp_path, {"velocity_umap_2d": vectors}))

    identity = json.loads((tmp_path / "dataset_identity.json").read_text(encoding="utf-8"))
    assert identity["vector_fields"] == {
        "default_field": "velocity_umap",
        "fields": {
            "velocity_umap": {
                "label": "velocity_umap",
                "basis": "umap",
                "available_dimensions": [2],
                "default_dimension": 2,
                "files": {"2d": "vectors/velocity_umap_2d.bin"},
            }
        },
    }
    exported = np.fromfile(
        tmp_path / "vectors" / "velocity_umap_2d.bin",
        dtype=np.float32,
    ).reshape(3, 2)
    np.testing.assert_allclose(exported, vectors * 2.0)


def test_vector_field_default_dimension_is_highest_available_dimension(
    tmp_path: Path,
) -> None:
    vectors = {
        "velocity_umap_2d": _vectors_2d(),
        "velocity_umap_3d": _embedding_3d(),
    }
    kwargs = _prepare_kwargs(tmp_path / "export", vectors)
    kwargs["X_umap_3d"] = _embedding_3d()
    prepare(**kwargs)
    prepared_identity = json.loads(
        (tmp_path / "export" / "dataset_identity.json").read_text(encoding="utf-8")
    )
    prepared_field = prepared_identity["vector_fields"]["fields"]["velocity_umap"]
    assert prepared_field["available_dimensions"] == [2, 3]
    assert prepared_field["default_dimension"] == 3

    adata = _adata(vectors)
    adata.obsm["X_umap_3d"] = _embedding_3d()
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Exact vector dimensions",
        dataset_id="exact-vector-dimensions",
    )
    direct_field = adapter.get_dataset_identity()["vector_fields"]["fields"]["velocity_umap"]
    assert direct_field["available_dimensions"] == [2, 3]
    assert direct_field["default_dimension"] == 3


def test_dataset_default_dimension_is_highest_declared_dimension(
    tmp_path: Path,
) -> None:
    embedding_1d = np.array([[0.0], [1.0], [2.0]], dtype=np.float32)
    kwargs = _prepare_kwargs(tmp_path / "export", {})
    kwargs["X_umap_1d"] = embedding_1d
    kwargs["X_umap_3d"] = _embedding_3d()
    prepare(**kwargs)
    prepared_identity = json.loads(
        (tmp_path / "export" / "dataset_identity.json").read_text(encoding="utf-8")
    )
    assert prepared_identity["embeddings"]["available_dimensions"] == [1, 2, 3]
    assert prepared_identity["embeddings"]["default_dimension"] == 3

    adata = _adata({})
    adata.obsm["X_umap_1d"] = embedding_1d
    adata.obsm["X_umap_3d"] = _embedding_3d()
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Exact dataset dimensions",
        dataset_id="exact-dataset-dimensions",
    )
    direct_embeddings = adapter.get_dataset_identity()["embeddings"]
    assert direct_embeddings["available_dimensions"] == [1, 2, 3]
    assert direct_embeddings["default_dimension"] == 3


def test_sparse_vector_values_follow_the_same_exact_contract(
    tmp_path: Path,
) -> None:
    vectors = _vectors_2d()
    sparse_vectors = sparse.csr_matrix(vectors)
    prepare(
        **_prepare_kwargs(
            tmp_path / "export",
            {"velocity_umap_2d": sparse_vectors},
        )
    )
    prepared = np.fromfile(
        tmp_path / "export" / "vectors" / "velocity_umap_2d.bin",
        dtype=np.float32,
    ).reshape(3, 2)
    np.testing.assert_allclose(prepared, vectors * 2.0)

    adapter = _adapter({"velocity_umap_2d": sparse_vectors})
    direct = np.frombuffer(
        adapter.get_vector_field_binary("velocity_umap", 2),
        dtype=np.float32,
    ).reshape(3, 2)
    np.testing.assert_allclose(direct, vectors * 2.0)


def test_prepare_multiple_vector_fields_require_an_exact_explicit_default(
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "export"
    fields = {
        "velocity_umap_2d": _vectors_2d(),
        "drift_umap_2d": _vectors_2d(0.1),
    }
    kwargs = _prepare_kwargs(out_dir, fields)

    with pytest.raises(ValueError, match="vector_field_default"):
        prepare(**kwargs)
    assert not out_dir.exists()

    kwargs["vector_field_default"] = "drift_umap"
    prepare(**kwargs)
    identity = json.loads((out_dir / "dataset_identity.json").read_text(encoding="utf-8"))
    assert identity["vector_fields"]["default_field"] == "drift_umap"


def test_vector_field_ids_must_be_case_insensitively_unique(
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    fields = {
        "Velocity_umap_2d": _vectors_2d(),
        "velocity_umap_2d": _vectors_2d(0.1),
    }
    kwargs = _prepare_kwargs(out_dir, fields)
    kwargs["vector_field_default"] = "Velocity_umap"

    with pytest.raises(ValueError, match="case-insensitively"):
        prepare(**kwargs)
    assert not out_dir.exists()

    with pytest.raises(ValueError, match="case-insensitively"):
        _adapter(fields, vector_field_default="Velocity_umap")


def test_vector_default_without_vector_fields_is_rejected_before_mutation(
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, {})
    kwargs["vector_field_default"] = "velocity_umap"

    with pytest.raises(ValueError, match="no vector fields"):
        prepare(**kwargs)
    assert not out_dir.exists()

    with pytest.raises(ValueError, match="no vector fields"):
        _adapter({}, vector_field_default="velocity_umap")


@pytest.mark.parametrize(
    "invalid_default",
    ["", "missing_umap", 7],
)
def test_prepare_rejects_invalid_vector_default_before_output_mutation(
    tmp_path: Path,
    invalid_default: object,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(
        out_dir,
        {"velocity_umap_2d": _vectors_2d()},
    )
    kwargs["vector_field_default"] = invalid_default

    with pytest.raises((TypeError, ValueError), match="vector_field_default"):
        prepare(**kwargs)
    assert not out_dir.exists()


@pytest.mark.parametrize(
    "invalid_key",
    [
        "",
        7,
        "velocity",
        "velocity_pca_2d",
        "velocity_umap_4d",
        "velocity_umap_2D",
        "velocity_umap_extra",
    ],
)
def test_prepare_rejects_unknown_or_non_string_vector_keys_before_mutation(
    tmp_path: Path,
    invalid_key: object,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, {invalid_key: _vectors_2d()})

    with pytest.raises((TypeError, ValueError), match="vector field key"):
        prepare(**kwargs)
    assert not out_dir.exists()


@pytest.mark.parametrize(
    ("key", "value", "message"),
    [
        ("velocity_umap_2d", None, "must not be None"),
        ("velocity_umap_2d", np.arange(3, dtype=np.float32), "2D array"),
        (
            "velocity_umap_2d",
            np.ones((2, 2), dtype=np.float32),
            "exactly 3 rows",
        ),
        (
            "velocity_umap_2d",
            np.ones((3, 3), dtype=np.float32),
            "exactly 2 columns",
        ),
        (
            "velocity_umap_2d",
            np.ones((3, 4), dtype=np.float32),
            "exactly 2 columns",
        ),
        (
            "velocity_umap_2d",
            np.ones((3, 2), dtype=np.complex64),
            "real numeric",
        ),
        (
            "velocity_umap_2d",
            np.ones((3, 2), dtype=np.bool_),
            "real numeric",
        ),
        (
            "velocity_umap_2d",
            np.array([[0.0, np.nan], [0.0, 0.0], [0.0, 0.0]]),
            "finite",
        ),
        (
            "velocity_umap_2d",
            np.full((3, 2), 1e100, dtype=np.float64),
            "float32 range",
        ),
    ],
)
def test_prepare_rejects_invalid_vector_values_before_output_mutation(
    tmp_path: Path,
    key: str,
    value: object,
    message: str,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, {key: value})

    with pytest.raises((TypeError, ValueError), match=message):
        prepare(**kwargs)
    assert not out_dir.exists()


def test_prepare_rejects_vector_dimension_without_matching_embedding(
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(
        out_dir,
        {"velocity_umap_3d": _embedding_3d()},
    )

    with pytest.raises(ValueError, match="matching 3D embedding"):
        prepare(**kwargs)
    assert not out_dir.exists()


def test_anndata_single_vector_field_metadata_and_binary_are_exact() -> None:
    vectors = _vectors_2d()
    adapter = _adapter({"velocity_umap_2d": vectors})

    identity = adapter.get_dataset_identity()
    assert identity["vector_fields"] == {
        "default_field": "velocity_umap",
        "fields": {
            "velocity_umap": {
                "label": "velocity_umap",
                "basis": "umap",
                "available_dimensions": [2],
                "default_dimension": 2,
                "files": {"2d": "vectors/velocity_umap_2d.bin"},
            }
        },
    }
    exported = np.frombuffer(
        adapter.get_vector_field_binary("velocity_umap", 2),
        dtype=np.float32,
    ).reshape(3, 2)
    np.testing.assert_allclose(exported, vectors * 2.0)


def test_anndata_multiple_vector_fields_require_explicit_default() -> None:
    fields = {
        "velocity_umap_2d": _vectors_2d(),
        "drift_umap_2d": _vectors_2d(0.1),
    }
    with pytest.raises(ValueError, match="vector_field_default"):
        _adapter(fields)

    adapter = _adapter(fields, vector_field_default="drift_umap")
    assert adapter.get_dataset_identity()["vector_fields"]["default_field"] == "drift_umap"


@pytest.mark.parametrize(
    "key",
    [
        "X_umap_harmony",
        "my_umap_history",
        "pca_umap_projection",
        "velocity_umap_extra",
        "velocity_umap_4d",
    ],
)
def test_anndata_ignores_obsm_keys_outside_the_exact_vector_declaration_grammar(
    key: str,
) -> None:
    adapter = _adapter({key: np.ones((3, 4), dtype=np.float32)})
    assert "vector_fields" not in adapter.get_dataset_identity()


def test_anndata_discovers_only_exact_vector_keys_among_unrelated_obsm_data() -> None:
    adapter = _adapter(
        {
            "velocity_umap_2d": _vectors_2d(),
            "velocity_umap_history": np.ones((3, 7), dtype=np.float32),
        }
    )
    assert list(adapter.get_dataset_identity()["vector_fields"]["fields"]) == ["velocity_umap"]


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        (
            {"velocity_umap_3d": _embedding_3d()},
            "matching 3D embedding",
        ),
        (
            {
                "velocity_umap_2d": np.array(
                    [[0.0, np.inf], [0.0, 0.0], [0.0, 0.0]],
                )
            },
            "finite",
        ),
        (
            {
                "velocity_umap_2d": np.full(
                    (3, 2),
                    1e100,
                    dtype=np.float64,
                )
            },
            "float32 range",
        ),
    ],
)
def test_anndata_rejects_malformed_vector_fields(
    fields: dict[str, object],
    message: str,
) -> None:
    with pytest.raises((TypeError, ValueError), match=message):
        _adapter(fields)


def test_anndata_revalidates_vector_values_after_initialization() -> None:
    adata = _adata({"velocity_umap_2d": _vectors_2d()})
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Mutable vector field",
        dataset_id="mutable-vector-field",
    )
    adata.obsm["velocity_umap_2d"] = np.full(
        (3, 2),
        np.nan,
        dtype=np.float32,
    )

    with pytest.raises(ValueError, match="finite"):
        adapter.get_vector_field_binary("velocity_umap", 2)


def test_anndata_vector_binary_requires_native_field_id_and_dimension() -> None:
    adapter = _adapter({"velocity_umap_2d": _vectors_2d()})

    with pytest.raises(TypeError, match="field_id"):
        adapter.get_vector_field_binary(7, 2)
    with pytest.raises(TypeError, match="dim"):
        adapter.get_vector_field_binary("velocity_umap", 2.0)


def test_vector_field_default_is_explicit_in_every_public_vector_input_api() -> None:
    from cellucid.anndata_server import AnnDataServer, serve_anndata
    from cellucid.jupyter import AnnDataViewer, show_anndata

    public_apis = [
        prepare,
        AnnDataAdapter,
        AnnDataAdapter.from_file,
        AnnDataServer,
        serve_anndata,
        AnnDataViewer,
        show_anndata,
    ]
    for api in public_apis:
        assert "vector_field_default" in inspect.signature(api).parameters


def test_prepare_exposes_no_unsupported_4d_api() -> None:
    assert "X_umap_" + "4d" not in inspect.signature(prepare).parameters
