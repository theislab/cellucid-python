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
from cellucid.prepare_data import Embedding, prepare
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
    # Range 4, so its scale factor (0.5) differs from the 2D one (2.0) and a
    # vector scaled by the wrong embedding cannot pass.
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [0.0, 2.0, 1.0],
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


def _vectors_3d() -> np.ndarray:
    return np.array(
        [
            [0.1, 0.2, 0.3],
            [0.4, 0.5, 0.6],
            [0.7, 0.8, 0.9],
        ],
        dtype=np.float32,
    )


def _prepare_kwargs(
    out_dir: Path,
    vector_fields: object,
) -> dict[str, object]:
    embedding = _embedding_2d()
    return {
        "latent_space": embedding.copy(),
        "obs": pd.DataFrame(index=["cell-1", "cell-2", "cell-3"]),
        "embeddings": {"X_umap": embedding},
        "out_dir": out_dir,
        "dataset_name": "Exact vector fields",
        "dataset_id": "exact-vector-fields",
        "obs_categorical_dtype": "uint16",
        "centroid_min_points": 1,
        "vector_fields": vector_fields,
    }


def _identity(out_dir: Path) -> dict:
    return json.loads((out_dir / "dataset_identity.json").read_text(encoding="utf-8"))


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
        serve_vector_fields=True,
    )


def test_vector_fields_bind_arrays_to_embedding_ids(tmp_path: Path) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, {"velocity": _vectors_2d()})

    with pytest.raises(TypeError, match="must map embedding ids to arrays"):
        prepare(**kwargs)
    assert not out_dir.exists()

    # In AnnData the binding is declared by the obsm key, and a key without a
    # dimension suffix declares none.
    adapter = _adapter({"velocity_umap": _vectors_2d()})
    assert "vector_fields" not in adapter.get_dataset_identity()


def test_transition_drift_exposes_no_unsuffixed_output_switch() -> None:
    assert "explicit_dim_suffix" not in inspect.signature(add_transition_drift_to_obsm).parameters


_VECTOR_FIELD_KEY_ORDER = ["label", "files"]


def test_a_vector_field_entry_uses_the_specified_key_order(tmp_path: Path) -> None:
    # Key order is the order the output format specification prints and the
    # order cellucid-r writes, so one input produces one file whichever writer
    # produced it. A dict comparison passes whatever the order is, so the order
    # is asserted directly.
    prepare(**_prepare_kwargs(tmp_path, {"velocity": {"X_umap": _vectors_2d()}}))
    entry = _identity(tmp_path)["vector_fields"]["fields"]["velocity"]
    assert list(entry) == _VECTOR_FIELD_KEY_ORDER

    adapter = _adapter({"velocity_umap_2d": _vectors_2d()})
    served = adapter.get_dataset_identity()["vector_fields"]["fields"]["velocity_umap"]
    assert list(served) == _VECTOR_FIELD_KEY_ORDER


def test_prepare_single_vector_field_is_unambiguous_and_metadata_is_exact(
    tmp_path: Path,
) -> None:
    vectors = _vectors_2d()
    prepare(**_prepare_kwargs(tmp_path, {"velocity": {"X_umap": vectors}}))

    assert _identity(tmp_path)["vector_fields"] == {
        "default_field": "velocity",
        "fields": {
            "velocity": {
                "label": "velocity",
                "files": {"X_umap": "vectors/0_0.bin"},
            }
        },
    }
    exported = np.fromfile(tmp_path / "vectors" / "0_0.bin", dtype="<f4").reshape(3, 2)
    np.testing.assert_allclose(exported, vectors * 2.0)


def test_a_field_bound_to_two_embeddings_names_both_in_item_order(
    tmp_path: Path,
) -> None:
    kwargs = _prepare_kwargs(
        tmp_path / "export",
        # Supplied 3D first: the files object still follows item order.
        {"velocity": {"X_umap_3d": _vectors_3d(), "X_umap": _vectors_2d()}},
    )
    kwargs["embeddings"] = {"X_umap": _embedding_2d(), "X_umap_3d": _embedding_3d()}
    prepare(**kwargs)

    field = _identity(tmp_path / "export")["vector_fields"]["fields"]["velocity"]
    assert list(field["files"].items()) == [
        ("X_umap", "vectors/0_0.bin"),
        ("X_umap_3d", "vectors/0_1.bin"),
    ]
    vectors_3d = np.fromfile(tmp_path / "export" / "vectors" / "0_1.bin", dtype="<f4")
    np.testing.assert_array_equal(
        vectors_3d,
        (_vectors_3d().astype(np.float64) * 0.5).astype(np.float32).ravel(),
    )

    adata = _adata({"velocity_umap_2d": _vectors_2d(), "velocity_umap_3d": _vectors_3d()})
    adata.obsm["X_umap_3d"] = _embedding_3d()
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Exact vector bindings",
        dataset_id="exact-vector-bindings",
        serve_vector_fields=True,
    )
    served = adapter.get_dataset_identity()["vector_fields"]["fields"]["velocity_umap"]
    assert list(served["files"].items()) == [
        ("X_umap_2d", "vectors/0_0.bin"),
        ("X_umap_3d", "vectors/0_1.bin"),
    ]


def test_a_vector_is_scaled_by_the_embedding_it_is_bound_to(tmp_path: Path) -> None:
    # X_tsne spans 8 units (scale 0.25), X_umap spans 1 (scale 2.0).
    tsne = _embedding_2d() * 8.0
    kwargs = _prepare_kwargs(tmp_path, {"velocity": {"X_tsne": _vectors_2d()}})
    kwargs["embeddings"] = {"X_umap": _embedding_2d(), "X_tsne": tsne}
    kwargs["embedding_defaults"] = {"2d": "X_umap"}
    prepare(**kwargs)

    identity = _identity(tmp_path)
    assert list(identity["embeddings"]["items"]) == ["X_tsne", "X_umap"]
    assert identity["vector_fields"]["fields"]["velocity"]["files"] == {
        "X_tsne": "vectors/0_0.bin"
    }
    exported = np.fromfile(tmp_path / "vectors" / "0_0.bin", dtype="<f4")
    np.testing.assert_array_equal(
        exported,
        (_vectors_2d().astype(np.float64) * 0.25).astype(np.float32).ravel(),
    )


def test_dataset_default_dimension_is_highest_declared_dimension(
    tmp_path: Path,
) -> None:
    embedding_1d = np.array([[0.0], [1.0], [2.0]], dtype=np.float32)
    kwargs = _prepare_kwargs(tmp_path / "export", {})
    kwargs["embeddings"] = {
        "X_umap_3d": _embedding_3d(),
        "X_umap_1d": embedding_1d,
        "X_umap": _embedding_2d(),
    }
    prepare(**kwargs)
    prepared = _identity(tmp_path / "export")["embeddings"]
    assert list(prepared["items"]) == ["X_umap_1d", "X_umap", "X_umap_3d"]
    assert prepared["default_dimension"] == 3

    adata = _adata({})
    adata.obsm["X_umap_1d"] = embedding_1d
    adata.obsm["X_umap_3d"] = _embedding_3d()
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Exact dataset dimensions",
        dataset_id="exact-dataset-dimensions",
        serve_vector_fields=True,
    )
    direct = adapter.get_dataset_identity()["embeddings"]
    assert list(direct["items"]) == ["X_umap_1d", "X_umap_2d", "X_umap_3d"]
    assert direct["default_dimension"] == 3


def test_sparse_vector_values_follow_the_same_exact_contract(
    tmp_path: Path,
) -> None:
    vectors = _vectors_2d()
    sparse_vectors = sparse.csr_matrix(vectors)
    prepare(**_prepare_kwargs(tmp_path / "export", {"velocity": {"X_umap": sparse_vectors}}))
    prepared = np.fromfile(
        tmp_path / "export" / "vectors" / "0_0.bin",
        dtype="<f4",
    ).reshape(3, 2)
    np.testing.assert_allclose(prepared, vectors * 2.0)

    adapter = _adapter({"velocity_umap_2d": sparse_vectors})
    direct = np.frombuffer(
        adapter.get_vector_field_binary("velocity_umap", "X_umap_2d"),
        dtype="<f4",
    ).reshape(3, 2)
    np.testing.assert_allclose(direct, vectors * 2.0)


def test_prepare_multiple_vector_fields_require_an_exact_explicit_default(
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "export"
    fields = {
        "velocity": {"X_umap": _vectors_2d()},
        "drift": {"X_umap": _vectors_2d(0.1)},
    }
    kwargs = _prepare_kwargs(out_dir, fields)

    with pytest.raises(ValueError, match="vector_field_default"):
        prepare(**kwargs)
    assert not out_dir.exists()

    kwargs["vector_field_default"] = "drift"
    prepare(**kwargs)
    identity = _identity(out_dir)
    assert identity["vector_fields"]["default_field"] == "drift"
    # Fields are indexed in code-point order, whatever order they were given in.
    assert list(identity["vector_fields"]["fields"]) == ["drift", "velocity"]


def test_case_distinct_vector_field_ids_are_two_indexed_fields(
    tmp_path: Path,
) -> None:
    """Vector payloads are named by index, so case no longer collides at a path."""
    out_dir = tmp_path / "export"
    fields = {
        "Velocity": {"X_umap": _vectors_2d()},
        "velocity": {"X_umap": _vectors_2d(0.1)},
    }
    kwargs = _prepare_kwargs(out_dir, fields)
    kwargs["vector_field_default"] = "Velocity"

    prepare(**kwargs)

    identity = _identity(out_dir)
    assert {
        field_id: field["files"] for field_id, field in identity["vector_fields"]["fields"].items()
    } == {
        "Velocity": {"X_umap": "vectors/0_0.bin"},
        "velocity": {"X_umap": "vectors/1_0.bin"},
    }
    assert sorted(path.name for path in (out_dir / "vectors").iterdir()) == [
        "0_0.bin",
        "1_0.bin",
    ]

    adapter = _adapter(
        {"Velocity_umap_2d": _vectors_2d(), "velocity_umap_2d": _vectors_2d(0.1)},
        vector_field_default="Velocity_umap",
    )
    assert adapter.get_vector_field_route("1_0.bin") == ("velocity_umap", "X_umap_2d")
    assert adapter.get_vector_field_route("2_0.bin") is None
    assert adapter.get_vector_field_route("0_1.bin") is None


def test_vector_default_without_vector_fields_is_rejected_before_mutation(
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, {})
    kwargs["vector_field_default"] = "velocity"

    with pytest.raises(ValueError, match="no vector fields"):
        prepare(**kwargs)
    assert not out_dir.exists()

    with pytest.raises(ValueError, match="no vector fields"):
        _adapter({}, vector_field_default="velocity_umap")


@pytest.mark.parametrize(
    "invalid_default",
    ["", "missing", 7],
)
def test_prepare_rejects_invalid_vector_default_before_output_mutation(
    tmp_path: Path,
    invalid_default: object,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, {"velocity": {"X_umap": _vectors_2d()}})
    kwargs["vector_field_default"] = invalid_default

    with pytest.raises((TypeError, ValueError), match="vector_field_default"):
        prepare(**kwargs)
    assert not out_dir.exists()


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        ({"": {"X_umap": _vectors_2d()}}, "vector field id"),
        ({7: {"X_umap": _vectors_2d()}}, "vector field id"),
        ({"velocity\n": {"X_umap": _vectors_2d()}}, "displayed verbatim"),
        ({"velocity": _vectors_2d()}, "must map embedding ids"),
        ({"velocity": {}}, "binds no embedding"),
        ({"velocity": {"X_tsne": _vectors_2d()}}, "not one of the exported embeddings"),
        ({"velocity": {7: _vectors_2d()}}, "embedding id"),
        ({"velocity": {"": _vectors_2d()}}, "embedding id"),
    ],
)
def test_prepare_rejects_invalid_field_ids_and_bindings_before_mutation(
    tmp_path: Path,
    fields: object,
    message: str,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, fields)

    with pytest.raises((TypeError, ValueError), match=message):
        prepare(**kwargs)
    assert not out_dir.exists()


@pytest.mark.parametrize("fields", ["velocity", [("velocity", {})]])
def test_prepare_rejects_vector_fields_that_are_not_a_mapping(
    tmp_path: Path,
    fields: object,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(TypeError, match="vector_fields must map each field id"):
        prepare(**_prepare_kwargs(out_dir, fields))
    assert not out_dir.exists()


@pytest.mark.parametrize(
    ("value", "message"),
    [
        (None, "must not be None"),
        (np.arange(3, dtype=np.float32), "2D array"),
        (np.ones((2, 2), dtype=np.float32), "exactly 3 rows"),
        (np.ones((3, 3), dtype=np.float32), "exactly 2 columns"),
        (np.ones((3, 4), dtype=np.float32), "exactly 2 columns"),
        (np.ones((3, 2), dtype=np.complex64), "real numeric"),
        (np.ones((3, 2), dtype=np.bool_), "real numeric"),
        (np.array([[0.0, np.nan], [0.0, 0.0], [0.0, 0.0]]), "finite"),
        (np.full((3, 2), 1e100, dtype=np.float64), "float32 range"),
    ],
)
def test_prepare_rejects_invalid_vector_values_before_output_mutation(
    tmp_path: Path,
    value: object,
    message: str,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, {"velocity": {"X_umap": value}})

    with pytest.raises(
        (TypeError, ValueError),
        match=rf"Vector field 'velocity' on embedding 'X_umap'.*{message}",
    ):
        prepare(**kwargs)
    assert not out_dir.exists()


def test_a_vector_has_the_width_of_the_embedding_it_is_bound_to(tmp_path: Path) -> None:
    # A wide source read at two columns is a 2D embedding, so its vectors are 2D.
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, {"velocity": {"X_pca": _vectors_3d()}})
    kwargs["embeddings"] = {
        "X_pca": Embedding(np.tile(_embedding_2d(), (1, 3)), columns=[0, 1]),
    }

    with pytest.raises(ValueError, match="exactly 2 columns"):
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
                "files": {"X_umap_2d": "vectors/0_0.bin"},
            }
        },
    }
    exported = np.frombuffer(
        adapter.get_vector_field_binary("velocity_umap", "X_umap_2d"),
        dtype="<f4",
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


def test_anndata_umap_vector_binds_only_to_a_slot_showing_that_umap() -> None:
    # No 3D slot: the 3D vector has nothing to bind to.
    assert "vector_fields" not in _adapter({"velocity_umap_3d": _vectors_3d()}).get_dataset_identity()

    # The 2D slot shows t-SNE, and a UMAP vector drawn over it would point the
    # wrong way, so it is not bound.
    adata = _adata({"velocity_umap_2d": _vectors_2d()})
    adata.obsm["X_tsne"] = _embedding_2d() * 3.0
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Exact vector slot",
        dataset_id="exact-vector-slot",
        embedding_keys={"2d": "X_tsne"},
        serve_vector_fields=True,
    )
    assert "vector_fields" not in adapter.get_dataset_identity()


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        (
            {"velocity_umap_2d": np.array([[0.0, np.inf], [0.0, 0.0], [0.0, 0.0]])},
            "finite",
        ),
        (
            {"velocity_umap_2d": np.full((3, 2), 1e100, dtype=np.float64)},
            "float32 range",
        ),
        (
            {"velocity_umap_2d": np.ones((3, 3), dtype=np.float32)},
            "exactly 2 columns",
        ),
    ],
)
def test_anndata_rejects_malformed_vector_fields(
    fields: dict[str, object],
    message: str,
) -> None:
    with pytest.raises((TypeError, ValueError), match=rf"'velocity_umap_2d'.*{message}"):
        _adapter(fields)


def test_anndata_revalidates_vector_values_after_initialization() -> None:
    adata = _adata({"velocity_umap_2d": _vectors_2d()})
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Mutable vector field",
        dataset_id="mutable-vector-field",
        serve_vector_fields=True,
    )
    adata.obsm["velocity_umap_2d"] = np.full(
        (3, 2),
        np.nan,
        dtype=np.float32,
    )

    with pytest.raises(ValueError, match="finite"):
        adapter.get_vector_field_binary("velocity_umap", "X_umap_2d")


def test_anndata_vector_binary_requires_native_field_and_embedding_ids() -> None:
    adapter = _adapter({"velocity_umap_2d": _vectors_2d()})

    with pytest.raises(TypeError, match="field_id"):
        adapter.get_vector_field_binary(7, "X_umap_2d")
    with pytest.raises(TypeError, match="embedding_id"):
        adapter.get_vector_field_binary("velocity_umap", 2)
    with pytest.raises(ValueError, match="not bound to embedding"):
        adapter.get_vector_field_binary("velocity_umap", "X_umap_3d")


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


def test_prepare_exposes_embeddings_and_no_per_dimension_parameters() -> None:
    parameters = inspect.signature(prepare).parameters
    assert "embeddings" in parameters
    assert "embedding_defaults" in parameters
    for dimension in (1, 2, 3, 4):
        assert f"X_umap_{dimension}d" not in parameters
