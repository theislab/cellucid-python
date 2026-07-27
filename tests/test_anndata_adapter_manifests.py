from __future__ import annotations

from unittest import mock

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid.anndata_adapter import (
    AnnDataAdapter,
    _categorical_storage,
)
from cellucid.prepare_data import (
    _assert_unique_filename_components,
    _require_portable_filename_component,
    _select_category_dtype,
)
from cellucid.prepare_data import (
    _json_category_values as _prepare_json_category_values,
)


def _adata_with_current_fields() -> ad.AnnData:
    obs = pd.DataFrame(
        {
            "score": np.array([0.25, 0.5, 0.75], dtype=np.float32),
            "cluster": pd.Categorical.from_codes(
                [0, 1, 0],
                categories=pd.Index(["A", "B"], dtype=object),
            ),
        },
        index=pd.Index(["cell-1", "cell-2", "cell-3"], dtype=object),
    )
    var = pd.DataFrame(index=pd.Index(["Gene_A", "Gene-B"], dtype=object))
    adata = ad.AnnData(
        X=np.array(
            [
                [0.0, 1.0],
                [2.0, 3.0],
                [4.0, 5.0],
            ],
            dtype=np.float32,
        ),
        obs=obs,
        var=var,
    )
    adata.obsm["X_umap_2d"] = np.array(
        [[0.0, 0.0], [1.0, 0.5], [0.5, 1.0]],
        dtype=np.float32,
    )
    adata.obsm["X_latent"] = np.array(
        [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]],
        dtype=np.float32,
    )
    return adata


def test_anndata_adapter_emits_exact_current_compact_manifests() -> None:
    adapter = AnnDataAdapter(
        _adata_with_current_fields(),
        latent_key="X_latent",
        dataset_name="Exact current",
        dataset_id="exact-current",
    )

    obs = adapter.get_obs_manifest()
    assert list(obs) == [
        "_format",
        "n_points",
        "centroid_outlier_quantile",
        "latent_key",
        "compression",
        "_obsSchemas",
        "_continuousFields",
        "_categoricalFields",
    ]
    assert obs["_format"] == "compact_v1"
    assert obs["_continuousFields"] == [["score"]]
    assert len(obs["_categoricalFields"]) == 1
    categorical = obs["_categoricalFields"][0]
    assert categorical[:4] == ["cluster", ["A", "B"], "uint8", 255]
    assert isinstance(categorical[4], dict)
    assert len(categorical) == 5
    assert obs["_obsSchemas"] == {
        "continuous": {
            "pathPattern": "obs/{key}.values.f32",
            "ext": "f32",
            "dtype": "float32",
            "quantized": False,
        },
        "categorical": {
            "codesPathPattern": "obs/{key}.codes.{ext}",
            "outlierPathPattern": "obs/{key}.outliers.f32",
            "outlierExt": "f32",
            "outlierDtype": "float32",
            "outlierQuantized": False,
        },
    }

    var = adapter.get_var_manifest()
    assert list(var) == [
        "_format",
        "n_points",
        "var_gene_id_column",
        "compression",
        "quantization",
        "_varSchema",
        "fields",
    ]
    assert var["_format"] == "compact_v1"
    assert var["var_gene_id_column"] is None
    assert var["quantization"] is None
    assert var["_varSchema"] == {
        "kind": "continuous",
        "pathPattern": "var/{key}.values.f32",
        "ext": "f32",
        "dtype": "float32",
        "quantized": False,
    }
    assert var["fields"] == [["Gene_A"], ["Gene-B"]]

    identity = adapter.get_dataset_identity()
    assert "_anndata_adapter" not in identity
    assert list(identity["source"]) == ["name"]
    assert identity["source"] == {"name": "In-memory AnnData"}
    assert identity["export_settings"] == {
        "compression": None,
        "var_quantization": None,
        "obs_continuous_quantization": None,
        "obs_categorical_dtype": "auto",
    }


def test_identity_obs_fields_follow_the_exact_compact_manifest_order() -> None:
    adata = _adata_with_current_fields()
    adata.obs = adata.obs[["cluster", "score"]]
    adapter = AnnDataAdapter(
        adata,
        latent_key="X_latent",
        dataset_name="Canonical observation order",
        dataset_id="canonical-observation-order",
    )

    manifest = adapter.get_obs_manifest()
    assert [field[0] for field in manifest["_continuousFields"]] == ["score"]
    assert [field[0] for field in manifest["_categoricalFields"]] == ["cluster"]
    identity = adapter.get_dataset_identity()
    assert identity["obs_fields"] == [
        {"key": "score", "kind": "continuous"},
        {"key": "cluster", "kind": "category", "n_categories": 2},
    ]
    assert identity["stats"]["n_obs_fields"] == 2
    assert identity["stats"]["n_continuous_fields"] == 1
    assert identity["stats"]["n_categorical_fields"] == 1


def test_anndata_adapter_emits_only_schemas_owned_by_present_fields() -> None:
    adata = _adata_with_current_fields()
    del adata.obs["cluster"]
    continuous_adapter = AnnDataAdapter(
        adata,
        latent_key="X_latent",
        dataset_name="Continuous only",
        dataset_id="continuous-only",
    )
    assert list(continuous_adapter.get_obs_manifest()["_obsSchemas"]) == ["continuous"]

    del adata.obs["score"]
    empty_adapter = AnnDataAdapter(
        adata,
        latent_key="X_latent",
        dataset_name="No fields",
        dataset_id="no-fields",
    )
    assert empty_adapter.get_obs_manifest()["_obsSchemas"] == {}


def test_anndata_adapter_marks_outliers_absent_without_an_explicit_latent_key() -> None:
    adata = _adata_with_current_fields()
    adapter = AnnDataAdapter(
        adata,
        dataset_name="No latent",
        dataset_id="no-latent",
    )
    categorical_schema = adapter.get_obs_manifest()["_obsSchemas"]["categorical"]
    assert categorical_schema == {
        "codesPathPattern": "obs/{key}.codes.{ext}",
        "outlierPathPattern": None,
        "outlierExt": None,
        "outlierDtype": None,
        "outlierQuantized": False,
    }
    with pytest.raises(ValueError, match="no declared latent space"):
        adapter.get_obs_outlier_quantiles("cluster")

    with pytest.raises(ValueError, match="was not found in adata.obsm"):
        AnnDataAdapter(
            adata,
            latent_key="missing",
            dataset_name="Missing latent",
            dataset_id="missing-latent",
        )


def test_anndata_adapter_does_not_advertise_or_invent_missing_expression() -> None:
    adata = _adata_with_current_fields()
    adata.X = None
    adapter = AnnDataAdapter(
        adata,
        latent_key="X_latent",
        dataset_name="No expression",
        dataset_id="no-expression",
    )

    assert adapter.get_var_manifest()["fields"] == []
    assert adapter.get_dataset_identity()["stats"]["n_genes"] == 0
    with pytest.raises(ValueError, match="no X expression matrix"):
        adapter.get_gene_expression("Gene_A")


def test_anndata_adapter_rejects_zero_observation_datasets() -> None:
    adata = ad.AnnData(
        X=np.empty((0, 1), dtype=np.float32),
        obs=pd.DataFrame(index=pd.Index([], dtype=object)),
        var=pd.DataFrame(index=pd.Index(["Gene_A"], dtype=object)),
    )
    with pytest.raises(ValueError, match="at least one observation"):
        AnnDataAdapter(
            adata,
            dataset_name="Empty",
            dataset_id="empty",
        )


def test_anndata_adapter_preserves_boolean_category_identity() -> None:
    adata = _adata_with_current_fields()
    adata.obs["flag"] = pd.Categorical([False, True, False])
    adapter = AnnDataAdapter(
        adata,
        latent_key="X_latent",
        dataset_name="Boolean categories",
        dataset_id="boolean-categories",
    )

    categorical_fields = adapter.get_obs_manifest()["_categoricalFields"]
    flag = next(field for field in categorical_fields if field[0] == "flag")
    assert flag[1] == [False, True]
    _data, categories, missing_value = adapter.get_obs_categorical_codes("flag")
    assert categories == [False, True]
    assert missing_value == 255


def test_python_producers_share_exact_categorical_storage_boundaries() -> None:
    assert _categorical_storage(255, field_key="group") == (
        np.uint8,
        "uint8",
        255,
    )
    assert _categorical_storage(256, field_key="group") == (
        np.uint16,
        "uint16",
        65_535,
    )
    assert _categorical_storage(65_535, field_key="group") == (
        np.uint16,
        "uint16",
        65_535,
    )
    with pytest.raises(ValueError, match="at most 65,535"):
        _categorical_storage(65_536, field_key="group")

    assert _select_category_dtype(255) == (np.uint8, 255)
    assert _select_category_dtype(256) == (np.uint16, 65_535)
    assert _select_category_dtype(65_535) == (np.uint16, 65_535)
    with pytest.raises(ValueError, match="at most 65,535"):
        _select_category_dtype(65_536)


def test_prepare_producer_preserves_exact_json_category_identity() -> None:
    assert _prepare_json_category_values(
        [np.bool_(False), np.bool_(True), np.int64(2), np.float64(2.5), "2"],
        field_key="group",
    ) == [False, True, 2, 2.5, "2"]
    with pytest.raises(ValueError, match="exact JSON representation"):
        _prepare_json_category_values([1, 1.0], field_key="group")
    with pytest.raises(ValueError, match="exact integer range"):
        _prepare_json_category_values([2**53], field_key="group")


@pytest.mark.parametrize(
    "component",
    [
        "",
        "_leading",
        ".leading",
        "-leading",
        "trailing.",
        "A B",
        "A/B",
        "ümlaut",
        "a" * 181,
        "CON",
        "con.txt",
        "PRN.csv",
        "AUX",
        "NUL.json",
        "COM1",
        "com9.bin",
        "LPT1",
        "lpt9.txt",
    ],
)
def test_portable_filename_component_rejects_every_cross_platform_hazard(
    component: str,
) -> None:
    with pytest.raises(ValueError, match="ASCII|reserved Windows|end with"):
        _require_portable_filename_component(component)


def test_prepared_components_are_portable_while_direct_routes_preserve_case() -> None:
    valid = [
        "A",
        "A.b_c-",
        "a" * 180,
    ]
    assert [_require_portable_filename_component(component) for component in valid] == valid

    with pytest.raises(ValueError, match="case-insensitively"):
        _assert_unique_filename_components(
            ["Field", "field"],
            label="Observation field",
        )

    adata = _adata_with_current_fields()
    adata.obs = pd.DataFrame(
        {
            "Field": [1.0, 2.0, 3.0],
            "field": [4.0, 5.0, 6.0],
        },
        index=adata.obs_names,
    )
    adapter = AnnDataAdapter(
        adata,
        latent_key="X_latent",
        dataset_name="Case-distinct obs",
        dataset_id="case-distinct-obs",
    )
    assert [field[0] for field in adapter.get_obs_manifest()["_continuousFields"]] == [
        "Field",
        "field",
    ]

    adata = _adata_with_current_fields()
    adata.var_names = ["Gene", "gene"]
    adapter = AnnDataAdapter(
        adata,
        latent_key="X_latent",
        dataset_name="Case-distinct var",
        dataset_id="case-distinct-var",
    )
    assert adapter.get_var_manifest()["fields"] == [["Gene"], ["gene"]]

    adata.var_names = ["A/B", "A B"]
    with pytest.raises(ValueError, match="collide.*A_B"):
        AnnDataAdapter(
            adata,
            latent_key="X_latent",
            dataset_name="Wire collision",
            dataset_id="wire-collision",
        )


def test_portable_filename_component_requires_a_native_string() -> None:
    with pytest.raises(TypeError, match="native string"):
        _require_portable_filename_component(7)


@pytest.mark.parametrize(
    "invalid_value",
    [np.nan, np.inf, -np.inf, np.finfo(np.float64).max],
)
def test_continuous_obs_requires_finite_float32_values_at_initialization(
    invalid_value: float,
) -> None:
    adata = _adata_with_current_fields()
    adata.obs["score"] = np.array([0.25, invalid_value, 0.75], dtype=np.float64)

    with pytest.raises(ValueError, match=r"score.*finite|score.*float32"):
        AnnDataAdapter(
            adata,
            latent_key="X_latent",
            dataset_name="Invalid continuous obs",
            dataset_id="invalid-continuous-obs",
        )


def test_continuous_obs_is_revalidated_after_adapter_initialization() -> None:
    adata = _adata_with_current_fields()
    adapter = AnnDataAdapter(
        adata,
        latent_key="X_latent",
        dataset_name="Mutable continuous obs",
        dataset_id="mutable-continuous-obs",
    )
    adata.obs["score"] = np.array([0.25, np.nan, 0.75], dtype=np.float64)

    with pytest.raises(ValueError, match=r"score.*finite"):
        adapter.get_obs_continuous_values("score")


@pytest.mark.parametrize(
    ("keyword", "value"),
    [
        ("centroid_outlier_quantile", 0.5),
        ("centroid_min_points", 0),
    ],
)
def test_anndata_adapter_rejects_invalid_centroid_configuration(
    keyword: str,
    value: float,
) -> None:
    with pytest.raises(ValueError):
        AnnDataAdapter(
            _adata_with_current_fields(),
            latent_key="X_latent",
            dataset_name="Invalid centroids",
            dataset_id="invalid-centroids",
            **{keyword: value},
        )


@pytest.mark.parametrize(
    ("metadata", "missing_key"),
    [
        ({"dataset_id": "declared-id"}, "dataset_name"),
        ({"dataset_name": "Declared name"}, "dataset_id"),
    ],
)
def test_anndata_adapter_requires_explicit_dataset_identity(
    metadata: dict[str, str],
    missing_key: str,
) -> None:
    with pytest.raises(TypeError, match=missing_key):
        AnnDataAdapter(_adata_with_current_fields(), **metadata)


@pytest.mark.parametrize(
    ("dataset_name", "dataset_id", "invalid_key"),
    [
        ("", "declared-id", "dataset_name"),
        ("Declared name", "", "dataset_id"),
        (7, "declared-id", "dataset_name"),
        ("Declared name", 7, "dataset_id"),
    ],
)
def test_anndata_adapter_dataset_identity_is_exact_nonempty_text(
    dataset_name: object,
    dataset_id: object,
    invalid_key: str,
) -> None:
    with pytest.raises((TypeError, ValueError), match=invalid_key):
        AnnDataAdapter(
            _adata_with_current_fields(),
            dataset_name=dataset_name,
            dataset_id=dataset_id,
        )


def test_malformed_explicit_umap_cannot_be_replaced_by_generic_umap() -> None:
    adata = _adata_with_current_fields()
    adata.obsm["X_umap_2d"] = np.ones((adata.n_obs, 3), dtype=np.float32)

    with pytest.raises(ValueError, match=r"X_umap_2d.*exactly 2 columns"):
        AnnDataAdapter(
            adata,
            dataset_name="Malformed explicit UMAP",
            dataset_id="malformed-explicit-umap",
        )


def test_generic_umap_is_ignored_when_an_explicit_embedding_exists() -> None:
    adata = _adata_with_current_fields()
    adata.obsm["X_umap"] = np.asarray(adata.obsm["X_umap_2d"]).copy()

    adapter = AnnDataAdapter(
        adata,
        dataset_name="Explicit UMAP",
        dataset_id="explicit-umap",
    )
    assert adapter.get_dataset_identity()["embeddings"]["available_dimensions"] == [2]


@pytest.mark.parametrize("invalid_value", [np.nan, np.inf, -np.inf])
def test_nonfinite_umap_rejects_during_adapter_initialization(
    invalid_value: float,
) -> None:
    adata = _adata_with_current_fields()
    embedding = np.asarray(adata.obsm["X_umap_2d"]).copy()
    embedding[1, 0] = invalid_value
    adata.obsm["X_umap_2d"] = embedding

    with pytest.raises(ValueError, match=r"X_umap.*finite"):
        AnnDataAdapter(
            adata,
            dataset_name="Nonfinite UMAP",
            dataset_id="nonfinite-umap",
        )


def test_degenerate_umap_rejects_instead_of_inventing_a_normalization_range() -> None:
    adata = _adata_with_current_fields()
    adata.obsm["X_umap_2d"] = np.ones((adata.n_obs, 2), dtype=np.float32)

    with pytest.raises(ValueError, match=r"X_umap.*no coordinate variation"):
        AnnDataAdapter(
            adata,
            dataset_name="Degenerate UMAP",
            dataset_id="degenerate-umap",
        )


def test_extreme_finite_umap_normalizes_without_float32_range_overflow() -> None:
    adata = _adata_with_current_fields()
    limit = np.finfo(np.float32).max
    adata.obsm["X_umap_2d"] = np.array(
        [[-limit, 0.0], [limit, 0.0], [0.0, limit]],
        dtype=np.float32,
    )
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Extreme finite UMAP",
        dataset_id="extreme-finite-umap",
    )

    normalized = adapter.get_embedding(2)
    assert np.isfinite(normalized).all()
    assert normalized.min() >= -1.0
    assert normalized.max() <= 1.0


def test_embedding_shape_mutation_after_initialization_is_rejected() -> None:
    adata = _adata_with_current_fields()
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Mutable UMAP",
        dataset_id="mutable-umap",
    )
    adata.obsm["X_umap_2d"] = np.ones((adata.n_obs, 3), dtype=np.float32)

    with pytest.raises(ValueError, match=r"changed after adapter initialization"):
        adapter.get_embedding(2)


def test_numeric_gene_identifier_column_is_not_string_coerced() -> None:
    adata = _adata_with_current_fields()
    adata.var["gene_id"] = [7, 8]

    with pytest.raises(
        TypeError,
        match=r"var column 'gene_id' identifier at position 0 must be a string",
    ):
        AnnDataAdapter(
            adata,
            gene_id_column="gene_id",
            dataset_name="Numeric gene identifiers",
            dataset_id="numeric-gene-identifiers",
        )


def test_zarr_source_is_not_marked_as_lazy_when_anndata_materializes_it(
    tmp_path,
) -> None:
    store = tmp_path / "fixture.zarr"
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")
    adata = _adata_with_current_fields()

    with mock.patch("anndata.read_zarr", return_value=adata):
        adapter = AnnDataAdapter.from_file(
            store,
            dataset_name="Materialized Zarr",
            dataset_id="materialized-zarr",
        )

    assert adapter.is_backed is False
    assert not hasattr(adapter, "_is_lazy")
    assert adapter.get_dataset_identity()["source"] == {
        "name": "Zarr store",
    }


def test_from_file_does_not_derive_identity_from_private_path(tmp_path) -> None:
    store = tmp_path / "must-not-become-a-dataset-name.zarr"
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")

    with (
        mock.patch("anndata.read_zarr", return_value=_adata_with_current_fields()),
        pytest.raises(TypeError, match="dataset_name"),
    ):
        AnnDataAdapter.from_file(
            store,
            dataset_id="declared-id",
        )


def test_real_h5ad_backed_loading_reports_actual_state_and_serves_values(
    tmp_path,
) -> None:
    source = tmp_path / "fixture.h5ad"
    _adata_with_current_fields().write_h5ad(source)

    with AnnDataAdapter.from_file(
        source,
        dataset_name="Backed H5AD",
        dataset_id="backed-h5ad",
    ) as adapter:
        assert adapter.is_backed is True
        assert not hasattr(adapter, "_is_lazy")
        assert adapter.get_dataset_identity()["source"] == {
            "name": "H5AD file",
        }
        values = np.frombuffer(
            adapter.get_gene_expression("Gene_A"),
            dtype=np.float32,
        )
        np.testing.assert_array_equal(values, np.array([0.0, 2.0, 4.0]))


def test_real_zarr_loading_reports_materialized_state(tmp_path) -> None:
    pytest.importorskip("zarr")
    store = tmp_path / "fixture.zarr"
    with pytest.warns(UserWarning, match="Writing zarr v2 data"):
        _adata_with_current_fields().write_zarr(store)

    with AnnDataAdapter.from_file(
        store,
        dataset_name="Materialized Zarr",
        dataset_id="materialized-zarr-real",
    ) as adapter:
        assert adapter.is_backed is False
        assert not hasattr(adapter, "_is_lazy")
        assert adapter.get_dataset_identity()["source"] == {
            "name": "Zarr store",
        }
        np.testing.assert_array_equal(
            adapter.get_gene_ids(),
            ["Gene_A", "Gene-B"],
        )
