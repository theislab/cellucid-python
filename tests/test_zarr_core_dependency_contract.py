from __future__ import annotations

from importlib.metadata import version
from pathlib import Path
from unittest import mock

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from packaging.requirements import Requirement
from packaging.version import Version

from cellucid.anndata_adapter import AnnDataAdapter, _classify_anndata_path
from cellucid.cli import _detect_data_format

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_core_metadata_declares_one_current_bounded_direct_data_runtime() -> None:
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    anndata_declaration = Requirement("anndata>=0.12.19,<0.13")
    zarr_declaration = Requirement("zarr>=3.1.4,<4")
    numcodecs_declaration = Requirement("numcodecs>=0.16.3,<0.17")
    assert pyproject.count('"anndata>=0.12.19,<0.13",') == 1
    assert pyproject.count('"zarr>=3.1.4,<4",') == 1
    assert pyproject.count('"numcodecs>=0.16.3,<0.17",') == 1
    assert Version(version("anndata")) in anndata_declaration.specifier
    assert Version(version("zarr")) in zarr_declaration.specifier
    assert Version(version("numcodecs")) in numcodecs_declaration.specifier


@pytest.mark.parametrize("suffix", [".zarr", ""])
def test_python_zarr_reader_rejects_every_non_v2_root(
    tmp_path: Path,
    suffix: str,
) -> None:
    store = tmp_path / f"format-three{suffix}"
    store.mkdir()
    (store / "zarr.json").write_text(
        '{"zarr_format":3,"node_type":"group"}',
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="exact Zarr v2 root contract"):
        _classify_anndata_path(store)
    assert _detect_data_format(store) == "unknown"


def test_zarr_read_failure_propagates_without_optional_dependency_translation(
    tmp_path: Path,
) -> None:
    store = tmp_path / "exact-v2.zarr"
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")
    failure = ModuleNotFoundError("synthetic exact read failure", name="zarr")

    with (
        mock.patch("anndata.read_zarr", side_effect=failure),
        pytest.raises(ModuleNotFoundError) as raised,
    ):
        AnnDataAdapter.from_file(
            store,
            dataset_name="Exact failure",
            dataset_id="exact-failure",
        )

    assert raised.value is failure


def test_h5ad_read_failure_propagates_without_dependency_translation(
    tmp_path: Path,
) -> None:
    source = tmp_path / "exact.h5ad"
    source.touch()
    failure = ModuleNotFoundError("synthetic exact read failure", name="h5py")

    with (
        mock.patch("anndata.read_h5ad", side_effect=failure),
        pytest.raises(ModuleNotFoundError) as raised,
    ):
        AnnDataAdapter.from_file(
            source,
            dataset_name="Exact failure",
            dataset_id="exact-failure",
        )

    assert raised.value is failure


def test_real_zarr_v2_store_runs_the_read_and_serve_preparation_path(
    tmp_path: Path,
) -> None:
    source = ad.AnnData(
        X=np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]], dtype=np.float32),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical.from_codes(
                    [0, 1, 0],
                    categories=pd.Index(["alpha", "beta"], dtype=object),
                ),
                "score": np.array([0.25, 0.5, 0.75], dtype=np.float32),
            },
            index=pd.Index(
                ["cell-1", "cell-2", "cell-3"],
                dtype=object,
            ),
        ),
        var=pd.DataFrame(index=pd.Index(["Gene_A", "Gene_B"], dtype=object)),
        obsm={
            "X_umap_2d": np.array(
                [[-2.0, 1.0], [0.0, 0.0], [2.0, -1.0]],
                dtype=np.float32,
            )
        },
    )
    store = tmp_path / "serve-ready.zarr"
    with pytest.warns(UserWarning, match="Writing zarr v2 data"):
        source.write_zarr(store)
    assert (store / ".zgroup").is_file()
    assert (store / ".zattrs").is_file()
    assert not (store / "zarr.json").exists()

    with AnnDataAdapter.from_file(
        store,
        dataset_name="Zarr v2 contract",
        dataset_id="zarr-v2-contract",
    ) as adapter:
        assert adapter.is_backed is False
        identity = adapter.get_dataset_identity()
        assert identity["source"] == {"name": "Zarr store"}
        assert identity["embeddings"] == {
            "available_dimensions": [2],
            "default_dimension": 2,
            "files": {"2d": "points_2d.bin"},
        }
        obs_manifest = adapter.get_obs_manifest()
        assert obs_manifest["_continuousFields"] == [["score"]]
        assert obs_manifest["_categoricalFields"] == [
            [
                "cell_type",
                ["alpha", "beta"],
                "uint8",
                255,
                {"2": []},
            ]
        ]
        assert adapter.get_var_manifest()["fields"] == [["Gene_A"], ["Gene_B"]]
        np.testing.assert_allclose(
            np.frombuffer(adapter.get_points_binary(2), dtype=np.float32).reshape(3, 2),
            np.array([[-1.0, 0.5], [0.0, 0.0], [1.0, -0.5]], dtype=np.float32),
        )
        np.testing.assert_array_equal(
            np.frombuffer(
                adapter.get_gene_expression("Gene_B"),
                dtype=np.float32,
            ),
            np.array([10.0, 20.0, 30.0], dtype=np.float32),
        )
