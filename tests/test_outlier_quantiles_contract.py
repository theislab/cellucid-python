"""
The latent-space outlier quantiles :func:`cellucid.prepare` stores.

- The rule (``cellucid.prepare_data._outliers``) gives, bit for bit, the bytes
  of the shared conformance fixture, which the web repository's numpy
  generator (``scripts/generate-latent-outlier-conformance.py``) writes
  independently of every implementation; the web app and cellucid-r are
  checked against the same file (pinned here by its SHA-256).
- ``prepare(latent_key=...)`` stores ``uns['cellucid']['outlierQuantiles']``
  with every categorical obs column's quantiles in that latent space -- a
  categorical, a text column, a boolean and a nullable boolean column -- and
  none for a numeric column or one with no category of ``minPoints`` cells;
  without ``latent_key`` nothing is stored.
- Each field is its own float32 array, chunked and compressed like every
  array the viewer reads whole, so the viewer reads one field without the
  others.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import zarr

from cellucid import prepare, restore_original_order
from cellucid.prepare_data._outliers import latent_outlier_quantiles

FIXTURE = Path(__file__).parent / "fixtures" / "latent_outlier_conformance.json"
FIXTURE_SHA256 = "a812db4dd25c9b5b44549c001a0e3927043042ca12ad7d8567c73f583ecf498d"

N_CELLS = 60


def test_the_rule_gives_the_conformance_bytes() -> None:
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256
    for case in json.loads(FIXTURE.read_text(encoding="utf-8"))["cases"]:
        # None (JSON null): a NaN row's value, a cell the latent space does not place.
        stored = [np.nan if value is None else value for value in case["latent"]]
        latent = np.asarray(stored, dtype=np.float64).astype(np.float32).reshape(-1, case["width"])
        quantiles = latent_outlier_quantiles(latent, np.asarray(case["codes"]), case["nCategories"], case["minPoints"])
        assert quantiles.astype("<f4").tobytes().hex() == case["quantilesHex"], case["name"]


def _arguments() -> dict:
    rng = np.random.default_rng(11)
    obs = pd.DataFrame(
        {
            "cell_type": rng.choice(["T", "B", "NK"], N_CELLS),
            "rare": ["a"] * 2 + ["b"] * (N_CELLS - 2),
            "barcode": [f"AAAC{i:04d}" for i in range(N_CELLS)],
            "donor": pd.Categorical([None, *rng.choice(["d2", "d1"], N_CELLS - 1)], categories=["d2", "d1"]),
            "tiny": pd.Categorical(["x", "y", "z"] * (N_CELLS // 3)),
            "is_doublet": rng.random(N_CELLS) > 0.7,
            "flag": pd.array([True, True, None] * (N_CELLS // 3), dtype="boolean"),
            "score": rng.normal(size=N_CELLS),
        },
        index=[f"cell{i}" for i in range(N_CELLS)],
    )
    return {
        "embeddings": {
            "X_umap": rng.normal(size=(N_CELLS, 2)),
            "X_scVI": rng.normal(size=(N_CELLS, 10)) * 1e3 + 7,
        },
        "dataset_id": "outliers",
        "dataset_name": "Outliers",
        "obs": obs,
        "latent_key": "X_scVI",
        "centroid_min_points": 21,
    }


@pytest.mark.parametrize(
    ("name", "output_format"), [("data.h5ad", "h5ad"), ("data.zarr", "zarr"), ("data.zarr.zip", "zarr.zip")]
)
def test_prepare_stores_every_categorical_fields_quantiles(tmp_path: Path, name: str, output_format: str) -> None:
    path = prepare(tmp_path / name, format=output_format, **_arguments())
    if output_format == "h5ad":
        adata = ad.read_h5ad(path)
    elif output_format == "zarr":
        adata = ad.read_zarr(path)
    else:
        store = zarr.storage.ZipStore(path, mode="r")
        try:
            adata = ad.read_zarr(store)
        finally:
            store.close()
    stored = adata.uns["cellucid"]["outlierQuantiles"]
    assert stored["latentKey"] == "X_scVI"
    assert stored["minPoints"] == 21
    # 'rare' (58 of 60 in one category), every categorical, the boolean and the
    # nullable boolean; not the per-cell text, the 20-cell categories of
    # 'tiny', or the number.
    assert sorted(stored["fields"]) == sorted(["cell_type", "rare", "donor", "is_doublet", "flag"])
    latent = np.asarray(adata.obsm["X_scVI"], dtype=np.float32)
    for key, quantiles in stored["fields"].items():
        assert quantiles.dtype == np.float32 and quantiles.shape == (N_CELLS,)
        column = adata.obs[key]
        if isinstance(column.dtype, pd.CategoricalDtype):
            codes, n = column.cat.codes.to_numpy(), len(column.cat.categories)
        else:
            values = column.astype("boolean")
            codes = np.where(values.isna(), -1, values.fillna(False).astype(int)).astype(np.int64)
            n = 2
        expected = latent_outlier_quantiles(latent, codes, n, 21)
        assert quantiles.tobytes() == expected.tobytes(), key
    # The nullable boolean's missing values and the small category have none
    # (the cells by their original positions: the file stores them in
    # priority order).
    fields = restore_original_order(adata).uns["cellucid"]["outlierQuantiles"]["fields"]
    assert np.isnan(fields["flag"][2::3]).all()
    assert np.isnan(fields["rare"][:2]).all()
    assert not np.isnan(fields["rare"][2:]).any()


def test_each_field_is_its_own_array_read_whole(tmp_path: Path) -> None:
    path = prepare(tmp_path / "data.zarr", format="zarr", **_arguments())
    metadata = json.loads((path / ".zmetadata").read_text())["metadata"]
    array = metadata["uns/cellucid/outlierQuantiles/fields/cell_type/.zarray"]
    assert array["dtype"] == "<f4"
    assert array["shape"] == [N_CELLS]
    assert array["chunks"] == [N_CELLS]
    assert array["filters"] == [{"id": "shuffle", "elementsize": 4}]
    assert metadata["uns/cellucid/outlierQuantiles/.zattrs"]["encoding-type"] == "dict"
    assert metadata["uns/cellucid/outlierQuantiles/fields/.zattrs"]["encoding-type"] == "dict"


def test_without_a_latent_space_nothing_is_stored(tmp_path: Path) -> None:
    arguments = _arguments()
    del arguments["latent_key"]
    adata = ad.read_h5ad(prepare(tmp_path / "data.h5ad", format="h5ad", **arguments))
    assert "outlierQuantiles" not in adata.uns["cellucid"]


def test_the_default_fewest_cells_is_ten(tmp_path: Path) -> None:
    arguments = _arguments()
    del arguments["centroid_min_points"]
    adata = ad.read_h5ad(prepare(tmp_path / "data.h5ad", format="h5ad", **arguments))
    stored = adata.uns["cellucid"]["outlierQuantiles"]
    assert stored["minPoints"] == 10
    assert "tiny" in stored["fields"]
    assert "defaults" not in adata.uns["cellucid"] or "centroidMinPoints" not in adata.uns["cellucid"]["defaults"]
