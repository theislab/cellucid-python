"""
:func:`cellucid.subsample_cells` follows the subsampling specification.

- Every case of the shared conformance fixture keeps exactly the cells it
  lists. The fixture is a byte-identical copy of the web app's
  ``tests/fixtures/subsample-conformance.json``, computed in plain Python by
  its ``scripts/generate-subsample-conformance.py`` independently of every
  implementation, and pinned here by its SHA-256 as in the web app and in
  cellucid-r.
- The cell-set digest and the session encoding of a cell set are the
  fixture's, and so is the priority order a subsample's cells are opened in.
- It takes an AnnData (a categorical obs column, an obsm key) as it takes
  arrays, keeps ``uniform`` and ``blocks`` nested, and refuses what the
  specification refuses.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import cellucid
from cellucid.subsample import (
    cell_set_digest,
    count_from_fraction,
    decode_cell_set,
    encode_cell_set,
)

FIXTURE = Path(__file__).parent / "fixtures" / "subsample_conformance.json"
FIXTURE_SHA256 = "a805a302a88e77e0c7e07e998811ce79825b5a34890db79a9462e801ab84a60d"
DOCUMENT = json.loads(FIXTURE.read_text("utf-8"))


def _run(case: dict) -> np.ndarray:
    method = case["method"]
    if method == "uniform":
        return cellucid.subsample_cells(case["nObs"], case["k"], seed=case["seed"])
    if method == "blocks":
        return cellucid.subsample_cells(case["nObs"], case["k"], method="blocks", seed=case["seed"], block_size=case["blockSize"])
    if method == "stratified":
        codes = np.asarray([-1 if code is None else code for code in case["codes"]])
        categories = [f"c{i}" for i in range(case["nCategories"])]
        groupby = pd.Categorical.from_codes(codes, categories=categories)
        return cellucid.subsample_cells(
            case["nObs"], case["k"], method="stratified", seed=case["seed"], groupby=groupby,
            allocation=case["allocation"], min_per_group=case["minPerCategory"],
        )
    values = np.asarray([np.nan if v is None else float(v) for v in case["coordinates"]], dtype=np.float32)
    coordinates = values.reshape(case["nObs"], case["dimensions"])
    return cellucid.subsample_cells(
        case["nObs"], case["k"], method="spatial", seed=case["seed"], basis=coordinates,
        columns=tuple(range(case["dimensions"])), grid_size=case["gridSize"],
    )


def test_the_fixture_is_the_pinned_one() -> None:
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256
    assert DOCUMENT["version"] == 1


def test_sizes_are_the_specifications() -> None:
    for item in DOCUMENT["refusals"]:
        with pytest.raises(ValueError, match="cannot hold"):
            _run(item)
    for item in DOCUMENT["fraction"]:
        assert count_from_fraction(item["fraction"], item["nObs"]) == item["k"]


@pytest.mark.parametrize("case", DOCUMENT["cases"], ids=lambda case: case["name"])
def test_every_case_keeps_the_fixtures_cells(case: dict) -> None:
    kept = _run(case)
    if case.get("kept", 0) is None:
        assert kept.tolist() == list(range(case["nObs"]))
        assert cell_set_digest(case["nObs"], None) == case["digest"]
        return
    if "kept" in case:
        assert kept.tolist() == case["kept"]
    else:
        assert kept.size == case["keptCount"]
        assert kept[:8].tolist() == case["keptHead"]
        assert kept[-8:].tolist() == case["keptTail"]
    assert cell_set_digest(case["nObs"], kept) == case["digest"]


def test_the_session_encoding_is_the_fixtures() -> None:
    for item in DOCUMENT["encoding"]:
        if item["kept"] is None:
            assert decode_cell_set(b"", item["nObs"]) is None
            assert cell_set_digest(item["nObs"], None) == item["digest"]
            continue
        assert encode_cell_set(item["kept"]).hex() == item["bytesHex"]
        assert decode_cell_set(bytes.fromhex(item["bytesHex"]), item["nObs"]).tolist() == item["kept"]
        assert cell_set_digest(item["nObs"], item["kept"]) == item["digest"]


def test_an_anndata_is_read_by_column_and_key() -> None:
    rng = np.random.default_rng(0)
    n = 400
    adata = ad.AnnData(
        obs=pd.DataFrame({"cluster": pd.Categorical(rng.choice(["a", "b", "c"], n))}, index=[f"c{i}" for i in range(n)]),
        obsm={"X_umap": rng.normal(size=(n, 2)).astype(np.float32)},
    )
    by_name = cellucid.subsample_cells(adata, 50, method="stratified", groupby="cluster", seed=3)
    by_values = cellucid.subsample_cells(n, 50, method="stratified", groupby=adata.obs["cluster"].values, seed=3)
    assert by_name.tolist() == by_values.tolist()
    # A text column by name: its values, made categories in code-point order
    # (as the web app and R take a text column), a missing value in none.
    text = rng.choice(np.array(["b", "a", "C", "\u00e9"], dtype=object), n)
    text[::17] = None
    adata.obs["label"] = text
    by_text = cellucid.subsample_cells(adata, 60, method="stratified", groupby="label", seed=5, allocation="balanced")
    ordered = pd.Categorical(text, categories=sorted({value for value in text if value is not None}))
    assert by_text.tolist() == cellucid.subsample_cells(n, 60, method="stratified", groupby=ordered, seed=5, allocation="balanced").tolist()
    spatial = cellucid.subsample_cells(adata, fraction=0.2, method="spatial", basis="X_umap", seed=1, grid_size=6)
    assert spatial.size == 80
    assert spatial.tolist() == cellucid.subsample_cells(n, 80, method="spatial", basis=adata.obsm["X_umap"], seed=1, grid_size=6).tolist()


def test_uniform_and_blocks_are_nested() -> None:
    for method, extra in (("uniform", {}), ("blocks", {"block_size": 37})):
        previous: set[int] | None = None
        for k in (3, 40, 900, 4999):
            kept = set(cellucid.subsample_cells(5000, k, method=method, seed=11, **extra).tolist())
            if previous is not None:
                assert previous <= kept
            previous = kept


def test_requests_outside_the_specification_are_refused() -> None:
    with pytest.raises(ValueError):
        cellucid.subsample_cells(10, 0)
    with pytest.raises(ValueError):
        cellucid.subsample_cells(10, 3, seed=-1)
    with pytest.raises(ValueError):
        cellucid.subsample_cells(10, 3, seed=2**32)
    with pytest.raises(ValueError):
        cellucid.subsample_cells(10, fraction=0)
    with pytest.raises(TypeError):
        cellucid.subsample_cells(10, 3, fraction=0.5)
    with pytest.raises(TypeError):
        cellucid.subsample_cells(10, 3, method="blocks")
    with pytest.raises(ValueError):
        cellucid.subsample_cells(4, 2, method="stratified", groupby=pd.Categorical(["a", "b", "c", "a"]), min_per_group=1)


def test_a_session_records_its_subsample_cells_as_the_web_app_writes_them() -> None:
    import base64

    from cellucid.session_bundle import _decode_cells

    kept = cellucid.subsample_cells(10_000, 1_234, method="blocks", seed=9, block_size=64)
    record = {
        "spec": {"method": "blocks", "nObs": 10_000, "k": 1_234, "seed": 9, "blockSize": 64},
        "kept": base64.b64encode(encode_cell_set(kept)).decode("ascii"),
        "digest": cell_set_digest(10_000, kept),
    }
    decoded = _decode_cells(record, 1_234)
    assert decoded["kept"].tolist() == kept.tolist()
    assert decoded["spec"]["blockSize"] == 64
    with pytest.raises(ValueError, match="digest"):
        _decode_cells({**record, "digest": "0" * 16}, 1_234)
    with pytest.raises(ValueError, match="cell count"):
        _decode_cells(record, 1_235)
    with pytest.raises(ValueError, match="exactly spec"):
        _decode_cells({**record, "extra": 1}, 1_234)
