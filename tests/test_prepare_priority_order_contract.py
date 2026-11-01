"""The priority order ``prepare()`` stores cells in, and ``restore_original_order()``.

The written contract is the Cellucid-ready AnnData specification
(``docs/user_guide/python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification.md``,
"Priority order"): rows by ascending ``uniform`` rank under the file's seed,
``obs['cellucid_original_position']``, ``uns['cellucid']['priority']`` and,
beside a sparse ``X``, ``varm['cellucid_prefix_counts']``. cellucid-r's
``tests/testthat/test-priority-order.R`` pins the same on its writer, and its
cross-writer test holds the two writers to one store.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import warnings
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import zarr
from scipy import sparse

from cellucid import prepare, restore_original_order, subsample_cells
from cellucid.prepare_data._priority import (
    LEVEL_RATIO,
    prefix_levels,
    stored_order,
    stored_order_of,
)
from cellucid.subsample import _ranks

REPOSITORY = Path(__file__).resolve().parents[1]
WEB_ROOT = REPOSITORY.parent / "cellucid"
N_CELLS = 300
N_GENES = 12

OUTPUTS = [
    pytest.param("data.h5ad", "h5ad", None, id="h5ad"),
    pytest.param("data.zarr", "zarr", 2, id="zarr-v2"),
    pytest.param("data.zarr", "zarr", 3, id="zarr-v3"),
    pytest.param("data.zarr.zip", "zarr.zip", 2, id="zarr-zip"),
]


def _inputs() -> dict[str, Any]:
    rng = np.random.default_rng(11)
    counts = sparse.random(N_CELLS, N_GENES, density=0.3, format="csc", random_state=4, dtype=np.float32)
    counts.data = rng.integers(1, 9, counts.nnz).astype(np.float32)
    graph = sparse.random(N_CELLS, N_CELLS, density=0.02, format="csr", random_state=5, dtype=np.float32)
    graph = (graph + graph.T).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(rng.choice(["B", "T", "NK"], N_CELLS)),
            "barcode": [f"BC{i:04d}" for i in range(N_CELLS)],
            "score": rng.normal(size=N_CELLS),
            "batch": pd.array([None if i % 7 == 0 else i % 3 for i in range(N_CELLS)], dtype="Int64"),
        },
        index=[f"cell{i}" for i in range(N_CELLS)],
    )
    return {
        "embeddings": {"X_umap": rng.normal(size=(N_CELLS, 2)), "X_pca": rng.normal(size=(N_CELLS, 5))},
        "obs": obs,
        "var": pd.DataFrame(index=[f"g{i}" for i in range(N_GENES)]),
        "gene_expression": counts,
        "connectivities": graph,
        "dataset_id": "priority",
        "dataset_name": "Priority",
        "latent_key": "X_pca",
    }


def _stored(path: Path) -> ad.AnnData:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        if path.suffix == ".h5ad":
            return ad.read_h5ad(path)
        if path.name.endswith(".zip"):
            store = zarr.storage.ZipStore(path, mode="r")
            try:
                return ad.read_zarr(store)
            finally:
                store.close()
        return ad.read_zarr(path)


@pytest.mark.parametrize(("name", "output_format", "zarr_format"), OUTPUTS)
@pytest.mark.parametrize("seed", [0, 7, 2**32 - 1])
def test_cells_are_stored_by_ascending_rank_with_every_per_cell_element(
    tmp_path: Path, name: str, output_format: str, zarr_format: int | None, seed: int
) -> None:
    given = _inputs()
    arguments = {**given, "priority_seed": seed}
    if zarr_format is not None:
        arguments["zarr_format"] = zarr_format
    stored = _stored(prepare(tmp_path / name, format=output_format, **arguments))

    original = stored.obs["cellucid_original_position"].to_numpy()
    assert original.dtype == np.uint32
    assert sorted(original.tolist()) == list(range(N_CELLS))
    assert np.all(np.diff(_ranks(original, seed).astype(np.int64)) > 0)
    # Every uniform subsample under the file's seed is a prefix of the rows.
    for k in (1, 17, 150, N_CELLS - 1):
        kept = subsample_cells(N_CELLS, n_cells=k, seed=seed)
        assert sorted(original[:k].tolist()) == kept.tolist()

    rows = original.astype(np.int64)
    assert stored.obs_names.tolist() == given["obs"].index[rows].tolist()
    assert stored.obs["barcode"].tolist() == given["obs"]["barcode"].iloc[rows].tolist()
    assert stored.obs["batch"].tolist() == given["obs"]["batch"].iloc[rows].tolist()
    np.testing.assert_array_equal(stored.obsm["X_umap"], given["embeddings"]["X_umap"][rows].astype(np.float32))
    assert (given["gene_expression"][rows] != stored.X).nnz == 0
    x = stored.X
    for gene in range(N_GENES):
        assert np.all(np.diff(x.indices[x.indptr[gene] : x.indptr[gene + 1]]) > 0)
    assert (stored.obsp["connectivities"] != given["connectivities"][rows][:, rows]).nnz == 0

    priority = stored.uns["cellucid"]["priority"]
    assert int(priority["seed"]) == seed and int(priority["specVersion"]) == 1
    levels = np.asarray(priority["levels"])
    np.testing.assert_array_equal(levels, prefix_levels(N_CELLS))
    counts = stored.varm["cellucid_prefix_counts"]
    assert counts.dtype == np.uint32 and counts.shape == (N_GENES, len(levels))
    for gene in range(N_GENES):
        np.testing.assert_array_equal(counts[gene], np.searchsorted(x.indices[x.indptr[gene] : x.indptr[gene + 1]], levels))

    # Back in the order given: everything as it was.
    restored = restore_original_order(stored)
    assert restored.obs_names.tolist() == given["obs"].index.tolist()
    assert "cellucid_original_position" not in restored.obs
    assert "priority" not in restored.uns["cellucid"] and "cellucid_prefix_counts" not in restored.varm
    assert (given["gene_expression"] != restored.X).nnz == 0
    assert (restored.obsp["connectivities"] != given["connectivities"]).nnz == 0
    np.testing.assert_array_equal(restored.obsm["X_pca"], given["embeddings"]["X_pca"].astype(np.float32))
    # The outlier quantiles follow the cells: as prepare() writes them for the
    # cells in their original order, up to the order the sums ran in.
    written_as_given = _stored(
        prepare(tmp_path / f"given-{name}", format=output_format, **{**arguments, "priority_seed": None})
    )
    for key, values in restored.uns["cellucid"]["outlierQuantiles"]["fields"].items():
        np.testing.assert_allclose(values, written_as_given.uns["cellucid"]["outlierQuantiles"]["fields"][key], rtol=1e-6)


def test_the_prefix_levels_step_down_by_three_quarters_as_cellucid_r_counts_them() -> None:
    assert LEVEL_RATIO == (3, 4)
    assert prefix_levels(1).tolist() == []
    assert prefix_levels(4).tolist() == [1, 2, 3]
    big = prefix_levels(2**32 - 1)
    assert big[-1] == (2**32 - 1) * 3 // 4 and np.all(np.diff(big) > 0)
    # The rule written out: each level floor(previous * 3 / 4), repeats kept once.
    level, expected = 10_000, []
    while (level := level * 3 // 4) >= 1:
        if not expected or expected[-1] != level:
            expected.append(level)
    assert prefix_levels(10_000).tolist() == expected[::-1]


def test_priority_seed_none_stores_the_cells_as_given_and_a_bad_seed_is_refused(tmp_path: Path) -> None:
    given = _inputs()
    stored = _stored(prepare(tmp_path / "as-given.zarr", format="zarr", priority_seed=None, **given))
    assert stored.obs_names.tolist() == given["obs"].index.tolist()
    assert "cellucid_original_position" not in stored.obs
    assert "priority" not in stored.uns["cellucid"] and "cellucid_prefix_counts" not in stored.varm
    assert restore_original_order(stored) is stored
    for seed in (-1, 2**32, 1.5, True, "0"):
        with pytest.raises((TypeError, ValueError), match="priority_seed"):
            prepare(tmp_path / f"bad-{seed!s}.zarr", format="zarr", priority_seed=seed, **given)


def test_a_file_prepared_again_keeps_its_cells_original_positions(tmp_path: Path) -> None:
    first = _stored(prepare(tmp_path / "first.zarr", format="zarr", priority_seed=3, **_inputs()))
    again = _stored(prepare(tmp_path / "again.zarr", format="zarr", data=first, priority_seed=3))
    assert again.obs_names.tolist() == first.obs_names.tolist()
    np.testing.assert_array_equal(again.obs["cellucid_original_position"], first.obs["cellucid_original_position"])
    # A column that is not the record's order: the record does not describe
    # the cells, so they are taken as given (never refused, never misread).
    broken = first.copy()
    broken.obs["cellucid_original_position"] = np.zeros(N_CELLS, dtype=np.uint32)
    with pytest.warns(UserWarning, match="is not its seed's order"):
        rewritten = _stored_quiet(prepare(tmp_path / "broken.zarr", format="zarr", data=broken, priority_seed=None))
    assert rewritten.obs_names.tolist() == first.obs_names.tolist()
    assert "cellucid_original_position" not in rewritten.obs
    with pytest.warns(UserWarning, match="is not its seed's order"):
        assert restore_original_order(broken) is broken


def _stored_quiet(path: Path) -> ad.AnnData:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ad.read_zarr(path)


def test_the_record_names_its_cell_count_and_is_verified_before_it_is_used(tmp_path: Path) -> None:
    """``uns['cellucid']['priority']`` records ``nObs``; an object whose record
    no longer describes it -- cells added or removed, reordered, a newer
    specification, other levels -- is taken as cells Cellucid did not write
    (each row its own original position), with a warning: the rule the web
    app reads such a file by (``priority-layout.js``)."""
    stored = _stored(prepare(tmp_path / "v.zarr", format="zarr", priority_seed=7, **_inputs()))
    record = stored.uns["cellucid"]["priority"]
    assert int(record["nObs"]) == N_CELLS and int(record["seed"]) == 7
    np.testing.assert_array_equal(record["levels"], prefix_levels(N_CELLS))
    verified = stored_order_of(stored)
    assert verified.notice is None
    np.testing.assert_array_equal(verified.positions, stored.obs["cellucid_original_position"])
    # The column is a check, not the source: a record without it verifies.
    np.testing.assert_array_equal(stored_order(record, N_CELLS, None).positions, verified.positions)
    column = stored.obs["cellucid_original_position"].to_numpy()
    swapped = column.copy()
    swapped[[3, 4]] = swapped[[4, 3]]
    cases = [
        (dict(record), N_CELLS - 1, column[:-1], "written for 300 cells, and there are 299"),
        (dict(record), N_CELLS, swapped, "is not its seed's order"),
        ({**record, "specVersion": 2}, N_CELLS, column, "version 2 of the subsampling specification"),
        ({key: value for key, value in record.items() if key != "nObs"}, N_CELLS, column, "records no cell count"),
        ({**record, "levels": prefix_levels(2 * N_CELLS)}, N_CELLS, column, "levels are not those of 300 cells"),
        ({**record, "seed": -1}, N_CELLS, column, "seed is not a whole number"),
        (None, N_CELLS, column, "has no uns\\['cellucid'\\]\\['priority'\\] beside it"),
    ]
    for given, n_obs, positions, why in cases:
        order = stored_order(given, n_obs, positions)
        assert order.positions is None, why
        assert order.notice is not None and __import__("re").search(why, order.notice), (why, order.notice)
    assert stored_order(None, N_CELLS, None) == stored_order(None, N_CELLS, None)
    assert stored_order(None, N_CELLS, None).notice is None

    # Objects changed after they were read: restore leaves them as they are.
    dropped = stored[np.arange(N_CELLS) != 150].copy()
    with pytest.warns(UserWarning, match="written for 300 cells, and there are 299"):
        assert restore_original_order(dropped) is dropped
    resorted = stored[np.argsort(stored.obs_names.to_numpy())].copy()
    with pytest.warns(UserWarning, match="is not its seed's order"):
        assert restore_original_order(resorted) is resorted
    # A record whose column was dropped still restores, from its seed.
    no_column = stored.copy()
    del no_column.obs["cellucid_original_position"]
    assert restore_original_order(no_column).obs_names.tolist() == [f"cell{i}" for i in range(N_CELLS)]


def test_restore_reads_a_backed_anndata_into_memory_in_original_order(tmp_path: Path) -> None:
    given = _inputs()
    path = prepare(tmp_path / "b.h5ad", format="h5ad", **given)
    backed = ad.read_h5ad(path, backed="r")
    try:
        restored = restore_original_order(backed)
    finally:
        backed.file.close()
    assert not restored.isbacked
    assert restored.obs_names.tolist() == given["obs"].index.tolist()
    np.testing.assert_array_equal(restored.X.toarray(), given["gene_expression"].toarray())
    np.testing.assert_array_equal(restored.obsm["X_umap"], given["embeddings"]["X_umap"].astype(np.float32))
    assert "cellucid_original_position" not in restored.obs and "priority" not in restored.uns["cellucid"]


def test_prepare_without_a_seed_puts_a_stored_order_back_first(tmp_path: Path) -> None:
    """An object read in its stored priority order, prepared with
    ``priority_seed=None``, is written in its original order: rows are then
    original positions, as the viewer names them, and no stale position
    column or record is left."""
    given = _inputs()
    stored = _stored(prepare(tmp_path / "p.zarr", format="zarr", **given))
    as_given = _stored(prepare(tmp_path / "g.zarr", format="zarr", data=stored, priority_seed=None))
    assert as_given.obs_names.tolist() == given["obs"].index.tolist()
    assert "cellucid_original_position" not in as_given.obs
    assert "priority" not in as_given.uns["cellucid"]
    np.testing.assert_array_equal(as_given.obsm["X_umap"], given["embeddings"]["X_umap"].astype(np.float32))
    # obs= given alone (no record beside the column): the column is not written.
    with pytest.warns(UserWarning, match="is not written"):
        plain = _stored_quiet(prepare(tmp_path / "o.zarr", format="zarr", priority_seed=None,
                                      **{**given, "obs": stored.obs, "embeddings": {"X_umap": stored.obsm["X_umap"]},
                                         "gene_expression": stored.X, "connectivities": None, "latent_key": None}))
    assert "cellucid_original_position" not in plain.obs


def test_a_mudata_orders_each_modality_by_its_cells_priority(tmp_path: Path) -> None:
    md = pytest.importorskip("mudata")
    rng = np.random.default_rng(2)
    names = [f"c{i}" for i in range(400)]
    rna = ad.AnnData(
        X=sparse.random(400, 20, density=0.2, format="csr", random_state=3, dtype=np.float32),
        obs=pd.DataFrame(index=names),
        var=pd.DataFrame(index=[f"g{i}" for i in range(20)]),
    )
    rna.obsm["X_umap"] = rng.normal(size=(400, 2)).astype(np.float32)
    measured = rng.permutation(400)[:250]
    prot = ad.AnnData(
        X=sparse.random(250, 6, density=0.5, format="csr", random_state=6, dtype=np.float32),
        obs=pd.DataFrame(index=[names[i] for i in measured]),
        var=pd.DataFrame(index=[f"p{i}" for i in range(6)]),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mdata = md.MuData({"rna": rna, "prot": prot})
        path = prepare(tmp_path / "m.h5mu", format="h5mu", data=mdata, dataset_id="m", dataset_name="M")
        stored = md.read_h5mu(path)
    levels = np.asarray(stored.uns["cellucid"]["priority"]["levels"])
    for name in ("rna", "prot"):
        obsmap = np.asarray(stored.obsmap[name]).astype(np.int64)
        rows = obsmap[obsmap > 0]
        assert np.all(np.diff(rows) > 0), name
        below = np.concatenate([[0], np.cumsum(obsmap > 0)])[levels]
        x = sparse.csc_matrix(stored.mod[name].X)
        counts = stored.mod[name].varm["cellucid_prefix_counts"]
        for gene in range(x.shape[1]):
            np.testing.assert_array_equal(counts[gene], np.searchsorted(x.indices[x.indptr[gene] : x.indptr[gene + 1]], below))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        restored = restore_original_order(stored)
    assert restored.obs_names.tolist() == names
    assert restored.mod["prot"].obs_names.tolist() == prot.obs_names.tolist()
    np.testing.assert_array_equal(restored.mod["prot"].X.toarray(), prot.X.toarray())
    np.testing.assert_array_equal(restored.mod["rna"].X.toarray(), rna.X.toarray())


_WEB_PREFIX = r"""
const [web, target, k] = process.argv.slice(1);
const { readFileSync, readdirSync, statSync } = await import('node:fs');
const path = await import('node:path');
const { pathToFileURL } = await import('node:url');
const load = rel => import(pathToFileURL(path.join(web, rel)).href);
const { openDatasetStore } = await load('assets/js/data/anndata/open.js');
const { AnnDataReader } = await load('assets/js/data/anndata/reader.js');
console.log = () => {};
const name = path.basename(target);
let source;
if (statSync(target).isDirectory()) {
  const files = new Map();
  const walk = dir => {
    for (const entry of readdirSync(dir)) {
      const full = path.join(dir, entry);
      if (statSync(full).isDirectory()) walk(full);
      else files.set(path.relative(target, full).split(path.sep).join('/'), new File([readFileSync(full)], entry));
    }
  };
  walk(target);
  source = { files, rootName: name };
} else {
  source = { file: new File([readFileSync(target)], name) };
}
const read = async cellSubset => {
  const opened = await openDatasetStore({ source }, { decodedChunkCacheBytes: 64 << 20 });
  const reader = new AnnDataReader();
  await reader.open(opened.root, { format: opened.readerFormat, name: opened.name, decodedChunkBytes: 64 << 20, link: opened.link, cellSubset });
  const genes = await reader.getMatrixColumns('X', Array.from({ length: reader.nVars }, (_, g) => g));
  const out = {
    nObs: reader.nObs,
    prefix: reader.cellSubset?.prefix ?? null,
    kept: reader.cellSubset?.kept ? Array.from(reader.cellSubset.kept) : null,
    genes: genes.map(column => Array.from(column)),
    umap: Array.from(await reader.readObsmColumns('X_umap', [0, 1])),
    cellType: Array.from((await reader.getObsField('cell_type')).codes),
  };
  reader.close();
  return out;
};
const all = await read(null);
const part = await read({ mode: 'subsample', method: 'uniform', k: Number(k) });
process.stdout.write(JSON.stringify({ all, part }));
"""


@pytest.mark.parametrize(("name", "output_format", "zarr_format"), OUTPUTS)
def test_the_web_app_opens_a_uniform_subsample_as_the_prefix_of_the_rows(
    tmp_path: Path, name: str, output_format: str, zarr_format: int | None
) -> None:
    node = shutil.which("node")
    if node is None or not (WEB_ROOT / "assets/js/data/anndata/priority-layout.js").is_file():
        pytest.skip("Node and the web app (with its priority reader) are needed")
    arguments = {**_inputs(), "priority_seed": 0}
    if zarr_format is not None:
        arguments["zarr_format"] = zarr_format
    path = prepare(tmp_path / name, format=output_format, **arguments)
    k = 120
    completed = subprocess.run(
        [node, "--input-type=module", "-e", _WEB_PREFIX, str(WEB_ROOT), str(path), str(k)],
        capture_output=True, text=True, check=False, timeout=300,
    )
    assert completed.returncode == 0, completed.stderr
    read = json.loads(completed.stdout)
    stored = _stored(path)
    original = stored.obs["cellucid_original_position"].to_numpy()
    part, whole = read["part"], read["all"]
    assert part["nObs"] == k and part["prefix"] == k
    # The cells kept are the uniform subsample's, named by original position.
    assert part["kept"] == sorted(original[:k].tolist())
    # Every value is the full read's at the first k rows.
    for gene in range(N_GENES):
        assert part["genes"][gene] == whole["genes"][gene][:k]
    assert part["umap"] == whole["umap"][: 2 * k]
    assert part["cellType"] == whole["cellType"][:k]


def test_a_session_is_applied_in_original_order_and_refuses_a_stored_order(tmp_path: Path) -> None:
    """A session names cells by original position: applied to a file read back
    in its stored order it refuses, and points at ``restore_original_order``."""
    from cellucid import apply_cellucid_session_to_anndata

    stored = _stored(prepare(tmp_path / "s.zarr", format="zarr", **_inputs()))
    fixture = Path(__file__).parent / "fixtures" / "session_format3" / "session.cellucid-session"
    with pytest.raises(ValueError, match="restore_original_order"):
        apply_cellucid_session_to_anndata(fixture, stored, expected_dataset_id="format3")
    # A record that no longer describes the object (rows reordered since):
    # each row is its own original position, as the web app read the file;
    # the warning says why, and the session's own checks apply.
    resorted = stored[np.argsort(stored.obs_names.to_numpy())[::-1]].copy()
    with pytest.warns(UserWarning, match="is not its seed's order"), pytest.raises(ValueError) as refused:
        apply_cellucid_session_to_anndata(fixture, resorted, expected_dataset_id="format3")
    assert "restore_original_order" not in str(refused.value)


_WEB_MUDATA_PREFIX = _WEB_PREFIX.replace(
    """  const genes = await reader.getMatrixColumns('X', Array.from({ length: reader.nVars }, (_, g) => g));
  const out = {
    nObs: reader.nObs,
    prefix: reader.cellSubset?.prefix ?? null,
    kept: reader.cellSubset?.kept ? Array.from(reader.cellSubset.kept) : null,
    genes: genes.map(column => Array.from(column)),
    umap: Array.from(await reader.readObsmColumns('X_umap', [0, 1])),
    cellType: Array.from((await reader.getObsField('cell_type')).codes),
  };""",
    """  const out = { nObs: reader.nObs, prefix: reader.cellSubset?.prefix ?? null, genes: {} };
  for (const name of ['rna', 'prot']) {
    const { index } = await reader.readVarNames(`mod/${name}`, null);
    const columns = await reader.getMatrixColumns(`mod/${name}/X`, index.map((_, g) => g));
    out.genes[name] = columns.map(column => Array.from(column, v => (Number.isNaN(v) ? null : v)));
  }""",
)


def test_the_web_app_opens_a_mudata_prefix_with_every_modalitys_first_rows(tmp_path: Path) -> None:
    md = pytest.importorskip("mudata")
    node = shutil.which("node")
    if node is None or not (WEB_ROOT / "assets/js/data/anndata/priority-layout.js").is_file():
        pytest.skip("Node and the web app (with its priority reader) are needed")
    rng = np.random.default_rng(8)
    names = [f"c{i}" for i in range(500)]
    rna = ad.AnnData(
        X=sparse.random(500, 8, density=0.3, format="csr", random_state=1, dtype=np.float32),
        obs=pd.DataFrame(index=names),
        var=pd.DataFrame(index=[f"g{i}" for i in range(8)]),
    )
    rna.obsm["X_umap"] = rng.normal(size=(500, 2)).astype(np.float32)
    measured = rng.permutation(500)[:320]
    prot = ad.AnnData(
        X=sparse.random(320, 4, density=0.5, format="csr", random_state=2, dtype=np.float32),
        obs=pd.DataFrame(index=[names[i] for i in measured]),
        var=pd.DataFrame(index=[f"p{i}" for i in range(4)]),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = prepare(tmp_path / "m.h5mu", format="h5mu", data=md.MuData({"rna": rna, "prot": prot}),
                       dataset_id="m", dataset_name="M")
    k = 90
    completed = subprocess.run(
        [node, "--input-type=module", "-e", _WEB_MUDATA_PREFIX, str(WEB_ROOT), str(path), str(k)],
        capture_output=True, text=True, check=False, timeout=300,
    )
    assert completed.returncode == 0, completed.stderr
    read = json.loads(completed.stdout)
    part, whole = read["part"], read["all"]
    assert part["prefix"] == k
    for name in ("rna", "prot"):
        for gene, column in enumerate(part["genes"][name]):
            assert column == whole["genes"][name][gene][:k], (name, gene)
