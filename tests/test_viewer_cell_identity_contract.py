"""The Python side of a viewer names cells as the served store does.

The viewer names a cell by its original position. What Python serves -- an
object's snapshot, a file read where it is -- holds the cells in priority
order with each one's original position, and the names the Python side
attaches to cells (``viewer.obs_names``, ``selection.obs_names``, a highlight
or a column given by names, ``pull``, ``subset``) come from the same served
store, never from the live object in its own row order:

- an object read from a Cellucid-ready file in priority order (anndata,
  mudata) is named in its original order, and refused by ``pull``/``subset``
  until restored;
- an object changed after it was shown keeps the names it was shown with;
- a file served as it is is named at the version served, read again when it
  is written again;
- an object whose rows are the viewer's cells in another order is refused
  with repeated names too (the order is checked by sequence, and by the
  verified priority record).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from _serverless_viewer import serverless_viewer
from scipy import sparse

import cellucid
from cellucid.anndata_server import AnnDataServer

N = 300


def _adata(names: list[str] | None = None) -> ad.AnnData:
    rng = np.random.default_rng(7)
    adata = ad.AnnData(
        X=sparse.random(N, 12, density=0.3, format="csr", dtype=np.float32, random_state=2),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(rng.integers(0, 4, N).astype(str))},
            index=names if names is not None else [f"cell{i}" for i in range(N)],
        ),
        var=pd.DataFrame(index=[f"g{j}" for j in range(12)]),
    )
    adata.obsm["X_umap"] = rng.normal(size=(N, 2)).astype(np.float32)
    return adata


def _shown(data: Any, monkeypatch: pytest.MonkeyPatch) -> tuple[Any, list[dict[str, Any]]]:
    """A viewer of ``data`` as show_anndata makes one, its server made (not
    bound), its commands kept instead of sent."""
    viewer = serverless_viewer()
    viewer.data = data
    viewer._server = AnnDataServer(data, port=0, quiet=True, serve_web_ui=False, dataset_id="ready", dataset_name="Ready")
    viewer._served_cell_count = viewer._server.server_info["n_cells"]
    sent: list[dict[str, Any]] = []
    monkeypatch.setattr(viewer, "send_message", sent.append)
    return viewer, sent


def _cells(message: dict[str, Any]) -> list[int]:
    cells = message["cells"]
    return np.frombuffer(bytes(cells), dtype="<u4").tolist() if isinstance(cells, bytes | bytearray) else list(cells)


def test_an_object_read_from_a_file_in_priority_order_is_named_in_original_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = _adata()
    path = cellucid.prepare(tmp_path / "ready.zarr", format="zarr", data=original, dataset_id="ready", dataset_name="Ready")
    stored = ad.read_zarr(path)
    assert not stored.obs_names.equals(original.obs_names), "prepare() stores the cells in priority order"
    viewer, sent = _shown(stored, monkeypatch)
    try:
        assert viewer.obs_names.equals(original.obs_names)
        # A selection of the viewer's cells 0..4 names cell0..cell4.
        viewer._handle_frontend_message(
            {"type": "selection", "cells": bytearray(np.arange(5, dtype="<u4").tobytes()), "source": "lasso"}
        )
        assert list(viewer.selection.obs_names) == [f"cell{i}" for i in range(5)]
        # Cells named are the viewer's cells of those names.
        viewer.highlight_cells(pd.Index(["cell17", "cell3"]))
        assert sorted(_cells(sent[-1])) == [3, 17]
        # A column given by names lands on the cells of those names.
        viewer.push(pd.Series([1.0, 2.0], index=["cell17", "cell3"]), key="score")
        assert _cells(sent[-1]) == [3, 17] and np.asarray(sent[-1]["values"]).tolist() == [2.0, 1.0]
        # The object as stored is refused; restored, it is the viewer's cells in order.
        with pytest.raises(ValueError, match="restore_original_order"):
            viewer.subset(stored)
        with pytest.raises(ValueError, match="restore_original_order"):
            viewer.pull(stored)
        restored = cellucid.restore_original_order(stored)
        viewer._require_viewer_order(restored)
        assert viewer.subset(restored).obs_names.tolist() == [f"cell{i}" for i in range(5)]
    finally:
        viewer._server.stop()


def test_an_object_changed_after_it_was_shown_keeps_the_names_it_was_shown_with(monkeypatch: pytest.MonkeyPatch) -> None:
    adata = _adata()
    viewer, _ = _shown(adata, monkeypatch)
    try:
        adata.obs_names = [f"renamed{i}" for i in range(N)]
        assert viewer.obs_names.tolist() == [f"cell{i}" for i in range(N)]
        stopped = viewer._server
        viewer._server = None
        stopped.stop()
        # A viewer that stopped keeps the last names.
        assert viewer.obs_names.tolist() == [f"cell{i}" for i in range(N)]
    finally:
        if viewer._server is not None:
            viewer._server.stop()


def test_a_file_served_as_it_is_is_named_at_the_version_served(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = cellucid.prepare(tmp_path / "ready.h5ad", format="h5ad", data=_adata(), dataset_id="ready", dataset_name="Ready")
    viewer, _ = _shown(str(path), monkeypatch)
    try:
        assert viewer.obs_names.tolist() == [f"cell{i}" for i in range(N)]
        # Written again to the same path (prepare() again, other names): the
        # browser reads the new version, and so do the names.
        cellucid.prepare(
            path, format="h5ad", data=_adata([f"other{i}" for i in range(N)]), dataset_id="ready", dataset_name="Ready", force=True
        )
        assert viewer.obs_names.tolist() == [f"other{i}" for i in range(N)]
    finally:
        viewer._server.stop()


def test_a_mudata_object_is_served_from_a_snapshot_named_in_original_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    md = pytest.importorskip("mudata")
    rna, prot = _adata(), _adata()
    prot = prot[:, :4].copy()
    prot.var_names = [f"p{j}" for j in range(4)]
    mdata = md.MuData({"rna": rna, "prot": prot})
    mdata.obsm["X_umap"] = rna.obsm["X_umap"]
    path = cellucid.prepare(tmp_path / "ready.zarr", format="zarr", data=mdata, dataset_id="ready", dataset_name="Ready")
    stored = md.read_zarr(path)
    assert not stored.obs_names.equals(mdata.obs_names)
    viewer, _ = _shown(stored, monkeypatch)
    try:
        assert viewer.obs_names.equals(mdata.obs_names)
        store = viewer._server._served.priority
        assert store is not None and store.layout == "mudata"
        store.laid_out()
        served = md.read_zarr(_dumped(store, tmp_path / "served.zarr"))
        # Served in priority order, each cell by its original position, every
        # modality's values the cells' own.
        positions = served.obs["cellucid_original_position"].to_numpy().astype(np.int64)
        assert served.obs_names.tolist() == [mdata.obs_names[p] for p in positions]
        for name in ("rna", "prot"):
            np.testing.assert_array_equal(served.mod[name].X.toarray(), mdata.mod[name].X.toarray()[positions])
        with pytest.raises(ValueError, match="restore_original_order"):
            viewer.pull(stored)
    finally:
        viewer._server.stop()


def _dumped(store: Any, folder: Path) -> Path:
    """Every key of a served store, written as a Zarr v2 folder."""
    for key in store.stored_keys():
        value = store.get(key)
        if value is None:
            continue
        target = folder / key
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(value)
    return folder


def test_an_object_reordered_is_refused_with_repeated_names_too(monkeypatch: pytest.MonkeyPatch) -> None:
    names = [f"cell{i // 2}" for i in range(N)]  # every name twice
    adata = _adata(names)
    viewer, _ = _shown(adata, monkeypatch)
    try:
        viewer._require_viewer_order(adata)
        reordered = adata[np.r_[np.arange(2, N), [0, 1]]].copy()
        with pytest.raises(ValueError, match="another order"):
            viewer._require_viewer_order(reordered)
    finally:
        viewer._server.stop()
