#!/usr/bin/env python3
"""Write the files the "open a format" captures pick, from public samples only.

Usage::

    python build-format-fixtures.py OUT_DIR

``capture.mjs`` runs this on demand (``formatFixture()`` in ``scenarios.mjs``)
with ``$CELLUCID_DOCS_PYTHON``, into a cache folder outside every repository,
so nothing here is committed but the recipe. Every file is derived from data
the workspace already ships:

- the Pancreas sample of ``cellucid-datasets/anndata/pancreas.zarr`` (scVelo's
  public E15.5 endocrinogenesis data, integer counts in ``X``), written as
  each tool writes it: an ``.h5ad`` (anndata), a Loom file, a Cell Ranger
  ``outs/`` folder, a MatrixMarket folder with a cell table, and -- through
  ``build-format-fixtures.R`` with ``$CELLUCID_DOCS_RSCRIPT`` -- a Seurat v5
  ``.rds`` and an ``.rda`` holding two Seurat objects;
- the CITE-seq MuData of the web repository's reader corpus
  (``cellucid/tests/fixtures/mudata/mudata-0.4.1.zip``,
  ``citeseq-subset.h5mu``: two modalities, proteins measured on 36 of 44
  cells), as it is.

The files stand in for the reader's own; the pictures show what Cellucid
does with each format, on cells a reader of the documentation recognises.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import warnings
import zipfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[4]
DATASETS = Path(os.environ.get("CELLUCID_DATASETS_REPO", WORKSPACE / "cellucid-datasets"))
WEB = Path(os.environ.get("CELLUCID_WEB_REPO", WORKSPACE / "cellucid"))
OBS = ["clusters_coarse", "clusters", "S_score", "G2M_score"]


def pancreas() -> ad.AnnData:
    adata = ad.read_zarr(DATASETS / "anndata" / "pancreas.zarr")
    adata.X = sp.csr_matrix(adata.X, dtype=np.float32)
    # The store's own record of its cells' order (prepare()'s position column,
    # beside uns['cellucid']['priority']) is not part of a file another tool
    # writes, and alone it names no order: the app would say so on every file.
    adata.obs = adata.obs.drop(columns=["cellucid_original_position"], errors="ignore")
    return adata


def write_loom(adata: ad.AnnData, target: Path) -> None:
    """Loom 3.0.0 as loompy writes it: genes x cells, obsm as cell attributes."""
    with h5py.File(target, "w") as h5:
        h5.attrs["LOOM_SPEC_VERSION"] = "3.0.0"
        h5.create_dataset(
            "matrix", data=adata.X.T.toarray(), chunks=(64, 64), compression="gzip"
        )
        h5.create_group("layers")
        rows = h5.create_group("row_attrs")
        rows.create_dataset("Gene", data=np.array(adata.var_names, dtype="S"))
        cols = h5.create_group("col_attrs")
        cols.create_dataset("CellID", data=np.array(adata.obs_names, dtype="S"))
        for key in OBS:
            values = adata.obs[key]
            if isinstance(values.dtype, pd.CategoricalDtype):
                cols.create_dataset(key, data=np.array(values.astype(str), dtype="S"))
            else:
                cols.create_dataset(key, data=values.to_numpy(np.float64))
        cols.create_dataset("X_umap", data=adata.obsm["X_umap"].astype(np.float32))
        cols.create_dataset("X_pca", data=adata.obsm["X_pca"][:, :30].astype(np.float32))
        h5.create_group("col_graphs")
        h5.create_group("row_graphs")
        attrs = h5.create_group("attrs")
        attrs.create_dataset("title", data=np.bytes_("Pancreatic endocrinogenesis (Loom)"))


def write_tenx_h5(adata: ad.AnnData, target: Path) -> None:
    """A Cell Ranger 3 ``filtered_feature_bc_matrix.h5`` (CSC, genes x cells)."""
    matrix = sp.csc_matrix(adata.X.T.astype(np.int32))
    barcodes = [f"{name}-1" for name in adata.obs_names]
    with h5py.File(target, "w") as h5:
        group = h5.create_group("matrix")
        group.create_dataset("data", data=matrix.data.astype(np.int32), compression="gzip")
        group.create_dataset("indices", data=matrix.indices.astype(np.int64), compression="gzip")
        group.create_dataset("indptr", data=matrix.indptr.astype(np.int64))
        group.create_dataset("shape", data=np.array(matrix.shape, dtype=np.int32))
        group.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        features = group.create_group("features")
        names = np.array(adata.var_names, dtype="S")
        features.create_dataset("id", data=names)
        features.create_dataset("name", data=names)
        features.create_dataset(
            "feature_type", data=np.array([b"Gene Expression"] * adata.n_vars)
        )
        features.create_dataset("genome", data=np.array([b"mm10"] * adata.n_vars))
        features.create_dataset("_all_tag_keys", data=np.array([b"genome"]))


def write_outs(adata: ad.AnnData, outs: Path) -> None:
    """A Cell Ranger ``outs/`` folder: the matrix plus ``analysis/``."""
    outs.mkdir(parents=True)
    write_tenx_h5(adata, outs / "filtered_feature_bc_matrix.h5")
    barcodes = [f"{name}-1" for name in adata.obs_names]
    umap = outs / "analysis" / "umap" / "gene_expression_2_components"
    umap.mkdir(parents=True)
    pd.DataFrame(
        {
            "Barcode": barcodes,
            "UMAP-1": adata.obsm["X_umap"][:, 0],
            "UMAP-2": adata.obsm["X_umap"][:, 1],
        }
    ).to_csv(umap / "projection.csv", index=False)
    clustering = outs / "analysis" / "clustering" / "gene_expression_clusters"
    clustering.mkdir(parents=True)
    pd.DataFrame(
        {"Barcode": barcodes, "Cluster": adata.obs["clusters"].astype(str).to_numpy()}
    ).to_csv(clustering / "clusters.csv", index=False)


def write_matrix_market(adata: ad.AnnData, folder: Path) -> None:
    """Cell Ranger 3 MatrixMarket files plus a cell table with a UMAP."""
    folder.mkdir(parents=True)
    with gzip.open(folder / "matrix.mtx.gz", "wb") as handle:
        scipy.io.mmwrite(handle, sp.coo_matrix(adata.X.T.astype(np.int32)), field="integer")
    names = list(adata.var_names)
    with gzip.open(folder / "features.tsv.gz", "wt") as handle:
        for name in names:
            handle.write(f"{name}\t{name}\tGene Expression\n")
    with gzip.open(folder / "barcodes.tsv.gz", "wt") as handle:
        for name in adata.obs_names:
            handle.write(f"{name}-1\n")
    table = adata.obs[OBS].copy()
    table.insert(0, "barcode", [f"{name}-1" for name in adata.obs_names])
    table["UMAP_1"] = adata.obsm["X_umap"][:, 0]
    table["UMAP_2"] = adata.obsm["X_umap"][:, 1]
    table.to_csv(folder / "cell_metadata.csv", index=False)


def export_for_r(adata: ad.AnnData, folder: Path) -> None:
    folder.mkdir(parents=True)
    scipy.io.mmwrite(str(folder / "X.mtx"), sp.coo_matrix(adata.X.T))
    (folder / "cells.txt").write_text("\n".join(adata.obs_names) + "\n")
    (folder / "genes.txt").write_text("\n".join(adata.var_names) + "\n")
    adata.obs[OBS].to_csv(folder / "obs.csv")
    pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names).to_csv(folder / "umap.csv")
    pd.DataFrame(adata.obsm["X_pca"][:, :30], index=adata.obs_names).to_csv(folder / "pca.csv")


def main() -> None:
    out = Path(sys.argv[1]).resolve()
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    adata = pancreas()
    plain = adata.copy()
    plain.uns = {"title": "Pancreatic endocrinogenesis"}
    plain.write_h5ad(out / "pancreas.h5ad", compression="gzip")
    # The published store itself, zipped uncompressed as `.zarr.zip` is.
    store = DATASETS / "anndata" / "pancreas.zarr"
    with zipfile.ZipFile(out / "pancreas.zarr.zip", "w", zipfile.ZIP_STORED) as archive:
        for path in sorted(store.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(store).as_posix())
    write_loom(adata, out / "pancreas.loom")
    write_outs(adata, out / "outs")
    write_matrix_market(adata, out / "pancreas_matrix")

    with zipfile.ZipFile(WEB / "tests" / "fixtures" / "mudata" / "mudata-0.4.1.zip") as corpus:
        (out / "citeseq-subset.h5mu").write_bytes(corpus.read("citeseq-subset.h5mu"))

    export = out / ".r-export"
    export_for_r(adata, export)
    rscript = os.environ.get("CELLUCID_DOCS_RSCRIPT", "Rscript")
    subprocess.run(
        [rscript, str(HERE / "build-format-fixtures.R"), str(export), str(out)],
        check=True,
    )
    shutil.rmtree(export)
    for path in sorted(out.rglob("*")):
        if path.is_file():
            print(f"{path.relative_to(out)}\t{path.stat().st_size}")


if __name__ == "__main__":
    main()
