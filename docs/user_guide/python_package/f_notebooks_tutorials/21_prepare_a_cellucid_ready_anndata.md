# Prepare a Cellucid-ready AnnData

**Audience:** computational users who share, host, or repeatedly open a dataset
**Time:** 20–30 minutes
**Goal:** write one standard AnnData the web app reads fast, and open it

Runnable version: {doc}`prepare_pancreas` (the checked-in Pancreas dataset,
2,531 cells).

## At a glance

- `prepare()` takes an AnnData (`data=adata`), or arrays and data frames, or
  both, and writes one standard AnnData (`.h5ad`, `.zarr` or `.zarr.zip`) —
  any AnnData tool reads it. A MuData is written as a Cellucid-ready MuData
  ({doc}`24_citeseq_mudata_prepare_and_show`).
- The layout is chosen for the browser: the expression matrix gene-major
  (CSC), arrays chunked for HTTP range reads, codecs browsers decode natively.
- `uns['cellucid']` names the dataset and carries its citation, source and
  launch defaults, so the file opens without `dataset_name`/`dataset_id`.
- Values are stored as given (float32). The one thing computed is the outlier
  filter's quantiles, when `latent_key` names a latent space.

## Why prepare (one paragraph)

The viewer reads any AnnData itself, in place. What `prepare()` changes is the
layout: in a row-major (CSR) matrix the browser must first scan every stored
value to find one gene, while in the gene-major matrix `prepare()` writes a gene
is one range. On a 561,947-cell dataset over a simulated 80 ms / 50 MB/s link
the first gene took 4,110 ms from the anndata-default file and 86 ms from the
prepared one ({doc}`../c_data_preparation_api/10_performance_and_storage`).

## Minimal call (copy/paste)

```python
import anndata as ad
from cellucid import prepare

adata = ad.read_h5ad("my_data.h5ad")

path = prepare(
    "my_dataset.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    dataset_id="my-dataset-v1",
    dataset_name="My dataset",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
)
```

## Step 1 — Choose the output and the identity

- **Container (`format`, required).** `"zarr"`: a folder of small files, for
  GitHub (no file over 100 MB, no LFS). `"h5ad"`: one file, for everything
  else — single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server),
  local files, `cellucid serve` and Jupyter. `"zarr.zip"`: the Zarr store in
  one file, wherever an `.h5ad` goes ({ref}`which-format`). The path ends in
  the matching `.h5ad`, `.zarr` or `.zarr.zip`.
  Writing to an existing path needs `force=True`.
- **`dataset_id`.** The stable key sessions and community annotations use:
  1–180 ASCII letters, digits, `.`, `_` or `-`. Version it
  (`pbmc-10k-v2`) whenever the cells change
  ({doc}`../b_concepts_mental_models/04_dataset_identity_and_reproducibility`).
- **`dataset_name`**, and optionally `dataset_description`, `citation_text` /
  `citation_url`, `source_name` / `source_url`.

## Step 2 — Gather the inputs

### Embeddings (required)

```python
embeddings = {
    "X_umap": adata.obsm["X_umap"],
    "X_pca": adata.obsm["X_pca"],  # also the latent space, named below
}
```

Every array is stored whole, in float32. Pass a latent space (scVI, PCA) for
the viewer's outlier filter as one of them and name it with `latent_key`:
`prepare()` then stores every categorical `obs` column's outlier quantiles in
the file, so the filter works at once
({doc}`../../web_app/e_filtering/02_outlier_filtering_per_active_field`).

### Expression and genes (for gene colouring)

`gene_expression` is cells × genes, sparse or dense; `var` has one row per
gene. Name the symbol column with `var_gene_id_column` when `var.index` holds
ids, and keep a subset with `gene_identifiers=[...]`
({doc}`../c_data_preparation_api/05_var_gene_metadata`).

### Optional

- `connectivities=adata.obsp["connectivities"]` (symmetric; a diagonal is left out)
- `vector_fields={"velocity": {"X_umap": adata.obsm["velocity_umap"]}}` (stored as `obsm['velocity_umap']`)
  ({doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`)
- launch defaults: `embedding_defaults`, `obs_keys`, `vector_field_default`,
  `normalize_embeddings`, `centroid_min_points` (the fewest cells a category
  needs to get outlier quantiles; 10 when None)

## Step 3 — Run `prepare(...)`

```python
from cellucid import prepare

path = prepare(
    "pbmc_10k_v1.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"], "X_pca": adata.obsm["X_pca"]},
    dataset_id="pbmc-10k-v1",
    dataset_name="PBMC 10k",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp["connectivities"],
    latent_key="X_pca",
    var_gene_id_column="gene_symbols",
    obs_keys=["cell_type", "leiden", "total_counts"],
    citation_text="Doe et al. 2026",
    citation_url="https://doi.org/10.0000/example",
)
```

Every input is checked first; a refusal names the argument and leaves nothing
on disk ({doc}`../c_data_preparation_api/11_troubleshooting_prepare`).

## Step 4 — Check and open it

```python
import anndata as ad

prepared = ad.read_h5ad(path)
print(prepared.uns["cellucid"])
```

In a notebook:

```python
from cellucid import show_anndata

viewer = show_anndata(path, height=650)
```

From a terminal: `cellucid serve pbmc_10k_v1.h5ad`. Or open the file with the
web app's **Open…** button, or host it and paste its address into the
**Link or server** bar.

## Interface reference

```{figure} ../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: The whole Cellucid window with the Pancreas sample coloured by cell_type: the sidebar shows Coloring and filtering with the categorical field selected and its eight-row legend, and the 3D embedding fills the canvas in the matching colours.
:width: 1440px

A loaded dataset in Cellucid: the sidebar controls the active view while the categorical legend maps directly to the colored points.
```

## Edge cases

### Shape mismatches

Every per-cell input needs one row per cell in the same order; the first
mismatch is named with both counts. Take every array from one AnnData.

### Duplicate names

Repeated cell or gene names are made unique (`MATR3`, `MATR3-1`) and the change
is printed.

### Very large counts

Integer counts above 2^24 keep the matrix in float64 so they read back exactly;
everything else is float32.

## Troubleshooting

### Symptom: “The file is huge”

Keep fewer genes (`gene_identifiers=[...]`), pass a sparse matrix, or subset
cells before the call.

### Symptom: “`prepare(...)` is slow”

Compressing the chunks takes most of the time, and it runs on several cores
for every container; a `.zarr.zip` adds packing the folder. Converting a large CSR
matrix to CSC also takes a few seconds.

### Symptom: “The viewer opens it but gene colouring is missing”

`gene_expression` and `var` were not passed, or the gene names are ids: set
`var_gene_id_column`.

## Next steps

- The contract the file follows: {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`
- Large datasets on a server: {doc}`22_large_dataset_server_mode_and_lazy_gene_expression`
