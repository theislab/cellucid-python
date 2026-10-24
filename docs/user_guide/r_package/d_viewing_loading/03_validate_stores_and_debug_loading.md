# Validate a Store & Debug Loading

**Audience:** computational users, maintainers
**Time:** 10–20 minutes
**Goal:** confirm a store is what you meant to write, and find out why a store
does not open.

A store written by `cellucid_prepare()` is a standard AnnData Zarr store, so the
best validation is the reference reader: Python's `anndata`. If anndata reads
the store without a warning and the values are yours, the web app reads it too.

## 1. Look at the store from R

`cellucid_prepare()` returns the store's path. Its metadata is plain JSON:

```r
library(cellucid)

set.seed(1)
umap <- matrix(rnorm(200), ncol = 2)
obs <- data.frame(cluster = factor(sample(c("A", "B"), 100, TRUE)))

store <- cellucid_prepare(
  out_path = file.path(tempdir(), "check.zarr"),
  format = "zarr",
  embeddings = list(X_umap = umap),
  dataset_id = "check",
  dataset_name = "Check",
  obs = obs,
  force = TRUE
)

list.files(store, all.files = TRUE, no.. = TRUE)
jsonlite::read_json(file.path(store, ".zattrs"))
jsonlite::read_json(file.path(store, "obs", ".zattrs"))
jsonlite::read_json(file.path(store, "obsm", "X_umap", ".zarray"))$shape
```

What to look for:

- the root `.zattrs` says `encoding-type: anndata`;
- `obs/.zattrs` lists every column you passed in `column-order`;
- each `obsm` array has one row per cell;
- `.zmetadata` exists (it always does: it is written last, and the store is
  moved into place only after it).

A store is written beside `out_path` under a temporary name and moved into place
only when complete, so a store at `out_path` is never half-written.

## 2. Read it back with Python's anndata

```python
import warnings
warnings.simplefilter("error")   # a warning is a failure

import anndata
adata = anndata.read_zarr("check.zarr")
print(adata)
print(adata.uns["cellucid"])
```

For a `.zarr.zip`:

```python
import anndata, zarr
adata = anndata.read_zarr(zarr.storage.ZipStore("check.zarr.zip", mode="r"))
```

Compare what you see with your R objects: `adata.obs_names` with the cell names,
`adata.var_names` with the gene names, `adata.obsm["X_umap"]` with the
embedding (as float32), `adata.X` with the expression matrix.

## 3. Open it in the web app

Open the `.zarr.zip` with **Open… › File…**, or the hosted store with the
**Link or server** bar ({doc}`01_open_stores_in_cellucid_web_app`). The browser console (open
the developer tools) logs what the reader opened, for example
`[AnnDataReader] Opened check.zarr: 100 cells, 0 genes`.

## Common failures

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| **Open… › File…** refuses the file | The file is not a ZIP, or it is a `.zarr` folder you zipped yourself without its hidden files | Write it with `format = "zarr.zip"` |
| A link reports nothing published at that address | Wrong path, branch, or a folder address without the trailing `/` store | Open the address in a browser tab; for a folder, check that `<address>/.zmetadata` downloads |
| A link reports CORS or ignored range requests | The host does not let other sites read its files | Use GitHub, Hugging Face or Zenodo, or fix the host's headers ({doc}`02_host_stores_for_sharing`) |
| The app refuses `uns['cellucid']` | A value in it was edited by hand into something the viewer does not accept | Write the store again with `cellucid_prepare()` |
| No vector field overlay | The field is not bound to the UMAP the view shows | {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement` |
| Fields or genes are missing | Wrong `obs_keys`, or `gene_identifiers` narrowed the genes | {doc}`../c_data_preparation_api/04_obs_cell_metadata`, {doc}`../c_data_preparation_api/05_var_gene_metadata` |

More: {doc}`../i_troubleshooting_index/03_store_format_and_validation_issues`
and {doc}`../i_troubleshooting_index/04_web_app_loading_issues`.
