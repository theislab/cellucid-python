# `cellucid_prepare()` overview

**Audience:** everyone using the R package
**Time:** 10 minutes
**Goal:** know what the function writes, which arguments it takes, and where each one ends up.

`cellucid_prepare()` takes a Seurat, SingleCellExperiment or
SummarizedExperiment object (`data =`), or plain R objects -- matrices, sparse
`Matrix` matrices, and data frames -- and writes them as a **Cellucid-ready
AnnData** (an object with several assays as a **Cellucid-ready MuData**):

- a standard AnnData store in **AnnData Zarr** (Zarr format 2, or 3 with
  `zarr_format = 3`), which Python's `anndata` reads without a warning;
- **gene-major** expression (`X` as a float32 `csc_matrix` from a sparse
  matrix), chunked so the web app reads one gene with one small range of the
  file;
- `uns['cellucid']`: the dataset's id, name, description, citation, source,
  and the launch defaults you chose;
- consolidated metadata (`.zmetadata`), so a store on a web host opens with
  one request.

The output is a **`.zarr` folder** (`format = "zarr"`), or **one `.zarr.zip`
file** holding the same files. Open either in the web app
({doc}`../d_viewing_loading/01_open_stores_in_cellucid_web_app`), or read it in
Python with `anndata.read_zarr()`.

---

## From an object

```r
cellucid_prepare("pbmc.zarr", "zarr", data = seu, dataset_id = "pbmc", dataset_name = "PBMC")
```

Every other argument defaults to what the object holds, and an argument you
give replaces that one component (`obs = seu[[c("seurat_clusters")]]` stores
only that column). What each part of the object becomes:
{doc}`../e_integrations_recipes/01_seurat_recipe`,
{doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`; several assays:
{doc}`../d_viewing_loading/07_multimodal_assays`. To look before writing,
`cellucid_view(seu)` ({doc}`../d_viewing_loading/04_view_in_rstudio_positron`).

## Minimal call

```r
library(cellucid)

cellucid_prepare(
  out_path = "pbmc.zarr",
  format = "zarr",
  embeddings = list(X_umap = umap),   # one row per cell
  dataset_id = "pbmc",
  dataset_name = "PBMC"
)
```

Without `data`, five arguments are required: `out_path`, `format`,
`embeddings`, `dataset_id`, and `dataset_name`; with `data`, only `out_path`
and `format` (the object supplies the rest, and its own record or the
arguments its identity). Leaving one out is reported as a missing argument,
with what a valid value is:

```text
cellucid_prepare() is missing 1 required argument:
  - dataset_name: one non-empty name for the dataset, shown to the reader verbatim.
```

## A typical call

```r
cellucid_prepare(
  out_path = "pbmc.zarr",
  format = "zarr",
  embeddings = list(X_umap = umap, X_pca = pca),
  dataset_id = "pbmc-3k",
  dataset_name = "PBMC 3k",
  obs = obs,                         # data.frame, one row per cell
  var = var,                         # data.frame, one row per gene
  gene_expression = expr,            # cells x genes, sparse or dense
  connectivities = knn,              # cells x cells, symmetric
  embedding_defaults = c("2d" = "X_umap"),
  latent_key = "X_pca",
  var_gene_id_column = "symbol"
)
```

## Signature

```r
cellucid_prepare(
  out_path,
  format,
  data = NULL,
  embeddings = NULL,
  dataset_id = NULL,
  dataset_name = NULL,
  obs = NULL,
  var = NULL,
  gene_expression = NULL,
  connectivities = NULL,
  vector_fields = NULL,
  modalities = NULL,
  gene_identifiers = NULL,
  dataset_description = NULL,
  citation_text = NULL,
  citation_url = NULL,
  source_name = NULL,
  source_url = NULL,
  embedding_defaults = NULL,
  latent_key = NULL,
  obs_keys = NULL,
  var_gene_id_column = NULL,
  vector_field_default = NULL,
  normalize_embeddings = NULL,
  centroid_min_points = NULL,
  modality = NULL,
  zarr_format = NULL,
  force = FALSE
)
```

It returns, invisibly, the absolute path of the written store.

## Arguments and where they go

### The data

| Argument | Required | What it is | In the store |
|---|---|---|---|
| `out_path` | yes | A folder path ending in `.zarr` for `format = "zarr"`, a file path ending in `.zarr.zip` for `format = "zarr.zip"`. Missing parent folders are created. | the store itself |
| `format` | yes | `"zarr"` (a `.zarr` folder) or `"zarr.zip"` (one `.zarr.zip` file). No default: a folder of many small files and one file suit different places, and only you know where the store goes. `"zarr"` writes a `.zarr` folder, for GitHub (no file over 100 MB, no LFS); `"zarr.zip"` writes the same store as one `.zarr.zip` file, for everything else: single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files, `cellucid serve` and Jupyter ({doc}`index`). There is no `"h5ad"`: this package writes no HDF5, so it needs no HDF5 library; Python's `cellucid.prepare(format = "h5ad")` writes `.h5ad`. | the store itself |
| `data` | no | A Seurat, SingleCellExperiment or SummarizedExperiment object; every argument below defaults to what it holds, and one you pass replaces that component. | every element; several assays or experiments as a Cellucid-ready MuData |
| `embeddings` | without `data` | Named list of numeric matrices, one row per cell. | `obsm[<name>]`, float32, every column |
| `obs` | no | `data.frame` of cell metadata, one row per cell. | `obs`, every column |
| `var` | no | `data.frame` of gene metadata, one row per gene; needs `gene_expression`. | `var`, every column |
| `gene_expression` | no | Cells x genes matrix, sparse (`Matrix`) or dense. | `X`, float32 (float64 for integer counts above 2^24), gene-major |
| `connectivities` | no | Cells x cells neighbour graph, exactly symmetric, zero diagonal. | `obsp['connectivities']`, float32 CSR |
| `vector_fields` | no | Per-cell vectors bound to an embedding `X_<basis>`, as many columns as it: `list(velocity = list(X_umap = v2))`. | `obsm['<field>_<basis>']` (`velocity_umap`) |
| `modalities` | no | The modalities of `data` to write, in this order (default: all, the Seurat active assay or the main experiment first); one name writes a plain AnnData. | `mod/<name>`, `mod-order` |
| `gene_identifiers` | no | The genes to store, in this order (default: all). | narrows `X` and `var` |

### The dataset's description (`uns['cellucid']`)

| Argument | Required | In the store |
|---|---|---|
| `dataset_id` | without `data`, or an object that names none | `uns['cellucid']['id']` (1-180 portable ASCII characters) |
| `dataset_name` | without `data`, or an object that names none | `uns['cellucid']['name']` and `uns['title']` |
| `dataset_description` | no | `uns['cellucid']['description']` |
| `citation_text`, `citation_url` | no (the URL needs the text) | `uns['cellucid']['citation']` = `{text, url}` |
| `source_name`, `source_url` | no (the URL needs the name) | `uns['cellucid']['source']` = `{name, url}` |

### Launch defaults (`uns['cellucid']['defaults']`)

Each one is stored only when you set it; anything you leave `NULL` stays the
viewer's own choice.

| Argument | Stored as | Meaning |
|---|---|---|
| `embedding_defaults` | `embedding1d`, `embedding2d`, `embedding3d` | The embedding each view dimension shows first: an obsm key, or `list(key = , columns = )` with 1-based columns (stored 0-based). |
| `latent_key` | `latentKey` | The latent space the viewer's outlier filter measures in, such as `"X_pca"`; never a 2D/3D layout. The outlier quantiles of every categorical obs column are computed in it and stored in `uns['cellucid']['outlierQuantiles']`, so the filter works as soon as the file opens. |
| `obs_keys` | `obsKeys` | The obs columns the viewer shows, in this order. Every column is stored either way. |
| `var_gene_id_column` | `geneIdColumn` | The var column whose values name genes in the viewer. |
| `vector_field_default` | `vectorFieldDefault` | The field selected first, by its `vector_fields` name (`"velocity"`), which is also the viewer's field id. |
| `normalize_embeddings` | `normalizeEmbeddings` | `FALSE` shows embeddings at their own coordinates. |
| `centroid_min_points` | `centroidMinPoints` | Fewest cells a category needs to get outlier quantiles (10 when `NULL`). |
| `modality` | `modality` | The modality a multimodal store opens on; without it the first named like gene expression (`rna`, `gex`, `gene expression`, `gene_expression`, `sct`, any case), else the first. Refused for a store of one feature space. |

### Writing

| Argument | Default | Meaning |
|---|---|---|
| `zarr_format` | `NULL` (2) | `2` (read by every `anndata` since 0.8 and zarr-python 2 and 3) or `3` (`anndata` 0.11 and later). The web app reads both equally fast. |
| `force` | `FALSE` | `TRUE` replaces an existing store at `out_path`. A folder that is not a Zarr store is never replaced. |

## How a call runs

1. Every argument is checked before anything is written: shapes, types,
   names, and the rules on each page of this section. A mistake stops the call
   with a message naming the argument ({doc}`11_troubleshooting_prepare`).
2. The store is written **beside** `out_path` under a temporary hidden name
   (`.<name>.<process id>.partial`).
3. Only when the store is complete is it moved to `out_path`. With
   `force = TRUE` the old store is moved aside first and deleted after the new
   one is in place.

A failed or interrupted call therefore leaves nothing at `out_path`, and never
a half-written store; the store it would have replaced stays as it was.

## Messages

`cellucid_prepare()` prints a message (not an error) for each part of an
object it leaves out (`From the Seurat object: ...`: a column of a type
AnnData does not store, a reduction that does not cover every cell, a graph
that is not symmetric), and when it makes repeated names unique, for example:

```text
Gene names: 1 duplicate name was made unique, e.g. 'MATR3' → 'MATR3-1'.
```

See {doc}`05_var_gene_metadata` and {doc}`02_input_requirements_global`.

## Next

- {doc}`02_input_requirements_global`: the rules every input shares
- {doc}`09_cellucid_ready_anndata_specification`: the exact layout of the store
- {doc}`../e_integrations_recipes/01_seurat_recipe` and
  {doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`: from Seurat
  and SingleCellExperiment objects
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/01_prepare_overview`
