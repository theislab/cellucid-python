# Public Functions

**Audience:** everyone who calls the package  
**Time:** 5 minutes  
**What you'll get:** the complete public surface of `cellucid-r` on one page

`cellucid-r` exports two functions. It has no classes, no S4 methods, and no
other exported helpers.

| Function | What it does | Returns |
|---|---|---|
| `cellucid_prepare()` | Writes a Cellucid-ready AnnData store: a `.zarr` folder (`format = "zarr"`) or one `.zarr.zip` file (`format = "zarr.zip"`) | Invisibly, the absolute path of the written store |
| `cellucid_transition_drift()` | Computes the drift of every cell in an embedding's coordinates from a cell-cell transition matrix, a vector field to pass to `cellucid_prepare(vector_fields = )` | A double matrix as wide as the embedding |

## Signature

```r
cellucid_prepare(
  out_path,
  format,
  embeddings,
  dataset_id,
  dataset_name,
  obs = NULL,
  var = NULL,
  gene_expression = NULL,
  connectivities = NULL,
  vector_fields = NULL,
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
  zarr_format = NULL,
  force = FALSE
)
```

`out_path`, `format`, `embeddings`, `dataset_id`, and `dataset_name` have no
default.
Leaving any of them out is reported, all at once, as
`cellucid_prepare() is missing N required argument(s):` followed by each name
and what a valid value is.

## Arguments by purpose

### Where the store goes

| Argument | Default | Meaning | Explained in |
|---|---|---|---|
| `out_path` | required | A folder ending in `.zarr` for `format = "zarr"`, a file ending in `.zarr.zip` for `format = "zarr.zip"`; a mismatch is refused with the path to write instead; missing parent folders are created | {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification` |
| `format` | required | `"zarr"` or `"zarr.zip"`. No default: a folder of many small files and one file suit different hosts, and only you know where the store goes. `"zarr"` writes a `.zarr` folder, for GitHub (no file over 100 MB, no LFS); `"zarr.zip"` writes the same store as one `.zarr.zip` file, for everything else: single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files, `cellucid serve` and Jupyter. There is no `"h5ad"`: this package writes no HDF5, so it needs no HDF5 library; Python's `cellucid.prepare(format = "h5ad")` writes `.h5ad`. | {doc}`../c_data_preparation_api/index` |
| `zarr_format` | `NULL` (2) | `2` or `3`: the Zarr format of the store; anndata 0.8 and later read format 2, anndata 0.11 and later format 3, the web app both equally fast | {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification` |
| `force` | `FALSE` | `TRUE` replaces an existing store once the new one is complete; never a folder that is not a Zarr store | {doc}`../c_data_preparation_api/01_cellucid_prepare_overview` |

### Data

| Argument | Default | Stored as | Explained in |
|---|---|---|---|
| `embeddings` | required | `obsm[<key>]`, float32, every column | {doc}`../c_data_preparation_api/03_embeddings_and_coordinates` |
| `obs` | `NULL` | `obs`, every column | {doc}`../c_data_preparation_api/04_obs_cell_metadata` |
| `var` | `NULL` | `var`, every column (requires `gene_expression`) | {doc}`../c_data_preparation_api/05_var_gene_metadata` |
| `gene_expression` | `NULL` | `X`, cells x genes, gene-major float32, float64 for integer counts above 2^24 (`csc_matrix` for a sparse matrix, stored zeros dropped) | {doc}`../c_data_preparation_api/06_gene_expression_matrix` |
| `connectivities` | `NULL` | `obsp['connectivities']`, float32 `csr_matrix` | {doc}`../c_data_preparation_api/07_connectivities_knn_graph` |
| `vector_fields` | `NULL` | `obsm['<field>_<basis>']` for each embedding `X_<basis>` (`velocity_umap`) | {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement` |
| `gene_identifiers` | `NULL` | Narrows `X` and `var` to these genes, in this order | {doc}`../c_data_preparation_api/05_var_gene_metadata` |

### What the dataset says about itself (`uns['cellucid']`)

| Argument | Default | Stored as |
|---|---|---|
| `dataset_id` | required | `uns['cellucid']['id']` |
| `dataset_name` | required | `uns['cellucid']['name']` and `uns['title']` |
| `dataset_description` | `NULL` | `uns['cellucid']['description']` |
| `citation_text`, `citation_url` | `NULL` | `uns['cellucid']['citation']` `text`, `url` (a URL needs its text) |
| `source_name`, `source_url` | `NULL` | `uns['cellucid']['source']` `name`, `url` (a URL needs its name) |

See {doc}`../b_concepts_mental_models/03_dataset_identity_and_reproducibility`.

### Launch defaults (`uns['cellucid']['defaults']`, only the ones you set)

| Argument | Stored key | Meaning |
|---|---|---|
| `embedding_defaults` | `embedding1d`, `embedding2d`, `embedding3d` | The embedding each view dimension shows first: an `obsm` key, or `list(key = , columns = )` with 1-based columns (stored 0-based) |
| `latent_key` | `latentKey` | The latent space (an `embeddings` key, such as `"X_pca"`) the viewer's outlier filter measures in; also stores `outlierQuantiles` (below) |
| `obs_keys` | `obsKeys` | The `obs` columns the viewer shows, in order (every column is stored either way) |
| `var_gene_id_column` | `geneIdColumn` | The `var` column that names genes in the viewer |
| `vector_field_default` | `vectorFieldDefault` | The vector field selected first, by its `vector_fields` name (`"velocity"`), which is also the viewer's field id |
| `normalize_embeddings` | `normalizeEmbeddings` | `FALSE` shows embeddings at their own coordinates |
| `centroid_min_points` | `centroidMinPoints` | The fewest cells a category needs to get outlier quantiles (10 when `NULL`) |

An argument left at `NULL` writes no key, so the viewer keeps its own choice.

### Outlier quantiles (`uns['cellucid']['outlierQuantiles']`, only with `latent_key`)

| Key | Type | Value |
|---|---|---|
| `latentKey` | text | `latent_key` |
| `minPoints` | int64 | `centroid_min_points`, or 10 |
| `fields` | dict | one float32 array of one value per cell for every `obs` column the viewer shows as categorical (a factor, a character column, a logical column with or without `NA`), in `obs` column order; only a column with a category of at least `minPoints` cells |

Each value is the share of the cell's category at its distance from the
category's centre in the latent space, or closer; `NaN` for `NA` and for a
category under `minPoints` cells. `?cellucid_prepare` states the rule, which
Python's `cellucid.prepare()` and the web app follow bit for bit.

## `cellucid_transition_drift()`

```r
cellucid_transition_drift(transition_matrix, embedding, normalize_rows)
```

All three arguments are required.

| Argument | Meaning |
|---|---|
| `transition_matrix` | Cells x cells, finite and non-negative: a `Matrix` sparse or dense matrix with numeric values, or a base numeric matrix; row `i` holds the weights of the cells cell `i` moves to (CellRank's transition matrix, for example) |
| `embedding` | The coordinates the drift is computed in, one row per cell, in any form `cellucid_prepare()` takes for an embedding; every column is used, every value finite and representable in float32 |
| `normalize_rows` | `TRUE` divides every row of `transition_matrix` by its sum (a zero row sum is refused); `FALSE` uses it as given |

It returns `T %*% E - E` (the product divided by the row sums first with
`normalize_rows = TRUE`) as a double matrix with the embedding's dimensions
and dimnames, the computation of Python's `cellucid.compute_transition_drift()`.
Store it with `vector_fields = list(T_fwd = list(X_umap = drift))`; see
{doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`.

## Messages

`cellucid_prepare()` reports with `message()` (not a warning) when it makes
repeated cell or gene names unique, for example
`Gene names: 1 duplicate name was made unique, e.g. 'A' → 'A-1'.` Every other
problem stops the call with an error; see
{doc}`02_error_messages_and_exceptions_document_patterns`.
