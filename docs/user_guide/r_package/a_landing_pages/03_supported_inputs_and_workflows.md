# Supported Inputs and Workflows

**Audience:** computational users and “power” wet lab users  
**Time:** 10 minutes  
**Goal:** know exactly what `cellucid_prepare()` expects and how to avoid the most common gotchas.

## The supported workflow (high level)

1) Look at a Seurat, SingleCellExperiment or SummarizedExperiment object with
   `cellucid_view(object)` (the Viewer pane or a browser; nothing written).
2) Write it with `cellucid_prepare(out_path, format, data = object, ...)`, or
   write matrices and data frames you extracted yourself.
3) Open the `.zarr` folder or `.zarr.zip` file in the Cellucid web app, or
   share it ({doc}`../d_viewing_loading/index`).

## What inputs are supported?

### Required inputs

Without `data`, these five carry no default (with `data`, only `out_path` and
`format`); leaving one out stops the call with
`cellucid_prepare() is missing 1 required argument:` and what a valid value is.

- **`out_path`**: a folder path ending in `.zarr` for `format = "zarr"`, or a file path ending in `.zarr.zip` for `format = "zarr.zip"`
- **`format`**: `"zarr"` (a folder of small files, for GitHub) or `"zarr.zip"` (one file, for everything else: single-file hosts such as Zenodo, Figshare, Hugging Face, S3 or a lab server, local files, `cellucid serve` and Jupyter); there is no default, because only you know where the store goes ({doc}`../c_data_preparation_api/index`)
- **`embeddings`**: a named list of matrices with one row per cell, any number
  of them, any width (`list(X_umap = umap, X_pca = pca)`); each name becomes an `obsm` key
- **`dataset_id`**: a portable identifier (ASCII letters, digits, `.`, `_`, `-`)
- **`dataset_name`**: the name shown in the viewer

### Optional inputs

- **`data`**: a Seurat, SingleCellExperiment or SummarizedExperiment object;
  every input below then defaults to what it holds, and one you pass replaces
  that component
- **`modalities`**: which assays (Seurat) or experiments (SingleCellExperiment)
  of `data` to write, in order; with several, the store is a Cellucid-ready
  MuData ({doc}`../d_viewing_loading/07_multimodal_assays`)
- **`obs`**: `data.frame` of cell metadata, one row per cell
- **`gene_expression`**: cells × genes, dense `matrix` or sparse `Matrix` (recommended)
- **`var`**: `data.frame`, one row per gene; only together with `gene_expression`
- **`gene_identifiers`**: the genes to store (default: all)
- **`connectivities`**: cells × cells neighbour graph
- **`vector_fields`**: per-cell vectors drawn over a UMAP, `list(velocity = list(X_umap = v2))`
- **Metadata**: `dataset_description`, `citation_text`, `citation_url`, `source_name`, `source_url`
- **Launch defaults**: `embedding_defaults`, `latent_key`, `obs_keys`,
  `var_gene_id_column`, `vector_field_default`, `normalize_embeddings`,
  `centroid_min_points`, `modality` (the modality a multimodal store opens
  on); with `latent_key`, the outlier quantiles of every categorical `obs`
  column are stored too
- **Output**: `zarr_format` (2, the default, or 3), `force` (replace an existing store)

## Accepted types (what “matrix-like” means)

| Input | Accepted |
|---|---|
| embeddings, vector fields | base numeric `matrix`, a `data.frame` whose every column is a plain numeric vector, a `Matrix` matrix, or a plain numeric vector (one column) |
| `gene_expression` | a `Matrix` sparse matrix (stored as CSC, never densified), or a base numeric `matrix`, numeric `data.frame`, or dense `Matrix` |
| `connectivities` | a `Matrix` matrix (recommended), or a base numeric or logical `matrix` (needs `Matrix` installed) |
| `obs`, `var` columns | `factor`, `character`, `logical`, `integer`, `numeric` |

```{warning}
A column or matrix carrying a class of its own — `bit64::integer64`, a `units`
measurement, a `Date` — is **refused**, not silently converted, because the
conversion would drop the attribute that says what the values mean. Convert it
deliberately (for example with `as.character()` or `as.numeric()`) before the call.
```

## Shapes and orientation (the #1 source of bugs)

Rows are **cells**; columns are dimensions, genes, or components.

| Input | Expected shape |
|---|---|
| an embedding | `(n_cells, k)`, any `k ≥ 1` |
| `obs` | `n_cells` rows |
| `gene_expression` | `(n_cells, n_genes)` |
| `var` | `n_genes` rows |
| `connectivities` | `(n_cells, n_cells)` |
| a vector field | `(n_cells, k)`, `k` = the column count of the embedding it is drawn over |
| `cellucid_transition_drift()`'s `transition_matrix` | `(n_cells, n_cells)` |

```{warning}
Many R containers store expression as **genes × cells**. With `data = ` the
object's own orientation is read and nothing needs transposing; when you pass
`gene_expression` yourself, you almost always need to transpose:

- Seurat: `GetAssayData(seu, layer = "data")` is genes × cells → `Matrix::t(...)`
- SingleCellExperiment: `assay(sce, "logcounts")` is genes × cells → `Matrix::t(...)`

`cellucid_prepare()` checks the dimensions but never guesses the orientation;
a square matrix can be the wrong way round and still match both axes.
```

## Cell identity = row order (and row names)

Row `i` of every input is the same cell. The store's cell names (`obs_names`)
come from the row names of `obs`, `gene_expression`, and the embeddings; where
several of them carry row names, they must be **the same names in the same
order**, or the call stops and names the first row where they disagree. With
no row names at all, cells are named `"0"` to `"n - 1"`, as anndata names them.

## Names that become paths

Every `obs`/`var` column name and every `obsm` key becomes a folder in the
store, so it must be a portable file name: no `/ \ : * ? " < > |` or control
character, no trailing space or dot, not a Windows device name, not `_index`
for a column, and distinct from its neighbours by more than letter case. See
{doc}`../c_data_preparation_api/02_input_requirements_global`.

## Recommended workflow decision tree

- **Seurat object** → {doc}`../e_integrations_recipes/01_seurat_recipe`
- **SingleCellExperiment** → {doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`
- **Raw matrices** → {doc}`../e_integrations_recipes/03_raw_matrices_and_data_frames_recipe`

Then:
- Write the store → {doc}`../c_data_preparation_api/index`
- Open it in the browser → {doc}`../d_viewing_loading/index`

## Next steps

- Copy/paste example: {doc}`04_quick_start_3_levels`
- Full API documentation: {doc}`../c_data_preparation_api/index`
