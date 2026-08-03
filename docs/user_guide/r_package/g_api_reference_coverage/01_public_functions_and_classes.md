# Public Functions and Objects

**Audience:** everyone looking an argument up  
**Time:** 5 minutes  
**What you'll get:** every argument of `cellucid_prepare()` with its default and the page that explains it

`cellucid-r` exports exactly one user-facing function. There are no classes, no
S4 generics, and no companion helpers: `cellucid_prepare()` writes an export
folder and returns `NULL` invisibly.

## The shortest call that works

Three arguments carry no default, and you also need one embedding, a latent
space, and cell metadata:

```r
cellucid_prepare(
  latent_space = pca,
  obs = obs,
  X_umap_2d = umap2,
  obs_categorical_dtype = "uint16",
  dataset_name = "My dataset",
  dataset_id = "my_dataset",
  out_dir = file.path(getwd(), "exports", "my_dataset")
)
```

Everything else on this page is optional. Start at
{doc}`../a_landing_pages/04_quick_start_3_levels` if you want a runnable
version, or {doc}`../c_data_preparation_api/01_cellucid_prepare_overview` for
what the exporter does with each input.

## `cellucid_prepare()` arguments

A dash in the **Default** column means the argument has no default at all:
leaving it out is reported as a missing argument, not as a value of the wrong
type.

### Required

| Argument | Default | What it is |
|---|---|---|
| `latent_space` | `NULL`, but required | Numeric matrix-like `(n_cells, n_dims)`. Measured, never written: the per-cell outlier quantiles are distances in it. |
| `obs` | `NULL`, but required | `data.frame` with `n_cells` rows. |
| `obs_categorical_dtype` | — | Exactly `"uint8"` or `"uint16"`. The storage width of every categorical code. |
| `dataset_name` | — | Non-empty name, shown to the reader verbatim. |
| `dataset_id` | — | Portable identifier: 1–180 ASCII letters, numbers, `.`, `_`, `-`. |

At least one of `X_umap_1d`, `X_umap_2d`, `X_umap_3d` is also required.

Details: {doc}`../c_data_preparation_api/02_input_requirements_global` ·
{doc}`../c_data_preparation_api/04_obs_cell_metadata`

### Embeddings and overlays

| Argument | Default | What it is |
|---|---|---|
| `X_umap_1d` | `NULL` | `(n_cells, 1)`, or a plain length-`n_cells` numeric vector. |
| `X_umap_2d` | `NULL` | `(n_cells, 2)`. |
| `X_umap_3d` | `NULL` | `(n_cells, 3)`. |
| `vector_fields` | `NULL` | Named list of per-cell displacements. Every name must match `<field>_umap_<1\|2\|3>d`. |
| `vector_field_default` | `NULL` | Field id to select on load. Required when `vector_fields` holds more than one field. |
| `connectivities` | `NULL` | `(n_cells, n_cells)` symmetric graph with a zero diagonal. |

Details: {doc}`../c_data_preparation_api/03_embeddings_and_coordinates` ·
{doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement` ·
{doc}`../c_data_preparation_api/07_connectivities_knn_graph`

### Gene expression

| Argument | Default | What it is |
|---|---|---|
| `var` | `NULL` | `data.frame` with `n_genes` rows. Required when `gene_expression` is supplied. |
| `gene_expression` | `NULL` | Matrix-like `(n_cells, n_genes)`. Cells are rows; the matrix is never transposed for you. |
| `var_gene_id_column` | `NULL` | Column of `var` holding the gene identifiers. `NULL` uses `rownames(var)`. |
| `gene_identifiers` | `NULL` | Character vector selecting which genes to export. `NULL` exports all of them. |

Details: {doc}`../c_data_preparation_api/05_var_gene_metadata` ·
{doc}`../c_data_preparation_api/06_gene_expression_matrix`

### Metadata selection and category summaries

| Argument | Default | What it is |
|---|---|---|
| `obs_keys` | `NULL` | Columns of `obs` to export. `NULL` exports all of them. |
| `centroid_outlier_quantile` | `0.95` | Distance quantile kept as inliers when placing category centroids. `NULL` disables centroids; per-cell outlier quantiles are still computed. |
| `centroid_min_points` | `10` | Cells a category needs before it gets a centroid and non-`NaN` outlier quantiles. |

Details: {doc}`../c_data_preparation_api/04_obs_cell_metadata`

### Size and encoding

| Argument | Default | What it is |
|---|---|---|
| `compression` | `NULL` | gzip level `1`–`9`, or `NULL` for uncompressed payloads. |
| `var_quantization` | `NULL` | `8` or `16` bits for gene expression, or `NULL` for float32. |
| `obs_continuous_quantization` | `NULL` | `8` or `16` bits for continuous obs fields *and* categorical outlier quantiles, or `NULL` for float32. |

Details: {doc}`../c_data_preparation_api/10_performance_tuning_prepare_export`

### Output and provenance

| Argument | Default | What it is |
|---|---|---|
| `out_dir` | `file.path(getwd(), "exports")` | Where the export goes. It must name a dedicated directory: the working directory and the home directory are refused. |
| `force` | `FALSE` | `TRUE` atomically replaces an existing generation. `FALSE` requires that `out_dir` does not already exist. |
| `dataset_description` | `NULL` | Optional description. `NULL` and `""` both publish an empty description. |
| `source_name` | `NULL` | Data source name. Required whenever `source_url` or `source_citation` is given. |
| `source_url` | `NULL` | Data source URL. |
| `source_citation` | `NULL` | Citation text. |
| `created_at` | `NULL` | Exact UTC timestamp, `YYYY-MM-DDTHH:MM:SSZ`. Defaults to now; set it for a reproducible `dataset_identity.json`. |

Details: {doc}`../b_concepts_mental_models/03_dataset_identity_and_reproducibility` ·
{doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`

## Return value

`invisible(NULL)`. The function is called for its side effect — the files it
writes to `out_dir` — and there is no object to keep.

## When a call fails

Every failure is a `stop()` raised before publication, with nothing written. The
error families and how to read them are on
{doc}`02_error_messages_and_exceptions_document_patterns`; the
symptom-by-symptom guide is
{doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`.
