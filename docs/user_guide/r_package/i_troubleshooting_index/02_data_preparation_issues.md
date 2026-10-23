# Data Preparation Issues

**Audience:** everyone calling `cellucid_prepare()`  
**Time:** 5 minutes  
**What you'll get:** the root causes behind most failures, in priority order

Use this page when `cellucid_prepare()` stops with an error, or writes a store
that looks wrong in the viewer.

```{important}
The full guide is {doc}`../c_data_preparation_api/11_troubleshooting_prepare`,
and every message is listed in
{doc}`../g_api_reference_coverage/02_error_messages_and_exceptions_document_patterns`.
```

## Most common root causes (in priority order)

### 1) Cells in different orders

Symptoms:
- `The row names of <a> and <b> disagree at row <i> (...). Every input must list the same cells in the same order.`
- or, when the inputs carry no row names, a store whose clusters or genes sit
  on the wrong cells

Fix:
- choose one cell order and reorder *every* input to it; keep row names on
  `obs`, `gene_expression` and the embeddings so the check can catch a mismatch
- see {doc}`../c_data_preparation_api/02_input_requirements_global`

### 2) Expression in the wrong orientation

Symptoms:
- `gene_expression has <n> rows, but the embeddings have <m> cells. gene_expression is cells x genes; transpose a genes x cells matrix with t().`
- `var has <n> rows, but gene_expression has <m> genes (columns).`

Fix:
- `gene_expression` is **cells x genes**; Seurat and SingleCellExperiment store
  genes x cells, so pass `Matrix::t(...)`
- see {doc}`../c_data_preparation_api/06_gene_expression_matrix`

### 3) Names that cannot become paths

Symptoms:
- `obs column '<name>' contains a character that cannot appear in a file name on every platform ...`
- `... differ only in letter case. ...`
- `obs column '_index' is reserved ...`

Fix:
- every `obs` / `var` column name and every `obsm` key becomes a folder in the
  store: rename the offending ones (for example
  `names(obs) <- make.names(names(obs), unique = TRUE)`)
- see {doc}`../c_data_preparation_api/04_obs_cell_metadata`

### 4) Columns with no AnnData encoding

Symptoms:
- `obs column '<name>' is a Date column, which has no AnnData encoding Cellucid writes. ...`
- `obs column '<name>' holds 1 infinite value ...`

Fix:
- convert dates, times, lists and other classed columns with `as.character()`
  or `as.numeric()`; replace `Inf` with `NA`
- a character column whose values repeat becomes a categorical; convert a
  numeric code you meant as a category with `factor()`

### 5) A vector field bound to an embedding not named `X_<basis>`

Symptom:
- `Vector field '<field>' on embedding '<key>': the viewer binds a vector field to the embedding 'X_<basis>' it is named for ('velocity_umap' to 'X_umap'), so the embedding needs an 'X_<basis>' key. ...`

Fix:
- name the embedding `X_<basis>` (`X_umap`, `X_pca`, `X_branches`), or store
  the vectors as an embedding and bind them in the viewer
- see {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`

## Quick “what do I do next?” links

- Seurat recipe: {doc}`../e_integrations_recipes/01_seurat_recipe`
- SingleCellExperiment recipe: {doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`
- The store looks wrong after writing: {doc}`03_store_format_and_validation_issues`
