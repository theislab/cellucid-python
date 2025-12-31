# Supported Inputs and Workflows

**Audience:** computational users and “power” wet lab users  
**Time:** 10–15 minutes  
**Goal:** know exactly what `cellucid_prepare()` expects and how to avoid the most common gotchas.

## The supported workflow (high level)

1) Extract arrays/data.frames from your R object(s) (Seurat, SingleCellExperiment, or raw matrices).
2) Call `cellucid::cellucid_prepare(...)` (or the alias `cellucid::prepare(...)`).
3) Open the output folder in the Cellucid web app.

```{note}
The web app is the viewer UI. The R package is for export only (for now).
If you want server mode and notebook embedding, see {doc}`../../python_package/index`.
```

## What inputs are supported?

### Required inputs

`cellucid_prepare()` requires:

- **At least one embedding**: one of `X_umap_1d`, `X_umap_2d`, `X_umap_3d`
  - shape must be `(n_cells, 1/2/3)`
- **`latent_space`**:
  - numeric matrix-like, shape `(n_cells, n_dims)`
- **`obs`**:
  - `data.frame` with `n_cells` rows

Why `latent_space` is required: it is used to compute **per-cell outlier quantiles** for categorical fields (see {doc}`../c_data_preparation_api/04_obs_cell_metadata`).

### Optional inputs

- **Gene expression (`gene_expression`)**:
  - matrix-like, shape `(n_cells, n_genes)`
  - can be dense `matrix` or sparse `Matrix` (recommended)
- **Gene metadata (`var`)**:
  - `data.frame`, shape `(n_genes, n_var_columns)`; required when `gene_expression` is provided
- **Subset genes (`gene_identifiers`)**:
  - character vector of gene IDs to export (reduces disk use)
- **Connectivity (`connectivities`)**:
  - matrix or `Matrix::Matrix`, shape `(n_cells, n_cells)`
  - exported as unique undirected edge pairs
- **Vector fields (`vector_fields`)**:
  - named list of per-cell vectors (1D/2D/3D), used for velocity/displacement overlays

## Accepted types (what “matrix-like” means)

Embeddings and latent space:
- base R `matrix`
- `data.frame` (converted via `as.matrix`)

Gene expression:
- base R `matrix`
- `Matrix` sparse matrices (e.g. `dgCMatrix`)

Connectivities:
- base R `matrix`
- `Matrix` sparse matrices (recommended)

## Shapes and orientation (the #1 source of bugs)

Cellucid’s export format assumes:

- **rows = cells**
- **columns = dimensions / genes / components**

Concretely:

| Input | Expected shape |
|---|---|
| `X_umap_2d` | `(n_cells, 2)` |
| `latent_space` | `(n_cells, n_latent_dims)` |
| `obs` | `n_cells` rows |
| `gene_expression` | `(n_cells, n_genes)` |
| `var` | `n_genes` rows |
| `connectivities` | `(n_cells, n_cells)` |
| vector field (2D) | `(n_cells, 2)` |

```{warning}
Many R containers store expression as **genes × cells**. For Cellucid export you almost always need to transpose:

- Seurat: `GetAssayData(seu, slot = "data")` is usually genes × cells → use `Matrix::t(...)`
- SingleCellExperiment: `assay(sce, "logcounts")` is genes × cells → use `Matrix::t(...)`
```

## Cell identity = row order

There is no separate “cell ID” file in the export folder. Instead:

> Cell `i` is “whatever is in row `i`”.

That means you must keep the **same cell order** across:
- embeddings
- `latent_space`
- `obs`
- gene expression rows
- connectivity matrix rows/cols
- vector fields

If your source object has cell IDs (Seurat/SCE), use them to explicitly align inputs before export.

## Common limitations (know these early)

- **4D embeddings are not supported** (`X_umap_4d` must be `NULL`).
- **Gene expression is “one file per gene”** in the export format; exporting 20k+ genes is possible but can be slow and creates many files.
- **Invalid characters in some keys become underscores** (obs column names and gene IDs are sanitized for filenames). This can cause filename collisions if two keys sanitize to the same string.
- **Vector field IDs are stricter**: they must already be filesystem-safe; otherwise export fails (see {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`).

## Recommended workflow decision tree

- If you have a **Seurat object** → start with {doc}`../e_integrations_recipes/01_seurat_recipe`
- If you have a **SingleCellExperiment** → start with {doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`
- If you have **raw matrices** → start with {doc}`../e_integrations_recipes/03_raw_matrices_and_data_frames_recipe`

Then:
- Export → {doc}`../c_data_preparation_api/index`
- Open in browser → {doc}`../d_viewing_loading/index`

## Next steps

- Want a copy/paste example? {doc}`04_quick_start_3_levels`
- Want the full API documentation? {doc}`../c_data_preparation_api/index`
