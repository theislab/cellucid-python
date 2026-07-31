# Input Requirements (Global)

**Audience:** computational users (recommended for everyone exporting real data)  
**Time:** 10–15 minutes  
**Goal:** prevent subtle exports that “work” but load incorrectly in the viewer.

This page documents the global rules that apply to *every* input you pass to `cellucid_prepare()`.

## Rule 1: cell identity is the row order

Cellucid’s export format does not store separate cell IDs. Instead:

> Cell `i` is the `i`-th row in every exported array.

This means you must keep a consistent row order across:

- `X_umap_1d` / `X_umap_2d` / `X_umap_3d`
- `latent_space`
- `obs`
- `gene_expression` rows (if provided)
- `connectivities` rows/cols (if provided)
- every vector field (if provided)

### Practical strategy (recommended)

If your source data has cell IDs (Seurat/SCE), pick a canonical ordering:

```r
cell_ids <- colnames(seu)  # Seurat
# or:
cell_ids <- colnames(sce)  # SingleCellExperiment
```

Then reorder everything to match `cell_ids` **explicitly** before exporting.

## Rule 2: required shapes

At least one embedding is required:

| Argument | Shape |
|---|---|
| `X_umap_1d` | `(n_cells, 1)` |
| `X_umap_2d` | `(n_cells, 2)` |
| `X_umap_3d` | `(n_cells, 3)` |

And:

| Argument | Shape |
|---|---|
| `latent_space` | `(n_cells, n_latent_dims)` |
| `obs` | `n_cells` rows |
| `gene_expression` | `(n_cells, n_genes)` |
| `var` | `n_genes` rows |
| `connectivities` | `(n_cells, n_cells)` |

```{warning}
If you pass `gene_expression` in the common “genes × cells” orientation, export will fail (shape mismatch) or silently produce nonsense.
Always ensure **cells × genes** for Cellucid export.
```

## Rule 3: no missing values in embeddings

Embeddings are validated before normalization. `NA`, `NaN`, infinities, and
nonpositive coordinate ranges reject the complete candidate.

Recommendation:
- remove cells with missing embedding coordinates or recompute the embedding
  from reviewed input.

## Rule 4: careful with non-numeric `obs` columns

`obs` is a data.frame. `cellucid-r` classifies columns as:

- **continuous**: `is.numeric(x)` is `TRUE`
- **categorical**: factors, logicals, and all other types (including character)

This means:
- a character column like `"sample_id"` becomes categorical (good),
- but a `Date` column becomes categorical (probably not what you want),
- and a numeric-looking character column becomes categorical unless you convert it.

Recommendation:
- explicitly coerce `obs` columns you care about (`as.numeric`, `factor`, etc.)

Details: {doc}`04_obs_cell_metadata`

## Rule 5: avoid filename collisions

Observation keys, gene IDs, dataset IDs, and vector-field IDs are used exactly
in manifests and paths. They must already be 1–180-byte portable ASCII
components: start with a letter or digit, contain only letters, digits, `.`,
`_`, or `-`, not end with `.`, not be a Windows device name, and be unique
under case-insensitive comparison.

This applies to the identifiers you **export**, because the rule is about the
file each one is written to. `obs_keys` and `gene_identifiers` narrow the export,
and they narrow this rule with it: a column or gene you leave out is written to
no path and is not checked here. Gene identifiers carry one extra rule that is
not about paths and so is *not* narrowed — every ID in `var` must be a non-empty
string and must be distinct, because `gene_identifiers` addresses `var` rows by
identifier. The `cellucid` Python package applies the identical scopes.

Recommendation:
- choose stable portable identifiers before export.
- `cellucid-r` aborts the complete candidate on any unsafe identifier or
  collision among the exported ones; it does not rewrite names.

## Quick “preflight” checks (copy/paste)

Use this pattern before you export real data.

```r
stopifnot(is.matrix(latent_space) || is.data.frame(latent_space))
stopifnot(is.data.frame(obs))
stopifnot(is.matrix(X_umap_2d) && ncol(X_umap_2d) == 2)

n_cells <- nrow(X_umap_2d)
stopifnot(nrow(latent_space) == n_cells)
stopifnot(nrow(obs) == n_cells)

if (!is.null(gene_expression)) {
  stopifnot(nrow(gene_expression) == n_cells)
  stopifnot(!is.null(var))
  stopifnot(nrow(var) == ncol(gene_expression))
}

if (!is.null(connectivities)) {
  stopifnot(nrow(connectivities) == n_cells)
  stopifnot(ncol(connectivities) == n_cells)
}
```

If you have cell IDs, add explicit alignment checks:

```r
stopifnot(!is.null(rownames(X_umap_2d)))
stopifnot(!is.null(rownames(latent_space)))
stopifnot(!is.null(rownames(obs)))
stopifnot(identical(rownames(X_umap_2d), rownames(obs)))
stopifnot(identical(rownames(latent_space), rownames(obs)))
```

## Troubleshooting pointers

- Export fails with “n_cells mismatch” → {doc}`11_troubleshooting_prepare_export`
- Viewer loads but fields look wrong (often a row-order bug) → {doc}`../i_troubleshooting_index/03_export_format_and_validation_issues`
