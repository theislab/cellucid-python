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

## Rule 5: identifiers name fields, not files

Observation keys, gene IDs, and vector-field IDs are recorded exactly in the
manifests. They are **not** filenames: every payload file is named by an integer
index, and the manifest is what says which identifier that index belongs to.

So an identifier carries no filename rule at all — no ASCII restriction, no
length limit, no ban on `/`, spaces, or non-ASCII characters, no
case-insensitive collision rule, and no Windows device-name rule. What survives
is what the identity is actually for. Every exported identifier must:

1. be a **non-empty string**;
2. be **distinct within its axis**, so a lookup by identifier resolves exactly
   one field; and
3. be **text the viewer can draw exactly as it is stored** — no control
   characters (`U+0001`-`U+001F`, `U+007F`-`U+009F`), no zero-width characters
   (`U+200B`, `U+2060`, `U+FEFF`), and no leading or trailing whitespace of any
   kind, including `U+00A0` NO-BREAK SPACE.

Rule 3 is the same rule every category label obeys, because a gene name and a
category label are drawn in the same legend. `cellucid_prepare()` rejects the
complete candidate rather than trimming, since trimming would rewrite an
annotation you did not ask to change. The `cellucid` Python package enforces the
identical rules.

Rules 1 and 3 are checked on what you actually **export**: `obs_keys` and
`gene_identifiers` narrow the export and narrow those checks with it. Gene
uniqueness is not narrowed — every ID in `var` must be distinct, because
`gene_identifiers` addresses `var` rows by identifier whether or not that row is
exported.

`dataset_id` is the one identifier that really is a path segment: it names the
dataset's directory and appears in URLs. It must still be a 1–180-byte portable
ASCII component — start with a letter or digit, otherwise only letters, digits,
`.`, `_`, or `-`, not ending with `.`, and not a Windows device name.

Recommendation:
- choose a stable, versioned `dataset_id` before export.
- `cellucid-r` aborts the complete candidate on a duplicate or undrawable
  identifier; it does not rewrite names.

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
