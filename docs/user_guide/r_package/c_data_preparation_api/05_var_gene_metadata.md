# `var`: Gene / Feature Metadata

**Audience:** computational users  
**Time:** 10–15 minutes  
**Goal:** ensure gene identifiers are correct, stable, and loadable in the viewer.

`var` is the gene/feature metadata table (analogous to `adata.var` in AnnData).

You only need `var` if you export gene expression (`gene_expression`).

## Required invariants

If `gene_expression` is provided:

- `var` must be provided
- `nrow(var) == ncol(gene_expression)`

```{warning}
`gene_expression` must be shaped `(n_cells, n_genes)`, so `var` must have `n_genes` rows.
This is the reverse of many R containers where expression is stored as genes × cells.
```

## Choosing gene identifiers (`var_gene_id_column`)

The exporter needs a **gene ID string** for each column in `gene_expression`.

### `NULL` uses `rownames(var)`

If `var_gene_id_column = NULL`:
- gene IDs are taken from `rownames(var)`
- `rownames(var)` must contain exact non-empty strings

Recommendation:
- always set `rownames(var)` to the gene ID you want users to search in the viewer

Common choices:
- gene symbols (`MS4A1`, `CD3D`, …)
- Ensembl IDs (`ENSG000001...`)

### Custom column

If your gene IDs live in a column:

```r
cellucid_prepare(..., var_gene_id_column = "gene_id")
```

Export fails if the column does not exist.

## Subsetting genes (`gene_identifiers`)

If you do not want to export all genes, pass a character vector:

```r
markers <- c("MS4A1", "CD3D", "NKG7")
cellucid_prepare(..., gene_identifiers = markers)
```

Behavior:
- every requested identifier must be present in the gene ID list
- a gene ID that is absent stops the export with
  `gene_identifiers contains identifiers not found in var: …`, which lists the
  first five and appends `...` when there are more
- a repeated identifier stops the export with
  `gene_identifiers must not contain duplicate identifiers.`
- nothing is skipped and no partial export is written
- narrowing the export narrows the filename contract with it: only the
  identifiers you actually export must be portable filename components (see
  below)

Intersect your panel with your gene IDs before the call when a marker list may
not apply to the dataset:

```r
markers <- intersect(markers, rownames(var))
```

This is the single best lever for reducing disk size for large datasets.

## Filename safety and collisions

Gene IDs are used in two places:

1) In the manifest (`var_manifest.json`) as the **true gene key**
2) As part of the output filename under `var/`

The exporter uses each gene ID exactly. Every **exported** ID must already
satisfy the portable identifier contract and be unique under case-insensitive
filesystem comparison.

The two identifier rules have different scopes, because they answer different
questions:

| Rule | Scope | Why |
| --- | --- | --- |
| IDs are non-empty strings and **distinct** | every row of `var` | `gene_identifiers` addresses rows by ID, so a repeated ID names no single row, and the export stops with `Gene identifiers must be unique.` |
| IDs are **portable filename components**, unique case-insensitively | only the genes actually exported | the rule is about the file the gene is written to, and a deselected gene is written to none |

So a `var` that carries `HLA-DRB1/2` still exports cleanly as long as
`gene_identifiers` leaves that gene out. This is the same scoping `obs_keys`
already has (see {doc}`04_obs_cell_metadata`). `prepare()` in the Python package
applies the identical scopes.

```{warning}
An unsafe ID or case-insensitive collision **among the exported genes** rejects
the complete candidate before publication.
```

Practical recommendation:
- avoid gene IDs containing `/`, `\\`, whitespace, trailing `.`, or non-ASCII
  characters
- check `anyDuplicated(tolower(gene_ids)) == 0` over the IDs you export

## Troubleshooting pointers

- “var has X rows but gene_expression has Y genes” → your expression orientation is wrong or var is mismatched.
- “My gene IDs are rejected” → check `var_gene_id_column`, portability,
  duplicates, and case-insensitive collisions.
- See also: {doc}`06_gene_expression_matrix` and {doc}`11_troubleshooting_prepare_export`
