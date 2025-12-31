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

### Default: `"index"` uses `rownames(var)`

If `var_gene_id_column = "index"` (default):
- gene IDs are taken from `rownames(var)`
- if `rownames(var)` is `NULL`, gene IDs fall back to `"0".."n_genes-1"`

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
- missing gene IDs trigger a warning and are skipped
- export proceeds with the intersection

This is the single best lever for reducing disk size for large datasets.

## Filename safety and collisions

Gene IDs are used in two places:

1) In the manifest (`var_manifest.json`) as the **true gene key**
2) As part of the output filename under `var/`

For filenames, the exporter sanitizes gene IDs:
- unsupported characters become underscores
- leading/trailing dots/underscores are removed

```{warning}
If two different gene IDs sanitize to the same filename, one will overwrite the other.
`cellucid-r` does not currently detect this collision for you.
```

Practical recommendation:
- ensure `unique(safe_gene_ids)` is true, where `safe_gene_ids` is your gene IDs after sanitization
- avoid gene IDs containing `/`, `\\`, whitespace, or trailing punctuation

## Troubleshooting pointers

- “var has X rows but gene_expression has Y genes” → your expression orientation is wrong or var is mismatched.
- “My gene search is missing genes” → check `var_gene_id_column`, duplicates, and sanitization collisions.
- See also: {doc}`06_gene_expression_matrix` and {doc}`11_troubleshooting_prepare_export`
