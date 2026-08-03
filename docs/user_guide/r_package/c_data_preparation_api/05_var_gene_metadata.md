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

```{warning}
A `data.frame` you never gave row names to is not row-name-less: R reports the
automatic sequence `"1"` … `"n_genes"`, which would publish an export whose genes
are named after their own row numbers. `cellucid_prepare()` refuses that rather
than let a reader search `CD8A` and find nothing:

    var has only automatic row names, so rownames(var) would name the genes '1'
    to '2000'. Set rownames(var) to the gene identifiers, or pass
    var_gene_id_column.
```

Common choices:
- gene symbols (`MS4A1`, `CD3D`, …)
- accessions (`ENSG00000156738`, …)

Whichever you choose is the **only** identity the export records, and it is
exactly what a reader types into the gene search box. An accession-keyed export
therefore matches nothing when a reader searches `MS4A1`. If you want symbols to
be searchable, resolve them into a `var` column yourself and select that column;
`cellucid_prepare()` never performs a symbol lookup and ships no mapping.

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
  `gene_identifiers contains identifiers not found in var: c("MS4A1").`, which
  lists the first five as an R vector and appends ` and N more` when there are
  more
- a repeated identifier stops the export with
  `gene_identifiers must not contain duplicate identifiers.`
- nothing is skipped and no partial export is written

Intersect your panel with your gene IDs before the call when a marker list may
not apply to the dataset:

```r
markers <- intersect(markers, rownames(var))
```

This is the single best lever for reducing disk size for large datasets.

## A gene ID is a name, not a filename

Gene IDs are used in exactly one place: `var_manifest.json`, as the **true gene
key**. They never become paths. Each gene's expression payload is written to
`var/<index>.values.*`, where `<index>` is the gene's position in the manifest's
`fields` array, so a gene name has only to be readable and unambiguous.

That leaves two rules, with different scopes because they answer different
questions:

| Rule | Scope | Why |
| --- | --- | --- |
| IDs are **distinct** | every row of `var` | `gene_identifiers` addresses rows by ID, so a repeated ID names no single row, and the export stops with `Gene key 'MS4A1' is duplicated.` |
| IDs are non-empty and **drawable verbatim** — no control or zero-width characters, no leading or trailing whitespace | only the genes actually exported | being drawable is a property of a name the viewer shows, and a deselected gene reaches no manifest and no legend |

A `var` that carries `HLA-DRB1/2`, `Wnt/β-catenin target`, or a gene name with an
interior space exports cleanly. There is nothing to rename and nothing to escape.
`prepare()` in the Python package applies the identical rules and scopes.

There is also only **one** identity per gene. Whatever `var_gene_id_column`
selects is the name recorded in the manifest and the name a reader searches for.
No accession is kept beside it, and there is no separate display-name argument:
if you want symbols in the viewer, put symbols in the column you select.

```{warning}
A duplicate gene ID anywhere in `var` rejects the complete candidate before
publication, even when `gene_identifiers` leaves that row out.
```

Practical recommendation:
- check `anyDuplicated(gene_ids) == 0` over every row of `var`
- check `identical(gene_ids, trimws(gene_ids))` over the IDs you export
- resolve accessions to symbols *before* the call if that is what you want
  searchable; `cellucid_prepare()` performs no lookup and ships no mapping

## Troubleshooting pointers

- “var has X rows but gene_expression has Y genes” → your expression orientation is wrong or var is mismatched.
- “My gene IDs are rejected” → check `var_gene_id_column`, duplicates, and
  invisible characters in the IDs you export.
- See also: {doc}`06_gene_expression_matrix` and {doc}`11_troubleshooting_prepare_export`
