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
- `rownames(var)` must contain exact non-empty strings (R row names are always
  unique, so this set needs no renaming)

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

Whichever you choose is the **default** name set: what a reader sees first and
types into the gene search box. To offer both, keep one as the default and add
the other with `var_gene_name_columns` (below). `cellucid_prepare()` never
performs a symbol lookup and ships no mapping: resolve accessions to symbols in
a `var` column yourself.

### Custom column

If your gene IDs live in a column:

```r
cellucid_prepare(..., var_gene_id_column = "gene_id")
```

Export fails if the column does not exist. An empty value in the column (`""`
or `NA`) takes the gene's name from `rownames(var)`, and a repeated value is made
unique (see below).

### Alternative name sets (`var_gene_name_columns`)

```r
cellucid_prepare(
  ...,
  var_gene_id_column = "gene_symbols",
  var_gene_name_columns = c("gene_ids")
)
```

Each column listed becomes one more set of names for the same genes, stored as
metadata in `var_manifest.json` under `gene_name_sets`. No expression payload is
duplicated: every set names the same `var/` files, so the viewer can switch sets
without reloading expression. A column must exist, may be listed once, and may
not be the `var_gene_id_column` itself:

```text
var_gene_name_columns lists 'sym', which var_gene_id_column already uses for the default gene names.
```

To offer the row names as an alternative while a column is the default, copy
them into a column first: `var$ensembl_id <- rownames(var)`.

### Empty and repeated names: the make-unique rule

Every name set is completed by one rule, shared by `cellucid_prepare()`,
`cellucid.prepare()` in Python, the Python AnnData server, and the viewer — the
rule anndata's `var_names_make_unique()` applies:

1. An empty value (`""` or `NA`) takes the row's name from `rownames(var)`.
2. The first occurrence of a name keeps it.
3. Each later occurrence becomes `name-1`, `name-2`, …, skipping any value that
   already exists, so `A, A, A-1, A` becomes `A, A-2, A-1, A-3`.

The rule runs over the **whole** `var`, so a gene keeps its name however you
narrow the export, and `gene_identifiers` is matched against the unique names. A
set that was already unique is unchanged. What the rule changed is recorded per
set in `var_manifest.json` (`filled_from_index`, `renamed`). The cases every
implementation must agree on are pinned in a conformance fixture that both
packages' test suites run.

## Subsetting genes (`gene_identifiers`)

If you do not want to export all genes, pass a character vector:

```r
markers <- c("MS4A1", "CD3D", "NKG7")
cellucid_prepare(..., gene_identifiers = markers)
```

Behavior:
- every requested identifier must be present among the default gene names,
  after repeated names were made unique (the second `MATR3` is `MATR3-1`)
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

Gene IDs are used in exactly one place: `var_manifest.json`. They never become
paths. Each gene's expression payload is written to `var/<index>.values.*`,
where `<index>` is the gene's position in the manifest's `fields` array, so a
gene name has only to be readable.

The exported names — in every name set — must be non-empty and **drawable
verbatim**: no control or zero-width characters, no leading or trailing
whitespace. A name that breaks this rejects the complete candidate before
publication. A gene `gene_identifiers` leaves out reaches no manifest and is not
held to the rule.

A `var` that carries `HLA-DRB1/2`, `Wnt/β-catenin target`, or a gene name with an
interior space exports cleanly. There is nothing to rename and nothing to escape.
`prepare()` in the Python package applies the identical rules and scopes.

Practical recommendation:
- check `identical(gene_ids, trimws(gene_ids))` over the IDs you export
- resolve accessions to symbols *before* the call if that is what you want
  searchable; `cellucid_prepare()` performs no lookup and ships no mapping

## Troubleshooting pointers

- “var has X rows but gene_expression has Y genes” → your expression orientation is wrong or var is mismatched.
- “My gene IDs are rejected” → check `var_gene_id_column` and invisible
  characters in the IDs you export.
- “A gene search finds `MATR3` but not my second `MATR3`” → it is exported as
  `MATR3-1`.
- See also: {doc}`06_gene_expression_matrix` and {doc}`11_troubleshooting_prepare_export`
