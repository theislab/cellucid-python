# `var`: gene metadata and gene names

**Audience:** everyone storing gene expression
**Time:** 10 minutes
**Goal:** get the genes named the way you want to search for them in the viewer.

## What to pass

`var` is an optional `data.frame` with **one row per column of
`gene_expression`**, in the same order. Every column is stored, with the same
encodings as `obs` columns ({doc}`04_obs_cell_metadata`), and column names
follow the same path rules ({doc}`02_input_requirements_global`).

`var` describes the genes of `gene_expression`, so passing `var` (or
`gene_identifiers`, or `var_gene_id_column`) without `gene_expression` is an
error.

## The var index (gene names in the store)

The var index comes from, in this order:

1. `rownames(var)`, when you set them (a data frame's automatic `"1"`, `"2"`,
   ... row names are row numbers, not gene names, and do not count);
2. else `colnames(gene_expression)`;
3. else the values of `var_gene_id_column`.

If both `rownames(var)` and `colnames(gene_expression)` are set they must be
the same names in the same order:

```text
rownames(var) and colnames(gene_expression) disagree at gene 1 ('g1' and 'ENSG1'). var must describe the columns of gene_expression in their order.
```

With none of the three the call stops:

```text
The genes have no names: set rownames(var) or colnames(gene_expression) to the gene identifiers, or pass var_gene_id_column.
```

Gene names must not be missing or empty. Repeated names are made unique the way
anndata's `var_names_make_unique()` does it -- the first occurrence keeps its
name, the next ones become `"<name>-1"`, `"<name>-2"`, skipping names already
taken -- and a message says what changed:

```text
Gene names: 1 duplicate name was made unique, e.g. 'MATR3' → 'MATR3-1'.
```

Without `var`, the store holds a `var` data frame with the gene names and no
columns.

## `var_gene_id_column`: the names the viewer shows

Often the row names are stable ids (Ensembl) and a column holds the symbols
people search for. Name that column:

```r
var <- data.frame(symbol = gene_symbols, row.names = ensembl_ids)

cellucid_prepare(..., var = var, gene_expression = expr,
                 var_gene_id_column = "symbol")
```

The store keeps the Ensembl ids as the var index and the `symbol` column as it
is; `uns['cellucid']['defaults']['geneIdColumn'] = "symbol"` tells the viewer to
name genes by that column. The column must exist in `var` and is checked by the
rule the viewer reads it with:

- text and factor values name genes as they are;
- integer values name genes by their digits (`7157L` is `"7157"`);
- a double or logical column cannot name genes and is refused
  (`var column 'score': floating-point values cannot name genes.`);
- an empty value (`NA` or `""`) takes the row's var index name;
- repeated names are made unique as above.

The column is stored unchanged, and the viewer applies the last two rules when
it reads it. When they change a name, a message says so ahead of time:

```text
Gene names from var column 'symbol': 1 gene has no 'symbol' value and uses its var index name. 1 duplicate name was made unique, e.g. 'CD3E' → 'CD3E-1'.
```

Other var columns stay available in the viewer under **Session › Data mapping**
({doc}`../../web_app/b_data_loading/13_data_mapping`).

## `gene_identifiers`: store a subset of genes

```r
gene_identifiers = c("CD3E", "MS4A1", "NKG7")
```

stores only those genes, in that order: `X` and `var` are narrowed together.
The identifiers are matched against the names the viewer shows --
`var_gene_id_column` when set, else the var index -- **after** repeated names
were made unique, so the second `"CD3E"` is addressed as `"CD3E-1"`. They must
be distinct, and every one must exist:

```text
gene_identifiers contains identifiers not found in var: c("NOPE"). Identifiers are matched against the gene names the viewer shows (var_gene_id_column, or the var row names), after repeated names were made unique.
```

## Next

- {doc}`06_gene_expression_matrix`
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/05_var_gene_metadata`
