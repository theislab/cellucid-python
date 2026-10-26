# `var` (gene metadata) and gene names

**Audience:** everyone with gene expression
**Time:** 10 minutes
**Goal:** name the genes, choose the names the viewer shows, and keep a subset

`var` is a pandas DataFrame with one row per column of `gene_expression`, in
the same order. Every column is stored, with the encodings of
{doc}`04_obs_cell_metadata`. `var` without `gene_expression` is refused.

(python_package-make-unique-gene-names)=
## Gene names

- `var.index` names the genes. Usually these are stable ids (Ensembl) or
  symbols; repeated names are made unique (`MATR3`, `MATR3-1`) and reported.
- With pandas' default `0 … n-1` index, `var_gene_id_column` names the genes;
  with neither, the call is refused:

  ```text
  The genes have no names: give var an index of gene identifiers, or pass var_gene_id_column naming the var column that holds them.
  ```

## The names the viewer shows: `var_gene_id_column=`

When the index holds ids and a column holds symbols, name that column:

```python
cellucid.prepare(
    "data.h5ad",
    format="h5ad",
    embeddings={"X_umap": umap},
    dataset_id="pbmc-3k",
    dataset_name="PBMC 3k",
    var=adata.var,                    # index: Ensembl ids, column "gene_symbols"
    gene_expression=adata.X,
    var_gene_id_column="gene_symbols",
)
```

The column is stored as it is and named in
`uns['cellucid']['defaults']['geneIdColumn']`. The viewer reads it, fills an
empty value with the row's index name, and makes repeated names unique; what
that changes is printed when the file is written:

```text
Gene names from var column 'gene_symbols': 3 duplicate names were made unique, e.g. 'MATR3' → 'MATR3-1'.
```

The column must exist in `var` and hold text, categories or integers (a float
or boolean column cannot name genes). The viewer can switch to any other `var`
column or to the index later, in the data mapping
({doc}`../../web_app/b_data_loading/13_data_mapping`).

## Keeping some genes: `gene_identifiers=`

```python
cellucid.prepare(
    "data.h5ad",
    format="h5ad",
    embeddings={"X_umap": umap},
    dataset_id="pbmc-3k",
    dataset_name="PBMC 3k",
    var=adata.var,
    gene_expression=adata.X,
    var_gene_id_column="gene_symbols",
    gene_identifiers=["CD3E", "MS4A1", "NKG7"],
)
```

Only those genes are written, in that order: the columns of `X` and the rows of
`var` together. Identifiers are matched against the names the viewer shows
(`var_gene_id_column`, else the index), after repeated names were made unique.
A duplicate or unknown identifier is refused, naming the first five unknown
ones.

## Next

- {doc}`06_gene_expression_matrix`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/05_var_gene_metadata`
