# Gene expression (`X`)

**Audience:** everyone who wants gene search and gene colouring
**Time:** 10 minutes
**Goal:** pass expression in the right orientation and know how it is stored.

## What to pass

`gene_expression` is an optional **cells x genes** matrix:

- a sparse `Matrix` matrix (`dgCMatrix` and any other sparse numeric class), or
- a dense numeric matrix, a dense `Matrix`, or a numeric `data.frame`.

It must have one row per cell, in the cell order of the embeddings. With
`data =`, an object's expression is read as it is (Seurat's `data` layer, else
`counts`; a SingleCellExperiment's `logcounts`, else `counts`, else the first
assay), with no transposing by you. Seurat and SingleCellExperiment keep
expression as genes x cells, so a matrix you take out of one yourself is
transposed first:

```r
expr <- Matrix::t(Seurat::GetAssayData(seu, layer = "data"))
```

A wrong row count is refused with a reminder:

```text
gene_expression has 2000 rows, but the embeddings have 5000 cells. gene_expression is cells x genes; transpose a genes x cells matrix with t().
```

A square matrix in the wrong orientation cannot be detected from its shape;
check `nrow(expr)` against your cell count yourself.

Which values to pass is your choice: the viewer shows what you store.
Normalized, log-transformed values (Seurat's `data` layer, a
SingleCellExperiment's `logcounts`) are the usual choice for colouring.

## How it is stored

Values must be finite and within the float32 range, and are stored as
float32 -- unless they are integer counts float32 cannot hold exactly: when
every value is a whole number and one exceeds 2^24 (16,777,216) in magnitude,
the whole `X` is stored as float64, so every count reads back as it was given.
The rule looks at the stored values (a sparse matrix's non-zeros, a dense
matrix's every value) of the genes that are stored, after `gene_identifiers`
has chosen them; Python's `cellucid.prepare()` applies the same rule. A pattern matrix (`ngCMatrix`, no values) is refused.

**Sparse input** is stored as `X`, an AnnData `csc_matrix`: column-compressed
with one column per gene, so each gene's values are one contiguous range.

- `X/data`: float32 values (float64 for counts above 2^24, as above)
- `X/indices`: int32 cell indices
- `X/indptr`: int32 offsets, one per gene plus one

(indices and offsets become int64 only if a matrix ever held 2^31 or more
stored values.) `data` and `indices` are chunked at 65,536 values, so reading
one gene fetches little more than that gene. A stored zero is no value and
is dropped (as `Matrix::drop0()` drops it); the matrix stays sparse.

**Dense input** is stored as a dense float32 (or float64, as above) `X` array in tall, gene-major
chunks: each chunk holds up to 262,144 rows of one or a few genes.

See {doc}`09_cellucid_ready_anndata_specification` for every detail.

## Gene names and subsets

Gene names come from `rownames(var)` or `colnames(gene_expression)`;
`gene_identifiers` stores a subset; `var_gene_id_column` picks the names the
viewer shows. See {doc}`05_var_gene_metadata`.

## Size

Expression is usually the largest part of a store. Each chunk is
zlib-compressed; on simulated data, 10 million stored values took about 33 MB
({doc}`10_performance_and_storage`). To make a store smaller, store fewer genes
(`gene_identifiers`, for example highly variable genes) -- the cell count is
what it is.

## Without expression

`gene_expression` is optional. Without it the store has no `X` and an empty
`var`; the viewer shows embeddings and `obs` fields but no genes.

## Next

- {doc}`07_connectivities_knn_graph`
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/06_gene_expression_matrix`
