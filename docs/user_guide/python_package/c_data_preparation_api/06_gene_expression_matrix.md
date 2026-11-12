# Gene expression matrix

**Audience:** everyone with gene expression
**Time:** 10 minutes
**Goal:** know what `gene_expression=` becomes and why it reads fast

`gene_expression` is the cells × genes matrix the viewer colors by: any scipy
sparse format, or a dense NumPy array. It is stored as `X`.

## What is stored

| Input | `X` |
|---|---|
| sparse (CSR, CSC, COO, sparse arrays) | `csc_matrix`: gene-major, sorted indices, duplicates summed, stored zeros dropped |
| dense | a dense array, chunked in tall strips of genes |

Values are stored as **float32**. Integer counts are exact in float32 up to
2^24 (16,777,216); when any integer count of the genes stored (after
`gene_identifiers`) is larger, the whole matrix is stored as float64 so every
count reads back as given. Fractional values
(normalized, log-transformed) are rounded to float32 once, when written. Every
value must be finite and within float32's range.

Pass the matrix you want to color by: raw counts, normalized values, or a
layer (`gene_expression=adata.layers["log1p"]`). Only one matrix is stored.
Raw counts keep every choice open: the web app's **Session › Data mapping ›
Expression values** normalizes them to 10,000 counts per cell and takes
`log1p` for coloring, filters and analyses when asked, and its pseudobulk
differential expression sums raw counts per sample.

(total-counts-prepare)=
## Each cell's total count

Normalizing a cell's counts divides by its total count (its library size).
When `X` holds counts (every value a whole number of at least 0),
`prepare()` writes each cell's sum over the genes stored as
`obs['total_counts']` (float64) and records it in
`uns['cellucid']['totalCounts']` as `{"field": "total_counts", "over":
"stored genes"}`, so the viewer reads one column instead of every gene.

Totals over a gene panel are not the cell's depth: after
`sc.pp.highly_variable_genes(..., subset=True)` the stored genes are a few
thousand of the transcriptome. Compute the totals **before** selecting genes
and name the column:

```python
sc.pp.calculate_qc_metrics(adata, inplace=True)   # obs['total_counts'], all genes
adata = adata[:, adata.var.highly_variable].copy()
cellucid.prepare("pbmc.zarr", format="zarr", data=adata, total_counts="total_counts")
```

A named column is recorded `"over": "all genes"` and must hold, for every
cell, a finite total at least the sum of its stored genes (a log of the
totals is refused). A column of the data's own named `total_counts` that the
argument does not name is kept as given, and recorded only when it is exactly
the stored genes' sums. With a MuData, `total_counts={"rna": "total_counts"}`
names a column of the MuData's `obs` per modality (or `"rna:<column>"`, one of
the modality's own `obs`); without it each counts modality gets its own
`mod/<name>/obs/total_counts`.

When the source itself measured a gene panel -- a targeted panel, or an
integration's shared genes (the built-in atlas samples' 8,192) -- the totals
cover that panel, not the transcriptome. Say so, and the viewer says so where
it names the totals:

```python
cellucid.prepare("atlas.zarr", format="zarr", data=adata,
                 total_counts="total_counts", total_counts_over="source panel")
```

The column is recorded `"over": "source panel"` and checked as an `"all
genes"` one is (with a MuData, one value, or `{"rna": "source panel"}`). The
contract: {ref}`total-counts-spec`.

## Why gene-major

The viewer reads one gene at a time (coloring by it) or a batch of genes (an
analysis). In a CSC matrix a gene is one contiguous range of `data` and
`indices`, read with a single range request; in a CSR matrix the browser must
first scan every stored entry to build a gene index. On a 561,947-cell dataset
the first gene took 86 ms from the prepared file and 4,110 ms from the same
data stored CSR, over a simulated 80 ms / 50 MB/s link
({doc}`10_performance_and_storage`).

## Sizes

`X/data` and `X/indices` are cut into chunks of 16,384 values and deflated —
the indices after a byte shuffle, the values as they are, which is how real
expression values compress best; on the dataset above a gene cost 8.6 KB
stored (2.0 KB of values, 6.6 KB of indices), against 13.1 KB in the CSR
float64 input. The full layout
is in {doc}`09_cellucid_ready_anndata_specification`.

## Next

- {doc}`07_connectivities_knn_graph`
- {doc}`10_performance_and_storage`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/06_gene_expression_matrix`
