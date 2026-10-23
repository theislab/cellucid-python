# Gene expression matrix

**Audience:** everyone with gene expression
**Time:** 10 minutes
**Goal:** know what `gene_expression=` becomes and why it reads fast

`gene_expression` is the cells × genes matrix the viewer colours by: any scipy
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

Pass the matrix you want to colour by: raw counts, normalized values, or a
layer (`gene_expression=adata.layers["log1p"]`). Only one matrix is stored.

## Why gene-major

The viewer reads one gene at a time (colouring by it) or a batch of genes (an
analysis). In a CSC matrix a gene is one contiguous range of `data` and
`indices`, read with a single range request; in a CSR matrix the browser must
first scan every stored entry to build a gene index. On a 561,947-cell dataset
the first gene took 86 ms from the prepared file and 4,110 ms from the same
data stored CSR, over a simulated 80 ms / 50 MB/s link
({doc}`10_performance_and_storage`).

## Sizes

`X/data` and `X/indices` are cut into chunks of 65,536 values and deflated —
the indices after a byte shuffle, the values as they are, which is how real
expression values compress best; on the dataset above a gene cost 8.6 KB
stored (2.0 KB of values, 6.6 KB of indices), against 13.1 KB in the CSR
float64 input. The full layout
is in {doc}`09_cellucid_ready_anndata_specification`.

## Next

- {doc}`07_connectivities_knn_graph`
- {doc}`10_performance_and_storage`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/06_gene_expression_matrix`
