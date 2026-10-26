# Connectivities (neighbour graph)

**Audience:** computational users
**Time:** 10 minutes
**Goal:** store a neighbour graph the viewer can draw and use.

## What to pass

`connectivities` is an optional **cells x cells** matrix of edge weights:

- a `Matrix` sparse matrix (recommended; any numeric or pattern class -- a
  pattern matrix gives every edge weight 1), or
- a dense numeric or logical base matrix (needs the `Matrix` package installed).

It is stored as `obsp['connectivities']`, the name scanpy uses, so the viewer
finds it without configuration.

## Rules

- **Shape:** exactly cells x cells.
  `connectivities must be 3 x 3 (one row and one column per cell), got 2 x 2.`
- **Weights:** finite and non-negative.
- **No self-loops:** the diagonal must be exactly zero.
  `Connectivity diagonal values must all be exactly zero.`
- **Exactly symmetric** in topology and weight (an undirected graph).
  `Connectivity matrix must be exactly symmetric.`
- A stored zero is not an edge and is dropped (`Matrix::drop0()`).

With `data =`, the object's graph is chosen and its diagonal removed for you
(Seurat's `<active assay>_snn`, a SingleCellExperiment `colPair` named
`connectivities` or like `snn`; {doc}`../e_integrations_recipes/index`).
Graphs you pass yourself from common tools often need one step first:

```r
# Seurat: the shared-nearest-neighbour graph links each cell to itself.
graph <- seu@graphs[["RNA_snn"]]
Matrix::diag(graph) <- 0

# A directed k-nearest-neighbour graph: make it symmetric (stays sparse).
graph <- (knn + Matrix::t(knn)) / 2
```

Avoid `pmax(knn, t(knn))` on large sparse matrices: `Matrix` turns it into a
dense matrix, which for 100,000 cells needs tens of gigabytes.

## How it is stored

`obsp['connectivities']` is an AnnData `csr_matrix` of shape cells x cells:

- `data`: float32 weights
- `indices`: int32 column indices
- `indptr`: int32 row offsets

For a symmetric matrix the column-compressed arrays of R's `dgCMatrix` are
exactly the row-compressed arrays of the same matrix, so they are written
without reordering. Python reads it as `scipy.sparse.csr_matrix`.

## What the viewer does with it

The viewer reads the graph's edges when they are first needed, not when the
dataset opens, so a large graph does not slow opening. `cellucid_prepare()`
stores one graph, under `connectivities`.

## Next

- {doc}`08_vector_fields_velocity_displacement`
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/07_connectivities_knn_graph`
