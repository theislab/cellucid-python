# Connectivities / KNN Graph

**Audience:** computational users  
**Time:** 10–20 minutes  
**Goal:** export a cell–cell graph safely (without shape/order bugs).

`connectivities` is an optional `(n_cells, n_cells)` adjacency/connectivity matrix.

In the Cellucid ecosystem, this is typically a KNN/SNN graph derived from:
- scRNA-seq neighbor graphs,
- integrated embeddings,
- or any graph you want to visualize/use for neighbor-based interactions.

## Requirements

- shape must be `(n_cells, n_cells)`
- rows/cols must correspond to the same cell ordering as embeddings/obs
- input can be:
  - a base R `matrix`, or
  - a `Matrix::Matrix` sparse matrix (recommended)
- the `Matrix` package must be installed to export connectivities

```{note}
If you do not need connectivity features, skip this entirely. Exports load fine without it.
```

## What `cellucid-r` does with the matrix

The exporter turns your connectivity matrix into an **undirected, unweighted edge list**:

1) Converts to a sparse matrix (if needed).
2) Symmetrizes: `conn_sym = conn + t(conn)`.
3) Binarizes: any non-zero becomes `1`.
4) Extracts edges and keeps only unique undirected pairs (`src < dst`).
5) Writes two parallel arrays:
   - `edges.src.bin`
   - `edges.dst.bin`

Indices are **zero-based** in the binary files.

## Output files

- `connectivity_manifest.json`
- `connectivity/edges.src.bin` (or `.gz`)
- `connectivity/edges.dst.bin` (or `.gz`)

The manifest records:
- `n_cells`
- `n_edges`
- `max_neighbors` (max degree after symmetrization)
- `index_dtype` / `index_bytes` (`uint16`, `uint32`, or `uint64` depending on `n_cells`)

## Index dtype and dataset size limits

Cell indices are stored using the smallest unsigned integer that fits:

| `n_cells` | dtype |
|---:|---|
| ≤ 65,535 | `uint16` |
| ≤ 4,294,967,295 | `uint32` |
| larger | `uint64` |

## Practical recipes (where to get a graph)

### Seurat

Seurat often stores graphs in `seu@graphs`, for example:
- `seu@graphs$RNA_snn`

You must ensure it matches the cell order you export:

```r
cells <- colnames(seu)
conn <- seu@graphs$RNA_snn
conn <- conn[cells, cells]
```

Then pass `connectivities = conn`.

### SingleCellExperiment

Graph storage varies by workflow. Common patterns include:
- adjacency matrices stored in `metadata(sce)`
- KNN graphs stored as edge lists / igraph objects

If you can obtain an adjacency matrix `conn` with shape `(n_cells, n_cells)`, you can export it.

## Edge cases

### Directed or asymmetric graphs

Export always symmetrizes. If you pass a directed graph, the output becomes undirected.

### Weighted graphs

Weights are discarded (non-zero becomes 1). If you need weights for analysis, keep them elsewhere.

### Dense matrices

A dense `(n_cells, n_cells)` matrix becomes huge quickly. Prefer sparse matrices.

## Troubleshooting pointers

- “Matrix package is required to export connectivity matrices” → install `Matrix`.
- “Connectivity matrix shape ... does not match number of cells” → you have an ordering or subset mismatch.
- Viewer is slow/crashes after adding connectivities → your graph is too dense (too many edges).
