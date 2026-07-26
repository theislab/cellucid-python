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
- values must be real numeric or logical, finite, and non-negative
- topology and weights must be exactly symmetric
- every diagonal value must be exactly `0`
- the `Matrix` package must be installed when the input is sparse
- sparse matrices must not contain stored zero entries or duplicate
  coordinates

```{note}
If you do not need connectivity features, skip this entirely. Exports load fine without it.
```

## Exact exported representation

Cellucid validates the graph without changing it. It does not symmetrize,
binarize, discard weights, coalesce coordinates, remove self-edges, or
otherwise reinterpret scientific input. Asymmetric, directed, negative,
non-finite, nonzero-diagonal, sparse-stored-zero, and duplicate-coordinate
graphs fail before graph artifacts are written. Logical `TRUE` is the explicit
unit edge and is written as the exact Float64 weight `1.0`.

For a valid graph, the exporter:

1) extracts each undirected edge once from the upper triangle (`src < dst`);
2) sorts the pairs lexicographically;
3) writes three aligned arrays:
   - `edges.src.bin`
   - `edges.dst.bin`
   - `edges.weights.f64.bin`

Indices are **zero-based** in the binary files.
Weights are exact little-endian Float64 values.
A valid graph with no edges remains present and reports `n_edges: 0` and
`max_neighbors: 0`; all three payloads are zero length.

## Output files

- `connectivity_manifest.json`
- `connectivity/edges.src.bin` (or `.gz`)
- `connectivity/edges.dst.bin` (or `.gz`)
- `connectivity/edges.weights.f64.bin` (or `.gz`)

The manifest records:
- `n_cells`
- `n_edges`
- `max_neighbors` (maximum degree in the validated graph)
- `index_dtype` / `index_bytes` (`uint16` or `uint32`)
- `weightsPath`, `weight_dtype: "float64"`, and `weight_bytes: 8`

## Index dtype and dataset size limits

Cell indices are stored using the smallest unsigned integer that fits:

| `n_cells` | dtype |
|---:|---|
| 1–65,536 | `uint16` |
| 65,537–4,294,967,296 | `uint32` |
| larger | rejected; no current `uint64` graph format |

R matrices have lower practical dimension and memory limits than the format's
mathematical `uint32` ceiling. Prefer sparse matrices and use server-backed
loading when a graph exceeds the browser working-set limit.

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

Rejected. Construct and review the intended symmetric weighted graph explicitly
before calling `cellucid_prepare()`.

### Weighted graphs

Accepted when weights are finite, non-negative, exactly symmetric, and exactly
representable as Float64. They are preserved in the aligned weight payload.

### Self-edges, negative values, or non-finite values

Rejected. The diagonal must be exactly zero and every value must be finite and
non-negative.

### Dense matrices

A dense `(n_cells, n_cells)` matrix becomes huge quickly. Prefer sparse matrices.

## Troubleshooting pointers

- “Matrix package is required to export connectivity matrices” → install `Matrix`.
- “Connectivity matrix shape ... does not match number of cells” → you have an ordering or subset mismatch.
- Viewer is slow/crashes after adding connectivities → your graph is too dense (too many edges).
