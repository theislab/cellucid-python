# Connectivities (KNN graph)

**Audience:** computational users (optional feature; most datasets can be explored without it)  
**Time:** 15–45 minutes depending on dataset size  
**Goal:** export a cell–cell neighborhood graph in a format the viewer can load efficiently

`connectivities` is an optional cell × cell matrix that encodes a KNN-style graph (neighbors between cells).

In AnnData, this is typically:
- `adata.obsp["connectivities"]` (often from Scanpy’s neighbor graph).

---

## Fast path (when to skip)

Skip `connectivities` if:
- you don’t need any neighbor-graph-driven features,
- you’re exporting a very large dataset and export time is a concern,
- you are still iterating on embeddings/obs and want a minimal export first.

Add it later once the rest of your export pipeline is stable.

---

## Practical path (computational users)

### Expected type and shape

`connectivities` must be a square matrix:

```text
(n_cells, n_cells)
```

It can be:
- a SciPy sparse matrix (recommended), or
- a dense real numeric or boolean array.

Row/column order must match your embeddings and `obs`.

The matrix is accepted only when all of these conditions hold:

- its shape is exactly `(n_cells, n_cells)`;
- every value is finite and non-negative;
- topology and weights are exactly symmetric;
- every diagonal value is exactly `0`.

Sparse inputs must not store explicit zero entries or duplicate coordinates.
Cellucid does not symmetrize, binarize, drop weights, coalesce coordinates,
remove self-edges, or otherwise reinterpret the graph. Asymmetric, directed,
negative, non-finite, nonzero-diagonal, and structurally ambiguous inputs fail
before export. Boolean `True` is the explicitly supported unit edge and is
written as the exact Float64 weight `1.0`.

### Exact exported representation

After validation, the exporter:

1) extracts each undirected edge once from the upper triangle (`src < dst`);
2) sorts edges lexicographically by `(src, dst)`;
3) writes three aligned arrays:
   - `edges.src.bin(.gz)`
   - `edges.dst.bin(.gz)`
   - `edges.weights.f64.bin(.gz)`

The weight payload is little-endian Float64 and preserves each accepted input
weight exactly. The source matrix is not modified. A valid graph with zero
edges remains a present graph: its manifest reports `n_edges: 0` and
`max_neighbors: 0`, and all three payloads are zero length.

### Dtype and scaling limits

Edge indices are stored using the smallest integer dtype that fits `n_cells`:

- `uint16` for `1 ≤ n_cells ≤ 65,536` (indices `0..65,535`);
- `uint32` for `65,537 ≤ n_cells ≤ 4,294,967,296`;
- larger cell axes are rejected because the current format has no `uint64`
  connectivity representation.

The mathematical `uint32` ceiling does not imply that a browser can hold a
graph of that size. Cellucid preflights the required working set and directs
oversized interactive workloads to server-backed loading instead of attempting
an unsafe allocation.

### Output files

When connectivities are exported, you get:

```text
out_dir/
├── connectivity_manifest.json
└── connectivity/
    ├── edges.src.bin(.gz)
    ├── edges.dst.bin(.gz)
    └── edges.weights.f64.bin(.gz)
```

Full format spec: {doc}`09_output_format_specification_exports_directory`

---

## Performance considerations

Export time scales roughly with:
- number of cells (`n_cells`),
- and number of nonzero neighbor entries (`nnz`).

Practical tips:
- Prefer a sparse matrix input.
- Keep K reasonably sized (typical KNN graphs have 10–50 neighbors).
- For huge datasets, consider skipping connectivities in static exports and relying on other workflows.

---

## Edge cases and common footguns

- **Wrong shape**: any non-square matrix fails.
- **Row/col order mismatch**: graph edges connect the wrong cells (hard to detect visually).
- **Negative weights, asymmetry, directionality, or self-edges**: rejected
  instead of transformed.
- **Sparse stored zeros or duplicate coordinates**: rejected instead of
  coalesced or discarded.
- **Very dense graphs**: huge `nnz` → export time and output size explode.
- **Disconnected or empty graphs**: valid when the exact matrix contract holds.

---

## Troubleshooting (connectivities)

### Symptom: export errors with “Connectivity matrix shape … does not match number of cells”

Meaning:
- Your connectivities matrix is not aligned to the embedding/obs row order or was computed on a different subset.

Fix:
- recompute connectivities on the exact same cell set,
- or subset/reorder connectivities to match `n_cells`.

### Symptom: features that depend on the graph are disabled in the viewer

Meaning:
- `connectivity_manifest.json` wasn’t exported or wasn’t loaded.

How to confirm:
- check for `<out_dir>/connectivity_manifest.json`
- check that the manifest-declared source, destination, and Float64 weight
  payloads all exist

Fix:
- pass a valid graph and export to a fresh output directory.

---

## Next steps

- Optional vector fields: {doc}`08_vector_fields_velocity_displacement`
- Output directory format and manifests: {doc}`09_output_format_specification_exports_directory`
