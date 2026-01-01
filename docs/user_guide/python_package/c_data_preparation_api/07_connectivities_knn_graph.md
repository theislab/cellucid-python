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
- a scipy sparse matrix (recommended), or
- any array-like that can be converted to CSR.

Row/column order must match your embeddings and `obs`.

### How the exporter interprets the matrix (important)

The current Python exporter:

1) Converts to CSR (if not already).
2) **Symmetrizes** the matrix: `C_sym = C + C.T`
3) **Binarizes** edges by setting all nonzero values to `1`
4) Extracts unique undirected edges by keeping only pairs where `src < dst`
5) Sorts edges by `(src, dst)` for better gzip compression
6) Writes two parallel arrays:
   - `edges.src.bin(.gz)`
   - `edges.dst.bin(.gz)`

Important implications:
- Any original edge **weights are discarded** (graph becomes unweighted).
- Directionality is discarded (graph becomes undirected).
- Self-edges are not exported.

### Dtype and scaling limits

Edge indices are stored using the smallest integer dtype that fits `n_cells`:
- `uint16` if `n_cells ≤ 65,535`
- `uint32` if `n_cells ≤ 4,294,967,295`
- `uint64` otherwise (practically unrealistic for interactive visualization)

The web app converts indices to `uint32` internally.

### Output files

When connectivities are exported, you get:

```text
out_dir/
├── connectivity_manifest.json
└── connectivity/
    ├── edges.src.bin(.gz)
    └── edges.dst.bin(.gz)
```

Full format spec: {doc}`09_output_format_specification_exports_directory`

---

## Performance considerations

Export time scales roughly with:
- number of cells (`n_cells`),
- and number of nonzero neighbor entries (`nnz`).

The exporter currently iterates through neighbors in Python to extract unique edges.

Practical tips:
- Prefer a sparse matrix input.
- Keep K reasonably sized (typical KNN graphs have 10–50 neighbors).
- For huge datasets, consider skipping connectivities in static exports and relying on other workflows.

---

## Edge cases and common footguns

- **Wrong shape**: any non-square matrix fails.
- **Row/col order mismatch**: graph edges connect the wrong cells (hard to detect visually).
- **Very dense graphs**: huge `nnz` → export time and output size explode.
- **Disconnected graphs**: valid export, but some features may behave unexpectedly.

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
- check that `connectivity/edges.src.bin(.gz)` exists

Fix:
- export connectivities and re-export with `force=True`.

---

## Next steps

- Optional vector fields: {doc}`08_vector_fields_velocity_displacement`
- Output directory format and manifests: {doc}`09_output_format_specification_exports_directory`
