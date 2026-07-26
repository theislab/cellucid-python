# AnnData mode (`show_anndata()` / `serve_anndata()`)

This page explains how Cellucid can visualize **AnnData directly**, without running `prepare()` first.

Use AnnData mode when:
- you’re iterating inside an analysis notebook,
- you want a quick “does this dataset look sane?” view,
- or you have a large `.h5ad` and want read-only-backed access without exporting.

If you want maximum performance + easiest sharing, prefer export mode:
{doc}`07_exported_directory_mode_show_and_serve`.

## At a glance

**Audience**
- Computational users working with Scanpy/AnnData
- Wet lab users who received a `.h5ad`/`.zarr` from a collaborator

**Prerequisites**
- `pip install cellucid`
- A `.h5ad` / `.zarr` / in-memory `AnnData`

## Quick start

### Notebook (recommended)

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
# A .zarr path or in-memory AnnData uses the same required identity arguments.
```

### Terminal

```bash
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1
```

### Python (script)

```python
from cellucid import serve_anndata

serve_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
)  # blocks
```

## Minimum AnnData requirements (the “must haves”)

### 1) UMAP embedding in `adata.obsm`

Cellucid currently expects UMAP coordinates under one of these keys:

- explicit (preferred):
  - `X_umap_1d` (shape `(n_cells, 1)`)
  - `X_umap_2d` (shape `(n_cells, 2)`)
  - `X_umap_3d` (shape `(n_cells, 3)`)

If no valid UMAP embedding is found, you’ll get an error like:
“No valid UMAP embeddings found in `adata.obsm` …”.

### 2) (Optional but common) Gene expression in `adata.X`

- If `adata.X` is missing, Cellucid can still show embeddings + obs fields.
- Gene search / expression coloring will not work without expression values.

### 3) (Optional) Cell metadata in `adata.obs`

Cellucid derives obs fields automatically:
- numeric → continuous
- categorical/string/bool → categorical

## UMAP key resolution rules (detailed, for debugging)

Cellucid reads only the explicitly dimensioned keys
`X_umap_1d`, `X_umap_2d`, and `X_umap_3d`.

### Common edge cases

- `X_umap_3d` exists but has shape `(n_cells, 2)` → construction fails.

### Fix pattern (recommended)

Store the result under the exact key that declares its dimension:

```python
adata.obsm["X_umap_2d"] = umap_coordinates_2d
```

## Loading and memory behavior (critical for large datasets)

### `.h5ad` (always read-only-backed)

When you pass a `.h5ad` path:
- Cellucid opens it with `anndata.read_h5ad(..., backed="r")`.
- This is usually the best option for large datasets if you don’t want to export.

### `.zarr` (eager)

Cellucid calls `anndata.read_zarr`, which materializes the AnnData arrays in
memory. Use a `.h5ad` path when read-only-backed access is required.

### In-memory `AnnData`

In-memory datasets are convenient but can be risky for large matrices:
- they can trigger expensive conversions (sparse format changes),
- and “accidentally load everything” becomes more likely.

Rule of thumb: if your dataset is big enough that you worry about RAM, view from a file.

## Gene IDs (`gene_id_column`) and duplicates

Cellucid needs a mapping from “what you type in the gene search box” → a column index in the matrix.

By default, `gene_id_column=None` uses `var.index`. Any non-blank string names
that exact `var` column, so `gene_id_column="index"` selects a column literally
named `"index"`.

If your gene IDs are in a column (e.g. `"gene_symbols"`), pass:

```python
viewer = show_anndata(
    "data.h5ad",
    gene_id_column="gene_symbols",
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

```{warning}
AnnData construction rejects duplicate gene IDs and gene IDs whose portable
filename components collide, including case-insensitive collisions. Choose a
unique, portable identifier in `var.index` or the selected
`gene_id_column`.
```

## Vector fields (velocity/drift overlays)

Cellucid can render per-cell displacement vectors as an animated overlay.

### Naming convention (UMAP basis)

Required dimension-suffixed keys in `adata.obsm`:

- `velocity_umap_2d` (shape `(n_cells, 2)`)
- `velocity_umap_3d` (shape `(n_cells, 3)`)
- `T_fwd_umap_2d`, `T_bwd_umap_2d` (CellRank-style drift)

An unsuffixed key such as `velocity_umap` is not discovered as a vector field.

Rules:
- vectors must match `n_cells` and the current dimension (2D vs 3D)
- vectors are scaled into the same normalized coordinate space as the embedding points
- when the AnnData object declares more than one field ID, pass the exact
  `vector_field_default` to `show_anndata(...)`, `serve_anndata(...)`, or
  `AnnDataViewer(...)`; construction raises `ValueError` if it is omitted or
  does not name an available field

For example, an object containing both `velocity_umap_2d` and
`T_fwd_umap_2d` needs:

```python
viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
    vector_field_default="velocity",
)
```

```{figure} ../../../_static/screenshots/vector_field_velocity/overlay-controls.png
:alt: A synthetic 2D dataset with the Vector Field Overlay enabled and animated particle-flow controls visible.
:width: 100%

A current-format 2D vector field rendered as particle flow with its field, density, speed, trail, size, opacity, palette, and LOD controls visible.
```

## Connectivity (KNN graph)

If `adata.obsp["connectivities"]` exists, Cellucid can expose connectivity edges.

Notes:
- the matrix must already be square, finite, non-negative, exactly symmetric
  in topology and weight, and zero on the diagonal;
- sparse inputs must not store explicit zeros or duplicate coordinates;
- Cellucid preserves exact Float64 weights and rejects asymmetric, directed,
  negative, non-finite, or self-edge-containing graphs rather than
  transforming them;
- a valid zero-edge matrix is still an explicitly present graph,
- large graphs can be expensive to render in the browser.

## Troubleshooting (AnnData mode)

### “No valid UMAP embeddings found in adata.obsm”

**Likely causes**
- you didn’t compute UMAP yet
- keys are present but have wrong shape
- `n_cells` mismatch between `adata.obsm[...]` and `adata.n_obs`

**How to confirm**
- print `adata.obsm.keys()` and the shapes of candidate arrays

**Fix**
- compute UMAP (Scanpy) and store explicit keys (`X_umap_2d` / `X_umap_3d`)

### “Gene ‘X’ not found in var”

**Likely causes**
- you’re searching by gene symbol but `var.index` contains Ensembl IDs (or vice versa)
- wrong `gene_id_column`

**Fix**
- set `gene_id_column=...` to match your gene IDs

For all other issues: {doc}`15_troubleshooting_viewing`.

## Next steps

- Export mode (best performance): {doc}`07_exported_directory_mode_show_and_serve`
- Server details + endpoints: {doc}`09_server_mode_advanced`
