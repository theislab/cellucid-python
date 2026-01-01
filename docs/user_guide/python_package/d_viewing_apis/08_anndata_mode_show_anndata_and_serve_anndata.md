# AnnData mode (`show_anndata()` / `serve_anndata()`)

This page explains how Cellucid can visualize **AnnData directly**, without running `prepare()` first.

Use AnnData mode when:
- you’re iterating inside an analysis notebook,
- you want a quick “does this dataset look sane?” view,
- or you have a large `.h5ad`/`.zarr` and want lazy loading without exporting.

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

viewer = show_anndata("data.h5ad", height=600)  # backed mode by default
# viewer = show_anndata("data.zarr", height=600)
# viewer = show_anndata(adata, height=600)
```

### Terminal

```bash
cellucid serve /path/to/data.h5ad
```

### Python (script)

```python
from cellucid import serve_anndata

serve_anndata("data.h5ad")  # blocks
```

## Minimum AnnData requirements (the “must haves”)

### 1) UMAP embedding in `adata.obsm`

Cellucid currently expects UMAP coordinates under one of these keys:

- explicit (preferred):
  - `X_umap_1d` (shape `(n_cells, 1)`)
  - `X_umap_2d` (shape `(n_cells, 2)`)
  - `X_umap_3d` (shape `(n_cells, 3)`)
- or a generic fallback:
  - `X_umap` (shape `(n_cells, 1|2|3)`) **only if no explicit key for that dimension exists**

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

Cellucid looks for embeddings in this order:

1) `X_umap_1d`, `X_umap_2d`, `X_umap_3d` (explicit)
2) `X_umap` (fallback) if:
   - it has shape `(n_cells, 1|2|3)`, and
   - there is no explicit key for that same dimension (explicit wins)

### Common edge cases

- `X_umap` is present but has the wrong shape → ignored.
- `X_umap_3d` exists but has shape `(n_cells, 2)` → ignored.
- both `X_umap_2d` and `X_umap` (2D) exist → `X_umap_2d` wins; `X_umap` is ignored.

### Fix pattern (recommended)

Always store explicit keys so there’s no ambiguity:

```python
# After computing UMAP:
X_umap = adata.obsm["X_umap"]
if X_umap.shape[1] == 1:
    adata.obsm["X_umap_1d"] = X_umap
elif X_umap.shape[1] == 2:
    adata.obsm["X_umap_2d"] = X_umap
elif X_umap.shape[1] == 3:
    adata.obsm["X_umap_3d"] = X_umap
```

## Lazy loading and memory behavior (critical for large datasets)

### `.h5ad` (default: backed mode)

When you pass a `.h5ad` path:
- Cellucid opens it in **backed** mode by default (lazy, memory-efficient).
- This is usually the best option for large datasets if you don’t want to export.

Disable backed mode only if you know you want full in-memory behavior:

```bash
cellucid serve data.h5ad --no-backed
```

### `.zarr` (naturally chunked)

Zarr is directory-based and chunked; access is naturally “lazy-ish”.

### In-memory `AnnData`

In-memory datasets are convenient but can be risky for large matrices:
- they can trigger expensive conversions (sparse format changes),
- and “accidentally load everything” becomes more likely.

Rule of thumb: if your dataset is big enough that you worry about RAM, view from a file.

## Gene IDs (`gene_id_column`) and duplicates

Cellucid needs a mapping from “what you type in the gene search box” → a column index in the matrix.

By default it uses:
- `var.index` (`gene_id_column="index"`)

If your gene IDs are in a column (e.g. `"gene_symbols"`), pass:

```python
viewer = show_anndata("data.h5ad", gene_id_column="gene_symbols")
```

```{warning}
Duplicate gene IDs are ambiguous.

In AnnData mode, Cellucid builds an internal dict `{gene_id: index}`. If gene IDs repeat, later entries overwrite earlier ones, and the viewer will effectively show one of the duplicates.

Best practice: ensure your chosen gene ID column is unique before viewing/exporting.
```

## Vector fields (velocity/drift overlays)

Cellucid can render per-cell displacement vectors as an animated overlay.

### Naming convention (UMAP basis)

Preferred explicit keys in `adata.obsm`:

- `velocity_umap_2d` (shape `(n_cells, 2)`)
- `velocity_umap_3d` (shape `(n_cells, 3)`)
- `T_fwd_umap_2d`, `T_bwd_umap_2d` (CellRank-style drift)

Implicit keys are also supported:

- `velocity_umap` with shape `(n_cells, 2)` or `(n_cells, 3)`

Rules:
- vectors must match `n_cells` and the current dimension (2D vs 3D)
- vectors are scaled into the same normalized coordinate space as the embedding points

<!-- SCREENSHOT PLACEHOLDER
ID: vector-field-overlay-controls
Suggested filename: vector_field_velocity/python_01_overlay-controls.png
Where it appears: Python Package Guide → Viewing APIs → AnnData mode → Vector fields
Capture:
  - UI location: Cellucid viewer UI (overlay controls panel)
  - State prerequisites: dataset includes at least one vector field; viewer set to matching dimension
  - Action to reach state: load AnnData with `velocity_umap_2d` present; enable overlay
Crop:
  - Include: the overlay toggle + vector field dropdown + dimension indicator (2D/3D)
Alt text:
  - Cellucid overlay panel showing vector field controls.
Caption:
  - Vector fields appear only when the dataset advertises UMAP-aligned displacement vectors for the current dimension.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for vector field overlay controls.
:width: 100%

Vector field (velocity/drift) overlay controls in the Cellucid viewer.
```

## Connectivity (KNN graph)

If `adata.obsp["connectivities"]` exists, Cellucid can expose connectivity edges.

Notes:
- the matrix is symmetrized and binarized for visualization (weights are not preserved),
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
