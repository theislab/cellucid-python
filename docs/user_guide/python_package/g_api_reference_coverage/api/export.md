# Export / Data Preparation (`prepare`)

```{eval-rst}
.. currentmodule:: cellucid
```

This page documents {func}`~cellucid.prepare`, which writes an **exported dataset directory** that can be:
- opened via the web app (file picker),
- served with `cellucid serve ./export`,
- embedded in notebooks with {func}`~cellucid.show`.

If you want the long-form, user-guide style walkthroughs, see:
- {doc}`../../c_data_preparation_api/index` (data prep API, shapes, edge cases)

---

## Audience + prerequisites

**Audience**
- Wet lab / beginner: use the copy/paste example and the troubleshooting section.
- Computational: read the parameter and performance sections before exporting large datasets.
- Developer: read the output format + determinism notes if you need stable exports for papers/CI.

**Prerequisites**
- `pip install cellucid`
- Typically: `numpy`, `pandas`, `scipy`
- If sourcing from AnnData: you’ll usually have `adata.obsm`, `adata.obs`, `adata.var`, `adata.X`

---

## Fast path (copy/paste)

```python
from cellucid import prepare

X_umap = adata.obsm["X_umap"]  # shape: (n_cells, 2) or (n_cells, 3)

prepare(
    latent_space=adata.obsm.get("X_pca", X_umap),
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp.get("connectivities"),

    # Provide at least one embedding (1D/2D/3D). 4D is reserved (not implemented yet).
    X_umap_2d=adata.obsm.get("X_umap_2d", X_umap if X_umap.shape[1] == 2 else None),
    X_umap_3d=adata.obsm.get("X_umap_3d", X_umap if X_umap.shape[1] == 3 else None),

    out_dir="./my_export",
    compression=6,
    var_quantization=8,
    obs_continuous_quantization=8,
)
```

---

## Practical path (what to decide before you export)

### 1) Do you want reproducibility or convenience?

- Convenience: use {func}`~cellucid.show_anndata` / {func}`~cellucid.serve_anndata` (no export).
- Reproducibility + shareability: use {func}`~cellucid.prepare` once, then reuse the folder.

### 2) Choose compression and quantization

These trade size vs speed vs fidelity:
- `compression=6` is a good default for gzip.
- `var_quantization=8` is usually enough for coloring by gene expression.
- `obs_continuous_quantization=8` is usually enough for QC metrics and scores.

If you need *exact* values preserved:
- set quantization options to `None` (writes float32).

### 3) Decide which obs/genes you will ship

- `obs_keys=None` exports all `obs` columns. For very wide `obs`, consider selecting a subset.
- `gene_identifiers=None` exports all genes. For huge `n_genes`, consider a curated list.

### 4) Vector fields (velocity / drift overlays)

Vector fields are optional, but powerful:
- They are **per-cell displacement vectors in embedding space** (not “3D arrows in physical space”).
- Naming convention:
  - Explicit: `<field>_umap_<dim>d` (recommended)
  - Implicit: `<field>_umap` (only if explicit keys aren’t provided)

See {doc}`vector_fields` for helper functions and naming conventions.

---

## Screenshot placeholder (optional)

<!-- SCREENSHOT PLACEHOLDER
ID: export-loaded-in-web-app
Where it appears: Export / Data Preparation → Fast path
Capture:
  - UI location: the Cellucid web app dataset loader after successfully loading an exported folder
  - State prerequisites: you have a working export folder created by `prepare(...)`
  - Action to reach state:
    - start a server: `cellucid serve ./my_export`
    - open the viewer URL in your browser
Crop:
  - Include: dataset name + point count + one visible embedding view
  - Exclude: private dataset path, browser bookmarks, user accounts
Redact:
  - Remove: sample IDs, lab-internal dataset names (if private)
Annotations:
  - Callouts:
    - (1) dataset name / identity indicator
    - (2) points visible (success state)
    - (3) field selector / legend (to confirm metadata loaded)
Alt text:
  - Cellucid web viewer showing an exported dataset loaded successfully, with points visible in an embedding.
Caption:
  - After exporting with `prepare(...)`, serving the folder and opening the viewer URL loads the dataset with metadata and embeddings ready for exploration.
-->
```{figure} ../../../../_static/screenshots/placeholder-screenshot.svg
:alt: Cellucid web viewer showing an exported dataset loaded successfully, with points visible in an embedding.
:width: 100%

After exporting with `prepare(...)`, serving the folder and opening the viewer URL loads the dataset with metadata and embeddings ready for exploration.
```

---

## Output directory layout (high-level)

An exported dataset directory typically contains:

```
my_export/
├── dataset_identity.json
├── obs_manifest.json
├── var_manifest.json                 # optional (gene expression)
├── connectivity_manifest.json        # optional (KNN edges)
├── points_1d.bin.gz                  # optional
├── points_2d.bin.gz                  # optional
├── points_3d.bin.gz                  # optional
├── obs/                              # obs field binaries
├── var/                              # gene expression binaries
├── connectivity/                     # KNN edge binaries
└── vectors/                          # optional (vector field binaries)
```

Notes:
- You must provide at least one of `points_1d/2d/3d`.
- 4D (`points_4d`) is reserved for future development.

---

## API reference

```{eval-rst}
.. autofunction:: prepare
```

---

## Edge cases (do not skip)

### Missing embeddings
- If none of `X_umap_1d`, `X_umap_2d`, `X_umap_3d` are provided, export fails.
- If you pass `X_umap_4d`, export raises `NotImplementedError` (reserved for future viewer support).

### Shape mismatches
- All inputs must agree on `n_cells` (rows).
- `gene_expression` must be `(n_cells, n_genes)` and `var` must describe those genes.

### Required inputs (common surprises)
- `latent_space` is required (used for outlier quantiles); if you don’t have PCA, you can often reuse an embedding as a fallback.
- `obs` is required and must be aligned to the same cell order as embeddings.
- If you provide `gene_expression`, you must also provide `var` (to name/describe genes).

### File overwrites vs “skip existing”
- By default, `prepare(...)` skips writing files that already exist (to avoid accidental overwrites).
- Use `force=True` if you intentionally want to regenerate files.

### Vector field validation
- `vector_fields` must be a dict of arrays (or sparse matrices).
- Each vector field must be 1D or 2D and have 1/2/3 components (for 1D/2D/3D overlays).
- Keys should follow naming conventions so the viewer can find the right dimension.

### NaN/Inf and constant-value fields
- Continuous quantization reserves a missing-value marker; NaN/Inf will be mapped to “missing”.
- Constant-value fields are allowed, but are not informative for coloring/filtering.

### Very large datasets
- Export size scales with:
  - number of cells × number of exported embeddings
  - number of obs fields
  - number of genes included
- For huge datasets, consider:
  - fewer genes,
  - quantization,
  - compression,
  - serving from a fast filesystem.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “At least one dimensional embedding must be provided”
Fix:
- Provide one of `X_umap_1d`, `X_umap_2d`, `X_umap_3d`.

### Symptom: “All embeddings must have the same number of cells”
Fix:
- Ensure every embedding array has exactly the same number of rows (cell order must match).

### Symptom: “obs has N rows, but embeddings have M cells”
Fix:
- Ensure `obs` row order corresponds to the embedding row order (and gene expression if present).

### Symptom: “Export folder is huge”
Fix:
- Enable quantization (`var_quantization`, `obs_continuous_quantization`).
- Enable gzip (`compression`).
- Export fewer genes (`gene_identifiers=`) and/or fewer obs columns (`obs_keys=`).

---

### Symptom: “latent_space is required for outlier quantile calculation”
Fix:
- Provide `latent_space=...` with shape `(n_cells, n_latent_dims)`.
- Common choices:
  - `adata.obsm["X_pca"]`
  - `adata.obsm["X_scvi"]`
  - as a fallback: reuse `X_umap_2d`/`X_umap_3d` (less ideal, but workable)

---

### Symptom: “var is required if gene_expression is provided”
Fix:
- Pass `var=adata.var` whenever you pass `gene_expression=adata.X`.

---

### Symptom: “4D visualization is not yet implemented”
Fix:
- Do not pass `X_umap_4d`.
- Use `X_umap_1d`, `X_umap_2d`, or `X_umap_3d`.

---

### Symptom: “⚠ Skipping … already exists”
What’s happening:
- `prepare(...)` avoids overwriting existing files by default.

Fix:
- If you want to overwrite, call `prepare(..., force=True)`.
- If you want a clean export, write to a new `out_dir`.

---

### Symptom: “Vector field '…' must have 1/2/3 components”
Fix:
- Ensure each vector field array is shaped `(n_cells, dim)` with `dim ∈ {1,2,3}` (or 1D arrays for 1D).
- Ensure the vector field’s `n_cells` matches the embedding `n_cells`.

---

### Symptom: “Export is extremely slow”
Likely causes:
- You are exporting many genes (large `n_cells × n_genes`).
- Compression level is very high.

Fix:
- Export fewer genes (`gene_identifiers=`) or skip gene expression entirely for metadata-only exports.
- Use quantization (`var_quantization=8`) and moderate compression (`compression=6`).

---

## See also

- {doc}`server` for serving exports (local + SSH tunnel)
- {doc}`jupyter` for embedding exports in notebooks
- {doc}`../../c_data_preparation_api/09_output_format_specification_exports_directory` for a deeper format spec (user guide section)
