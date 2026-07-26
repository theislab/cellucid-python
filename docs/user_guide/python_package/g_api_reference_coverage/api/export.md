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

X_umap_2d = adata.obsm["X_umap_2d"]

prepare(
    latent_space=adata.obsm["X_pca"],
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp.get("connectivities"),

    X_umap_2d=X_umap_2d,

    out_dir="./my_export",
    dataset_name="My study",
    dataset_id="my-study-v1",
    obs_categorical_dtype="uint16",
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
- Naming convention: `<field>_umap_<dim>d`.

See {doc}`vector_fields` for helper functions and naming conventions.

---

## Interface reference

```{figure} ../../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: Cellucid web app with the sidebar open and a single-cell embedding colored by cell type.
:width: 100%

A loaded dataset in Cellucid: the sidebar controls the active view while the categorical legend maps directly to the colored points.
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

---

## API reference

```{eval-rst}
.. autofunction:: prepare
```

---

## Edge cases (do not skip)

### Missing embeddings
- If none of `X_umap_1d`, `X_umap_2d`, `X_umap_3d` are provided, export fails.
- Only `X_umap_1d`, `X_umap_2d`, and `X_umap_3d` are accepted keyword names.

### Shape mismatches
- All inputs must agree on `n_cells` (rows).
- `gene_expression` must be `(n_cells, n_genes)` and `var` must describe those genes.

### Required inputs (common surprises)
- `latent_space` is required and must be supplied explicitly.
- `obs` is required and must be aligned to the same cell order as embeddings.
- If you provide `gene_expression`, you must also provide `var` (to name/describe genes).

### Existing export directories
- A non-empty target directory raises `FileExistsError` when `force=False`.
- `force=True` publishes one complete replacement generation.

### Vector field validation
- `vector_fields` must be a dict of arrays (or sparse matrices).
- Each vector field must be 1D or 2D and have 1/2/3 components (for 1D/2D/3D overlays).
- Keys must match `<field>_umap_1d`, `<field>_umap_2d`, or
  `<field>_umap_3d`, and the matching embedding dimension must exist.

### NaN/Inf and constant-value fields
- Scientific arrays must contain real finite values representable as
  `float32`; invalid values reject before publication.
- Constant continuous fields are valid without quantization. Quantization
  rejects a constant field because the current compact contract requires
  `minValue < maxValue`.

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

---

### Symptom: “var is required if gene_expression is provided”
Fix:
- Pass `var=adata.var` whenever you pass `gene_expression=adata.X`.

---

### Symptom: an unsupported embedding keyword is rejected
Fix:
- Pass only `X_umap_1d`, `X_umap_2d`, or `X_umap_3d`.

---

### Symptom: “non-empty export directory”
What’s happening:
- `prepare(...)` refuses to mix a new generation with existing artifacts.

Fix:
- If you want one complete replacement, call `prepare(..., force=True)`.
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
- Export fewer genes (`gene_identifiers=`) or omit gene expression for a
  metadata-only export.
- Use quantization (`var_quantization=8`) and moderate compression (`compression=6`).

---

## See also

- {doc}`server` for serving exports (local + SSH tunnel)
- {doc}`jupyter` for embedding exports in notebooks
- {doc}`../../c_data_preparation_api/09_output_format_specification_exports_directory` for a deeper format spec (user guide section)
