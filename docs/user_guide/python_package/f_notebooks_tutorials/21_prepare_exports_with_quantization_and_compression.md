# Prepare exports with quantization and compression

This tutorial shows how to create a **reproducible, shareable export folder** for Cellucid.

Exports are the recommended way to use Cellucid when you care about:
- speed (fast loading in browser / fewer crashes)
- reproducibility (stable dataset identity; fewer hidden notebook-side assumptions)
- sharing (hand a folder to a collaborator; host it on a server)

You will learn:
- how to call `cellucid.prepare(...)` from an `AnnData`
- what quantization and gzip compression actually do
- how to choose safe defaults
- how to validate the export folder and debug common failures

---

## At a glance

**Audience**
- Computational users: primary audience (you control preprocessing/export settings).
- Beginners: you can still follow, but start with {doc}`11_open_a_dataset_and_color_by_clusters` first.

**Time**
- Small datasets (<50k cells): ~10–30 minutes
- Large datasets: depends on I/O and compression settings

**Prerequisites**
- `pip install cellucid`
- An `AnnData` with:
  - at least one embedding (recommended: UMAP 2D or 3D)
  - an expression matrix in `adata.X` (sparse is fine)

---

## Why exports exist (one paragraph mental model)

When you call `show_anndata(...)`, the viewer can work, but it is not optimized for repeated use:
- data is served dynamically from an `AnnData` object
- your notebook state becomes part of “the dataset”

When you call `prepare(...)`, Cellucid writes a **dataset directory** with:
- embedding point clouds (`points_2d.bin(.gz)`, `points_3d.bin(.gz)`, ...)
- obs metadata manifest + binary columns
- var (gene) manifest + per-gene expression vectors
- optional connectivities and vector fields
- a `dataset_identity.json` so the dataset can be recognized and sessions can be applied

After that, you view the export folder with:
- `show("./exports/my_dataset")` (in notebooks), or
- `cellucid serve ./exports/my_dataset` (CLI), or
- any static hosting setup (advanced; out of scope here)

---

## Minimal export (copy/paste)

```python
from pathlib import Path
import numpy as np

from cellucid import prepare, show

out_dir = Path("exports") / "my_dataset"

X_umap_2d = np.asarray(adata.obsm["X_umap"])  # adjust if your key differs
connectivities = adata.obsp["connectivities"] if "connectivities" in adata.obsp else None

prepare(
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    var_gene_id_column="index",
    connectivities=connectivities,
    X_umap_2d=X_umap_2d,
    out_dir=out_dir,
    # Practical defaults (tune later):
    var_quantization=8,
    obs_continuous_quantization=8,
    obs_categorical_dtype="auto",
    compression=6,
    dataset_name="My dataset",
    dataset_id="my_dataset",
)

viewer = show(out_dir, height=650)
viewer
```

If any key names differ (`X_umap`, connectivities location, gene IDs), continue with the step-by-step sections below.

---

## Step 1 — Choose a stable export directory and dataset identity

### Export directory (`out_dir`)

Pick a directory that you can:
- keep stable across runs
- share with collaborators (zip it, put it on a server, commit to a data repo)

Examples:
- `exports/pbmc_demo`
- `exports/he_developmental_v1`

### Dataset identity (`dataset_id`)

`dataset_id` is how Cellucid decides “is this the same dataset?”.
It matters for:
- session bundles (restore highlights/filters/etc)
- collaboration/sharing

Best practice:
- set `dataset_id` explicitly
- keep it stable if the underlying data is the same

---

## Step 2 — Extract the required inputs from AnnData

`prepare(...)` does not take an `AnnData` object directly; it expects the pieces.

### Required: at least one embedding

Most datasets have 2D UMAP in:
- `adata.obsm["X_umap"]` with shape `(n_cells, 2)`

If your UMAP is 3D, you might have:
- `adata.obsm["X_umap_3d"]` with shape `(n_cells, 3)`

Extract safely:

```python
import numpy as np

X_umap_2d = np.asarray(adata.obsm["X_umap"])  # (n_cells, 2)
```

```{note}
Embeddings are normalized by `prepare(...)` into `[-1, 1]` independently per dimension. This preserves geometry but discards absolute scale (which is usually fine for visualization).
```

### Required (for gene coloring): expression matrix + var table

```python
obs = adata.obs
var = adata.var
X = adata.X
```

If you want the viewer’s gene search to use gene symbols stored in a column:

```python
var_gene_id_column = "gene_symbols"  # example; use your actual column name
```

Otherwise:

```python
var_gene_id_column = "index"  # use var.index (often equals adata.var_names)
```

### Optional: connectivity graph

If you have a KNN graph:

```python
connectivities = adata.obsp["connectivities"] if "connectivities" in adata.obsp else None
```

### Optional: latent space for outlier quantiles/centroids

If you have a latent space (PCA, scVI, etc.), you can pass it:

```python
latent_space = adata.obsm.get("X_pca")  # example
```

This is used for:
- category centroid computation
- per-cell outlier quantiles (useful for some UI behaviors)

---

## Step 3 — Choose quantization + compression settings (and what they mean)

### Quantization (8-bit vs 16-bit vs none)

Quantization reduces file size by storing continuous values as integers with a scale:

- `var_quantization=8` means per-gene expression is stored as `uint8`
  - ~4× smaller than float32
  - usually visually indistinguishable for colormaps
- `var_quantization=16` is larger but preserves more precision
- `None` stores full float32 (largest)

Important details:
- reserved “missing” codes are used for NaN/Inf:
  - 8-bit reserves `255`
  - 16-bit reserves `65535`
- constant-valued fields are handled (scale is adjusted to avoid divide-by-zero)

### Compression (gzip)

Compression reduces disk size and (often) network transfer size:
- `compression=6` is a common balance (smaller files, still reasonable CPU time)
- `compression=1` is faster, slightly larger files
- `compression=9` is slowest, smallest files

Compression adds a `.gz` suffix to the written files.

### Recommended defaults

For most users:
- `var_quantization=8`
- `obs_continuous_quantization=8`
- `obs_categorical_dtype="auto"`
- `compression=6`

If you are doing extremely subtle continuous gradients and care about smoothness:
- try `var_quantization=16`

---

## Step 4 — Run `prepare(...)`

```python
from cellucid import prepare

prepare(
    latent_space=latent_space,
    obs=obs,
    var=var,
    gene_expression=X,
    var_gene_id_column=var_gene_id_column,
    connectivities=connectivities,
    X_umap_2d=X_umap_2d,
    out_dir=out_dir,
    force=False,
    var_quantization=8,
    obs_continuous_quantization=8,
    obs_categorical_dtype="auto",
    compression=6,
    dataset_name="My dataset",
    dataset_description="Exported from AnnData for Cellucid viewing",
    dataset_id="my_dataset",
)
```

---

## Step 5 — Validate the export folder (before sharing it)

### Check key files exist

In a shell:

```bash
ls -la exports/my_dataset
```

You should see (at minimum):
- `dataset_identity.json`
- `obs_manifest.json`
- at least one `points_*d.bin` (or `.bin.gz`)

### Open in a notebook

```python
from cellucid import show

viewer = show("exports/my_dataset", height=650)
viewer
```

### Or open via CLI

```bash
cellucid serve exports/my_dataset
```

---

## Screenshot placeholder (optional: show export loaded successfully)

<!-- SCREENSHOT PLACEHOLDER
ID: python-notebooks-export-opened-success
Suggested filename: data_loading/01_export-folder-opened-success.png
Where it appears: User Guide → Python Package → Notebooks/Tutorials → 21_prepare_exports_with_quantization_and_compression.md
Capture:
  - UI location: browser (or notebook) showing a successfully loaded export folder
  - State prerequisites: exported dataset opened; points visible; a field selected
  - Action to reach state: run `show("./exports/<id>")` or `cellucid serve ./exports/<id>` and open the URL
Crop:
  - Include: dataset name + point count area (if shown) and a visible embedding
  - Exclude: personal paths/usernames
Alt text:
  - Cellucid viewer showing a loaded exported dataset with points visible.
Caption:
  - After `prepare(...)`, you can open the export folder directly; this is the fastest and most reproducible viewing mode.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for opening an exported dataset successfully.
:width: 100%

An exported dataset opened successfully (points visible and interactive).
```

---

## Edge cases (and how to handle them)

### Shape mismatches (the most common export failure)

All per-cell arrays must have the same first dimension (`n_cells`):
- `obs` rows
- embeddings (`X_umap_*d`)
- `gene_expression` rows
- `connectivities` shape `(n_cells, n_cells)`

Confirm:

```python
print("n_cells:", adata.n_obs)
print("obs rows:", adata.obs.shape[0])
print("X shape:", adata.X.shape)
print("umap shape:", X_umap_2d.shape)
```

### Too many categories for uint8

If you force `obs_categorical_dtype="uint8"`, you can represent at most **254** categories (plus one missing sentinel).

If you have more categories, use:
- `obs_categorical_dtype="auto"` (safe)
- or `obs_categorical_dtype="uint16"`

### Gene IDs not searchable / confusing

If the viewer’s gene search feels wrong:
- you probably picked the wrong `var_gene_id_column`
- or you have duplicates/aliases

Confirm:
```python
var_gene_id_column
var[var_gene_id_column].head()
```

### “prepare skipped files” surprises

By default, `force=False` means:
- if a file already exists, it is **not overwritten**

This is safe, but can surprise you during iteration.

If you are intentionally re-exporting:
- set `force=True`
- or delete the export directory and re-run

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Export folder is huge”

Likely causes:
- `compression=None` (no gzip)
- `var_quantization=None` (float32 expression)
- exporting all genes for a very large matrix

Fix:
- enable gzip (`compression=6`)
- quantize (`var_quantization=8`)
- consider exporting a gene subset (`gene_identifiers=...`) for interactive workflows

### Symptom: “prepare(...) is extremely slow”

Likely causes:
- high gzip level (`compression=9`)
- extremely large sparse matrix being densified somewhere (avoid)
- disk is slow or network-mounted

Fix:
- use `compression=1` or `compression=6`
- export on fast local disk when possible
- export fewer genes first, validate workflow, then scale up

### Symptom: “The viewer loads the export but gene coloring is missing”

Likely causes:
- `gene_expression` was not provided (or was empty)
- `var` was not provided (required to define gene IDs)
- `var_gene_id_column` points to a missing column

Fix:
- verify inputs before running `prepare(...)`
- open `var_manifest.json` to confirm genes are present

---

## Next steps

- {doc}`04_real_world_dataset_recipes_gallery` (full dataset export notebooks you can adapt)
- {doc}`22_large_dataset_server_mode_and_lazy_gene_expression` (scaling patterns)
- {doc}`32_session_persistence_and_restoring_analysis_artifacts` (advanced reproducibility via sessions)

If you want the full export format spec:
- {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`
