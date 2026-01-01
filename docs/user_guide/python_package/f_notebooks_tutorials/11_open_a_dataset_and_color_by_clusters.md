# Open a dataset and color by clusters

This tutorial gets you to the first “real success state”:
- a dataset loads in Cellucid
- you can color points by a cluster/label column
- you understand what to do if that column is missing or looks wrong

It assumes **zero knowledge** of Cellucid, and only minimal Python.

---

## At a glance

**Audience**
- Wet lab / beginner: follow the “Minimal cells” and the UI steps.
- Computational users: scan the “Data requirements” and “Edge cases”.

**Time**
- If you already have UMAP + clusters: ~5–10 minutes
- If you need to compute them: ~20–40 minutes (dataset-dependent)

**Prerequisites**
- `pip install cellucid`
- A dataset in one of these forms:
  - an in-memory `AnnData`
  - a `.h5ad` file
  - a `.zarr` store directory

---

## Minimal cells (copy/paste)

If you already have `adata`:

```python
from cellucid import show_anndata

viewer = show_anndata(adata, height=650)
viewer
```

If you have a file:

```python
from cellucid import show_anndata

viewer = show_anndata("my_dataset.h5ad", height=650)  # or "my_dataset.zarr"
viewer
```

Then inside the viewer UI:
- find the **Field** (color-by) selector
- choose your cluster column (often called `leiden`, `louvain`, `cell_type`, `seurat_clusters`, etc.)

If you don’t have a cluster column or UMAP, continue with the step-by-step sections below.

---

## Step 0 — Get an `AnnData` (three safe options)

### Option A: you already have `adata` in memory

You can skip this step.

### Option B: load a `.h5ad` from disk

```python
import anndata as ad

adata = ad.read_h5ad("my_dataset.h5ad")
adata
```

### Option C: use a small demo dataset (downloads data)

```{warning}
Some demo datasets download files from the internet. If your environment has restricted network access, use Option A or B instead.
```

Example (Scanpy PBMC3k):

```python
import scanpy as sc

adata = sc.datasets.pbmc3k()
adata
```

---

## Step 1 — Sanity-check the dataset *before* opening the viewer

Cellucid assumes that:
- each row in `adata` is one cell (`adata.n_obs`)
- each column in `adata` is one feature/gene (`adata.n_vars`)
- per-cell metadata lives in `adata.obs`
- embeddings live in `adata.obsm`

Run:

```python
print("cells:", adata.n_obs)
print("genes:", adata.n_vars)
print("obs columns:", list(adata.obs.columns)[:10], "...")
print("obsm keys:", list(adata.obsm.keys()))
```

### What you need for this tutorial

You need:
- **one embedding** (2D or 3D recommended) in `adata.obsm`
- **one “cluster labels” column** in `adata.obs` (categorical/string is fine)

---

## Step 2 — Make sure you have an embedding (UMAP is the common case)

### If you already have a UMAP embedding

Look for a key like:
- `X_umap` (Scanpy convention; usually 2D)
- `X_umap_3d` (if you explicitly stored 3D)
- sometimes project-specific names like `umap`, `X_umap_harmony`, etc.

If you see `X_umap` in `adata.obsm`, you can proceed.

### If you do NOT have UMAP (compute it with Scanpy)

```{note}
This is the “standard recipe” for many single-cell workflows. If you already computed UMAP elsewhere (Seurat, scVI, etc.), you do not need to recompute—just ensure the embedding is present in `adata.obsm`.
```

Minimal Scanpy pipeline:

```python
import scanpy as sc

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var["highly_variable"]].copy()

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

Now you should have:

```python
list(adata.obsm.keys())
```

and `X_umap` should appear.

---

## Step 3 — Make sure you have a cluster/label column

### If you already have a cluster column

Common names:
- `leiden`, `louvain`
- `cell_type`, `celltype`
- `seurat_clusters`
- `cluster`, `clusters`

You can list categorical-looking fields quickly:

```python
import pandas as pd

categorical_like = [
    c for c in adata.obs.columns
    if pd.api.types.is_categorical_dtype(adata.obs[c])
    or pd.api.types.is_object_dtype(adata.obs[c])
    or pd.api.types.is_string_dtype(adata.obs[c])
]
categorical_like[:25]
```

### If you do NOT have a cluster column (compute Leiden)

```{note}
Clustering is not “required” for Cellucid, but it is one of the easiest ways to start exploring. If you already have labels from another tool, prefer those.
```

```python
import scanpy as sc

# Requires that neighbors graph exists (sc.pp.neighbors)
sc.tl.leiden(adata, key_added="leiden")
adata.obs["leiden"].value_counts().head()
```

---

## Step 4 — Open the viewer from Python

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    height=650,
    dataset_name="My dataset (tutorial)",
)
viewer
```

```{tip}
If you want Cellucid to use a specific gene identifier column (e.g. gene symbols stored in `adata.var["gene_symbols"]`), pass:

`show_anndata(adata, gene_id_column="gene_symbols")`
```

---

## Step 5 — Color by clusters (UI steps)

In the viewer:

1) Open the **Field selector** (color-by)
2) Find your cluster column (e.g. `leiden`)
3) Select it

### What success looks like

- points change from neutral gray to a category palette
- a legend appears showing cluster names and colors
- the number of categories matches what you expect (roughly)

### Optional: set color-by from Python (only if you want automation)

This is useful in teaching notebooks where you want the viewer to land in a known state.

```python
@viewer.on_ready
def _set_default_color(_event):
    viewer.set_color_by("leiden")
```

```{note}
If the field name is wrong or missing, nothing catastrophic happens; the command will simply not do what you expect. Use the UI field selector to confirm the real field name.
```

---

## Screenshot placeholders (optional but helpful)

<!-- SCREENSHOT PLACEHOLDER
ID: python-notebooks-color-by-clusters-field-selector
Suggested filename: web_app/00_field-selector-color-by-clusters.png
Where it appears: User Guide → Python Package → Notebooks/Tutorials → 11_open_a_dataset_and_color_by_clusters.md
Capture:
  - UI location: field selector open (obs fields visible)
  - State prerequisites: dataset loaded; at least one categorical obs field present (e.g. leiden)
  - Action to reach state: open field selector and hover the cluster field to show tooltip (if any)
Crop:
  - Include: field selector UI + search box + the selected field highlighted
  - Include: enough of the canvas to show categorical coloring applied
Redact:
  - Remove: any private dataset name/sample IDs visible in the UI
Alt text:
  - Field selector open with a cluster field selected, and points colored categorically.
Caption:
  - Use the field selector to color points by a cluster/label column (categorical obs field).
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for coloring by a cluster field.
:width: 100%

Field selector open with a cluster field selected; points are colored categorically and a legend is visible.
```

---

## Edge cases (common “gotchas”)

### “I don’t see my cluster column in the field list”

Possible causes:
- the column is not in `adata.obs` (it might be in a separate table)
- the column contains non-scalar objects (lists/dicts), which don’t convert cleanly to categories
- the column name is different than you think (case/underscores)

Confirm in Python:
```python
print("leiden" in adata.obs.columns)
print(adata.obs.dtypes.head(20))
```

### “My clusters are numeric (0/1/2/...) and the viewer treats them as continuous”

If a column is numeric, Cellucid treats it as a **continuous** field.

Fix: convert to categorical in Python:
```python
adata.obs["leiden"] = adata.obs["leiden"].astype("category")
```

### “I have too many categories”

If your field has hundreds/thousands of categories (e.g. raw barcodes, per-cell IDs), categorical coloring becomes unusable.

Best practice:
- don’t use unique-per-cell identifiers as “cluster” fields
- aggregate or filter categories before visualizing

### “Some cells are missing / invisible”

This is usually filters, not coloring.
See web-app docs:
- {doc}`../../web_app/e_filtering/index`

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “The viewer iframe is blank or says it can’t load”

Likely causes (most common first):
1) the Python server failed to start (port conflict, exception)
2) web UI cache is missing/corrupt (first run offline, blocked network)
3) remote kernel without proxy/tunnel

How to confirm:
```python
viewer.debug_connection()
```

Fix:
- pick a different port (advanced: construct viewer with `port=...`)
- clear the web cache and retry (see {doc}`../../web_app/b_data_loading/05_jupyter_tutorial`)
- use Jupyter Server Proxy or SSH tunneling (see {doc}`22_large_dataset_server_mode_and_lazy_gene_expression`)

### Symptom: “Color-by works, but categories look wrong”

Likely causes:
- labels are not what you think (e.g. mixed dtypes, trailing spaces)
- you’re using a different column (similar name)
- there are missing values (`NaN`) that show up as “missing”

How to confirm:
```python
adata.obs["leiden"].astype(str).value_counts().head(20)
```

Fix:
- clean labels (strip whitespace, unify casing)
- convert to categorical explicitly

---

## Next steps

- {doc}`12_find_marker_genes_and_export_a_figure` (go from clusters → marker genes → figure)
- {doc}`13_compare_two_groups_and_interpret_results` (compare groups safely)
- {doc}`21_prepare_exports_with_quantization_and_compression` (reproducible, shareable exports)
