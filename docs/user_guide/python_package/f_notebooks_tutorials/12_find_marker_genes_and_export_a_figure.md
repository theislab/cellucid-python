# Find marker genes and export a figure

This tutorial connects three things that often feel disconnected to new users:

1) **Clusters/labels** (a categorical obs field)
2) **Marker genes** (a ranked list you compute in Python)
3) **A publication-grade figure** (exported from the Cellucid viewer)

It is intentionally explicit about:
- gene identifier pitfalls (symbols vs Ensembl IDs, duplicates)
- what “marker gene” actually means (and what it does *not* prove)
- how to export figures without guessing

---

## At a glance

**Audience**
- Wet lab / beginner: follow the “Fast path” and the UI instructions.
- Computational users: focus on the gene-ID and edge-case sections.

**Time**
- ~20–45 minutes (depends on dataset size and whether you already have clusters)

**Prerequisites**
- You completed {doc}`11_open_a_dataset_and_color_by_clusters` (or you already have a cluster column)
- Optional but recommended: `pip install scanpy`

---

## Fast path (minimal workflow)

1) Compute marker genes for your clusters in Python.
2) In the viewer, search one of those genes and color by expression.
3) Export a figure (PNG/SVG).

---

## Step 0 — Confirm you have a cluster column

Pick a cluster/label field in `adata.obs`, for example:
- `leiden` (Scanpy)
- `seurat_clusters` (Seurat exports)
- `cell_type` (manual annotation)

```python
cluster_key = "leiden"  # change if needed
cluster_key in adata.obs.columns
```

If you don’t have a cluster column yet, go back to:
- {doc}`11_open_a_dataset_and_color_by_clusters`

---

## Step 1 — Compute marker genes (Scanpy example)

```{note}
This tutorial uses Scanpy because it is widely used and easy to reproduce, but “marker genes” is a general concept: any differential expression approach can produce a ranked list.
```

Minimal Scanpy workflow:

```python
import scanpy as sc

# If your dataset is already normalized/log-transformed, do NOT repeat preprocessing blindly.
# The minimal safe assumption for rank_genes_groups is that adata.X represents comparable values across cells.

sc.tl.rank_genes_groups(
    adata,
    groupby=cluster_key,
    method="wilcoxon",
)
```

Inspect results:

```python
sc.get.rank_genes_groups_df(adata, group=adata.obs[cluster_key].cat.categories[0]).head(10)
```

If `cluster_key` is not categorical, convert it:

```python
adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
```

---

## Step 2 — Pick candidate markers (and sanity-check them)

Markers are only useful if you check basic sanity:
- Are they expressed in the cluster?
- Are they absent/low elsewhere?
- Are you accidentally picking mitochondrial/ribosomal “QC markers”?

Example: pick top 5 markers per cluster:

```python
clusters = list(adata.obs[cluster_key].cat.categories)
top_markers_by_cluster = {}
for c in clusters:
    df = sc.get.rank_genes_groups_df(adata, group=c)
    top_markers_by_cluster[c] = df["names"].head(5).tolist()

top_markers_by_cluster
```

---

## Step 3 — The gene ID trap (read this before you search genes in the viewer)

When you search a gene in the Cellucid UI, the viewer uses the **gene IDs** it was given.

In AnnData mode (`show_anndata`), those gene IDs are typically:
- `adata.var_names` (default)
- OR a column you choose via `gene_id_column=...`

### Confirm what your gene IDs are

```python
print("var_names example:", list(adata.var_names[:5]))
print("var columns:", list(adata.var.columns)[:10])
```

Common situations:
- `adata.var_names` are Ensembl IDs (`ENSG...`) but you want to search by symbols (`GATA3`) → you need `gene_id_column`
- `adata.var_names` are symbols but contain duplicates → some tools accept this; viewers often do not

### If you have a gene-symbol column, use it

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    height=650,
    gene_id_column="gene_symbols",
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer
```

```{important}
Gene names must be distinct across `var`, and the ones that reach a manifest
must read as the value they store. Names are never filenames, so `HLA-DRB1/2`
is fine. Direct AnnData viewing serves every gene in `var`, so it holds all of
them to the display rule; `prepare()` holds only the genes `gene_identifiers`
actually exports. Both reject duplicate names anywhere in `var`.
```

---

## Step 4 — Color by a marker gene in the viewer (UI steps)

In the viewer:

1) Open the **Field selector**
2) Switch to **gene expression** (var fields) if the UI distinguishes obs vs genes
3) Search for a marker gene you picked (e.g. one from `top_markers_by_cluster`)
4) Select it and confirm the expression pattern is biologically plausible

```{note}
This tutorial does not assume a specific UI layout for the field selector (it evolves). If you get stuck, use the web-app docs for the exact UI walkthrough:
- {doc}`../../web_app/d_fields_coloring_legends/index`
```



---

## Step 5 — Export a figure (UI steps)

Exporting figures is a web-app feature. In the viewer:

1) Open the **Figure Export** panel
2) Choose output format:
   - **SVG** for vector editing (Illustrator/Inkscape) and publication
   - **PNG** for quick slides/docs
3) Confirm that:
   - the legend is visible and readable
   - the colormap matches your interpretation (and is colorblind-friendly when needed)
4) Export

For a complete figure-export walkthrough (including edge cases), see:
- {doc}`../../web_app/k_figure_export/index`

---

## Edge cases (gene markers + export)

### “My marker gene does not exist in the viewer search”

Likely causes:
- gene ID mismatch (you computed markers using symbols but the viewer uses Ensembl IDs)
- your gene list contains aliases (e.g. old symbols) not present in `adata.var_names`
- duplicates or `var_names` normalization changed the displayed IDs

How to confirm:
```python
gene = top_markers_by_cluster[clusters[0]][0]
gene in set(map(str, adata.var_names))
```

Fix:
- use `gene_id_column=...` when calling `show_anndata(...)`
- map/standardize gene IDs upstream (Ensembl↔symbol mapping)

### “Everything is gray / expression looks all-zero”

Likely causes:
- expression matrix is missing/empty
- you are viewing a layer that does not contain expression you expect
- values are all zeros for that gene in that subset (real biological or preprocessing artifact)

Confirm:
```python
gene = top_markers_by_cluster[clusters[0]][0]
idx = list(adata.var_names).index(gene) if gene in adata.var_names else None
idx
```

If you have Scanpy:
```python
import numpy as np
import scipy.sparse as sp

X = adata.X
v = X[:, idx].toarray().ravel() if sp.issparse(X) else X[:, idx]
float(np.nanmin(v)), float(np.nanmax(v))
```

### “My SVG export is huge / slow”

This is common when:
- you export very large point clouds
- you export at very high resolution or with many overlays

Use the figure-export troubleshooting:
- {doc}`../../web_app/k_figure_export/07_troubleshooting_figure_export`

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: rank_genes_groups crashes or returns nonsense

Likely causes:
- the data is not normalized/logged the way you think
- clusters are tiny (few cells) → unstable statistics
- batch effects dominate and clusters are not biological groups

How to confirm:
- check group sizes: `adata.obs[cluster_key].value_counts()`
- confirm preprocessing state (e.g., `adata.raw` usage, `adata.layers`)

Fix:
- use a pipeline consistent with your analysis goal
- filter clusters with too few cells
- consider batch-aware workflows (outside scope of this tutorial)

### Symptom: figure export looks different from what you saw

Likely causes:
- the export panel uses explicit settings (size, background, legend placement)
- you changed the active field or camera between “look” and “export”

Fix:
- open export panel, set options, and export immediately
- if you need reproducibility, save a session bundle (advanced; see {doc}`32_session_persistence_and_restoring_analysis_artifacts`)

---

## Next steps

- {doc}`13_compare_two_groups_and_interpret_results` (turn “markers” into a comparison story)
- {doc}`21_prepare_exports_with_quantization_and_compression` (export for speed + sharing)
- Web-app deep dive: {doc}`../../web_app/k_figure_export/index`
