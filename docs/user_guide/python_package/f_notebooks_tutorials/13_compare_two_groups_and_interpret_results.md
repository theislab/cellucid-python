# Compare two groups and interpret results

This tutorial shows how to compare two groups of cells *without fooling yourself*.

You will learn:
- how to define “Group A” and “Group B” in a way you can reproduce
- how to perform a basic comparison in Python
- how to interpret results and avoid the most common traps (batch effects, tiny groups, “p<1e-50” illusions)
- how to use Cellucid selections as inputs to analysis (optional)

---

## At a glance

**Audience**
- Wet lab / beginner: focus on Option A (cluster-based groups) and the interpretation checklist.
- Computational users: also read the “selection-based groups” and edge-case sections.

**Time**
- ~30–60 minutes

**Prerequisites**
- A dataset loaded as `AnnData`
- A cluster/label field (from {doc}`11_open_a_dataset_and_color_by_clusters`)
- Optional: `scanpy` for differential expression convenience

---

## Choose how you want to define groups

There are two common ways to define groups:

### Option A (recommended for beginners): compare clusters/labels

Group A = all cells with `adata.obs[cluster_key] == "X"`

Pros:
- reproducible
- easy to communicate
- stable across notebooks and sessions

Cons:
- depends on cluster quality and preprocessing choices

### Option B (optional): compare manual selections from the viewer

Group A = cells you lasso/select in the UI

Pros:
- extremely intuitive (“these cells here vs those cells there”)
- great for exploratory hypotheses

Cons:
- easy to introduce selection bias
- harder to reproduce unless you save a session bundle

---

## Option A — Compare two clusters (fast path)

### Step 1 — Pick your cluster key and two groups

```python
cluster_key = "leiden"  # change if needed
adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

counts = adata.obs[cluster_key].value_counts()
counts
```

Pick two labels:

```python
group_a_label = counts.index[0]
group_b_label = counts.index[1]
group_a_label, group_b_label
```

### Step 2 — Create boolean masks

```python
mask_a = (adata.obs[cluster_key].astype(str) == str(group_a_label)).to_numpy()
mask_b = (adata.obs[cluster_key].astype(str) == str(group_b_label)).to_numpy()

print("Group A cells:", int(mask_a.sum()))
print("Group B cells:", int(mask_b.sum()))
```

```{important}
If either group has very few cells (e.g. < 30–50), most “marker” results will be unstable and highly sensitive to preprocessing choices.
```

### Step 3 — A simple comparison you can trust (mean expression for a gene list)

If you already have a short gene list (e.g. markers), you can compare means/medians:

```python
genes = list(adata.var_names[:10])  # replace with genes you care about

# Dense-safe example:
import numpy as np
import scipy.sparse as sp

X = adata.X
gene_ix = [list(adata.var_names).index(g) for g in genes]

XA = X[mask_a][:, gene_ix]
XB = X[mask_b][:, gene_ix]

XA = XA.toarray() if sp.issparse(XA) else np.asarray(XA)
XB = XB.toarray() if sp.issparse(XB) else np.asarray(XB)

mean_a = XA.mean(axis=0)
mean_b = XB.mean(axis=0)

list(zip(genes, mean_a, mean_b))
```

### Step 4 — Differential expression (Scanpy convenience)

If you want a ranked list for just these two groups:

```python
import scanpy as sc

tmp = adata.copy()
tmp.obs["comparison_group"] = "other"
tmp.obs.loc[mask_a, "comparison_group"] = "A"
tmp.obs.loc[mask_b, "comparison_group"] = "B"
tmp.obs["comparison_group"] = tmp.obs["comparison_group"].astype("category")

tmp = tmp[tmp.obs["comparison_group"].isin(["A", "B"])].copy()

sc.tl.rank_genes_groups(tmp, groupby="comparison_group", method="wilcoxon")
sc.get.rank_genes_groups_df(tmp, group="A").head(15)
```

---

## Option B — Compare two manual selections (viewer → Python hooks)

This is optional, but very powerful when you want “this region vs that region”.

### Step 1 — Show the dataset

```python
from cellucid import show_anndata

viewer = show_anndata(adata, height=650)
viewer
```

### Step 2 — Capture two selections from the UI

```{warning}
This pattern is intentionally simple, not “production quality”. It stores the first selection as Group A and the second selection as Group B.

For robust event handling (timeouts, multiple viewers, deduping), see {doc}`23_programmatic_highlighting_and_selection_callbacks`.
```

```python
group_a_cells = None
group_b_cells = None

@viewer.on_selection
def _capture_selection(event):
    global group_a_cells, group_b_cells
    cells = event.get("cells", [])
    if group_a_cells is None:
        group_a_cells = cells
        print(f"Stored Group A: {len(cells)} cells")
    else:
        group_b_cells = cells
        print(f"Stored Group B: {len(cells)} cells")
```

In the viewer UI:
1) make a selection for Group A (lasso)
2) make a second selection for Group B

Confirm in Python:

```python
len(group_a_cells or []), len(group_b_cells or [])
```

### Step 3 — Use those indices to subset `adata`

```python
import numpy as np

mask_a = np.zeros(adata.n_obs, dtype=bool)
mask_b = np.zeros(adata.n_obs, dtype=bool)
mask_a[group_a_cells] = True
mask_b[group_b_cells] = True

subset_a = adata[mask_a].copy()
subset_b = adata[mask_b].copy()
subset_a, subset_b
```

---

## Interpretation checklist (read this even if you’re experienced)

When you get “top marker genes” or “differential expression”:

1) **Effect size beats p-value.**
   - Large datasets produce extremely small p-values even for tiny effects.
2) **Check group sizes.**
   - Tiny groups are noisy and often reflect selection artifacts.
3) **Watch for confounders.**
   - Batch, donor, sample prep, cell cycle, mitochondrial content can dominate.
4) **QC markers are not biology (usually).**
   - Ribosomal/mitochondrial genes often indicate technical differences.
5) **Selections can be biased.**
   - If you manually select “interesting-looking” regions, you are conditioning on the visualization.

```{note}
This tutorial is not a full statistical DE guide. The goal is to help you avoid the most common interpretation failures in exploratory workflows.
```

---

## Screenshot placeholders (optional but helpful)

<!-- SCREENSHOT PLACEHOLDER
ID: python-notebooks-two-selections-captured
Suggested filename: highlighting_selection/00_two-selections-and-python-callback-output.png
Where it appears: User Guide → Python Package → Notebooks/Tutorials → 13_compare_two_groups_and_interpret_results.md
Capture:
  - UI location: notebook output with viewer + Python output below
  - State prerequisites: two selections have been made; Python printed “Stored Group A/B”
  - Action to reach state: make two selections in the viewer and watch callback prints in notebook
Crop:
  - Include: viewer canvas showing selected cells (highlighted/outlined)
  - Include: the notebook output area showing the printed “Stored Group A/B” lines
Redact:
  - Remove: private dataset names, file paths, usernames
Alt text:
  - Viewer selection made in notebook and Python output confirming captured groups.
Caption:
  - Selection hooks send the selected cell indices back to Python so you can run analysis on exactly what you selected in the UI.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for capturing two selections via hooks.
:width: 100%

Two selections captured from the viewer, with Python output confirming the stored group sizes.
```

---

## Edge cases (common problems)

### “My two groups overlap”

This can happen if you select overlapping regions, or if your UI selection tool unions selections.

How to confirm:
```python
import numpy as np
overlap = np.logical_and(mask_a, mask_b).sum()
overlap
```

Fix:
- decide whether overlap is meaningful; otherwise remove overlap explicitly:
  - `mask_b = np.logical_and(mask_b, ~mask_a)`

### “Selection indices don’t match my AnnData ordering”

Selection indices always refer to the **row order the viewer is currently using**.

In normal usage, that is the same as `adata.obs_names` order at the time you called `show_anndata(...)`.

If you subset/reorder `adata` after opening the viewer, you can desynchronize your mental model.

Best practice:
- treat `adata` as immutable after creating `viewer`, or recreate the viewer after reordering

### “I can’t reproduce my selection later”

Manual selections are inherently ephemeral unless you capture state.

If you need reproducibility, use:
- cluster labels (Option A)
- or save a session bundle (advanced; see {doc}`32_session_persistence_and_restoring_analysis_artifacts`)

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Hooks don’t fire (I select cells but Python prints nothing)”

Likely causes:
- the viewer never reached “ready”
- browser cannot reach the kernel-side server port (remote notebook without proxy/tunnel)
- corporate network/proxy strips requests

How to confirm:
```python
viewer.debug_connection()
```

Fix:
- use Jupyter Server Proxy or SSH port forwarding (see {doc}`22_large_dataset_server_mode_and_lazy_gene_expression`)
- verify server health endpoint is reachable
- for robust patterns, see {doc}`23_programmatic_highlighting_and_selection_callbacks`

### Symptom: “DE results are dominated by mitochondrial/ribosomal genes”

Likely causes:
- you are comparing groups with different QC profiles (not biology)

Fix:
- regress/correct QC covariates (analysis-specific)
- filter QC genes from marker lists when telling a biological story

---

## Next steps

- {doc}`23_programmatic_highlighting_and_selection_callbacks` (robust hooks patterns)
- {doc}`21_prepare_exports_with_quantization_and_compression` (shareable exports)
- Web app selection details: {doc}`../../web_app/f_highlighting_selection/index`
