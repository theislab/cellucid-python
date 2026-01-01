# Beginner notebooks (wet lab friendly)

These tutorials are optimized for:
- wet lab scientists and non-technical users who want a clear “do this → see that” path
- computational users who want a quick refresher on the UI and common pitfalls

They use the simplest possible approach first (`show_anndata(...)`), then explain the “why” behind it.

```{important}
If you are running Cellucid in a notebook for the first time, read {doc}`00_start_here` first. It explains the Python-server + iframe mental model and the fastest troubleshooting checklist.
```

## Recommended order

1) {doc}`11_open_a_dataset_and_color_by_clusters`
2) {doc}`12_find_marker_genes_and_export_a_figure`
3) {doc}`13_compare_two_groups_and_interpret_results`

## What you will learn

- How to open a dataset in the Cellucid viewer from Python
- How to color by “clusters” / “cell type” / “batch” (metadata fields)
- How to select cells and interpret the selection
- How to export a figure without guessing (and how to report issues with screenshots)

## If you want the fastest path (copy/paste)

```python
from cellucid import show_anndata

# adata = ...  # your AnnData
viewer = show_anndata(adata, height=650)
viewer
```

Then in the viewer:
- choose a field to color by (clusters/cell type)
- select some cells (lasso)
- export a figure (SVG/PNG) if you want to share

For the complete, step-by-step version, start with {doc}`11_open_a_dataset_and_color_by_clusters`.

```{toctree}
:maxdepth: 1
:hidden:

11_open_a_dataset_and_color_by_clusters
12_find_marker_genes_and_export_a_figure
13_compare_two_groups_and_interpret_results
```
