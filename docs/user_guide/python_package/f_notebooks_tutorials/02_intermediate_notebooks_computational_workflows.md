# Intermediate notebooks (computational workflows)

These tutorials assume you can load an `AnnData` and are comfortable in Python, but want workflows that are:
- reproducible (stable dataset identity)
- scalable (large datasets, lazy loading, fewer browser crashes)
- shareable (Cellucid-ready AnnData files you can hand to collaborators)

## Recommended order

1) {doc}`21_prepare_a_cellucid_ready_anndata`
2) {doc}`22_large_dataset_server_mode_and_lazy_gene_expression`
3) {doc}`23_programmatic_highlighting_and_selection_callbacks`
4) For a runnable real dataset recipe: {doc}`04_real_world_dataset_recipes_gallery`

## What you will learn

- When to open the AnnData you have and when to `prepare(...)` it first
- What a Cellucid-ready AnnData changes for size and speed
- How to run server mode safely (local/remote) and what “lazy gene expression” means
- How to react to selection events and drive the UI programmatically

```{tip}
If you are preparing data for colleagues, prefer:

`prepare(...)` → share the file → `show_anndata("my_dataset.h5ad")`

This avoids “my notebook had hidden state” problems and makes the dataset identity stable.
```

```{toctree}
:maxdepth: 1
:hidden:

21_prepare_a_cellucid_ready_anndata
22_large_dataset_server_mode_and_lazy_gene_expression
23_programmatic_highlighting_and_selection_callbacks
```
