# Integrations & Recipes

This section is for practical “how do I get my data out of *this* object?” workflows.

`cellucid-r` intentionally has minimal dependencies, so it does **not** accept Seurat/SCE objects directly.
Instead, you extract:
- embeddings (UMAP 2D/3D)
- a latent space (often PCA)
- `obs` (cell metadata)
- optional gene expression + gene metadata
- optional connectivities / vector fields

Then you call `cellucid_prepare()`.

If you want the deep API docs (shapes, dtypes, file formats), see {doc}`../c_data_preparation_api/index`.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`versions;1.5em;sd-mr-1` Seurat Recipe
:link: 01_seurat_recipe
:link-type: doc

End-to-end export from a Seurat object (UMAP + PCA + metadata + expression + optional SNN graph).
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` SingleCellExperiment Recipe
:link: 02_singlecellexperiment_recipe
:link-type: doc

End-to-end export from an SCE object (reducedDims + colData + assays).
:::

:::{grid-item-card} {octicon}`database;1.5em;sd-mr-1` Raw Matrices + Data Frames
:link: 03_raw_matrices_and_data_frames_recipe
:link-type: doc

Minimal-dependency workflow: just matrices/data.frames, no Seurat/SCE required.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
