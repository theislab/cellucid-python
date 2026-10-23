# Integrations & Recipes

This section is for practical “how do I get my data out of *this* object?”
workflows.

`cellucid-r` has minimal dependencies, so `cellucid_prepare()` takes plain R
inputs rather than a Seurat or SingleCellExperiment object. You extract:

- the embeddings (UMAP, t-SNE, PCA, …) — each becomes `obsm[<key>]`;
- the cell metadata, a `data.frame` — `obs`;
- optionally gene expression, **cells x genes**, and gene metadata — `X`, `var`;
- optionally a symmetric neighbour graph — `obsp['connectivities']`;
- optionally vector fields — `obsm['<field>_<basis>']` for an embedding `X_<basis>` (`velocity_umap`).

Then one `cellucid_prepare()` call writes a Cellucid-ready AnnData store.

For every argument in detail, see {doc}`../c_data_preparation_api/index`.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`versions;1.5em;sd-mr-1` Seurat Recipe
:link: 01_seurat_recipe
:link-type: doc

From a Seurat v5 object: reductions, metadata, expression, and the SNN graph.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` SingleCellExperiment Recipe
:link: 02_singlecellexperiment_recipe
:link-type: doc

From an SCE object: reducedDims, colData, rowData, and an assay.
:::

:::{grid-item-card} {octicon}`database;1.5em;sd-mr-1` Raw Matrices + Data Frames
:link: 03_raw_matrices_and_data_frames_recipe
:link-type: doc

No Seurat or SCE: matrices and data frames only.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
