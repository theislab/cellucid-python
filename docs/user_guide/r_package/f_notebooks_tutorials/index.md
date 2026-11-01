# Notebook Tutorials

Short, runnable walkthroughs. Paste the code into an R session (RStudio,
Positron, an R Markdown notebook) or open the notebooks (`.ipynb`, R kernel)
in Jupyter or Google Colab, and run them as they are.

## Cellucid in your pipeline (notebooks, each with a Python twin)

Show an object in one line, work in the viewer, and bring what you made back
onto the Seurat or SingleCellExperiment object in one call — or push what you
computed into the viewer. Each opens in Google Colab from its badge.

| Topic | R | Python twin |
|---|---|---|
| View an object in a notebook (Jupyter, Colab, RStudio) | {doc}`07_view_in_a_notebook` | {doc}`../../python_package/f_notebooks_tutorials/41_show_a_dataset_in_any_notebook` |
| Pull highlights, labels and scores back | {doc}`08_pull_everything_back` | {doc}`../../python_package/f_notebooks_tutorials/42_pull_everything_back` |
| Push clusters, refine, pull back (`FindClusters`) | {doc}`09_push_clusters_refine_and_pull` | {doc}`../../python_package/f_notebooks_tutorials/43_push_clusters_refine_and_pull` |
| No notebook: the web app's exports in one call | {doc}`10_round_trip_without_a_notebook` | {doc}`../../python_package/f_notebooks_tutorials/44_round_trip_without_a_notebook` |
| Gene sets: GMT in, signature scores out | {doc}`11_gene_sets_round_trip` | {doc}`../../python_package/f_notebooks_tutorials/45_gene_sets_round_trip` |
| A selection, re-analysed | {doc}`12_subset_and_reanalyse` | {doc}`../../python_package/f_notebooks_tutorials/46_export_an_h5ad_subset_and_reanalyse` |
| Subsampled views | {doc}`13_subsampled_views` | {doc}`../../python_package/f_notebooks_tutorials/47_subsampled_views` |

In a Jupyter notebook R hears the viewer between cells: show the viewer in one
cell, work in it, and call `cellucid_pull()` in the next
({doc}`../../python_package/f_notebooks_tutorials/06_cellucid_in_your_pipeline_python_and_r`
explains the two-cell pattern and lists the calls side by side).

## Walkthroughs

Three start from a Seurat or SingleCellExperiment object (the first three
cards below); three simulate their data from matrices (Beginner, Intermediate
and Advanced).

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`device-desktop;1.5em;sd-mr-1` Explore a Seurat Object
:link: 04_explore_a_seurat_object_from_r
:link-type: doc

`cellucid_view(seu)`, colour and highlight from the console, and the cells you
lasso back in R.
:::

:::{grid-item-card} {octicon}`mark-github;1.5em;sd-mr-1` Share a SingleCellExperiment
:link: 05_share_a_singlecellexperiment_on_github
:link-type: doc

`cellucid_prepare(data = sce)` to a Zarr folder, published on GitHub and
opened from one link.
:::

:::{grid-item-card} {octicon}`stack;1.5em;sd-mr-1` CITE-seq: RNA and ADT
:link: 06_citeseq_seurat_to_mudata
:link-type: doc

A Seurat object of two assays: the Modality select, and one Cellucid-ready
MuData for others.
:::

:::{grid-item-card} {octicon}`person;1.5em;sd-mr-1` Beginner (Wet Lab Friendly)
:link: 01_beginner_notebooks_wet_lab_friendly
:link-type: doc

“I just want to see my cells.” One embedding, cell metadata, one file to open.
:::

:::{grid-item-card} {octicon}`graph;1.5em;sd-mr-1` Intermediate (Computational)
:link: 02_intermediate_notebooks_computational_workflows
:link-type: doc

Gene expression, gene symbols, a gene panel, launch defaults, and citation metadata.
:::

:::{grid-item-card} {octicon}`tools;1.5em;sd-mr-1` Advanced (Expert / Developer)
:link: 03_advanced_notebooks_expert_developer
:link-type: doc

Vector fields, a neighbour graph, and reading the store back value for value.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
