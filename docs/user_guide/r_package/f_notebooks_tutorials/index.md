# Notebook Tutorials

Short, runnable walkthroughs. Three start from a Seurat or
SingleCellExperiment object (the first three cards below); three simulate
their data from matrices (Beginner, Intermediate and Advanced).
Paste the code into an R session (or an R Markdown notebook) and run it as it
is.

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
