# Notebook Tutorials

Three short, runnable walkthroughs. Each one simulates its data, so you can paste
the code into an R session (or an R Markdown notebook) and run it as it is.

If you start from a Seurat or SingleCellExperiment object, the recipes in
{doc}`../e_integrations_recipes/index` are shorter.

::::{grid} 1 1 2 2
:gutter: 3

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
