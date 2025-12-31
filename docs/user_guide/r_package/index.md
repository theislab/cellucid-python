# R Package Guide

This guide documents the **`cellucid` R package** (repo: `cellucid-r`): how to export single-cell data from R into the **Cellucid viewer export format** so you can explore it in the **Cellucid web app**.

```{note}
`cellucid-r` is an export-focused package. It does **not** embed the viewer in R yet (no RStudio widget / Shiny integration yet).

If you want notebook embedding, servers, and programmatic control, see the {doc}`../python_package/index` docs instead.

For the web UI experience after export (loading, filtering, analysis, sessions, community annotation), see {doc}`../web_app/index`.
```

## Where to start (pick one)

- “I have a Seurat object” → {doc}`e_integrations_recipes/01_seurat_recipe`
- “I have a SingleCellExperiment” → {doc}`e_integrations_recipes/02_singlecellexperiment_recipe`
- “I just have matrices/data.frames” → {doc}`e_integrations_recipes/03_raw_matrices_and_data_frames_recipe`
- “I want to understand what gets written to disk” → {doc}`b_concepts_mental_models/01_what_is_an_export_folder`
- “I want to export vector fields (velocity/drift overlays)” → {doc}`c_data_preparation_api/08_vector_fields_velocity_displacement`
- “I want a full end-to-end walkthrough” → {doc}`f_notebooks_tutorials/index`
- “Something failed” → {doc}`i_troubleshooting_index/index`

## Chapters (what each folder contains)

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`rocket;1.5em;sd-mr-1` Landing Pages
:link: a_landing_pages/index
:link-type: doc

Quick starts and “what is this?” explanations for mixed audiences.
:::

:::{grid-item-card} {octicon}`book;1.5em;sd-mr-1` Concepts & Mental Models
:link: b_concepts_mental_models/index
:link-type: doc

How R exports map to the browser viewer, plus dataset identity and reproducibility.
:::

:::{grid-item-card} {octicon}`package-dependencies;1.5em;sd-mr-1` Data Preparation API
:link: c_data_preparation_api/index
:link-type: doc

Deep documentation of `cellucid_prepare()` / `prepare()` including shapes, dtypes, edge cases, and file outputs.
:::

:::{grid-item-card} {octicon}`browser;1.5em;sd-mr-1` Viewing & Loading
:link: d_viewing_loading/index
:link-type: doc

How to open exported folders in the Cellucid web app, and how to host them for sharing.
:::

:::{grid-item-card} {octicon}`versions;1.5em;sd-mr-1` Integrations & Recipes
:link: e_integrations_recipes/index
:link-type: doc

Practical recipes for Seurat and SingleCellExperiment (and “raw matrix” workflows).
:::

:::{grid-item-card} {octicon}`device-desktop;1.5em;sd-mr-1` Notebook Tutorials
:link: f_notebooks_tutorials/index
:link-type: doc

Very detailed, notebook-style walkthroughs with lots of “why”, edge cases, and troubleshooting.
:::

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` API Reference Coverage
:link: g_api_reference_coverage/index
:link-type: doc

A map of what parts of the R API are documented where (useful for maintainers and contributors).
:::

:::{grid-item-card} {octicon}`tools;1.5em;sd-mr-1` Developer Docs
:link: h_developer_docs/index
:link-type: doc

Implementation notes, testing/release guidance, and codebase architecture for `cellucid-r` contributors.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting Index
:link: i_troubleshooting_index/index
:link-type: doc

Symptom → diagnosis → fix across installation, exports, and web app loading.
:::

::::

```{toctree}
:maxdepth: 2
:hidden:

a_landing_pages/index
b_concepts_mental_models/index
c_data_preparation_api/index
d_viewing_loading/index
e_integrations_recipes/index
f_notebooks_tutorials/index
g_api_reference_coverage/index
h_developer_docs/index
i_troubleshooting_index/index
```
