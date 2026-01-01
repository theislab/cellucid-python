# Python Package Guide

This guide documents the **`cellucid` Python package**: how to prepare datasets, run servers, embed the viewer in notebooks, and integrate Python with the web UI (hooks/events).

If you’re looking for the end-user web UI docs (loading data, filtering, analysis, annotation), start with {doc}`../web_app/index`.

## Where to start (pick one)

- “I have an AnnData and want to view it” → {doc}`d_viewing_apis/index` (and {doc}`../web_app/b_data_loading/05_jupyter_tutorial`)
- “I want the fastest web viewing” → {doc}`c_data_preparation_api/index` (exports via `prepare`)
- “I want notebook embedding + programmatic control” → {doc}`d_viewing_apis/index` + {doc}`e_jupyter_hooks/index`
- “I’m new and want a guided path” → {doc}`a_landing_pages/index`
- “Something failed” → {doc}`i_troubleshooting_index/index`

## Chapters (what each folder contains)

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`rocket;1.5em;sd-mr-1` Landing Pages
:link: a_landing_pages/index
:link-type: doc

Quickstarts and entry points: installation, first run, and “what should I do next?” pages.
:::

:::{grid-item-card} {octicon}`book;1.5em;sd-mr-1` Concepts & Mental Models
:link: b_concepts_mental_models/index
:link-type: doc

How the Python package maps data → exports/servers → the web viewer, plus core terminology you’ll see across docs.
:::

:::{grid-item-card} {octicon}`package-dependencies;1.5em;sd-mr-1` Data Preparation API
:link: c_data_preparation_api/index
:link-type: doc

`prepare(...)` and export-time options (fields, gene IDs, embeddings, vector fields) for reproducible, high-performance viewing.
:::

:::{grid-item-card} {octicon}`globe;1.5em;sd-mr-1` Viewing APIs
:link: d_viewing_apis/index
:link-type: doc

How to view data: `serve(...)`, `serve_anndata(...)`, `show(...)`, and `show_anndata(...)` (Jupyter embedding and server mode).
:::

:::{grid-item-card} {octicon}`plug;1.5em;sd-mr-1` Jupyter Hooks (Python ↔ Frontend)
:link: e_jupyter_hooks/index
:link-type: doc

React to selection/hover/click events and send commands back to the viewer (highlight, color-by, visibility).
:::

:::{grid-item-card} {octicon}`device-desktop;1.5em;sd-mr-1` Notebook Tutorials
:link: f_notebooks_tutorials/index
:link-type: doc

Long-form notebook-style tutorials and “copy/paste” workflows for common analysis and visualization patterns.
:::

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` API Reference Coverage
:link: g_api_reference_coverage/index
:link-type: doc

A map of what parts of the Python API are documented where (useful for maintainers and contributors).
:::

:::{grid-item-card} {octicon}`tools;1.5em;sd-mr-1` Developer Docs
:link: h_developer_docs/index
:link-type: doc

Implementation notes, architectural decisions, and debugging patterns for the Python package.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting Index
:link: i_troubleshooting_index/index
:link-type: doc

Symptom → diagnosis → fix for installation, exports, servers, Jupyter embedding, and performance issues.
:::

:::{grid-item-card} {octicon}`code;1.5em;sd-mr-1` API Reference
:link: g_api_reference_coverage/api/index
:link-type: doc

Complete reference for all functions and classes in the Cellucid Python package.
:::

::::

```{toctree}
:maxdepth: 2
:hidden:

a_landing_pages/index
b_concepts_mental_models/index
c_data_preparation_api/index
d_viewing_apis/index
e_jupyter_hooks/index
examples/index
f_notebooks_tutorials/index
g_api_reference_coverage/index
h_developer_docs/index
i_troubleshooting_index/index
```
