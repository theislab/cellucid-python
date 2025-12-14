# User Guide

Step-by-step tutorials to get you started with Cellucid visualization.

---

## Getting Started

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` Loading Options Overview
:link: 01_loading_options_overview
:link-type: doc

Understand all 14 loading options and choose the best approach for your workflow: in-memory, h5ad, zarr, or pre-exported.
:::

:::{grid-item-card} {octicon}`file-directory;1.5em;sd-mr-1` Local Demo
:link: 02_local_demo_tutorial
:link-type: doc

Export your data and open it directly in a browser with no server required.
:::

::::

---

## Viewing Methods

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`upload;1.5em;sd-mr-1` Browser File Picker
:link: 03_browser_file_picker_tutorial
:link-type: doc

Load datasets directly from your file system using the browser's native file picker API.
:::

:::{grid-item-card} {octicon}`server;1.5em;sd-mr-1` Server Mode
:link: 04_server_tutorial
:link-type: doc

Run a local HTTP server for larger datasets with efficient lazy loading of gene expression.
:::

:::{grid-item-card} {octicon}`code;1.5em;sd-mr-1` Jupyter Integration
:link: 05_jupyter_tutorial
:link-type: doc

Visualize AnnData objects directly within Jupyter notebooks with embedded interactive widgets.
:::

::::

---

## Quick Comparison

| Method | Best For | Requires Server |
|--------|----------|-----------------|
| Local Demo | Quick visualization, sharing static files | No |
| File Picker | Desktop viewing without setup | No |
| Server Mode | Large datasets, lazy gene loading | Yes |
| Jupyter | Interactive analysis workflows | No |

---

```{toctree}
:maxdepth: 1
:hidden:

01_loading_options_overview
02_local_demo_tutorial
03_browser_file_picker_tutorial
04_server_tutorial
05_jupyter_tutorial
```
