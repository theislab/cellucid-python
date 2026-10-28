# Viewing and Loading

This section answers: “How do I look at my object from R?” and “I wrote a store
with `cellucid_prepare()` — how do I view it, and how do I share it?”

`cellucid_view(object)` opens a Seurat, SingleCellExperiment or
SummarizedExperiment object (or any file the web app reads) in the RStudio or
Positron Viewer pane, a browser, or -- in a Jupyter notebook -- the cell's
output, with no file written and no Python; `cellucid_serve()` serves it to a
browser elsewhere.

`cellucid_prepare()` writes a **Cellucid-ready AnnData**: a standard AnnData Zarr
store, either a `.zarr` folder (`format = "zarr"`) or one `.zarr.zip` file (`format = "zarr.zip"`). The
web app opens it like any other AnnData store, so there is no separate format to
load.

- **On your own computer:** write a `.zarr.zip` and open it with **Open… ›
  File…**.
- **Shared with others:** put the `.zarr` folder or the `.zarr.zip` file on a web
  host (GitHub, Hugging Face, Zenodo, a bucket) and open it with the **Link or
  server** bar.

```{note}
The web app's own loading documentation is in {doc}`../../web_app/index`; every
way to open data from a web address is in
{doc}`../../web_app/b_data_loading/14_open_from_link`.
```

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`upload;1.5em;sd-mr-1` Open a Store in the Web App
:link: 01_open_stores_in_cellucid_web_app
:link-type: doc

Open a `.zarr.zip` from your computer, or a store from a web address.
:::

:::{grid-item-card} {octicon}`share;1.5em;sd-mr-1` Host a Store for Sharing
:link: 02_host_stores_for_sharing
:link-type: doc

Publish a store on GitHub, Hugging Face, Zenodo or any web host, and share one link.
:::

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Validate a Store & Debug Loading
:link: 03_validate_stores_and_debug_loading
:link-type: doc

Read the store back in R and in Python's anndata, and diagnose loading failures.
:::



:::{grid-item-card} {octicon}`device-desktop;1.5em;sd-mr-1` View in RStudio or Positron
:link: 04_view_in_rstudio_positron
:link-type: doc

`cellucid_view(seu)`: the object in the Viewer pane, from memory.
:::

:::{grid-item-card} {octicon}`code-square;1.5em;sd-mr-1` View in a Jupyter Notebook
:link: 08_view_in_jupyter_notebooks
:link-type: doc

`cellucid_view(seu)` in IRkernel: the viewer in the cell's output, no port.
:::

:::{grid-item-card} {octicon}`server;1.5em;sd-mr-1` Serve from a Cluster or over SSH
:link: 05_serve_from_hpc_and_ssh
:link-type: doc

`cellucid_serve()`, SSH tunnels, RStudio Server and Workbench proxies.
:::

:::{grid-item-card} {octicon}`zap;1.5em;sd-mr-1` Hooks, Selections and Sessions
:link: 06_hooks_selection_and_sessions
:link-type: doc

Commands from R, events back, a selection or the session in R.
:::

:::{grid-item-card} {octicon}`stack;1.5em;sd-mr-1` Multimodal Assays
:link: 07_multimodal_assays
:link-type: doc

RNA + ADT or ATAC as modalities of a Cellucid-ready MuData.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
