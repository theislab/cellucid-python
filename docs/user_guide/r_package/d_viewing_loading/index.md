# Viewing and Loading Stores

This section answers: “I wrote a store with `cellucid_prepare()` — how do I view
it, and how do I share it?”

`cellucid_prepare()` writes a **Cellucid-ready AnnData**: a standard AnnData Zarr
store, either a `.zarr` folder (`format = "zarr"`) or one `.zarr.zip` file (`format = "zarr.zip"`). The
web app opens it like any other AnnData store, so there is no separate format to
load.

- **On your own computer:** write a `.zarr.zip` and open it with the **Zarr ZIP**
  button.
- **Shared with others:** put the `.zarr` folder or the `.zarr.zip` file on a web
  host (GitHub, Hugging Face, Zenodo, a bucket) and open it with **Open from
  link**.

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

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
