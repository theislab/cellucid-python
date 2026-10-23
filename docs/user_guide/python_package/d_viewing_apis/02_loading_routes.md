# Every loading route

This page lists every way data reaches the Cellucid viewer, grouped by **who
reads the data**: the browser itself, or a Cellucid server that Python starts.
Every route reads AnnData — an `.h5ad` file, a Zarr store (folder or `.zip`),
or, through Python, an in-memory `AnnData` — with the same reader in the
browser.

If you want a decision tree instead of a table, go to {doc}`03_choose_your_workflow_decision_tree`.

## Definitions

- **AnnData**: `.h5ad`, a Zarr store (folder or `.zip`), or an in-memory `AnnData` object.
- **Cellucid-ready AnnData**: an AnnData written by `cellucid.prepare(...)`
  (or R's `cellucid_prepare()`): gene-major, chunked for range reads, and
  naming itself in `uns['cellucid']`
  ({doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`).
- **Lazy genes**: gene expression is fetched **on demand**, gene by gene,
  instead of loading the whole matrix. Every route below reads lazily.

## The browser reads it (no Python)

| Route | How you point Cellucid at the data | Formats | Best for |
|---|---|---|---|
| **Sample dataset** | the dataset menu, Cellucid's own catalog | Cellucid-ready AnnData | learning the UI with known-good data |
| **Open from link** | a link or id in **Open from link**, or `?data=…` ({doc}`../../web_app/b_data_loading/14_open_from_link`) | `.h5ad`, `.zarr.zip`, a Zarr store folder, a dataset catalog folder | opening data where it is already published; sharing it as one link |
| **Files** | the **H5AD**, **Zarr ZIP** and **Zarr folder** buttons, or a dropped file or folder | `.h5ad`, `.zarr.zip`, `.zarr` folder | viewing a file on your computer without Python |

## A Cellucid server serves it

| Route | Command | Formats | Best for |
|---|---|---|---|
| **Terminal** | `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset` (or a Zarr store) | `.h5ad`, Zarr | reliable viewing of large data; a remote machine over an SSH tunnel |
| **Python** | `cellucid.serve_anndata(<data.h5ad, a Zarr store or an AnnData>, dataset_name="My dataset", dataset_id="my-dataset")` | `.h5ad`, Zarr, in-memory | scripted server startup |
| **Jupyter** | `cellucid.show_anndata(<data.h5ad, a Zarr store or an AnnData>, dataset_name="My dataset", dataset_id="my-dataset")` | `.h5ad`, Zarr, in-memory | notebook exploration |

A Cellucid-ready AnnData declares its own name and id, so the identity
arguments are needed only for other AnnData.

The web app reaches a server through its **Server** field, or the server opens
it: `cellucid serve` opens your browser and a notebook embeds the viewer. The
viewer reads the store itself, in place, with range requests, exactly as it
reads a linked file; the server serves the store's bytes (an in-memory object
is first written once to a temporary `.h5ad`).

## What you should actually use (recommended defaults)

- **Local machine, terminal:** `cellucid serve …`
- **Notebook-based analysis:** `show_anndata`
- **Sharing with collaborators:** write a Cellucid-ready AnnData once
  (`prepare`), then publish the file anywhere **Open from link** reads, or
  serve it

## How the CLI maps to the Python API

`cellucid serve <path>` detects the format and serves it with
`serve_anndata(...)`:

- a file with the exact `.h5ad` suffix, or a `.zip` holding a Zarr store →
  served in place as one file
- a Zarr store folder (`zarr.json` for Zarr v3, `.zgroup` for Zarr v2, at its
  root) → served in place, key by key; the directory name does not need a
  `.zarr` suffix

Anything else is refused: `Unable to detect format for: <path>. Expected an
.h5ad file, a Zarr store folder, or a zipped Zarr store (.zip).`

Details: {doc}`04_cli_cellucid_serve_quickstart`.

## Where to read next (by route)

- The browser routes: {doc}`../../web_app/b_data_loading/index`
- Terminal: {doc}`04_cli_cellucid_serve_quickstart`
- Python servers: {doc}`05_python_serve_anndata_quickstart`
- Jupyter: {doc}`06_jupyter_show_anndata_quickstart`

## Edge cases (high-signal)

- Remote notebooks (JupyterHub, cloud) often require a proxy URL (HTTPS) instead of direct `http://127.0.0.1:<port>`: {doc}`10_notebook_widget_mode_advanced`.

## Troubleshooting

If anything in this table “should work” but doesn’t, start with:
- {doc}`15_troubleshooting_viewing` (symptom → diagnosis → fix)
