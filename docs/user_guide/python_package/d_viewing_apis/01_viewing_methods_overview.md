# Viewing methods overview

This page is the “map” of how you can open Cellucid against your data using `cellucid-python`.

If you just want to pick the right command, go to {doc}`03_choose_your_workflow_decision_tree`.

## The mental model (1 minute)

Cellucid always needs two things:

1) the **viewer UI** (HTML/JS) that runs in your browser, and
2) the **dataset's bytes**: an AnnData file or store the viewer reads itself, with range requests.

`cellucid-python` establishes the exact viewer generation and starts a **local
HTTP server** that serves both UI and data from the same origin. Open the exact
**Viewer URL** printed by that server:

```text
http://127.0.0.1:8765/?anndata=true
```

The `anndata=true` marker declares the user-served source before the first
paint, and the UI can load the data without cross-origin issues.

## At a glance

**Audience**
- Wet lab / non-technical: use the CLI quickstart and stop at “What success looks like”.
- Computational: pay attention to data formats (`.h5ad` vs Zarr store vs in-memory, prepared or not).
- Power users: pay attention to SSH tunnels, HTTPS notebooks, and caching.

**Prerequisites**
- `pip install cellucid`
- One data input: `.h5ad`, Zarr store (folder or `.zip`), any other format
  Cellucid reads ({doc}`02_loading_routes`), or an in-memory `AnnData`/`MuData`

## The three core entry points

| You want to… | Best entry point | You pass… | Typical place you run it |
|---|---|---|---|
| Open the viewer from the terminal | `cellucid serve …` | a path | terminal |
| Start a server from a Python script | `serve_anndata(…)` | a path, an `AnnData` or a `MuData` | Python (notebook optional) |
| Embed the viewer in a notebook | `show_anndata(…)` | a path, an `AnnData` or a `MuData` | Jupyter / JupyterLab (VS Code and Colab untested) |

### Server mode vs notebook mode

- **Server mode** (`cellucid serve`, `serve_anndata`) is a normal “keep this terminal running” server.
- **Notebook mode** (`show_anndata`) starts the same kind of server, but also renders an **iframe** in the notebook and gives you a `viewer` object for hooks/commands.

## Data inputs

| Data input | How you get it | Recommended APIs | Why you’d choose it |
|---|---|---|---|
| Cellucid-ready AnnData | `cellucid.prepare(...)` | `show_anndata("pbmc_cellucid.h5ad")` / `serve_anndata(...)` / `cellucid serve pbmc_cellucid.h5ad` | fastest; names itself; easiest to share and host |
| `.h5ad` path | AnnData HDF5 file | `show_anndata("…h5ad", dataset_name="My dataset", dataset_id="my-dataset")` / `serve_anndata("…h5ad", dataset_name="My dataset", dataset_id="my-dataset")` / `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset` | convenient + scalable: the viewer reads the file in place, with range requests |
| Zarr store path (folder or `.zip`) | AnnData Zarr store | `show_anndata("…zarr", dataset_name="My dataset", dataset_id="my-dataset")` / `serve_anndata("…zarr", dataset_name="My dataset", dataset_id="my-dataset")` / `cellucid serve data.zarr --dataset-name "My dataset" --dataset-id my-dataset` | convenient + scalable: read in place, key by key |
| In-memory `AnnData` | you already loaded it | `show_anndata(adata, dataset_name="My dataset", dataset_id="my-dataset")` / `serve_anndata(adata, dataset_name="My dataset", dataset_id="my-dataset")` | convenient; written once to a temporary `.h5ad` first |
| MuData (`.h5mu`, MuData Zarr, or in memory) | `mudata`, or `prepare(data=mdata, ...)` | the same calls; `modality="rna"` (`--modality rna`) picks the modality it opens on | several modalities over the same cells; an object is written once to a temporary `.h5mu` |
| Loom, 10x `.h5` / `outs/`, MatrixMarket folder, `.rds`/`.rda` | other tools | a path, with `dataset_name`/`dataset_id` (Python reads no identity from these) | no conversion: the browser reads them itself ({doc}`02_loading_routes`) |

## Where these fit among every loading route

Data reaches the viewer either through the browser itself (a sample, a link,
a file) or through a Cellucid server, which these APIs start.

- Every route: {doc}`02_loading_routes`
- Web-app perspective: {doc}`../../web_app/b_data_loading/01_loading_options_overview`

## Common “gotchas” (read once)

- **Large dataset?** Pass a path, which the viewer reads in place; an in-memory
  object is written once to a temporary `.h5ad` first. A file written by
  `prepare()` is gene-major, so every gene is one small read
  ({doc}`../c_data_preparation_api/10_performance_and_storage`).
- **Remote machine?** Don’t bind to `0.0.0.0` unless you mean to. Use an SSH tunnel first: {doc}`12_remote_servers_ssh_tunneling_and_cloud`.
- **HTTPS notebook?** Direct `http://127.0.0.1:<port>` iframes can be blocked.
  Expose the port through HTTPS and pass its exact base as
  `client_server_url=`: {doc}`10_notebook_widget_mode_advanced`.
- **Viewer source blocked?** Viewer-serving startup stops instead of using a
  previous generation. See {doc}`09_server_mode_advanced` and
  {doc}`15_troubleshooting_viewing`.

## Next steps

- Pick a workflow: {doc}`03_choose_your_workflow_decision_tree`
- Terminal flow: {doc}`04_cli_cellucid_serve_quickstart`
- Notebook flow: {doc}`06_jupyter_show_anndata_quickstart`
