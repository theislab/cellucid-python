# Viewing APIs (cellucid serve / serve_anndata / show_anndata + loading options)

Cellucid is a **web app** (WebGL viewer UI). `cellucid-python` gives you **Python + CLI entry points** to open that viewer against *your data*.

There are two big “modes”:

- **Server mode** (terminal/browser): `cellucid serve …`, `serve_anndata(…)`
- **Notebook mode** (a widget in the cell's output, no server): `show_anndata(…)`

All of them take any file or folder Cellucid reads — AnnData (`.h5ad`, a Zarr
store folder or `.zip`), MuData (`.h5mu`, MuData Zarr), Loom, a 10x `.h5` file
or Cell Ranger `outs/` folder, a MatrixMarket folder, an R object
(`.rds`/`.rda`/`.RData`, gzipped or not) — and `serve_anndata`/`show_anndata`
also an in-memory `AnnData` or `MuData` ({doc}`02_loading_routes`). The viewer
reads a file in place, with range requests where the format allows it; an
in-memory AnnData, and an AnnData file not in priority order, are served in
priority order, made as the viewer reads them, so a subsample sends only its
cells (a MuData object is written once to a temporary `.h5mu` first). A file
written by {func}`~cellucid.prepare` (a Cellucid-ready AnnData) reads fastest
and names itself, so it needs no `dataset_name`/`dataset_id`.

## At a glance

**Audience**
- Wet lab / non-technical: copy/paste the quickstarts and use the decision tree.
- Computational: focus on data format tradeoffs: `.h5ad` files, Zarr stores,
  and in-memory objects, prepared or not.
- Power user / developer: focus on networking (SSH tunnels), caching, and debug endpoints.

**Time**
- First successful view: ~5 minutes
- Full read (with edge cases): ~30–45 minutes

**Prerequisites**
- `pip install cellucid`
- One of: a `.h5ad`, a Zarr store (folder or `.zip`), any other file or folder
  Cellucid reads ({doc}`02_loading_routes`), or an in-memory `AnnData`/`MuData`

## Where to start

- Not sure what to do → {doc}`03_choose_your_workflow_decision_tree`
- Want a terminal “just show me the viewer” flow → {doc}`04_cli_cellucid_serve_quickstart`
- Want a scripted Python server (non-notebook) → {doc}`05_python_serve_anndata_quickstart`
- Want the viewer in a notebook cell (JupyterLab is tested; Notebook 7, VS Code,
  Colab, marimo and Databricks are expected to work, see
  {doc}`../e_jupyter_hooks/03_supported_environments_matrix`) →
  {doc}`06_jupyter_show_anndata_quickstart`
- Something failed → {doc}`15_troubleshooting_viewing`

## Quick start (copy/paste)

### Terminal (recommended “default”)

```bash
pip install cellucid
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1
# A Zarr store (folder or .zip) takes the same flags; a file written by
# cellucid.prepare() names itself and needs neither.
```

Open the printed **Viewer URL** in your browser. Stop with **Ctrl+C**.

### Notebook (fastest analysis loop)

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

### A Cellucid-ready AnnData (fastest + names itself)

```python
from cellucid import show_anndata

viewer = show_anndata("pbmc_cellucid.h5ad", height=600)  # written by cellucid.prepare()
```

## Read once: how the exact viewer generation is served

```{important}
By design, the Python server serves **both**:

1) your data: the AnnData store's bytes under `/_cellucid/anndata/`, which the viewer reads itself with range requests, and
2) the **viewer UI** (HTML/JS/CSS) by downloading it from `https://www.cellucid.com` and caching it locally.

This verified same-origin generation approach:
- avoids HTTPS→HTTP mixed-content problems in notebooks and in strict browsers,
- keeps the UI and data on the **same origin** (simplifies networking),
- and requires the configured viewer source at every server startup.

If the source generation cannot be fetched and verified, startup reports the
exact failure and does not bind the viewer server. Configure the generation
directory with `web_cache_dir=...`.
```

## Relationship to the web app + other repos

- **Cellucid (web app)** is the UI you interact with.
- **cellucid-python** (this documentation) is how you run servers, write Cellucid-ready AnnData files, and embed the viewer in notebooks.
- **cellucid-annotation** is a separate repo for community annotation workflows (publish/vote/audit).
- **cellucid-r** writes the same Cellucid-ready AnnData from R
  (`cellucid_prepare()`) and has its own viewer and server (`cellucid_view()`,
  `cellucid_serve()`) on the same viewer channel; this chapter documents the
  Python ones. See {doc}`../../r_package/index` for R workflows.

## Serving from somewhere that is not your laptop

- {doc}`12_remote_servers_ssh_tunneling_and_cloud` -- lab servers, cloud VMs,
  containers, and the SSH tunnel that covers most of them.
- {doc}`17_hpc_slurm_and_compute_node_serving` -- an HPC cluster, where the
  machine holding the data usually refuses your SSH connection and the tunnel
  has to terminate on a login node instead.

## Related docs

- Web-app-centric loading docs: {doc}`../../web_app/b_data_loading/01_loading_options_overview`,
  with the server ({doc}`../../web_app/b_data_loading/04_server_tutorial`) and
  Jupyter ({doc}`../../web_app/b_data_loading/05_jupyter_tutorial`) walkthroughs
- Preparing an AnnData (`prepare`): {doc}`../c_data_preparation_api/index`
- Hooks (Python ↔ frontend events/commands): {doc}`../e_jupyter_hooks/index`

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
