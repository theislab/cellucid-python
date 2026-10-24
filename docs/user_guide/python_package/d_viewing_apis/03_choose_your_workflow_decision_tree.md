# Choose your workflow (decision tree)

This page helps you choose **one workflow** that fits:
- your environment (terminal vs notebook vs remote),
- your data format (`.h5ad` vs Zarr store vs in-memory),
- and your goal (quick look vs scaling vs sharing).

If you already know what you want, jump straight to:
- CLI: {doc}`04_cli_cellucid_serve_quickstart`
- Jupyter: {doc}`06_jupyter_show_anndata_quickstart`

## Fast path (for beginners / wet lab)

If you just want to see your dataset with the least thinking:

1) Install:

```bash
pip install cellucid
```

2) Run (terminal):

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
```

3) Open the printed **Viewer URL** in your browser.

That’s it. If you get stuck, go to {doc}`15_troubleshooting_viewing`.

## Decision tree

### Step 1 — Are you in a notebook?

#### Yes: notebook (Jupyter/VSCode/Colab)

1) If you have a file written by `prepare()`:
   - Use `show_anndata("pbmc_cellucid.h5ad")`; it names itself →
     {doc}`06_jupyter_show_anndata_quickstart`
2) If you have any other `.h5ad` or `.zarr`:
   - Use `show_anndata("data.h5ad", dataset_name="My dataset", dataset_id="my-dataset")`
     or the same exact identity arguments with `"data.zarr"` →
     {doc}`06_jupyter_show_anndata_quickstart`
3) If your notebook is served over **HTTPS** (JupyterHub, cloud):
   - Read the proxy section before debugging: {doc}`10_notebook_widget_mode_advanced`

#### No: terminal / script

1) Easiest:
   - `cellucid serve <path>` (auto-detects format) → {doc}`04_cli_cellucid_serve_quickstart`
2) If you want to start a server from Python:
   - `serve_anndata(data, dataset_name=..., dataset_id=...)` →
     {doc}`05_python_serve_anndata_quickstart`

### Step 2 — What data format do you have?

#### I have a file written by `prepare()` (a Cellucid-ready AnnData)

It declares its own name and id, so open it as it is:
- Terminal: `cellucid serve /path/to/pbmc_cellucid.h5ad`
- Notebook: `show_anndata("pbmc_cellucid.h5ad")`
- Browser, no Python: the web app's **Open…** button
  ({doc}`../../web_app/b_data_loading/03_browser_file_picker_tutorial`).

#### I have a `.h5ad`

Prefer:
- Terminal: `cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`
- Notebook: `show_anndata("data.h5ad", dataset_name="My dataset", dataset_id="my-dataset")`

The viewer reads the file in place, with range requests; nothing is converted first.

#### I have a Zarr store (folder or `.zip`)

Prefer:
- Terminal: `cellucid serve /path/to/data.zarr --dataset-name "My dataset" --dataset-id my-dataset`
- Notebook: `show_anndata("data.zarr", dataset_name="My dataset", dataset_id="my-dataset")`

The viewer reads the store in place, key by key (a zipped store, as one file).

#### I only have an in-memory `AnnData`

Use:
- Notebook: `show_anndata(adata, dataset_name="My dataset", dataset_id="my-dataset")`
- Script: `serve_anndata(adata, dataset_name="My dataset", dataset_id="my-dataset")`

The object is written once to a temporary `.h5ad`, removed when the server
stops. If it is already on disk, pass the path instead, or write a
Cellucid-ready AnnData with `prepare()`.

### Step 3 — Are you on a remote machine (HPC / cloud)?

If yes, do **not** start by binding to `0.0.0.0`.

Recommended:
- run `cellucid serve … --no-browser` on the remote machine,
- create an **SSH tunnel** from your laptop,
- open the exact **Viewer URL** printed by the server through the local tunnel
  (it ends in `/?anndata=true`).

Step-by-step: {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

### Step 4 — Do you need to share with collaborators?

- Public sharing (no server): write a Cellucid-ready AnnData once (`prepare`),
  publish it, and share its `?data=` link:
  {doc}`../../web_app/b_data_loading/11_custom_dataset_repository`
- Private sharing: send the file, or serve it on an internal machine with an
  SSH tunnel.

## “I chose a path — what next?”

- CLI workflow: {doc}`04_cli_cellucid_serve_quickstart`
- Python server workflow: {doc}`05_python_serve_anndata_quickstart`
- Notebook workflow: {doc}`06_jupyter_show_anndata_quickstart`
- Preparing a file: {doc}`../c_data_preparation_api/index`
- AnnData requirements (UMAP keys, vector fields, gene IDs): {doc}`08_anndata_mode_show_anndata_and_serve_anndata`
- Public catalog reference repository:
  {doc}`../../web_app/b_data_loading/11_custom_dataset_repository`
