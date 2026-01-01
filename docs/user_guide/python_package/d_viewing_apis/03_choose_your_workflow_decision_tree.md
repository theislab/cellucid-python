# Choose your workflow (decision tree)

This page helps you choose **one workflow** that fits:
- your environment (terminal vs notebook vs remote),
- your data format (`exported` vs `.h5ad` vs `.zarr` vs in-memory),
- and your goal (quick look vs scaling vs sharing).

If you already know what you want, jump straight to:
- CLI: {doc}`04_cli_cellucid_serve_quickstart`
- Jupyter: {doc}`06_jupyter_show_and_show_anndata_quickstart`

## Fast path (for beginners / wet lab)

If you just want to see your dataset with the least thinking:

1) Install:

```bash
pip install cellucid
```

2) Run (terminal):

```bash
cellucid serve /path/to/data.h5ad
```

3) Open the printed **Viewer URL** in your browser.

That’s it. If you get stuck, go to {doc}`15_troubleshooting_viewing`.

## Decision tree

### Step 1 — Are you in a notebook?

#### Yes: notebook (Jupyter/VSCode/Colab)

1) If you already have an **export directory** (from `prepare`):
   - Use `show(export_dir)` → {doc}`06_jupyter_show_and_show_anndata_quickstart`
2) If you have a `.h5ad` or `.zarr`:
   - Use `show_anndata("data.h5ad")` / `show_anndata("data.zarr")` → {doc}`06_jupyter_show_and_show_anndata_quickstart`
3) If your notebook is served over **HTTPS** (JupyterHub, cloud):
   - Read the proxy section before debugging: {doc}`10_notebook_widget_mode_advanced`

#### No: terminal / script

1) Easiest:
   - `cellucid serve <path>` (auto-detects format) → {doc}`04_cli_cellucid_serve_quickstart`
2) If you want to start a server from Python:
   - `serve(export_dir)` or `serve_anndata(h5ad/zarr/AnnData)` → {doc}`05_python_serve_and_serve_anndata_quickstart`

### Step 2 — What data format do you have?

#### I have an **export directory** (recommended)

Use one of:
- Terminal: `cellucid serve /path/to/export_dir` → {doc}`04_cli_cellucid_serve_quickstart`
- Python: `serve("/path/to/export_dir")` → {doc}`05_python_serve_and_serve_anndata_quickstart`
- Notebook: `show("/path/to/export_dir")` → {doc}`06_jupyter_show_and_show_anndata_quickstart`

You can also load exports directly in the web app with a folder picker:
{doc}`../../web_app/b_data_loading/03_browser_file_picker_tutorial`.

#### I have a `.h5ad`

Prefer:
- Terminal: `cellucid serve /path/to/data.h5ad`
- Notebook: `show_anndata("data.h5ad")`

This uses **backed mode** (lazy loading) by default. Avoid `--no-backed` unless you know why.

#### I have a `.zarr`

Prefer:
- Terminal: `cellucid serve /path/to/data.zarr`
- Notebook: `show_anndata("data.zarr")`

Zarr is naturally chunked and works well for large datasets.

#### I only have an in-memory `AnnData`

Use:
- Notebook: `show_anndata(adata)` (best experience)
- Script: `serve_anndata(adata)` (works, but think about memory)

If your dataset is large, it’s often better to view from a file (`.h5ad` or `.zarr`) so you can use lazy loading.

### Step 3 — Are you on a remote machine (HPC / cloud)?

If yes, do **not** start by binding to `0.0.0.0`.

Recommended:
- run `cellucid serve … --no-browser` on the remote machine,
- create an **SSH tunnel** from your laptop,
- open `http://127.0.0.1:<port>/` locally.

Step-by-step: {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

### Step 4 — Do you need to share with collaborators?

- Public sharing (no server): export once (`prepare`) then use the GitHub workflow (#2):
  {doc}`../../web_app/b_data_loading/02_local_demo_tutorial`
- Private sharing: export once, then share the folder (or serve on an internal machine with an SSH tunnel).

## “I chose a path — what next?”

- CLI workflow: {doc}`04_cli_cellucid_serve_quickstart`
- Python server workflow: {doc}`05_python_serve_and_serve_anndata_quickstart`
- Notebook workflow: {doc}`06_jupyter_show_and_show_anndata_quickstart`
- Export format details + validation: {doc}`07_exported_directory_mode_show_and_serve`
- AnnData requirements (UMAP keys, vector fields, gene IDs): {doc}`08_anndata_mode_show_anndata_and_serve_anndata`
