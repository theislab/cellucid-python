# Viewing APIs (serve / serve_anndata / show / show_anndata + loading options)

Cellucid is a **web app** (WebGL viewer UI). `cellucid-python` gives you **Python + CLI entry points** to open that viewer against *your data*.

There are two big “modes”:

- **Server mode** (terminal/browser): `cellucid serve …`, `serve(…)`, `serve_anndata(…)`
- **Notebook mode** (Jupyter iframe): `show(…)`, `show_anndata(…)`

And two big “data shapes”:

- **Exported directory** (recommended): a folder produced by `cellucid.prepare(…)`
- **AnnData direct** (convenient): `.h5ad`, `.zarr`, or an in-memory `AnnData`

## At a glance

**Audience**
- Wet lab / non-technical: copy/paste the quickstarts and use the decision tree.
- Computational: focus on data format tradeoffs: exported directories,
  read-only-backed `.h5ad`, and eagerly loaded `.zarr`.
- Power user / developer: focus on networking (SSH tunnels), caching, and debug endpoints.

**Time**
- First successful view: ~5 minutes
- Full read (with edge cases): ~30–45 minutes

**Prerequisites**
- `pip install cellucid`
- One of: an export directory, a `.h5ad`, a `.zarr`, or an in-memory `AnnData`

## Where to start

- Not sure what to do → {doc}`03_choose_your_workflow_decision_tree`
- Want a terminal “just show me the viewer” flow → {doc}`04_cli_cellucid_serve_quickstart`
- Want a scripted Python server (non-notebook) → {doc}`05_python_serve_and_serve_anndata_quickstart`
- Want a notebook embed (Jupyter/VSCode/Colab) → {doc}`06_jupyter_show_and_show_anndata_quickstart`
- Something failed → {doc}`15_troubleshooting_viewing`

## Quick start (copy/paste)

### Terminal (recommended “default”)

```bash
pip install cellucid
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1
# For .zarr, use the same required identity flags.
# or: cellucid serve /path/to/export_dir
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

### Exported dataset (fastest + most reproducible)

```python
from cellucid import show

viewer = show("./export_dir", height=600)
```

## Read once: how the exact viewer generation is served

```{important}
By design, the Python server serves **both**:

1) your dataset API (e.g. `/dataset_identity.json`, `/points_3d.bin`, `/var/<index>.values.f32`), and
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
- **cellucid-python** (this documentation) is how you run servers, export data, and embed the viewer in notebooks.
- **cellucid-annotation** is a separate repo for community annotation workflows (publish/vote/audit).
- **cellucid-r** writes export folders that load in the web app. Python owns the
  server, notebook-embedding, and programmatic viewer APIs documented in this
  chapter; see {doc}`../../r_package/index` for R workflows.

## Related docs

- Web-app-centric loading docs: {doc}`../../web_app/b_data_loading/01_loading_options_overview`
- Exporting (`prepare`) docs: {doc}`../c_data_preparation_api/index`
- Hooks (Python ↔ frontend events/commands): {doc}`../e_jupyter_hooks/index`

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
