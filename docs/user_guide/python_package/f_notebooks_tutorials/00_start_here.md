# Start here (how to use these tutorials)

This section is the “slow, safe, and complete” path through Cellucid.

It is intentionally written for **mixed audiences** and uses “layered” writing:

1) **Fast path (wet lab / non-technical / new to Python)** — copy/paste, minimal choices, clear success criteria.
2) **Practical path (computational users)** — reproducible exports, scale/performance, workflow integration.
3) **Deep path (expert / developer)** — event schemas, session bundles, debugging, edge-case behavior.

## What Cellucid is (and what it is not)

**Cellucid** is a **web app** (the interactive viewer UI).

**`cellucid-python`** (this repo) is the **helper package** that:
- prepares datasets into a viewer-friendly export format (`cellucid.prepare(...)`)
- serves datasets over HTTP (`cellucid.serve(...)`, `cellucid.serve_anndata(...)`)
- embeds the web app in notebooks (`cellucid.show(...)`, `cellucid.show_anndata(...)`)
- enables Python↔frontend bidirectional communication (hooks/events and commands)

`cellucid-annotation` is a separate repo for community annotation workflows.

`cellucid-r` provides the R export workflow. These tutorials focus on Python
servers, notebook embedding, hooks, and Python-side preparation.

## At a glance

**Time**
- Minimal: 5–10 minutes (get a viewer visible and responsive)
- Typical tutorial: 20–60 minutes

**Prerequisites**
- Python 3.10+ (recommended)
- A Jupyter environment (JupyterLab, classic notebook, VSCode notebooks, or Colab)
- `pip install cellucid`

Common optional dependencies (depending on tutorial):
- `anndata`, `scanpy`, `numpy`, `scipy`, `pandas`

## Folder map (where everything lives)

All notebook/tutorial docs in this section live here:
- `cellucid-python/docs/user_guide/python_package/f_notebooks_tutorials/`

This folder contains:
- narrative tutorial pages (`.md`)
- runnable notebooks (`prepare_*.ipynb`, `jupyter_*.ipynb`)
- a small supporting dataset folder (`datasets/`)

## The one mental model that explains 90% of issues

When you call `show(...)` or `show_anndata(...)` in a notebook, Cellucid does two things:

1) Starts a **local HTTP server** in Python (serving data + viewer UI assets)
2) Displays an **iframe** in your notebook pointing to that server

If “nothing shows up”, it’s almost always one of these:
- the server didn’t start (exception, port conflict)
- the browser can’t reach the server (remote kernel without proxy/tunnel)
- the exact viewer generation could not be established from its configured
  source

## Install (minimal)

```bash
pip install cellucid
```

```{tip}
If you plan to follow the workflows using `AnnData`, install the common ecosystem packages too:

`pip install anndata scanpy`
```

## Network requirement

At every viewer/server startup, the package fetches the exact web generation
declared by:
- `https://www.cellucid.com`

It downloads every inventory object to a staging directory, verifies the full
file set, MIME types, byte lengths, hashes, and build identity, then publishes
that generation atomically. The local server serves it from the same origin as
the dataset API. Source or verification failure stops startup; a previous
generation is not substituted.

For details (generation directory, source configuration, and firewalls), see:
- {doc}`../../web_app/b_data_loading/05_jupyter_tutorial`

## How to run these tutorials

### Option A — Read as docs and copy/paste snippets

Start with:
- {doc}`index`
- {doc}`01_beginner_notebooks_wet_lab_friendly`

### Option B — Run the actual notebooks (recommended)

From the repo root (`cellucid-python/`):

```bash
jupyter lab docs/user_guide/python_package/f_notebooks_tutorials/
```

```{note}
The docs site does not execute notebooks during the build (`nb_execution_mode = "off"`), so “expected output” sections describe what you should see rather than embedding executed outputs.
```

### Option C — Remote/HPC notebooks

If your kernel runs on a remote machine but your browser is on your laptop, you need a browser-reachable URL for the Python server:
- Jupyter Server Proxy (recommended in JupyterHub/remote Jupyter)
- Colab’s proxy URL (obtain it explicitly)
- SSH port forwarding (classic HPC)

In every remote case, pass the resulting browser-reachable HTTP(S) base as
`client_server_url=`. This is covered in:
- {doc}`22_large_dataset_server_mode_and_lazy_gene_expression`

## What success looks like

## Troubleshooting (quick triage)

### 1) Verify the server is alive

In Python:
```python
print(viewer.server_url)
print(viewer.viewer_url)
```

Then open in a browser tab:
- `http://127.0.0.1:<port>/_cellucid/health`

### 2) Ask Cellucid for a structured debug report

```python
report = viewer.debug_connection()
report
```

This report is designed to answer:
- “Did the server respond?”
- “Is the local web generation exact and complete?”
- “Are frontend→Python events arriving?”

### 3) Most common fixes (in order)

1) **Port conflict** → choose another explicit port, or use `port=0` in Python.
2) **Web generation** → correct the reported source, inventory, object, or
   destination failure.
3) **Remote kernel** → expose the port, then pass its exact base as
   `client_server_url=`.
4) **Corporate network** → allow access to the configured viewer source.

## Next steps

- If you are new: {doc}`01_beginner_notebooks_wet_lab_friendly`
- If you want reproducible exports: {doc}`21_prepare_exports_with_quantization_and_compression`
- If you want hooks/callbacks: {doc}`23_programmatic_highlighting_and_selection_callbacks`
