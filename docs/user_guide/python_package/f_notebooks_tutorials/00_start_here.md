# Start here (how to use these tutorials)

This section is the “slow, safe, and complete” path through Cellucid.

It is intentionally written for **mixed audiences** and uses “layered” writing:

1) **Fast path (wet lab / non-technical / new to Python)** — copy/paste, minimal choices, clear success criteria.
2) **Practical path (computational users)** — reproducible datasets, scale/performance, workflow integration.
3) **Deep path (expert / developer)** — event schemas, session bundles, debugging, edge-case behavior.

## What Cellucid is (and what it is not)

**Cellucid** is a **web app** (the interactive viewer UI).

**`cellucid-python`** (this repo) is the **helper package** that:
- writes Cellucid-ready AnnData files (`cellucid.prepare(...)`)
- serves AnnData, MuData and the other formats the web app reads (Loom, 10x,
  MatrixMarket folders, R objects) over HTTP (`cellucid serve`,
  `cellucid.serve_anndata(...)`)
- shows the web app in notebooks as a widget (`cellucid.show_anndata(...)`)
- enables Python↔frontend bidirectional communication (hooks/events and commands)
- brings what you make in the viewer back onto your AnnData in one call
  (`viewer.pull(adata)`, `cellucid.pull(adata, exported_files)`) and pushes
  columns and gene sets into it (`viewer.push(...)`)

`cellucid-annotation` is a separate repo for community annotation workflows.

`cellucid-r` does the same from R (`cellucid_view()`, `cellucid_pull()`,
`cellucid_prepare()`); {doc}`06_cellucid_in_your_pipeline_python_and_r` pairs
every round-trip notebook here with its R twin.

## At a glance

**Time**
- Minimal: 5–10 minutes (get a viewer visible and responsive)
- Typical tutorial: 20–60 minutes

**Prerequisites**
- Python 3.11–3.14
- A notebook that shows widgets: JupyterLab, Jupyter Notebook 7, VS Code,
  Google Colab, marimo or Databricks (JupyterLab is tested by the project; see
  {doc}`../e_jupyter_hooks/03_supported_environments_matrix`)
- `pip install cellucid`

Common optional dependencies (depending on tutorial):
- `anndata`, `scanpy`, `numpy`, `scipy`, `pandas`

## Folder map (where everything lives)

All notebook/tutorial docs in this section live here:
- `cellucid-python/docs/user_guide/python_package/f_notebooks_tutorials/`

This folder contains:
- narrative tutorial pages (`.md`)
- {doc}`prepare_pancreas`, which writes and checks a compact Cellucid-ready AnnData
- {doc}`jupyter_pancreas`, which exercises embedding, hooks, and sessions
- a checked-in real Pancreas H5AD in `datasets/`

## The one mental model that explains 90% of issues

When you call `show_anndata(...)` in a notebook, Cellucid:

1) checks the dataset (an in-memory AnnData is served in priority order,
   made as the viewer reads it), and
2) shows a notebook widget in the cell's output. Your browser loads the
   Cellucid app from `https://www.cellucid.com` inside it, and the app reads
   the dataset from the kernel through the notebook's own connection. No port
   is opened and no server address is involved.

The notebook hands the widget's messages to the kernel **between cells**: the
viewer appears when the cell finishes, what you make in the viewer reaches the
kernel while it is idle, and the next cell reads it. If “nothing shows up”,
it is almost always one of these:
- the notebook cannot show widgets (JupyterLab or Notebook 7 without
  `anywidget` where the notebook server runs),
- the browser cannot load the app from www.cellucid.com (a blocking network:
  host a copy and pass `web_source_url=`),
- the code is not running in a notebook (a script: use
  `cellucid.AnnDataViewer(...)` for a browser tab).

## Install (minimal)

```bash
pip install cellucid
```

```{tip}
If you plan to follow the workflows using `AnnData`, install the common ecosystem packages too:

`pip install anndata scanpy`
```

## Network requirement

In a notebook, the **browser** loads the Cellucid app from
`https://www.cellucid.com` (or the copy you name with `web_source_url=`); the
kernel needs no internet, and your data never leaves it except to that page.
A browser-tab viewer (`AnnDataViewer`, `serve_anndata`, `cellucid serve`)
instead serves a verified copy of the app from this machine; see
{doc}`../../web_app/b_data_loading/05_jupyter_tutorial`.

## How to run these tutorials

### Option A — Read as docs and copy/paste snippets

Start with:
- {doc}`index`
- {doc}`01_beginner_notebooks_wet_lab_friendly`

### Option B — Run the actual Pancreas notebooks (recommended)

From the repo root (`cellucid-python/`):

```bash
python -m pip install -e . jupyterlab
jupyter lab docs/user_guide/python_package/f_notebooks_tutorials/jupyter_pancreas.ipynb
```

```{note}
The notebook resolves its checked-in input from any working directory inside
the clone. You do not need to edit a data path.
```

### Option C — Google Colab, VS Code, remote notebooks

Every round-trip notebook ({doc}`06_cellucid_in_your_pipeline_python_and_r`)
has an "Open in Colab" badge, and runs in VS Code or on a remote JupyterHub
as it is: the widget needs no port, proxy or tunnel. Outside a clone of this
repository the notebooks download the Pancreas file once.

## What success looks like

- {doc}`prepare_pancreas` writes the Cellucid-ready AnnData, reads it back with
  its `uns['cellucid']` identity and gene-major matrix, opens it, then removes
  the temporary file cleanly.
- {doc}`jupyter_pancreas` shows the viewer in the output, renders the real
  `clusters` field, pulls a viewer highlight back into AnnData with
  `viewer.pull()`, and its connection report says the channel is connected
  and the viewer opened all 2,531 cells.
- The notebooks of {doc}`06_cellucid_in_your_pipeline_python_and_r` bring
  highlights, labels, scores and gene sets back in one call, and push
  clusters and gene sets in.

## Troubleshooting (quick triage)

### 1) Ask Cellucid for a structured report

In a new cell:

```python
report = viewer.debug_connection()
report
```

It waits for nothing and answers:
- “Did the notebook draw the widget?” (`connected_outputs`)
- “Did the app load and open its channel?” (`channel_connected`,
  `web_source_url`)
- “What did the viewer last say, and which commands did it refuse?”
  (`state`, `recent_events`)

### 2) Most common fixes (in order)

1) **No output at all** → install `anywidget` in the notebook server's
   environment and reload the page.
2) **An output that stays empty** → the browser cannot load
   www.cellucid.com; host a copy and pass `web_source_url=`.
3) **A call raises "waits for the viewer"** → run it in a later cell; the
   viewer is heard between cells ({doc}`06_cellucid_in_your_pipeline_python_and_r`).
4) **A plain script** → `cellucid.AnnDataViewer(...)` opens the viewer in a
   browser tab instead.

## Next steps

- If you are new: {doc}`01_beginner_notebooks_wet_lab_friendly`
- If you want a shareable Cellucid-ready AnnData: {doc}`21_prepare_a_cellucid_ready_anndata`
- If you want hooks/callbacks: {doc}`23_programmatic_highlighting_and_selection_callbacks`
