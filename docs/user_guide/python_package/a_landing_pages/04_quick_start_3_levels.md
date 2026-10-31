# Quick start (3 levels)

**Audience:** everyone (choose your depth)  
**Time:** 5–45+ minutes  
**What you’ll get:** a working Cellucid viewer (in a notebook or browser) + a path to scaling/sharing

This page gives you three depth levels. Stop as soon as you have what you need:

- **Level 1**: see an AnnData immediately (minimal setup)
- **Level 2**: write a Cellucid-ready AnnData (fast, shareable, names itself)
- **Level 3**: treat the viewer as an interactive UI in your notebook (hooks + commands)

```{note}
If you prefer a notebook-style, highly verbose walkthrough (with many edge cases and troubleshooting), start with {doc}`../f_notebooks_tutorials/index`.
```

---

## Level 1 — Fast path: direct AnnData viewing

This is the fastest way to prove “Cellucid works” and to start exploring a dataset.

### Prerequisites

- A notebook that shows widgets: JupyterLab (tested), Jupyter Notebook 7,
  VS Code, Google Colab, marimo or Databricks (expected, not tested by the
  project; see {doc}`03_compatibility_matrix_must_be_explicit`)
- An AnnData object with at least one embedding in `obsm` — the `X_umap` that
  `sc.tl.umap()` writes is read as written

```{tip}
If you don’t want notebooks at all, skip to Level 2 and use the CLI: `cellucid serve ...`.
```

### Step 1: Load (or create) an AnnData

**Option A: load a real dataset**

```python
import anndata as ad

adata = ad.read_h5ad("my_dataset.h5ad")
```

**Option B: create a tiny toy dataset (recommended for a first-run smoke test)**

```python
import numpy as np
import pandas as pd
import anndata as ad

X = np.array([[0, 1],
              [2, 3],
              [4, 5]], dtype=np.float32)  # cells × genes

obs = pd.DataFrame(
    {"cluster": ["A", "A", "B"], "score": [0.1, 0.2, 0.3]},
    index=[f"cell_{i}" for i in range(X.shape[0])],
)
var = pd.DataFrame(index=["G1", "G2"])

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.obsm["X_umap"] = np.array([[0, 0],
                                 [1, 0],
                                 [0, 1]], dtype=np.float32)
```

### Step 2: Confirm you have an embedding

```python
list(adata.obsm.keys())
```

**Coming from Scanpy?** `sc.tl.umap(adata)` writes `adata.obsm["X_umap"]`, and
Cellucid reads it at the dimension its own column count states — two columns as
2D, three as 3D — with nothing to rename. With no `X_umap` at all, Cellucid
picks another `obsm` array by a deterministic best guess (t-SNE before PCA, for
example) and says which one; pass `embedding_keys={"2d": "X_tsne"}` to choose
it yourself.

### Step 3: Show the viewer

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    height=600,
    dataset_name="Small example",
    dataset_id="small-example",
)
```

### Step 4 (optional): Confirm hooks are working

```python
@viewer.on_ready
def _ready(ev):
    print("Viewer ready:", ev)
```

### Step 5 (optional): Programmatic control

```python
viewer.set_color_by("cluster")
viewer.highlight_cells([0, 2], color="#00cc66")
```

### What this does (mental model)

`show_anndata(...)`:
- shows a notebook widget whose output holds the Cellucid web app, which the
  browser loads from `https://www.cellucid.com`,
- hands the data to the viewer through the notebook's own connection; the
  viewer reads it itself, with range requests (an in-memory AnnData as a
  store in priority order made as the viewer reads it, so a subsample reads
  only its cells; a MuData the same way),
- carries commands to the viewer and UI events back to Python on the same
  connection. No port is opened and no server is started.

### What success looks like

- An embedded viewer appears in the notebook output cell.
- You see points (cells) in the embedding.
- Coloring by `cluster` changes point colors.

### Common pitfalls (Level 1)

- **Not in a notebook**: `show_anndata(...)` raises and names the browser-tab
  alternatives (`AnnDataViewer`, `serve_anndata`, `cellucid serve`).
- **Nothing to plot**: the error lists every `obsm` entry and why it was not
  used; compute an embedding or name one with `embedding_keys`.
- **The viewer appears only after the cell ends**: the kernel handles the
  widget's messages between cells, so put commands and hooks in later cells
  (commands sent earlier are delivered once it connects).
- **Which notebooks**: see {doc}`03_compatibility_matrix_must_be_explicit`.
- **Viewer source unavailable**: viewer-serving startup stops; see
  {doc}`../installation`.
- **A gene that will not colour**: `show_anndata(...)` starts fine on an object
  whose `adata.X` holds infinities, because nothing reads a gene column until
  you select one. Selecting it shows a notification titled
  `Gene cannot be shown` that names the gene, counts the offending values, and
  lists the first offending cells (`NaN` is drawn, as the neutral grey for "no
  value here"). Everything else in the session keeps working. Repair the
  matrix and reload — see
  {doc}`../d_viewing_apis/15_troubleshooting_viewing`.

---

## Level 2 — Practical path: `prepare(...)`, then open the file

This is the “real workflow” for large datasets and collaboration: you write a
**Cellucid-ready AnnData** — a standard AnnData laid out for fast browser
reads, naming itself in `uns['cellucid']` — that you can archive, share, or
host.

### Inputs

`prepare(...)` requires an output path, its `format`, a `dataset_id` and a
`dataset_name`, and at least one embedding: pass an AnnData (or a MuData) as
`data=`, which fills every component it holds, or pass the components one by
one, starting with `embeddings={...}`. An argument given explicitly replaces
that one component of `data`. The identity may be omitted when `data` declares
it in `uns['cellucid']`. Everything else is optional:
- `obs` (cell metadata)
- `var` + `gene_expression` (gene colouring and search)
- `connectivities` (the neighbour graph)
- `vector_fields` (velocity / drift overlays)
- a latent space for outlier measurement: pass it as one of the `embeddings`
  and name it with `latent_key`

```{warning}
**Row order is identity.** Every array you pass must be aligned so that row
`i` always refers to the same cell across embeddings, `obs`, expression,
connectivities, and vector fields. Take them all from one AnnData.
```

### Step 1: Prepare

```python
from cellucid import prepare

path = prepare(
    "pbmc_demo.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"], "X_pca": adata.obsm["X_pca"]},
    dataset_id="pbmc-demo",
    dataset_name="PBMC demo",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    latent_key="X_pca",
    force=True,  # replace an earlier output at the same path
)
```

`format` is required: `"zarr"` (a folder of small files) for GitHub, which
refuses files over 100 MB; `"h5ad"` (or `"zarr.zip"`, the Zarr store in one
file) for everything else — single-file hosts such as Zenodo, Figshare,
Hugging Face, S3 or a lab server, local files, `cellucid serve` and Jupyter
({ref}`which-format`). The path must end in the matching `.h5ad`, `.zarr` or
`.zarr.zip`. A MuData passed as `data=` is written as `format="h5mu"` (one
`.h5mu` file), or as MuData Zarr with `"zarr"` / `"zarr.zip"`; see
{doc}`../c_data_preparation_api/12_prepare_from_anndata_or_mudata`.

Or let the AnnData fill the components (every numeric `obsm` array is then
written, not only the two named above):

```python
path = prepare(
    "pbmc_demo.h5ad",
    format="h5ad",
    data=adata,  # obs, X with var, every numeric obsm array, obsp['connectivities']
    dataset_id="pbmc-demo",
    dataset_name="PBMC demo",
    latent_key="X_pca",
    force=True,
)
```

### Step 2: Open it (pick one)

**Option A: embed in a notebook**

```python
from cellucid import show_anndata

viewer = show_anndata(path, height=700)
```

**Option B: open in a browser via CLI (recommended for non-notebook users)**

`cellucid serve pbmc_demo.h5ad` — the file names itself, so no identity flags
are needed.

**Option C: open the file in the hosted web app**

If you prefer not to run a local server, open the web app and load the file
with **Open…**, or host it and paste its address into the **Link or server**
bar. For
click-by-click instructions, see {doc}`../../web_app/b_data_loading/index`.

### What you get on disk

One standard AnnData: `anndata.read_h5ad("pbmc_demo.h5ad")` reads it. Its
expression matrix is gene-major, its arrays are chunked for range reads, and
`uns['cellucid']` holds the id, name and any launch defaults you set. The full
contract: {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`.

### What success looks like

- The viewer opens the dataset under the name you gave it, and you can:
  - color by `obs` fields,
  - (optionally) search genes and color by expression.

### Common pitfalls (Level 2)

- **Wrong shapes**: `obs` length, embedding row counts, expression matrix shape must all match.
- **Expression orientation**: `gene_expression` must be `(n_cells, n_genes)` with `len(var) == n_genes`.
- **Existing output**: writing to the same path again needs `force=True`, or pick a fresh path.

Every refusal and its fix: {doc}`../c_data_preparation_api/11_troubleshooting_prepare`.

---

## Level 3 — Deep path: custom server + hooks + message handling

This level is for notebook users who want the viewer as a live UI:
- select cells in the browser → run Python analysis on that subset,
- compute results → send highlights/colors back to the viewer.

### Step 1: Create a viewer

A prepared file opens as it is (`show_anndata(path, height=700)`); any other
AnnData takes its identity:

```python
from cellucid import show_anndata
viewer = show_anndata(
    adata,
    dataset_name="PBMC demo",
    dataset_id="pbmc-demo",
)
```

### Step 2: Register hooks (frontend → Python)

```python
@viewer.on_ready
def _on_ready(ev):
    print("ready:", ev)

@viewer.on_selection
def _on_selection(ev):
    cells = ev.get("cells", [])
    print("selected:", len(cells), "cells")
```

### Step 3: Send commands (Python → frontend)

```python
# highlight first 10 selected cells in green
@viewer.on_selection
def _highlight(ev):
    cells = (ev.get("cells") or [])[:10]
    viewer.highlight_cells(cells, color="#00cc66")
```

### Step 4: Make callbacks robust (recommended patterns)

- Always guard against empty selections.
- Let analysis failures surface, or catch only errors you can report with
  explicit context; do not swallow callback exceptions. A hook runs on the
  viewer channel's thread, not in your cell: an exception it raises is logged
  (logger `cellucid.viewer_channel`) and the channel stays open, so record
  failures where your notebook will see them.
- If your analysis is slow, consider:
  - running it asynchronously,
  - debouncing rapid selection events.

### Step 5: Debug connectivity when things feel “haunted”

Cellucid provides a structured report:

```python
report = viewer.debug_connection()
report
```

Look for:
- `connected_outputs` (is the widget shown in an output that reached the kernel?)
- `channel_connected` (did the viewer open its channel?)
- `state.ready` (did the dataset open?)
- `recent_events` (what arrived, and the commands the viewer refused)
- `web_source_url` (where the browser loads the app from)

A browser-tab viewer (`AnnDataViewer`) reports its server instead
(`server_health`, `anndata_store_probes`, `frontend_roundtrip`).

### Step 6: Clean up

```python
viewer.stop()
```

If you created multiple viewers:

```python
from cellucid.jupyter import cleanup_all
cleanup_all()
```

### What success looks like

- Selecting cells triggers your Python callbacks.
- Python can send a highlight that visibly updates the viewer.

---

## Edge cases (read this before scaling up)

- **Huge datasets**: the viewer transfers only what a view needs; for a dataset
  you will share or host, write it once with `prepare()` so every gene is one
  small read.
- **Missing embeddings**: no embedding → no viewer; compute one first.
- **NaN/Inf values**: `prepare(...)` refuses an expression matrix or graph
  holding NaN or an infinity before writing anything, with a count of what it
  found and the first affected positions (`obs` columns may hold NaN, as a
  missing value, but not infinities). A row of NaN in an embedding or vector
  field is a cell that array does not place: written as given, left out of
  that embedding's views (counted and named in a short notice), kept
  everywhere else. The viewer reads the
  values itself and refuses one field at a time instead: a gene or continuous
  column holding an infinity is refused with a counted diagnosis, and the rest
  of the session keeps working. In the viewer `NaN` and infinity are not alike
  — `NaN` is drawable and renders as the neutral grey meaning "not measured
  here", an infinity has no position on a colour scale.
- **Category explosions**: categorical fields with very high cardinality can become unusable in the UI; consider grouping or filtering.
- **Remote kernels**: nothing to forward for a notebook widget; a browser-tab
  viewer needs an SSH tunnel or proxy (see
  {doc}`03_compatibility_matrix_must_be_explicit`).
- **Viewer-source access**: the browser loads the app from
  `https://www.cellucid.com` (a browser-tab server downloads and verifies it
  at startup; see {doc}`../installation`).

## Troubleshooting (quick start)

Organize by symptom. Each entry points to the most likely root causes.

### Symptom: “Nothing shows up in my notebook”

**Likely causes**
- The notebook front end does not show widgets (`anywidget` missing where the
  notebook server runs); the output then says so.
- The browser cannot load the app from `https://www.cellucid.com`; the output
  says it did not load.

**How to confirm**
- Run, in a new cell:
  ```python
  viewer.debug_connection()
  ```

**Fix**
- `pip install anywidget` where the notebook server runs, and reload the page.
- Allow `https://www.cellucid.com`, or host a copy and pass `web_source_url=`.

---

### Symptom: `prepare(...)` succeeds, but the viewer can’t find fields / genes

**Likely causes**
- You didn’t pass `gene_expression` + `var`.
- `obs_keys` lists only some `obs` columns, so the viewer lists only those.

**How to confirm**
- `anndata.read_h5ad(path)` shows what the file holds; `uns['cellucid']['defaults']`
  shows `obsKeys` when it was set.

**Fix**
- Prepare again with `var=adata.var` and `gene_expression=adata.X` (and
  `force=True`), or change `obs_keys`.

---

### Symptom: hooks don’t fire / selection doesn’t reach Python

**Likely causes**
- Viewer is not fully loaded yet (wait for `on_ready`).
- You have multiple viewers and are listening on the wrong one.
- The cell that registered the hook is still running: the widget's events
  arrive between cells.

**Fix**
- Add an `on_message` logger, call `viewer.debug_connection()`, and stop old viewers (`viewer.stop()`).

## Next steps

- Detailed `prepare()` inputs/edge cases: {doc}`../c_data_preparation_api/index`
- All viewing modes (CLI/server/notebook): {doc}`../d_viewing_apis/index`
- Hooks system deep dive: {doc}`../e_jupyter_hooks/index`
- Notebook-style tutorials: {doc}`../f_notebooks_tutorials/index`
- Symptom-based troubleshooting index: {doc}`../i_troubleshooting_index/index`
