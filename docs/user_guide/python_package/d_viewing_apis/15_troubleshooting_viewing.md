# Troubleshooting (viewing)

This page is symptom-driven troubleshooting for:
- `cellucid serve …` (CLI server mode),
- `serve_anndata(...)` (Python server mode),
- `show_anndata(...)` (a notebook widget), and `AnnDataViewer(...)` (a browser-tab viewer).

If you’re new, skim the “Quick triage” section first.

---

## Quick triage (do these in order)

### 1) What workflow are you using?

- Terminal/browser → you’re in **server mode** (`cellucid serve …`, `serve_anndata`, `AnnDataViewer`)
- Viewer inside a notebook cell → you’re in **notebook mode** (`show_anndata`)

This matters because notebook mode has no server: steps 2 and 3 do not apply
to it, and its failure modes are the notebook's (widget support, the browser
loading the app, messages arriving between cells).

### 2) Is the server alive? (server mode)

Open (replace `<port>`):

```text
http://127.0.0.1:<port>/_cellucid/health
```

If this doesn’t load, the viewer can’t load either.

### 3) Is the exact viewer generation being served? (server mode)

Open the exact **Viewer URL** printed by the server:

```text
http://127.0.0.1:<port>/?anndata=true
```

A wildcard bind prints the loopback origin plus the machine's own name; see
{doc}`12_remote_servers_ssh_tunneling_and_cloud`.

If this fails, the server may have stopped.

### 4) Did one field fail, or did everything fail?

A viewer that loads, draws points, and refuses exactly one gene or one
continuous `obs` column is not a broken server. The viewer checks each
continuous column as it reads it, and a column holding an infinity is refused
with a notification titled `Gene cannot be shown` that counts the offending
values and names the first affected cells. See
[a gene or column that cannot be drawn](#symptom-a-gene-or-column-is-found-but-cannot-be-drawn).

### 5) If you are in a notebook: run a structured debug report (in a new cell)

```python
viewer.debug_connection()
```

---

## Symptom: “Port already in use” / server won’t start

### Likely causes (ordered)
- another Cellucid server is still running
- another program is using the port (often `8765`)

### How to confirm
- `cellucid serve` ends with `Error: Port 8765 is already in use on 127.0.0.1,
  …`; `serve_anndata(...)` / `AnnDataServer.start()` raise an address-in-use
  `OSError`

### Fix
- choose a port explicitly: `cellucid serve ... --port 9000`
- in Python, pass `port=0` to request one operating-system-assigned port and
  then use the actual `server.url` / `server.viewer_url`
- an `AnnDataViewer` created without `port=` already takes the next free port
  after 8765 and prints a note; one created with an explicit `port=` raises
  instead (a notebook widget, `show_anndata`, uses no port)
- in notebooks, call `viewer.stop()` on old viewers (see {doc}`11_viewer_lifecycle_cleanup_ports_and_multiple_viewers`)

### Prevention
- stop servers when done
- use fixed ports only when you need them (SSH tunnels)

---

## Symptom: viewer startup cannot establish the web generation

This failure occurs before the viewer server is published.

### Likely causes (ordered)
- your environment cannot reach the configured `web_source_url`
- the inventory or a declared object has an invalid MIME type, length, hash,
  file set, or build identity
- the selected `web_cache_dir` is not writable

### How to confirm
- read the exception raised by `cellucid serve`, `serve_anndata(...)` or
  `AnnDataViewer(...)` (a notebook widget keeps no web build: the browser loads
  the app)
- exercise the same establishment contract directly with
  `cellucid.web_cache.ensure_web_ui_cached(force=True)` (or
  `viewer.ensure_web_ui_cached()`)

### Fix
- allow access to the configured source and correct the exact response failure
- pass a writable generation location:

```bash
cellucid serve /path/to/data --web-cache-dir /path/to/cache
```

### Prevention
- verify source access and destination permissions before a workshop or job
- do not plan on a previous local generation being substituted when the source
  is unavailable

---

## Symptom: Viewer loads but dataset is blank / stuck loading

### Likely causes (ordered)
- no plottable embedding in `obsm`, or a shape mismatch
- a Zarr store folder whose consolidated metadata is stale (written before
  later edits to the store)

### How to confirm
- open `/_cellucid/datasets` and see what it lists
- check `adata.obsm.keys()` and shapes; `/_cellucid/datasets` names the
  store's route and the launch choices the server passed

### Fix
- confirm which `obsm` array each slot shows (**Session › Data mapping** names
  it), or name one with `embedding_keys` (see
  {doc}`08_anndata_mode_show_anndata_and_serve_anndata`)

---

## Symptom: the viewer shows a different embedding than I expected

### What decides it
- a slot you name with `embedding_keys` (Python) or `--embedding-<n>d` (CLI)
  shows that array, at the 0-based columns you give
- every other slot follows the default rule: `X_umap_1d` / `X_umap_2d` /
  `X_umap_3d` fill their own slots, and with none of them present a bare
  `X_umap` of 1, 2, or 3 columns fills the slot of its own width
- only when neither fills any slot is a deterministic best guess made from the
  other `obsm` arrays (UMAP-like names first, latent spaces such as PCA last),
  and a wide array is then shown in 2D at columns 0 and 1
- nothing in your object is renamed, assigned, or written back — resolution is a
  read, and the server only reads the store's bytes

### How to confirm
- **Session › Data mapping** names the array each slot shows, and marks a
  guessed one
- a file written by `prepare()`: `adata.uns["cellucid"]["defaults"]` holds the
  `embedding1d`/`embedding2d`/`embedding3d` it opens on

### Fix
- name the array you want: `embedding_keys={"2d": "X_umap_harmony"}` or
  `--embedding-2d X_umap_harmony`; `embedding_keys={"2d": ("X_pca", [0, 1])}` or
  `--embedding-2d X_pca:0,1` for chosen columns of a wide array
- or write the choice into the file once: `prepare(..., embedding_defaults=...)`

Full rules: {doc}`08_anndata_mode_show_anndata_and_serve_anndata`.

---

## Symptom: “Cellucid found nothing to plot”

### Likely causes (ordered)
- no embedding was ever computed, so `obsm` holds only tables, sparse arrays, or
  nothing at all
- every candidate array holds a non-finite value (or no spread) in the columns
  that would be plotted
- arrays do not have one row per cell (`adata.n_obs`)

### How to confirm

The message lists every `obsm` entry and says why it was not used:

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell (3,696) and at least one numeric column. obsm holds: X_df (table), spatial (3,696 × 2, NaN at cell 17). Add a 2-D embedding to obsm (for example with sc.tl.umap in Scanpy), then open the data again.
```

### Fix
- compute an embedding (`sc.tl.umap(adata)`), or repair the array the message
  names (the first non-finite cell is given)
- or name a usable array explicitly with `embedding_keys` / `--embedding-<n>d`
- a dimensional key whose array does not have exactly its own column count is an
  error, not a reinterpretation: `X_umap_3d` holding `(n_obs, 2)` is refused
  rather than shown as 2D
- serving a `.h5ad` or Zarr path? Re-save the file after a repair

---

## Symptom: “Observation field '…' has unsupported dtype”

### Likely causes (ordered)
- `adata.obs` carries a `datetime64`, timedelta, or period column; Cellucid
  shows categorical, boolean, numeric, and string columns only

### How to confirm

```python
print(adata.obs.dtypes)
```

### Fix
- convert the column: `adata.obs["collection_date"] = adata.obs["collection_date"].astype(str)`
- or leave it out by naming the columns to show:

```python
show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
    obs_keys=["cell_type", "leiden", "total_counts"],
)
```

- from the CLI, repeat `--obs-key` once per column

---

## Symptom: Gene search returns nothing / “Gene 'X' not found”

### Likely causes (ordered)
- your gene IDs are not in `var.index` (default)
- wrong `gene_id_column`
- the name is repeated in the chosen gene axis and is shown as `NAME-1`

### How to confirm

```python
print(adata.var.index[:5])
print(adata.var.columns)
```

### Fix
- pass `gene_id_column="..."` to `show_anndata(...)` / `serve_anndata(...)`
- ensure the chosen gene ID field is unique

This symptom is about a gene the viewer never found. A gene it *did* find and
then refused is the next section.

---

## Symptom: a gene or column is found but cannot be drawn

The gene is in `var.index`, the search offers it, you select it — and the viewer
answers with a notification titled `Gene cannot be shown`. Everything else keeps
working.

### What happened

A colour scale needs a finite minimum and maximum. `NaN` is drawn as the
neutral grey meaning “no value here”, but one `+Inf` makes a field's range
infinite and collapses every other cell onto one colour, so the viewer refuses a
gene or continuous `obs` column holding an infinity. It checks each column as it
reads it, one column at a time, as you select it.

### How to confirm

The notification is the diagnosis. It names the gene or field, counts the
offending values by kind out of all cells, lists the first affected cell
indices, and gives a repair one-liner:

```text
Gene "A2ML1" contains 12 infinite values, so it has no colour scale.
Of 18,142,044 cells: 12 +Infinity.
First affected cells: 4,127, 88,301, 250,994, 1,004,112, 7,730,015, ….
NaN is drawn as missing; an infinity is not a value on any scale.
Repair the expression matrix before serving or preparing it — for example `adata.X.data[~np.isfinite(adata.X.data)] = 0` for a sparse matrix, or `np.nan_to_num(adata.X, copy=False)` for a dense one — then reload.
```

Read the counts before acting: twelve cells out of eighteen million and every
cell in the column want opposite responses. A column whose every value is `NaN`
is reported separately (`has no value in any of … cells`); it usually means the
column was never computed.

### Fix

Repair the values in the object, then reload:

```python
import numpy as np

adata.X.data[~np.isfinite(adata.X.data)] = 0          # sparse expression
np.nan_to_num(adata.X, copy=False)                    # dense expression
adata.obs["score"] = adata.obs["score"].replace([np.inf, -np.inf], np.nan).fillna(0)
```

Serving a `.h5ad` or Zarr path rather than an in-memory object? Re-save after
the repair.

:::{note}
`prepare()` refuses the same values with a counted diagnosis. Running it once
is also how to find every affected value at once — it scans the whole matrix,
while the viewer refuses one column at a time as you select it.
:::

---

## Symptom (notebook): the output is empty, or says the app did not load

### Likely causes (ordered)
- the notebook front end does not show widgets (`anywidget` missing where the
  notebook *server* runs)
- the browser cannot load the app from `https://www.cellucid.com` (offline,
  blocked address); the output says “The page at … has not answered as a
  Cellucid viewer”

### Fix
- `pip install anywidget` where the notebook server runs, and reload the page
- allow the address, or host a copy of the app and pass `web_source_url=`

Deep dive: {doc}`10_notebook_widget_mode_advanced`.

---

## Symptom (notebook): the viewer loads, but hooks and commands do nothing

### Likely causes (ordered)
- a cell is still running: the widget's messages reach Python between cells
- the output belongs to a stopped or older viewer (after `viewer.stop()` or a
  kernel restart it keeps its last view and says so)
- in a browser tab (`AnnDataViewer`) behind a proxy: the channel (one
  WebSocket to `/_cellucid/ws`) cannot open because the proxy does not pass
  WebSockets; the viewer says “… this viewer cannot open its channel (…)”

### How to confirm
- `viewer.debug_connection()` in a new cell reports `channel_connected` and
  `recent_events`

### Fix
- let the cell finish; re-run the cell that creates the viewer for a stale output
- behind a proxy, configure it to pass WebSocket upgrades for `/_cellucid/ws`
- commands sent meanwhile wait in Python and are delivered, in order, once the
  viewer connects

---

## Symptom: Vector field / velocity overlay is missing

### Likely causes (ordered)
- vectors not present (no `<field>_<basis>` key beside an `X_<basis>` in `obsm`)
- naming mismatch: the key names no embedding the object holds (e.g.
  `velocity_umap_2d` beside a bare `X_umap` with no `X_umap_2d` — this no
  longer binds; scVelo's `velocity_umap` does)
- width mismatch: the vector does not have as many columns as its embedding
- the view shows another embedding (you have vectors for `X_umap` but are
  viewing `X_umap_3d`)

### How to confirm
- inspect `adata.obsm.keys()` for a `<field>_<basis>` key beside each
  `X_<basis>` (`velocity_umap` beside `X_umap`) and compare their shapes
- the view must show the embedding the field is named for: a vector is drawn
  only over its own embedding

### Fix
- rename keys to scVelo's `<field>_<basis>` convention (`velocity_umap_2d` →
  `velocity_umap` beside `X_umap`), or bind the array under
  **Session › Data mapping › Vector fields**
- ensure you have vectors for the embedding you’re viewing

---

## Symptom: the neighbor-graph overlay is missing

### Likely causes (ordered)
- the store has no `obsp['connectivities']`, and no `knn_key` named another
  graph (a file written by `prepare()` without `connectivities=` has none)
- the named graph is not symmetric (such as `distances`), and the viewer
  refused it when it read it

### How to confirm
- `print(list(adata.obsp.keys()))`; `/_cellucid/datasets` shows `"knn"` under
  `defaults` when the server was given one

### Fix
- in Python, name the graph with `knn_key`:

```python
viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    knn_key="connectivities",
)
```

- from the terminal, add `--knn connectivities`:

```bash
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1 \
  --knn connectivities
```

- or choose a graph in the viewer under **Session › Data mapping**
- the keyword is accepted by `show_anndata(...)`, `serve_anndata(...)`,
  `AnnDataViewer(...)`, and `AnnDataServer(...)`
- naming a graph the store does not carry stops the server before it binds:
  `knn_key names obsp['connectivities'], which this store does not have.`
  Compute one with `sc.pp.neighbors(adata)` first

The graph is read only when a view draws it, so it costs nothing at startup.

---

## Symptom: Performance is extremely slow / browser crashes

### Likely causes (ordered)
- the first gene of a large row-major (CSR) matrix: the browser builds a gene
  index once, and later genes are fast
- a large in-memory AnnData whose cell-major `X` is still being laid out
  gene-major (in the background from the start; the first gene waits for it)
- huge categorical fields or heavy overlays

### Fix
- serve a path rather than an in-memory object:
  `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`
- write a Cellucid-ready AnnData once (`prepare`) and open that: its matrix is
  gene-major, so no gene index is built
- reduce or remove pathological fields (very high cardinality)

Performance guide: {doc}`14_performance_scaling_and_lazy_loading`.

---

## Reporting a bug (what to include)

Please capture:

- the exact command you ran (`cellucid serve ...`), including flags
- whether the input was `.h5ad` vs Zarr store vs in-memory, and whether
  `prepare()` wrote it
- the full startup transcript: the numbered steps and the banner
- the last `Error:` line verbatim — it is written after everything else, and it
  names the one condition to correct
- the output of:
  - `/_cellucid/health`
  - `/_cellucid/info`
  - `/_cellucid/datasets`
- if one gene or one continuous field failed: the full text of its
  `Gene cannot be shown` notification — it already names the field, counts the
  offending values, and lists the first affected cell indices
- in notebooks: `viewer.debug_connection()` output
- if safe: `adata.uns["cellucid"]` of the file

This will save you (and maintainers) a lot of back-and-forth.

## Server-side symptoms

A banner URL that will not open, a startup that prints nothing for minutes, and
an import Cellucid says is missing are covered in
{doc}`../i_troubleshooting_index/04_server_mode_issues`. Cluster-specific
failures are in {doc}`../i_troubleshooting_index/07_hpc_and_remote_access_issues`.
