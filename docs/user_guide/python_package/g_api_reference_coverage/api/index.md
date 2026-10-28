# API Reference (cellucid-python)

This section is the **canonical reference** for the public `cellucid` Python API (and the `cellucid` CLI).

**Choose your path**
- **Wet lab / beginner:** start with {doc}`jupyter` (notebooks) or {doc}`cli` (terminal).
- **Computational user:** start with {doc}`prepare` (Cellucid-ready AnnData) + {doc}`server` (serving AnnData).
- **Power user / developer:** start with {doc}`sessions` (session bundles) + {doc}`vector_fields` (velocity/drift overlays).

---

## Which function should I call?

| You have… | You want… | Use… |
|---|---|---|
| `AnnData` or `MuData` (in memory), or any file or folder Cellucid reads (`.h5ad`, AnnData or MuData Zarr, `.h5mu`, `.loom`, a 10x `.h5` or `outs/` folder, a MatrixMarket folder, `.rds`/`.rda`/`.RData`) | The fastest “just show me” notebook view | {func}`~cellucid.show_anndata` |
| Arrays, an `AnnData` or a `MuData` to share or host | A Cellucid-ready AnnData (or MuData) file | {func}`~cellucid.prepare` |
| Any file or folder Cellucid reads (`.h5ad`, AnnData or MuData Zarr, `.h5mu`, `.loom`, a 10x `.h5` or `outs/` folder, a MatrixMarket folder, `.rds`/`.rda`/`.RData`) | Browser tab (local HTTP server) | {func}`~cellucid.serve_anndata` or `cellucid serve …` (with `--dataset-name` and `--dataset-id` unless the AnnData or MuData declares `uns['cellucid']['name']` and `['id']`, as every file `prepare` writes does) |
| A viewer, or the web app's exports (CSV, GMT, session) | Everything made in Cellucid (highlights, labels, scores, gene sets) as `obs` columns, in one call | {meth}`viewer.pull(adata) <cellucid.AnnDataViewer.pull>` or {func}`~cellucid.pull` |
| An `obs` column or a gene set | Shown in the viewer | {meth}`viewer.push(...) <cellucid.AnnDataViewer.push>`, {meth}`viewer.push_gene_set(...) <cellucid.AnnDataViewer.push_gene_set>` |
| A session captured in the web app | Apply highlights/fields back onto `AnnData` | {class}`~cellucid.CellucidSessionBundle` + {func}`~cellucid.apply_cellucid_session_to_anndata` |
| A transition matrix + embedding | Vector field overlay in Cellucid | {func}`~cellucid.compute_transition_drift` / {func}`~cellucid.add_transition_drift_to_obsm` |

---

## API map

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {octicon}`play;1.5em;sd-mr-1` Jupyter
:link: jupyter
:link-type: doc

Notebook embedding, hooks/events, and controlling the viewer from Python.
:::

:::{grid-item-card} {octicon}`server;1.5em;sd-mr-1` Server
:link: server
:link-type: doc

Run a local server for any file or folder Cellucid reads, or an AnnData or MuData object (including SSH-tunnel workflows).
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Data Preparation
:link: prepare
:link-type: doc

`prepare(...)`: write a Cellucid-ready AnnData (`.h5ad`, `.zarr`, `.zarr.zip`) or MuData (`.h5mu`, MuData Zarr).
:::

:::{grid-item-card} {octicon}`device-desktop;1.5em;sd-mr-1` Viewers
:link: viewers
:link-type: doc

The `AnnDataViewer` class (Jupyter widget + hooks).
:::

:::{grid-item-card} {octicon}`sync;1.5em;sd-mr-1` Round trip
:link: round_trip
:link-type: doc

Pull highlights, labels, scores and gene sets into `obs` in one call; push columns and gene sets; GMT.
:::

:::{grid-item-card} {octicon}`archive;1.5em;sd-mr-1` Sessions
:link: sessions
:link-type: doc

`.cellucid-session` bundles: read, inspect, and apply back to AnnData.
:::

:::{grid-item-card} {octicon}`graph;1.5em;sd-mr-1` Vector Fields
:link: vector_fields
:link-type: doc

Compute/store vector overlays (e.g., CellRank drift) for visualization.
:::

:::{grid-item-card} {octicon}`terminal;1.5em;sd-mr-1` CLI
:link: cli
:link-type: doc

`cellucid serve`: detects and serves any file or folder Cellucid reads.
:::

::::

---

## Reference conventions

- “**Cellucid-ready AnnData**” means a standard AnnData written by {func}`~cellucid.prepare` (or R's `cellucid_prepare()`): gene-major, chunked for range reads, naming itself in `uns['cellucid']`.
- Every viewing API hands a file or folder, or an in-memory `AnnData` or `MuData`, to the viewer, which reads it itself: in place with HTTP range requests where the format allows it (AnnData, MuData, Loom, 10x), once from the start for an R object or a MatrixMarket folder. An in-memory AnnData, and an AnnData, MuData, Loom or 10x file (or `outs/` folder) not in priority order, are served as a store with their cells in priority order, made as the viewer reads it and holding what the viewer presents from the file itself; an in-memory MuData is first written once to a temporary `.h5mu`.
- Many APIs are **lazy-imported** from `cellucid.__init__` to keep CLI startup fast; importing `cellucid` is intentionally lightweight.

```{toctree}
:maxdepth: 2
:hidden:

jupyter
server
prepare
viewers
round_trip
sessions
subsample
vector_fields
cli
```
