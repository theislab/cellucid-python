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
| `AnnData` (in memory / `.h5ad` / Zarr store) | The fastest “just show me” notebook view | {func}`~cellucid.show_anndata` |
| Arrays / an `AnnData` to share or host | A Cellucid-ready AnnData file | {func}`~cellucid.prepare` |
| `.h5ad` / Zarr store | Browser tab (local HTTP server) | {func}`~cellucid.serve_anndata` or `cellucid serve …` (with `--dataset-name` and `--dataset-id` unless the store declares `uns['cellucid']['name']` and `['id']`, as every file `prepare` writes does) |
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

Run a local server for an AnnData file, store or object (including SSH-tunnel workflows).
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Data Preparation
:link: prepare
:link-type: doc

`prepare(...)`: write a Cellucid-ready AnnData (`.h5ad`, `.zarr`, `.zarr.zip`).
:::

:::{grid-item-card} {octicon}`device-desktop;1.5em;sd-mr-1` Viewers
:link: viewers
:link-type: doc

The `AnnDataViewer` class (Jupyter widget + hooks).
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

The CLI serve subcommand detects an `.h5ad`, a Zarr store folder, or a zipped Zarr store.
:::

::::

---

## Reference conventions

- “**Cellucid-ready AnnData**” means a standard AnnData written by {func}`~cellucid.prepare` (or R's `cellucid_prepare()`): gene-major, chunked for range reads, naming itself in `uns['cellucid']`.
- Every viewing API hands an `.h5ad` file, a Zarr store, or an in-memory `AnnData` to the viewer, which reads the store itself, in place, with HTTP range requests (an in-memory object is first written once to a temporary `.h5ad`).
- Many APIs are **lazy-imported** from `cellucid.__init__` to keep CLI startup fast; importing `cellucid` is intentionally lightweight.

```{toctree}
:maxdepth: 2
:hidden:

jupyter
server
prepare
viewers
sessions
vector_fields
cli
```
