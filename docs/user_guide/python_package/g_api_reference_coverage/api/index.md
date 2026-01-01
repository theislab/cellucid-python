# API Reference (cellucid-python)

This section is the **canonical reference** for the public `cellucid` Python API (and the `cellucid` CLI).

**Choose your path**
- **Wet lab / beginner:** start with {doc}`jupyter` (notebooks) or {doc}`cli` (terminal).
- **Computational user:** start with {doc}`export` (reproducible exports) + {doc}`server` (sharing + lazy loading).
- **Power user / developer:** start with {doc}`sessions` (session bundles) + {doc}`vector_fields` (velocity/drift overlays).

---

## Which function should I call?

| You have… | You want… | Use… |
|---|---|---|
| `AnnData` (in memory / `.h5ad` / `.zarr`) | The fastest “just show me” notebook view | {func}`~cellucid.show_anndata` |
| Exported dataset directory | Notebook embedding | {func}`~cellucid.show` |
| Exported dataset directory | Browser tab (local HTTP server) | {func}`~cellucid.serve` or `cellucid serve …` |
| `.h5ad` / `.zarr` | Browser tab without exporting first | {func}`~cellucid.serve_anndata` or `cellucid serve …` |
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

Run a local server for exported datasets or AnnData (including SSH-tunnel workflows).
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Export / Data Preparation
:link: export
:link-type: doc

`prepare(...)` and the exported on-disk format used by the web viewer.
:::

:::{grid-item-card} {octicon}`device-desktop;1.5em;sd-mr-1` Viewers
:link: viewers
:link-type: doc

`CellucidViewer` and `AnnDataViewer` classes (Jupyter widgets + hooks).
:::

:::{grid-item-card} {octicon}`database;1.5em;sd-mr-1` Adapters
:link: adapters
:link-type: doc

`AnnDataAdapter` (server-side adapter for `.h5ad`/`.zarr`/in-memory AnnData).
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

`cellucid serve` auto-detects `.h5ad`, `.zarr`, or exported directories.
:::

::::

---

## Reference conventions

- “**Exported dataset directory**” means a folder created by {func}`~cellucid.prepare` (contains `dataset_identity.json`, `obs_manifest.json`, and `points_*d.bin(.gz)`).
- “**AnnData mode**” means serving data directly from an `AnnData` object or `.h5ad`/`.zarr` without running `prepare` first (convenient, but typically slower).
- Many APIs are **lazy-imported** from `cellucid.__init__` to keep CLI startup fast; importing `cellucid` is intentionally lightweight.

```{toctree}
:maxdepth: 2
:hidden:

jupyter
server
export
viewers
adapters
sessions
vector_fields
cli
```
