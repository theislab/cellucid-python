# What is cellucid-python?

**Audience:** wet lab scientists, computational users, developers  
**Time:** 5–15 minutes  
**Goal:** understand what `cellucid` (cellucid-python) does and which workflow to use

`cellucid-python` is the repository folder name in this monorepo. The **Python package you install** is named **`cellucid`**:

```bash
pip install cellucid
```

If you remember one sentence:

> `cellucid` opens your AnnData in the **Cellucid web app** — served from your machine, or embedded in a notebook — and writes **Cellucid-ready AnnData** files that load fastest and can be shared or hosted anywhere.

```{note}
Cellucid itself is a **web app** (the viewer UI). `cellucid-python`,
`cellucid-r`, and `cellucid-annotation` provide language and workflow
integrations around that viewer.
```

## What `cellucid` (Python) does (today)

At a high level, the Python package supports three related workflows:

1) **View immediately: AnnData → viewer**  
   - Notebook: `show_anndata(adata, dataset_name="My dataset", dataset_id="my-dataset")`
     (or the same identity arguments with `"data.h5ad"`; both may be omitted
     when the data declares `uns['cellucid']['name']` and `['id']`)
   - Browser: `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`

2) **Prepare for speed + sharing: AnnData, MuData or arrays → Cellucid-ready AnnData (or MuData) → viewer**  
   - Python: `prepare("my_dataset.h5ad", format="h5ad", data=adata, dataset_id=..., dataset_name=...)`
     (or the components one by one: `embeddings=...`, `obs=...`, ...; a
     MuData is written as `.h5mu` or MuData Zarr)
   - Then view in a browser: `cellucid serve my_dataset.h5ad`  
     (or embed in a notebook: `show_anndata("my_dataset.h5ad")`), or host the
     file and open it from a link

3) **Notebook integration (bidirectional): viewer ↔ Python**  
   - Send commands to the UI (highlight, color-by, visibility, reset camera,
     data mapping)
   - Receive events from the UI (selection/hover/click/ready/data mapping) as
     Python callbacks

### The other Cellucid repos (how they relate)

- **Cellucid web app**: the browser UI you click around in (rendering, filters, highlights, analysis, figure export, sessions).
- **`cellucid` (cellucid-python)**: prepare/serve/embed AnnData + hooks/events for notebooks.
- **`cellucid-annotation`**: community annotation workflows (multi-user, GitHub-backed collaboration).
- **`cellucid-r`**: writes the same Cellucid-ready AnnData from Seurat,
  SingleCellExperiment, matrices, and data frames, and shows R objects in the
  viewer from an R session (`cellucid_view()`, `cellucid_serve()`); see
  {doc}`../../r_package/index`.

## Fast path (wet lab / beginner / non-technical)

If you’re a wet lab scientist or a non-technical collaborator, it helps to think in terms of “what do I open?”:

- You usually open **a dataset in the browser** (Cellucid web app).
- That dataset is either:
  - an **AnnData file** (`.h5ad`, `.zarr`, `.zarr.zip`) someone gave you or published, or
  - a dataset served by a collaborator’s machine/server.

**What you typically do next (in the web app)**
1) Load the dataset
2) Color by a field (clusters, condition, QC metric, or a gene)
3) Select a population
4) Compare groups / export a figure / save a session

```{figure} ../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: The whole Cellucid window with the Pancreas sample coloured by cell_type: the sidebar shows Coloring and filtering with the categorical field selected and its eight-row legend, and the 3D embedding fills the canvas in the matching colours.
:width: 1440px

A loaded dataset in Cellucid: the sidebar controls the active view while the categorical legend maps directly to the colored points.
```

## Practical path (computational users)

For computational users, the main question is: **Do I open my AnnData as it is, or prepare it first?**

### Option A — open it as it is (fast iteration, nothing converted)

Use this when:
- you want to **inspect a dataset quickly**,
- you’re working interactively and don’t care about a shareable on-disk artifact *yet*.

Typical entry points:
- Notebook: direct AnnData viewing with `dataset_name` and `dataset_id`
  (taken from `uns['cellucid']` when the data declares them)
- Browser: `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`

Tradeoffs:
- Very little setup
- The viewer reads the file or store itself, in place, with HTTP range
  requests where the format allows it: only the bytes a view needs are
  transferred, and nothing is converted first. An in-memory AnnData is served
  with its cells in priority order, each chunk made from the object when the
  viewer asks for it; a MuData object is written once to a temporary `.h5mu`.
- Not a deterministic “artifact” you can hand to collaborators

### Option B — prepare first (fastest viewing; one more step)

Use this when:
- you want the **fastest** experience in the web app (a gene-major matrix makes
  every gene one small read),
- you want a **shareable file** for collaborators or papers that opens with a
  name, a citation, and chosen defaults,
- you will host the dataset on the web.

Typical entry points:
- Python: `prepare("my_dataset.h5ad", format="h5ad", embeddings=..., dataset_id=..., dataset_name=..., ...)`
- Then view: `cellucid serve my_dataset.h5ad` (or open the file in the web app)

Tradeoffs:
- Fast loading and consistent performance
- Still a standard AnnData: every AnnData tool reads it
- Can be hosted (static) or served locally
- Requires you to choose inputs (embeddings, latent space, obs, etc.)

### Supported “starting points”

Cellucid-python supports viewing (with `cellucid serve`, `serve_anndata`
and `show_anndata`) from:
- **In-memory AnnData or MuData** objects (an AnnData served in priority
  order as the viewer reads it, a MuData written once to a temporary `.h5mu`)
- **`.h5ad`** files and **AnnData Zarr** stores, a folder or zipped (`.zip`)
- **`.h5mu`** files and **MuData Zarr** stores
- **`.loom`** files, **10x Cell Ranger `.h5`** files and **`outs/`** folders
- **MatrixMarket folders** (`matrix.mtx.gz` with barcodes and features, also
  GEO's `GSM..._matrix.mtx.gz` naming) and **R objects** (`.rds`, `.rda`,
  `.RData`, gzipped or not), which the viewer reads once from the start
- **Cellucid-ready AnnData or MuData** written by `prepare(...)`

A file `prepare()` wrote, an R object and a MatrixMarket folder are served as
they are on disk, and the browser reads them. Any other AnnData, MuData, Loom
or 10x file is served as a store with its cells in priority order, made from
the file as the browser reads it, so a subsample the viewer opens transfers
only its cells.

## Deep path (developer / maintainer)

The architecture is intentionally simple:

- Python answers the viewer's reads of one dataset's bytes (an `.h5ad`, a
  Zarr store, an `.h5mu`, a Loom or 10x file, a MatrixMarket folder, an R
  object, ...); the viewer reads it itself, with range requests where the
  format allows it.
- The **viewer UI** is the Cellucid web app.
  - In a notebook (`show_anndata`) it is a widget: the browser loads the app
    from `https://www.cellucid.com` into the cell's output, and the app's
    reads and messages travel to the kernel through the notebook's own
    connection (a Jupyter comm). No port, no server.
  - In a browser tab (`cellucid serve`, `serve_anndata`, `AnnDataViewer`) a
    local HTTP server establishes the exact verified app build and serves it
    from the **same origin** as the dataset.
- Bidirectional communication: one channel, the viewer channel (the widget's
  comm in a notebook, a WebSocket at `/_cellucid/ws` for a browser tab),
  authenticated with the viewer's token:
  - **Python → frontend**: commands on it (`viewer.send_message(...)` + convenience methods)
  - **Frontend → Python**: events on it (hooks), selections as binary cell indices

If you maintain or extend `cellucid`, the “mental model” to keep in mind is:

> A viewer is “just” a dataset server + a browser UI + a small message protocol to coordinate state.

## Key terms (used throughout this guide)

- **Cellucid-ready AnnData**: a standard AnnData written by `prepare()` for fast browser reads.
- **Dataset identity**: `uns['cellucid']['id']` and `['name']` (or `dataset_id`/`dataset_name` passed when serving); sessions are keyed by the id.
- **Embedding**: an `obsm` array the viewer plots, any number per dimension.
- **`obs` fields**: per-cell metadata (categorical or continuous).
- **`var` / gene expression**: per-gene metadata + expression values (optional).
- **Viewer**: the Python object that shows the UI and exposes hooks/commands
  (`cellucid.Widget` in a notebook, `AnnDataViewer` for a browser tab).
- **Hooks / events**: selection/hover/click/ready callbacks from the UI back into Python.

## Edge cases and limitations (high-level)

- **Web UI availability**: a notebook widget needs the *browser* to reach
  `https://www.cellucid.com` (or the copy named by `web_source_url=`); a
  browser-tab server establishes the complete exact build published by the
  configured web source at startup, and a source failure is raised before the
  server binds.
- **Notebooks without widget support** cannot show the viewer; see
  {doc}`03_compatibility_matrix_must_be_explicit`.

## Troubleshooting (common misconceptions)

### “Do I have to upload my data to a server?”

No by default:
- When you run `cellucid serve ...`, your data is served from your machine (default: `127.0.0.1` only).
- When you use `show_anndata(...)` in a notebook, the viewer reads the dataset
  from your kernel through the notebook's own connection; nothing is uploaded.

The viewer UI itself (HTML/JavaScript/CSS) comes from `https://www.cellucid.com`:
loaded by the browser for a notebook widget, or downloaded, verified and
served locally by a browser-tab server.

### “I installed `cellucid` but I don’t see anything”

Most common causes:
- you ran `show_anndata(...)` **outside** a notebook: it raises and names the
  browser-tab alternatives (`AnnDataViewer`, `serve_anndata`, `cellucid serve`),
- the notebook front end does not show widgets (`anywidget` missing where the
  notebook server runs),
- the browser cannot load the app from `https://www.cellucid.com`.

Start with:
- {doc}`../installation`
- {doc}`03_compatibility_matrix_must_be_explicit`
- {doc}`04_quick_start_3_levels`

## Next steps

- Install and verify the CLI/notebook prerequisites: {doc}`../installation`
- Check environment constraints early (VS Code, Colab, remote notebooks): {doc}`03_compatibility_matrix_must_be_explicit`
- Get your first view running: {doc}`04_quick_start_3_levels`
