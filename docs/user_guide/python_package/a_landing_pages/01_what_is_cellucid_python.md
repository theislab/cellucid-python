# What is cellucid-python?

**Audience:** wet lab scientists, computational users, developers  
**Time:** 5–15 minutes  
**Goal:** understand what `cellucid` (cellucid-python) does and which workflow to use

`cellucid-python` is the repository folder name in this monorepo. The **Python package you install** is named **`cellucid`**:

```bash
pip install cellucid
```

If you remember one sentence:

> `cellucid` turns your data (often AnnData) into something the **Cellucid web app** can load, and it can also **serve/embed** the web app so you can use it from the CLI and notebooks.

```{note}
Cellucid itself is a **web app** (the viewer UI). `cellucid-python` and `cellucid-annotation` are helper repos. `cellucid-r` is planned but not ready yet.
```

## What `cellucid` (Python) does (today)

At a high level, the Python package supports three related workflows:

1) **View immediately (no export): AnnData → viewer**  
   - Notebook: `show_anndata(adata)` (or `show_anndata("data.h5ad")`)
   - Browser: `cellucid serve data.h5ad` (auto-detects `.h5ad` / `.zarr`)

2) **Export for speed + sharing: arrays → export folder → viewer**  
   - Python: `prepare(..., out_dir="exports/my_dataset")`
   - Then view in a browser: `cellucid serve exports/my_dataset`  
     (or embed in a notebook: `show("exports/my_dataset")`)

3) **Notebook integration (bidirectional): viewer ↔ Python**  
   - Send commands to the UI (highlight, color-by, visibility, reset camera)
   - Receive events from the UI (selection/hover/click/ready) as Python callbacks

### The other Cellucid repos (how they relate)

- **Cellucid web app**: the browser UI you click around in (rendering, filters, highlights, analysis, figure export, sessions).
- **`cellucid` (cellucid-python)**: export/serve/embed data + hooks/events for notebooks.
- **`cellucid-annotation`**: community annotation workflows (multi-user, GitHub-backed collaboration).
- **`cellucid-r`**: planned export helper for R; not ready yet.

## Fast path (wet lab / beginner / non-technical)

If you’re a wet lab scientist or a non-technical collaborator, it helps to think in terms of “what do I open?”:

- You usually open **a dataset in the browser** (Cellucid web app).
- That dataset is either:
  - an **export folder** someone generated for you, or
  - a dataset served by a collaborator’s machine/server.

**What you typically do next (in the web app)**
1) Load the dataset
2) Color by a field (clusters, condition, QC metric, or a gene)
3) Select a population
4) Compare groups / export a figure / save a session

<!-- SCREENSHOT PLACEHOLDER
ID: python-landing-what-is-webapp-orientation
Where it appears: What is cellucid-python? → Fast path
Capture:
  - UI location: Cellucid web app (browser)
  - State prerequisites: any small-ish dataset loaded; a categorical field selected for color-by; left sidebar visible
  - Action to reach state: load a dataset → pick a field → optionally create a small selection
Crop:
  - Include: embedding canvas + left sidebar panels + dataset name/point count
  - Exclude: browser bookmarks, personal accounts, private dataset/sample IDs
Redact:
  - Remove: any private filenames/paths, sample IDs, patient IDs, emails
Annotations:
  - Callouts: (1) dataset load location, (2) color-by control, (3) selection tool, (4) “what success looks like” area
  - Style: 2–4 callouts max; consistent color; readable at 100% width
Alt text:
  - Cellucid web app with a dataset loaded and colored by a categorical field.
Caption:
  - Orientation: where to load a dataset, choose a field for coloring, and confirm the dataset is loaded.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Cellucid web app with a dataset loaded and colored by a categorical field.
:width: 100%

Orientation: where to load a dataset, choose a field for coloring, and confirm the dataset is loaded.
```

## Practical path (computational users)

For computational users, the main question is: **Do I export first, or do I serve AnnData directly?**

### Option A — AnnData direct (fast iteration; slower viewing)

Use this when:
- you want to **inspect a dataset quickly** without deciding export options yet,
- you’re working interactively and don’t care about a shareable on-disk artifact *yet*.

Typical entry points:
- Notebook: `show_anndata(adata)` or `show_anndata("data.h5ad")`
- Browser: `cellucid serve data.h5ad`

Tradeoffs:
- Very little setup
- Can be lazy-loading for `.h5ad` / `.zarr`
- Slower than exports for repeated use and large datasets
- Not a deterministic “artifact” you can hand to collaborators

### Option B — Export-first (reproducible + fast viewing; more setup)

Use this when:
- you want the **fastest** experience in the web app,
- you want a **shareable export folder** for collaborators or papers,
- you want explicit control over **compression/quantization** and output size.

Typical entry points:
- Python: `prepare(..., out_dir="exports/my_dataset")`
- Then view: `cellucid serve exports/my_dataset` (or open the folder in the web app)

Tradeoffs:
- Fast loading and consistent performance
- Reproducible artifact you can archive/share
- Can be hosted (static) or served locally
- Requires you to choose inputs (embeddings, latent space, obs, etc.)

### Supported “starting points”

Cellucid-python supports viewing from:
- **In-memory AnnData** (notebook/server)
- **`.h5ad` file** (server and notebook; default: lazy loading via backed mode)
- **`.zarr` store** (server and notebook; array access is chunked/lazy)
- **Pre-exported directory** created by `prepare(...)`

## Deep path (developer / maintainer)

The architecture is intentionally simple:

- A **local HTTP server** serves either:
  - files in an **export folder**, or
  - “virtual files” generated from **AnnData** on demand.
- The **viewer UI** is a web app.
  - In CLI/notebook/server workflows, the Python server runs in **hosted-asset proxy** mode:
    it fetches `index.html` + `/assets/*` from `https://www.cellucid.com` and caches them locally,
    then serves them from the **same origin** as the dataset.
- Notebook embedding uses an **iframe**.  
  This matters because browser security rules differ across Jupyter/VSCode/Colab.
- Bidirectional communication:
  - **Python → frontend**: `postMessage` to the iframe (`viewer.send_message(...)` + convenience methods)
  - **Frontend → Python**: HTTP `POST` to `/_cellucid/events` on the local server (hooks)

If you maintain or extend `cellucid`, the “mental model” to keep in mind is:

> A viewer is “just” a dataset server + a browser UI + a small message protocol to coordinate state.

## Key terms (used throughout this guide)

- **Export folder**: a directory of `.json` manifests + binary files that the web app loads.
- **Dataset identity**: `dataset_identity.json` (name/id/metadata; helps reproducibility).
- **Embedding**: `points_1d.bin`, `points_2d.bin`, `points_3d.bin` (and future `4d`).
- **`obs` fields**: per-cell metadata (categorical or continuous).
- **`var` / gene expression**: per-gene metadata + expression values (optional).
- **Viewer**: a notebook object that embeds the UI and exposes hooks/commands (`CellucidViewer` / `AnnDataViewer`).
- **Hooks / events**: selection/hover/click/ready callbacks from the UI back into Python.

## Edge cases and limitations (high-level)

- **Web UI availability**: the notebook/server UI is downloaded from `https://www.cellucid.com` on first use and cached.
  If you are offline the first time, you’ll see an explicit “Viewer UI unavailable” page.
- **HTTPS notebooks**: if your notebook is served from `https://...`, your browser may block an `http://127.0.0.1:<port>` iframe (mixed content).
  The notebook embed tries to use `jupyter-server-proxy` when available; see {doc}`03_compatibility_matrix_must_be_explicit`.
- **4D is reserved**: exports can write `points_4d.bin`, but 4D visualization is not supported in the web app yet.

## Troubleshooting (common misconceptions)

### “Do I have to upload my data to a server?”

No by default:
- When you run `cellucid serve ...`, your data is served from your machine (default: `127.0.0.1` only).
- When you use `show_anndata(...)` in a notebook, the dataset is served locally by a background server.

The one thing that may use the public internet is the **viewer UI assets** (HTML/JS/CSS), which are fetched from `https://www.cellucid.com` and cached locally.

### “I installed `cellucid` but I don’t see anything”

Most common causes:
- you ran `show(...)` or `show_anndata(...)` **outside** a notebook environment (it will print a URL instead),
- the viewer UI could not be fetched (offline / firewall), so the iframe shows an error page,
- your notebook is served from HTTPS and blocks HTTP loopback iframes.

Start with:
- {doc}`02_installation`
- {doc}`03_compatibility_matrix_must_be_explicit`
- {doc}`04_quick_start_3_levels`

## Next steps

- Install and verify the CLI/notebook prerequisites: {doc}`02_installation`
- Check environment constraints early (VSCode/Colab/HTTPS notebooks): {doc}`03_compatibility_matrix_must_be_explicit`
- Get your first view running: {doc}`04_quick_start_3_levels`
