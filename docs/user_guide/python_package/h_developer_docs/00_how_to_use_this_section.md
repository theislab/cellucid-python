# How to use this section

This chapter is a **map + runbook** for people who want to *work on* `cellucid-python` (the Python package and CLI).

It is written for mixed audiences:

- **Wet lab / beginner / non-scientific**: you may still land here when something “technical” breaks (Jupyter embedding, servers, ports). Use the *Fast path* boxes and the Troubleshooting sections; skip the implementation details.
- **Computational users**: use this section when you need to understand performance, formats, or reproducibility beyond “copy/paste the quickstart”.
- **Developers / maintainers**: use this section to make changes without accidentally breaking the Cellucid-ready AnnData contract, security assumptions, or notebook environments.

:::{important}
This is *not* the web app developer documentation. The web app code lives in `cellucid/`.
If you’re debugging UI state/rendering/modules, start here instead: {doc}`../../web_app/p_developer_docs/index`.
:::

---

## How these pages are meant to be read

Cellucid’s docs are intentionally layered (see {doc}`15_docs_development_and_style_guide`):

1) **Fast path**: minimal steps to do the thing safely.
2) **Practical path**: parameters, performance, and reproducibility.
3) **Deep path**: implementation details, schemas, edge cases, and extension points.

For most work, you will bounce between:
- an architecture page (what the system *is*),
- a protocol/format page (what must not break),
- and a debugging page (how to prove what’s happening).

---

## The three things you must keep straight

Almost every developer confusion in Cellucid comes from mixing these up:

### 1) **Writing a Cellucid-ready AnnData** (`prepare()`)

- Python checks an AnnData or MuData (`data=`), or arrays and data frames,
  and writes one standard AnnData (`.h5ad`, `.zarr`, `.zarr.zip`) — or, from a
  MuData, a standard MuData (`.h5mu`, MuData Zarr) — laid out for browser
  range reads, with `uns['cellucid']` naming the dataset.
- Nothing is served; the result opens anywhere an AnnData opens.

Start reading: {doc}`07_prepare_writer_architecture` →
{doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`

### 2) **Serving a dataset** (nothing converted)

- You point Cellucid at an in-memory AnnData or MuData, or at any file or
  folder of the input format registry (`src/cellucid/input_formats.json`):
  `.h5ad`, AnnData or MuData Zarr (folder or `.zip`), `.h5mu`, `.loom`, a 10x
  `.h5` file or `outs/` folder, a MatrixMarket folder, or an R
  object (`.rds`, `.rda`/`.RData`, gzipped or not).
- The Python server computes nothing: it names the dataset, serves the store's
  bytes under `/_cellucid/anndata/<dataset_id>…` with exact HTTP range
  answers, and states the launch choices in `/_cellucid/datasets`.
- The web app's own readers read the data in place (with range requests where
  the format allows it, once from the start for an R object or a MatrixMarket
  folder). An in-memory object is first written once to a temporary `.h5ad` or
  `.h5mu`.

Start reading: {doc}`09_server_mode_architecture_endpoints_and_security` → {doc}`01_codebase_architecture`

### 3) **Notebook embedding mode** (Jupyter iframe + hooks/events)

- A server runs (serving one dataset).
- The notebook displays an iframe pointed at that server.
- The iframe opens one WebSocket to the server (the viewer channel): Python's
  commands go to the viewer and the viewer's events come back on it.

Start reading: {doc}`10_jupyter_embedding_architecture` → {doc}`11_hooks_events_protocol_and_schema`

---

## When to add screenshots (and how)

Screenshots are often more helpful than logs for:
- “Where in the UI do I click?”
- “What does success look like?”
- “What error message did you see exactly?”

In this doc site (Sphinx + MyST), commit only captures from a reproduced,
visually inspected state.

Rule of thumb:

- For click-by-click steps, include an orientation or success-state capture when
  it materially clarifies the workflow.
- For a common failure mode, include the exact visible failure only when it has
  been reproduced safely with public or synthetic data.
- Write alt text for the visible controls and a caption that states what the
  reader should verify.

---

## If you are reporting a bug (what maintainers need)

To get a bug fixed quickly, include:

- **Which mode**: `prepare()` vs `cellucid serve` / `serve_anndata` vs Jupyter embedding.
- **Version info**:
  - `python -c "import cellucid; print(cellucid.__version__)"`
  - Python version + OS
- **Data scale**: `n_cells`, `n_genes`, and which embeddings (1D/2D/3D).
- **Repro steps**:
  - code (for Python/Jupyter) or click-by-click (for UI)
  - expected vs actual result
- **Logs**:
  - server logs (terminal where `cellucid serve` runs)
  - browser console logs (especially for UI/network failures)

For “I can’t embed in Jupyter / hooks don’t work / viewer blank”, include `viewer.debug_connection()` output:
{doc}`12_debugging_playbook`.
