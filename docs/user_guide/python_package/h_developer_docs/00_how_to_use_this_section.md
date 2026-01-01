# How to use this section

This chapter is a **map + runbook** for people who want to *work on* `cellucid-python` (the Python package and CLI).

It is written for mixed audiences:

- **Wet lab / beginner / non-scientific**: you may still land here when something “technical” breaks (Jupyter embedding, servers, ports). Use the *Fast path* boxes and the Troubleshooting sections; skip the implementation details.
- **Computational users**: use this section when you need to understand performance, formats, or reproducibility beyond “copy/paste the quickstart”.
- **Developers / maintainers**: use this section to make changes without accidentally breaking export compatibility, security assumptions, or notebook environments.

:::{important}
This is *not* the web app developer documentation. The web app code lives in `cellucid/`.
If you’re debugging UI state/rendering/modules, start here instead: {doc}`../../web_app/p_developer_docs/index`.
:::

---

## How these pages are meant to be read

Cellucid’s docs are intentionally layered (see `cellucid/markdown/DOCUMENTATION_MASTER_GUIDE.md`):

1) **Fast path**: minimal steps to do the thing safely.
2) **Practical path**: parameters, performance, and reproducibility.
3) **Deep path**: implementation details, schemas, edge cases, and extension points.

For most work, you will bounce between:
- an architecture page (what the system *is*),
- a protocol/format page (what must not break),
- and a debugging page (how to prove what’s happening).

---

## The three “modes” you must keep straight

Almost every developer confusion in Cellucid comes from mixing these up:

### 1) **Pre-exported dataset mode** (recommended for sharing + performance)

- Python writes an export folder via `cellucid.prepare(...)`.
- The viewer loads static files like `points_3d.bin(.gz)`, `obs_manifest.json`, `var_manifest.json`.
- A lightweight HTTP server (`cellucid.server.CellucidServer`) serves those files.

Start reading: {doc}`07_prepare_export_pipeline_architecture` → {doc}`08_export_format_spec_and_invariants` → {doc}`09_server_mode_architecture_endpoints_and_security`

### 2) **AnnData server mode** (convenient, slower, great for iteration)

- You point Cellucid at an in-memory AnnData, an `.h5ad`, or a `.zarr`.
- The Python server synthesizes “virtual export files” on demand.
- The browser viewer still thinks it’s loading the export format, but there are no on-disk export files.

Start reading: {doc}`09_server_mode_architecture_endpoints_and_security` → {doc}`01_codebase_architecture`

### 3) **Notebook embedding mode** (Jupyter iframe + hooks/events)

- A server runs (exported or AnnData).
- The notebook displays an iframe pointed at that server.
- Python controls the iframe via `postMessage` (Python → frontend), and receives events via HTTP POST (frontend → Python).

Start reading: {doc}`10_jupyter_embedding_architecture` → {doc}`11_hooks_events_protocol_and_schema`

---

## When to add screenshots (and how)

Screenshots are often more helpful than logs for:
- “Where in the UI do I click?”
- “What does success look like?”
- “What error message did you see exactly?”

In this doc site (Sphinx + MyST), use the screenshot placeholder style from:
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

Rule of thumb:
- If a page contains *click-by-click steps*, add at least one screenshot placeholder (“orientation” + “success state”).
- If a page contains *a common failure mode*, add a screenshot placeholder for that failure.

You will see blocks like:

````md
<!-- SCREENSHOT PLACEHOLDER
ID: unique-id
Capture: exact steps to reach state
Alt text: one sentence
Caption: one sentence
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot.
:width: 100%

<caption text>
```
````

When you’re ready, replace the placeholder SVG path with your real PNG/SVG and keep the comment block as the “production spec”.

---

## If you are reporting a bug (what maintainers need)

To get a bug fixed quickly, include:

- **Which mode**: export folder vs AnnData server vs Jupyter embedding.
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
