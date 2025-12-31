# “Which workflow is for me?” decision tree

**Audience:** everyone (especially mixed teams)  
**Time:** 5–10 minutes  
**Goal:** pick the simplest workflow that matches your data + sharing needs

---

## The short answer (most common choices)

- If you are a **wet-lab scientist** and someone already exported data for you → load the **Prepared** folder in the web app.
- If you are a **computational user with AnnData** → use `cellucid-python` (`show_anndata(...)` or `prepare(...); show(...)`).
- If you need **multi-user labeling/voting** → use **Community Annotation** (GitHub-backed).

<!-- SCREENSHOT PLACEHOLDER
ID: workflow-decision-tree-diagram
Where it appears: Which workflow is for me? → Decision tree
Capture:
  - Asset type: diagram (SVG preferred)
  - Content: a minimal decision tree that matches this page (Prepared vs AnnData vs Hosting vs Community annotation)
  - Goal: let a non-technical reader choose a path in <30 seconds
Crop:
  - Include only the diagram
Redact:
  - Use generic dataset labels (e.g., `pbmc_demo`)
Alt text:
  - Decision tree for choosing a Cellucid workflow.
Caption:
  - Decision tree for selecting the simplest Cellucid workflow based on what data you have and how you need to share it.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder diagram for the workflow decision tree.
:width: 100%

Decision tree for selecting the simplest Cellucid workflow based on what data you have and how you need to share it.
```

---

## Decision tree (plain language first)

::::{tab-set}

:::{tab-item} Wet‑Lab / New to Coding

1) Do you already have a folder called something like `export/` from a colleague?
   - Yes → use **Local data → Prepared** in the web app.
   - No → ask a computational colleague to export one with `cellucid-python` (this is the normal path).

2) Do you need to *share* what you did?
   - Yes → use **Save State** to download a `.cellucid-session` file, then send it.
   - If you need many people to vote on labels → use **Community Annotation** instead (it’s designed for that).

3) Does the dataset feel too big / slow?
   - Ask for a “server mode” workflow (your colleague can run a small local server so the browser doesn’t have to load everything eagerly).

:::

:::{tab-item} Computational / Power User

Start here:

```
Do you already have a prepared export folder?
  ├─ Yes → open in the web app (local Prepared or hosted HTTP server)
  └─ No →
       Do you have the dataset in Python (AnnData / arrays)?
         ├─ Yes → use cellucid-python:
         │        - fastest: show_anndata(adata)
         │        - reproducible/sharable: prepare(..., out_dir=...); show(out_dir)
         └─ No → you need to obtain a usable export or a server endpoint first
```

Then decide on sharing:

- One-to-one sharing / reproducibility → `.cellucid-session` bundles
- Many-to-many annotation / voting → Community Annotation repo

:::

::::

---

## Workflow table (when to use what)

| You want to… | Recommended workflow | Why | Common gotcha |
|---|---|---|---|
| Try Cellucid quickly | Demo dataset | Zero setup | Not your data |
| Explore your own data locally | Web app: Local data → H5AD/Zarr/Prepared | No server needed | Browser memory limits on huge datasets |
| Share a dataset with collaborators | Host the prepared export folder (HTTP) | Everyone opens the same dataset URL | CORS / hosting config |
| Work in notebooks | `cellucid-python` (embedded viewer + hooks) | Tight Python ↔ UI loop | Notebook iframe restrictions (fullscreen/pointer lock) |
| Compare multiple hypotheses | Multiview snapshots (“Keep view”) | Side-by-side views | Smoke mode is disabled when snapshots exist |
| Run a multi-user labeling round | Community Annotation | Conflict-free collaboration | Requires GitHub repo + app installation |

---

## Edge cases and “which path wins?”

### “I can load H5AD directly, why export?”

Direct H5AD/Zarr loading is convenient, but exporting often wins when you need:

- fast repeated loads,
- deterministic sharing (“this exact export corresponds to this paper figure”),
- server mode / lazy loading patterns.

### “My dataset is huge (millions of cells). What should I do?”

Prefer:

- a prepared export (optimized layout), and/or
- server mode (so the browser doesn’t need to ingest everything eagerly).

Also plan for performance from day one:

- keep category counts reasonable,
- avoid exploding the number of views,
- be cautious with smoke mode.

### “I have sensitive data”

Choose workflows based on your constraints:

- local prepared folder (no upload) is simplest,
- remote hosting requires careful access control,
- community annotation uses GitHub (be sure your repo is private if needed).

---

## Troubleshooting (choosing a workflow)

### Symptom: “I’m not sure what I have”

**How to confirm**

- If you have a folder with files like `dataset_identity.json` and embedding/field binaries → you likely have a prepared export.
- If you have a `.h5ad` file → use H5AD loading or `cellucid-python`.
- If you have a `.zarr` directory → use Zarr loading or `cellucid-python`.

### Symptom: “My colleague sent me a folder but the app can’t load it”

Likely you need the data-loading troubleshooting section (`b_data_loading/index`) and/or the system requirements page (`a_orientation/02_system_requirements`).

---

## Next steps

- If you picked a web workflow: continue to `b_data_loading/index`
- If you picked the notebook/Python workflow: go to `python_package/index`
