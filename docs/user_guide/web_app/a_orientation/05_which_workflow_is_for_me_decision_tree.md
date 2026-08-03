# “Which workflow is for me?” decision tree

**Audience:** everyone (especially mixed teams)  
**Time:** 5–10 minutes  
**Goal:** pick the simplest workflow that matches your data + sharing needs

---

## The short answer (most common choices)

- If you are a **wet-lab scientist** and someone already exported data for you → load the **Prepared** folder in the web app.
- If you have your own `.h5ad` **under 512 MiB** and no Python → open it with the
  **H5AD** button in the web app, nothing installed
  ({doc}`../b_data_loading/03_browser_file_picker_tutorial`).
- If you are a **computational user with AnnData** → use `cellucid-python` (`show_anndata(...)` or `prepare(...); show(...)`).
- If you work in **R** → `cellucid_prepare()` writes the same export folder
  ({doc}`../../r_package/index`).
- If you need **multi-user labeling/voting** → use **Community Annotation** (GitHub-backed).

---

## Decision tree (plain language first)

::::{tab-set}

:::{tab-item} Wet‑Lab / New to Coding

1) Do you already have a folder called something like `export/` from a colleague?
   - Yes → use **Local data → Prepared** in the web app.
   - No, but you have a `.h5ad` under 512 MiB → use **Local data → H5AD**. You
     do not need Python for this.
   - No → ask a computational colleague to export one, with `cellucid-python`
     or with `cellucid_prepare()` in R. Both write the same folder, so it does
     not matter which they use.

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
         │        - fastest: show_anndata(adata, dataset_name="My dataset", dataset_id="my-dataset")
         │        - reproducible/sharable: prepare(..., out_dir=...); show(out_dir)
         └─ No → you need to obtain a usable export or a server endpoint first
```

If the dataset lives in R rather than Python, `cellucid_prepare()` writes the
same export folder; the viewer reads either one identically
({doc}`../../r_package/index`).

Then decide on sharing:

- One-to-one sharing / reproducibility → `.cellucid-session` bundles
- A public or intranet catalog on a static host → serve the exports root over
  HTTP and open it with `?exportsBaseUrl=https://host/exports/`
- Many-to-many annotation / voting → Community Annotation repo

:::

::::

---

## Workflow table (when to use what)

| You want to… | Recommended workflow | Why | Common gotcha |
|---|---|---|---|
| Try Cellucid quickly | Demo dataset | Zero setup | Not your data |
| Explore your own data locally | Web app: Local data → Prepared, H5AD, or Zarr ZIP | No server needed | Browser memory limits on large direct files and archives |
| Share a dataset with collaborators | Host the prepared export folder (HTTP) and open it with `?exportsBaseUrl=…` | Everyone opens the same dataset URL | CORS / hosting config; the URL must end in `/` |
| Work in notebooks | `cellucid-python` (embedded viewer + hooks) | Tight Python ↔ UI loop | Notebook iframe restrictions (fullscreen/pointer lock) |
| Prepare an export without Python | `cellucid-r` (`cellucid_prepare()`) | Same export contract, read by the same viewer | Arguments are named differently from `prepare()` |
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
- If you have a `.h5ad` file → use the browser **H5AD** control when it is under
  512 MiB, otherwise `cellucid serve` or `show_anndata()`
  ({doc}`../b_data_loading/03_browser_file_picker_tutorial`,
  {doc}`../b_data_loading/04_server_tutorial`).
- If you have a `.zarr` directory → the browser reads a *packaged* store only,
  so either package it as one `.zarr.zip` for the **Zarr ZIP** control, or point
  `cellucid serve` / `show_anndata()` at the directory itself. There is no
  browser control that accepts a `.zarr` directory
  ({doc}`../b_data_loading/03_browser_file_picker_tutorial`,
  {doc}`../b_data_loading/04_server_tutorial`).

### Symptom: “My colleague sent me a folder but the app can’t load it”

Likely you need the data-loading troubleshooting section ({doc}`../b_data_loading/index`) and/or the system requirements page ({doc}`02_system_requirements`).

---

## Next steps

- If you picked a web workflow: continue to {doc}`../b_data_loading/index`
- If you picked the notebook/Python workflow: go to {doc}`../../python_package/index`
- If your data lives in R: go to {doc}`../../r_package/index`
