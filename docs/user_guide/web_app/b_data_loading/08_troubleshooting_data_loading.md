# Troubleshooting (data loading)

**Audience:** everyone (quick checks for beginners; deeper diagnostics for computational users)  
**Time:** 5–20 minutes (depending on the symptom)  
**What you’ll learn:** how to quickly localize whether a failure is (a) browser/UI, (b) network/CORS, or (c) your data format  
**Prerequisites:** none

---

## First: identify your loading path (so you debug the right thing)

Before deep-diving, answer:

1) Are you loading:
   - an **export folder** (recommended),
   - a **`.h5ad` / `.zarr`** file,
   - a **remote server** (`?remote=...`),
   - or a **GitHub exports root** (`?github=...`)?
2) Is it **browser-only** (no Python running) or **server-backed** (Python process running)?

If you’re unsure, start with {doc}`01_loading_options_overview`.

---

## 2-minute triage checklist (do this before anything else)

### Step 1 — Confirm the web app itself works

- Load a built-in demo dataset (if available).
- If demos fail too, this is likely an environment/browser issue:
  - see system requirements: {doc}`../a_orientation/02_system_requirements`

### Step 2 — Confirm you’re not fighting the wrong workflow

Common mismatch:
- Trying to load a very large `.h5ad` directly in the browser (often fails due to memory).

If your `.h5ad` is large:
- switch to server mode: {doc}`04_server_tutorial`

### Step 3 — Capture the “shape” of the failure

- “Nothing happens when I click” → file picker permissions / browser issues
- “It connects but then fails to fetch something” → network/CORS/URL issue
- “It loads but is blank / missing fields” → data missing embeddings/fields

---

## Symptom → diagnosis → fix (common issues)

### Symptom: “I clicked Folder / .h5ad / .zarr, but nothing happens”

**Likely causes (ordered)**
1) Browser does not support the required file picker capability (common on Safari / embedded browsers).
2) You denied the permission prompt (folder access).
3) The app is embedded in a context that blocks file access (some notebook/iframe setups).

**How to confirm**
- Try the same action in Chrome or Edge.
- Try a very small folder/file to rule out “large file = slow UI”.

**Fix**
- Use Chrome/Edge.
- If folder picking is blocked, use server mode instead: {doc}`04_server_tutorial`
- If you must stay browser-only, prefer pre-exported folders over `.h5ad`.

**Prevention**
- For workshops/collaborators, standardize on a supported browser + desktop environment.

---

### Symptom: “It loads forever / spinner never ends”

**Likely causes (ordered)**
1) Very large `.h5ad` loaded directly in the browser (browser memory pressure).
2) Remote server URL is wrong or unreachable.
3) GitHub raw fetch blocked (corporate firewall/content blocker).

**How to confirm**
- Try loading a demo dataset.
- If using server mode, open the health endpoint:

  ```text
  http://127.0.0.1:8765/health
  ```

- If using GitHub, open the raw manifest URL in a new tab (replace placeholders):

  ```text
  https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/datasets.json
  ```

**Fix**
- Large `.h5ad` → use server mode or export first:
  - server: {doc}`04_server_tutorial`
  - export + folder picker: {doc}`03_browser_file_picker_tutorial`
- Server mode → ensure URL/port matches the terminal output.
- GitHub → try a different network, disable blockers, or switch to server mode.

---

### Symptom: “No embedding / it says no UMAP”

**Likely causes**
- AnnData is missing required `obsm` keys.
- Export folder is missing `points_2d.bin(.gz)` / `points_3d.bin(.gz)`.

**How to confirm**
- AnnData (Python):

  ```python
  print(adata.obsm.keys())
  ```

- Exports: check that your export folder has at least one `points_*d.bin` file and that `dataset_identity.json["embeddings"]` lists it.

**Fix**
- Compute UMAP and store it in `obsm` under a supported key, or re-export with embeddings.
- Prefer exports for the smoothest experience.

---

### Symptom: “Dataset loads, but I see no points (blank canvas)”

**Likely causes (ordered)**
1) WebGL/GPU issue (context lost, unsupported browser/GPU).
2) Embedding contains NaN/Inf or extreme values.
3) All points are filtered/hidden (rare on first load, but possible with a saved session).

**How to confirm**
- Try a demo dataset.
- Open browser console and look for WebGL/context errors.
- For data sanity, check embeddings for NaN/Inf (Python):

  ```python
  import numpy as np
  X = adata.obsm.get("X_umap")
  print(np.isfinite(X).all(), X.min(), X.max())
  ```

**Fix**
- WebGL issues → see {doc}`../a_orientation/02_system_requirements`
- NaN/Inf → clean/recompute embeddings and re-export.

---

### Symptom: “Fields list is empty / missing expected metadata”

**Likely causes**
- `obs` is empty or was not exported.
- Field names changed (session/expectations mismatch).

**How to confirm**
- AnnData: `print(adata.obs.head())`
- Exports: verify `obs_manifest.json` exists and `dataset_identity.json["obs_fields"]` is non-empty.

**Fix**
- Ensure `obs` is populated before exporting.
- Re-export with `obs=adata.obs` (and appropriate categorical/continuous handling).

---

### Symptom: “Gene search returns nothing / genes are missing”

**Likely causes (ordered)**
1) No gene expression provided (`adata.X` missing/empty), so there is nothing to query.
2) Gene identifiers are not where Cellucid expects (e.g., you want symbols but stored elsewhere).
3) Data is extremely large and gene fetch is slow (looks like “nothing”).

**How to confirm**
- AnnData: `print(adata.X is None)`, `print(adata.n_vars)`, `print(adata.var.head())`
- Exports: ensure `var_manifest.json` and `var/` exist (gene expression exported).

**Fix**
- Provide a valid gene expression matrix and gene ids; if needed, pick the right gene id column in server/Jupyter mode (`gene_id_column=...`).
- For very large datasets, use server mode (best reliability): {doc}`04_server_tutorial`

---

### Symptom: “Vector field overlay toggle is disabled / no fields appear”

**Likely causes (ordered)**
1) The dataset truly has no vector fields (common).
2) Vector fields exist, but they weren’t exported / served (missing `vectors/` or wrong server path).
3) Naming mismatch: vectors are present but not discoverable (e.g., not `*_umap_2d` / `*_umap_3d`).
4) Dimension mismatch: you only have 2D vectors but you’re looking at 3D (or vice versa).

**How to confirm**
- **Exports**:
  - Open `<export_dir>/dataset_identity.json` and search for `vector_fields`.
  - Confirm `<export_dir>/vectors/` exists and contains `*_2d.bin(.gz)` / `*_3d.bin(.gz)` files.
- **AnnData**:
  - Print `adata.obsm.keys()` and look for `velocity_umap_2d`-style entries.

**Fix**
- Provide/export vector fields using the supported conventions and re-load the dataset.
- Switch to the dimension that actually has vectors (2D vs 3D).

For overlay behavior and deeper diagnostics, see:
- {doc}`../i_vector_field_velocity/index`
- {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

---

### Symptom: “GitHub connect says `datasets.json not found`”

**Likely causes (ordered)**
1) Wrong repo path (pointing at the dataset folder instead of the exports root).
2) Missing `datasets.json` (you didn’t generate/commit it).
3) Corporate network blocks `raw.githubusercontent.com`.

**How to confirm**
- Open the raw URL in a browser (replace placeholders):

  ```text
  https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/datasets.json
  ```

**Fix**
- Generate and commit `datasets.json` at the exports root.
- Confirm you’re pointing Cellucid at the exports root, not at `.../pbmc_demo/`.
- If blocked, use server mode or a different network.

See {doc}`02_local_demo_tutorial`.

---

### Symptom: “Remote server connect fails / CORS blocked / mixed content”

**Likely causes**
- Wrong URL (host/port mismatch).
- Server is bound to `127.0.0.1` on a remote machine, but you’re trying to access it directly.
- Browser blocks non-localhost HTTP in some cases (mixed content).

**How to confirm**
- Check server banner output and copy the exact viewer URL printed by the server.
- Visit the health endpoint.

**Fix**
- For remote machines: use an SSH tunnel (recommended): {doc}`04_server_tutorial`
- Keep the server bound to `127.0.0.1` and access via `http://localhost:<port>` through the tunnel.

---

## When to stop debugging and switch workflows (high leverage)

If you hit one of these, switching workflows saves time:

- **Large `.h5ad` in browser** → switch to server mode or exports.
- **GitHub too large / blocked** → switch to server mode or host exports on an accessible HTTP server.
- **File picker blocked** → switch to server mode.

---

## What to include when asking for help (so others can debug quickly)

Copy/paste:
- your loading path (export folder / `.h5ad` / `.zarr` / server / GitHub)
- approximate dataset size: `n_cells`, `n_genes`
- browser + OS
- any console error text
- for server mode: the exact `cellucid serve ...` command and the printed URL

---

## Related pages

- Loading options map: {doc}`01_loading_options_overview`
- Local file picker: {doc}`03_browser_file_picker_tutorial`
- Server mode: {doc}`04_server_tutorial`
- GitHub exports workflow: {doc}`02_local_demo_tutorial`
- Dataset identity: {doc}`06_dataset_identity_why_it_matters`
- Vector fields overlay: {doc}`../i_vector_field_velocity/index`
