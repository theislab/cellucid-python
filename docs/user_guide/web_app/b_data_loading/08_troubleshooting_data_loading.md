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
   - a browser-selected **`.h5ad`** or **Zarr ZIP** file,
   - a Python-served **`.h5ad` / `.zarr`** path,
   - a **server-backed Viewer URL** (the exact printed
     `/?source=remote` for prepared data or `/?anndata=true` for direct
     AnnData),
   - a **remote server** (`?remote=...`, only when the viewer origin can fetch that scheme),
   - a **GitHub exports root** (`?github=...`),
   - or a **self-hosted exports root** on a static host or CDN
     (`?exportsBaseUrl=...`, which drives the `Sample datasets:` dropdown)?
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

## Read the message before anything else

Cellucid's loading failures are written to be diagnostic, and the wording tells
you which of the three layers failed. Connection failures for the remote-server
and GitHub paths take the shape:

> *`<subject>` `<what happened>`: `<cause>` `<what to do>`*

and the middle clause is the one that identifies the layer:

| The message says | It means | Look at |
|---|---|---|
| nothing is published at that address (the server answered "not found") | The address resolved, but there is nothing there | Your path, folder, or branch |
| the server refused access | Reached, but not public | Repository visibility |
| the server rejected the request | Reached, but the request was malformed | Your address |
| the server reported a problem of its own | Reached; the far side failed | Wait, then retry |
| what is published there is larger than Cellucid opens in a browser | Reached, readable, valid — and over a browser ceiling | The export size. Serve it instead ({doc}`04_server_tutorial`), or publish a smaller export. The ceilings are tabulated in {doc}`../q_troubleshooting_index/index` |
| what is published there is not in a format Cellucid can read | Reached and readable, but not a Cellucid export | Re-export with `prepare()` |
| nothing readable came back — the connection may have failed, or that address may not hold Cellucid data | Ambiguous: no usable response arrived | Your network *and* your address |

:::{warning}
Read the size row and the format row as opposites. *Larger than Cellucid opens
in a browser* means the export is perfectly valid and merely too big, so
re-exporting it produces the same bytes and meets the same wall. That row's fix
is a different loading path, not a different file.
:::

There is one message that does **not** carry a cause clause, and its punctuation
is the tell — a full stop where the others have a colon:

> *`<subject>` `<what happened>`. Try again; if it keeps failing, check the
> address and your connection. Details: `<raw error text>`*

That is what Cellucid says when it could not classify the failure at all, and
naming a cause anyway would be a guess. It is rare, because the ordinary
outcomes above are all classified — a connection you cancel yourself, and one
that runs past its timeout, both read as *nothing readable came back*. What is
left is a transfer the browser abandoned from underneath the request, or a
defect inside Cellucid. That is why the trailing `Details:` exists at all, and
why it appears in this case only: it is the raw error text, and it is the only
evidence anyone has. Quote it if you report the problem.

Local file failures use neither shape. They quote the specific defect instead —
for example `The selected file is not a valid HDF5/H5AD file…` or
`No UMAP embedding this viewer can read was found in obsm…`.

---

## Symptom → diagnosis → fix (common issues)

### Symptom: “I clicked Prepared / H5AD / Zarr ZIP, but nothing happens”

**Likely causes (ordered)**
1) **You cancelled the file chooser.** This is by far the most common cause and
   it is deliberately silent: Cellucid shows no error, starts no spinner, and
   leaves the current dataset in place.
2) You denied the file or directory permission prompt.
3) The app is embedded in a context that blocks file access.
4) A managed-browser policy disabled local file access.
5) The selected object does not match the control: directory for **Prepared**,
   `.h5ad` for **H5AD**, or one `.zarr.zip`/`.zip` for **Zarr ZIP**.

**How to confirm**
- Click the same button again and confirm the dialog. Silence after a *confirmed*
  selection is a real problem; silence after a dismissed one is expected.
- Open the standalone Cellucid page in a supported current desktop browser.
- Use a small known-current fixture with the matching control.

**Fix**
- Allow the picker permission and use the standalone page.
- If organizational policy forbids local selection, choose server mode
  explicitly: {doc}`04_server_tutorial`.

**Prevention**
- For workshops, validate the three controls with a small fixture on the
  managed browser image before distributing data.

---

### Symptom: “It loads forever / spinner never ends”

:::{note}
An `.h5ad` **above 512 MiB never gets this far**. It is refused immediately,
before a byte is read, with `H5AD direct browser files must have a positive safe
size no larger than 512 MiB; use the Cellucid server or prepared format`. An
endless spinner therefore means a file under the ceiling that is still too big
once decompressed — or one of the network causes below.
:::

**Likely causes (ordered)**
1) Large `.h5ad` under 512 MiB loaded directly in the browser (memory pressure
   from the decompressed matrix).
2) Remote server URL is wrong or unreachable.
3) GitHub raw fetch blocked (corporate firewall/content blocker).

**How to confirm**
- Try loading a demo dataset.
- If using server mode, open the health endpoint:

  ```text
  http://127.0.0.1:8765/_cellucid/health
  ```

- If using GitHub, open the raw manifest URL in a new tab after substituting
  your exact owner, repository, branch, and path:

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

The two browser readers word this differently, and only one of them tells you
what your file actually contains.

A **browser-selected `.h5ad`** reports:

```text
No UMAP embedding this viewer can read was found in obsm. Expected one or more
of X_umap_1d, X_umap_2d, or X_umap_3d, or a plain X_umap of 1, 2, or 3 columns.
Available obsm keys: <what your file has>.
```

A file with an empty `obsm` reads `Available obsm keys: (none).`

A **browser-selected Zarr ZIP** reports the same requirement without the key
list:

```text
AnnData Zarr requires a UMAP embedding in obsm: X_umap_1d, X_umap_2d, or
X_umap_3d, or a plain X_umap of 1, 2, or 3 columns
```

So for a Zarr ZIP, list the keys yourself in Python (`print(adata.obsm.keys())`)
rather than looking for them in the message.

```{figure} ../../../_static/screenshots/data_loading/fail-missing-umap-embedding.png
:alt: A Cellucid notification listing the expected X_umap_1d, X_umap_2d and X_umap_3d keys and then the single obsm key the selected file actually contained, X_pca.
:width: 760px

The trailing list is the useful half: it names the `obsm` keys your file does
have, so you can see immediately whether the embedding is missing or merely
named differently. This capture predates the current wording, which also names
the plain `X_umap` among the keys it accepts; the key list at the end is
unchanged.
```

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
  X = adata.obsm["X_umap_2d"]
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
2) The genes were named from the wrong `var` column, so the names you search for are not the ones published.
3) Data is extremely large and gene fetch is slow (looks like “nothing”).

**How to confirm**
- AnnData: `print(adata.X is None)`, `print(adata.n_vars)`, `print(adata.var.head())`
- Exports: ensure `var_manifest.json` and `var/` exist (gene expression exported).
  The files under `var/` are numbered, so read the names out of
  `var_manifest.json` rather than out of the directory listing.

**Fix**
- Provide a valid gene expression matrix, and name the genes from the `var` column your readers will search — `gene_id_column=...` in server/Jupyter mode, `var_gene_id_column=...` when exporting.
- For very large datasets, use server mode (best reliability): {doc}`04_server_tutorial`

---

### Symptom: “Vector field overlay toggle is disabled / no fields appear”

**Likely causes (ordered)**
1) You are on a Python path (`cellucid serve`, `serve_anndata()`,
   `show_anndata()`) and did not pass `--vector-fields` /
   `serve_vector_fields=True`. The overlay is **off by default** there; the
   three browser pickers need no such flag.
2) The dataset truly has no vector fields (common).
3) Vector fields exist, but they weren’t exported / served (missing `vectors/` or wrong server path).
4) Naming mismatch: vectors are present but not discoverable (e.g., not `*_umap_2d` / `*_umap_3d`).
5) Dimension mismatch: you only have 2D vectors but you’re looking at 3D (or vice versa).

**How to confirm**
- **Python server or notebook**: read step 3 of the terminal startup report, or
  check `dataset_identity.json` for a `vector_fields` block. See
  {doc}`04_server_tutorial` and {doc}`05_jupyter_tutorial`.
- **Exports**:
  - Open `<export_dir>/dataset_identity.json` and search for `vector_fields`.
  - Confirm `<export_dir>/vectors/` exists and contains `*_2d.bin(.gz)` / `*_3d.bin(.gz)` files.
- **AnnData**:
  - Print `adata.obsm.keys()` and look for `velocity_umap_2d`-style entries.

**Fix**
- Restart the server or viewer with `--vector-fields` / `serve_vector_fields=True`.
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
3) The files are not on `main`, but the shorthand omitted `@branch`.
4) Corporate network blocks `raw.githubusercontent.com`.

**How to confirm**
- Cellucid raises two notifications, and the first one names the exact address
  it requested:

  ```text
  Resource not found: https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/datasets.json
  ```

- Open that URL in a browser tab. If it 404s there, the path or branch is wrong.

```{figure} ../../../_static/screenshots/data_loading/fail-github-catalog-not-found.png
:alt: Two stacked Cellucid notifications, the upper naming the raw.githubusercontent.com URL that returned not found, the lower explaining that the GitHub repository could not be reached and asking the reader to check for a typo.
:width: 776px

Both notifications for one wrong path. Debug with the upper one; it is the exact
request that failed.
```

**Fix**
- Generate and commit `datasets.json` at the exports root.
- Confirm you’re pointing Cellucid at the exports root, not at `.../pbmc_demo/`.
- Shorthand without `@branch` always resolves to `main`; use
  `owner/repo@branch/path` for another branch.
- If blocked, use server mode or a different network.

Compare your repository with the known-good
`theislab/cellucid-demo-custom-datasets/exports` connection in
{doc}`11_custom_dataset_repository`, then see {doc}`02_local_demo_tutorial`.

---

### Symptom: “Remote server connect fails / CORS blocked / mixed content”

```{figure} ../../../_static/screenshots/data_loading/fail-remote-server-unreachable.png
:alt: A Cellucid notification reading The remote server could not be reached, followed by the explanation that nothing readable came back, that the connection may have failed or the address may not hold Cellucid data, and advice to check the network and the address.
:width: 776px

Connecting to a port with nothing listening on it. Cellucid deliberately names
*both* possibilities, because a failed connection and a wrong address are
indistinguishable from the browser's side.
```

**Likely causes**
- Wrong URL (host/port mismatch) — this and a dead server produce the same
  message.
- You pasted the **Viewer URL** into the `Remote server:` field. That field
  wants the bare origin the banner prints as **Local URL**
  (`http://127.0.0.1:8765`) and refuses a trailing slash, a `?…` query, or a
  fragment.
- Server is bound to `127.0.0.1` on a remote machine, but you’re trying to access it directly.
- The Cellucid page is HTTPS and the address is `http://`, which Cellucid
  refuses before the browser does:
  `An HTTPS Cellucid page requires an explicit HTTPS remote server URL`.

**How to confirm**
- Check server banner output and copy the exact viewer URL printed by the server.
- Visit the health endpoint — it separates the two causes the message cannot:

  ```text
  http://127.0.0.1:8765/_cellucid/health
  ```

  JSON back means the address was wrong; nothing back means the server is not
  running or not reachable.

**Fix**
- For remote machines: use an SSH tunnel (recommended): {doc}`04_server_tutorial`
- Keep the server bound to `127.0.0.1` and access via `http://localhost:<port>` through the tunnel.

**Reassurance**

A failed connection does not throw away what you already had. The dataset that
was open stays open, and the control returns to its unconnected state:

```{figure} ../../../_static/screenshots/data_loading/keep-dataset-after-failed-connect.png
:alt: The full Cellucid window after a failed connection attempt: the previously loaded dataset is still rendered in the canvas and still summarised in the Session panel, while two error notifications sit in the lower right corner.
:width: 1440px

After a refused connection the previous dataset is untouched — same points, same
summary, same fields. Only the notifications changed. You can correct the
address and try again without reloading the page.
```

---

## When to stop debugging and switch workflows (high leverage)

If you hit one of these, switching workflows saves time:

- **Large `.h5ad` in browser** → switch to server mode or exports.
- **GitHub too large / blocked** → switch to server mode or host exports on an accessible HTTP server.
- **File picker blocked** → switch to server mode.

---

## What to include when asking for help (so others can debug quickly)

Copy/paste:
- your loading path (Prepared / H5AD / Zarr ZIP / server / GitHub)
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
- Complete custom repository reference: {doc}`11_custom_dataset_repository`
- Dataset identity: {doc}`06_dataset_identity_why_it_matters`
- Vector fields overlay: {doc}`../i_vector_field_velocity/index`
