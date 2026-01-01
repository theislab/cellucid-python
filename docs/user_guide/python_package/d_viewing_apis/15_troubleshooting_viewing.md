# Troubleshooting (viewing)

This page is symptom-driven troubleshooting for:
- `cellucid serve …` (CLI server mode),
- `serve(...)` / `serve_anndata(...)` (Python server mode),
- `show(...)` / `show_anndata(...)` (notebook embedding).

If you’re new, skim the “Quick triage” section first.

---

## Quick triage (do these in order)

### 1) What workflow are you using?

- Terminal/browser → you’re in **server mode** (`cellucid serve …`)
- Notebook iframe → you’re in **notebook mode** (`show` / `show_anndata`)

This matters because notebook mode adds extra failure modes (HTTPS mixed content, notebook proxies).

### 2) Is the server alive?

Open (replace `<port>`):

```text
http://127.0.0.1:<port>/_cellucid/health
```

If this doesn’t load, the viewer can’t load either.

### 3) Is the viewer UI available (hosted-asset proxy)?

Open:

```text
http://127.0.0.1:<port>/
```

If you see a “viewer UI could not be loaded” page, it’s a UI cache/network issue (see below).

### 4) Is the dataset metadata reachable?

Open:

```text
http://127.0.0.1:<port>/dataset_identity.json
```

If this fails:
- exported folders may be incomplete, or
- AnnData server may have crashed while generating metadata.

### 5) If you are in a notebook: run a structured debug report

```python
viewer.debug_connection()
```

---

## Symptom: “Port already in use” / server won’t start

### Likely causes (ordered)
- another Cellucid server is still running
- another program is using the port (often `8765`)

### How to confirm
- the CLI prints “Port X in use, using Y”
- `/_cellucid/health` works on a different port than you expected

### Fix
- use the printed Viewer URL (don’t assume `8765`)
- choose a port explicitly: `cellucid serve ... --port 9000`
- in notebooks, call `viewer.stop()` on old viewers (see {doc}`11_viewer_lifecycle_cleanup_ports_and_multiple_viewers`)

### Prevention
- stop servers when done
- use fixed ports only when you need them (SSH tunnels)

---

## Symptom: “Cellucid viewer UI could not be loaded”

This is the hosted-asset proxy error page.

### Likely causes (ordered)
- you are offline and have no cached UI assets yet
- your environment cannot reach `https://www.cellucid.com` (DNS/firewall/proxy)
- the cache directory is not writable or is being cleared (ephemeral temp dirs)

### How to confirm
- open the server root (`http://127.0.0.1:<port>/`) and read the reason text
- in notebooks: `viewer.debug_connection()` → look at `web_ui.cache` and `web_ui.prefetch`

### Fix
- run once while online (starting a server triggers a best-effort UI prefetch)
- set a persistent cache dir:

```bash
export CELLUCID_WEB_PROXY_CACHE_DIR=$HOME/.cache/cellucid-web
```

- clear and re-download if you suspect a corrupted cache:

```python
from cellucid import clear_web_cache
clear_web_cache()
```

### Prevention
- keep the cache in a persistent, writable directory
- if you plan to work offline later, run `cellucid serve ...` once while online ahead of time

<!-- SCREENSHOT PLACEHOLDER
ID: viewer-ui-unavailable-page
Suggested filename: web_app/viewing_01_ui-unavailable.png
Where it appears: Python Package Guide → Viewing APIs → Troubleshooting → UI unavailable
Capture:
  - UI location: browser tab at `http://127.0.0.1:<port>/`
  - State prerequisites: UI fetch fails and no cached copy exists
  - Action to reach state: set cache dir to empty path, disconnect internet, start server, open viewer URL
Crop:
  - Include: the error title + the “What to do” bullet list
  - Exclude: browser bookmarks/personal info
Alt text:
  - Cellucid error page stating the viewer UI could not be loaded.
Caption:
  - When the server cannot fetch the viewer UI and no cached copy exists, it serves an explanatory error page with next steps (network + cache directory).
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the “viewer UI unavailable” error page.
:width: 100%

Hosted-asset proxy error page shown when the viewer UI cannot be fetched and no cached copy exists.
```

---

## Symptom: Viewer loads but dataset is blank / stuck loading

### Likely causes (ordered)
- export directory is missing required files (manifests or points binaries)
- you passed the wrong directory (parent folder vs dataset folder)
- AnnData mode: missing UMAP keys or shape mismatch

### How to confirm
- open `/_cellucid/datasets` and see what it lists
- open `/dataset_identity.json` and check:
  - `stats.n_cells`
  - `embeddings.available_dimensions`
- in AnnData mode: check `adata.obsm.keys()` and shapes

### Fix
- for exports: validate your folder layout (see {doc}`07_exported_directory_mode_show_and_serve`)
- for AnnData: add explicit `X_umap_2d` / `X_umap_3d` keys (see {doc}`08_anndata_mode_show_anndata_and_serve_anndata`)

---

## Symptom: “No valid UMAP embeddings found in adata.obsm”

### Likely causes (ordered)
- UMAP was never computed
- keys exist but have wrong shape
- the embedding array has a row mismatch vs `adata.n_obs`

### How to confirm

```python
print(list(adata.obsm.keys()))
for k in ["X_umap_1d", "X_umap_2d", "X_umap_3d", "X_umap"]:
    if k in adata.obsm:
        print(k, getattr(adata.obsm[k], "shape", None))
print("n_cells:", adata.n_obs)
```

### Fix
- compute UMAP and store explicit dimension keys (2D/3D)

---

## Symptom: Gene search returns nothing / “Gene 'X' not found”

### Likely causes (ordered)
- your gene IDs are not in `var.index` (default)
- wrong `gene_id_column`
- duplicate gene IDs (ambiguous)

### How to confirm

```python
print(adata.var.index[:5])
print(adata.var.columns)
```

### Fix
- pass `gene_id_column="..."` to `show_anndata(...)` / `serve_anndata(...)`
- ensure the chosen gene ID field is unique

---

## Symptom (notebook): iframe shows “proxy required” / mixed-content type errors

### Likely causes (ordered)
- notebook served over HTTPS
- kernel/server is remote (JupyterHub/cloud)
- `jupyter-server-proxy` is missing/disabled

### How to confirm
- the iframe message explicitly suggests `jupyter-server-proxy`
- `viewer.debug_connection()` reports `jupyter_server_proxy.installed: False`

### Fix
- install/enable `jupyter-server-proxy` (recommended)
- or set `CELLUCID_CLIENT_SERVER_URL` to an HTTPS-reachable URL (advanced)

Deep dive: {doc}`10_notebook_widget_mode_advanced`.

---

## Symptom: Vector field / velocity overlay is missing

### Likely causes (ordered)
- vectors not present (exports: missing `vector_fields` metadata; AnnData: missing `obsm` keys)
- naming mismatch (e.g. `velocity_umap2d` vs `velocity_umap_2d`)
- dimension mismatch (you have 2D vectors but are viewing 3D)

### How to confirm
- exports: open `/dataset_identity.json` and search for `vector_fields`
- AnnData: inspect `adata.obsm.keys()` for `*_umap_2d` / `*_umap_3d`

### Fix
- rename keys to the supported convention
- ensure you have vectors for the dimension you’re viewing

---

## Symptom: Performance is extremely slow / browser crashes

### Likely causes (ordered)
- using browser `.h5ad` loading instead of Python-backed mode
- forcing in-memory mode (`--no-backed`) on a large dataset
- huge categorical fields or heavy overlays

### Fix
- use `cellucid serve data.h5ad` (backed) or `data.zarr`
- export once (`prepare`) and use export mode
- reduce or remove pathological fields (very high cardinality)

Performance guide: {doc}`14_performance_scaling_and_lazy_loading`.

---

## Reporting a bug (what to include)

Please capture:

- the exact command you ran (`cellucid serve ...`), including flags
- whether the input was export vs `.h5ad` vs `.zarr`
- the output of:
  - `/_cellucid/health`
  - `/_cellucid/info`
- in notebooks: `viewer.debug_connection()` output
- if safe: `dataset_identity.json` (or at least its `stats` block)

This will save you (and maintainers) a lot of back-and-forth.
