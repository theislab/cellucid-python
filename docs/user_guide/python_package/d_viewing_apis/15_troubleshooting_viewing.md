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

### 3) Is the exact viewer generation being served?

Open the exact **Viewer URL** printed by the server:

```text
Prepared export: http://127.0.0.1:<port>/?source=remote
Direct AnnData:  http://127.0.0.1:<port>/?anndata=true
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
- startup raises an address-in-use `OSError`

### Fix
- choose a port explicitly: `cellucid serve ... --port 9000`
- in Python, pass `port=0` to request one operating-system-assigned port and
  then use the actual `viewer.server_url`
- in notebooks, call `viewer.stop()` on old viewers (see {doc}`11_viewer_lifecycle_cleanup_ports_and_multiple_viewers`)

### Prevention
- stop servers when done
- use fixed ports only when you need them (SSH tunnels)

---

## Symptom: viewer startup cannot establish the web generation

This failure occurs before the viewer server is published.

### Likely causes (ordered)
- your environment cannot reach the configured `web_source_url`
- the inventory or a declared object has an invalid MIME type, length, hash,
  file set, or build identity
- the selected `web_cache_dir` is not writable

### How to confirm
- read the exception raised by `cellucid serve`, `show(...)`, or
  `show_anndata(...)`
- exercise the same establishment contract directly with
  `ensure_web_ui_cached(force=True)`

### Fix
- allow access to the configured source and correct the exact response failure
- pass a writable generation location:

```bash
cellucid serve /path/to/data --web-cache-dir /path/to/cache
```

### Prevention
- verify source access and destination permissions before a workshop or job
- do not plan on a previous local generation being substituted when the source
  is unavailable

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
for k in ["X_umap_1d", "X_umap_2d", "X_umap_3d"]:
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
- the chosen gene IDs are duplicated or have portable filename collisions;
  AnnData construction rejects them

### How to confirm

```python
print(adata.var.index[:5])
print(adata.var.columns)
```

### Fix
- pass `gene_id_column="..."` to `show_anndata(...)` / `serve_anndata(...)`
- ensure the chosen gene ID field is unique

---

## Symptom (notebook): iframe is blank, unreachable, or mixed-content blocked

### Likely causes (ordered)
- notebook served over HTTPS
- kernel/server is remote (JupyterHub/cloud)
- the browser-facing proxy/tunnel base was not supplied as `client_server_url=`

### How to confirm
- compare `viewer.viewer_url` with the URL that the browser can actually reach
- inspect the browser console for mixed-content or connection errors

### Fix
- expose the selected port through an HTTPS proxy or tunnel
- pass that exact browser-reachable base as `client_server_url=` when creating
  the viewer

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
- loading a large in-memory AnnData or Zarr store
- huge categorical fields or heavy overlays

### Fix
- use `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`
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
