# Server (browser tab + local HTTP server)

```{eval-rst}
.. currentmodule:: cellucid
```

This page documents **server-mode** usage: open Cellucid in a browser tab while a local server serves the data.

You have three equivalent entry points:
- **CLI:** `cellucid serve …` (recommended for most users; auto-detects format)
- **Python (exported data):** {func}`~cellucid.serve` / {class}`~cellucid.CellucidServer`
- **Python (AnnData):** {func}`~cellucid.serve_anndata` / {class}`~cellucid.AnnDataServer`

---

## Fast path (beginner-friendly)

### A) Serve anything (auto-detected) from the terminal

```bash
# Exported folder, .h5ad, or .zarr are all supported:
cellucid serve /path/to/data
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
cellucid serve /path/to/data.zarr --dataset-name "My dataset" --dataset-id my-dataset
```

Then open the printed URL (or let it auto-open).

### B) Serve an exported dataset from Python

```python
from cellucid import serve

serve("./my_export")  # blocks until Ctrl+C
```

---

## Practical path (deployment patterns)

### 1) Exported dataset vs AnnData mode

**Exported dataset** (recommended for speed + reproducibility):
- Create with {func}`~cellucid.prepare` (see {doc}`export`)
- Serve with {func}`~cellucid.serve` or `cellucid serve ./my_export`

**AnnData mode** (convenient, often slower):
- Serve with {func}`~cellucid.serve_anndata` or
  `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`
- Useful during exploratory analysis when you don’t want to write an export yet

### 2) Local machine vs remote server (HPC / cloud VM)

**Local machine**:
- Keep default host `127.0.0.1` (only accessible on your machine).

**Remote server + SSH tunnel** (recommended for HPC):
1. On the remote machine (where the data lives):
   ```bash
   cellucid serve /path/to/data --host 127.0.0.1 --port 8765 --no-browser
   ```
2. On your laptop:
   ```bash
   ssh -L 8765:127.0.0.1:8765 user@remote
   ```
3. In your laptop browser:
   - Open the printed prepared-export Viewer URL:
     `http://127.0.0.1:8765/?source=remote`

**Direct LAN access** (use carefully):
- Bind to `0.0.0.0` to make the server accessible to other machines on the network:
  ```bash
  cellucid serve /path/to/data --host 0.0.0.0
  ```
- Do this only if you understand your network/security posture.

---

## Deep path (how server mode works)

### What the server provides

- Static or virtual endpoints for dataset files (exported data or AnnData-backed “virtual files”)
- CORS headers so the viewer UI can request data
- `Host` validation ahead of routing: only one well-formed authority naming
  `localhost` or a loopback IP literal on the bound port is served, which is what
  defeats DNS rebinding. Behind a reverse proxy, declare its host name with
  `--allowed-host HOST` or `allowed_hosts=[...]`; see
  {doc}`../../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`
- Health/info endpoints:
  - `/_cellucid/health`
  - `/_cellucid/info`
  - `/_cellucid/datasets`
  - `/_cellucid/protocol` (which wire capabilities this `cellucid` accepts)

### Viewer UI hosting and source requirements

By default, Cellucid establishes and serves the exact verified **viewer UI**
generation from the same server origin. This avoids:
- mixed-content issues (HTTPS notebook + HTTP localhost)
- cross-origin iframe restrictions

If the runtime cannot fetch and verify the exact source generation, startup
raises before binding the server. Use `--web-cache-dir PATH` or the
`web_cache_dir=` argument to select its publication directory.

---

## API reference

### Functions

```{eval-rst}
.. autofunction:: serve
```

```{eval-rst}
.. autofunction:: serve_anndata
```

### Classes

```{eval-rst}
.. autoclass:: CellucidServer
   :members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: AnnDataServer
   :members:
   :show-inheritance:
```

---

## Edge cases (do not skip)

### “It works locally but not on the remote server”
- If the server is on a remote machine, your browser cannot directly reach its localhost.
- Use SSH tunneling (preferred) or bind to an external interface + open firewall rules (higher risk).

### “I bound to 0.0.0.0 and now anyone can access it”
- Binding `--host 0.0.0.0` exposes the server to your network.
- Cellucid servers are designed for local/private use (not hardened for public internet exposure).

### “The dataset directory contains multiple datasets”
- {class}`~cellucid.CellucidServer` supports serving a directory that contains multiple exported datasets as subfolders.
- `server.viewer_url` opens that served catalog without embedding an arbitrary
  dataset id. One unique entry is selected automatically; multiple entries
  require an exact dataset-id selection.
- Confirm the exact ids and paths by visiting `/_cellucid/datasets`.
- Once any immediate child is recognized as a dataset candidate, **every**
  immediate subdirectory must be one complete current export. A stray or
  partial subdirectory rejects the root instead of being silently omitted.

### “Exported folder missing some files”
- `CellucidServer` validates the complete declared prepared-artifact inventory
  during construction.
- Missing required metadata, point files, or any artifact declared by a
  current manifest raises before a socket is bound.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Port already in use”
Fix:
- Use `--port 0` if supported by your workflow (otherwise choose a free port).
- Or pick another port manually: `--port 9000`.

---

### Symptom: “Data directory not found”
Fix:
- Confirm the path exists on the machine where you ran the server.
- If you’re on a remote server, remember your local machine path is different.

---

### Symptom: “I can open `/_cellucid/health` but the dataset won’t load”
Likely causes:
- You served a directory that is not a valid exported dataset (missing manifests/points files).
- You served a parent folder with multiple datasets but opened the wrong path.

How to confirm:
- Visit `/_cellucid/datasets` and confirm the dataset path you intended exists.
- Try fetching `dataset_identity.json` and `obs_manifest.json` in the browser.

Fix:
- Serve the correct dataset directory (the folder produced by `prepare(...)`).
- Or re-export with {func}`~cellucid.prepare`.

---

### Symptom: “CORS blocked” (browser console error)
Likely causes:
- A corporate environment / browser extension blocks cross-origin requests.

Fix:
- Prefer opening the viewer URL served by the same origin (default).
- Disable the conflicting extension for the site, or verify the request in a
  clean profile with extensions disabled.

---

### Symptom: web-generation startup failure
Likely causes:
- The configured source inventory or one of its declared objects could not be
  fetched or verified.

Fix:
- Ensure source access at startup and pass a writable `web_cache_dir`.

---

### Symptom: “Browser can’t connect / connection refused”
Fix:
- Confirm server is running and you are opening the correct URL.
- If remote: confirm your SSH tunnel is active and mapped to the correct port.

---

### Symptom: “It works for small datasets but crashes/gets slow for large ones”
Likely causes:
- AnnData mode is generating virtual files on demand and may hit memory/CPU limits.

Fix:
- Export with {func}`~cellucid.prepare` (quantization + compression) and serve the exported directory.
- Use read-only-backed `.h5ad` input when staying in AnnData mode; Zarr input
  is loaded eagerly.

---

### Symptom: direct AnnData startup says “Importing dependencies…” and feels stuck
What’s happening:
- The first import may be slow in fresh environments because large scientific dependencies are loaded.

Fix:
- Wait once (subsequent runs are faster in the same process).
- If it truly hangs, run with `--verbose` and check for import errors.

---

## See also

- {doc}`cli` for CLI options and examples
- {doc}`export` for `prepare(...)` and export layout
- {doc}`jupyter` for notebook embedding
