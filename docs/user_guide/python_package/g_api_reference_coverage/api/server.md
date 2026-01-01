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
cellucid serve /path/to/data.h5ad
cellucid serve /path/to/data.zarr
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
- Serve with {func}`~cellucid.serve_anndata` or `cellucid serve data.h5ad`
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
   - Open `http://127.0.0.1:8765/`

**Direct LAN access** (use carefully):
- Bind to `0.0.0.0` to make the server accessible to other machines on the network:
  ```bash
  cellucid serve /path/to/data --host 0.0.0.0
  ```
- Do this only if you understand your network/security posture.

---

## Screenshot placeholders (optional but recommended)

<!-- SCREENSHOT PLACEHOLDER
ID: server-terminal-banner
Where it appears: Server (browser tab + local HTTP server) → Fast path
Capture:
  - UI location: terminal output after starting the server
  - State prerequisites: server started successfully
  - Action to reach state: run `cellucid serve ./my_export` and wait for the banner
Crop:
  - Include: the final “SERVER RUNNING” banner + the Local URL/Viewer URL lines
  - Exclude: full terminal prompt with username/hostname, any private paths
Redact:
  - Remove: usernames, hostnames, private dataset paths
Annotations:
  - Callouts:
    - (1) Local URL (API/data origin)
    - (2) Viewer URL (what to open in browser)
Alt text:
  - Terminal output showing the Cellucid server banner with local and viewer URLs.
Caption:
  - The server prints a local URL (data origin) and a viewer URL; open the viewer URL in your browser to load the dataset.
-->
```{figure} ../../../../_static/screenshots/placeholder-screenshot.svg
:alt: Terminal output showing the Cellucid server banner with local and viewer URLs.
:width: 100%

The server prints a local URL (data origin) and a viewer URL; open the viewer URL in your browser to load the dataset.
```

---

## Deep path (how server mode works)

### What the server provides

- Static or virtual endpoints for dataset files (exported data or AnnData-backed “virtual files”)
- CORS headers so the viewer UI can request data
- Health/info endpoints:
  - `/_cellucid/health`
  - `/_cellucid/info`
  - `/_cellucid/datasets`

### Viewer UI hosting and offline behavior

By default, Cellucid serves the **viewer UI** from the same server origin via a hosted-asset proxy. This avoids:
- mixed-content issues (HTTPS notebook + HTTP localhost)
- cross-origin iframe restrictions

If the runtime cannot reach the hosted viewer assets, you may see a “viewer UI unavailable” error page.
Use the environment variable:
- `CELLUCID_WEB_PROXY_CACHE_DIR` to control where the UI cache is stored.

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
- Confirm by visiting `/_cellucid/datasets`.

### “Exported folder missing some files”
- The viewer expects a minimum set of files (`obs_manifest.json` and at least one `points_*d.bin(.gz)`).
- If those are missing, the server may still start, but the viewer will fail to load the dataset correctly.

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
- Disable conflicting extensions for the session, or use a different browser profile.

---

### Symptom: “Viewer UI unavailable” / blank page
Likely causes:
- The hosted viewer assets could not be fetched and are not cached.

Fix:
- Ensure HTTPS access to the hosted UI once (to populate cache), or set a persistent `CELLUCID_WEB_PROXY_CACHE_DIR`.

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
- Use `.h5ad` backed mode (default) or `.zarr` for lazy loading if staying in AnnData mode.

---

### Symptom: “`cellucid serve data.h5ad` says ‘Importing dependencies…’ and feels stuck”
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
