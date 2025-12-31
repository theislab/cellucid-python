# Server Mode (CLI + Python) — Recommended for Large Datasets

Server Mode is the most reliable way to view **large** `.h5ad` / `.zarr` datasets in Cellucid.

You run a small local Python server that:
- reads your data efficiently (often lazily)
- serves only what the viewer needs, on demand
- adds the right CORS headers so the web viewer can fetch from it

Then you open Cellucid with a `?remote=...` URL.

This tutorial covers options **#6–#11** from the “14 loading options” list.

## At A Glance

**Audience**
- Wet lab / beginner: follow the copy/paste commands; you do not need to write code.
- Computational users: focus on lazy loading (`backed` mode), memory, and SSH tunnel workflows.

**Time**
- Local server (same machine): ~5–10 minutes
- Remote server (SSH tunnel): ~15–30 minutes

**Prerequisites**
- `pip install cellucid`
- Your data in one of these forms:
  - pre-exported folder from `prepare()`
  - `.h5ad` file
  - `.zarr` directory

## Security Model (Read Once)

- By default, the server binds to `127.0.0.1` (localhost).
  - This means **only your machine** can access it.
- If you bind to `0.0.0.0`, other machines on your network may be able to access it.
  - Do this only if you understand the risk and you trust your network.

**Best practice for remote machines:**
- Keep the server bound to `127.0.0.1` on the remote machine.
- Use an **SSH tunnel** so you still access it via `http://localhost:...` on your laptop.

This is safer and also avoids browser mixed-content issues.

## Fast Path (CLI)

1) Start the server:

```bash
cellucid serve /path/to/data.h5ad
```

2) Open the viewer:

```text
https://www.cellucid.com?remote=http://127.0.0.1:8765
```

3) Keep the terminal running while you use the viewer.
4) Stop the server with **Ctrl+C**.

<!-- SCREENSHOT PLACEHOLDER
ID: data-loading-server-terminal-banner
Suggested filename: data_loading/08_server-terminal-banner.png
Where it appears: Data Loading → 04_server_tutorial → Fast Path
Capture:
  - UI location: terminal window
  - State prerequisites: server started successfully
  - Action to reach state: run `cellucid serve ...`
Crop:
  - Include: the printed server URL and the suggested viewer URL
  - Exclude: your username/hostnames if sensitive
Redact:
  - Replace: file paths with generic `/path/to/data.h5ad`
Alt text:
  - Terminal output showing the Cellucid server URL.
Caption:
  - Explain that the viewer URL is what you open in the browser.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Cellucid server banner.
:width: 100%

When the server starts, it prints a URL you open in the browser.
```

## Option #6 — Serve a Pre-exported Folder (Best Performance)

Use this if you already ran `prepare()`.

```bash
cellucid serve /path/to/export_dir
```

Notes:
- This is usually the fastest experience.
- If the directory contains `dataset_identity.json`, the CLI will detect it as an export.

## Option #7/#8 — Serve `.h5ad` or `.zarr` Directly (Auto-detected)

```bash
# h5ad
cellucid serve /path/to/data.h5ad

# zarr (directory)
cellucid serve /path/to/data.zarr
```

Why this is good:
- For `.h5ad`, Cellucid can use **backed mode** (true lazy loading) instead of loading everything.
- For `.zarr`, chunked storage works well with lazy access.

If you see slow loads, confirm you are not forcing in-memory mode.

### Vector fields (velocity/drift) in server mode

Server mode supports the vector field / velocity overlay (if your dataset includes vectors).

Where vectors come from depends on your data:
- **Pre-exported folders**: vectors live in `vectors/`, and `dataset_identity.json` contains a `vector_fields` block.
- **AnnData (`.h5ad` / `.zarr`)**: vectors are discovered from `obsm` keys like `velocity_umap_2d`, `velocity_umap_3d`, `T_fwd_umap_2d`, etc.

Quick verification (high signal):
- Open `http://127.0.0.1:8765/dataset_identity.json` and search for `vector_fields`.

If the overlay toggle is disabled or the dropdown is empty, it’s usually:
- naming mismatch, or
- dimension mismatch (e.g. you only have 2D vectors but you’re viewing 3D).

See:
- {doc}`../i_vector_field_velocity/index`
- `cellucid/markdown/VECTOR_FIELD_OVERLAY_CONVENTIONS.md`

## CLI Options (What They Mean)

Run this anytime:

```bash
cellucid serve --help
```

Key flags:

- `--port, -p`:
  - Change port if `8765` is in use.

- `--host, -H`:
  - Default is `127.0.0.1` (local only).
  - Use `0.0.0.0` only if you need LAN access.

- `--no-browser`:
  - Don’t auto-open a browser tab.

- `--no-backed`:
  - Forces loading the entire AnnData into memory.
  - This is rarely what you want for large datasets.

- `--latent-key`:
  - Selects which `obsm` key to use as “latent space” (used for certain derived values).
  - If you don’t know, leave it alone.

## Option #9/#10/#11 — Python API Server Mode

You can start servers from Python (useful for scripted workflows).

- `serve(export_dir)` serves a pre-exported folder.
- `serve_anndata(data)` serves `.h5ad`, `.zarr`, or an in-memory `AnnData`.

```python
# Serve a pre-exported folder
from cellucid import serve

# serve("/path/to/export_dir", port=8765, host="127.0.0.1", open_browser=True)
```

```python
# Serve a .h5ad, .zarr, or AnnData
from cellucid import serve_anndata

# serve_anndata("/path/to/data.h5ad", port=8765, host="127.0.0.1", open_browser=True)
# serve_anndata("/path/to/data.zarr", port=8765, host="127.0.0.1", open_browser=True)
```

### Stopping the server from Python

If you use the class-based API (advanced), you can stop the server programmatically:

```python
from cellucid.anndata_server import AnnDataServer

server = AnnDataServer("data.h5ad", open_browser=False)
server.start_background()

# ... interact in the browser ...

server.stop()
```

For most users, **Ctrl+C** in the terminal is simplest.

## Remote Server Access (SSH Tunnel Workflow)

This is the recommended way to use Cellucid when your data lives on a remote machine (HPC, lab server).

### Step 1 — Start the server on the remote machine

On the remote machine:

```bash
cellucid serve /path/to/data.h5ad --no-browser
```

Keep it bound to `127.0.0.1` (default).

### Step 2 — Create an SSH tunnel from your laptop

On your laptop:

```bash
ssh -L 8765:localhost:8765 user@remote-host
```

Leave that SSH session open.

### Step 3 — Open Cellucid locally

```text
https://www.cellucid.com?remote=http://localhost:8765
```

Why this works well:
- Your browser talks only to `localhost`.
- You avoid exposing the server to the public internet.
- You avoid mixed-content blocks for non-localhost HTTP.

## Edge Cases

- **Port already in use**:
  - The server may automatically pick a new port.
  - Always copy the printed URL.

- **Windows firewall prompt**:
  - If you allow public network access accidentally, other machines may reach your server.

- **Large `.h5ad` in memory mode** (`--no-backed`):
  - Can cause huge RAM use.

- **Mixed content**:
  - `https://www.cellucid.com` loading from `http://localhost:...` is usually OK.
  - Loading from `http://remote-host:...` may be blocked by the browser.
  - Prefer an SSH tunnel.

- **Vector fields only exist in one dimension**:
  - If you have `*_umap_2d` vectors but no `*_umap_3d`, the overlay dropdown will be empty in 3D.
  - Switch to 2D, or provide the missing dimension.

## Troubleshooting (Massive)

### Symptom: “Port already in use”

**Likely causes**
- Another service is using that port.

**How to confirm**
- Try a different port: `--port 9000`.

**Fix**
- Run:

  ```bash
  cellucid serve /path/to/data.h5ad --port 9000
  ```

- Then open:

  ```text
  https://www.cellucid.com?remote=http://127.0.0.1:9000
  ```

---

### Symptom: “Cellucid says it can’t connect to the remote server”

**Likely causes**
1) You typed the wrong URL (port mismatch).
2) The server is bound to `127.0.0.1` but you are trying to access it from another machine.
3) Your browser blocked mixed-content requests.

**How to confirm**
- Open the server URL directly in your browser:

  ```text
  http://127.0.0.1:8765/health
  ```

  (If the server is running locally, you should get a small JSON response.)

**Fix**
- Use an SSH tunnel for remote machines.
- Ensure the server is actually running and the URL matches.

---

### Symptom: “It connects, but genes are missing / gene search returns nothing”

**Likely causes**
- Your AnnData has no expression matrix (`adata.X` empty) or no `var`.
- Your gene IDs are stored under a different `var` column.

**How to confirm**
- In Python: `print(adata.X is None)`, `print(adata.var.head())`.

**Fix**
- Use `show_anndata(..., gene_id_column="...")` in Jupyter.
- Or export via `prepare(var_gene_id_column="...")`.

---

### Symptom: “Vector field overlay toggle is disabled / no fields appear”

**Likely causes (ordered)**
1) The dataset contains no vectors (common).
2) Vector fields exist, but are not named using the expected convention (`*_umap_2d`, `*_umap_3d`).
3) Dimension mismatch: vectors exist for 2D but you’re viewing 3D (or vice versa).

**How to confirm**
- Open `http://127.0.0.1:8765/dataset_identity.json` and check for a `vector_fields` block.
- If using AnnData, list `obsm` keys and look for `velocity_umap_2d`-style entries.

**Fix**
- Rename/regenerate vector fields to follow the convention in `cellucid/markdown/VECTOR_FIELD_OVERLAY_CONVENTIONS.md`.
- Switch the viewer to the dimension that has vectors.

For overlay UI behavior and deeper debugging, see:
- {doc}`../i_vector_field_velocity/index`
- {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

---

### Symptom: “It’s extremely slow”

**Likely causes**
- Large dataset + dense X
- Running over a high-latency SSH tunnel

**Fix**
- Prefer pre-exported data for best performance.
- If using SSH, run the server close to the data (on the same machine as the file) and tunnel to it.

## Next Steps

- Dataset identity (sessions/sharing): {doc}`06_dataset_identity_why_it_matters`
- Full troubleshooting matrix: {doc}`08_troubleshooting_data_loading`
- Want notebook embedding + programmatic control? → {doc}`05_jupyter_tutorial`
- Want browser-only loading without any server? → {doc}`03_browser_file_picker_tutorial`
