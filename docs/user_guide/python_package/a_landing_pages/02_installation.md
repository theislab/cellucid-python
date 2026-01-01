# Installation

**Audience:** everyone  
**Time:** 5–20 minutes  
**Goal:** install `cellucid`, verify the CLI, and avoid the most common environment pitfalls

## Requirements

- **Python:** `>= 3.10` (see `requires-python` in `pyproject.toml`)
- **OS:** macOS / Linux / Windows (development is primarily macOS/Linux; please report Windows issues)
- **Browser:** a modern Chromium-based browser or Firefox (for the web app UI)
- **Network (for first run):** the viewer UI is fetched from `https://www.cellucid.com` and cached locally  
  - If you’re offline the very first time, you can still serve data, but the embedded UI may show a “Viewer UI unavailable” page until cached.

```{note}
Your **data is not uploaded** to `cellucid.com` by the Python package. The network requirement is for fetching the **viewer UI assets** (HTML/JS/CSS) the first time, unless you already have them cached.
```

## Fast path (recommended)

### 1) Create a clean environment

Pick one of these. If you’re not sure, use **venv**.

**Option A: `venv` (works everywhere)**

```bash
python -m venv .venv
source .venv/bin/activate  # macOS/Linux
python -m pip install -U pip
```

**Option B: conda/mamba**

```bash
conda create -n cellucid python=3.11 -y
conda activate cellucid
python -m pip install -U pip
```

### 2) Install Cellucid (Python)

```bash
pip install cellucid
```

### 3) Verify the installation

```bash
python -c "import cellucid; print(cellucid.__version__)"
cellucid --version
cellucid serve --help
```

If these commands run without errors, your install is working.

## Practical path (computational users)

### Upgrading / pinning

Upgrade:

```bash
pip install -U cellucid
```

Pin a specific version (recommended for paper pipelines):

```bash
pip install "cellucid==0.0.1a2"
```

### Optional dependencies (what you might need next)

Cellucid installs core scientific dependencies, but common workflows often require additional packages:

- **Jupyter notebooks** (recommended for `show_anndata(...)` / hooks):  
  `pip install jupyterlab` (or `pip install notebook`)
- **Remote/HTTPS notebooks** (JupyterHub, some corporate setups):  
  `pip install jupyter-server-proxy` (helps avoid HTTPS→HTTP mixed-content blocking)
- **Zarr-backed AnnData** (`.zarr` stores):  
  `pip install zarr`
- **Single-cell analysis stack** (compute embeddings / neighbors / QC):  
  `pip install scanpy` (and whatever your lab uses)

### (Recommended) Prefetch the viewer UI cache

If you are:
- behind a firewall,
- frequently offline (air-gapped machines),
- or teaching/workshopping and want predictable first-load behavior,

prefetch the viewer UI assets once while online.

1) Check where the cache will live:

```bash
python -c "from cellucid import get_web_cache_dir; print(get_web_cache_dir())"
```

2) Prefetch the UI (best-effort):

```bash
python -c "from cellucid.web_cache import ensure_web_ui_cached; ensure_web_ui_cached()"
```

3) (Optional) Make the cache persistent

By default, Cellucid may use a temp directory for cached web assets. To make the cache persist across reboots, set:

```bash
export CELLUCID_WEB_PROXY_CACHE_DIR=\"$HOME/.cache/cellucid-web\"
```

On Windows (PowerShell), you can set it per-session like:

```powershell
$env:CELLUCID_WEB_PROXY_CACHE_DIR = \"$env:USERPROFILE\\.cache\\cellucid-web\"
```

To clear the cache:

```bash
python -c "from cellucid import clear_web_cache; print(clear_web_cache())"
```

## Deep path (developers/maintainers)

- Editable install (repo development):
  ```bash
  pip install -e ".[dev,docs]"
  ```
- Docs build dependencies (Sphinx/MyST):
  ```bash
  pip install -e ".[docs]"
  ```

## Edge cases and common footguns

- **Multiple Python environments**: “installed but can’t import” usually means you installed into one env and ran Python in another.
- **Conda + pip mixing**: it can work, but if you hit binary/ABI issues, prefer either “all conda” or “all pip” for core scientific deps.
- **HTTPS notebook pages**: if your notebook is served over HTTPS, the browser can block `http://127.0.0.1:<port>` iframes (mixed content). See {doc}`03_compatibility_matrix_must_be_explicit`.
- **First-run offline**: the viewer UI might not render until it has been cached once (see prefetch instructions above).

## Troubleshooting (installation)

### Symptom: `pip install cellucid` fails with a Python version error

**Likely causes**
- You’re on Python `< 3.10`.

**How to confirm**
```bash
python --version
```

**Fix**
- Create a new environment with Python 3.10+ (see “Fast path”).

---

### Symptom: `ModuleNotFoundError: No module named 'cellucid'`

**Likely causes**
- Installed into a different environment than the one you’re running.

**How to confirm**
```bash
which python
python -m pip show cellucid
```

**Fix**
- Activate the intended environment, then reinstall:
  ```bash
  pip install -U cellucid
  ```

---

### Symptom: `cellucid: command not found`

**Likely causes**
- The environment’s `bin/` directory is not on your PATH.
- You installed into a venv but didn’t activate it.

**Fix**
- Activate your environment (`source .venv/bin/activate`), then retry.
- Or run via `python -m cellucid.cli --help` (advanced / debugging).

---

### Symptom: the notebook viewer shows “Cellucid viewer UI could not be loaded”

**Likely causes**
- No network access to `https://www.cellucid.com` and no cached copy exists.
- Cache directory is not writable.

**How to confirm**
- Run:
  ```bash
  python -c "from cellucid import get_web_cache_dir; print(get_web_cache_dir())"
  ```

**Fix**
- Get online once and prefetch the UI cache:
  ```bash
  python -c "from cellucid.web_cache import ensure_web_ui_cached; ensure_web_ui_cached()"
  ```
- If you suspect a permissions issue, set `CELLUCID_WEB_PROXY_CACHE_DIR` to a writable location and retry.

---

### Symptom: the embedded viewer is blank in JupyterLab / JupyterHub

**Likely causes**
- Your notebook is served from HTTPS or a remote origin, so the browser blocks an HTTP loopback iframe.

**Fix**
- Install `jupyter-server-proxy` and restart the notebook server:
  ```bash
  pip install jupyter-server-proxy
  ```
- Or use the browser workflow instead (`cellucid serve ...`) and open the URL directly.

---

### Symptom: browser opens but shows “Port already in use”

**Likely causes**
- Another process is using the port (default: 8765).

**Fix**
- Pick a different port:
  ```bash
  cellucid serve /path/to/data --port 9000
  ```
- Or stop old viewers/servers (in notebooks, use `viewer.stop()` or restart the kernel).

## Next steps

- Environment gotchas first: {doc}`03_compatibility_matrix_must_be_explicit`
- Then get something running end-to-end: {doc}`04_quick_start_3_levels`
