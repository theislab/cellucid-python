# Installation

**Audience:** everyone  
**Time:** 5–20 minutes  
**Goal:** install `cellucid`, verify the CLI, and avoid the most common environment pitfalls

:::{note}
**Active package version: 0.9.1.** Version 0.9.1 is the PyPI submission
release. If the [PyPI package page](https://pypi.org/project/cellucid/) lists
0.9.1, the normal install and pin commands below retrieve it. Until then, they <!-- CELLUCID_VERSION -->
resolve only to versions available from PyPI; install the active source
straight from the official GitHub repository (see **From the official source
repository** below) when validating 0.9.1 before publication. <!-- CELLUCID_VERSION -->
:::

## Requirements

- **Python:** 3.11, 3.12, 3.13, or 3.14 (see `requires-python` in `pyproject.toml`)
- **OS:** macOS, Linux, or Windows
- **Browser:** a current Chrome, Edge, Firefox, or Safari release
- **Network (at viewer startup):** the Python server establishes the exact web
  generation published by `https://www.cellucid.com` before it starts serving
  the viewer. The configured source must be reachable for every
  `cellucid serve`, `show(...)`, or `show_anndata(...)` startup that serves the
  web UI.

```{note}
Your **data is not uploaded** to `cellucid.com` by the Python package. The
network request retrieves only the declared viewer files (HTML, JavaScript,
CSS, icons, and related static assets). Cellucid verifies the source inventory,
byte lengths, hashes, MIME types, build identity, and complete file set before
publishing that generation locally. An existing cached generation is not used
as a substitute when the configured source is unavailable.
```

## Fast path (recommended)

### 1) Create a clean environment

Pick one of these. If you’re not sure, use **venv**.

**Option A: `venv` (works everywhere)**

macOS/Linux:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
```

Windows PowerShell:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
py -m pip install -U pip
```

Windows Command Prompt:

```batch
py -m venv .venv
.\.venv\Scripts\activate.bat
py -m pip install -U pip
```

**Option B: conda/mamba**

```bash
conda create -n cellucid python=3.11 -y
conda activate cellucid
python -m pip install -U pip
```

### 2) Install Cellucid (Python)

From PyPI:

```bash
pip install cellucid
```

#### From the official source repository

The [PyPI package page](https://pypi.org/project/cellucid/) is authoritative for
registry availability. Until it lists the active version, `pip install cellucid`
resolves only to the versions PyPI already holds, and those are earlier
generations whose API differs from the one this guide describes. To install the
active source directly, without cloning:

```bash
pip install "cellucid @ git+https://github.com/theislab/cellucid-python"
```

This is the same one-liner the repository `README.md` gives. Clone the
repository instead only when you intend to edit it — see **Deep path** for the
editable install.

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
pip install "cellucid==0.9.1"  # CELLUCID_VERSION
```

### Core scientific runtime

The core installation includes the coupled direct-data runtime:

- `anndata>=0.12.19,<0.13`
- `zarr>=3.1.4,<4`
- `numcodecs>=0.16.3,<0.17`

Python `.zarr` input is one complete Zarr v2 root containing both `.zgroup` and
`.zattrs`. The core installation supplies the complete reader runtime,
including the compiled numcodecs wheel used on Python 3.14.

### Additional tools you might need

- **Jupyter notebooks** (recommended for `show_anndata(...)` / hooks):  
  `pip install jupyterlab` (or `pip install notebook`)
- **Single-cell analysis stack** (compute embeddings / neighbors / QC):  
  `pip install scanpy` (and whatever your lab uses)

Cellucid already installs `jupyter-server-proxy`. For remote/HTTPS notebooks,
the Jupyter deployment must configure a route for the selected Cellucid port,
and the caller must pass that exact browser base as `client_server_url=`.
Package installation alone does not configure routing.

### Verify viewer-source access before a workshop or deployment

`ensure_web_ui_cached()` exercises the same exact source-generation contract as
viewer startup. Run it ahead of time to confirm that the configured source is
reachable and the selected cache directory is writable. It does not enable
later offline startup.

1) Check where the cache will live:

```bash
python -c "from cellucid import get_web_cache_dir; print(get_web_cache_dir())"
```

2) Establish and verify the source generation:

```bash
python -c "from cellucid.web_cache import ensure_web_ui_cached; print(ensure_web_ui_cached())"
```

3) (Optional) Select an explicit cache directory

The default cache is the process-independent temporary directory reported by
`get_web_cache_dir()`. For a different location, pass
`--web-cache-dir PATH` to `cellucid serve`, or pass `web_cache_dir=PATH` to the
Python viewing API. The selected path must be a real directory or a missing
path; symbolic links and Windows reparse points are rejected.

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
- **HTTPS notebook pages**: if your notebook is served over HTTPS, the browser can block `http://127.0.0.1:<port>` iframes (mixed content). See {doc}`a_landing_pages/03_compatibility_matrix_must_be_explicit`.
- **Blocked viewer source**: startup stops if the configured source inventory or
  any declared asset cannot be fetched and verified. Use `--no-web-ui` only
  when you intentionally want the scientific data endpoints without a viewer.

## Troubleshooting (installation)

### Symptom: `pip install cellucid` fails with a Python version error

**Likely causes**
- You’re outside the supported Python 3.11–3.14 range.

**How to confirm**
```bash
python --version
```

**Fix**
- Create a new environment with Python 3.11–3.14 (see “Fast path”).

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

### Symptom: viewer startup reports that the web generation could not be established

**Likely causes**
- The configured web source is unavailable or blocked.
- A response has the wrong MIME type, byte length, hash, or build identity.
- Cache directory is not writable.

**How to confirm**
- Run:
  ```bash
  python -c "from cellucid import get_web_cache_dir; print(get_web_cache_dir())"
  ```

**Fix**
- Confirm the exact source contract directly:
  ```bash
  python -c "from cellucid.web_cache import ensure_web_ui_cached; print(ensure_web_ui_cached())"
  ```
- If the source is blocked, allow access to the configured
  `--web-source-url`; Cellucid does not substitute an older cached build.
- If the destination is not writable, retry with
  `--web-cache-dir /path/to/writable/directory`.

---

### Symptom: the embedded viewer is blank in JupyterLab / JupyterHub

**Likely causes**
- Your notebook is served from HTTPS or a remote origin, so the browser blocks an HTTP loopback iframe.

**Fix**
- Configure a browser-reachable HTTPS route (or an SSH tunnel) for the selected
  Cellucid port and pass its exact base as `client_server_url=`.
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

- Environment gotchas first: {doc}`a_landing_pages/03_compatibility_matrix_must_be_explicit`
- Then get something running end-to-end: {doc}`a_landing_pages/04_quick_start_3_levels`
