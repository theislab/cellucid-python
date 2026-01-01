# Local development setup

This page walks you through **setting up a development environment** for `cellucid-python` (Python package + CLI + docs).

It is intentionally detailed so that:
- first-time contributors don’t get stuck on environment issues,
- computational users can reproduce bugs reliably,
- and maintainers have a repeatable workflow.

---

## Prerequisites

Required:
- Python `>=3.10`
- Git

Recommended:
- a fresh virtual environment (`venv`, `conda`, `uv`, etc.)

Optional (only if you work on notebook embedding):
- JupyterLab / VSCode notebooks / classic notebook

---

## Step 0 — Choose which repo(s) you’re working on

Cellucid is split across repositories:

- **Python package** (this page): `cellucid-python/`
- **Web app** (UI/rendering): `cellucid/` (see {doc}`../../web_app/p_developer_docs/index`)

If you are changing the export format or hooks protocol, you will likely touch both repos.

---

## Step 1 — Create an isolated environment

Using `venv`:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
python -m pip install --upgrade pip
```

If you use conda, the key requirement is “Python 3.10+ and a clean environment”; the rest is the same.

---

## Step 2 — Install in editable mode (dev + docs extras)

From the `cellucid-python/` folder:

```bash
python -m pip install -e ".[dev,docs]"
```

What this does:
- installs runtime dependencies,
- installs developer tooling (`pytest`, `ruff`, `mypy`, `pre-commit`),
- installs documentation tooling (Sphinx + MyST).

---

## Step 3 — Quick sanity checks

### 3.1 Import + version

```bash
python -c "import cellucid; print(cellucid.__version__)"
```

### 3.2 CLI help

```bash
cellucid --help
cellucid serve --help
```

If the CLI is slow or crashes, jump to: {doc}`05_cli_architecture_and_commands`.

### 3.3 Run tests

```bash
pytest
```

If tests are missing for what you changed, see: {doc}`13_testing_and_ci`.

---

## Step 4 — Build the docs locally

From the `cellucid-python/` folder:

```bash
make -C docs html
```

Then open:
- `cellucid-python/docs/_build/html/index.html`

Common variants:

```bash
make -C docs clean
make -C docs linkcheck
```

Doc-writing conventions live here:
{doc}`15_docs_development_and_style_guide`

---

## Step 5 — Run a viewer locally (for manual testing)

You can test with either:

### Option A: exported folder (fastest + most representative)

1) Export your dataset:

```python
from cellucid import prepare
# ... call prepare(...) to write an export folder ...
```

2) Serve it:

```bash
cellucid serve ./my_export
```

### Option B: AnnData server (convenient during iteration)

```bash
cellucid serve ./data.h5ad
```

Both modes should open the viewer in a browser automatically (unless `--no-browser`).

---

## Step 6 — Notebook embedding (optional, but common)

Notebook embedding is the “hardest environment” because of proxies, mixed-content, and remote kernels.

Minimal test:

```python
from cellucid import show_anndata
viewer = show_anndata("path/to/data.h5ad")
viewer  # show the iframe output
```

If the iframe is blank or hooks don’t work:
1) run `viewer.debug_connection()`,
2) follow {doc}`12_debugging_playbook`,
3) see notebook architecture details in {doc}`10_jupyter_embedding_architecture`.

---

## Offline / airgapped development notes

The Python servers run in **hosted-asset proxy** mode by default:
- they fetch the web UI from `https://www.cellucid.com`,
- cache it locally,
- and serve it from the same origin as the dataset server (avoids mixed-content).

If you are offline:
- run once while online to populate the cache,
- then reuse that cache.

You can control where the cache lives:
- `CELLUCID_WEB_PROXY_CACHE_DIR` (see {doc}`06_configuration_env_vars_and_logging`)

---

## Troubleshooting

### Symptom: “`pip install -e ".[dev,docs]"` fails”

Likely causes:
- conflicting packages already installed in the environment,
- old pip/setuptools,
- corporate proxy restrictions.

Fix checklist:
1) ensure you activated the correct environment,
2) upgrade pip: `python -m pip install --upgrade pip`,
3) try a fresh environment.

### Symptom: “Docs build fails with missing dependencies”

Confirm you installed the docs extras:

```bash
python -m pip show sphinx myst-nb
```

Then rebuild:

```bash
make -C docs clean html
```

### Symptom: “Browser opens but shows ‘viewer UI unavailable’”

This usually means the hosted UI could not be fetched and no cached copy exists.

Fix:
- ensure network access to `https://www.cellucid.com`, or
- pre-populate the cache once while online, or
- set `CELLUCID_WEB_PROXY_CACHE_DIR` to a writable, persistent path.
