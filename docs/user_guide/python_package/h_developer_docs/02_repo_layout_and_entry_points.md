# Repo layout and entry points

This page is a **practical map of the `cellucid-python` repository**: where the important code lives, which files are “entry points”, and how to navigate changes without getting lost.

:::{admonition} Audience
:class: note

- Beginners: use the “Where do I change X?” table.
- Experienced devs: use the module map + entry points + invariants.
:::

---

## Repo layout (top level)

The important top-level directories are:

- `src/`: the Python package (`src/cellucid/`).
- `docs/`: the Sphinx + MyST documentation site (this is what ReadTheDocs builds).
- `tests/`: pytest tests.
- `.github/workflows/`: CI (docs build, PyPI publish, RTD trigger).
- `scripts/`: project-specific dev scripts (may be environment-specific).

If you like a “tree view”:

```text
cellucid-python/
  .github/workflows/
  docs/
    conf.py
    index.md
    user_guide/
  src/
    cellucid/
      __init__.py
      cli.py
      prepare_data.py
      server.py
      anndata_adapter.py
      anndata_server.py
      jupyter.py
      _server_base.py
      session_bundle.py
      session_codecs.py
      anndata_session.py
      vector_fields.py
      web_cache.py
  tests/
  pyproject.toml
  README.md
  CONTRIBUTING.md
  CHANGELOG.md
```

---

## The entry points (where execution “starts”)

### 1) CLI entry point

The installed command `cellucid` is defined in `cellucid-python/pyproject.toml`:

- `cellucid = "cellucid.cli:main"`

Implementation:
- `cellucid-python/src/cellucid/cli.py`

Right now the CLI is intentionally minimal and centers around:
- `cellucid serve <path>` (auto-detect `.h5ad`, `.zarr`, or export directory)

See: {doc}`05_cli_architecture_and_commands`

### 2) Public Python API entry point

The public import surface is defined in:
- `cellucid-python/src/cellucid/__init__.py`

Important design choice:
- it uses `__getattr__` to **lazy-import** heavy dependencies (numpy/pandas/scipy/anndata) so `cellucid --help` stays fast.

If you add a new public API:
1) implement it in a module under `src/cellucid/`,
2) expose it via `__getattr__`,
3) add it to `__all__`,
4) document it in the API docs or user guide.

### 3) Server “entry points”

There are two server implementations:

- Exported/static: `cellucid-python/src/cellucid/server.py`
  - public function: `cellucid.serve(...)`
- AnnData/dynamic: `cellucid-python/src/cellucid/anndata_server.py`
  - public function: `cellucid.serve_anndata(...)`

Both share CORS + web-asset proxy + event/session upload handlers from:
- `cellucid-python/src/cellucid/_server_base.py`

### 4) Jupyter entry point

Notebook embedding lives in:
- `cellucid-python/src/cellucid/jupyter.py`

Public functions/classes (via `cellucid.__init__`):
- `show(...)`, `show_anndata(...)`
- `CellucidViewer`, `AnnDataViewer`

See: {doc}`10_jupyter_embedding_architecture`

---

## “Where do I change X?” (common developer tasks)

| You want to change… | Start in | Notes |
|---|---|---|
| Export filenames/manifests | `src/cellucid/prepare_data.py` | Coordinate with web app; see {doc}`08_export_format_spec_and_invariants` |
| Quantization behavior | `src/cellucid/prepare_data.py` | Reserved NaN/Inf markers must match viewer expectations |
| AnnData lazy loading and caching | `src/cellucid/anndata_adapter.py` | Watch memory blowups (CSR→CSC), LRU cache |
| HTTP routes for AnnData | `src/cellucid/anndata_server.py` | Keep behavior consistent with export mode |
| CORS / origin rules | `src/cellucid/_server_base.py` | Security-sensitive |
| Hosted UI caching/offline | `src/cellucid/web_cache.py` | Be careful with cache invalidation |
| Notebook embedding issues | `src/cellucid/jupyter.py` | Remote notebook environments are tricky |
| Session bundle parsing | `src/cellucid/session_bundle.py` | Treat bundles as untrusted input |
| Applying sessions to AnnData | `src/cellucid/anndata_session.py` | Mismatch policy + column conflict policy matter |
| Vector field helpers | `src/cellucid/vector_fields.py` | Keep naming conventions stable (`*_umap_2d`, etc.) |

---

## Search tips (fast navigation for real work)

When you’re trying to find where something is implemented, search the repo before guessing.

Examples (run from the `cellucid-python/` folder):

```bash
rg "dataset_identity\\.json" src/cellucid
rg "CELLUCID_WEB_PROXY_CACHE_DIR" -S src/cellucid
rg "requestSessionBundle" src/cellucid/jupyter.py
rg "obs_manifest" src/cellucid/prepare_data.py
```

---

## Troubleshooting

### Symptom: “I can’t find where `cellucid.prepare` is defined”

`cellucid.prepare` is a lazy import from `cellucid.__init__`.

Look here:
- `cellucid-python/src/cellucid/__init__.py` (`__getattr__`)
- `cellucid-python/src/cellucid/prepare_data.py` (`def prepare(...)`)

### Symptom: “I changed a module but the CLI got slow”

Likely causes:
- importing heavy deps at module import time (numpy/pandas/scipy/anndata),
- importing `cellucid` at CLI import time in a way that triggers heavy imports.

Rule: keep `cellucid.cli` lightweight; push heavy imports into the code path that needs them.
