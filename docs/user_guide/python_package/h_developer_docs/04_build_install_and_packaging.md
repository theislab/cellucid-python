# Build, install, and packaging

This page documents how `cellucid-python` is packaged and shipped: how dependencies are defined, how versioning works, and how to build/install distributions.

If you are cutting a release, also read: {doc}`14_release_process`.

---

## Packaging model (what’s in `pyproject.toml`)

`cellucid-python` uses:
- **setuptools** as the build backend (`setuptools.build_meta`)
- a `src/` layout (`src/cellucid/`)
- extras for developer tooling and docs

Key sections:

- `[project]`
  - `name = "cellucid"`
  - `version = "0.0.9"` (alpha) <!-- CELLUCID_VERSION -->
  - runtime dependencies (`numpy`, `pandas`, `scipy`, `anndata`, `ipython`, `jupyter-server-proxy`)
- `[project.optional-dependencies]`
  - `dev`: `pytest`, `ruff`, `mypy`, `pre-commit`
  - `docs`: `sphinx`, `myst-nb`, theme, etc.
- `[project.scripts]`
  - `cellucid = "cellucid.cli:main"`

---

## Versioning and `cellucid.__version__`

`cellucid.__version__` is loaded from installed package metadata:
- `importlib.metadata.version("cellucid")`

This means:
- in a normal install, `__version__` matches the package version,
- in some editable/dev situations, metadata might be missing and `__version__` falls back to `"0.0.0"` (dev guard).

If you are debugging a version mismatch, confirm:

```bash
python -c "import cellucid; print(cellucid.__version__); import importlib.metadata as m; print(m.version('cellucid'))"
python -m pip show cellucid
```

---

## Installing (recommended patterns)

### Editable install for development

```bash
python -m pip install -e ".[dev,docs]"
```

### “User-like” install for testing packaging

In a clean environment:

```bash
python -m pip install .
```

---

## Building distributions (wheel + sdist)

The CI release workflow uses `python -m build`.

Locally:

```bash
python -m pip install --upgrade build
python -m build
```

Outputs appear under:
- `cellucid-python/dist/`

To test the wheel:

```bash
python -m pip install dist/*.whl
python -c "import cellucid; print(cellucid.__version__)"
cellucid --help
```

---

## Dependency philosophy (why “extras” exist)

Cellucid has a split personality:

- As a library/CLI, it should be lightweight to install and import.
- As a docs site, it needs Sphinx + MyST + theme + notebook tooling.

So:
- runtime dependencies live in `[project.dependencies]`,
- docs tooling lives in `[project.optional-dependencies].docs`,
- dev tooling lives in `[project.optional-dependencies].dev`.

This keeps `pip install cellucid` usable for typical users while allowing maintainers to build docs and run QA.

---

## Troubleshooting

### Symptom: “`cellucid.__version__` is `0.0.0`”

Likely cause:
- installed from source without metadata (some editable workflows) or a broken install.

Fix:
1) reinstall: `python -m pip install -e .`
2) confirm only one `cellucid` is on `sys.path` (`python -c "import cellucid; print(cellucid.__file__)"`)

### Symptom: “`cellucid` command runs the wrong code”

Likely cause:
- multiple environments/installs on the machine,
- old `cellucid` in PATH.

Confirm:

```bash
which cellucid
python -c "import cellucid; print(cellucid.__file__)"
```
