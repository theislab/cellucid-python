# Testing and CI

This page documents how to test `cellucid-python` locally and what CI currently enforces.

---

## What tests exist

Tests live under:
- `cellucid-python/tests/`

The suite covers input-format classification and serving every registry
format (`test_input_formats.py`, `test_serve_input_formats.py`), the
Cellucid-ready AnnData and MuData writers, exact scientific value and identity
contracts, HTTP routes and server lifecycle, the viewer channel
(`test_viewer_channel_contract.py`) and Jupyter messaging, sessions, vector
fields, cache and path confinement, documentation, cross-platform repository
paths, and release artifacts. The MuData tests need `mudata`, which the `dev`
extra brings with it (they are skipped without it). On Python 3.11 pip
installs mudata 0.3 (mudata 0.4 needs Python 3.12); `tests/conftest.py` sets its
`pull_on_update=False`, so every MuData the suite builds or reads keeps each
modality's columns in the modality, as mudata 0.4 does.
`test_minimum_python_grammar.py` parses every Python file with the grammar of
the oldest Python `requires-python` admits, including the f-strings only Python
3.12 accepts. Changes to a public contract require a focused regression test at
the owning boundary.

---

## Running tests locally

From the `cellucid-python/` folder:

```bash
python -m pytest
```

Useful variants:

```bash
python -m pytest -k sessions
python -m pytest -k vector_fields
```

Coverage (if installed):

```bash
python -m pytest --cov=cellucid
```

---

## Linting, formatting, and types

### Ruff

Format the Python files you intentionally changed, then lint the surfaces CI
lints (`src` and `tests`; the release workflow adds
`scripts/validate_release.py` and `scripts/normalize_sdist.py`):

```bash
python -m ruff format path/to/changed_file.py
python -m ruff check src tests
```

To check an intentionally formatted file without changing it:

```bash
python -m ruff format --check path/to/changed_file.py
```

### mypy

```bash
python -m mypy src/cellucid
```

Note: mypy is configured pragmatically (`ignore_missing_imports = true`).

---

## CI workflows (GitHub Actions)

Workflows live under:
- `cellucid-python/.github/workflows/`

Current workflows:

- `docs-check.yml`: builds Sphinx docs with warnings as errors
  (`sphinx-build -W --keep-going`) on PRs and pushes to `main`
- `test.yml`: installs `.[dev]` (which brings `mudata`) and runs the complete pytest suite,
  Ruff, mypy, and distribution checks on Python 3.11–3.14 on Ubuntu, and on
  3.11, 3.12 and 3.14 on macOS and on Windows
- `pypi-publish.yml`: proves tag provenance, validates source and docs,
  normalizes and verifies the source-distribution digest, installs both
  artifacts on all three operating systems, and publishes through trusted
  identity on tags `v*`
- `readthedocs.yml`: triggers a ReadTheDocs build via API

---

## Writing new tests (recommended patterns)

General principles:

- Prefer small, synthetic inputs over real datasets (speed + privacy).
- Test behavior and error messages for edge cases (shape mismatch, NaN handling, dtype selection).
- Avoid tests that require network access (hosted UI proxy) or browsers.

Core test categories:

1) The Cellucid-ready AnnData writer (`tests/test_prepare_anndata_contract.py`,
   and `tests/test_prepare_data_and_mudata_contract.py` for `data=` and MuData):
   - every output (`.h5ad`, Zarr v2/v3, `.zarr.zip`, `.h5mu`) reads back exactly with anndata or mudata
   - chunks, codecs, consolidated metadata, and `uns['cellucid']` match the specification; Zarr v3 metadata names only codecs of the v3 specification
   - every refusal happens before anything is written
2) Server routes:
   - health/info endpoints
   - the server answers store routes with exact `Range`/`206`/`HEAD`/`ETag` bytes
3) The viewer channel and notebook event plumbing:
   - the WebSocket upgrade, `Origin` and token checks, `hello`/`welcome`,
     the outbox and binary frames (`tests/_viewer_channel_client.py` is a
     minimal channel client)
   - `_handle_frontend_message` correctly updates `viewer.state`
   - hook registry triggers expected callbacks

---

## Troubleshooting

### Symptom: “Tests pass locally but docs-check fails”

Docs-check installs only the docs extra and builds Sphinx.

Common causes:
- broken MyST syntax (unclosed fences),
- bad `{doc}` links,
- missing files referenced by toctrees.

Fix: run the exact CI command, which treats warnings as errors:

```bash
sphinx-build -W --keep-going -b html docs docs/_build/html
```

### Symptom: “Ruff complains but formatting looks fine”

Run both checks on the files you changed:

```bash
python -m ruff format path/to/changed_file.py
python -m ruff check path/to/changed_file.py
```

Any rule suppression must be narrow and justified in the review description.
