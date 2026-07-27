# Testing and CI

This page documents how to test `cellucid-python` locally and what CI currently enforces.

---

## What tests exist

Tests live under:
- `cellucid-python/tests/`

The suite covers AnnData and Zarr loading, prepared exports, exact scientific
value and identity contracts, HTTP routes and server lifecycle, Jupyter
messaging, sessions, vector fields, cache and path confinement, documentation,
cross-platform repository paths, and release artifacts. Changes to a public
contract require a focused regression test at the owning boundary.

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

Format the Python files you intentionally changed, then lint the maintained
source, test, and release-script surfaces:

```bash
python -m ruff format path/to/changed_file.py
python -m ruff check src tests scripts
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

- `docs-check.yml`: builds Sphinx docs on PRs and pushes to `main`
- `test.yml`: runs the complete pytest suite, Ruff, mypy, and distribution
  checks across Python 3.11–3.14 on Ubuntu, representative Python versions on
  macOS, and representative Python versions on Windows
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

1) Export format:
   - `prepare(...)` produces required files
   - manifests reference existing files
   - quantization markers and min/max are correct
2) Server routes:
   - health/info endpoints
   - AnnData server route handlers return bytes with correct length/dtype
3) Notebook event plumbing:
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

Fix:

```bash
make -C docs clean html
```

### Symptom: “Ruff complains but formatting looks fine”

Run both checks on the files you changed:

```bash
python -m ruff format path/to/changed_file.py
python -m ruff check path/to/changed_file.py
```

Any rule suppression must be narrow and justified in the review description.
