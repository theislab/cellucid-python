# Testing and CI

This page documents how to test `cellucid-python` locally and what CI currently enforces.

---

## What tests exist today

Tests live under:
- `cellucid-python/tests/`

Current coverage focuses on:
- session bundle decoding + “apply session to AnnData”
- vector field utilities

If you change servers or export format behavior, you should add tests even if none exist yet for that area.

---

## Running tests locally

From the `cellucid-python/` folder:

```bash
pytest
```

Useful variants:

```bash
pytest -v --tb=short
pytest -k sessions
pytest -k vector_fields
```

Coverage (if installed):

```bash
pytest --cov=cellucid
```

---

## Linting, formatting, and types

### Ruff (format + lint)

```bash
ruff format .
ruff check .
```

Auto-fix what’s safe:

```bash
ruff check . --fix
```

### mypy

```bash
mypy src/cellucid
```

Note: mypy is configured pragmatically (`ignore_missing_imports = true`).

---

## CI workflows (GitHub Actions)

Workflows live under:
- `cellucid-python/.github/workflows/`

Current workflows:

- `docs-check.yml`: builds Sphinx docs on PRs and pushes to `main`
- `pypi-publish.yml`: builds and publishes distributions to PyPI on tags `v*`
- `readthedocs.yml`: triggers a ReadTheDocs build via API

Important current limitation:
- There is not (yet) a dedicated “run pytest in CI” workflow. If you add tests for critical behavior, consider adding CI coverage in a future change.

---

## Writing new tests (recommended patterns)

General principles:

- Prefer small, synthetic inputs over real datasets (speed + privacy).
- Test behavior and error messages for edge cases (shape mismatch, NaN handling, dtype selection).
- Avoid tests that require network access (hosted UI proxy) or browsers.

Suggested “test categories” for future expansion:

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

Run:

```bash
ruff format .
ruff check . --fix
```

If you need to suppress a lint rule, do it surgically and justify it in the PR description.
