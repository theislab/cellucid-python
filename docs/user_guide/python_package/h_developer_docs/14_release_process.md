# Release process

This page describes how to ship a new `cellucid` Python package release to PyPI and keep ReadTheDocs in sync.

It assumes you are a maintainer with access to:
- push tags to the repo,
- the `PYPI_API_TOKEN` secret in GitHub Actions,
- and the ReadTheDocs token if you use the API trigger.

---

## Versioning policy (current)

`cellucid` is currently in alpha (`0.0.9`). <!-- CELLUCID_VERSION -->

Practical guidance:
- keep versions monotonic and tagged as `v<version>` (e.g. `v0.0.9`) <!-- CELLUCID_VERSION -->
- document user-facing changes in `CHANGELOG.md`
- avoid breaking export format compatibility without a coordinated web app change

---

## Release checklist (step-by-step)

### Step 1 — Update version

Update:
- `cellucid-python/pyproject.toml` → `[project].version`

Confirm locally:

```bash
python -c "import cellucid; print(cellucid.__version__)"
```

### Step 2 — Update changelog

Update:
- `cellucid-python/CHANGELOG.md`

Include:
- Added/Changed/Fixed/Removed sections as appropriate
- a date
- links at the bottom (`[Unreleased]`, `[<version>]`)

### Step 3 — Run quality gates locally

Recommended:

```bash
ruff format .
ruff check .
mypy src/cellucid
pytest
make -C docs clean html
```

### Step 4 — Tag the release

Create a git tag with a `v` prefix (required for CI trigger):

```bash
git tag v0.0.9  # CELLUCID_VERSION
git push origin v0.0.9  # CELLUCID_VERSION
```

### Step 5 — Watch CI publish to PyPI

The workflow:
- `cellucid-python/.github/workflows/pypi-publish.yml`

Trigger:
- push a tag matching `v*`

What it does:
1) checks out code
2) installs `build`
3) runs `python -m build` (wheel + sdist)
4) publishes via `pypa/gh-action-pypi-publish`

### Step 6 — Ensure docs are built

Docs are built in two ways:

- PR/branch validation:
  - `docs-check.yml` builds docs on PRs and pushes to main
- ReadTheDocs:
  - `.readthedocs.yaml` installs the `docs` extra and builds with `docs/conf.py`
  - `readthedocs.yml` can trigger an RTD build via API (if configured)

---

## Release pitfalls (common)

### Tag/version mismatch

If you tag `v0.0.9` but `pyproject.toml` still says a different version: <!-- CELLUCID_VERSION -->
- your published package will have the old version,
- and PyPI publishing may fail if that version already exists.

Always update `pyproject.toml` first.

### Broken docs links

Docs-check will fail if:
- `{doc}` targets don’t exist,
- MyST fences are malformed.

Run:

```bash
make -C docs clean html
```

### Network-dependent docs behavior

The docs build should not require network access.
Avoid adding steps that fetch remote resources at build time.

---

## Troubleshooting CI publishing

### Symptom: “PyPI publish failed: invalid token”

Confirm:
- `PYPI_API_TOKEN` secret exists and is valid,
- the workflow is running on the correct repository (forks won’t have secrets).

### Symptom: “Build failed”

Reproduce locally:

```bash
python -m pip install --upgrade build
python -m build
```

Then fix packaging errors before re-tagging (you may need a new version).
