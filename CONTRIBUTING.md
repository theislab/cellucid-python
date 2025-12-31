# Contributing to Cellucid (Python)

Contributions are welcome — code, docs, bug reports, design feedback, and reproducible examples all matter.

This file focuses on `cellucid-python` (the Python package + CLI), but also helps you route issues/PRs to the correct repo in the Cellucid ecosystem.

---

## Which repo should I contribute to?

Cellucid is split by responsibility:

| Repo | What it is | Contribute here when you… |
|---|---|---|
| `cellucid` | Web app (UI + state + WebGL rendering) | are fixing UI bugs, rendering/performance, figure export, sessions, or community annotation frontend |
| `cellucid-python` (this repo) | Python package + CLI (`prepare`, `serve`, `show_anndata`, hooks) + Sphinx docs | are fixing Python/CLI bugs, data prep/export, server endpoints, Jupyter hooks, or docs on ReadTheDocs |
| `cellucid-r` | R package for exporting data to the Cellucid viewer format | are changing the R exporter (`cellucid_prepare`) or adding R-side docs/tests |
| `cellucid-annotation` | GitHub repo template for community annotation | are changing the repo schema/validation/workflows |

If you’re not sure where a bug belongs, open an issue in the repo you’re currently using and include:
- how you loaded data (exports vs h5ad/zarr vs remote server vs Jupyter),
- the UI environment (hosted app vs local app vs Jupyter iframe),
- and the smallest reproduction you can share.

---

## Fast paths (pick your contribution type)

### I want to report a bug

Please include:
- Cellucid version (`python -c "import cellucid; print(cellucid.__version__)"` if available)
- Python version + OS
- how you loaded data (exports / `.h5ad` / `.zarr` / remote server / Jupyter)
- exact steps to reproduce (click-by-click or code)
- expected vs actual behavior
- logs:
  - Python server logs (if using `cellucid serve`)
  - browser console logs (if UI bug; enable debug with `localStorage.setItem('CELLUCID_DEBUG','true')` and reload)

### I want to contribute docs only

Docs live under `cellucid-python/docs/` and use Sphinx + MyST.

Fast workflow:
1) Edit the relevant `.md` (or notebook) under `docs/`
2) Build docs locally (see “Build docs” below)
3) Submit a PR

### I want to add/modify Python code

Fast workflow:
1) Create a dev environment (see “Development setup”)
2) Make a small, focused change
3) Add/adjust tests when behavior changes
4) Run `ruff`, `mypy`, and `pytest`
5) Submit a PR with a clear “what/why/how to verify”

---

## Development setup (cellucid-python)

### Prerequisites

- Python `>=3.10`
- Git

Recommended:
- a clean virtual environment (venv/conda)

### Clone

```bash
git clone https://github.com/theislab/cellucid-python.git
cd cellucid-python
```

### Install (editable + dev extras)

Using `venv`:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev,docs]"
```

This installs:
- runtime deps (numpy/pandas/scipy/anndata/etc.)
- dev tooling (pytest/ruff/mypy/pre-commit)
- docs tooling (sphinx/myst-nb/pydata-sphinx-theme/etc.)

---

## Code quality checks

### Formatting + lint (Ruff)

This repo uses Ruff for formatting and linting.

```bash
ruff format .
ruff check .
```

If you want Ruff to auto-fix what it can:

```bash
ruff check . --fix
```

### Type checking (mypy)

```bash
mypy src/cellucid
```

Notes:
- The mypy config is intentionally pragmatic (`ignore_missing_imports = true`).
- If you add new public APIs, prefer adding type hints even if the internal implementation stays flexible.

### Pre-commit (optional but recommended)

If you want checks to run automatically before you commit:

```bash
pre-commit install
pre-commit run --all-files
```

---

## Tests

Run the test suite:

```bash
pytest
```

Useful variants:

```bash
pytest -v --tb=short
pytest --cov=cellucid
```

Guidelines:
- Add tests when you change behavior (especially for edge cases and error messages).
- Prefer small synthetic inputs over large real datasets (privacy + speed).

---

## Build docs (Sphinx + MyST)

Docs live in `cellucid-python/docs/`.

Build HTML:

```bash
make -C docs html
```

Or explicitly:

```bash
sphinx-build -b html docs docs/_build/html
```

Open:
- `cellucid-python/docs/_build/html/index.html`

Doc-writing expectations (especially important for Cellucid):
- Write for mixed audiences: wet lab users, computational users, and developers.
- Prefer “layered” pages:
  - fast path (click-by-click / copy-paste)
  - practical path (parameters, performance)
  - deep path (design rationale, edge cases)
- Always include edge cases + troubleshooting for anything that can fail.

---

## PR guidelines (what makes a good PR here)

### Scope and clarity

- Keep PRs small and focused (one feature/bugfix at a time).
- Include a short “Why?” and “How to verify?” in the PR description.

### User-facing changes

If your change affects user-visible behavior:
- update docs and/or CLI help text
- add/adjust tests
- update `CHANGELOG.md` if it’s a user-facing change worth calling out

### Data privacy

Do not attach private datasets or patient-derived data to issues/PRs.
If you need a reproduction:
- reduce to a tiny synthetic dataset, or
- provide a script that generates a minimal failing example.

---

## Troubleshooting (common contributor problems)

### `pip install -e ".[dev,docs]"` fails

Common causes:
- old `pip` (upgrade first)
- environment mixing (conda + system Python)

Fix:
- create a fresh environment and retry
- run `python -m pip install --upgrade pip` before installing

### `ruff`/`mypy`/`pytest` aren’t found

Cause:
- dev extras not installed (or you forgot to activate your env).

Fix:
- activate your env and reinstall dev extras:
  - `python -m pip install -e ".[dev]"`

### Docs build fails on notebooks

Common causes:
- missing optional doc dependencies
- notebook execution errors (if notebooks are executed during build)

Fix:
- install docs extras: `python -m pip install -e ".[docs]"`
- run the failing page’s code cells locally and fix the root error

---

## Thanks

Thank you for improving Cellucid — contributions here directly improve how quickly people can explore and trust their single-cell results.
