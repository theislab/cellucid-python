# Contributing to Cellucid (Python)

Contributions are welcome — code, docs, bug reports, design feedback, and reproducible examples all matter.

This file focuses on `cellucid-python` (the Python package + CLI), but also helps you route issues/PRs to the correct repo in the Cellucid ecosystem.

By participating, you agree to follow the project’s Code of Conduct:
- `CODE_OF_CONDUCT.md`

If you’re reporting a security issue, please follow:
- `SECURITY.md`

---

## Which repo should I contribute to?

Cellucid is split by responsibility:

| Repo | What it is | Contribute here when you… |
|---|---|---|
| `cellucid` | Web app (UI + state + WebGL rendering) | are fixing UI bugs, rendering/performance, figure export, sessions, or community annotation frontend |
| `cellucid-python` (this repo) | Python package + CLI (`prepare`, `serve`, `show_anndata`, hooks) + Sphinx docs | are fixing Python/CLI bugs, `prepare()` (Cellucid-ready AnnData), server endpoints, Jupyter hooks, or docs on ReadTheDocs |
| `cellucid-r` | R package writing the same Cellucid-ready AnnData from R | are changing the R writer (`cellucid_prepare`) or adding R-side docs/tests |
| `cellucid-annotation` | GitHub repo template for community annotation | are changing the repo schema/validation/workflows |
| `cellucid-datasets` | The published demo catalog | are correcting a published generation or its catalog entry |
| `cellucid-demo-custom-datasets` | Worked example of publishing your own datasets | are changing that guide or its synthetic examples |

If you’re not sure where a bug belongs, open an issue in the repo you’re currently using and include:
- how you loaded data (`.h5ad`/`.zarr`, prepared with `prepare()` or not, vs remote server vs Jupyter),
- the UI environment (hosted app vs local app vs Jupyter iframe),
- and the smallest reproduction you can share.

---

## Fast paths (pick your contribution type)

### I want to report a bug

Please include:
- Cellucid version (`python -c "import cellucid; print(cellucid.__version__)"` if available)
- Python version + OS
- how you loaded data (`.h5ad` / `.zarr`, prepared or not / remote server / Jupyter)
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

- Python 3.11 through 3.14
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

## Comments say what the code does; history goes to the decision records

A comment or a docstring says what the code does and why it must be so: the
invariant it keeps, the failure it prevents, the unit of a number, the module
that owns the other half. It is written in the present tense and stays true
for as long as the code below it is unchanged. A docstring is held to this
first: Sphinx publishes it on the API pages.

A comment or a docstring never holds:

1. a CI run number or address;
2. a finding, audit or plan id (`F123`, `L25-3`, `R9W-5`, "round 9", "bug
   hunt", a section of a plan), or a path into a repository the reader cannot
   open;
3. a date of our own work (a date of an outside version the code depends on,
   such as a protocol revision or a Python release, stays);
4. who decided, or their words;
5. what the code used to do, what was tried, what was found or fixed;
6. numbers measured on a machine, with the machine, engine or dataset.

Each of these goes, in the same change, to the page that records decisions
for that part of the package (all under `docs/`, the documentation site):

| Code | Where its history goes |
|---|---|
| `prepare()` and the writer (`src/cellucid/prepare_data/`): codecs, shuffles, chunk sizes, sharding | The specification, `docs/user_guide/python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification.md` |
| The servers, the served store and the notebook widget (`anndata_server.py`, `_served_store.py`, `jupyter/`) | `docs/user_guide/web_app/p_developer_docs/26_decisions_served_store_and_widgets.md` |
| A test that pins a behavior of the web app | The decision record of that subsystem; `docs/user_guide/web_app/p_developer_docs/index.md` lists them |

In a decision record, a decision goes to "Decisions and why", numbers to
"Measurements that justify it", an attempt to "Tried and reverted". The
record tells the story in plain words, with its date and its numbers; it
holds no finding, audit or plan id either. Where a constant or a design rests
on one, the comment keeps the reason in one clause and names the page:
`(measured: 26_decisions_served_store_and_widgets.md)`.

What stays in a comment: the failure the code prevents, told as a condition
("on Windows a file held open cannot have another moved over it, so the file
is closed before the answer is sent"); pointers to tests, modules and pages
that exist; facts about outside versions.

Tests follow the same rule. A test's name says what it pins and holds no id.
A test module's docstring says what the file pins and what breaks without it;
the incident that led to it goes to the record.

`tests/test_comment_history_contract.py` checks the comments, docstrings and test names of `src`, `tests` and `scripts`, and the comments of `.github/workflows`, for finding ids, rounds, CI run addresses and plan pointers.

**Before and after: a reason told through an incident**
(`src/cellucid/anndata_server.py`, `_send_file`).

```python
# Before
"""...
On Windows a file held open cannot have another moved over it (a person's
``write_h5ad`` to a temporary name and a rename; ``open_shared`` lets it be
removed, not replaced), and a handle still held after the answer arrived
failed that save with "Access is denied" (Python package run on
windows-latest, 2026-11-07).
"""

# After
"""...
On Windows a file held open cannot have another moved over it (a person's
``write_h5ad`` to a temporary name and a rename; ``open_shared`` lets it be
removed, not replaced), so a handle held after the answer would fail that
save.
"""
```

Record 26, "Decisions and why", gains the row: 2026-11-07; the file is held
open only while its body is sent; why: on Windows a handle still held after
the answer arrived failed a person's save with "Access is denied" (the
package's CI); rules out: holding the file for the length of the request.

**Before and after: a measurement that justifies a design**
(`src/cellucid/prepare_data/_store.py`, the module docstring).

```python
# Before
"""...
Measured on the 561,947-cell, 8,192-gene Suo atlas with the web app's own
reader (``scripts/benchmark-anndata.mjs`` in the web repository; the numbers
are in the specification page): shuffling cut the Zarr store from 243.4 MB to
198 MB -- the gene indices by a third, embeddings and graph weights by 7-9% --
while shuffled expression values of real data grew 14-45% ...
"""

# After
"""...
The codec, the shuffle and the chunk sizes are measured choices, taken with
the web app's own reader (``scripts/benchmark-anndata.mjs`` in the web
repository; measured: 09_cellucid_ready_anndata_specification.md).
"""
```

The specification already holds those numbers with their dataset and date; a
number it lacks is added there in the same change.

---

## Tests

Run the test suite:

```bash
pytest
```

After every change, run the core tests first: each part of the package once,
in about half a minute on an idle machine. The tests marked `core` say what they protect.

```bash
python -m pytest -m core
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
sphinx-build -W --keep-going -b html docs docs/_build/html
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
