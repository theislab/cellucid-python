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
(`test_viewer_channel_contract.py`) and Jupyter messaging, session files,
the history and the apply (`test_sessions.py`, against the golden history in
`tests/fixtures/history_golden/` that the web app and cellucid-r read too;
`test_session_layout.py`, the chunk layout the web app's writer makes: parts
shared, cyclic, expanding past what a page holds, unreferenced or renamed,
steps and values in many pieces, chunks stored as they are, parameters naming
lazy arrays;
`test_history_words.py`, the methods text and pseudocode rules the golden
branches do not reach, string for string with cellucid-r;
`test_history_model.py`, the generated model against
`history_model.schema.json`), the local MCP server (`test_mcp.py`: the
shared `mcp_tools.json` and `local_agents.json`, the agent
messages of the viewer channel, a viewer's publication and loopback endpoint
and its refusals, and `cellucid mcp` over stdio against a fake page:
initialize, lists, calls, errors, cancellation, elicitation, a list that never changes, a
viewer stopping), vector fields, cache and path confinement, documentation, cross-platform repository
paths, and release artifacts. The MuData tests need `mudata`, which the `dev`
extra brings with it (they are skipped without it). On Python 3.11 pip
installs mudata 0.3 (mudata 0.4 needs Python 3.12); `tests/conftest.py` sets its
`pull_on_update=False`, so every MuData the suite builds or reads keeps each
modality's columns in the modality, as mudata 0.4 does.
`test_minimum_python_grammar.py` parses every Python file with the grammar of
the oldest Python `requires-python` admits, including the f-strings only Python
3.12 accepts. Changes to a public contract require a focused regression test at
the owning boundary.

Some tests compare with the web app or run its readers. Every file this
package shares with the web app is compared byte for byte with the web
repository's file, its home (`assert_shared_file()` in `tests/_web_root.py`):
the action catalog, `mcp_tools.json`, the history model schema, the input
formats, the viewer protocol, the data-mapping, subsampling and latent-outlier
conformance fixtures, the category colors, the golden
history, the round-trip fixtures and the older-layout AnnData files. Others
read the session README's and codecs' constants (`test_sessions.py`,
`test_session_layout.py`) and run the web app's AnnData and priority readers
on what `prepare()` writes (`test_prepare_data_and_mudata_contract.py`,
`test_prepare_priority_order_contract.py`). Others compare with cellucid-r's
sources: every public name's R twin, `prepare()`'s arguments, the launch
defaults and the offline texts (`test_r_python_parity_contract.py`), the R
tutorials' calls (`test_round_trip_tutorials_contract.py`), the MuData
cross-writer stores (`test_prepare_data_and_mudata_contract.py`), and
cellucid-r's copy of each file only the two packages share, whose home is
here: the notebook widget's front end, `local_agents.json`, the
transition-drift fixture and the MuData cross-writer harness. No test holds a
digest of a shared file: a difference names the one command that makes every
copy agree, `node scripts/sync-shared-files.mjs` in the web repository
({ref}`the shared files <shared-files>`).

They find the web app at `CELLUCID_WEB_ROOT` (a checkout or worktree of the
`cellucid` repository) and cellucid-r at `CELLUCID_R_ROOT`, else in the
checkouts beside this repository (`tests/_web_root.py`). A checkout that is
there must hold the file a test reads: one it lacks fails the test, and so does
a named folder that is not a checkout of that repository. Only where a sibling
is not there at all (this repository alone, an sdist, an installed package, a
variable naming a folder that does not exist) are its tests skipped, saying
where they looked, and only when
`CELLUCID_REQUIRE_SIBLINGS` is not `1`: this repository's CI checks out
cellucid-r and the web app's code (each at `main`, after the package is
installed, so the day's pushes of the three repositories have landed) beside
this one, names them with the two variables, installs Node, and sets it, so
there a missing checkout or Node fails.

```bash
CELLUCID_WEB_ROOT=/path/to/cellucid CELLUCID_R_ROOT=/path/to/cellucid-r python -m pytest
```

---

## Running tests locally

From the `cellucid-python/` folder:

```bash
python -m pytest
```

After a change, run the core tests first: the tests marked `core` cover each
part of the package once (334 tests, about 21 s). The full run remains the
check before a change is done.

```bash
python -m pytest -m core
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
- `test.yml`: installs `.[dev]` (which brings `mudata`), checks out cellucid-r and
  the web app's code beside this repository, and runs the complete pytest suite
  (the cross-repository tests required, `CELLUCID_REQUIRE_SIBLINGS=1`), Ruff,
  mypy, and distribution checks on Python 3.11–3.14 on Ubuntu, and on 3.11,
  3.12 and 3.14 on macOS and on Windows
- `pypi-publish.yml`: runs when a GitHub Release is published (and on
  nothing else); proves the release and tag provenance, validates source and
  docs, normalizes the source distribution, installs both
  artifacts on all three operating systems, and publishes through trusted
  identity ({doc}`14_release_process`)
- `readthedocs.yml`: triggers the Read the Docs `latest` build via API on a
  push to `main`

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

## Decisions and why

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | CI checks out every sibling repository the tests compare with, at `main`, names it and requires it (`test_ci_checks_out_every_sibling_the_tests_compare_with`) | A test reading cellucid-r without its checkout failed every job of a run | A comparison that skips in CI |
| by 2026-11-12 | One stand-in server for every test that opens a viewer without serving (`tests/_stand_in_server.py`), checked against what the viewers read | The server gained `declared_dataset_id`, and three of four hand-written stand-ins lacked it | A stand-in written by hand in each test |
