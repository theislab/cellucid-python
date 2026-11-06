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
- `.github/workflows/`: CI (tests, docs build, PyPI publish, RTD trigger).
- `scripts/`: release and maintenance scripts (`validate_release.py`, `sweep_versions.py`, `normalize_sdist.py`, `publishing/`, conformance-fixture generators).
- `examples/`: example scripts.

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
      _contracts.py
      _console.py
      _launch_defaults.py
      _non_finite.py
      _plural.py
      data_mapping.py
      subsample.py
      prepare_data/
        __init__.py
        _prepare.py
        _data.py
        _inputs.py
        _gene_major.py
        _mudata.py
        _outliers.py
        _priority.py
        _store.py
      input_formats.json
      history_model.schema.json
      _anndata_files.py
      _scratch.py
      _source_arrays.py
      _served_store.py
      _served_formats.py
      anndata_server.py
      jupyter/
        __init__.py
        _wire.py
        _hooks.py
        _history.py
        _viewer.py
        _widget.py
        static/
          widget.js
          widget.css
      _server_base.py
      _viewer_channel.py
      _websocket.py
      session.py
      session_codecs.py
      history.py
      _history_model.py
      _history_rules.py
      anndata_session/
        __init__.py
        _schema.py
        _records.py
        _primitives.py
        _cells.py
        _labels.py
        _highlights.py
        _fields.py
        _mudata.py
        _apply.py
        _columns.py
        _results.py
      vector_fields.py
      web_cache.py
      py.typed
  examples/
  scripts/
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
- the `serve` subcommand (serves any file or folder Cellucid reads, by the input format registry `src/cellucid/input_formats.json`: AnnData, MuData, Loom, 10x, MatrixMarket, R objects)

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

There is one server implementation:

- the dataset's bytes (the browser reads them): `cellucid-python/src/cellucid/anndata_server.py`
  - public function: `cellucid.serve_anndata(...)`, class `cellucid.AnnDataServer`

It builds on `cellucid-python/src/cellucid/_server_base.py` (CORS, `Host`
checks, byte ranges, verified web-asset serving) and hands `/_cellucid/ws` to
the viewer channel in `cellucid-python/src/cellucid/_viewer_channel.py`.

### 4) Jupyter entry point

Notebook embedding lives in:
- `cellucid-python/src/cellucid/jupyter/`

Public functions/classes (via `cellucid.__init__`):
- `show_anndata(...)`
- `AnnDataViewer`

See: {doc}`10_jupyter_embedding_architecture`

---

## “Where do I change X?” (common developer tasks)

| You want to change… | Start in | Notes |
|---|---|---|
| What `prepare()` accepts and how it converts inputs | `src/cellucid/prepare_data/_inputs.py` | Keep the rules equal to cellucid-r's; see {doc}`07_prepare_writer_architecture` |
| Chunks, codecs, encodings of a prepared file | `src/cellucid/prepare_data/_store.py` | Change only with measurements; the web reader and cellucid-r must agree |
| Launch defaults (`uns['cellucid']['defaults']`, `/_cellucid/datasets`) | `src/cellucid/_launch_defaults.py` | Key names are the web app's `requireAnnDataLaunchDefaults` |
| Accepted input formats (by `input_formats.json`), undecodable codecs | `src/cellucid/_anndata_files.py` | Metadata only, never an array; the browser's decodable codec sets must match the web app's HDF5/Zarr reader |
| Priority order of a prepared file (`priority_seed`, the record, original positions, prefix counts, `restore_original_order()`, `stored_order()`) | `src/cellucid/prepare_data/_priority.py` | Rule for rule cellucid-r's `R/priority.R`; the web app's `data/anndata/priority-layout.js` verifies what it writes |
| Which cells a subsample keeps (`subsample_cells()`) | `src/cellucid/subsample.py` | The subsampling specification; the shared conformance fixture pins it in all three repositories |
| The store in priority order served for an object or a file not in priority order | `src/cellucid/_served_store.py` | Chunks made on request from the source (a `Range` makes only its bytes); the arrays `prepare()` writes, value for value, chunked alike but for an expression matrix's entries (2^18 values a chunk: a request on this machine, at the server and in the browser, costs more than 2^14 values move); uncompressed; expression matrices laid out gene-major by a background worker, a chunk asked for before made from the source at once (`503` + `Retry-After` only for a file's matrix not counted yet) |
| The gene-major layout of an expression matrix (`prepare()`, the backed path, the served store) | `src/cellucid/prepare_data/_gene_major.py` | One pass for every source; values as `prepare()` stores them; bytes pinned against the earlier writer |
| MuData, Loom and 10x files served in priority order by path | `src/cellucid/_served_formats.py` | Mirrors the web app's layout views (`views/loom.js`, `views/tenx-*.js`) |
| The notebook widget (`Widget`, `show_anndata()`) | `src/cellucid/jupyter/_widget.py`, `jupyter/static/widget.js` | The page's reads and the viewer channel over the widget's comm; `widget.js` is byte for byte cellucid-r's `inst/widget/widget.js` |
| HTTP routes for the served data | `src/cellucid/anndata_server.py` | Exact `Range`/`206`/`HEAD`/`ETag` answers; the web app's reader depends on them |
| The viewer channel (commands, events, session bundles) | `src/cellucid/_viewer_channel.py`, `src/cellucid/jupyter/_wire.py` | Cross-language rules in the web app's `data/viewer-protocol.js`; see {doc}`11_hooks_events_protocol_and_schema` |
| CORS / origin rules | `src/cellucid/_server_base.py` | Security-sensitive |
| Verified web-generation establishment | `src/cellucid/web_cache.py` | Preserve inventory, byte-verification, staging, and atomic-publication invariants |
| Notebook embedding issues | `src/cellucid/jupyter/_viewer.py` | Remote notebook environments are tricky |
| Reading session files | `src/cellucid/session.py` | Treat session files as untrusted input; mirror every check of the web app's `history-file.js` |
| Reading a history (trackers, branches, steps, notes, turns, diffs, methods, pseudocode) | `src/cellucid/history.py` over `_history_model.py` (generated by `scripts/generate_history_model.py`) and `_history_rules.py` | Every record the history model's; the web app's words; the same reading as cellucid-r, held to the golden history |
| Applying one step of a session to AnnData | `src/cellucid/anndata_session/_apply.py` (labels `_labels.py`, fields `_fields.py`, highlights `_highlights.py`, cells `_cells.py`) | Check only what the step needs, plan every row before writing, write all at once, never replace a column the person made; the same rows as cellucid-r (`history_golden.expected.json`) |
| The live history in a notebook (turns, `load_session`, `history_changed`) | `src/cellucid/jupyter/_history.py`, `_viewer.py` | One protocol for Python and R (`viewer_protocol.json`) |
| Pull, push, subset and gene sets (the round trip) | `src/cellucid/anndata_session/_results.py`, `_columns.py` (the one column writer sessions and pulls share) | Rule for rule cellucid-r's `R/round-trip.R` |
| Vector field helpers | `src/cellucid/vector_fields.py` | Write keys by the vector-field key rule (`<field>_<basis>` for `X_<basis>`, `cellucid.data_mapping.parse_vector_field_key`) |

---

## Search tips (fast navigation for real work)

When you’re trying to find where something is implemented, search the repo before guessing.

Examples (run from the `cellucid-python/` folder):

```bash
rg "uns\\[.cellucid.\\]" src/cellucid
rg "web_cache_dir" src/cellucid
rg "requestSessionBundle" src/cellucid/jupyter
rg "SLICED_CHUNK_ELEMENTS" src/cellucid/prepare_data
```

---

## Troubleshooting

### Symptom: “I can’t find where `cellucid.prepare` is defined”

`cellucid.prepare` is a lazy import from `cellucid.__init__`.

Look here:
- `cellucid-python/src/cellucid/__init__.py` (`__getattr__`)
- `cellucid-python/src/cellucid/prepare_data/_prepare.py` (`def prepare(...)`)

### Symptom: “I changed a module but the CLI got slow”

Likely causes:
- importing heavy deps at module import time (numpy/pandas/scipy/anndata),
- importing `cellucid` at CLI import time in a way that triggers heavy imports.

Rule: keep `cellucid.cli` lightweight; push heavy imports into the code path that needs them.
