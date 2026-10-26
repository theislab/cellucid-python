# Extension points overview

This page is a map of **safe ways to extend Cellucid Python** and a decision guide for when changes require coordination with the web app.

If you are planning a change, start here before coding:

- Cellucid-ready AnnData changes → {doc}`07_prepare_writer_architecture`
- hooks/protocol changes → {doc}`20_extension_point_add_new_hook_event_or_command`

---

## The first decision: “Python-only” vs “contract change”

### Python-only changes (low coordination)

Examples:
- add a helper function (e.g., vector field utilities)
- add a CLI subcommand that doesn’t change the Cellucid-ready AnnData contract
- add docs pages or notebooks
- improve error messages, validation, performance inside existing contracts

Usually touches:
- `src/cellucid/*.py`
- `docs/`
- `tests/`

### Contract changes (requires coordination)

Examples:
- new elements or `uns['cellucid']` keys in the Cellucid-ready AnnData
- new launch defaults
- new hook event types or new frontend command types
- changes to session bundle chunk codecs

These must be coordinated with:
- the web app loader and UI behavior (`cellucid/` repo)
- cellucid-r, which writes the same Cellucid-ready AnnData and speaks the same
  viewer channel (the shared `viewer-protocol.json`)

Start with:
- {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`
- {doc}`11_hooks_events_protocol_and_schema`
- web app developer docs: {doc}`../../web_app/p_developer_docs/index`

---

## Extension point catalog

### A) Add a new public Python API function/class

Checklist:
1) implement in `src/cellucid/<module>.py`
2) export via `src/cellucid/__init__.py` (`__getattr__` + `__all__`)
3) add docs (user guide or API reference)
4) add tests for edge cases

Watch out for:
- importing heavy dependencies at module import time (CLI startup)

### B) Add a new CLI command

Checklist:
1) add a subparser in `src/cellucid/cli.py`
2) keep parser lightweight (no heavy imports)
3) import heavy deps inside command handler
4) document in {doc}`05_cli_architecture_and_commands`

### C) Add something to the Cellucid-ready AnnData (a new element or `uns['cellucid']` key)

This is the most coordination-heavy extension point.

Checklist:
- read {doc}`07_prepare_writer_architecture`
- update together:
  - `prepare_data/` (the Python writer) and cellucid-r's `cellucid_prepare()`
  - the web app's AnnData reader (`cellucid/assets/js/data/anndata/`, and
    `anndata-cellucid-metadata.js` / `requireAnnDataLaunchDefaults` for a key)
  - `_launch_defaults.py` for a launch default, so the server accepts it too
  - {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`

### D) Add a new hook event or command

Checklist:
- follow {doc}`20_extension_point_add_new_hook_event_or_command`
- update:
  - frontend event/command handlers (web app)
  - Python viewer helpers (`jupyter/`) if adding a new public method
  - docs and examples

### E) Extend session bundle parsing

Session bundles are untrusted input and must remain robust.

Checklist:
- define one exact schema for every new chunk or codec
- keep hard caps on decompression and allocations
- update tests in `tests/test_sessions.py` (add synthetic bundles)

---

## Contract-change strategy

When adding new fields/files:
- update every writer, server adapter, loader, validator, test, and format document
  in the same change;
- reject unknown keys and unsupported schemas with a clear error;
- do not retain alternate schema readers or redirect paths.

---

## Troubleshooting

### Symptom: “I added a key to `uns['cellucid']` and the viewer refuses the file”

Cause:
- the web app refuses any key it does not know, so a misspelt default never
  passes unnoticed.

Fix:
- add the key to the web app's `anndata-cellucid-metadata.js` (or
  `requireAnnDataLaunchDefaults`) in the same change, and to cellucid-r.

### Symptom: “I added an element but the viewer never reads it”

Cause:
- the web app's AnnData reader is not hooked up to it.

Fix:
- update the web app's reader in the same change.
