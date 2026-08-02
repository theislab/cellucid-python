# Extension points overview

This page is a map of **safe ways to extend Cellucid Python** and a decision guide for when changes require coordination with the web app.

If you are planning a change, start here before coding:

- export format changes → {doc}`19_extension_point_add_new_export_feature`
- hooks/protocol changes → {doc}`20_extension_point_add_new_hook_event_or_command`

---

## The first decision: “Python-only” vs “contract change”

### Python-only changes (low coordination)

Examples:
- add a helper function (e.g., vector field utilities)
- add a CLI subcommand that doesn’t change the export format
- add docs pages or notebooks
- improve error messages, validation, performance inside existing contracts

Usually touches:
- `src/cellucid/*.py`
- `docs/`
- `tests/`

### Contract changes (requires coordination)

Examples:
- new files in the export folder
- changes to manifest schemas
- new hook event types or new frontend command types
- changes to session bundle chunk codecs

These must be coordinated with:
- the web app loader and UI behavior (`cellucid/` repo)

Start with:
- {doc}`08_export_format_spec_and_invariants`
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

### C) Add a new export feature (new file or new manifest fields)

This is the most coordination-heavy extension point.

Checklist:
- follow {doc}`19_extension_point_add_new_export_feature`
- update both:
  - `prepare_data/` (export writer)
  - `anndata_adapter.py` + `anndata_server.py` (dynamic mode parity), if applicable
- update the web app loader to consume it
- update the export spec docs

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

### Symptom: “I added a feature but it works only in AnnData mode”

Cause:
- you updated `AnnDataAdapter` but not `prepare`, or vice versa.

Fix:
- decide if the feature should exist in both modes,
- implement in both code paths, and update docs accordingly.

### Symptom: “I added a new file but the viewer never requests it”

Cause:
- the web app doesn’t know it exists (manifest/schema missing),
- or it is not hooked into the loader path.

Fix:
- update the export manifest/identity to advertise it,
- update the web app loader logic.
