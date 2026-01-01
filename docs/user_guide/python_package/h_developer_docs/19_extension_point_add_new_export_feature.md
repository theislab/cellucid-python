# Extension point: add a new export feature

This page is a **step-by-step guide** for adding a new capability to the Cellucid export format.

Examples of “export features”:
- a new per-cell binary field type
- a new manifest section
- a new kind of overlay (beyond the current vectors/connectivity)
- new metadata that the UI should display

This work almost always requires coordination with the web app.

Prerequisites:
- read {doc}`07_prepare_export_pipeline_architecture`
- read {doc}`08_export_format_spec_and_invariants`
- understand server behavior: {doc}`09_server_mode_architecture_endpoints_and_security`

---

## Step 1 — Define the user-visible behavior (before touching code)

Write down:

- What does the user get in the UI?
- What is the minimal example dataset where it works?
- What happens when the feature is missing (degrades gracefully vs hard error)?
- What are the scale limits (n_cells, n_features, file size)?

If you can’t describe the behavior, you can’t design the format.

---

## Step 2 — Decide where the data lives in the export folder

Choose:

- a new top-level file (e.g., `new_feature_manifest.json`), or
- a new directory (e.g., `new_feature/`), or
- an additive extension to an existing manifest (`obs_manifest.json`, `dataset_identity.json`, etc.).

Recommendations:
- Keep binaries headerless and describable by manifest context (consistent with existing format).
- Use a manifest to describe:
  - dtype
  - shape or how to infer shape
  - compression semantics (`.gz` raw gzip vs HTTP gzip)
  - how to map keys to filenames

---

## Step 3 — Implement the export writer (`prepare_data.py`)

Implementation file:
- `cellucid-python/src/cellucid/prepare_data.py`

Checklist:

1) Add input parameter(s) to `prepare(...)` (if user-facing).
2) Validate shapes early (fail fast with clear messages).
3) Write binaries using `_write_binary(...)` so compression behavior matches the rest of the format.
4) Update the relevant manifest(s) to advertise the feature.
5) Update `dataset_identity.json` to include high-level metadata if it improves UX.

Avoid:
- writing files without manifest references (they become “dead data”)
- silently overwriting existing files unless `force=True`

---

## Step 4 — Implement AnnData server parity (if applicable)

If the feature should work when serving AnnData directly:

1) Update `cellucid-python/src/cellucid/anndata_adapter.py` to generate the same bytes.
2) Update `cellucid-python/src/cellucid/anndata_server.py` to serve the route(s).

Goal:
- the browser viewer should not need to know whether the data came from disk or AnnData.

If parity is not required, document it explicitly (and ideally warn in UI or docs).

---

## Step 5 — Update the web app loader (required for new contracts)

In the web app repo (`cellucid/`), you typically must:

1) update the dataset loader to parse the new manifest fields,
2) fetch the new file(s),
3) integrate data into UI/rendering/state.

Start with web app developer docs:
- {doc}`../../web_app/p_developer_docs/index`

---

## Step 6 — Update documentation

Minimum doc updates:

- update the spec: {doc}`08_export_format_spec_and_invariants`
- update user guide pages that mention exports:
  - {doc}`../c_data_preparation_api/index`
  - {doc}`../d_viewing_apis/index`
- add troubleshooting entries for common failure modes

If the feature is UI-visible:
- add screenshot placeholders per the screenshots guide.

---

## Step 7 — Add tests (even if minimal)

Suggested tests:

1) Export writer test:
   - create a small synthetic dataset
   - run `prepare` into a temp folder
   - assert new files exist and manifests reference them
2) AnnData server test (if parity):
   - instantiate adapter on synthetic AnnData
   - call the new `get_*_binary(...)` method
   - assert bytes length and dtype

If you touch session bundles or codecs:
- add tests in `tests/test_sessions.py` using synthetic bundles (no real datasets).

---

## Step 8 — Compatibility review

Before merging:

- Does the change break older exports?
- Does it break older viewers?
- Is there a clear error message when a mismatch occurs?
- Is the change additive and optional where possible?

If you must break compatibility:
- bump a version marker (`dataset_identity.json.version` or a manifest `_format` field),
- and coordinate rollout with the web app.

---

## Troubleshooting

### Symptom: “I wrote new files but the UI ignores them”

Cause:
- the UI only knows about what manifests advertise.

Fix:
- add the feature to a manifest schema,
- update the web app loader to consume it.

### Symptom: “Works in export mode, not in AnnData mode”

Cause:
- missing parity implementation in `AnnDataAdapter` or `AnnDataRequestHandler`.

Fix:
- either implement parity or document that the feature is export-only.
