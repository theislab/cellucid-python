# Session bundles: apply to AnnData (`apply_session_to_anndata`)

This page documents how to take UI state captured from the viewer (a `.cellucid-session` bundle) and apply parts of it back onto an `AnnData`.

The main workflow is:

1. Capture a session bundle (no download): {doc}`10_session_bundles_get_session_bundle`
2. Apply it to AnnData:
   - `viewer.apply_session_to_anndata(...)` (one-liner), or
   - `bundle.apply_to_anndata(...)` (explicit bundle object), or
   - `cellucid.apply_cellucid_session_to_anndata(...)` (function-level control)

## At a glance

**Audience**
- Wet lab / beginner: this is optional unless you’re collaborating or exporting curated labels back into analysis.
- Computational users: this is the bridge from “interactive UI edits” → “analysis-ready AnnData”.
- Developers: read mismatch policies and column naming rules.

## Quickstart (copy/paste)

### One-liner (recommended)

```python
adata2 = viewer.apply_session_to_anndata(adata, inplace=False)
```

### Explicit two-step

```python
bundle = viewer.get_session_bundle(timeout=60)
adata2 = bundle.apply_to_anndata(adata, inplace=False)
```

```{important}
Session application is currently **index-based**:

The session references cells by **row index** (0-based). Only apply a session to an `AnnData` whose row order matches the dataset that produced the session.
```

## What gets applied (current behavior)

As implemented in `cellucid-python/src/cellucid/anndata_session.py`, applying a session can add:

### 1) Highlight groups → boolean `adata.obs` columns

If the bundle contains highlight metadata:
- For each highlight group, a boolean mask column is created in `adata.obs`.
- Default prefix: `cellucid_highlight__`

Example output:

```text
adata.obs["cellucid_highlight__highlight_12"]  # True/False per cell
```

### 2) User-defined categorical fields → `adata.obs` / `adata.var` columns

If the bundle contains user-defined categorical fields:
- Fields sourced from `obs` become `adata.obs` categorical columns.
- Fields sourced from `var` become `adata.var` categorical columns.

### 3) Provenance metadata → `adata.uns["cellucid"]`

By default, application stores:
- the session manifest,
- dataset fingerprint,
- and “applied settings”
under:

```text
adata.uns["cellucid"]["session"]
```

## Dataset mismatch policies (avoid accidental corruption)

Session bundles include a dataset fingerprint with (at least):
- `cellCount`
- `varCount`
- `datasetId` (if available)

When applying:
- if counts/IDs mismatch, you choose what happens via `dataset_mismatch`:
  - `"error"`: raise and stop
  - `"warn_skip"` (default): warn and skip dataset-dependent chunks
  - `"skip"`: silently skip dataset-dependent chunks

In practice, if you subset/shuffle your AnnData after the session was created, you will almost always have a mismatch.

## Column naming and conflicts

Applying a session may add columns that collide with existing names.

Controls:
- `column_conflict="suffix"` (default): create `name__2`, `name__3`, ...
- `column_conflict="overwrite"`: replace existing columns
- `column_conflict="error"`: raise on conflict

Controls for prefixes:
- `highlights_prefix="cellucid_highlight__"`
- `user_defined_prefix=""`

## Full API (most important knobs)

The function-level API is:

```python
from cellucid import apply_cellucid_session_to_anndata

adata2, summary = apply_cellucid_session_to_anndata(
    bundle,
    adata,
    inplace=False,
    dataset_mismatch="warn_skip",
    expected_dataset_id=None,
    add_highlights=True,
    add_user_defined_fields=True,
    store_uns=True,
    column_conflict="suffix",
    return_summary=True,
)
```

This returns an `ApplySummary` with:
- which columns were added,
- whether chunks were skipped due to mismatch,
- mismatch reasons.

## Edge cases

- Applying to backed AnnData may be slow if it triggers copies; plan memory accordingly.
- If pandas is not installed, applying a session raises an ImportError (it needs pandas to write categorical/boolean columns).
- If a highlight group contains out-of-range indices, they are bounds-checked and safely ignored.

## Troubleshooting

### Symptom: “Dataset fingerprint mismatch”

This is almost always because:
- you filtered/subset/reordered `adata` after opening the viewer, or
- you’re applying the session to a different dataset.

Fix options:
- apply to the exact same AnnData (same row order),
- set `dataset_mismatch="error"` to catch mistakes early,
- or accept skipping dataset-dependent chunks if you only want provenance stored.

### Symptom: expected columns weren’t added

Confirm:
- the session actually contains the relevant chunks:
  ```python
  bundle.list_chunk_ids()
  ```
- you didn’t skip due to mismatch:
  - set `return_summary=True` and inspect `summary`.

## Next steps

- Session bundle capture: {doc}`10_session_bundles_get_session_bundle`
- Session ↔ AnnData conceptual bridge: {doc}`../b_concepts_mental_models/05_sessions_to_anndata_bridge`
