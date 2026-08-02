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
- Developers: read the identity checks and column naming rules.

## Quickstart (copy/paste)

### One-liner (recommended)

```python
adata2 = viewer.apply_session_to_anndata(adata, inplace=False)
```

### Explicit two-step

```python
bundle = viewer.get_session_bundle(timeout=60)
adata2 = bundle.apply_to_anndata(
    adata,
    expected_dataset_id="my-study-v1",
    inplace=False,
)
```

```{important}
Session application is currently **index-based**:

The session references cells by **row index** (0-based). Only apply a session to an `AnnData` whose row order matches the dataset that produced the session.

Nothing in Python can verify that for you — see
[`cellOrder`: recorded here, enforced in the viewer](#cellorder-recorded-here-enforced-in-the-viewer).
```

## What gets applied (current behavior)

As implemented in `cellucid-python/src/cellucid/anndata_session/`, applying a session can add:

### 1) Highlight groups → boolean `adata.obs` columns

If the bundle contains highlight metadata:
- For each highlight group, a boolean mask column is created in `adata.obs`.
- Default prefix: `cellucid_highlight__`

Example output:

```text
adata.obs["cellucid_highlight__highlight_12"]  # True/False per cell
```

### 2) User-defined fields → `adata.obs` columns

If the bundle contains user-defined fields:
- categorical fields decode exact labels and cell-aligned codes;
- continuous fields copy one exact cell-aligned source;
- a field declared with source `var` copies that gene's values across cells.

All materialized fields are cell-aligned `adata.obs` columns.

### 3) Provenance metadata → `adata.uns["cellucid"]`

By default, application stores:
- the session manifest,
- dataset fingerprint,
- and “applied settings”
under:

```text
adata.uns["cellucid"]["session"]
```

## Exact dataset identity (avoid accidental corruption)

A session bundle's `datasetFingerprint` has four required fields and one
optional one:

| Field | Required | Checked by `apply_cellucid_session_to_anndata` |
|---|---|---|
| `sourceType` | yes | must be a non-empty string |
| `datasetId` | yes | must equal `expected_dataset_id` |
| `cellCount` | yes | must equal `adata.n_obs` |
| `varCount` | yes | must equal `adata.n_vars` |
| `cellOrder` | no | shape only, when present |

Application requires `expected_dataset_id`. Cell/variable counts and the dataset
ID must match exactly; a mismatch raises `ValueError` naming every field that
differs, before the target is mutated.

In practice, if you subset/shuffle your AnnData after the session was created, you will almost always have a mismatch.

### `cellOrder`: recorded here, enforced in the viewer

`cellOrder` is the viewer's record of which cell ordering a session was saved
against — a digest of the coordinates on screen at save time, plus the
embedding dimension (`1`, `2`, or `3`) the digest was taken in. It is what stops
a dataset republished under the same name with the same cell and gene counts but
a different row order from silently re-pointing every saved selection at
different cells.

This reader treats it asymmetrically, and the asymmetry is deliberate:

- **It is accepted, not required.** A bundle validates with or without
  `cellOrder`.
- **Its shape is checked strictly whenever it is present.** Exactly `dimension`
  and `digest`; `dimension` exactly the integer 1, 2, or 3; `digest` exactly 16
  lowercase hexadecimal characters. Anything else raises before payload
  decoding.
- **Its value is never verified.** The digest is taken over exported
  coordinates, which an `AnnData` need not contain, and
  `apply_session_to_anndata` writes into an `AnnData` whose row order the caller
  controls. There is nothing here to check it against.

Enforcing cell-order identity is therefore the viewer's job, not this reader's.

```{important}
A bundle this reader accepts is not necessarily one the web app will open.

The web app **requires** `cellOrder` and refuses a session file saved before the
record existed, on the reasoning that a file recording selections as bare row
indices can never be shown to still mark the same cells. If you hold a bundle
without `cellOrder`, Python will still apply it — under the row-order assumption
below, which is now yours alone to uphold — but the viewer will not load it.
Re-create the state in the app and save a fresh session file.
```

When `store_uns=True`, the fingerprint is copied verbatim into
`adata.uns["cellucid"]["session"]["dataset_fingerprint"]`, `cellOrder` included,
so the provenance of an applied session survives even though this reader did not
act on it.

## Column naming and conflicts

Applying a session may add columns that collide with existing names. Any
existing target column raises `ValueError`; Cellucid does not rename or
overwrite it.

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
    expected_dataset_id="my-study-v1",
    inplace=False,
    add_highlights=True,
    add_user_defined_fields=True,
    store_uns=True,
    return_summary=True,
)
```

This returns an `ApplySummary` with:
- `added_obs_columns`: the exact ordered names materialized by the successful
  application.

## Edge cases

- The default `inplace=False` path rejects backed AnnData before calling
  `adata.copy()`. Explicitly materialize the intended target with
  `adata.to_memory()` first.
- `inplace=True` preserves backed `X` and changes the live object's in-memory
  `obs`/`uns`; it does not write those changes into the source H5AD.
- Pandas is a core `cellucid` dependency and writes the requested categorical
  and boolean columns.
- If a highlight group contains out-of-range indices, application raises.
- Every manifest chunk must belong to the closed current contributor/chunk
  inventory documented in
  {doc}`../g_api_reference_coverage/api/sessions`; unknown or mismatched entries
  raise before payload decoding.

## Troubleshooting

### Symptom: “Dataset fingerprint mismatch”

This is almost always because:
- you filtered/subset/reordered `adata` after opening the viewer, or
- you’re applying the session to a different dataset.

Fix:
- apply to the exact same AnnData (same row order), and
- pass the exact dataset ID used to create the viewer.

### Symptom: expected columns weren’t added

Confirm:
- the session actually contains the relevant chunks:
  ```python
  bundle.list_chunk_ids()
  ```
- the dataset ID and dimensions match exactly.

## Next steps

- Session bundle capture: {doc}`10_session_bundles_get_session_bundle`
- Session ↔ AnnData conceptual bridge: {doc}`../b_concepts_mental_models/05_sessions_to_anndata_bridge`
