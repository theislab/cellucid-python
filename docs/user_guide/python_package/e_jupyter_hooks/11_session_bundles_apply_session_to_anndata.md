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
the `cellOrder` section below.
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
- a field declared with source `var` copies that gene's values across cells:
  the gene is found by name among the gene names the session's `dataMapping`
  shows (made unique with the shared rule), and its values are read from the
  mapping's expression matrix — `X`, `raw.X`, or a layer.

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

A session bundle's `datasetFingerprint` has exactly six required fields:

| Field | Checked by `apply_cellucid_session_to_anndata` |
|---|---|
| `sourceType` | must be a non-empty string |
| `datasetId` | must equal `expected_dataset_id` |
| `cellCount` | must equal `adata.n_obs` |
| `dataMapping` | every key it names must exist on `adata` |
| `varCount` | must equal the gene count of the mapped axis: `adata.raw.n_vars` for `expression: "raw"`, `0` when the mapping reads no expression, else `adata.n_vars` |
| `cellOrder` | shape only |

Application requires `expected_dataset_id`. Counts, the dataset ID, and the
mapping must match exactly; a mismatch raises `ValueError` before the target is
mutated.

In practice, if you subset/shuffle your AnnData after the session was created, you will almost always have a mismatch.

### `dataMapping`: which arrays the session showed

`dataMapping` records which stored arrays the session was showing. Its `kind`
is `"export"` — exactly `kind`, `embeddings` (the `1d`/`2d`/`3d` slots),
`geneNames` (`"index"` or `"column:<name>"`), and `vectorFields` — or
`"anndata"`, which adds `knn`, `expression` (`"X"`, `"raw"`, `"layer:<name>"`,
or null), and `latent`. For an AnnData mapping, Python requires each named `obsm`
embedding to exist and be wide enough for its columns, and the named `obsp`
graph, latent space, layer or `raw`, and gene-name column to exist. An export
mapping names exported embedding ids, which an `AnnData` need not carry, so only
its gene-name column is checked.

### `cellOrder`: shape-checked here, enforced in the viewer

`cellOrder` is exactly `embedding`, `dimension` (1, 2, or 3), and a
16-character lowercase hexadecimal `digest` of the coordinates of that
embedding. `embedding` is `{"key", "columns"}`: an exported embedding id with
`columns: null`, or an AnnData `obsm` key with the 0-based columns that were
read, one per dimension. It is what stops a dataset republished under the same
name with the same cell and gene counts but a different row order from silently
re-pointing every saved selection at different cells.

- **It is required.** A bundle without it is refused.
- **Its shape is checked strictly.** Anything other than the three keys above,
  a `dimension` outside 1–3, a malformed `digest`, or an `embedding` whose
  columns do not fit its kind raises before payload decoding.
- **Its digest is never verified.** It is taken over the normalized coordinates
  the viewer drew, and `apply_session_to_anndata` writes into an `AnnData` whose
  row order the caller controls. There is nothing here to check it against.

Enforcing cell-order identity is therefore the viewer's job, not this reader's,
so a bundle this reader accepts is not necessarily one the web app will open.

When `store_uns=True`, the fingerprint is copied verbatim into
`adata.uns["cellucid"]["session"]["dataset_fingerprint"]`, `cellOrder` and
`dataMapping` included.

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
