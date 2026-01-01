# Sessions → AnnData (No-Download Bridge)

**Audience:** computational users + power users (wet lab users: skim the “Fast path”)  
**Time:** 20–40 minutes  
**What you’ll learn:**
- How to capture the current viewer state into Python **without** a manual browser download
- What parts of a session can be applied back to AnnData today
- The safety model (dataset mismatch policies + index-based identity)
- How to debug “bundle capture is stuck” and “apply did nothing”

**Prerequisites:**
- A viewer running in a notebook (`show(...)` or `show_anndata(...)`)
- An `AnnData` object you want to annotate (in memory)

For deeper design rationale and file-format notes, see {doc}`05_sessions_to_anndata_design`.

---

## Mental model: “live state” vs “durable state”

Cellucid exposes two complementary layers of state in notebooks:

1) **`viewer.state` (live snapshot; volatile)**  
   A small “latest event” snapshot updated from UI events (ready/selection/hover/click/…).

2) **Session bundle (`.cellucid-session`; durable artifact)**  
   A file-format snapshot of meaningful UI state (highlights, user-defined fields, overlay registries, layout state).
   In notebooks you can request it as a `CellucidSessionBundle` object.

If you want to **reproduce** or **transfer** interactive state, use the **session bundle**.

---

## Fast path (wet lab / non-technical): save your work into the data object

If you don’t want to think about file formats:

1) Open your data in Cellucid (in Jupyter).
2) Create highlights/labels in the viewer (select cells, label groups, etc.).
3) In Python, run:

   ```python
   adata2 = viewer.apply_session_to_anndata(adata, inplace=False)
   ```

4) Use `adata2` for downstream analysis (Scanpy plots, marker genes, etc.).

This turns interactive choices (highlights/labels) into explicit columns you can store and version.

---

## Step-by-step workflow: capture → (optionally save) → apply

### Step 1: create and load a viewer

```python
from cellucid import show_anndata

viewer = show_anndata("data.h5ad", height=650)
viewer.wait_for_ready(timeout=60)
```

### Step 2: interact in the UI

In the viewer:
- color by a field,
- select cells,
- create highlight groups,
- create user-defined categorical fields (if you use that UI).

### Step 3: capture the session bundle into Python (no browser download)

```python
bundle = viewer.get_session_bundle(timeout=60)
bundle
```

The returned object is a handle to a `.cellucid-session` file written to a **temporary file** on disk.

### Step 4 (optional but recommended): persist the artifact

If you want to keep a copy of “what I saw in the UI”:

```python
bundle.save("./sessions/my_state.cellucid-session")
```

### Step 5: apply to AnnData (non-destructive by default)

```python
adata2 = bundle.apply_to_anndata(adata, inplace=False)
```

Or use the one-liner (captures + applies + cleans up the temp bundle file):

```python
adata2 = viewer.apply_session_to_anndata(adata, inplace=False)
```

---

## What gets applied to AnnData (today)

`apply_cellucid_session_to_anndata(...)` currently materializes:

### 1) Highlights → `adata.obs` boolean columns

- One boolean `.obs` column per highlight group.
- Default column naming: `cellucid_highlight__<groupId>` (sanitized to safe characters).

This is the “most useful” mapping because it makes interactive groups explicit and scriptable:
- `adata2.obs["cellucid_highlight__my_group"]` becomes a mask you can use in analysis.

### 2) User-defined categorical fields → `adata.obs` or `adata.var`

- Categorical columns decoded from session chunks like `user-defined/codes/<fieldId>`.
- Category labels come from `core/field-overlays`.
- The target is controlled by the field’s `source` (`obs` vs `var`).

### 3) Metadata stored under `adata.uns`

If `store_uns=True` (default), the apply step stores:
- `adata.uns["cellucid"]["session"]["manifest"]`
- `adata.uns["cellucid"]["session"]["dataset_fingerprint"]`
- `adata.uns["cellucid"]["session"]["applied"]` (apply options)
- some decoded JSON chunks for debugging

```{important}
Not everything in the session bundle becomes a clean AnnData mutation.
The goal is to materialize the parts that make sense as dataset annotations (highlights, categorical labels).
UI-only state (camera, panel layout, view geometry) generally stays in the session bundle, not in AnnData.
```

---

## Safety model: dataset mismatch + index-based identity

### Sessions are treated as untrusted input

Even if the session came from “your own browser”, treat it like external input:
- bounds checks on indices and code arrays,
- file framing validation (MAGIC header),
- size caps to avoid runaway memory usage.

### Dataset mismatch policy (shape + optional ID)

Before applying dataset-dependent chunks, Python checks a lightweight fingerprint stored in the bundle:
- `cellCount` vs `adata.n_obs`
- `varCount` vs `adata.n_vars`
- optionally `datasetId` vs `expected_dataset_id`

Policy (`dataset_mismatch`):
- `"error"`: raise on mismatch (strictest; best for pipelines)
- `"warn_skip"`: warn and skip dataset-dependent chunks (default; safest interactive behavior)
- `"skip"`: silently skip dataset-dependent chunks

### Critical constraint: cell identity is index-based (today)

Highlights and user-defined codes are aligned to **cell indices** (row positions).

Only apply a session to an `AnnData` whose row order matches the dataset that produced the session.

Recommended extra guard:
- pass `expected_dataset_id=...` when you have a stable dataset ID string (see {doc}`04_dataset_identity_and_reproducibility`).

---

## Advanced usage

### Apply from a saved `.cellucid-session` file (web app or notebook)

```python
from cellucid import CellucidSessionBundle

bundle = CellucidSessionBundle("./sessions/my_state.cellucid-session")
adata2 = bundle.apply_to_anndata(adata, inplace=False)
```

### Control column naming + conflict behavior

If you apply multiple sessions, you may create naming conflicts.

Key knobs:
- `highlights_prefix=...`
- `user_defined_prefix=...`
- `column_conflict="error" | "overwrite" | "suffix"` (default: `"suffix"`)

Example:

```python
adata2 = bundle.apply_to_anndata(
    adata,
    inplace=False,
    highlights_prefix="hl__",
    user_defined_prefix="udf__",
    column_conflict="suffix",
)
```

### Get a structured apply summary

```python
adata2, summary = bundle.apply_to_anndata(adata, inplace=False, return_summary=True)
summary
```

---

## Troubleshooting (large, symptom-based)

### Symptom: `viewer.get_session_bundle()` times out

Likely causes (ordered):
1) The viewer never reached “ready” (dataset load stuck).
2) The notebook iframe is not displayed (no frontend to respond).
3) The browser cannot POST to the server endpoint `/_cellucid/session_bundle` (proxy/mixed content/remote kernel).

How to confirm:
- Run `viewer.debug_connection()` and inspect:
  - `state.ready`
  - `frontend_roundtrip`
  - `viewer_index_probe`
  - `frontend_console`
- In the browser devtools network tab:
  - look for a POST to `/_cellucid/session_bundle?...`

Fix:
- Call `viewer.display()` to ensure the iframe exists.
- Call `viewer.wait_for_ready(timeout=60)` and only then request a bundle.
- If remote: use SSH tunneling or `jupyter-server-proxy` (see {doc}`../../web_app/b_data_loading/05_jupyter_tutorial`).

### Symptom: “No pending session bundle request (viewerId/requestId)”

This usually means:
- the upload arrived late (after the request TTL expired), or
- the viewer tab is from an old kernel/viewerId, or
- you copied a URL across environments and the viewerId no longer matches.

Fix:
- refresh the viewer (or recreate it),
- request a new bundle again.

### Symptom: `apply_to_anndata(...)` applies nothing / skips due to mismatch

Likely causes:
- `adata.n_obs` / `adata.n_vars` differ from the session’s fingerprint,
- `expected_dataset_id` mismatched (or inferred ID differs),
- the AnnData was subset/reordered since the session was created.

How to confirm:
- Use `return_summary=True` and inspect `summary.mismatch_reasons`.
- Inspect `bundle.dataset_fingerprint`.

Fix:
- apply to the matching AnnData version (same row order),
- or treat the mismatch as intentional and accept that dataset-dependent chunks are skipped.

### Symptom: “My highlight groups didn’t become columns”

Likely causes:
- highlights are missing from the session (never created, or not saved by that feature yet),
- dataset mismatch caused highlights to be skipped,
- the session chunk IDs differ across versions.

How to confirm:
- `bundle.list_chunk_ids()` and look for:
  - `highlights/meta`
  - `highlights/cells/<groupId>`

Fix:
- recreate highlights in the viewer and capture a fresh bundle,
- ensure dataset identity matches.

### Symptom: “Applying user-defined fields fails or gives weird categories”

Likely causes:
- the field is not categorical (only categorical is applied today),
- the categories list is missing,
- codes contain out-of-range values (the apply code clamps invalid codes to NA).

Fix:
- confirm the field is categorical in the UI,
- apply with `return_summary=True` and inspect which columns were added.

### Symptom: ImportError about pandas

`apply_cellucid_session_to_anndata` requires `pandas` because it writes `.obs`/`.var` columns.

Fix:
- `pip install pandas`

---

## Next steps

- Deep design and file format: {doc}`05_sessions_to_anndata_design`
- Dataset identity and versioning: {doc}`04_dataset_identity_and_reproducibility`
