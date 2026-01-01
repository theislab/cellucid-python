# Session persistence and restoring analysis artifacts

This tutorial explains Cellucid **session bundles** and how to use them as a bridge:

Viewer state → `.cellucid-session` file → Python → `AnnData`

This is the workflow you want when:
- you made careful selections/highlight groups in the UI and want them in your analysis code
- you want to reproduce an interactive exploration later
- you want to share “exactly what I did in the viewer” with a collaborator

---

## At a glance

**Audience**
- Computational users (power users)
- Developers (protocol-level details)

**Time**
- Minimal “capture a session bundle”: ~10–15 minutes
- Full “apply to AnnData and analyze”: ~30–60 minutes

**Prerequisites**
- A notebook environment
- A viewer (`show(...)` or `show_anndata(...)`)

---

## What a session bundle is

A session bundle is a single file with extension:
- `.cellucid-session`

It contains:
- a manifest (`datasetFingerprint`, chunk list, codecs)
- a sequence of chunks (JSON or binary; optionally gzip-compressed)

In Python, it is represented by:
- `cellucid.CellucidSessionBundle`

```{note}
Session bundles are treated as **untrusted input** in Python: the reader and AnnData-applier do bounds checks and have dataset mismatch policies.
```

---

## Two ways to capture a session

### Option A: UI download (web app workflow)

Use the viewer’s UI to save/download a `.cellucid-session`.

See web-app docs:
- {doc}`../../web_app/l_sessions_sharing/index`

### Option B (recommended in notebooks): capture from Python (no browser download)

Cellucid can request the session from the iframe and receive it back over HTTP.

---

## Step-by-step: capture a session bundle in Python

### Step 1 — Open a viewer

```python
from cellucid import show_anndata

viewer = show_anndata(adata, height=650)
viewer
```

### Step 2 — Create some meaningful state in the UI

Examples:
- create highlight groups (select cells, then “confirm” into a group)
- create user-defined categorical fields (if supported by the UI)

For the highlight mental model and UI walkthrough:
- {doc}`../../web_app/f_highlighting_selection/index`

### Step 3 — Capture the bundle from Python

```python
bundle = viewer.get_session_bundle(timeout=60)
bundle
```

If you want to block explicitly until the viewer is ready:

```python
viewer.wait_for_ready(timeout=30)
bundle = viewer.get_session_bundle(timeout=60)
```

### Step 4 — Inspect what’s inside

```python
bundle.manifest.keys()
bundle.dataset_fingerprint
bundle.list_chunk_ids()
```

---

## Save the bundle as an artifact (recommended)

The session bundle returned by `get_session_bundle()` may be stored in a temporary location.
Save it to a stable path:

```python
saved_path = bundle.save("sessions/my_analysis.cellucid-session")
saved_path
```

---

## Apply a session bundle to an AnnData (the “session → analysis” bridge)

### What gets applied today (important)

The current Python applier focuses on data that is meaningful in analysis contexts:

- **Highlights** → new boolean `adata.obs` columns
  - one column per highlight group
- **User-defined categorical fields** → new categorical columns in `adata.obs` and/or `adata.var`

Other session state (camera position, UI layout, etc.) is not applied to AnnData (it’s UI state, not analysis state).

### Apply with the convenience API

```python
adata2 = viewer.apply_session_to_anndata(adata, inplace=False)
adata2
```

### Apply with full control

```python
from cellucid import apply_cellucid_session_to_anndata

adata2, summary = apply_cellucid_session_to_anndata(
    bundle,
    adata,
    inplace=False,
    dataset_mismatch="warn_skip",          # "error" | "warn_skip" | "skip"
    expected_dataset_id=None,             # set if you want strict dataset-id matching
    add_highlights=True,
    highlights_prefix="cellucid_highlight__",
    add_user_defined_fields=True,
    user_defined_prefix="",
    include_deleted_user_defined_fields=False,
    store_uns=True,                       # store manifest + chunk metadata under adata.uns['cellucid']
    column_conflict="suffix",             # "error" | "overwrite" | "suffix"
    return_summary=True,
)

summary
```

### How to use the result in analysis

List added columns:

```python
summary.added_obs_columns
```

Subset by a highlight group:

```python
highlight_col = summary.added_obs_columns[0]
selected = adata2[adata2.obs[highlight_col].astype(bool)].copy()
selected
```

---

## Dataset mismatch policies (read this!)

Session bundles include a dataset fingerprint (cell count, var count, dataset id).

If you apply a session to the “wrong” dataset, one of these should happen:
- **error**: stop immediately (strict)
- **warn_skip** (default): warn and skip dataset-dependent chunks
- **skip**: silently skip dataset-dependent chunks

Best practice:
- when you care about reproducibility, set `expected_dataset_id` and use `dataset_mismatch="error"`

---

## Screenshot placeholders (optional but helpful)

<!-- SCREENSHOT PLACEHOLDER
ID: python-notebooks-session-highlight-groups-before-capture
Suggested filename: sessions_sharing/00_highlight-groups-before-session-capture.png
Where it appears: User Guide → Python Package → Notebooks/Tutorials → 32_session_persistence_and_restoring_analysis_artifacts.md
Capture:
  - UI location: highlight UI panel (pages/groups visible)
  - State prerequisites: at least two highlight groups exist (created from selections)
  - Action to reach state: create two highlight groups and open the highlight panel
Crop:
  - Include: highlight page/group list and group names/colors
  - Include: enough canvas to show highlighted cells
Redact:
  - Remove: private dataset names
Alt text:
  - Highlight panel showing multiple named highlight groups before capturing a session bundle.
Caption:
  - Session bundles can include highlight groups; capturing the session lets you bring those groups back into Python as `adata.obs` columns.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for highlight groups before session capture.
:width: 100%

Highlight groups created in the UI (these can be captured in a session bundle).
```

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: `viewer.get_session_bundle()` times out

Likely causes:
- viewer never reached “ready”
- iframe cannot communicate to Python server (proxy/tunnel issue)

How to confirm:
```python
viewer.debug_connection()
```

Fix:
- wait for ready: `viewer.wait_for_ready()`
- ensure remote notebooks use a proxy/tunnel (see {doc}`22_large_dataset_server_mode_and_lazy_gene_expression`)

### Symptom: “Applying session adds no columns”

Likely causes:
- the session did not contain highlight groups or user-defined fields
- dataset mismatch caused dataset-dependent chunks to be skipped

How to confirm:
```python
bundle.list_chunk_ids()
bundle.dataset_fingerprint
```

Fix:
- create highlight groups in the UI before capturing
- use `dataset_mismatch="error"` while developing to catch mismatch early

### Symptom: “Column names look weird / got suffixed”

This is expected when a column name already exists:
- policy `column_conflict="suffix"` creates `name__2`, `name__3`, ...

If you want strict behavior:
- use `column_conflict="error"`

---

## Next steps

- {doc}`05_jupyter_embedding_hooks_sessions_gallery` (runnable session notebooks)
- {doc}`23_programmatic_highlighting_and_selection_callbacks` (robust hooks patterns)
- Web-app sessions: {doc}`../../web_app/l_sessions_sharing/index`
