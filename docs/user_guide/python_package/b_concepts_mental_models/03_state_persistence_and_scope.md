# State, persistence, and scope

**Audience:** everyone (wet lab → developer)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- The three persistence layers (live state vs session bundle vs export folder)
- What lives in Python vs in the browser vs on disk
- What tends to reset when you reload/re-run
- The safest way to make interactive work reproducible

---

## Mental model (one sentence)

Cellucid has **live state** (changes while you interact), **durable state** (session bundles), and **reproducible data artifacts** (export folders) — and you should choose the right one depending on whether you’re exploring, sharing, or publishing.

---

## The three persistence layers (this prevents most confusion)

### 1) Live, in-memory state (volatile)

This is “what you see right now” in the viewer:
- camera position,
- active color-by field,
- filters,
- highlights,
- analysis panel contents,
- view layout (live view + snapshots),
- etc.

Live state is lost if:
- the tab crashes,
- you refresh the page,
- you close the browser,
- your notebook output iframe is destroyed.

### 2) Session bundle (`.cellucid-session`) (durable, shareable)

A session bundle is a single file that captures **a portable snapshot of the app state** (not the dataset itself).

You can create it:
- in the web app via **Save State** (downloads a file), or
- in notebooks via `viewer.get_session_bundle()` (no-download capture).

You can use it to:
- reopen your work later,
- send “this exact view/highlight/filter setup” to a collaborator,
- apply selected parts of state back onto an AnnData in Python (see {doc}`05_sessions_to_anndata_bridge`).

### 3) Export folder (reproducible dataset artifact)

An export folder (from `prepare(...)`) is a versioned dataset representation designed for:
- fast loading,
- stable behavior across machines,
- long-term sharing/hosting,
- paper-ready reproducibility.

The export folder is where you put stable identity metadata like `dataset_identity.json` (see {doc}`04_dataset_identity_and_reproducibility`).

---

## Where state lives (Python vs browser vs disk)

### Python-side state (what the notebook “knows”)

The Python `viewer` object owns:
- connection identifiers (`viewerId`, `viewerToken`),
- hook registrations (`@viewer.on_selection`, etc.),
- a small latest-event snapshot (`viewer.state`),
- the lifetime of the local server (`viewer.stop()`).

Python does **not** automatically have:
- the full viewer UI state (filters, highlights, camera),
- your snapshot layout,
- your analysis panel state.

If you need durable UI state in Python, you must capture it via a **session bundle**.

### Browser-side state (what the web app “knows”)

The web app owns the interactive UI state and may also keep some small preferences in browser storage.

Practical implication:
- if you refresh/reload, don’t assume “everything comes back” unless you saved a session.

### Disk-side artifacts (what survives across machines)

- Export folder: stable dataset artifact, shareable.
- Session bundle: stable state artifact, shareable.
- Notebook code + environment: the “how to reproduce” part.

---

## Scope: dataset-dependent vs dataset-agnostic state

This distinction matters both for **restoring sessions** and for **applying sessions to AnnData**:

- **Dataset-dependent state** (must match the same dataset):
  - highlights (membership lists),
  - filters based on obs/var fields,
  - user-defined codes aligned to cell indices,
  - anything that references specific cell indices or gene indices.

- **Dataset-agnostic state** (can sometimes apply even if dataset differs):
  - some UI layout/panel geometry,
  - some global preferences,
  - some non-data-specific settings.

If a session detects a dataset mismatch, Cellucid intentionally restores only the safer “agnostic” subset.

---

## A reproducible workflow (recommended)

This is a practical, high-signal pattern that works for both beginners and experts.

### Step 1: export (once)

```python
from cellucid import prepare
prepare(..., out_dir="./exports/pbmc_v1")
```

### Step 2: explore (many times)

```python
from cellucid import show
viewer = show("./exports/pbmc_v1")
viewer.wait_for_ready(timeout=60)
```

### Step 3: save interactive state (when it matters)

In notebooks:

```python
bundle = viewer.get_session_bundle(timeout=60)
bundle.save("./sessions/pbmc_v1_my_findings.cellucid-session")
```

In the web app: use **Save State** (downloads a `.cellucid-session` file).

### Step 4: materialize key decisions back into AnnData (optional)

For example, turn highlights into `adata.obs` columns so downstream analysis can be scripted and versioned:

```python
adata2 = bundle.apply_to_anndata(adata, inplace=False)
```

See: {doc}`05_sessions_to_anndata_bridge`.

---

## Screenshot: where “Save State” lives in the web app

If you want to make this chapter friendlier for wet lab users, capture a screenshot of the session controls.

<!-- SCREENSHOT PLACEHOLDER
ID: state-persistence-save-state-button
Suggested filename: sessions_sharing/save-state-button-location.png
Where it appears: Python Package → Concepts → State, persistence, and scope → Session bundle section
Capture:
  - UI location: Cellucid web app, “Session state” panel
  - State prerequisites: any dataset loaded
  - Action to reach state: open the “User data” area and locate “Session state”; show the “Save State” button
Crop:
  - Include: the panel header + the Save/Load controls (so the reader can find it)
  - Exclude: private dataset identifiers
Redact:
  - Remove: dataset names if sensitive
Annotations:
  - Callouts: #1 Save State button, #2 Load State button
Alt text:
  - Session state panel with Save State and Load State controls.
Caption:
  - Session bundles are the durable way to preserve interactive UI state; you can save/load them in the web app or capture them into Python in notebooks.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Session state panel.
:width: 100%

Session bundles are the durable way to preserve interactive UI state across reloads and across machines.
```

---

## Edge cases and caveats (read this if you share sessions)

### Index-based identity (today)

Selections/highlights/codes are aligned to **cell indices** (row positions). If you apply a session to an AnnData whose row order differs, your highlights can silently point to different cells.

Mitigations:
- keep exports versioned and immutable once shared,
- treat cell order as part of dataset identity,
- use `expected_dataset_id=...` when applying sessions in Python.

### “Why did my session restore only part of the state?”

Most commonly:
- dataset mismatch (different cell count / var count / dataset id),
- fields renamed/deleted since the session was created,
- the session was created on an older/newer version with different feature coverage.

### “Why did my notebook re-run ‘forget’ the viewer?”

The notebook code is not the viewer. Re-running cells creates new viewer instances unless you capture state as a session bundle.

---

## Troubleshooting

### Symptom: “I refreshed and lost my highlights/filters”

Likely causes:
- you did not save a session bundle,
- you saved a session but did not reload it after refresh.

Fix:
- save state in the web app (Save State) or in Python (`viewer.get_session_bundle()`),
- reload the session after reloading the dataset.

### Symptom: “I loaded a session but nothing changed”

Likely causes:
- you loaded the session before the dataset finished loading,
- the session was for a different dataset and was skipped,
- the session restore was still in progress (large sessions can restore progressively).

How to confirm:
- watch for session-related warnings/notifications in the UI,
- verify dataset ID and counts match (see {doc}`04_dataset_identity_and_reproducibility`).

---

## See also

- Web app session mental model: {doc}`../../web_app/l_sessions_sharing/01_session_mental_model`
- What gets saved/restored: {doc}`../../web_app/l_sessions_sharing/02_what_gets_saved_and_restored`

## Next steps

- Stable IDs and versioning: {doc}`04_dataset_identity_and_reproducibility`
- Pull sessions into Python: {doc}`05_sessions_to_anndata_bridge`
