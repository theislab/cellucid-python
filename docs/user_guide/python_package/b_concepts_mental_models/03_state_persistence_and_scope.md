# State, persistence, and scope

**Audience:** everyone (wet lab → developer)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- The three persistence layers (live state vs session file vs dataset file)
- What lives in Python vs in the browser vs on disk
- What tends to reset when you reload/re-run
- The safest way to make interactive work reproducible

---

## Mental model (one sentence)

Cellucid has **live state** (changes while you interact), **durable state** (session files), and **reproducible data artifacts** (Cellucid-ready AnnData files) — and you should choose the right one depending on whether you’re exploring, sharing, or publishing.

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
- the notebook output that shows the viewer is removed or the notebook reloads.

### 2) Session file (`.cellucid-session`) (durable, shareable)

A session file is a single file that holds **the history of the app state**
(not the dataset itself): every step taken, by a person, an agent or a
notebook, and the state each step left behind, in one or more trackers.

You can create it:
- in the web app via **Share…** under **State › Session** (downloads a file of
  the history, the current state only, or the steps only, as chosen), or
- in notebooks via `viewer.request_session()` and, in the next cell,
  its `result()` (no-download capture).

You can use it to:
- reopen your work later, and go back to any step,
- send “this exact view/highlight/filter setup” to a collaborator,
- apply the state at any step back onto an AnnData in Python (see {doc}`05_sessions_to_anndata_bridge`).

### 3) Dataset file (reproducible dataset artifact)

A Cellucid-ready AnnData (from `prepare(...)`) is a versioned dataset file designed for:
- fast loading,
- stable behavior across machines,
- long-term sharing/hosting,
- paper-ready reproducibility.

The file carries its stable identity in `uns['cellucid']` (see {doc}`04_dataset_identity_and_reproducibility`).

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

If you need durable UI state in Python, you must capture it via a **session file**.

### Browser-side state (what the web app “knows”)

The web app owns the interactive UI state and may also keep some small preferences in browser storage.

Practical implication:
- if you refresh/reload, don’t assume “everything comes back” unless you saved a session.

### Disk-side artifacts (what survives across machines)

- Dataset file: stable dataset artifact, shareable.
- Session file: stable state artifact (its whole history), shareable.
- Notebook code + environment: the “how to reproduce” part.

---

## Scope: dataset-dependent vs dataset-agnostic state

This distinction matters both for **restoring sessions** and for **applying sessions to AnnData**:

- **Dataset-dependent state** (must match the same dataset):
  - highlights (membership lists),
  - filters based on obs/var fields,
  - user-defined codes aligned to cell indices,
  - anything that references specific cell indices or gene indices.

- **Dataset-agnostic state** (no reference to cells or genes):
  - some UI layout/panel geometry,
  - some global preferences,
  - some non-data-specific settings.

A session restores only onto the dataset it was saved from. If its dataset
fingerprint does not match the loaded dataset, the whole restore is rejected
and rolled back; no dataset-agnostic subset is salvaged. The obs fields and
genes the session uses are checked, each missing one named, before anything
applies.

---

## A reproducible workflow (recommended)

This is a practical, high-signal pattern that works for both beginners and experts.

### Step 1: prepare (once)

```python
from cellucid import prepare
prepare(
    "pbmc_v1.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    dataset_id="pbmc-v1",
    dataset_name="PBMC v1",
    obs=adata.obs,
)
```

### Step 2: explore (many times)

```python
from cellucid import show_anndata
viewer = show_anndata("pbmc_v1.h5ad")
```

### Step 3: save interactive state (when it matters)

In notebooks (two cells: the widget hears the viewer between cells):

```python
request = viewer.request_session()
```

```python
session = request.result()
session.save("./sessions/pbmc_v1_my_findings.cellucid-session")
```

In the web app: use **Share…** (downloads a `.cellucid-session` file).

### Step 4: materialize key decisions back into AnnData (optional)

For example, turn highlights into `adata.obs` columns so downstream analysis can be scripted and versioned:

```python
plan = session.apply(adata)      # what would change, one row each; nothing written yet
report = plan.apply()            # writes the rows into adata and records the apply
```

See: {doc}`05_sessions_to_anndata_bridge`.

---

## Screenshot: where Share… and Load… live in the web app

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The State section with the Pancreatic endocrinogenesis sample open: Session with Share… and Load…, and the History with the recording switch, the branch shown and the graph of the steps taken so far, the step shown marked.
:width: 540px

The **State** section of the sidebar keeps **Share…** and **Load…** under **Session**, above the history they write and read.
```

---

## Edge cases and caveats (read this if you share sessions)

### Index-based identity (today)

Selections, highlights and field values name cells by **original position** (row positions). If you apply a session to an AnnData whose row order differs, your highlights would point to different cells. Applying a session refuses such an object when it has the embedding the session's cell-order digest was taken on; without that embedding, the order is not checked.

Mitigations:
- keep dataset files versioned and immutable once shared,
- treat cell order as part of dataset identity,
- give each dataset an id (`uns['cellucid']['id']`): applying a session compares it with the session's.

### “Why was my session refused?”

A session restores whole or not at all. Most commonly:
- dataset mismatch (different cell count / var count / dataset id / cell order),
- fields or genes the session uses were renamed or deleted since it was
  created (the error names each missing one),
- the session was created by a version with a different session format.

### “Why did my notebook re-run ‘forget’ the viewer?”

The notebook code is not the viewer. Re-running cells creates new viewer instances unless you capture state as a session file.

---

## Troubleshooting

### Symptom: “I refreshed and lost my highlights/filters”

Likely causes:
- the page opened its default dataset after the refresh, not yours;
- **Autosave in this browser** (State › History settings) is off, or the browser keeps no data
  for the site (some private windows);
- a notebook viewer: a new viewer starts a new history.

Fix:
- open the same dataset again and choose **Continue** when Cellucid offers to
  continue where you left off (or **State › Session › Load…**):
  every step is kept in this browser as you make it
  ({doc}`../../web_app/m_state_history/03_saved_states_and_autosave`);
- keep a session file for anything you need on another machine: **Share…**
  under **State › Session** in the web app, or `viewer.request_session()`
  in Python, then load it after opening the dataset.

### Symptom: “I loaded a session but nothing changed”

Likely causes:
- you loaded the session before the dataset finished loading,
- the session was rejected because its dataset identity did not match,
- the session load reported another validation error.

How to confirm:
- read the session-load error shown by the UI,
- verify dataset ID and counts match (see {doc}`04_dataset_identity_and_reproducibility`).

---

## See also

- Web app session mental model: {doc}`../../web_app/l_sessions_sharing/01_session_mental_model`
- What gets saved/restored: {doc}`../../web_app/l_sessions_sharing/02_what_gets_saved_and_restored`

## Next steps

- Stable IDs and versioning: {doc}`04_dataset_identity_and_reproducibility`
- Pull sessions into Python: {doc}`05_sessions_to_anndata_bridge`
