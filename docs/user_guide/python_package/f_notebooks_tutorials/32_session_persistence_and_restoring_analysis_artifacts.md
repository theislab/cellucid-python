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
- A viewer (`show_anndata(...)`)

---

## What a session bundle is

A session bundle is a single file with extension:
- `.cellucid-session`

It contains:
- a manifest (`datasetFingerprint`, `source`, `requires` — the obs columns and
  genes the session uses — and the chunk list with codecs)
- a sequence of chunks (JSON or binary; optionally gzip-compressed)

In Python, it is represented by:
- `cellucid.CellucidSessionBundle`

```{note}
Session bundles are treated as **untrusted input** in Python: the reader and
AnnData applier validate bounds, shapes, and exact dataset identity before
mutation.
```

---

## Two ways to capture a session

### Option A: UI download (web app workflow)

Use the viewer’s UI to save/download a `.cellucid-session`.

See web-app docs:
- {doc}`../../web_app/l_sessions_sharing/index`

### Option B (recommended in notebooks): capture from Python (no browser download)

Cellucid asks the viewer for the session and receives it back on the viewer
channel (in a notebook, the notebook's own connection; in a browser tab, the
WebSocket the viewer opened to this machine), written to a temporary file;
nothing is downloaded in the browser.

```{tip}
To bring highlight groups, label fields, scores and gene sets back onto
`adata.obs` without a session file, one call does it:
`viewer.pull(adata)` ({doc}`42_pull_everything_back`). A session also keeps
the view itself, and can be reopened in the web app.
```

---

## Step-by-step: capture a session bundle in Python

### Step 1 — Open a viewer

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    height=650,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer
```

### Step 2 — Create some meaningful state in the UI

Examples:
- create highlight groups (select cells, then “confirm” into a group)
- create user-defined categorical fields (if supported by the UI)

For the highlight mental model and UI walkthrough:
- {doc}`../../web_app/f_highlighting_selection/index`

### Step 3 — Capture the bundle from Python

Ask in one cell:

```python
request = viewer.request_session_bundle()
```

and collect it in the next:

```python
bundle = request.result()
bundle
```

A notebook widget hears the viewer between cells, so the bundle arrives after
the cell that asked has finished; `request.done` says whether it is there. (A
browser-tab viewer, `AnnDataViewer`, answers while a cell runs: there
`request.result(timeout=60)` may follow in the same cell.)

### Step 4 — Inspect what’s inside

```python
bundle.manifest.keys()
bundle.dataset_fingerprint
bundle.list_chunk_ids()
```

---

## Save the bundle as an artifact (recommended)

The session bundle `request.result()` returns is stored in a temporary location.
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
- **User-defined categorical and continuous fields** → new cell-aligned columns in `adata.obs`

Other session state (camera position, UI layout, etc.) is not applied to AnnData (it’s UI state, not analysis state).

### Apply with the convenience API

```python
adata2 = viewer.apply_session_to_anndata(adata, request, inplace=False)
adata2
```

### Apply with full control

```python
from cellucid import apply_cellucid_session_to_anndata

adata2, summary = apply_cellucid_session_to_anndata(
    bundle,
    adata,
    expected_dataset_id="my-study-v1",
    inplace=False,
    add_highlights=True,
    highlights_prefix="cellucid_highlight__",
    add_user_defined_fields=True,
    user_defined_prefix="",
    include_deleted_user_defined_fields=False,
    store_uns=True,                       # store manifest + chunk metadata under adata.uns['cellucid']
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

## Exact dataset identity (read this!)

A bundle's `datasetFingerprint` carries exactly `datasetId`, `name`,
`cellCount`, `varCount`, a `cellOrder` record, and a `dataMapping` record — the
current (format 2) session fingerprint:

```python
bundle.dataset_fingerprint
```

```text
{'datasetId': 'my-study-v1',
 'name': 'My study',
 'cellCount': 3696,
 'varCount': 3753,
 'cellOrder': {'embedding': {'key': 'X_umap', 'columns': [0, 1]},
               'dimension': 2,
               'digest': '9f2c41a0b7d3e185'},
 'dataMapping': {'kind': 'anndata',
                 'embeddings': {'enabled': [{'key': 'X_umap', 'columns': [0, 1]}],
                                'defaults': {'1d': None,
                                             '2d': {'key': 'X_umap', 'columns': [0, 1]},
                                             '3d': None}},
                 'knn': None,
                 'geneNames': 'index',
                 'expression': 'X',
                 'latent': None,
                 'vectorFields': []}}
```

`expected_dataset_id` is required. If the ID or the cell count does not match,
if the `dataMapping` names an `obsm`, `obsp`, layer, `raw`, or `var` column the
target does not have, or if `varCount` disagrees with the gene axis that mapping
reads (`adata.raw.n_vars` for `expression: "raw"`), application raises before
mutating the target. A gene-sourced user-defined field is copied from that same
axis: the gene is found among the mapping's gene names, made unique with the
shared rule, and read from the mapping's expression matrix.

### What `cellOrder` is, and why Python does not check its digest

Everything a session records about your cells is a **row number**. That is safe
only for as long as row 417 keeps meaning the same cell. Rewrite a dataset
from a re-sorted `AnnData` and the id, the cell count, and the gene count are
all still identical, while every stored row number now denotes something else.

`cellOrder` closes that hole in the viewer. It is a digest of the coordinates
of one embedding — the one `cellOrder.embedding` names, an `obsm` key with its
0-based columns — plus that embedding's dimension.
Every embedding is normalized on its own, so a digest only means anything for
the embedding it was taken over.

Python requires `cellOrder`, validates its shape strictly (`embedding` a
`{key, columns}` record that fits the mapping's kind, `dimension` exactly 1, 2,
or 3, `digest` exactly 16 lowercase hex characters), and never checks the digest
value. There is nothing here to check it against: the digest covers coordinates
the viewer normalized and drew, and you — not Cellucid — chose the row order of
the object you are applying the session to. Verifying cell-order identity is the
viewer's job.

The practical consequence: **apply a session only to the `AnnData` that produced
it, in its original cell order** (a session names cells by original position;
a Cellucid-ready file read with anndata comes back in its stored priority order:
`cellucid.restore_original_order(adata)` first, and applying to it before is
refused). In the notebook that is your responsibility: a file object subset or
re-sorted since it was read only gets a warning (its priority record no
longer describes it, so its rows are taken as they are), and an object of
your own gets none: go back to the unmodified object.

```python
# Safe: the object the viewer was opened on, untouched.
adata2 = viewer.apply_session_to_anndata(adata, request, inplace=False)

# Not safe: a re-sorted view. The counts still match, so nothing raises.
resorted = adata[adata.obs.sort_values("cell_type").index].copy()
```

With `store_uns=True` the whole fingerprint, `cellOrder` and `dataMapping`
included, is written to `adata.uns["cellucid"]["session"]["dataset_fingerprint"]`,
so the applied result carries a record of the ordering and arrays the session
was saved on.

---

## Interface reference

```{figure} ../../../_static/screenshots/highlighting_selection/highlighting-selected-page.png
:alt: The Highlighting panel with one page renamed "Left lobe" carrying 826 cells, the summary "826 cells highlighted", and one group row labelled "Lasso (826 cells) (826)".
:width: 532px

Highlighted cells are stored on a named page that can be used by analysis and session workflows.
```

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: `request.result()` says the bundle has not arrived

Likely causes:
- it was called in the same cell as `request_session_bundle()`: the viewer
  answers between cells
- the viewer has not opened its dataset yet, or its channel is not open

How to confirm (in a new cell):
```python
request.done
viewer.debug_connection()
```

Fix:
- run `request.result()` again in a later cell
- check `channel_connected` in the report

### Symptom: “Applying session adds no columns”

Likely causes:
- the session did not contain highlight groups or user-defined fields
- exact dataset identity or dimensions did not match

How to confirm:
```python
bundle.list_chunk_ids()
bundle.dataset_fingerprint
```

Fix:
- create highlight groups in the UI before capturing
- pass the exact dataset ID used to create the viewer

### Symptom: “Column already exists”

Cellucid rejects an existing target column. Choose distinct
`highlights_prefix` and `user_defined_prefix` values or remove the conflicting
column explicitly before applying.

### Symptom: Python applies a bundle the web app refuses to open

The two readers are not symmetric: only the web app compares the `cellOrder`
digest against the coordinates it loads, so a re-ordered dataset is caught there
and not here.

How to confirm:
```python
bundle.dataset_fingerprint["cellOrder"]
```

Fix: open the dataset generation the session was saved on, or re-create the
state in the viewer and save a new session file. A digest is only meaningful if
it was taken at the moment the state was captured.

### Symptom: applying raises that the data mapping does not match

The session was saved while the viewer showed an array — an embedding, a graph,
a layer, or a gene-name column — that this `AnnData` does not have. The message
names each one.

```python
bundle.dataset_fingerprint["dataMapping"]
```

Fix: apply the session to the object it was saved on.

---

## Next steps

- {doc}`05_jupyter_embedding_hooks_sessions_gallery` (runnable session notebooks)
- {doc}`23_programmatic_highlighting_and_selection_callbacks` (robust hooks patterns)
- Web-app sessions: {doc}`../../web_app/l_sessions_sharing/index`
