# Session persistence and restoring analysis artifacts

This tutorial explains Cellucid **sessions** -- the viewer's whole history in
one file -- and how to use them as a bridge:

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
- Minimal “capture a session”: ~10–15 minutes
- Full “apply to AnnData and analyze”: ~30–60 minutes

**Prerequisites**
- A notebook environment
- A viewer (`show_anndata(...)`)

---

## What a session is

A session is a single file with extension:
- `.cellucid-session` (session format 4)

It holds the dataset's **history**:
- a manifest (what wrote it, `datasetFingerprint`, `source`, `requires` — the
  obs columns and genes the shown state uses —, the subsample it was saved
  on, what the share kept, and the trackers)
- each tracker's steps: which action each applied, with which parameters, who
  took it (the person, an agent, a notebook) and its label, in a tree with
  named branches, with notes, bookmarks, turns and the journal of the
  history's own edits
- the values of the viewer's state the steps show (labels, user-defined
  fields, highlights, gene sets, filters, colours, views, camera, layout),
  each stored once

In Python, it is read with {func}`cellucid.read_session` as a
{class}`cellucid.CellucidSession`; its trackers, branches and steps are the
reading model of {mod}`cellucid.history`.

```{note}
Sessions are treated as **untrusted input** in Python: `read_session` reads
the whole file and checks every chunk, array, value and record before
anything is returned, and the AnnData applier checks the dataset before it
plans anything. A session file of format 3, saved by an earlier Cellucid, is
refused with the web app's sentence; one written by a newer Cellucid than
this package is refused as newer: *This session was written by a newer
Cellucid (X) than this package reads (Y); update the cellucid package to read
it.*
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

## Step-by-step: capture a session in Python

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

### Step 3 — Capture the session from Python

Ask in one cell:

```python
request = viewer.request_session()
```

and collect it in the next:

```python
session = request.result()
session
```

A notebook widget hears the viewer between cells, so the session arrives after
the cell that asked has finished; `request.done` says whether it is there. (A
browser-tab viewer, `AnnDataViewer`, answers while a cell runs: there
`request.result(timeout=60)` may follow in the same cell.)

### Step 4 — Inspect what’s inside

```python
session.dataset_fingerprint
tracker = session.tracker()              # the tracker the viewer shows
tracker.steps[["id", "branch", "label", "origin"]]
tracker.branches                         # each branch: its fork step, its tip
tracker.step(tracker.head_id).state()["highlights"]   # the highlight pages at the step shown
```

---

## Save the session as an artifact (recommended)

The session `request.result()` returns is stored in a temporary location.
Save it to a stable path:

```python
saved_path = session.save("sessions/my_analysis.cellucid-session")
saved_path
```

---

## Apply a session to an AnnData (the “session → analysis” bridge)

### What gets applied (important)

The Python applier plans the state of one step -- by default the step the
viewer shows -- as one row per change, and writes the rows that are
meaningful in analysis contexts:

- **Labels** → renamed `obs` columns and categories
- **Highlights** → new boolean `adata.obs` columns
  - one column per highlight group, on every page,
    `cellucid_highlight__<page>__<group>`
- **User-defined categorical and continuous fields** → new cell-aligned
  columns in `adata.obs`, each named by its key
- **Gene sets** → `adata.uns['cellucid']['gene_sets']`
- **A record of the apply** → appended to
  `adata.uns['cellucid']['history_applied']` (JSON text; read it with
  `cellucid.history_applied(adata)`)

Presentation (filters, colours, views, the camera, analyses) is reported as
one row per domain the step changed and never applied: it is how the data was
shown, not the data. Any other step: `step="s3"`, a branch's tip with
`branch=`, another tracker with `tracker=`; `tracker.steps` lists the ids.

### Plan, then apply

```python
plan = viewer.apply_session_to_anndata(adata, request)   # dry run: nothing written
plan.table
```

Each row has `kind`, `domain`, `name`, `column`, `from`, `cells`, `status`
and `note`. The status is `write` (the object lacks it), `unchanged` (the
object holds it already), `replace` (a column or gene set an earlier apply
wrote, with `overwrite=True`) or `not-applied` (the note says why: a column
the person made is never replaced; a column or category the object does not
have; a deletion; a presentation domain).

```python
report = plan.apply()          # writes every row at once into adata
report.applied                 # False when there was nothing to write
cellucid.history_applied(adata)[-1]
```

Applying the same step again writes nothing: every row is `unchanged` or
`not-applied`.

### Apply with full control

```python
import cellucid

plan = cellucid.apply_cellucid_session_to_anndata(
    adata,
    session,
    tracker=None,                 # the tracker the viewer shows
    branch=None,                  # its active branch
    step=None,                    # the step the viewer shows
    dry_run=True,                 # plan only; False applies at once
    overwrite=False,              # replace what an earlier apply wrote with other values
    keep_session=False,           # True keeps the session file (base64) in the apply's record
)
report = plan.apply()
```

### How to use the result in analysis

List the columns written:

```python
written = report.provenance.columns
written
```

Subset by a highlight group:

```python
highlight_col = next(name for name in written if name.startswith("cellucid_highlight__"))
selected = adata[adata.obs[highlight_col].fillna(False).astype(bool)].copy()
selected
```

---

## Dataset identity (read this!)

A session's `datasetFingerprint` carries exactly `datasetId`, `name`,
`cellCount`, `varCount`, a `cellOrder` record, and a `dataMapping` record:

```python
session.dataset_fingerprint
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

Only what the chosen step needs is checked, before anything is planned: the
dataset id when `adata.uns['cellucid']['id']` names one, the cell count, the
keys the `dataMapping` names (an `obsm`, `obsp`, layer, `raw`, or `var`
column), the gene count of the axis that mapping reads (`adata.raw.n_vars`
for `expression: "raw"`), the cell order when `adata` has the embedding
`cellOrder` names, and the obs columns and genes the step's state uses. A
mismatch raises before the target is changed. A gene-sourced user-defined
field is copied from that same axis: the gene is found among the mapping's
gene names, made unique with the shared rule, and read from the mapping's
expression matrix.

### What `cellOrder` is, and when it is checked

Everything a session records about your cells is a **position**: the cell's row
in the dataset's original order. That is safe only for as long as position 417
keeps meaning the same cell. Rewrite a dataset
from a re-sorted `AnnData` and the id, the cell count, and the gene count are
all still identical, while every stored position now denotes something else.

`cellOrder` closes that hole. It is a digest of the stored coordinates of one
embedding — the one `cellOrder.embedding` names, an `obsm` key with its
0-based columns — each cell named by its original position, plus that
embedding's dimension. A digest only means anything for the embedding it was
taken over; how the viewer opened the cells (all of them, a share, grown step
by step) does not change it.

The viewer compares it on the cells it opens. Applying a session compares it
when `adata` has that embedding: the digest of its stored values over the
session's cells must be the session's, or the apply is refused (a digest is a
sum over the cells, so a mismatch says the cells differ, not which). An
object without that embedding is not checked for cell order.

The practical consequence: **apply a session only to the `AnnData` that produced
it, in its original cell order** (a session names cells by original position;
a Cellucid-ready file read with anndata comes back in its stored priority order:
`cellucid.restore_original_order(adata)` first, and applying to it before is
refused). A file object subset or re-sorted since it was read gets a warning
(its priority record no longer describes it, so its rows are taken as they
are); the cell-order digest refuses it when it keeps the session's embedding.

```python
# Safe: the object the viewer was opened on, untouched.
plan = viewer.apply_session_to_anndata(adata, request)

# Refused: a re-sorted copy that keeps X_umap gives another cell-order digest.
resorted = adata[adata.obs.sort_values("cell_type").index].copy()
```

`report.provenance.session_sha256` names the session file an apply came from;
with `keep_session=True` the record keeps the file itself, its fingerprint
included.

---

## Interface reference

```{figure} ../../../_static/screenshots/highlighting_selection/highlighting-selected-page.png
:alt: The Highlighting panel with one page renamed "Left lobe" carrying 826 cells, the summary "826 cells highlighted", and one group row labelled "Lasso (826 cells) (826)".
:width: 532px

Highlighted cells are stored on a named page that can be used by analysis and session workflows.
```

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: `request.result()` says the session has not arrived

Likely causes:
- it was called in the same cell as `request_session()`: the viewer
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
- the step applied holds no highlight groups or user-defined fields (an
  earlier step, another branch, or another tracker may)
- a column of that name exists already and Cellucid did not write it (its
  row is `not-applied`)
- the same step was applied before (every row is `unchanged`)

How to confirm:
```python
plan.table[["kind", "name", "status", "note"]]
tracker.steps[["id", "branch", "label"]]
tracker.step(tracker.head_id).state()["field-overlays"]["userFields"]
```

Fix:
- create highlight groups in the UI before capturing
- apply the step, branch or tracker that holds them

### Symptom: “This object's cells are not the cells the session was saved on, in the same order”

The cell-order digest of the session's embedding differs: `adata`'s rows are
other cells, or in another order.

How to confirm:
```python
session.dataset_fingerprint["cellOrder"]
```

Fix: apply the session to the object it was saved on, in its original cell
order. A digest is only meaningful if it was taken at the moment the state
was captured.

### Symptom: applying raises that the data mapping does not match

The session was saved while the viewer showed an array — an embedding, a graph,
a layer, or a gene-name column — that this `AnnData` does not have. The message
names each one.

```python
session.dataset_fingerprint["dataMapping"]
```

Fix: apply the session to the object it was saved on.

---

## Next steps

- {doc}`05_jupyter_embedding_hooks_sessions_gallery` (runnable session notebooks)
- {doc}`23_programmatic_highlighting_and_selection_callbacks` (robust hooks patterns)
- Web-app sessions: {doc}`../../web_app/l_sessions_sharing/index`
