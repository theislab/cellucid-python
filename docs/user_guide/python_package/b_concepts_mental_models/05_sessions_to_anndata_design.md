# Sessions → AnnData (Design notes)

**Audience:** developer / power user  
**Time:** 20–30 minutes  
**Goal:** understand the protocol and guardrails behind `viewer.request_session()` and `apply_cellucid_session_to_anndata(...)`.

```{note}
This page is a “deep path” companion to {doc}`05_sessions_to_anndata_bridge`.
Most users should start there.
```

---

## What’s implemented (today)

### 1) Jupyter: “no download” session capture

In a notebook, Python can request the viewer's session (its whole history) from the frontend and receive it as a {class}`~cellucid.CellucidSession`.

High-level flow:
1) Python sends `{"type": "requestSessionBundle", "requestId": ...}` to the viewer on its channel (the notebook widget's comm, or `/_cellucid/ws` for a browser tab).
2) The web app writes the dataset's history as raw `.cellucid-session` bytes (session format 4).
3) The web app sends those bytes back on the same channel as one binary frame
   (header `{"type": "session_bundle", "requestId": ...}`, then the bundle).
4) Python streams the frame to a temp file (bounded memory) and validates it.
5) A `type="session_bundle"` event completes the request, so
   `request.result()` returns (in a notebook, in a later cell).

Why this matters:
- Users never have to find/download a file from the browser.
- The artifact is the same `.cellucid-session` format used by the web app.

### 2) Python: a `.cellucid-session` reader that checks everything first

`cellucid.read_session(path)` reads the file chunk by chunk, in order, and
checks all of it before it returns: every chunk inflated to exactly its
declared size, every array against its digest, every value against its hash,
every record against the history model, every tracker's tree whole. A damaged
file is refused whole, never read in part; a file written by a newer
Cellucid than this package is refused as newer, never as damaged. It then
answers from memory through the reading model ({mod}`cellucid.history`):
`session.tracker(name_or_id)` is a `Tracker` (its steps as a table, its
branches, turns, notes and journal), `tracker.step(id).state()` the state
after a step, `step.diff(other)` what differs between two steps, and
`session.conversations(tracker)` the conversations kept with a tracker.
cellucid-r reads the same records function for function
(`cellucid_read_session()`, `cellucid_tracker()`, `cellucid_steps()`, ...).

### 3) AnnData mutation (safe, minimal, DRY)

`apply_cellucid_session_to_anndata(...)` plans the state at one step
(default: the active tracker's head) as rows, and writes them only when asked
(`dry_run=True` is the default; `plan.apply()` or `dry_run=False` writes):
- **Labels** → renamed `.obs` columns and categories
- **Highlights** → one boolean `.obs` column per highlight group
- **User-defined categorical and continuous fields** → cell-aligned `.obs` columns
- **Colors** → each written, renamed or relabeled column's
  `uns["<column>_colors"]`, as Scanpy keeps them: a categorical field the
  colors it is shown in (in the view in focus, a category's color the person
  or an agent chose, else its `palette`, else Cellucid's own category
  colors), any other written column none; a dataset column recolored in the
  view in focus, the colors it is shown in ({ref}`apply-column-colors`)
- **Gene sets** → `adata.uns["cellucid"]["gene_sets"]`, and one boolean
  `var` column each, `cellucid_gene_set__<name>`, on the gene axis the data
  mapping reads
- **Presentation** (filters, colors, views, the camera, analyses, the
  expression values) → one row per domain the step changed, reported and
  never applied; the expression values (normalized, log1p) are recorded in
  the apply's provenance (`expression`), and `X` keeps its values

Each row's status is `write`, `unchanged`, `replace` (with `overwrite=True`,
only a column or gene set an earlier apply wrote) or `not-applied` (with the
reason). Every apply that writes appends one provenance record to
`adata.uns["cellucid"]["history_applied"]` (JSON text holding a list).

The rules are the same in cellucid-python, cellucid-r
(`cellucid_apply_session()`) and the web app (whose **Save AnnData with this
state…** writes a new `.zarr.zip` of the open dataset with the state shown
applied, {ref}`history-save-anndata`), and the three are held to one shared
expected plan (its rows and its gene-set `var` columns): the golden history
`history_golden.cellucid-session` and its `history_golden.expected.json`,
written by the web repository's `scripts/generate-history-golden.mjs` and
pinned by digest in all three test suites.

---

## Why it’s implemented this way

### “Session file first”: one artifact everywhere

Everything is built around the `.cellucid-session` file:
- the web app's **Share…** and **Load…** buttons (**State › Session**)
- Jupyter no-download capture
- Python and R application to an object
- pushing a saved history into a viewer (`viewer.load_session(path)`; its
  trackers are added beside the person's own)

Benefits:
- reproducibility: the session is a sharable file artifact, and its history
  says who did what, step by step,
- one shape: every state domain's value is stored the same way (a skeleton
  and its arrays, by content hash), so a new feature's state needs no
  encoder of its own in any reader.

### Treat session files as untrusted input

Even if the session “came from your own browser”, treat it as external input:
- validate framing (the signature, every frame length),
- bound every decompression by the size its chunk declares (no gzip “zip
  bomb” can inflate past it),
- check every array against its digest and every value against its hash,
- check every record against the history model and every tracker's tree whole,
- bounds-check cells and codes when a state is applied.

This is critical for:
- robustness (sessions should not brick restore),
- safety (avoid accidental data corruption when applying).

---

## Protocol details (what happens on the wire)

### Frontend → Python events (hooks)

Everything travels on the viewer channel: one WebSocket the viewer opens to
`/_cellucid/ws` and authenticates once with its `viewerId` and `viewerToken`
({doc}`../e_jupyter_hooks/05_architecture_viewer_channel`).

The viewer sends UI events as:
- JSON text frames (ready, hover, click, data mapping, command errors)
- binary frames for selections (the selected cell indices as little-endian
  `uint32`, strictly ascending) and session bundles

The channel belongs to one viewer, so events carry no `viewerId`.

### Python → frontend commands

Python sends commands as JSON text frames on the same channel, and a command
that names cells (a highlight, a visibility change) as one binary frame of
`uint32` cell indices. A command sent before the viewer has connected waits
and is sent, in order, when it does.

### Notebook turns and the live history

An `action` command names the notebook turn it belongs to: `turn` (its id,
`python-` and 12 hexadecimal digits) and `turnLabel`, both text inside
`with viewer.turn("..."):` and both null outside one. The viewer groups a
turn's steps and undoes them together. Every change of the viewer's history
reaches Python as a `history_changed` event (`kind`: `step`, `moved`,
`edited`, `tracker` or `loaded`), and `viewer.load_session(path)` sends a
session file as one binary `load_session` command, answered by
`session_loaded` (the tracker shown and its head) or `command_error`.

### Session bundle frames

Guardrails:
- Python registers an expected `(viewerId, requestId)` pair until the request
  is collected or canceled.
- Python accepts exactly one bundle for that pair (consumes the pending
  request); a late one is dropped.
- The frame is streamed to disk, so the temporary directory's free space is
  the only size bound; an invalid bundle makes `request.result()` raise.

---

## File format (high level)

A `.cellucid-session` file is:

1) the signature `CELLUCID_SESSION/4\n` (session format 4; a session of another format is refused naming it, any other signature as not a session)  
2) `manifestByteLength` (u32 LE)  
3) `manifest` JSON bytes (UTF‑8)  
4) repeated chunks:
   - `chunkByteLength` (u32 LE)
   - `chunkBytes` (gzip'd, or as they are where gzip saves little)

The manifest contains:
- what wrote it (the Cellucid version and the digest of its action
  catalog), the action types its steps use, and what the share kept
  (`scope`)
- dataset fingerprint (id, display name, counts, cell order, data mapping)
- the dataset source (sample, link, server, or local file name)
- `requires`: the obs columns and genes the shown state uses
- the trackers (each the history model's `Tracker`: its branches, turns,
  notes and head; its steps, journal and conversations as chunk ids) and the
  active one
- the chunks: each tracker's steps, the parts (sub-values several values
  share, stored once), every array once by its digest, the values (by
  content hash; the ones the active head shows first; none in a share of the
  steps only), each tracker's journal, and the conversations when they were
  included; long steps, parts and values are cut into several chunks

In Python, you can inspect:

```python
import cellucid

session = cellucid.read_session("my.cellucid-session")
session.dataset_fingerprint
session.requires
session.trackers
tracker = session.tracker()
tracker.steps.head()
tracker.step(tracker.head_id).state()["field-overlays"]
```

The full format is on {doc}`../g_api_reference_coverage/api/sessions`.

---

## Identity and mismatch guardrails (critical constraint)

### Cell identity is index-based (today)

State values name cells by **original position** (the row in the order the
dataset gives its cells):
- a highlight group's cells are an ascending list of original positions
- a user-defined field's codes or values are a column: one value per cell, or
  per cell of an ascending list of original positions

This is why fingerprint rejection exists.

### Fingerprint rejection checks what the chosen step needs

Python compares the session's `datasetId` with `adata.uns["cellucid"]["id"]`
when the object names one, compares `cellCount` with the target AnnData (the
subsample's cells, for a session saved on one), requires every key the
`dataMapping` names to exist on it, compares `varCount` with the gene count
of the axis that mapping reads, takes the `cellOrder` digest when the object
has that embedding, and requires every obs column and gene the chosen step's
state uses. Any mismatch raises `ValueError` before anything is written.

This prevents the worst failure mode:
> applying “correct-looking” highlights to the wrong biological cells.

### `cellOrder` and `dataMapping`: recorded by the viewer, checked on apply

The manifest fingerprint carries two records beside those three scalars, and both
are required. `cellOrder` is exactly `embedding`, `dimension` (1, 2, or 3), and a
16-character lowercase hexadecimal `digest` of that embedding's stored
coordinates, each cell named by its original position. `embedding` is `{"key", "columns"}`: an `obsm` key with the 0-based
columns that were read, one per dimension. It exists because the three scalars
all survive a rewrite of the same dataset from re-sorted input, and index-based
selections would then silently denote different cells.

`dataMapping` records which stored arrays the session was showing: `kind`
(`"anndata"`), `embeddings` (`enabled`, the embeddings that are turned on, and
`defaults`, the `1d`/`2d`/`3d` default of each dimension), `knn`, `geneNames`
(`"index"` or `"column:<name>"`), `expression` (`"X"`, `"raw"`,
`"layer:<name>"`, or null), `log1p`, `highlyVariable`, `latent`, and
`vectorFields`, plus `modality` for a dataset of several modalities (a MuData).

The two readers check them as follows:

| | web app | Python (apply) |
|---|---|---|
| `cellOrder` and `dataMapping` present | required | required |
| record shapes | validated | validated |
| `dataMapping` keys exist on the dataset | checked | checked against the `AnnData` |
| `cellOrder` digest value | compared against the stored embedding of the cells it opens | compared against the `AnnData`'s stored embedding over the session's cells, when it has that embedding |

A digest is a sum over the cells: a mismatch says the cells differ (other
cells, or another order), not which. An `AnnData` without the session's
embedding is not checked for cell order.

## Troubleshooting (protocol-level)

### Symptom: “Dropped a session bundle nobody is waiting for” (log)

Meaning:
- the bundle arrived after its request was canceled (`request.cancel()`,
  or the viewer stopped).

Fix:
- request a new bundle again from the *current* viewer instance.

### Symptom: “This file is not a Cellucid session file of format 4 …”

Meaning:
- the bytes the viewer sent on its channel for a `requestSessionBundle` are
  not a `.cellucid-session` file of the current format (or a file passed to
  `cellucid.read_session(...)` is not one). On the channel, the frame is
  dropped and `request.result()` raises `RuntimeError` ("The viewer
  did not capture the session bundle: The viewer sent an invalid session
  bundle: …") instead of waiting for its timeout.

Fix:
- make sure the viewer page and the Python package are current versions of
  the same release (reload the viewer after upgrading), then request a new
  bundle.

---

## Related docs

- User-facing workflow: {doc}`05_sessions_to_anndata_bridge`
- Web app session model: {doc}`../../web_app/l_sessions_sharing/01_session_mental_model`
