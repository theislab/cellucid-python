# Sessions → AnnData (No-Download Bridge)

**Audience:** computational users + power users (wet lab users: skim the “Fast path”)  
**Time:** 20–40 minutes  
**What you’ll learn:**
- How to capture the viewer's history into Python **without** a manual browser download
- What a session holds (trackers, branches, steps, the state at every step)
  and what of it can be applied back to AnnData (a plan first, then the change)
- The safety model (dataset identity, cell order and index-based cell identity)
- How to debug “session capture is stuck” and “apply did nothing”

**Prerequisites:**
- A viewer running in a notebook (`show_anndata(...)`)
- An `AnnData` (or `MuData`) object you want to annotate, in the same cell
  order as the dataset the viewer showed

For deeper design rationale and file-format notes, see {doc}`05_sessions_to_anndata_design`.

---

## Mental model: “live state” vs “durable state”

Cellucid exposes two complementary layers of state in notebooks:

1) **`viewer.state` (live snapshot; volatile)**  
   A small “latest event” snapshot updated from UI events (ready/selection/hover/click/…).

2) **The viewer's results (`viewer.results`, `viewer.pull(adata)`)**  
   Every highlight group, label, cluster, score and gene set the viewer holds,
   sent to Python as it is made. One call writes them onto your AnnData.

3) **Session (`.cellucid-session`; durable artifact)**  
   The dataset's history: every step taken in the viewer (by a person, an
   agent or a notebook), grouped in trackers with named branches, and the
   state every step left behind (labels, highlights, fields, views, camera,
   layout...). In notebooks you can request it as a
   {class}`~cellucid.CellucidSession`.

If you want your **results as columns**, use `viewer.pull(adata)`. If you want
to **reproduce**, **transfer** or **go back to** an interactive state, use the
**session**.

---

## Fast path (wet lab / non-technical): save your work into the data object

If you don’t want to think about file formats:

1) Open your data in Cellucid (in Jupyter).
2) Create highlights/labels in the viewer (select cells, label groups, etc.).
3) In Python, in a new cell, run:

   ```python
   summary = viewer.pull(adata)
   summary.columns
   ```

4) Use `adata` for downstream analysis (Scanpy plots, marker genes, etc.):
   each highlight group is a `cellucid_highlight__<id>` column, each label
   field a categorical column ({doc}`../g_api_reference_coverage/api/round_trip`).

This turns interactive choices (highlights/labels) into explicit columns you can store and version.

---

## Step-by-step workflow: capture → (optionally save) → apply

### Step 1: create and load a viewer

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    height=650,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

### Step 2: interact in the UI

In the viewer:
- color by a field,
- select cells,
- create highlight groups,
- create user-defined categorical fields (if you use that UI).

### Step 3: capture the session into Python (no browser download)

A notebook widget hears the viewer only between cells, so ask in one cell:

```python
request = viewer.request_session()
```

and collect it in the next:

```python
session = request.result()
tracker = session.tracker()                       # the tracker the viewer showed
tracker.steps[["id", "branch", "label", "origin"]]
```

The returned {class}`~cellucid.CellucidSession` was read and checked whole
from a `.cellucid-session` file written to a **temporary file** on disk.

### Step 4 (optional but recommended): persist the artifact

If you want to keep a copy of “what I did in the UI”:

```python
session.save("./sessions/my_state.cellucid-session")
```

### Step 5: plan, then apply to AnnData

```python
plan = session.apply(adata)      # dry run: one row per change, nothing written
plan.table
report = plan.apply()            # writes the rows into adata and records the apply
```

That applies the state the viewer showed (the active tracker's head); pass
`step="s4"` (an id from `tracker.steps`) for the state after an earlier step,
`branch=` for a branch's tip, or `tracker=` for another tracker's head.

Or let the viewer apply what the request brings -- in a notebook in a later
cell, in a browser-tab viewer (`AnnDataViewer`) at once:

```python
plan = viewer.apply_session_to_anndata(adata, request)
report = plan.apply()
```

With `dry_run=False` the viewer applies the plan at once and deletes the
temporary session file (`cleanup_file=True`).

---

## What gets applied to AnnData (today)

`apply_cellucid_session_to_anndata(...)` plans the state at one step as rows,
then writes them. Each row says what happens to it: `write` (the object lacks
it), `unchanged` (the object holds it already), `replace` (a column or gene
set an earlier apply wrote, with `overwrite=True`) or `not-applied` (with the
reason).

### 1) Labels → renamed `adata.obs` columns and categories

- An obs column renamed in the session is renamed, keeping its place.
- A category renamed in a categorical (or text) obs column is renamed; two
  categories renamed to one label become one.
- A gene's display name and a field deleted in the session are recorded as
  `not-applied` rows, never applied to data. A rename naming a column or
  category `adata` does not have is `not-applied`, and its note says which.

### 2) Highlights → `adata.obs` boolean columns

- One boolean `.obs` column per highlight group, of every page.
- The column is `cellucid_highlight__<page>__<group>`, from the page's name
  and the group's label (ASCII letters and digits, every other run one
  underscore, `_2`, `_3`... when two read alike).

This is the “most useful” mapping because it makes interactive groups explicit and scriptable:
- `adata.obs["cellucid_highlight__QC__Doublets"]` becomes a mask you can use in analysis.

### 3) User-defined fields → cell-aligned `adata.obs` columns

- Each column is named by the field's key.
- A categorical field (labels, a merge of categories, an annotation) is a
  `category` column with its categories in order; a cell with no category is
  missing.
- A continuous field a notebook pushed is `float32`, `NaN` for no value.
- A continuous copy of another field copies that field: an obs column, or a
  gene's values across cells (it does not create a gene-aligned `adata.var`
  column).

### 4) Gene sets and the apply's record under `adata.uns`

- The step's gene sets go to `adata.uns["cellucid"]["gene_sets"]`, as
  `viewer.pull(adata)` writes them, and each set written is also a boolean
  `var` column, `cellucid_gene_set__<name>`, on the genes the data mapping
  reads (`X` or a layer; a MuData: the mapping's modality), `True` on the
  set's genes
  ({ref}`apply-gene-set-var-columns`).
- Every apply that writes something appends one record to
  `adata.uns["cellucid"]["history_applied"]` (JSON text holding a list; read
  it with `cellucid.history_applied(adata)`): the session's SHA-256, the
  tracker, branch and step, when, the package and its version, the columns it
  wrote, the renames it applied and what it did not apply. With
  `keep_session=True` it also keeps the session file (base64), so the history
  can be rebuilt from the object alone.

### What is never written

- A column or gene set the person made is never replaced: its row is
  `not-applied` and says so, and the rest applies.
- `obs['cellucid_original_position']`, the priority-order column, is never
  written, renamed or replaced: a plan that would is refused.
- Applying the same step again writes nothing: every row is `unchanged` or
  `not-applied`, and `report.applied` is False.

```{important}
Not everything in a session becomes a clean AnnData mutation.
The goal is to materialize the parts that make sense as dataset annotations (labels, highlights, fields, gene sets).
Presentation (filters, colours, views, the camera, analyses) is reported as one
`presentation` row per domain the step changed and never applied; it stays in
the session (`tracker.step(step).state()` reads it).
```

---

## Safety model: dataset mismatch + index-based identity

### Sessions are treated as untrusted input

Even if the session came from “your own browser”, treat it like external input.
`cellucid.read_session(...)` reads and checks the whole file before anything
is applied:
- file framing validation (the signature, every frame length),
- every decompression bounded by the size its chunk declares,
- every array checked against its digest and every value against its hash,
- every record checked against the history model, every tracker's tree whole,
- bounds checks on cells and codes when a state is applied.

A file written by a newer Cellucid than this package is refused as newer:
*This session was written by a newer Cellucid (X) than this package reads
(Y); update the cellucid package to read it.*

### Dataset match checks (only what the chosen step needs)

Before planning anything, Python checks the object against the session:
- `datasetId` vs `adata.uns["cellucid"]["id"]`, when the object names one
- `cellCount` vs `adata.n_obs` (the subsample's cells, for a session saved on
  one)
- every key the `dataMapping` names vs what `adata` holds
- `varCount` vs the gene count of the mapped axis (`adata.raw.n_vars` when the
  mapping reads expression from `raw`, `0` when it reads none, `adata.n_vars`
  otherwise)
- the `cellOrder` digest, when `adata` has the embedding it names
- every obs column and gene the chosen step's state uses vs what `adata` holds

Any mismatch raises `ValueError` before Cellucid changes the target object.

### The two records: `cellOrder` and `dataMapping`

The current session fingerprint (format 4) has exactly six required fields:
the three scalars (`datasetId`, `cellCount`, `varCount`), the dataset's display
`name` (or null, informational), `cellOrder`, and `dataMapping`.

`cellOrder` is exactly `embedding`, `dimension` (1, 2, or 3), and a
16-character lowercase hexadecimal `digest` of that embedding's stored
coordinates, each cell named by its original position. `embedding` is `{"key", "columns"}`: an `obsm` key with the 0-based
columns that were read, one per dimension. It ties the session to one cell
ordering, which the scalars cannot do: they all survive a rewrite of the
same dataset from re-sorted input.

`dataMapping` records which stored arrays the session was showing: `kind`
(`"anndata"`), `embeddings` (`enabled`, the embeddings that are turned on, and
`defaults`, the `1d`/`2d`/`3d` default of each dimension), `knn`, `geneNames`
(`"index"` or `"column:<name>"`), `expression` (`"X"`, `"raw"`,
`"layer:<name>"`, or null), `log1p`, `highlyVariable`, `latent`, and
`vectorFields`, plus `modality` for a MuData. A user-defined field copied from a gene is resolved by name
on that axis — through the mapping's gene names, made unique with the shared
rule — and its values are read from the mapping's expression matrix.

When `adata` has the embedding `cellOrder` names, the apply takes the same
digest of its stored values over the session's cells and refuses the object
when it differs: the rows are other cells, or in another order. A digest is a
sum over the cells, so a mismatch says the cells differ, not which.

### Critical constraint: cells are named by original position

Highlights and user-defined fields name each cell by its **original
position**: its row in the order the dataset's cells were given, whatever order
a file stores them in (a Cellucid-ready file stores them in priority order) or
the viewer opened them in. A session saved on a subsample (or a step taken
while one was open) applies to the cells that were open: the others are
`<NA>` in a highlight column and missing in a field's.

Apply a session to an `AnnData` in that original order: the object you
prepared or served, or a Cellucid-ready file read back through
{func}`~cellucid.restore_original_order` (an `AnnData` still in a stored
priority order -- its `uns['cellucid']['priority']` describes its cells -- is
refused with that advice; one whose record no longer describes its cells,
because it was subset, combined or reordered since it was read, is taken row
by row as it is, with a warning, as the web app reads such a file). The
cell-order digest catches other cells or another order only when `adata` has
the session's embedding.

Dataset ids are compared when the object names one in
`uns['cellucid']['id']` (see {doc}`04_dataset_identity_and_reproducibility`).

---

## Advanced usage

### Apply from a saved `.cellucid-session` file (web app or notebook)

```python
import cellucid

session = cellucid.read_session("./sessions/my_state.cellucid-session")
plan = session.apply(adata)
report = plan.apply()
```

### Apply another step, branch or tracker

```python
tracker = session.tracker("Annotation")
tracker.branch("Without score").steps[["id", "label"]]

plan = session.apply(adata, tracker="Annotation", branch="Without score")
plan = session.apply(adata, step="s5")              # the state after step s5
```

### Replace what an earlier apply wrote

A column or gene set an earlier apply wrote, holding other values now, is
`not-applied` until `overwrite=True` makes it `replace`:

```python
plan = session.apply(adata, step="s13", overwrite=True)
plan.table[plan.table["status"] == "replace"]
```

A column the person made is never replaced, whatever `overwrite` says.

### Backed AnnData targets

Applying changes the object's `obs` and `uns` in memory. A backed AnnData
keeps `X` on disk, and its source file is not rewritten: write the object to
keep the result.

---

## Troubleshooting (large, symptom-based)

### Symptom: `request.result()` raises at once in a notebook

It says a notebook widget hears the viewer only between cells: collect
`request.result()` in the next cell.

### Symptom: `request.result()` says the session has not arrived

Likely causes (ordered):
1) The viewer never reached “ready” (dataset load stuck).
2) The viewer is not shown in any output (no frontend to respond).
3) The viewer's channel is not connected.

How to confirm:
- Run `viewer.debug_connection()` and inspect `state.ready`,
  `connected_outputs` and `channel_connected`.

Fix:
- Show the viewer (end a cell with `viewer`), wait until it shows the data,
  then ask again.

### Symptom: “Dropped a session bundle nobody is waiting for” (log)

This means the session arrived after its request was cancelled
(`request.cancel()`, or the viewer stopped). A viewer tab from an old kernel cannot send one at all:
its channel is refused and the view freezes.

Fix:
- request the session again;
- after a kernel restart, re-run the cell that shows the viewer.

### Symptom: the apply raises a dataset mismatch

Likely causes:
- `adata.n_obs` / `adata.n_vars` differ from the session’s fingerprint,
- `adata.uns["cellucid"]["id"]` names another dataset than the session's,
- the AnnData was subset/reordered since the session was created (the
  cell-order digest differs).

How to confirm:
- Inspect `session.dataset_fingerprint`.

Fix:
- apply to the matching AnnData version, in its original order (a
  Cellucid-ready file read back: `cellucid.restore_original_order(adata)`).

### Symptom: “My highlight groups didn’t become columns”

Likely causes:
- the state at the step applied holds no highlight group (they were made
  after it, on another branch, or on another tracker),
- a column of that name exists already and Cellucid did not write it (its row
  is `not-applied`),
- dataset identity or shape validation failed.

How to confirm:
- `plan.table` lists every highlight row and its status;
- `session.tracker(tracker).step(step).state()["highlights"]["pages"]` lists the groups.

Fix:
- apply the step (`step=`), branch (`branch=`) or tracker (`tracker=`) that holds them,
- ensure dataset identity matches.

### Symptom: “Applying user-defined fields fails or gives weird categories”

Likely causes:
- a categorical field's codes name a category it does not list,
- a continuous field names a source column or gene `adata` does not hold.

Each of these raises `ValueError` naming the field before anything is changed.

Fix:
- apply to the dataset the session was saved on,
- save the session again from the viewer.

## Next steps

- Deep design and file format: {doc}`05_sessions_to_anndata_design`
- Dataset identity and versioning: {doc}`04_dataset_identity_and_reproducibility`
