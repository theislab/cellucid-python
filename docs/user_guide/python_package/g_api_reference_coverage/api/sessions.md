# Sessions (`.cellucid-session` files)

```{eval-rst}
.. currentmodule:: cellucid
```

A Cellucid session file holds *the history of what was done in the viewer*:
every step a person, an agent or a notebook took on a dataset (renaming a
category, highlighting cells, merging categories, pushing a column, saving a
gene set...), grouped in **trackers** (independent lines of work, each a tree
of steps with named branches), with notes, bookmarks, turns, the journal of
the history's own edits, and the state every step left behind. Python reads
it, and applies the state at any step to an AnnData: a plan first, then the
change.

- File extension: **`.cellucid-session`** (session format 4)
- Use cases:
  - “I selected and annotated cells in the UI; now I want those labels back in Python.”
  - “I want the labels as they were *before* that last merge, or on another branch.”
  - “I want a reproducible trail of what was done, and by whom.”

This page documents:
- {func}`~cellucid.read_session` and {class}`~cellucid.CellucidSession` (read a session: its trackers, branches, steps, notes and turns)
- {mod}`cellucid.history` (the reading model: `Tracker`, `Branch`, `Step`, `Diff`)
- {func}`~cellucid.apply_cellucid_session_to_anndata`, {class}`~cellucid.ApplyPlan`, {class}`~cellucid.ApplyReport` and {func}`~cellucid.history_applied` (plan and apply the state at one step)

---

## Audience + prerequisites

**Audience**
- Wet lab / beginner: use the “Apply to AnnData” workflow and the troubleshooting section.
- Computational: pay attention to exact dataset identity and cell ordering assumptions.
- Developer: read the “Format + safety” section before handling session files from untrusted sources.

**Prerequisites**
- If applying to AnnData: an `AnnData` (or `MuData`) object. Pandas is part
  of the Cellucid core installation.
- A session file, saved from the web app (**Share…**) or asked of a
  notebook viewer.

---

## Fast path (apply what was done in the viewer to an AnnData)

```python
import cellucid

session = cellucid.read_session("my-session.cellucid-session")

plan = session.apply(adata)      # what applying the step the viewer showed would change
plan.table                       # one row per change: kind, domain, name, column, from, cells, status, note
report = plan.apply()            # writes the rows into adata, all at once, and records the apply
```

Now inspect:
- `adata.obs`: renamed columns and categories, one boolean column per
  highlight group, one column per user-defined field.
- `adata.uns["cellucid"]["gene_sets"]`: the gene sets of the step.
- `cellucid.history_applied(adata)`: every apply recorded on `adata`.

---

## Practical path (common workflows)

### 1) How to obtain a session file

#### A) From the web app (download a `.cellucid-session` file)

Use this when you explore in the standalone browser viewer (notebook or not):
**Share…** writes the dataset's history (the trackers and branches ticked),
the current state only, or the steps only to one file. A file of the steps
only holds no values: it can be read, not applied.

```{figure} ../../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The State section with the Pancreatic endocrinogenesis sample open: Session with Share… and Load…, and the History with the recording switch, the branch shown and the graph of the steps taken so far, the step shown marked.
:width: 540px

The **State** section of the sidebar keeps **Share…** and **Load…** under **Session**, above the history they write and read.
```

#### B) From a notebook viewer without downloading

A notebook or browser-tab viewer hands Python the session directly
(`viewer.request_session()`, then `result()`; {doc}`jupyter`). To bring
results onto an AnnData without a session, `viewer.pull(adata)`
({doc}`round_trip`).

### 2) Read the history

```python
session = cellucid.read_session("my-session.cellucid-session")

session.dataset_fingerprint   # dataset id, name, counts, cell order, data mapping
session.source                # where it was opened from, or None
session.requires              # Requires(obs=(...), var=(...)): the fields the shown state uses
session.trackers              # every tracker, in the file's order

tracker = session.tracker()                   # the tracker the viewer showed
tracker = session.tracker("Annotation")       # or one by its name or id
tracker.steps[["id", "branch", "label", "origin"]]
tracker.branch("Without score").steps         # a branch: its path from the root to its tip
tracker.notes_table()                         # every note: on the tracker, a branch, a group, a step or a turn
tracker.group("QC").steps                     # a group: a named run of consecutive steps of one branch
```

A tracker (`cellucid.history.Tracker`) holds every field of the history
model's `Tracker` in snake case (`head_id`, `active_branch`, `main_branch`,
`paused`, `start`, `removed`, `branches`, `groups`, `step_records`, `turns`,
`notes`, `journal`, `chats`). `paused` is the tracker's recording switch:
`None` while it records, else a `cellucid.history.Pause` with `at` (since
when), `by` (who paused it, then everyone who acted while it was paused) and
`draft` (always `None` in a file: the web app squashes what differs into one
step before it saves); the changes made while paused became one squashed step
when it resumed. `tracker.steps` is
every step as a DataFrame, one row per step in step order (parents before
children), with the columns `cellucid.history.STEP_COLUMNS`: `id`, `branch`,
`parent`, `kind` (`root`, `action`, `outside` -- a change made outside an
action -- or `squashed`), `type` (the action), `v` (the version of the
action's parameters it was recorded with), `label` (what it did, in a
sentence), `code` (its line of pseudocode, as recorded), `origin` (`user`,
`agent`, `host` -- a Python or R notebook --, `session`, `replay`, `startup`
or `outside`), `agent`, `model`, `host`, `turn`, `at` (a UTC timestamp),
`status` (`applied` or `unapplicable`), `stale` (its result no longer matches
its inputs at the head), `bookmark`, `group` (the group it is in) and
`n_notes`.

A **squashed** step has no action of its own (`type`, `params`, `v` and
`code` are `None`): it is the net change of several steps merged into one
("Merge into one step", a squashed group or turn, the history cleared keeping
the state) or of what changed while recording was paused, and
`step.squashed` says why (`cause`: `steps`, `clear`, `paused` or `share`),
how many steps it replaced (`count`, `None` when none was recorded), who
acted in it (`by`) and what the replaced steps carried that it keeps
(`carried`, a `cellucid.history.Carried`: the ids of their notes, now its own,
and their bookmarks; empty unless the cause is `steps`). A squash in the
middle of a branch keeps the id of the last step it replaced, so the steps
and branches after it keep their parent. Its state is exact, as any step's. A step's `code` is its
pseudocode as it was recorded: parameters of an earlier version (`v` below
the action's version in `cellucid.action_catalog.json`) are shown and applied
as they were; the web app upgrades them only to replay them.

| Call | What it returns |
|---|---|
| `tracker.step("s5")` | one step (`cellucid.history.Step`): `by` (who took it), `params`, `outcome`, `touches`, `notes`, `bookmark`, ... |
| `tracker.branch(name_or_id)` | one branch (`cellucid.history.Branch`): `main`, `path`, `steps` (its path from the root), `own_steps` (the steps made on it); `None` is the active branch |
| `tracker.group(name_or_id)` | one group (`cellucid.history.Group`): `name`, `notes`, `step_ids`, `first`, `last`, `branch`, `steps` (its steps as a table); each step names its group (`step.group`) |
| `tracker.turn(turn_id)` | one turn: a notebook's or an agent's group of steps (`tracker.turns` lists them; `viewer.turn()` makes one) |
| `tracker.path_to("s8")` | the steps from the root to a step, root first |
| `tracker.stale()` | the results stale at the head: `{"step", "cause", "domains"}` each |
| `tracker.methods(branch)` / `tracker.pseudocode(branch)` | a branch (the main one when left out) as methods text, or as pseudocode, as the web app writes them: the tracker's first note opens the methods; groups are headings in the methods and `# --- QC ---` … `# --- end of QC ---` in the pseudocode; a bookmark is `# Bookmark: <name>` (one a squashed step kept from a merged step adds "(on a step merged into this one)"); a squashed step is said as what it stands for (changes made with recording paused, steps merged, a fresh start, a shared state) and holds only values |
| `session.conversations(tracker)` | the agent conversations kept with a tracker, when the person chose to include them |

A name is compared as the web app compares names (case and extra white space
aside). An id or name the session does not hold is refused, naming the ones
it does.

### 3) The state at any step

```python
step = tracker.step("s5")
state = step.state()                            # every state domain's value after step s5

state["field-overlays"]["renamedCategories"]    # [{"source", "field", "category", "name"}]
state["highlights"]["pages"]                    # pages of highlight groups, cells by original position
state["gene-sets"]["sets"]                      # [{"name", "description", "genes"}]

step.diff("s13").table                          # what differs between two steps, domain by domain
```

`step.state()` returns every state domain's value as the web app captures it
(`cells`, `mapping`, `field-overlays`, `gene-sets`, `views`, `filters`,
`colors`, `highlights`, `camera`...). Cells are **original positions** (rows
of `obs` in the order the dataset gives them), `Uint32Array` cells come back
as `uint32` arrays and `Float32Array` values as `float32` (`NaN`: no value).
A session shared as its steps only (`session.scope.content == "steps"`)
holds no values, and its states are refused. A session shared as the current
state only (`"state"`) holds one tracker: the dataset as opened and one
squashed step, "Shared state"; it reads and applies as any session. No
picture of a step is ever stored: the web app draws one from the live view
when it is asked for.

### 4) Apply to AnnData (a plan, then the change)

```python
plan = cellucid.apply_cellucid_session_to_anndata(
    adata,
    "my-session.cellucid-session",
    tracker="Annotation",         # default: the tracker the viewer showed
    branch="Without score",       # default: the tracker's head on its active branch
)
plan.table
plan.writes_anything             # whether a row would write or replace
report = plan.apply()
report.applied, report.provenance
```

`session.apply(adata, ...)` is the same call with that session. With
`dry_run=True` (the default) nothing is written: the {class}`~cellucid.ApplyPlan`
lists every change, one row each. `plan.apply()`, or `dry_run=False`, writes
it into `adata` in place and returns the {class}`~cellucid.ApplyReport`.
`step="s5"` applies the state after that step (on `branch`, when both are
given).

The rows, in domain order:

| `kind` | What it is |
|---|---|
| `column-renamed` | an obs column renamed in place; a gene's display name is recorded, never applied |
| `categories-renamed` / `categories-merged` | a category renamed, or renamed into one the column has: the two become one |
| `field-deleted` | a deletion, recorded and never applied: the column stays |
| `user-field` | one column per field the session made, named by its key: a categorical as `category` (no category: missing), a column a notebook pushed as `float32` (no value: `NaN`), a continuous copy of another field copied from that obs column or gene (under the session's data mapping); a field deleted in the session is not applied |
| `gene-set` | each gene set, written to `adata.uns["cellucid"]["gene_sets"]` as `cellucid.pull` writes them, and as one `bool` column `cellucid_gene_set__<name>` of the `var` of the gene axis the data mapping reads, `True` on the set's genes ({ref}`apply-gene-set-var-columns`) |
| `highlight` | one boolean column per group of every page, `cellucid_highlight__<page>__<group>` from the page's name and the group's label |
| `column-colors` | a dataset column whose categories were recolored in the view in focus: its `uns['<column>_colors']` written with the colors it is shown in ({ref}`apply-column-colors`) |
| `presentation` | one row per other domain the step changed (filters, colors, views, the camera, analyses, the expression values): how the data was shown, reported and never applied |

Each row's `status`:
- `write`: the object lacks it;
- `unchanged`: the object holds it already. Applying the same step again
  writes nothing (`report.applied` is False) and is never a failure;
- `replace`: a column or gene set an earlier apply wrote, holding other
  values, with `overwrite=True`;
- `not-applied`, with the reason in `note`: a column or gene set the person
  made is never replaced, whatever `overwrite` says; a column or category the
  object does not have; a deletion; a presentation domain.

Every apply that writes something appends one record to
`adata.uns["cellucid"]["history_applied"]`, JSON text holding a list (JSON
text keeps `.h5ad` and Zarr files writable). Read it with
`cellucid.history_applied(adata)`: each record names the session's SHA-256,
the tracker, branch and step, `applied_at`, the package and its version, the
columns it wrote (which a later apply replaces only with `overwrite=True`),
the renames it applied, what it did not apply, the expression values the
genes were shown in at the step (`expression`, `{"normalize", "log1p"}`, None
when the session holds none: recorded, never applied, so `X` keeps its
values), and, with `keep_session=True`, the session file itself (base64), so the history can be
rebuilt from the object alone.

Only what the chosen step needs is checked, before anything is written: the
cell count (of the subsample the session was saved on, when it was), the
dataset id (the one `adata` was shown under, `dataset_id=`, when given: a
viewer given another id than the one the data declares saves its session under
it; else `adata.uns["cellucid"]["id"]` when it names one), the data mapping's
keys and gene count, the cell order (when `adata` has the embedding the
session names, the cell-order digest of its stored values over the session's
cells must be the session's), and the obs columns and genes the step's state
uses (every missing one named at once). A plan that would write, rename or
replace `obs['cellucid_original_position']` (the priority-order column) is
refused, dry runs included. Everything is then written at once,
so a failed apply leaves `adata` as it was. `plan.apply()` makes the plan
again first and refuses an object changed since.

The target may also be a MuData: the session is checked against what the
viewer showed (the cells, and the modality its data mapping names), labels
and columns apply to the global `mdata.obs`, and `mdata.uns` records the
apply.

(apply-column-colors)=
#### Colors go with the columns

Colors go with the columns, as Scanpy keeps them (`uns['<column>_colors']`,
one color a category in the column's order, which `sc.pl.*` and the viewer
open the column in). The same rule in Python, R and the web app's "Save
AnnData with this state…":

| What the apply does to a column | Its `uns['<column>_colors']` |
|---|---|
| writes or replaces a categorical field | the colors it is shown in, one `#rrggbb` a category, every category, used or not: in the view in focus at the step, the color the person or an agent gave a category (the session's `colors` value), else the color it opens with — its `palette` (a field made from highlight pages, in Create categorical: its pages' colors), else Cellucid's own category colors by position |
| writes or replaces any other column (a pushed or copied continuous field, a highlight group) | none: colors held under its name described another column, and are removed |
| renames it | carried along to the new name (a list of text; colors of another form are removed); stale colors under the new name are replaced |
| relabels a categorical in place, no categories merged | kept as they are |
| merges categories of a categorical | each category keeps the color of the first one it was made from (it takes that one's place); colors that did not cover the column's categories (fewer than its categories, or not text) are removed |
| relabels a text column | removed (Scanpy colors a text column by its sorted labels, which a rename reorders) |
| leaves it as it is (`unchanged`, `not-applied`) | kept as they are, unless the session recolored it (below) |
| a column that came with the dataset, recolored in the view in focus at the step (one `column-colors` row each) | the colors it is shown in, one `#rrggbb` a category in the column's order once the apply's renames are made: the color the person or an agent gave a category, else the color it opens with — the column's own `uns` colors as the viewer draws them (`#rgb` widened, the alpha of `#rrggbbaa` dropped), when there is a hex color for every category, else Cellucid's own category colors by position. `write` when that differs from what the column holds, `unchanged` when it holds them already; `not-applied` for a column that is not categorical, a column of a MuData modality, or a column with none of the categories recolored (a recolored category the column lacks is named in the note) |

The plan's `colors` row (when the step's colors differ from the tracker's
opening) stays `not-applied`: every other color the step holds — a color map
or range, a gene's colors, another view's colors — is how the data was shown.
Its note says so when a field's or a dataset column's colors went with its
column.

A dry run writes nothing. For a MuData the colors go to `mdata.uns`, beside
the global `mdata.obs` columns they color.

---

## Deep path (format + safety)

### Treat session files as untrusted input

Session files may come from collaborators, community annotation workflows,
or downloads of unknown origin. {func}`~cellucid.read_session` reads and
checks the whole file, as the web app does before it applies anything:
- the signature and every frame length against the bytes the file holds;
- every chunk inflated to exactly its declared size (zip-bomb protection);
- every array against its digest, every value against its hash;
- every record against the history model (`history_model.schema.json`, byte
  for byte the web app's and cellucid-r's);
- every tracker's tree whole: one root, branches that are lines from their
  fork step to their tip, groups that are runs of consecutive steps holding
  whole turns, unique names, the head on the active branch's path;
- what the share left out is absent (`scope`: its `content`, notes,
  conversations, Recently deleted, times); a value or array no step refers
  to, one a step needs that the file lacks, and trailing bytes are refused.

Each refusal is one sentence starting *This session file is damaged:* and
saying what. A file written by a newer Cellucid than this package is refused
as newer, never as damaged: *This session was written by a newer Cellucid
(X) than this package reads (Y); update the cellucid package to read it.*

### Format overview (for developers)

The `.cellucid-session` framing (session format 4) is:
1. the signature `CELLUCID_SESSION/4\n` -- the format version is part of it
2. `manifestByteLength` (u32 little-endian)
3. UTF-8 manifest JSON bytes
4. repeated chunks: `[chunkByteLength (u32 LE), chunkBytes...]`

The manifest holds `createdAt`, `app` (the Cellucid version that wrote it and
the SHA-256 of its action catalog), `actions` (one stamp per action type
and parameters version its steps use, `{type, version, params}`, sorted by
type then version, `params` the digest of that version's parameter schema), `scope` (what the share
chose to include: the trackers and branches, `content` -- `history`, `state`
or `steps` --, `results` -- `include` or `recompute` --, notes,
conversations, Recently deleted, times, the dataset's source),
`datasetFingerprint` (`datasetId`, `name`, `cellCount`,
`varCount`, `cellOrder`, `dataMapping`), `source` (where the dataset was
opened from: a sample, a link, a server, or the name of a local file or
folder; `null` for a notebook, or a share with the dataset's name only),
`requires` (the obs columns and genes the shown state uses), `chrome` (page
chrome kept outside the history, which the web app puts back after a load:
`{"pressed": {<group id>: <button id>}}`, the pressed button of each
pressed-button group, such as the selection tool; `{"pressed": {}}` from a
file written without a page; read-only as `session.chrome`, never applied), `cells` (a subsample's), `genes` (a gene
sample's), `activeTracker`, `trackers` and `chunks`.

| Chunk (id) | What it holds |
|---|---|
| `steps` (`steps/<n>/<k>`, eager) | a run of one tracker's steps as JSON, in id order; a tracker lists its steps chunks, read in order as one list; `before`, `after` and `reads` name state values by their SHA-256 hash, `params` and `outcome` are skeletons |
| `parts` (`parts/<priority>/<k>`) | `{<name>: <skeleton>}`: sub-values stored once (below); eager for what the shown values name, else lazy |
| `array` (`array/<digest>`) | one typed array a value or step holds, stored once by its digest `<type>:<length>:<sha256 of its little-endian bytes>`: `raw` elements, or `cell-set` (a uvarint bound, then the cell-set codec's smallest encoding); eager for what the shown values name, else lazy |
| `values` (`values/<priority>/<k>`) | `{<hash>: <skeleton>}`: eager, what the active head shows; lazy, every other value a step or a journal entry refers to; at least one eager chunk, none in a share of the steps only |
| `journal` (`journal/<n>`, lazy) | one tracker's edits of the history itself (a note, a bookmark, a group, a squash, a branch locked, a step deleted) |
| `chat` (`chat/<n>/<m>`, lazy) | one conversation kept with a tracker, when the person chose to include them |

Every eager chunk comes before any lazy one, and within a priority the kinds
come in the order of the table. `steps`, `parts` and `values` are JSON cut
into pieces of about 16 MiB of text. Each chunk is stored gzip'd (`codec:
"gzip"`) where that saves at least a sixteenth of its bytes, else as it is
(`codec: "none"`, its `storedBytes` equal to its `uncompressedBytes`).

A skeleton is a value as JSON with each typed array as `{"$array": <digest>}`,
each number JSON cannot hold as `{"$number": "NaN" | "Infinity" |
"-Infinity" | "-0"}`, and each object or list whose skeleton text is at least
512 characters as `{"$part": <name>}`, stored once in a `parts` chunk under
the SHA-256 of its text (a part names the parts inside it). Values that share
a sub-value (the gene names of every analysis result, a user field in each
later value) store it once. A value's hash is the SHA-256 of the value
itself, joined, as canonical text: keys sorted by UTF-16 code units, numbers
and strings written as JavaScript's `JSON.stringify` writes them. A part is
storage only: `read_session` joins each part once, never trusts its name (the
joined value is checked against its own hash), and refuses a part that names
itself through others, a value whose text would pass 2^29 − 24 characters
(refused before it is expanded), and a part, array or value nothing refers to.
A step's parameters may name arrays and parts only the lazy chunks hold;
`read_session` reads the whole file and joins them. The cell-set decoder is
`cellucid.session_codecs.decode_cell_set`.

A session of another format version is refused with the web app's sentence
naming it (`other_session_format_message()`): *This session was saved by an
earlier version of Cellucid (session format 3), which this version no longer
reads: a session file now holds the history of every step (session format 4).
It was not opened. Open it in the version of Cellucid that saved it, or
re-create the session in this version.* Any other file that does not start
with `CELLUCID_SESSION/4\n` is refused with: *This file is not a Cellucid
session file of format 4 (it does not start with CELLUCID_SESSION/4).*

### Closed value shapes

Applying reads four state domains, each a closed shape
(`app/state/domains/README.md` in the web app): `cells` (`{open, genes}`),
`field-overlays` (`{renamedFields, renamedCategories, deletedFields,
userFields}`), `highlights` (`{active, pages}`) and `gene-sets` (`{sets}`). A
key this version does not know is refused before anything is written: a
value it does not understand is never applied in part.

---

## API reference

### `read_session`

```{eval-rst}
.. autofunction:: read_session
```

### `CellucidSession`

```{eval-rst}
.. autoclass:: CellucidSession
   :members:
   :show-inheritance:
```

### The reading model (`cellucid.history`)

```{eval-rst}
.. autoclass:: cellucid.history.Tracker
   :members:
   :no-special-members:

.. autoclass:: cellucid.history.Branch
   :members:
   :no-special-members:

.. autoclass:: cellucid.history.Step
   :members:
   :no-special-members:

.. autoclass:: cellucid.history.Diff
   :members:
   :no-special-members:
```

### `apply_cellucid_session_to_anndata`

```{eval-rst}
.. autofunction:: apply_cellucid_session_to_anndata
```

### `ApplyPlan`, `ApplyReport` and `history_applied`

```{eval-rst}
.. py:class:: ApplyPlan

   What applying one step to an object would change, written nowhere yet: the
   history model's ``ApplyPlan`` (``session_sha256``, ``tracker``, ``branch``,
   ``step``, ``requires``, ``rows``), each row an ``ApplyPlanRow`` (``kind``,
   ``domain``, ``name``, ``column``, ``from_``, ``cells``, ``status``,
   ``note``).

.. autoproperty:: cellucid.ApplyPlan.table

.. autoproperty:: cellucid.ApplyPlan.writes_anything

.. automethod:: cellucid.ApplyPlan.apply

.. py:class:: ApplyReport

   What an apply did: its ``plan``, whether it wrote anything (``applied``),
   and the record it wrote (``provenance``, the history model's
   ``Provenance``; ``None`` when nothing was written).

.. autoproperty:: cellucid.ApplyReport.adata

.. autofunction:: history_applied
```

---

## Edge cases (do not skip)

### Dataset mismatch (most common)

A session is only meaningfully applicable to the dataset it was saved on.
Common mismatch causes:
- you filtered or reordered cells before saving the session,
- you regenerated the AnnData with different preprocessing,
- you are applying a session to a subset.

Any mismatch raises before the target is changed.

### Cell ordering assumptions

Highlights and field values are stored by **original cell position** (the
row in the order the dataset gives the cells). If your AnnData row order
differs, labels would land on the wrong cells: when `adata` has the embedding
the session names, the cell-order digest refuses it (a digest is a sum over
the cells, so a mismatch says the cells differ, not which). An AnnData read
from a Cellucid-ready file in its stored priority order is refused until
`cellucid.restore_original_order` puts it back.

A session saved on a **subsample** records which cells were open (the
manifest's `cells`: the specification, the kept cells and their digest), and
a step taken while a subsample was open records it in its `cells` state.
Applying it checks that `adata` has the dataset's full cell count; the state
then applies to the cells open both when it was saved and when the step was
taken. A cell that was not open is `<NA>` in a highlight column and missing
in a field's -- unknown, not "not selected".

A session saved while genes were read for a sample of the open cells (a
large dataset stored cell by cell opens every cell on the map and reads a
gene for the first of them) records that sample:
`CellucidSession.genes` is `{"cells", "order", "digest"}` (how many of the
open cells, the file's cell order they were open in -- every cell's, or a
share's spread through a file in its own order -- or `None` for a
subsample's own order, and their cell-set digest), or `None`. A session
holds no gene values, so applying it is the same either way; the record says
which cells the viewer's gene results (a signature score, a marker table)
were computed on.

### Columns the object has already

A highlight or user-defined field column whose name the object already uses
is never replaced when Cellucid did not write it: its row is `not-applied`
and says so, and the rest applies. A column an earlier apply wrote is
`unchanged` when it holds the same values, and is replaced only with
`overwrite=True`. Gene sets follow the same rule.

### Backed AnnData

Applying changes the object's `obs` and `uns` in memory. A backed AnnData
keeps `X` on disk, and its source file is not rewritten: write the object to
keep the result.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “This session was saved on …, and this object has …”
Likely causes:
- You are applying the session to a different AnnData (another cell count,
  gene count, or dataset id in `uns["cellucid"]["id"]`).
- The viewer was given another `dataset_id=` than the one the object declares
  (“saved on the dataset 'x', and this object is 'y'”): the session is saved
  under the id the viewer was given.

Fix:
- Apply to the original AnnData the session was saved on.
- For a viewer given another id: apply through it
  (`viewer.apply_session_to_anndata(adata, request)` names the id itself), or
  name the id: `apply_cellucid_session_to_anndata(adata, path, dataset_id="x")`.

### Symptom: “This object's cells are not the cells the session was saved on, in the same order …”
The cell-order digest of the session's embedding differs: the rows are other
cells, or in another order. Put the cells in their original order, or apply
the session to the object it was saved on.

### Symptom: “This session was written by a newer Cellucid (…) than this package reads (…)”
The web app that saved it is newer than this package. Update the `cellucid`
package.

### Symptom: “This session was saved by an earlier version of Cellucid (session format 3) …”
A session saved before Cellucid kept the history of every step is not read.
Open it in the version of Cellucid that saved it, or re-create the session.

### Symptom: “This file is not a Cellucid session file of format 4 …”
Fix:
- Confirm the file is actually a `.cellucid-session` file and not (for
  example) a `.json` file or a renamed file.

### Symptom: “This session file is damaged: …”
The file breaks the format where the sentence says (a chunk of another size
than declared, a value that does not match its hash, an array no value refers
to, trailing bytes...). Treat it as corrupted or tampered with, and save the
session again from the viewer.

### Symptom: “This session was shared as its steps only …”
The share left the values out, so there is no state to apply. Replay it in
Cellucid on its dataset, or ask for a share that keeps the values.

### Symptom: “The state at step … uses the obs column … which this dataset does not have”
The chosen step's state colors by, filters on or binds a highlight group to
fields the AnnData lacks: every missing obs column and gene is listed. Apply
the session to the AnnData it was saved on (or with the same data mapping and
gene names).

### Symptom: “There is no step …”, “There is no branch …” or “There is no tracker …”
The `step`, `branch` or `tracker` is not in the session; `session.trackers`,
`tracker.branches` and `tracker.steps` list what it holds.

### Symptom: a row is `not-applied`
Its `note` says why: a column or gene set the person made, a column or
category the object does not have (a column dropped since, a modality's
column of a MuData), a deletion, a presentation domain. Everything else was
applied; `cellucid.history_applied(adata)` lists what each apply skipped.

### Symptom: “The object changed since this plan was made …”
Something changed `adata` between the dry run and `plan.apply()`. Make the
plan again and apply that.

### Symptom: “No highlight columns were added”
Likely causes:
- The state at that step holds no highlight group.

How to confirm:
- `tracker.step(step).state()["highlights"]["pages"]`.

---

## See also

- {doc}`jupyter` for asking a notebook viewer for its session
- {doc}`round_trip` for pulling everything (highlights, labels, scores, gene sets) from a viewer or the web app's exports in one call
- {doc}`../../e_jupyter_hooks/index` for deeper hook/event documentation
