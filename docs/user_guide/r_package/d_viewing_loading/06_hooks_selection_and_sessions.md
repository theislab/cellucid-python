# Hooks, Selections and Sessions

**Audience:** analysts who want the viewer and R to work together
**What you'll do:** act on an open viewer, react to what you do in it, and
bring a selection or a whole session back to R.

Every viewer has a channel (a WebSocket) to the R session that opened it.
Actions sent before the page has connected wait and go out in order as soon as
it does.

## Actions

Everything a person can do in the viewer -- and everything an agent reads
from it -- is an action of one catalog, and every action is a function of
this package, generated from the catalog (`inst/action_catalog.json`, the
web app's own, with its help page: `?cellucid_highlight_cells`):

```r
viewer <- cellucid_view(seu)

cellucid_color_by(viewer, field = list(source = "obs", key = "seurat_clusters"))
cellucid_highlight_cells(viewer, cells = WhichCells(seu, idents = "B"), color = "#377eb8")
cellucid_cells_visibility_set(viewer, cells = seu$nCount_RNA < 500, visible = FALSE)
cellucid_highlight_clear_all(viewer)
cellucid_camera_reset(viewer, view = "live")
cellucid_mapping_set(viewer, modality = "ADT")   # a Seurat object with an ADT assay
```

Each is one step of the viewer's history, labeled as R's, undone like a
person's click. `cellucid_act(viewer, "color.by", field = ...)` sends any
action named in data (its arguments are written `.viewer` and `.action`, so
no action parameter is partially matched to them). Fields, genes and
categories are named, never positions. Cells are given by name, by 1-based
position, or as `TRUE`/`FALSE` per cell, and travel as original positions in
one binary frame. Every action returns a request, invisibly, at once:

```r
request <- cellucid_fields_list(viewer)
fields <- cellucid_action_result(request)   # in a notebook: in a later cell
request$step                                # the history step a change made
request$unrecorded                          # or why a change made none
```

`cellucid_action_result(request, timeout = 60)` returns the outcome (a read's
result, as lists) or stops with the viewer's reason when it did not apply the
action. A change the viewer applied but recorded no history step for keeps
why in `request$unrecorded` -- `"unchanged"` (it was so already),
`"retracted"` (it put back the step just before it, which was taken out),
`"off"` (the History settings do not record this kind of change),
`"paused"` (recording is paused), `"elsewhere"` (another tracker was shown
meanwhile) or `"no-history"` (the dataset's history has not started) -- and
`cellucid_action_result()` returns its outcome with a warning of class
`cellucid_action_not_recorded` ("cellucid_color_by(): it changed nothing (it
was so already), so no step was recorded."); a read records no step by
nature and warns nothing. A refusal nobody collects arrives as a `command_error` naming its
request (a warning unless a hook handles it). `cellucid_mapping_set()` changes
only the groups it names; the mapping the viewer then shows becomes
`viewer$data_mapping`. An object is served with one
matrix per assay (its `data` layer, else `counts`) and no layers, so the
mapping's `expression` stays `"X"`; change the modality, the graph (`knn`) or
the embeddings instead. A file given by path may hold layers
(`expression = "layer:counts"`).

```{figure} ../../../_static/screenshots/r_package/cellucid-view-driven-from-r.png
:alt: Cellucid opened by cellucid_view() on a Seurat object in an R session, after R called cellucid_color_by(viewer, ...) and cellucid_highlight_cells(viewer, ...) on the Ductal cells: the points are colored by clusters with the Ductal cells highlighted, and the Dataset box names the object.
:width: 1440px

`cellucid_view()` on a Seurat object of the Pancreas sample, opened in a
browser tab, after R ran `cellucid_color_by()` on `clusters` and
`cellucid_highlight_cells()` on the Ductal cells. The Dataset box reads *Zarr v2
folder · R*: the object is served from memory by this R session.
```

## The history

Each action R sends is one step of the viewer's history, labeled **R** in the
viewer's **State › History** list, beside the person's own steps and an
agent's. R reads the history and moves it as the viewer's own buttons do:

```r
history <- cellucid_action_result(cellucid_history_get(viewer))
history$head                     # the step shown now
length(history$steps)            # the steps from the first to the one shown
cellucid_history_undo(viewer)    # back one step, as Undo does
cellucid_history_redo(viewer)    # forward again
cellucid_history_jump(viewer, step = "s3")   # show any step history_get() lists
```

`cellucid_history_get(viewer, all_branches = TRUE)` lists every branch, not
only the path to the step shown; each step says what it did (`label`), who
took it (`origin`, with `agent` or `host`), its pseudocode line (`code`), when,
what it changed (`touches`) and whether it no longer applies (`status`,
`note`). Undo,
redo and jumps record no step; their outcome is `head` (and `moved`, false
when there was nothing to move to). The history is changed through the
generated `cellucid_history_*()` functions (`cellucid_history_note_add()`,
`cellucid_history_branch_switch()`, `cellucid_history_turn_undo()`, …), the
same actions the State panel and agents use. How the history works in the
viewer: {doc}`../../web_app/m_state_history/index`.

### Turns

```r
turn <- cellucid_turn(viewer, "Mark the B cells", {
  cellucid_color_by(viewer, field = list(source = "obs", key = "seurat_clusters"))
  cellucid_highlight_cells(viewer, cells = WhichCells(seu, idents = "B"))
})
cellucid_history_note_add(viewer, target = list(turn = turn), text = "B cells by cluster.")
cellucid_history_turn_undo(viewer, turn = turn)   # every step of the turn, at once
```

`cellucid_turn(viewer, label, code)` makes every action sent while `code`
runs one turn of the history: the State panel shows its steps as one group,
labeled `label` and by R, and undoes them together, as it does an agent's
turn. It returns the turn's id (`"r-"` and 12 hexadecimal digits),
invisibly. Turns do not nest; a turn whose code stops on an error closes all
the same.

### What changed, and a session pushed in

```r
cellucid_on(viewer, "history_changed", function(change) message(change$kind, " ", change$head))
cellucid_last_history_change(viewer)   # list(kind, tracker, head, step, branch)

request <- cellucid_load_session(viewer, "pbmc.cellucid-session")
cellucid_load_session_result(request)  # list(tracker = list(id, name), head)
```

The viewer tells every change of its history as a `history_changed` event:
`kind` is `"step"` (a step recorded; `step` is its id), `"moved"` (undo,
redo, a jump, a branch switched to), `"edited"` (a note, a bookmark, a branch
renamed, locked, archived or deleted), `"tracker"` (a tracker started, used,
renamed or removed) or `"loaded"` (a session loaded); `tracker` is the
tracker in use, `head` the step shown and `branch` the active branch.
`cellucid_last_history_change()` returns the latest (`NULL` before the
first).

`cellucid_load_session(viewer, path)` sends a `.cellucid-session` file saved
on the dataset the viewer shows: its trackers are added beside the person's
own, and the viewer shows the file's active tracker, as opening the file in
the viewer does. `cellucid_load_session_result(request, timeout = 60)`
returns what the viewer shows then, or stops with the viewer's reason
(another dataset's session, a damaged file); in a Jupyter notebook, collect
it in a later cell.

## Events

```r
cellucid_on(viewer, "selection", function(event) {
  message(length(event$cells), " cells selected with the ", event$source)
})
cellucid_on(viewer, "click", function(event) str(event))
```

The events are `selection`, `hover`, `click`, `ready`, `data_mapping`,
`command_error`, `action_result` (an action's answer, which also completes its
request), `pong`, `debug_snapshot`, `session_bundle`, `result`,
`result_removed` and `results` (what the viewer made, as
`cellucid_results()` lists it), `color_by` (the field it colors by),
`history_changed` (every change of its history) and `session_loaded` (a
session `cellucid_load_session()` pushed in, loaded); `message` catches every
one of them, before its own handlers, with its kind as `event` (Python's
`on_message`). Any other name is refused, as Python's hooks refuse it: a
misspelt event would never be heard. A handler runs
whenever R is idle; one that fails gives a warning and the others still run.
`cellucid_off(viewer, "hover")` removes a kind's handlers.
`cellucid_wait_event(viewer, "click", timeout = 30, predicate = NULL)` waits
for the next event of a kind (one the predicate accepts) and returns it, as
Python's `viewer.wait_for_event()`; like every wait, it cannot work inside a
Jupyter notebook cell, where a handler does.

Cells are 1-based positions in every event, as the actions take them: a
`selection` gives `cells` and `cell_names`, a `hover` or `click` gives `cell`
and `cell_name`. (Python's events are 0-based, as Python indexes.) `pong` and
`debug_snapshot` answer `cellucid_debug_connection(viewer)`, the report to read
first when a viewer stays blank or does not react; the channel's own
keep-alive pings, every 20 seconds, raise no event.

## A selection back in R

```r
selected <- cellucid_selection(viewer)
selected$cells        # 1-based positions
selected$cell_names   # the object's cell names
sub <- subset(seu, cells = selected$cell_names)
```

A selection of any size arrives at once (millions of cells decode in a fraction
of a second).

## The session

```r
request <- cellucid_request_session(viewer)
session <- cellucid_session_result(request, file = "pbmc.cellucid-session")
cellucid_trackers(session)
tracker <- cellucid_tracker(session)
cellucid_steps(tracker)[, c("id", "branch", "label", "origin")]
```

`cellucid_request_session()` asks the viewer for its whole session -- its
history: every step a person, an agent or R took, each tracker and branch,
its notes, bookmarks and turns, and the state each step leaves (highlights,
colors, filters, annotations, the camera) -- and `cellucid_session_result()`
saves it as a `.cellucid-session` file (session format 4: the whole history,
every tracker and branch with its notes and Recently deleted, without the
conversations; the viewer's **Share…** writes a chosen part of the same) and
returns it read and checked whole: the
one way, in an R console and in a Jupyter notebook (there, collect it in a
later cell: the viewer answers between cells). Open it again in the web app
with the same dataset, under the same dataset id.
`cellucid_session_result()` stops with what to do when the viewer is not open
or does not answer within `timeout` seconds; the request stays open, so it can
be collected again.

A saved file is read with `cellucid_read_session(path)`, which checks it
whole first: a damaged file is refused with "This session file is damaged: "
and what is wrong; a file of format 3, saved by an earlier Cellucid, is
refused with the viewer's own sentence; and a file written by a newer
Cellucid than this package with "This session was written by a newer
Cellucid (X) than this package reads (Y); update the cellucid package to
read it." The history reads as data frames and records, function for
function cellucid-python's `cellucid.history`:

| Function | What it returns |
|---|---|
| `cellucid_trackers(session)`, `cellucid_tracker(session, tracker = NULL)` | Every tracker, one row each (`id`, `name`, `active`, `removed`, `start`, `head`, `active_branch`, `main_branch`, `recording`, `paused_at` (both derived from the tracker's `paused`: `NULL` while it records, else `at`, `by` and `draft`), `steps`, `branches`, `groups`, `created_at`); one tracker by name or id (the active one by default) |
| `cellucid_branches(tracker)`, `cellucid_branch(tracker, branch = NULL)` | Every branch (`id`, `name`, `base`, `tip`, `main`, `locked`, `archived`, `n_notes`, `created_at`); one branch's record |
| `cellucid_groups(tracker)`, `cellucid_group(tracker, group)` | Every group, a named run of consecutive steps of one branch (`id`, `name`, `branch`, `first`, `last`, `steps`, `n_notes`, `created_at`); one group's record by name or id, with its steps' ids |
| `cellucid_steps(tracker, branch = NULL)` | Every step in step order, or a branch's path from the root to its tip: `id`, `branch`, `parent`, `kind` (`root`, `action`, `outside`, `squashed`), `type`, `v` (the version of the action's parameters), `label`, `code`, `origin`, `agent`, `model`, `host`, `turn`, `at`, `status`, `stale`, `bookmark`, `group`, `n_notes` -- the columns of Python's `Tracker.steps` |
| `cellucid_step(tracker, step)`, `cellucid_step_state(tracker, step)` | One step's record (a squashed step's `squashed` says why it was made, who acted and what the replaced steps carried: `carried`, their note ids and bookmarks; a squash in mid-branch keeps the last replaced step's id); the state after it, a named list by state domain |
| `cellucid_step_diff(tracker, a, b, other = NULL)` | What differs between two steps, one row per change (`domain`, `name`, `line`, in the web app's words) |
| `cellucid_stale(tracker)` | The results stale at the head: the step, the later step that changed an input it read, and those inputs |
| `cellucid_notes(tracker)`, `cellucid_turns(tracker)`, `cellucid_journal(tracker)` | Every note; every turn (a notebook's or an agent's group of steps); every edit of the history itself |
| `cellucid_methods(tracker, branch = NULL)`, `cellucid_pseudocode(tracker, branch = NULL)` | A branch (main by default) as methods text and as pseudocode, as the web app writes them: the tracker's first note opens the methods; groups are headings in the methods and `# --- QC ---` … `# --- end of QC ---` in the pseudocode; a bookmark is `# Bookmark: <name>` (one a squashed step kept from a merged step adds "(on a step merged into this one)"); a squashed step is said as what it stands for and holds only values |
| `cellucid_conversations(session, tracker = NULL)` | The agent conversations the file kept with a tracker |

`?cellucid_trackers` and `?cellucid_read_session` describe every column and
field. A file shared as its steps only (without the values they leave
behind) reads, but has no states: `cellucid_step_state()` refuses it.

## Writing a step back onto the object

```r
session <- cellucid_read_session("pbmc.cellucid-session")
plan <- cellucid_apply_session(seu, session)   # dry run: one row per change
plan
seu <- cellucid_apply_plan(plan)
cellucid_history_applied(seu)

# Another step, a branch's last step, or the viewer's session, applied at once:
seu <- cellucid_apply_session(seu, session, step = "s5", dry_run = FALSE)
seu <- cellucid_apply_session(seu, session, branch = "Without score", dry_run = FALSE)
seu <- cellucid_apply_session(seu, cellucid_request_session(viewer), dry_run = FALSE)
```

`cellucid_apply_session(data, session, tracker = NULL, branch = NULL,
step = NULL, ..., dry_run = TRUE, overwrite = FALSE, keep_session = FALSE,
timeout = 60)` applies the state one step left -- by default the step the
viewer showed (the active tracker's head); `branch` takes a branch's last
step, `step` any step (on `branch` when both are given) -- to the object it
was saved on. `session` is a session `cellucid_read_session()` read, a file,
or a `cellucid_request_session()` request (collected first; in a notebook, in
a later cell). The rules are cellucid-python's
`apply_cellucid_session_to_anndata()`, row for row: both are held to one
shared expected plan.

**A plan first.** With `dry_run = TRUE` (the default) nothing is written:
the `cellucid_apply_plan` returned lists every change, one row each, in
domain order -- `column-renamed`, `categories-renamed`, `categories-merged`,
`field-deleted` (recorded, never applied: the column stays), `user-field` (a
column per field the session made: a factor, or numbers), `gene-set`,
`highlight` (one `TRUE`/`FALSE` column per group,
`cellucid_highlight__<page>__<group>` from the page's name and the group's
label) and `presentation` -- with its `status`:

| Status | Meaning |
|---|---|
| `write` | The object lacks it; applying writes it |
| `unchanged` | The object holds exactly it already |
| `replace` | A column or gene set an earlier apply wrote holds other values, and `overwrite = TRUE` |
| `not-applied` | Not written; `note` says why (a column the person made, an earlier apply's column without `overwrite`, a deletion, presentation) |

`cellucid_apply_plan(plan)`, or `dry_run = FALSE`, writes the plan all at
once: labels renamed in place, each user-defined field and highlight group a
column of the cell metadata (Seurat `meta.data`, `colData`), the gene sets in
`cellucid$gene_sets` (origin `"session"`) and, each set written or replaced,
one `TRUE`/`FALSE` column of the feature metadata (below). R copies on change: the object
passed in is not changed, and the one returned carries the attribute
`"cellucid_apply_report"`, `list(plan, applied, provenance)`. A plan whose
object changed since it was made is refused: make it again.

- **Applying again changes nothing.** The same step applied twice gives
  every row `unchanged` or `not-applied`, and `applied` is `FALSE`.
- **A column the person made is never replaced**, whatever `overwrite` says:
  its row says so and the rest applies.
- **Presentation is reported, never applied.** Each other domain whose value
  at the step differs from the tracker's opening (filters, colors, views,
  the camera, analyses, the expression values) is one `presentation` row:
  how the data was shown, not the data.
- **Each apply that writes is recorded.** One provenance record is appended
  to `cellucid$history_applied` (Seurat `misc$cellucid`;
  `metadata(x)$cellucid`): the session's SHA-256, the tracker, branch and
  step, when, the package and its version, the columns it wrote, what it
  renamed, what it did not apply, and the expression values the genes were
  shown in at the step (`expression`: `list(normalize, log1p)`, `NULL` when
  the session holds none) -- what the person looked at, recorded and never
  applied: the object's expression values stay as stored, and
  `NormalizeData()` (`scale.factor = 1e4`) gives normalized and log1p
  together. `cellucid_history_applied(seu)` reads
  them, oldest first; `keep_session = TRUE` also keeps the session file
  itself (base64), so the history can be rebuilt from the object alone.
  cellucid-python's `uns['cellucid']['history_applied']` holds the same
  records.
- **Gene sets mark their genes.** Each set written or replaced is one
  `TRUE`/`FALSE` column, `cellucid_gene_set__<name>`, of the feature
  metadata of the gene axis the data mapping reads: the assay's for a Seurat
  object (`meta.features`, or an Assay5's `meta.data`), `rowData()` of the
  main experiment, or of the `altExp()` the mapping names, for a
  SingleCellExperiment or SummarizedExperiment (a data frame has none).
  `TRUE` on each row whose gene, as the mapping names genes, the set lists;
  a gene the object does not have is left out, with a warning. The name is
  cellucid-python's: the set's name with each run of characters other than
  ASCII letters and digits one `_`, trimmed (`gene_set` when nothing is
  left), `_2`, `_3` on a clash, given over every set of the step in order. A
  mapping that reads the raw matrix or no expression gets no column (the
  sets stay in the record, with a warning). A set an earlier apply wrote and
  that is replaced (`overwrite = TRUE`) has its column replaced in place; a
  column of that name beside a set no apply wrote is the object's own and is
  never replaced (a warning says so). A set the step no longer holds gets no
  row: its record and column stay. The columns are not rows of the plan; the
  shared fixture pins them (`apply.sequence[].var_columns`), as it pins the
  web app's **Save AnnData with this state…**.
- **Colors go with the columns**, where the object keeps what an AnnData
  keeps in `uns` (Seurat `misc`, `metadata(x)` -- as zellkonverter carries
  `uns` -- or a data frame's attributes), under Scanpy's name,
  `<column>_colors`: a categorical field gets the colors it is shown in, one
  `#rrggbb` a level in the factor's order, every level, used or not (in the
  view in focus at the step, a level's color the person or an agent chose,
  else its `palette` -- a field made from highlight pages opens in its pages'
  colors -- else Cellucid's own category colors); any other column written
  gets none (colors left under its name described
  another column and are removed). A column of the object (a factor, such as
  `seurat_clusters`) recolored in the view in focus gets the colors it is
  shown in too, under its name after the step's renames and in its levels'
  order after the step's relabels: a level's color the person or an agent
  chose, else the color it opens with -- its own `<column>_colors` as the
  viewer draws them (`#rrggbb` in lowercase), when there is one a level,
  else Cellucid's -- one `column-colors` row of the plan each (`unchanged`
  when the object holds them already; `not-applied`, saying why, for a
  column that is not a factor or one with none of the levels recolored). A
  renamed column carries its colors (a
  character vector; another form is removed), merged levels keep their first
  level's color (colors that did not cover the levels are removed), a text
  column relabeled loses them, and a column left as it is keeps them. The
  rule is cellucid-python's ({ref}`apply-column-colors`); `cellucid_pull()`
  of a session file writes a field's colors the same way, under the column
  name it writes, and `cellucid_pull()` from a live viewer writes a
  categorical field with the colors its view in focus shows it in (the
  viewer's result carries them).
- **The priority-order column is never touched.** A plan that would write,
  rename or replace the cell-metadata column `cellucid_original_position`
  (each cell's original position in a store written in priority order) is
  refused before anything is written, dry runs included.
- **Only what the chosen step needs is checked**, before anything is written:
  the cell count, the dataset id (the one the object was shown under:
  `dataset_id =` when a viewer was given another id than the one the object
  declares -- `cellucid_apply_session(object, request)` names it itself for
  an object declaring the id the viewer's data declared -- else the one its
  cellucid record names, when it names one), the data mapping the session was saved with and its gene count, the
  cells' order when the object has the embedding the session names (a
  cell-order digest of its stored values), and the obs columns and genes the
  step's state uses.

A step's state knows the cells the viewer had open then: a cell outside them
is `NA` in a highlight column and in a field's, never "not selected". A
session shared as its steps only cannot be applied. An object read from a
Cellucid-ready store in its stored priority order is refused until
`cellucid_restore_original_order()` puts its cells back. If `cellucid_view()`
was given arguments that change what is shown (`embeddings`, `obs`,
`modalities`, …), pass the same ones after the session.
