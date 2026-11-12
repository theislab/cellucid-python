# Sessions: apply to AnnData (`apply_session_to_anndata`)

This page documents how to take what was done in the viewer (a
`.cellucid-session` file: the history of every step) and apply the state at
one of its steps back onto an `AnnData`.

```{tip}
To bring results from a **live viewer** onto your AnnData, `viewer.pull(adata)`
is simpler: one call, no session, every highlight, label, score and gene set
({doc}`../g_api_reference_coverage/api/round_trip`). `cellucid.pull(adata,
"exports/")` reads a saved `.cellucid-session` (the state it shows) together
with the web app's CSV and GMT exports. This page is about applying a session
itself: any step of it, with its renames.
```

The main workflow is:

1. Get a session (no download): {doc}`10_session_bundles_request_session_bundle`,
   or save one from the web app (**State › Session › Share…**)
2. Plan the apply: what applying one step would change, one row each, written
   nowhere yet:
   - `session.apply(adata)` (a {class}`~cellucid.CellucidSession`), or
   - `cellucid.apply_cellucid_session_to_anndata(adata, session_or_path)`, or
   - `viewer.apply_session_to_anndata(adata, request)` (collect what a
     request brought and plan it; in a notebook, in a later cell than the
     request)
3. Apply the plan: `plan.apply()` writes it into `adata` and records the
   apply.

## At a glance

**Audience**
- Wet lab / beginner: this is optional unless you’re collaborating or exporting curated labels back into analysis.
- Computational users: this is the bridge from “interactive UI edits” → “analysis-ready AnnData”.
- Developers: read the identity checks, the plan's rows and the provenance record.

## Quickstart (copy/paste)

### In a notebook (two cells)

```python
request = viewer.request_session()
```

```python
import cellucid

plan = viewer.apply_session_to_anndata(adata, request)   # dry run: one row per change
plan.table
report = plan.apply()                                     # writes into adata, records the apply
cellucid.history_applied(adata)
```

`viewer.apply_session_to_anndata` collects the session (`request.result()`)
and plans it. Its keywords are the apply's (`tracker`, `branch`, `step`,
`dry_run`, `overwrite`, `keep_session`), and two of its own: `timeout` (the
wait for the session) and `cleanup_file` (delete the session's temporary file
once applied with `dry_run=False`; a plan leaves it in place, and
`request.result().save(...)` keeps a copy).

### From a session file

```python
session = cellucid.read_session("my-session.cellucid-session")
plan = session.apply(adata)
report = plan.apply()
```

`dry_run=False` applies at once and returns the report:
`session.apply(adata, dry_run=False)`.

### Another step, another branch, another tracker

```python
tracker = session.tracker()                       # the tracker the viewer showed
tracker.steps[["id", "branch", "label", "origin"]]

before_merge = session.apply(adata, step="s4")
without_score = session.apply(adata, branch="Without score")
exploration = session.apply(adata, tracker="Exploration")
```

By default the state applied is the one the viewer showed: the active
tracker's head, on its active branch. `branch` applies a branch's tip, or
`step` when that step is on the branch's path; `step` alone applies that step
on the branch it was made on. Trackers and branches are named by their name
or their id; an unknown one is refused, naming the ones there are.
`session.trackers` lists every tracker.

```{important}
Session application is **index-based**:

The session references cells by **original position** (0-based): the row of a
cell in the order the dataset's cells were given, whatever order a file stores
them in. Apply a session to an `AnnData` in that order: the object you prepared
or served, or a Cellucid-ready file read back through
{func}`~cellucid.restore_original_order` (an `AnnData` still in a
Cellucid-ready file's stored order is refused). When `adata` has the
embedding the session names, its cell order is checked (the `cellOrder`
section below).
```

```{note}
The target may also be a `mudata.MuData` (for a viewer opened on one). The
session is then checked against the MuData as the viewer showed it — the
modality its `dataMapping` names, with every modality's own `obs` columns as
`<modality>:<column>` — and the labels and columns it applies go to the
MuData's global `obs`; `mdata.uns` records the apply.
```

## The plan

A plan ({class}`~cellucid.ApplyPlan`) holds the ids it chose
(`session_sha256`, `tracker`, `branch`, `step`), `requires` (the obs columns
and genes the step's state uses) and `rows`. `plan.table` is the rows as a
DataFrame, one row per change, in domain order:

| Column | Meaning |
|---|---|
| `kind` | `column-renamed`, `categories-renamed`, `categories-merged`, `field-deleted`, `user-field`, `gene-set`, `highlight`, `column-colors` or `presentation` |
| `domain` | the state domain it comes from (`field-overlays`, `gene-sets`, `highlights`, `colors`, or a presentation domain) |
| `name` | the column, category or gene set it writes (a presentation row: the domain's name) |
| `column` | the column a category is in |
| `from` | the old name or category |
| `cells` | the cells it labels |
| `status` | `write`, `unchanged`, `replace` or `not-applied` |
| `note` | why, in a sentence |

The statuses:

- `write`: the object lacks it; applying writes it.
- `unchanged`: the object holds exactly it already. Applying the same step
  again changes nothing and says so: every row `unchanged` or `not-applied`,
  `plan.writes_anything` false, and the report's `applied` false.
- `replace`: a column or gene set an earlier apply wrote, holding other
  values, with `overwrite=True`.
- `not-applied`, with the reason in `note`: a column or gene set the person
  made is never replaced, whatever `overwrite` says; a column or category the
  object does not have; a deletion; a presentation domain; a column or gene
  set an earlier apply wrote with other values, without `overwrite`.

`plan.apply()` makes the plan again first and refuses an object that changed
since ("make the plan again"), then writes every row to write or replace at
once. It returns an {class}`~cellucid.ApplyReport`: `plan`, `applied` (whether
anything was written), `provenance` (the record written, or `None`) and
`adata`. Every check runs before anything is written, so a failed apply
leaves the object as it was. The object is written in place.

## What the rows write

As implemented in `cellucid-python/src/cellucid/anndata_session/`:

### 1) Labels → renamed `adata.obs` columns and categories

- `column-renamed`: an obs column the session renamed is renamed, keeping its
  place. A gene's display name is reported, never applied.
- `categories-renamed`: a category of a categorical obs column is renamed in
  place (a text column has its values renamed).
- `categories-merged`: a category renamed to a label the column already has,
  or to the label another rename gives: the two become one category.
- `field-deleted`: a field deleted in the session is reported, never applied
  (the column stays).

A rename naming a column or category `adata` does not have, or a new name the
object already uses for a column it did not get from Cellucid, is
`not-applied` with the reason.

### 2) User-defined fields → `adata.obs` columns

One column per field the session made, named by the field's key:

- a categorical field: a `category` column with the field's categories in
  order; a cell with no category is missing;
- a continuous field a notebook pushed: `float32`, `NaN` for no value;
- a continuous copy of another field copies that field: an obs column, or a
  gene found by name among the gene names the session's `dataMapping` shows
  (made unique with the shared rule), read from the mapping's expression
  matrix — `X`, `raw.X`, or a layer.

A field deleted in the session, and a purged one (which holds no values), is
`not-applied`.

A categorical field is written with the colors it is shown in, where Scanpy
keeps a column's colors, so `sc.pl.umap(adata, color=key)` and the viewer show
the same colors:

```text
adata.uns["<key>_colors"]  # ['#808080', '#0077c8', ...], one a category, in the column's order
```

Each category's color is the one it has in the view in focus at the applied
step: the color the person or an agent gave it there (Colors, `color.category`),
else the color the field opens with — its `palette` (a field made from
highlight pages, in Create categorical, opens in its pages' colors), else
Cellucid's own category colors by position. Every category has one, used or
not. A column that came with the dataset and was recolored in the view in
focus (`clusters`, `cell_type`) is written the same way, one `column-colors`
row each: the categories recolored, the rest in the colors it opens with (its
own `uns` colors, else Cellucid's). The plan's Colors row stays `not-applied`
(every other color — a color map, a gene's, another view's — is how the data
was shown) and says so when colors went with a column. A column the apply writes
without categories (a highlight group, a continuous field) has no colors:
colors left under its name described another column and are removed. Renamed
columns carry their colors; merged categories keep their first category's
({ref}`the full rule <apply-column-colors>`).

### 3) Highlight groups → boolean `adata.obs` columns

One `bool` column per group of every page, named
`cellucid_highlight__<page>__<group>` from the page's name and the group's
label: ASCII letters and digits kept, every other run one underscore, none at
either end; `_2`, `_3`, … when two read alike.

```text
adata.obs["cellucid_highlight__Highlights__T_cells"]  # True/False per cell
```

(apply-gene-set-var-columns)=
### 4) Gene sets → `adata.uns["cellucid"]["gene_sets"]` and boolean `var` columns

Each gene set of the step, written as `cellucid.pull` writes them (`genes`,
`description`, `origin` `"session"`, `pulled_from`: the session file's name,
`id`). A set of the same name that no session wrote is never replaced.

Each set written or replaced is also one `bool` column of `var`, so the
genes it lists can be selected with the rest of the gene table:

```text
adata.var["cellucid_gene_set__T_cell_markers"]  # True on the set's genes
```

- **The name** is `cellucid_gene_set__` and the set's name by the highlight
  columns' rule: ASCII letters and digits kept, every other run one
  underscore, none at either end (`gene_set` when nothing is left); `_2`,
  `_3`, … on a clash. The names are given over every set of the step, in
  its order, so a set's column is the same whichever of them an apply
  writes.
- **The `var`** is that of the gene axis the data mapping reads: `adata.var`
  when it reads `X` or a layer; for a MuData, the var of the mapping's
  modality (`mdata.mod[modality].var`; the global `mdata.var` is not
  written). A mapping that reads `raw` or no expression has no `var` to mark:
  the sets stay in `uns` only, with a warning.
- **The values**: `True` on each row whose gene, as the mapping names genes
  (`geneNames`), the set lists. A gene the dataset does not have is left out,
  with a warning.
- **Replacing.** A new column goes after `var`'s columns. A set an earlier
  apply wrote and that is replaced (`overwrite=True`) has its column replaced
  in place. A `var` column of that name beside a set no apply wrote is the
  dataset's own: it is never replaced (a warning says so).
- **A set the step no longer holds** gets no row: its record and its column
  stay, as for a user-defined field.

The `var` columns are not rows of the plan and are not in the provenance's
`columns` (which lists obs columns); the plan's rows are the same as without
them. The web app's **Save AnnData with this state…** writes the same columns
({ref}`history-save-anndata`), and the three are pinned by the shared
fixture `history_golden.expected.json` (`apply.sequence[].var_columns`).

### 5) Presentation → reported, not applied

The data mapping, the expression values, views, the active field, filters,
colors, hidden cells, the camera and its path, settings, analyses and the
layout are how the data was shown, not the data. Each such domain whose
value at the step differs from the tracker's opening is one `presentation`
row, `not-applied`. The session file keeps them; `step.state()` reads them.

The **expression values** (Session › Data mapping › Expression values:
normalized to 10,000 counts per cell, log1p, both, or as stored) are a way of
looking at the stored values, so an apply never changes `X` or a layer.
What the person looked at is recorded instead: each apply's record (below)
holds the step's `expression`, `{"normalize": bool, "log1p": bool}`, as the
person asked it (a matrix that does not hold counts was shown as stored; None
for a session that holds no expression values), so a notebook can follow the
same transform, for example with `sc.pp.normalize_total(adata,
target_sum=1e4)` and `sc.pp.log1p(adata)` on a copy.

### 6) Provenance → `adata.uns["cellucid"]["history_applied"]`

Each apply that writes something appends one record to
`adata.uns["cellucid"]["history_applied"]`, a JSON text (a `str`, which
`.h5ad` and Zarr files hold as it is) of a list of records, oldest first.
`cellucid.history_applied(adata)` reads it:

```python
cellucid.history_applied(adata)[-1]
# {"session_sha256": ..., "tracker": "t1", "branch": "b1", "step": "s7",
#  "applied_at": "2026-11-04T09:12:44.512Z", "cellucid_version": ..., "package": "python",
#  "columns": [...], "renamed": [...], "skipped": [...],
#  "expression": {"normalize": True, "log1p": True}, "session": None}
```

`columns` are the columns it wrote (an earlier apply's columns are the ones
a later apply replaces only with `overwrite=True`), `renamed` each rename it
applied (`kind`, `column`, `from`, `to`), `skipped` each row not applied
(`kind`, `name`, `note`), `expression` the expression values the genes were
shown in at the step (recorded, never applied; None when the session holds
none), and `session` the session file itself (base64) with
`keep_session=True`, so the history can be rebuilt from the object alone.
cellucid-r's `cellucid_history_applied()` reads the same records.

## What is checked

Only what the chosen step needs, before anything is written:

| Check | Rule |
|---|---|
| The share | a session shared as its steps only (without the values they leave behind) cannot be applied |
| The step | the tracker, branch and step must be in the session; a step named with a branch must be on its path |
| Cell order on file | an `AnnData` in a Cellucid-ready file's stored order is refused |
| `datasetId` | must equal the id the object was shown under, `dataset_id=` (a viewer given another id than the one the data declares saves its session under that id; `viewer.apply_session_to_anndata` names it itself for an object declaring the id the viewer's data declared), else `uns['cellucid']['id']` when the object names one; checked when the session names one |
| `cellCount` | must equal `adata.n_obs`; for a session saved on a subsample, `adata.n_obs` must be the cells it was drawn from |
| `dataMapping` | every key it names must exist on `adata` (for a MuData, its modality too) |
| `varCount` | must equal the gene count of the mapped axis: `adata.raw.n_vars` for `expression: "raw"`, `0` when the mapping reads no expression, else `adata.n_vars` |
| `cellOrder` | when `adata` has the embedding it names, the cell-order digest of its stored values over the session's cells must be the session's |
| The step's fields | the obs columns and genes the step's views color by or filter on, and its highlight groups are bound to, must exist; all missing ones are named in one error |
| The priority-order column | no row may write, rename or replace `obs['cellucid_original_position']`, nor rename another column onto it (it holds each row's original position in a file stored in priority order); a plan that would is refused, dry runs included |

The fingerprint's `name` is not checked: it is the dataset's display name (or
null), for messages. The manifest's `requires` (`session.requires`) is the
shown state's fields; an apply checks the chosen step's own (`plan.requires`).
A value of the state this version does not read is refused, never applied in
part. A mismatch raises `ValueError` before the target is changed.

A file written by a newer Cellucid than the package is refused with: "This
session was written by a newer Cellucid (X) than this package reads (Y);
update the cellucid package to read it."

In practice, if you subset/shuffle your AnnData after the session was created, you will almost always have a mismatch.

### `dataMapping`: which arrays the session showed

`dataMapping` records which stored arrays the session was showing: `kind`
(`"anndata"`), `embeddings` (the enabled embeddings and the `1d`/`2d`/`3d`
defaults), `knn`, `geneNames` (`"index"` or `"column:<name>"`), `expression`
(`"X"`, `"raw"`, `"layer:<name>"`, or null), `log1p` (the expression matrices
holding `log1p` values, each with its `base`, null for e), `highlyVariable`
(the true/false var column marking the highly variable genes, or null),
`latent`, and `vectorFields`. Python requires each named `obsm` embedding to
exist and be wide enough for its columns, the named `obsp` graph, latent
space, layer or `raw`, gene-name column and every matrix `log1p` names to
exist, and the highly variable column to be a true/false column of the
mapped gene axis.

### `cellOrder`: the cells, in their order

`cellOrder` is exactly `embedding`, `dimension` (1, 2, or 3), and a
16-character lowercase hexadecimal `digest` of that embedding's stored
coordinates, each cell named by its original position. `embedding` is
`{"key", "columns"}`: an `obsm` key with the 0-based columns that were read,
one per dimension. It is what stops a dataset republished under the same name
with the same cell and gene counts but a different row order from silently
re-pointing every saved selection at different cells.

- **It is required.** A session without it is refused when the file is read,
  as is anything other than the three keys above, a `dimension` outside 1–3,
  a malformed `digest`, or an `embedding` whose columns do not fit its kind.
- **Its digest is checked when `adata` has that embedding**, over the cells
  the session was saved on: a mismatch means the rows are other cells, or in
  another order, and the apply is refused. A digest is a sum over the cells,
  so a mismatch says that the cells differ, not which.
- **An object without that embedding** (or with fewer columns) is not checked
  on cell order.

## Edge cases

- A backed `AnnData` is changed in memory (`obs`, `var`, `uns`); the file it is
  backed by is not written.
- A highlight group or field naming a cell beyond `adata`'s cells raises.
- A session saved on a subsample, or a step taken while a subsample was open,
  applies to the cells that were open: a cell that was not is `<NA>` in a
  highlight column (the nullable `boolean` dtype) and missing in a field's --
  unknown, not "not selected".
- Every value applied is a closed shape documented in
  {doc}`../g_api_reference_coverage/api/sessions`; a key this version does not
  know is refused before anything is written.

## Troubleshooting

### Symptom: “This session was saved on …, and this object has …”

This is almost always because:
- you filtered/subset/reordered `adata` after opening the viewer, or
- you’re applying the session to a different dataset.

Fix:
- apply to the exact same AnnData, in its original cell order (`cellucid.restore_original_order` for a Cellucid-ready file read with anndata).

### Symptom: “This session was saved on the dataset 'x', and this object is 'y'”

The viewer was given another id than the one the object declares
(`show_anndata(adata, dataset_id="x")`, `cellucid serve --dataset-id x`), and
saved its session under it. Applied through that viewer
(`viewer.apply_session_to_anndata(adata, request)`), the session goes back to
an object declaring the id its data declared; applied from a file, name the
id it was shown under: `apply_cellucid_session_to_anndata(adata, path,
dataset_id="x")`. The cell count, data mapping and cell order are checked as
always.

### Symptom: a row is `not-applied`

Read its `note`. A column the object already has that Cellucid did not write
is never replaced: rename or drop that column, or leave the row. A column an
earlier apply wrote with other values is replaced with `overwrite=True`.

### Symptom: expected columns weren’t added

Confirm the state at that step holds them:

```python
state = session.tracker().step("s4").state()
state["highlights"]["pages"], state["field-overlays"]["userFields"]
```

## Next steps

- Getting a session from a viewer: {doc}`10_session_bundles_request_session_bundle`
- Session ↔ AnnData conceptual bridge: {doc}`../b_concepts_mental_models/05_sessions_to_anndata_bridge`
