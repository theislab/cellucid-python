# Round trip (results back from the viewer, columns and gene sets into it)

```{eval-rst}
.. currentmodule:: cellucid
```

What a person makes in Cellucid -- highlight groups (lasso, k-NN, proximity,
category and range picks), labels and annotations, clusters made from pages,
gene-signature scores, gene sets -- comes back into your `AnnData` (or
`MuData`) in one call, and what you compute in Python goes into the viewer
the same way. R's `cellucid_pull()`, `cellucid_push()` and friends do the same
for Seurat, SingleCellExperiment and data frames
({doc}`../../../r_package/g_api_reference_coverage/01_public_functions_and_classes`).

| You want… | Call… |
|---|---|
| everything made in a live viewer, as `obs` columns and gene sets | {meth}`AnnDataViewer.pull` (or {func}`pull` with the viewer) |
| the same from the web app's exports, without a viewer | {func}`pull` with the exported files or their folder |
| the cells of a highlight group or the last selection, as an AnnData | {meth}`AnnDataViewer.subset` |
| the last selection, with obs names | {attr}`AnnDataViewer.selection` |
| every gene set the viewer holds | {meth}`AnnDataViewer.gene_sets` |
| a column (labels or a score) shown in the viewer | {meth}`AnnDataViewer.push` |
| a gene set offered in the viewer | {meth}`AnnDataViewer.push_gene_set` |
| a GMT file | {func}`write_gmt` / {func}`read_gmt` |
| cells highlighted from Python, by obs name | {meth}`AnnDataViewer.highlight_cells` |

---

## Nothing waits

The viewer sends what is made as it is made: once its dataset opens (and
whenever its channel opens again) every result and then an inventory, and
each change after. {attr}`AnnDataViewer.results` is that copy, and `pull`,
`subset`, `gene_sets` and `selection` read it at once, with no request. In a
notebook the widget hears the viewer between cells, so a lasso drawn and
`viewer.pull(adata)` run in the next cell finds it; a browser-tab
{class}`AnnDataViewer` works the same way.

- **One page, one dataset.** The copy is of one *epoch*: one page's account
  of one dataset it opened. A reloaded page, or the viewer shown again in
  another output, replaces what an older page sent, whatever its revisions.
  While the person has another dataset open in the viewer (one picked in it
  instead of the notebook's), `results` and everything reading it raise,
  naming the dataset to open again: that dataset's cells are not `adata`'s.
- **Pushes and mappings wait for the viewer, not the cell.** `push`,
  `push_gene_set` and {meth}`AnnDataViewer.set_data_mapping` send at once and
  return. A column of a format only the browser reads goes once the viewer
  has said how many cells it has, and a data mapping once the viewer has
  reported its own -- so a notebook run top to bottom ("Run All") works. The
  field or set is listed in `viewer.results`, and the mapping is
  {attr}`AnnDataViewer.data_mapping`, from a later cell; anything refused
  arrives as a `command_error` (printed unless a hook handles it).
- **Sessions in two steps.** A session bundle is asked for, then collected
  -- the one way, in a notebook (in a later cell) and a browser tab alike:
  `request = viewer.request_session_bundle()`, then `request.result()` or
  `viewer.apply_session_to_anndata(adata, request)` (R:
  `cellucid_request_session()`, then `cellucid_session_result()` or
  `cellucid_apply_session(object, request)`).

## Fast path

```python
import cellucid

viewer = cellucid.show_anndata(adata, dataset_name="My study", dataset_id="my-study-v1")
# … lasso cells, label clusters, score a gene signature in the viewer, then:
summary = viewer.pull(adata)
summary.columns                       # every obs column written, and what it holds
adata.uns["cellucid"]["gene_sets"]    # every gene set, by name

viewer.push(adata.obs["leiden"])      # the clusters, as a field in the viewer
viewer.push_gene_set("T cells", ["CD3E", "CD3D", "IL7R"])
t_cells = viewer.subset(adata, "Lasso selection")
```

Without a notebook, export from the web app (Highlights: **CSV** or
**.h5ad**; a field's **Export** button; a gene set editor's **Export**; Save
session), then:

```python
cellucid.pull(adata, "~/Downloads/cellucid-exports")   # a folder, a file, or a list of files
```

---

## What `pull` writes

| Result | `obs` column | dtype |
|---|---|---|
| a highlight group | `cellucid_highlight__<group id>` (`highlight_prefix`) | `bool` |
| a user-defined categorical field (labels, annotations, clusters) | `<key>` (`field_prefix`) | `category`, the viewer's categories in order |
| a continuous field, a score (a gene signature's: `<set> score`) | `<key>` (`field_prefix`) | `float32` |
| the cells the filters leave visible (`visibility=True`) | `cellucid_visible` | `bool` |
| a gene set | -- | `uns["cellucid"]["gene_sets"][name]`: `genes`, `symbols` (from a viewer), `description`, `origin`, `pulled_from`, `id` |

- **Genes by identifier.** A gene set's `genes` are the var index, as
  `adata.var_names` holds them -- Ensembl ids for a CELLxGENE dataset -- and
  `symbols` the names the viewer shows for them (a symbol column such as
  `feature_name`). {meth}`AnnDataViewer.gene_sets` returns the identifiers,
  or the symbols with `symbols=True`. A gene set name used twice is
  `"name (2)"`, `"name (3)"`…, as everywhere in Cellucid.
- **Missing values.** A cell with no category or no value is missing (`NaN`).
- **Subsamples.** When the viewer opened a subsample, every result names
  cells by their row in `obs` (original positions). A cell the viewer did not
  open is `<NA>` in a highlight column (the nullable `boolean` dtype) and
  missing in a field's: unknown, not "not selected".
- **The cells in the viewer's order.** A viewer names a cell by its
  original position: an object's own row, or, for a file `prepare()` wrote
  (cells stored in priority order, `uns["cellucid"]["priority"]` verified
  against its cells), the original position its record gives the row. `pull` and `subset` from a
  viewer refuse an object holding the viewer's cells in another order --
  that file read with `anndata.read_zarr` as stored (put it back with
  {func}`~cellucid.restore_original_order`), or an object sorted since -- and an
  object whose cells the viewer's dataset does not have; nothing is written.
  R's `cellucid_pull()` and `cellucid_subset()` refuse the same.
- **Again and again.** `uns["cellucid"]["pulled"]` records every column a
  pull wrote: where it came from (`pulled_from`: `"viewer"`, or the file's
  name) and the result it holds (`kind`, `id`). A pull of the same result
  writes the same column again, whatever suffix it took. With `prune=True`
  (the default) a viewer pull removes the columns -- and the gene sets -- of
  an earlier viewer pull whose result the viewer no longer holds. A result
  left out of a pull (`highlights=False`, a pushed field without
  `include_pushed`) keeps its column; a file's columns are never removed.
- **Conflicts.** A column of the same name that no pull from the same source
  wrote is a conflict: `on_conflict="error"` (the default) raises and writes
  nothing, `"replace"` overwrites it, `"suffix"` writes `<name>_2` beside it
  (and the next pull writes `<name>_2` again). Session files follow the same
  rule.
- **Pushed columns.** A field this notebook pushed comes back only with
  `include_pushed=True`: `adata` holds it already.
- **One call, one write.** Every column is built first and `obs` changes
  once: a pull that fails leaves `adata` as it was. A MuData takes the
  columns on its global `obs`.

{func}`apply_cellucid_session_to_anndata` writes a saved session's highlight
groups and fields with the same column rules; `pull` of a
`.cellucid-session` file plans its columns the same way, names them with
every other column of the pull, and adds the session's notebook gene sets.

---

## The web app's exports

| Export | File | Read back as |
|---|---|---|
| Highlights › **CSV** | a `highlight` column `cellucid_highlight__<group id>` per group of the page (`true`/`false`), a row per open cell | one `bool` column per group |
| a field's **Export** | the field under its name, of kind `categorical` or `continuous`, a row per open cell | a `category` or `float` column |
| Highlights › **.h5ad** | the highlighted cells as an AnnData, in the dataset's own order (ascending original position): obs names, every obs field (continuous ones as `float32` or `float64`, as the viewer holds them), each group as a `bool` column, the genes by identifier (the var index; the shown names in `var['symbol']` when they differ), the expression rows (CSR, `float32`) of the matrix the data mapping reads, the embeddings shown (`float32` or `float64`, as stored). A key an HDF5 file cannot name -- one holding `/`, or `.` -- is written renamed (`/` as `_`), and the app says which | `anndata.read_h5ad` |
| a gene set editor's **Export**, Marker genes › GMT | GMT, genes by identifier | {func}`read_gmt`, `uns["cellucid"]["gene_sets"]` |
| Save session | `.cellucid-session` | {func}`apply_cellucid_session_to_anndata` |

An obs-column CSV declares itself, so nothing is guessed from its values:

```text
obs_name,leiden,cellucid_highlight__highlight_1,score
kind,categorical,highlight,continuous
AAACCTG-1,0,true,0.25
AAACGGG-1,10,false,
```

The header row names the columns (`obs_name` first); the kind row (`kind`
first) gives each column's kind: `highlight` (`true`, `false`, or empty for
unknown), `categorical` (labels -- `0`, `1`, `10` or `True` stay labels;
empty: no category) or `continuous` (numbers; empty: no value). A CSV
without a kind row is refused, saying how to add one. A cell the file does
not list is unknown (`<NA>`, or missing), as an unopened cell of a subsample;
a cell the file names that `adata` does not have is refused, naming it.

GMT files follow one rule in the web app, Python and R: cells are
tab-separated and trimmed; each run of tabs and line breaks inside a name,
description or gene is written as one space; a line without a name reads as
`Set <n>`; a name used again is `name (2)`, `name (3)`…; each gene is kept
once. {func}`write_gmt` writes a name that repeats an earlier one, once
cleaned, as {func}`read_gmt` will name it.

---

## Pushing into the viewer

{meth}`AnnDataViewer.push` takes a `pandas.Series` aligned by its index
(obs names: every cell, or some, the rest without a value) or an array of one
value per cell. Categorical, string, object and `bool` values become a
categorical field (a `bool` column's labels are `"false"` and `"true"`, as R
pushes a logical one); numbers a continuous one. Pushing the same key again
replaces it; a field of that key the notebook did not push is never replaced
(the viewer refuses it with a `command_error`). A pushed column is saved with the
session (`user-defined/values/<id>` holds a continuous one's values).

{meth}`AnnDataViewer.push_gene_set` adds a gene set (or replaces the one of
that name), its genes by identifier (`adata.var_names`) or by the symbols the
viewer shows (identifiers first): Gene signature and Dot plot offer it under **From notebook**, and
it is saved with the session (`gene-sets/notebook`).

---

## The wire (for developers)

Results travel on the viewer channel ({doc}`../../e_jupyter_hooks/05_architecture_viewer_channel`,
protocol schema 2, `viewer-protocol.json`):

- the `result` event (binary): a header `{type, requestId, epoch, revision,
  kind, id, nObs, meta, cells, values}` and a body of `cells` original
  positions (`uint32`, strictly ascending, none when `cells` is `null`) then
  `values` (one per named cell, or per cell of the dataset when `cells` is
  `null`; a code type's largest value is "no category", `NaN` "no value").
  Kinds: `highlight`, `categorical`, `continuous`, `geneSet` (`genes` by
  identifier, the var index, and `symbols`, the names shown, in the same
  order), `visibility` (`visible`), `subset` (`open`, the open cells);
- `result_removed {epoch, revision, kind, id}` and `results {epoch, dataset,
  revision, nObs, items}`, the inventory ending a snapshot: after a dataset
  opens (a page reload included) and whenever the channel opens again the
  viewer sends every result and the inventory, so a notebook's copy is
  replaced whole. `dataset` is the id of the notebook's dataset the results
  describe, or `null` while another one is open in the viewer (then `items`
  is empty and no result, selection, hover or click is sent);
- an `epoch` (16 hexadecimal characters, drawn at random) is one page's
  account of one dataset: a new page, and every dataset a page opens, start
  one, with revisions from 1. The notebook keeps one epoch's results; an
  event of another epoch starts afresh. Within an epoch revisions grow with
  every result and removal: a frame whose revision the notebook already holds
  is one it applied;
- `setField` (binary: a pushed column) and `setGeneSet` are answered by the
  result they made (carrying their `requestId`), or by `command_error`; so
  are `highlight`, `setVisibility`, `setField` and `setGeneSet` refused while
  the viewer shows another dataset;
- `color_by {field, source, label}`: a gene by its identifier, and the name
  the viewer shows for it.
  {attr}`AnnDataViewer.results` is the live copy; `viewer.register_hook("result", …)` runs on each.
  Nothing in the round trip is a request the notebook waits on.

---

## API reference

```{eval-rst}
.. autofunction:: pull

.. autofunction:: write_gmt

.. autofunction:: read_gmt

.. autoclass:: PullSummary
   :members:

.. autoclass:: Selection
   :members:

.. autoclass:: ViewerResult
   :members:
```

The viewer methods ({meth}`AnnDataViewer.pull`,
{attr}`AnnDataViewer.results`, {attr}`AnnDataViewer.selection`,
{attr}`AnnDataViewer.obs_names`, {meth}`AnnDataViewer.subset`,
{meth}`AnnDataViewer.gene_sets`, {meth}`AnnDataViewer.push`,
{meth}`AnnDataViewer.push_gene_set`) are documented with the class in
{doc}`viewers`.
