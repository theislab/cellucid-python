# Data mapping

**Audience:** everyone who opens a dataset in Cellucid
**Time:** 5 minutes
**What you'll learn:**
- Which stored arrays Cellucid can show, and how to turn more on.
- How several embeddings of one dimension (UMAP, t-SNE, PCA) are shown side by side.
- What happens when a file has no UMAP.
- How to choose the neighbour graph, the gene names and the expression matrix,
  the latent space of outlier filters, and the vector fields.
- What a session records about it.

---

## At a glance

A dataset can hold several embeddings: a UMAP, a t-SNE, a force-directed
layout, a PCA, spatial coordinates, in 1D, 2D or 3D. The **data mapping** says
which of them views can show — any number, of any mix of dimensions — and
which one is each dimension's **default**. It also names the neighbour graph,
gene names and expression matrix Cellucid reads, the latent space outlier
filters measure in, and the vector fields it draws on which embedding.

It is the last group of the **Session** panel, **Data mapping**, below
**Session file**. Its **Embeddings** list is a checklist; **Read from file**
lists the rest.

| What you opened | What you can change |
|---|---|
| An `.h5ad` file or a Zarr ZIP, read in the browser or served from Python (a notebook or the command-line server) | Turn any numeric `obsm` array on or off, at any 1 to 3 of its columns, as many column selections as you like; choose any `obsp` graph, `X`, a layer or `raw.X` named by the var index or any var column, any dense `obsm` array as the latent space, and any `obsm` array that fits an embedding as its vector field |
| A prepared dataset (sample, folder, GitHub, self-hosted) | Turn the embeddings it was exported with on or off, choose each dimension's default, and choose among the gene-name sets and vector fields it was exported with |

A served file opens with the server's defaults (`--embedding-1d/-2d/-3d`, or
its best guess) and every other narrow `obsm` array turned on; the server reads
any array the checklist turns on when a view first shows it, and so for every
other choice.

Turning an embedding on loads nothing. Its coordinates are read the first time
a view shows it. After Embeddings, one disclosure per further group follows:
[Genes and expression](#genes-and-expression),
[Neighbour graph](#neighbour-graph), [Latent space](#latent-space) and
[Vector fields](#vector-fields). Every group is staged until **Apply**, and
applied with the embeddings as one step.

---

## Choose the embeddings views can show

Open **Session › Data mapping › Embeddings**. Each row is one embedding:

- the checkbox turns it on (views can show it) or off;
- its name and dimension, `UMAP 3D`, `PCA (unintegrated) 2D`;
- **Default** marks the one embedding per dimension that a view switched to
  that dimension shows first, and that keys 1, 2 and 3 show first (see
  {doc}`../c_core_interactions/05_dimension_switching_1d_2d_3d`).

What is listed:

- For a file read in the browser, every `obsm` array except vector fields
  (`<field>_umap_<d>d`). An array 1 to 3 columns wide is offered at all of its
  columns; a wider one at its first columns (see below). Entries that cannot be
  drawn are listed with the reason, and cannot be turned on: `obsm table, not
  an array`, `sparse obsm arrays cannot be plotted`.
- For a prepared dataset, every embedding it was exported with.

When a file opens, every embedding of a prepared dataset is on. For a file
read in the browser, the defaults are on, and so is every other array that is
1 to 3 columns wide; a wider array is on only when it is a default.

Changes are staged: nothing loads until you press **Apply**, and **Revert**
drops them. Turning the default of a dimension off makes the first other
embedding of that dimension that is on, as listed, its default; turning on the
first embedding of a dimension makes it the default. At least one embedding stays on, and an embedding a vector field is
drawn on cannot be turned off.

On Apply, a view that shows an embedding you turned off switches to its
dimension's default (or, with none left, to the dataset's opening embedding).
A view that follows its dimension's default — one you did not switch to a
particular embedding — follows a new default. Every other view keeps what it
shows. The embeddings those views switch to are loaded and checked first; if
one is refused, nothing changes.

Cells, highlights, selections, filters and colours are kept: they belong to
cells, not to coordinates. Every embedding is scaled to the same box, so
cameras stay where they were.

## Columns of a wide array

An array wider than three columns, such as a 50-column `X_pca`, is listed at
its first columns with a column menu per axis, **x**, **y** and (in 3D) **z**,
and a **1D/2D/3D** menu that sets how many axes it has. Columns are shown
counted from 1 (the column 0 of `obsm['X_pca']` is shown as 1, PC1).

**+ Columns of PCA** adds another selection of the same array, for example
components 3 and 4, as its own row; each selection is one embedding views can
show, named by its components (`PCA 3,4`). Each axis needs a different column,
and two rows cannot select the same columns; Apply stays disabled until that is
fixed. A selection you added can be removed with its ×.

The embedding is recorded as `X_pca` with its 0-based columns wherever it is
named exactly: in a view's hover text, in a session, on the Jupyter wire.

## Files without a UMAP: the best guess

A file opened in the browser normally shows the UMAP it carries:
`obsm['X_umap_1d']`, `obsm['X_umap_2d']` and `obsm['X_umap_3d']` are the
defaults of their dimensions, or, when none of those exist, the plain `X_umap`
that `sc.tl.umap()` writes is the default of the dimension of its width. A
present UMAP is always used, even when it cannot be drawn: a NaN in `X_umap` is
an error, not a reason to show something else.

When the file has no UMAP at all, Cellucid guesses the defaults instead of
refusing it:

1. It considers every `obsm` array that can be drawn, except vector fields.
2. It ranks their names: `umap`, then `tsne`, then `draw_graph` or `force`,
   then `phate`, `pacmap` or `trimap`, `spatial`, `diffmap`, any other name,
   and last `pca`, `scvi`, `latent`, `harmony` or `X_emb`.
3. Each dimension's default is the best-ranked array exactly that wide.
4. When no array is 1 to 3 columns wide, only 2D gets a default, from the first
   two columns of the best-ranked wide array (for example `X_pca`).
5. An array that fails the check below is skipped for the next one, and stays
   off.

The guess is the same for the same file every time, in the browser and in
`cellucid serve`. It is marked **guessed** in the checklist, and a notification
says what is shown: "This file has no default UMAP (X_umap_2d, X_umap_3d or
X_umap). Showing 'X_tsne' as 2D. Change it under Session › Data mapping." Its
**Show** button opens the checklist. On a desktop the checklist also opens by
itself; on a phone it is only marked open, so the plot is not covered as the
file opens.

When nothing in `obsm` can be drawn, the file is refused with a list of what
`obsm` holds and why each entry was not used:

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell
(10,000) and at least one numeric column. obsm holds: X_df (table),
spatial (10,000 × 2, NaN at cell 17).
```

## Several embeddings at once

Every view shows one of the embeddings that are on, of any dimension: two 3D
embeddings side by side, a 3D one beside a 2D one, three 2D ones. To compare
t-SNE and UMAP:

1. Open **Compare views** and press **Keep view**: a second view, a copy of the
   selected one, is added and selected.
2. In **View 2 shows**, choose `t-SNE · 2D`.

Each row of **Views** says what its view shows (`t-SNE · 2D`) and, under it,
the field it is coloured by. Select a row to edit that view: the menu under the
list is labelled with the view it edits.

Cells and highlight pages are the same in every view, so a highlight or
selection made in one view marks the same cells in all of them, whatever each
shows. Each view keeps its own embedding, colouring, filters, camera,
navigation mode, centroid labels and vector field. Views that show the same
embedding share its coordinates and GPU buffers, so a second view of an
embedding costs no extra memory for positions.

A vector field is bound to the embedding it was computed on. It is drawn in
every view that shows that embedding, and is unavailable in views that show
another.

## Genes and expression

The **Genes** disclosure chooses the gene list: which matrix gene values are
read from, and what the genes are called.

**Expression** lists `X`, every layer as `layers['counts']`, and `raw.X` with
its gene count, as `raw.X · 30,000 genes`, because `raw` has its own, often
longer, gene list. **None** shows no genes. A file opens on `X`. Without `X`
it opens on `raw.X`; without either, on its only layer, or, among several, on
one whose name contains `norm` or `log`, then any other, with `counts` and
`raw_counts` last (alphabetically within each). The line under the menus then
says "Expression is guessed: this file has no X." A served file opens on
`cellucid serve --expression` (`show_anndata(..., expression_source=...)`)
when it names one, and on the same rule otherwise.

**Gene names** lists the var index and every column of the var table the
chosen matrix uses (`raw.var` for `raw.X`). Text, categorical and integer
columns can name genes; integers are written as text (`7`, not `7.0`). A
floating-point or true/false column is listed disabled with the reason, as
`dispersions · floating-point values cannot name genes`. A gene whose column
value is empty or missing takes its var-index name, and repeated names are
made unique as AnnData's `var_names_make_unique()` does: the second `MATR3`
becomes `MATR3-1`. The line under the menus says what that changed:

```text
1 gene has no 'gene_symbols' value and uses its var index name.
1 duplicate name was made unique, e.g. 'LYZ' → 'LYZ-1'.
```

When the names are Ensembl ids (at least 90% of a sample look like
`ENSG00000139618`) and the var table has a symbol column, the line says so:
"Gene names are Ensembl ids; 'gene_symbols' holds symbols."

Apply reads the new gene list, and the values of every gene a view is
coloured by, before anything changes. Then:

- A view coloured by a gene keeps that gene, with its colour range, log scale
  and filter. Other names for the same matrix (the var index and a symbol
  column, or `X` and a layer) keep the same column: `ENSG00000090382` becomes
  `LYZ`. Between `raw.X` and the others the gene is found by name. A range or
  filter outside the new values is narrowed to them.
- A gene missing from the new list loses its colouring, and Apply says so:
  "Applied. 'PPBP' is not on the new gene list; its colouring was removed."
- Cells, highlights, selections, obs filters and colours are kept.
- Gene search, the gene menus and **Dataset › Genes** use the new list; any
  other gene loads when it is first used.
- Analysis windows close, since what they computed came from the previous
  values.

A prepared dataset offers the gene-name sets `prepare()` exported, such as
`symbol` and `gene_id`; its expression matrix is fixed.

## Latent space

The **Latent space** disclosure chooses where a categorical field's outlier
filter measures distance: each cell's distance to its category's centre in
that space, as the share of its category at that distance or closer (see
{doc}`../e_filtering/02_outlier_filtering_per_active_field`). Its menu lists
**None** and every dense `obsm` array with its width, as
`X_scVI · 10 columns`; an array that is sparse, a table, or larger than
512 MiB as 32-bit floats is listed disabled with the reason. A file opened in
the browser starts with **None**; a served file starts with
`cellucid serve --latent-key` (`show_anndata(..., latent_key=...)`), or None.

Apply reads the space, checks every value (finite, representable as a 32-bit
float) and computes the quantiles of every categorical field a view has
loaded before anything changes; a space with a bad value is refused in the
server's words ("NaN at cell 5") and stays unavailable. Then each field's
outlier filter uses the new quantiles in every view, and its threshold is
kept: "hide the farthest 5%" means the same in any space. **None** removes the
outlier filters. The browser computes exactly what `prepare()` and
`cellucid serve` compute, off the main thread.

## Vector fields

The **Vector fields** disclosure chooses which `obsm` arrays are drawn as
vectors (velocity, drift) on which embedding. It lists every embedding that is
turned on, each with a checkbox per array that fits it:

- an array exactly as wide as the embedding's dimension, read whole
  (`velocity_umap_2d` on `X_umap`);
- or an array as wide as the embedding's own array, read at the same columns
  (`velocity_pca` on `X_pca[2,3]` reads columns 2 and 3).

Arrays named for the embedding come first, marked **suggested**: the default
rule's `<field>_umap_<d>d` for a UMAP, and `<name>_<basis>` for `X_<basis>`
(`velocity_umap`, as scVelo writes it). `X_…` arrays are embeddings and are
not offered. A file opened in the browser starts with the default rule:
`<field>_umap_<d>d` is bound to the UMAP of that dimension. A served file
starts with that rule when `cellucid serve --vector-fields`
(`show_anndata(..., serve_vector_fields=True)`) asks for it, and with none
otherwise. A prepared dataset lists the bindings it was exported with, and any
of them can be turned off.

An embedding a field is bound to stays on in the checklist until the field is
unbound. Apply replaces the bindings and turns the overlay off; nothing is read until
the overlay shows a field. A view draws a field only while it shows an
embedding the field is bound to (see
{doc}`../i_vector_field_velocity/02_enabling_overlay_and_selecting_field`).
An array whose values cannot be drawn is reported when it is first shown, and
is then listed disabled with the reason.

## Neighbour graph

For an AnnData file, the **Neighbour graph** disclosure chooses the graph that
connectivity edges and the KNN drag tool use. Its menu lists **None** and
every `obsp` matrix, as `connectivities · sparse 10,000 × 10,000`: the
`connectivities` that `sc.pp.neighbors()` writes, a `<key>_connectivities`
from `sc.pp.neighbors(key_added=…)`, and any other square matrix with one row
per cell. A matrix that is not square, or not as large as the dataset, is
listed disabled with the reason.

A file opened in the browser starts with `obsp['connectivities']` when it has
one. A served file starts with the graph `cellucid serve --knn KEY` (or
`show_anndata(..., knn_key=...)`) named, or none.

Apply reads the graph and checks it before anything changes. It must be
symmetric, with a zero diagonal and no negative, infinite or missing weight.
A graph that is not, such as the `distances` that `sc.pp.neighbors()` also
writes, is refused, never repaired:

```text
obsp['distances'] is not symmetric — choose a connectivities matrix.
```

The previous graph stays, and the refused one is listed disabled with that
reason until the dataset is opened again. After a graph changes, shown
connectivity edges are redrawn from it, the KNN drag tool uses it on its next
drag, and **Dataset › Connectivity** counts its edges.

## What stays fixed

**Read from file** lists what a dataset is read with but cannot change.
Every other group has its own disclosure:

| Group | Prepared dataset | File read in the browser | Served file |
|---|---|---|---|
| Neighbour graph | fixed by `prepare()` (not listed) | menu; starts on `connectivities` | menu; starts on `--knn` |
| Gene names | menu of the exported name sets (fixed with one) | menu; starts on the var index | menu; starts on `--gene-id-column` |
| Expression | fixed by `prepare()` (not listed) | menu; starts on `X` | menu; starts on `--expression` |
| Latent space | fixed by `prepare()` (not listed) | menu; starts on None | menu; starts on `--latent-key` |
| Vector fields | checkboxes of the exported bindings (fixed with none) | checkboxes; `<field>_umap_<d>d` bound to `X_umap_<d>d` or `X_umap` | checkboxes; the same rule with `--vector-fields`, else none |

## Sessions

A session file records the data mapping — the embeddings that are on and each
dimension's default — in its dataset fingerprint, and with each view the
embedding it shows. Loading the session restores both:

- If the session's embeddings are ones this dataset can show, they are turned
  on before the views are restored, whatever the checklist shows now, and the
  session's checklist is applied exactly once the views are in place.
- If it names an embedding this dataset does not have, the session is refused
  and names it: "This session shows obsm['X_tsne'] columns [0, 1], which this
  dataset does not have."
- If it was saved with another neighbour graph, gene names, expression
  matrix, latent space or vector fields, that choice is restored when this
  dataset has it and can use it; otherwise the session is refused and names
  it: "This session uses obsp['scvi_connectivities'] as its neighbour graph,
  which this dataset does not have." A gene a view was coloured by is
  restored on the saved gene list, and the gene count is checked against that
  list (`raw.X` has its own). Outlier thresholds are restored with the
  latent space they were set in.

See {doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

## When an embedding is refused

Apply, a view's choice and a session all check an embedding's values before it
is drawn. Every value must be a finite number that fits in 32-bit floating
point, and the columns must not all be constant. A refusal names the first
offending cell:

| Message | Cause |
|---|---|
| `NaN at cell 17` | a missing value in row 17 |
| `Infinity at cell 3` / `-Infinity at cell 3` | an infinite value |
| `not representable as float32 at cell 5` | a value too large (or too small, but not zero) for 32-bit floats |
| `no coordinate variation` | every cell has the same coordinates |

The previous mapping stays in place, the status line shows the reason, and the
array is listed with that reason, and cannot be turned on, until the dataset is
opened again. Because a single cell has no coordinate variation, a one-cell
file cannot be opened.

## Memory

An embedding that is on but shown by no view costs nothing until a view shows
it. Once shown, its coordinates stay in memory while a view shows it, and
afterwards in a cache bounded by a quarter of the device memory the browser
reports (between 512 MiB and 2 GiB; 512 MiB where the browser does not report
it), so switching back is immediate. The least recently shown is released
first, and everything no view shows is released when memory runs short or the
page has been hidden for ten minutes; it is read again when a view shows it.

---

## Next steps

- {doc}`03_browser_file_picker_tutorial` — open an `.h5ad` or Zarr ZIP.
- {doc}`08_troubleshooting_data_loading` — when a file does not open.
- {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples` — multiview.
