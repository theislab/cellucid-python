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
**Session file**. Its **Embeddings** list is a checklist; each further group
has its own disclosure.

Every dataset is read as AnnData (an `.h5ad` or a Zarr store, or any other
format Cellucid reads as AnnData: MuData, Loom, 10x, an R object, a
MatrixMarket folder; opened in the browser, from a link or the samples, or
served from Python or R), and for every one
you can turn any numeric `obsm` array on or off, at any 1 to 3 of its columns,
as many column selections as you like; and choose any `obsp` graph, `X`, a
layer or `raw.X` named by the var index or any var column, any dense `obsm`
array as the latent space, and any `obsm` array that fits an embedding as its
vector field.

A served file is read by the same reader in the browser as a file you open
yourself: the server only hands over its bytes. Everything this page says about
a file read in the browser holds for it; only the choices it opens with can
differ (see [Launch defaults](#launch-defaults)).

Turning an embedding on loads nothing. Its coordinates are read the first time
a view shows it. After Embeddings, one disclosure per further group follows:
[Genes and expression](#genes-and-expression),
[Neighbour graph](#neighbour-graph), [Latent space](#latent-space) and
[Vector fields](#vector-fields). Every group is staged until **Apply**, and
applied with the embeddings as one step.

---

## Choose the embeddings views can show

Open **Session › Data mapping › Embeddings**. The rows are one flat list,
and each row is one embedding:

- the checkbox turns it on (views can show it) or off;
- its name (`UMAP`, `t-SNE`, `PCA 1,2,3`; a long one takes a
  second line, and hovering it shows the array and columns) and its dimension
  (`3D`). The dimension is a button with a chevron: it opens the row's editor, described under
  [Dimension and columns](#dimension-and-columns-of-an-array);
- for an embedding that is on, a second line: **Default** marks the one
  embedding per dimension that a view switched to that dimension shows first,
  and that keys 1, 2 and 3 show first (see
  {doc}`../c_core_interactions/05_dimension_switching_1d_2d_3d`), and
  **guessed** marks a default the [best guess](#files-without-a-umap-the-best-guess)
  chose. An embedding that is off is one line.

What is listed:

- Every `obsm` array except vector fields (a key `<field>_<basis>` beside an
  `X_<basis>` in the same file, such as `velocity_umap` beside `X_umap`; see
  [Vector fields](#vector-fields)). An array 1 to 3 columns wide is an
  embedding and is offered at all of its columns. A wider one (`X_pca`, a
  latent space thousands of columns wide, a 4-column `coordinates`) is a matrix
  of components, not an embedding: it is listed after the embeddings and
  offered at its first two columns, in 2D, the way such a matrix is plotted
  (PC1 against PC2). Either can be changed (see below). Entries that cannot be
  drawn are listed with the reason, and cannot be turned on: `obsm table, not
  an array`, `sparse obsm arrays cannot be plotted`.
- Arrays that are off and wider than three columns, and entries that cannot be
  drawn, are listed last, under a collapsed **Other arrays (N)** disclosure; a
  wide one there shows its width (`3,136 columns`). An array stays where it
  was listed while you edit, even after you turn it on.
When a file opens, the defaults are on, and so is every other array that is 1
to 3 columns wide; a wider array is on only when it is a default.

```{figure} ../../../_static/screenshots/data_loading/data-mapping-other-arrays.png
:alt: The Data mapping group of the Garcia sample with its Embeddings checklist open as one flat list of rows, each with its checkbox, name, dimension and Default mark, and at its end Other arrays (1) open, listing scVI 1,2, which is off, with its width, 24 columns.
:width: 488px

The Garcia sample's Embeddings: one flat list, each row with its checkbox,
name, dimension and **Default** mark, and **Other arrays** open at the end,
where the 24-column `scVI` array waits, off, with its width.
```

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
cameras stay where they were. A dataset you open next starts from its own
view: the default camera of the dimension it opens in, and the point size its
cell count calls for, whatever the dataset before it showed; its session or
published preset, when it has one, then sets its own.

## Dimension and columns of an array

Every `obsm` array two or more columns wide can be shown in any dimension its
columns allow, at any of its columns: a 2-column `X_umap` in 1D or 2D, a
3-column one in 1D, 2D or 3D, a 50-column `X_pca` in 1D, 2D or 3D. The
dimension a file opens on is only where it starts.

Click a row's dimension (`2D ⌄`) to open its editor; click it again to close
it. Editors start closed; one you opened stays open while you edit and after Apply.
The editor holds:

- **Dimension**, a **1D/2D/3D** menu offering the dimensions the array's
  columns allow;
- a column menu per axis, **x**, **y** and (in 3D) **z**. Columns are shown
  counted from 1 (the column 0 of `obsm['X_pca']` is shown as 1, PC1);
- for an array wider than three columns, **Add columns**, which adds another
  selection of the same array, for example components 3 and 4, as its own row,
  on and with its editor open.

A selection other than all of an array's columns in order is named by its
columns: `UMAP 1` is the first column of a 2-column UMAP in 1D, `Spatial 1,3`
the first and third columns of a 3-column array in 2D, `PCA 3,4` components 3
and 4. Each selection is one embedding views can show. Each axis needs a
different column, and two rows cannot select the same columns; Apply stays
disabled until that is fixed. A selection you added can be removed with its ×.

The embedding is recorded as `X_pca` with its 0-based columns wherever it is
named exactly: in a view's hover text, in a session, on the Jupyter wire.

## Files without a UMAP: the best guess

A file opened in the browser normally shows the UMAP it carries:
`obsm['X_umap_1d']`, `obsm['X_umap_2d']` and `obsm['X_umap_3d']` are the
defaults of their dimensions, or, when none of those exist, one UMAP under its
conventional name is the default of the dimension of its width: the first 1- to
3-column array of `X_umap` (what `sc.tl.umap()` writes, and a Seurat object's
`umap` reduction), `X_UMAP` (a SingleCellExperiment's `UMAP` reducedDim, as
scater's `runUMAP()` names it), `X_wnn.umap` and `X_ref.umap` (a Seurat
object's weighted-nearest-neighbour and Azimuth reference UMAPs). An R
object's reductions and reducedDims keep their names, case and all, as
`obsm['X_<name>']`. A present UMAP is always used, even when it cannot be
drawn: a NaN in `X_umap` is an error, not a reason to show something else.

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

A view opens in the dimension whose default ranks best by the same names, the
higher dimension when they rank the same. A segmentation export with three
2-column `X_umap_scale*_latents`, a 3-column `X_tsne_scale0_latents` and a
3-column `spatial` gets the t-SNE as its 3D default, and opens on
`X_umap_scale0_latents` in 2D. A file whose defaults are all UMAPs
(`X_umap_1d` to `X_umap_3d`, or the defaults `prepare()` names) opens in its
highest dimension. Keys 1, 2 and 3 switch to the others.

The guess is the same for the same file every time, in the browser and in
`cellucid serve`. It is marked **guessed** in the checklist, and a notification
says what is shown: "This file has no default UMAP (X_umap_1d, X_umap_2d,
X_umap_3d, X_umap, X_UMAP, X_wnn.umap or X_ref.umap). Showing 'X_tsne' as 2D. Change it under Session › Data mapping." Its
**Show** button opens the checklist. Nothing opens on its own: every dataset
starts with **Data mapping** collapsed, so the plot is never covered as the
file opens.

```{figure} ../../../_static/screenshots/data_loading/best-guess-notification.png
:alt: A Cellucid notification reading "This file has no default UMAP (X_umap_1d, X_umap_2d, X_umap_3d, X_umap, X_UMAP, X_wnn.umap or X_ref.umap). Showing 'X_pca' as 2D. Change it under Session › Data mapping.", with a Show button under the text and a dismiss button at its right.
:width: 760px

An `.h5ad` whose only `obsm` array is `X_pca`: the notification names every
key the default rule looks for, what is shown instead, and where to change it.
```

```{figure} ../../../_static/screenshots/data_loading/best-guess-show-checklist.png
:alt: The Data mapping group after Show on the best-guess notification: Embeddings open and summarised as 1 on · PCA 2D · guessed, its one row, PCA in 2D, on and marked Default and guessed, above the collapsed Genes, Neighbour graph, Latent space and Vector fields groups and the Revert and Apply buttons.
:width: 496px

**Show** opens **Data mapping › Embeddings** on the guess: `PCA` in 2D, on,
marked **Default** and **guessed**.
```

When nothing in `obsm` can be drawn, the file is refused with a list of what
`obsm` holds and why each entry was not used:

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell
(10,000) and at least one numeric column. obsm holds: X_df (table),
spatial (10,000 × 2, NaN at cell 17). Add a 2-D embedding to obsm (for example
with sc.tl.umap in Scanpy), then open the data again.
```

## Launch defaults

What an AnnData file opens with — each dimension's embedding, the neighbour
graph, the expression matrix and gene names, the latent space, the vector field
the overlay opens on, the obs columns shown, embedding normalization and the
fewest cells a category needs for outlier quantiles — is decided per choice,
in this order:

1. The launch option, when the file is served and the option was set:
   `cellucid serve --embedding-2d`, `--knn`, `--expression`,
   `--gene-id-column`, `--latent-key`, `--vector-field-default`, `--obs-key`,
   `--modality`, or the matching `serve_anndata()` / `show_anndata()` keyword
   (`embedding_keys`, `knn_key`, `expression_source`, `gene_id_column`,
   `latent_key`, `vector_field_default`, `obs_keys`, `modality`,
   `normalize_embeddings`, `centroid_min_points`). The server checks each
   name against the file's metadata before it starts, and passes only the
   options that were set.
2. The file's own `uns['cellucid']['defaults']`, when it declares them.
3. The rule this page describes for that choice (the UMAP, else the best guess;
   `obsp['connectivities']`; `X`; the var index, or its column of symbols
   when the index holds identifiers; no latent space; each
   `<field>_<basis>` vector field bound to its `X_<basis>`; every obs column).

Every choice can still be changed under **Data mapping** once the file is open.

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
read from, what the genes are called, whether each matrix holds `log1p`
values, and which genes are highly variable.

**Modality** comes first, for a dataset of several modalities only: a
MuData's modalities, or the feature types of a 10x matrix (Gene Expression,
Antibody Capture). It chooses whose features the rest of the disclosure lists;
Expression, Gene names, Values and Highly variable genes are that modality's,
and changing it sets them to its defaults. A dataset opens on its
gene-expression modality. A cell the modality did not measure has no value
(NaN), not zero. See {doc}`19_multimodal_mudata`.

**Expression** lists `X`, every layer as `layers['counts']`, and `raw.X` with
its gene count, as `raw.X · 30,000 genes`, because `raw` has its own, often
longer, gene list. **None** shows no genes. A file opens on `X`. Without `X`
it opens on `raw.X`; without either, on its only layer, or, among several, on
one whose name contains `norm` or `log`, then any other, with `counts` and
`raw_counts` last (alphabetically within each). The line under the menus then
says "Guessed: expression (no X)." A served file opens on
`cellucid serve --expression` (`show_anndata(..., expression_source=...)`)
when it names one (see [Launch defaults](#launch-defaults)).

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

A file opens named by its var index, unless the index holds database
identifiers while a var column holds their symbols, as every CELLxGENE file
does (Ensembl ids in the index, symbols in the categorical `feature_name`).
Then that column names the genes, marked **Guessed**, and each gene keeps its
id: the gene search finds `SFTPC` and `ENSG00000168484` alike, and lists the
id beside the symbol. The choice is made by content, not by name alone, for
each matrix's genes on their own (`X` and its layers, `raw.X`, each
modality):

- the index holds identifiers when at least 90% of an even sample of its
  names look like one: letters then at least five digits, optionally
  versioned (`ENSG00000168484`, `ENSMUSG00000064370.1`, `FBgn0000008`,
  `WBGene00000001`), a RefSeq (`NM_000546`) or UniProt (`P04637`) accession,
  `HGNC:5`, or a bare number (an Entrez id). A few symbols look like one too
  (`LINC01128`, `Gm12345`), which is why it takes 90%;
- only then are var columns read, one at a time: `feature_name`,
  `gene_symbols`, `gene_symbol`, `symbols`, `symbol`, `gene_name`,
  `gene_names`, `gene_short_name`, `feature_names`, `genes`, `gene`, `names`,
  `name` (any letter case) first, the first that holds symbols chosen; else
  every other text column, the one most like symbols chosen;
- a column holds symbols when at least half of the genes have a value, at
  least half of those are distinct (a biotype or a species column is not), at
  most half look like identifiers, at least 90% are one word starting with a
  letter and not a genomic interval (`chr1:1000-1050`), and it is not the
  index again. Integer columns are never symbols;
- a column that cannot be read is passed over (choosing it says why).

The first view never waits for gene names. A dataset is drawn first, and the
gene names of the matrix it opens with are read right after: its var index
and, when that holds identifiers, the columns above. Only that matrix's axis
is guessed then; another axis (`raw.X`, another modality) is guessed when you
first choose it under **Expression** or **Modality**. While names are read,
**Gene names** says **Reading gene names…**, the gene search says the same
and cannot be typed in, and analyses that need genes open once they are
read; a session, an export with genes and a change of the data mapping wait
for them. The var index is never used in their place meanwhile. A var index
that cannot be read leaves the genes unnamed and says why, in the gene search
and in **Gene names**, where another name source can still be chosen.

Whatever names them in the viewer, genes leave it by their identifier (the
var index name, as in `adata.var_names`): the gene the viewer colours by and
the gene sets a notebook receives, and gene-set GMT files. A gene set pushed
from a notebook may name genes by identifier or by the name shown
(identifiers first; a name two genes answer is reported, not guessed).

A file that names its genes itself (`uns['cellucid']['defaults']['geneIdColumn']`,
`--gene-id-column`) opens on that column for the matrix it opens with; the
other matrices and modalities keep their own default. Switching **Expression**
to a matrix on another var axis (to or from `raw.X`) takes that axis's
default. When you choose the var index and it holds identifiers while the
dataset found a column of symbols, the line says so:
"Gene names are identifiers; 'gene_symbols' holds symbols."

**Values** says how the expression matrix chosen above is stored: **As stored
(not log-transformed)**, **log1p (base e)**, **log1p (base 2)**, **log1p
(base 10)**, and any other base the file records. Each matrix keeps its own
choice: switch **Expression** to `layers['counts']` and **Values** shows that
layer's. Analyses read it: fold changes are computed on the values with the
`log1p` undone (`expm1`, or `base^x - 1`), as Scanpy's `rank_genes_groups`
does. The default follows Scanpy, which records `sc.pp.log1p` in
`uns['log1p']`:

- with `uns['log1p']`, `X` and every layer are `log1p`, in the base it
  records (`base`; none, as Scanpy writes `base=None`, is e). That is how
  Scanpy reads the record: `rank_genes_groups` applies it to whichever matrix
  it is given. A `counts` layer is guessed `log1p` too; set it to **As stored**;
- `raw.X` is what `adata.raw = adata` saved, conventionally the
  log-normalised matrix before scaling: it is `log1p` too, unless it is stored
  as integers (counts);
- without `uns['log1p']`, every matrix is as stored.

A default read from `uns['log1p']` is marked guessed (the line under the
menus says "Guessed: log1p values.") until you change it.

**Highly variable genes** lists **None** and every true/false column of the
chosen matrix's var table (`raw.var` for `raw.X`): `highly_variable` first,
then names that say so (`highly_variable_*`, `hvg`, `*_hvg`), then any other.
It starts on `highly_variable`, the column `sc.pp.highly_variable_genes`
writes and Scanpy uses as its `mask_var`, when that column is true/false
(the line under the menus says "Guessed: highly variable genes."); otherwise
on **None**. Choosing it reads nothing: analyses that keep to the highly
variable genes read the column when they first ask for it. A missing value in
a nullable column is not marked.

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
- Gene search, the gene menus and the **Genes** row of the **Dataset** box's
  **Details** use the new list; any
  other gene loads when it is first used.
- Analysis windows close, since what they computed came from the previous
  values.

## Latent space

The **Latent space** disclosure (**Outliers in**) chooses where a
categorical field's outlier filter measures distance: each cell's distance to
its category's centre in that space, as the share of its category at that
distance or closer (see
{doc}`../e_filtering/02_outlier_filtering_per_active_field`). It is never the
2D/3D layout, and Cellucid never picks one for you. Its menu lists **None**
and every dense `obsm` array that is not a vector field, with its width, as
`X_scVI · 10 columns`; the latent space the file stores outlier quantiles for
is marked `outliers stored`, and filters at once. An array that is sparse or
a table is listed disabled with the reason.

The dataset starts on the file's own `uns['cellucid']['defaults']['latentKey']`
(written by `prepare(latent_key=...)` or `cellucid_prepare(latent_key = ...)`),
or on the launch option (`cellucid serve --latent-key`,
`serve_anndata(..., latent_key=...)`, `show_anndata(..., latent_key=...)`);
otherwise on **None**, and the outlier control is one **Choose latent
space** button, with no slider, that opens this setting
([Turning it on](../e_filtering/02_outlier_filtering_per_active_field.md#turning-it-on-choose-a-latent-space)).

**Apply** reads nothing. Every outlier threshold is kept — "hide the farthest
5%" means the same in any space — the previous space's quantiles are dropped,
and a field whose outlier filter is in use gets its quantiles in the new
space at once, with a progress bar and **Cancel**. Any other field gets them
on its slider's first move: read from the file when it stores them for this
space, otherwise computed in the browser's reader worker from the latent
space, which is read once and kept in memory within Cellucid's budget
([Stored or computed](../e_filtering/02_outlier_filtering_per_active_field.md#stored-or-computed)).
The browser computes exactly the bytes `prepare()` stores. A latent space of
any size works: it is streamed from the file in row slabs. A space holding a
missing value (NaN) or a value beyond float32 fails at its first read with a
message naming the cell, and is marked unusable. **None** turns the outlier
filters off until a latent space is chosen again.

## Vector fields

The **Vector fields** disclosure chooses which `obsm` arrays are drawn as
vectors (velocity, drift) on which embedding. It lists every embedding that is
turned on, each with a checkbox per array that fits it:

- an array as wide as the embedding's own array, read at the same columns
  (`velocity_umap` on `X_umap`; `velocity_pca` on `X_pca[2,3]` reads columns 2
  and 3);
- or an array exactly as wide as the embedding's dimension, read whole (a
  2-column `drift` on `X_pca[0,1]`).

Arrays named for the embedding come first, marked **suggested**. The name
follows scVelo's convention, the vector-field key rule: an `obsm` key
`<field>_<basis>` is the field `<field>` drawn over `X_<basis>`, for any basis
(`velocity_umap` over `X_umap`, `velocity_umap_3d` over `X_umap_3d`,
`velocity_branches` over `X_branches`, `T_fwd_umap` over `X_umap`). A key that
starts with `X_` names an embedding, never a vector, and is not offered. When
the key splits at more than one `_` into a field and a basis whose `X_<basis>`
the file holds, the longest basis wins: with both `X_umap_3d` and `X_3d`,
`velocity_umap_3d` is `velocity` over `X_umap_3d`. Matching is exact and
case-sensitive.

A file starts with the default rule: each `<field>_<basis>` is bound to every
enabled view of `X_<basis>` when both are usable arrays with the same number of
columns, read at that view's columns (`velocity_pca`, 50 columns, over
`X_pca`, 50 columns, shown at columns `[0, 1]`). Every field of one basis is
bound (`velocity_umap` and `T_fwd_umap` both on `X_umap`). A vector named for
an embedding of another width, or one that cannot be read, is still not
offered as an embedding, but stays unbound until you bind it here. A
`velocity_umap_2d` beside a bare `X_umap`, with no `X_umap_2d`, names no
embedding in the file, so it is no vector field at all but an ordinary `obsm`
array (listed among the embeddings): rename it `velocity_umap`, as scVelo
writes it, or bind it here. A key with no such split is likewise an ordinary
array.

The overlay lists a field bound to the embedding it is named for under its
`<field>` part, so the arrays of one field on every embedding are one field:
`velocity_umap`, `velocity_umap_3d` and `velocity_branches` are all the field
`velocity`. An array you bind to another embedding is listed under its own
`obsm` key. `uns['cellucid']['defaults']['vectorFieldDefault']` or
`cellucid serve --vector-field-default`
(`show_anndata(..., vector_field_default=...)`) names the field the overlay
opens on, by that field id (`velocity`).

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

A file starts with `obsp['connectivities']` when it has one; a served file
starts with the graph `cellucid serve --knn KEY` (or
`show_anndata(..., knn_key=...)`) names instead. Either way the graph is read
only when edges are drawn or the KNN drag tool is used.

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
drag, and the **Neighbour graph** row of the **Dataset** box's **Details** counts
its edges.

## Where each group starts

Every group can be changed; this is what each one shows when the dataset
opens:

| Group | File read in the browser | Served file |
|---|---|---|
| Neighbour graph | `connectivities` | `--knn`, else as in the browser |
| Modality (several only) | the first gene-expression modality (`rna`, `gex`, `gene expression`, `gene_expression`, `sct`), else the first | as in the browser |
| Gene names | the var index, or the var column of symbols when the index holds identifiers (guessed) | `--gene-id-column`, else as in the browser |
| Expression | `X` | `--expression`, else as in the browser |
| Values (log1p) | `X`, every layer and a non-integer `raw.X` are `log1p` when `uns['log1p']` exists (its `base`, e when none); else as stored | as in the browser |
| Highly variable genes | `var['highly_variable']` (`raw.var` for `raw.X`) when it is true/false, else None | as in the browser |
| Latent space | the file's `latentKey`, else None | `--latent-key`, else as in the browser |
| Vector fields | each `<field>_<basis>` bound to `X_<basis>` when both have the same number of columns | the same rule |

A file that declares `uns['cellucid']['defaults']` starts on those (`knn`,
`modality`, `geneIdColumn`, `expression`, `latentKey`, `vectorFieldDefault`, and the
embeddings of each dimension); a choice a server names on its command line
takes precedence over the file's.

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
- If it was saved with another neighbour graph, modality, gene names, expression
  matrix, values (`log1p`) of a matrix, highly variable genes column, latent
  space or vector fields, that choice is restored when this
  dataset has it and can use it; otherwise the session is refused and names
  it: "This session uses obsp['scvi_connectivities'] as its neighbour graph,
  which this dataset does not have." A gene a view was coloured by is
  restored on the saved gene list, and the gene count is checked against that
  list (`raw.X` has its own). Outlier thresholds are restored with the
  latent space they were set in, by getting their quantiles (stored or
  computed, with progress and **Cancel**); a notification names each threshold
  that could not be restored.

See {doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

(embedding-refused)=
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
file cannot be opened, and a subsample opens at least two cells
({ref}`fewest-cells-open`).

## Memory

An embedding that is on but shown by no view costs nothing until a view shows
it. Once shown, its coordinates stay in memory while a view shows it, and
afterwards in a cache bounded by a sixteenth of the device memory (1 GiB on a
16 GiB machine; see {ref}`memory-budgets`), so switching back is immediate. The least recently shown is released
first, and everything no view shows is released when memory runs short or the
page has been hidden for ten minutes; it is read again when a view shows it.

---

## Next steps

- {doc}`03_browser_file_picker_tutorial` — open an `.h5ad`, a Zarr folder or a Zarr ZIP.
- {doc}`08_troubleshooting_data_loading` — when a file does not open.
- {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples` — multiview.
- {doc}`../../python_package/g_api_reference_coverage/api/viewers` — read and change the mapping from a notebook; `cellucid serve`'s mapping flags are in {doc}`../../python_package/g_api_reference_coverage/api/cli`.
