# Data mapping

**Audience:** everyone who opens a dataset in Cellucid
**Time:** 5 minutes
**What you'll learn:**
- Which stored arrays Cellucid plots, and how to choose others.
- What happens when a file has no UMAP.
- How one view can show t-SNE while another shows UMAP.
- What a session records about it.

---

## At a glance

A dataset can hold several embeddings: a UMAP, a t-SNE, a force-directed
layout, a PCA, spatial coordinates. The **data mapping** says which of them
Cellucid draws in each dimension (1D, 2D, 3D), and, read-only, which neighbour
graph, gene names, expression matrix, latent space and vector fields it reads.

It lives in **Session › Dataset › Data mapping**, a collapsible section under
the dataset facts.

| What you opened | What you can change |
|---|---|
| An `.h5ad` file or a Zarr ZIP, in the browser | Any numeric `obsm` array, for any dimension, at any columns |
| A prepared dataset (sample, folder, GitHub, self-hosted) | Switch between the embeddings it was exported with |
| A file served by `cellucid serve` or a notebook | Nothing in the panel; start the server with `--embedding-1d/-2d/-3d` |

The neighbour graph, gene names, expression matrix, latent space and vector
fields are listed under **Fixed**. They are not changed in the panel.

---

## Choose an embedding

Open **Session › Dataset › Data mapping**. Under **Embeddings**, each dimension
has one menu:

- For a file read in the browser, the menu lists every `obsm` entry wide enough
  for that dimension, as `X_pca · 50 columns`. Entries that cannot be drawn are
  listed disabled with the reason: `obsm_df · obsm table, not an array`,
  `X_sparse · sparse obsm arrays cannot be plotted`. Vector fields
  (`<field>_umap_<d>d`) are not offered as embeddings.
- **Not used** removes a dimension. It is accepted only when no view shows that
  dimension and another dimension remains.
- For a prepared dataset, the menu lists its embeddings of that dimension, as
  `UMAP (source file) · X_umap`. With only one embedding the menu is disabled.

Choices are staged: nothing loads until you press **Apply**, and **Revert**
drops them. Apply loads the new embedding and checks every value first. Then
the dimension's default changes, and every view that follows the default is
redrawn with it. A view that shows its own embedding (see
[One embedding per view](#one-embedding-per-view)) keeps it.

Cells, highlights, selections, filters and colours are kept: they belong to
cells, not to coordinates. Every embedding is scaled to the same box, so the
camera stays where it was. Centroid labels, level of detail and picking are
rebuilt for the new layout.

## Columns of a wide array

An array wider than the dimension, such as a 50-column `X_pca`, shows one
column menu per axis, **x**, **y** and (in 3D) **z**. Columns are counted from
0, so `X_pca` columns 0 and 1 are PC1 and PC2. Each axis needs a different
column; choosing the same column twice shows "Each axis needs a different
column." and Apply stays disabled.

The embedding is then known as `X_pca[2,3]` (the key and its columns) wherever
a layout is named: in a view's badge, in the per-view menu and in a session.

## Files without a UMAP: the best guess

A file opened in the browser normally shows the UMAP it carries:
`obsm['X_umap_1d']`, `obsm['X_umap_2d']` and `obsm['X_umap_3d']` fill their
dimensions, or, when none of those exist, the plain `X_umap` that
`sc.tl.umap()` writes fills the dimension of its width. A present UMAP is
always used, even when it cannot be drawn: a NaN in `X_umap` is an error, not
a reason to show something else.

When the file has no UMAP at all, Cellucid guesses instead of refusing it:

1. It considers every `obsm` array that can be drawn, except vector fields.
2. It ranks their names: `umap`, then `tsne`, then `draw_graph` or `force`,
   then `phate`, `pacmap` or `trimap`, `spatial`, `diffmap`, any other name,
   and last `pca`, `scvi`, `latent`, `harmony` or `X_emb`.
3. Each dimension takes the best-ranked array exactly that wide.
4. When no array is 1 to 3 columns wide, only 2D is filled, from the first two
   columns of the best-ranked wide array (for example `X_pca`).
5. An array that fails the check below is skipped for the next one.

The guess is the same for the same file every time, in the browser and in
`cellucid serve`. It is marked **guessed** in the panel, and a notification says
what is shown: "This file has no default UMAP (X_umap_2d, X_umap_3d or X_umap).
Showing 'X_tsne' as 2D. Change it under Session › Data mapping." Its **Show**
button opens the panel. On a desktop the panel also opens by itself; on a phone
it is only marked open, so the plot is not covered as the file opens.

When nothing in `obsm` can be drawn, the file is refused with a list of what
`obsm` holds and why each entry was not used:

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell
(10,000) and at least one numeric column. obsm holds: X_df (table),
spatial (10,000 × 2, NaN at cell 17).
```

## One embedding per view

With several embeddings in one dimension, **Compare views › Dimension** shows
an **Embedding** menu for the selected view. Its first entry is the dataset
default (`UMAP (dataset default)`); the others are every embedding of that
dimension.

To compare t-SNE and UMAP side by side:

1. Press **Keep view** to freeze the current view.
2. Select the kept view, open **Embedding** and choose `X_tsne`.

The view's badge then reads `2D · X_tsne`. Cells and highlight pages are the
same in every view, so a selection made in one panel marks the same cells in
the other; each view keeps its own field, filters, camera, centroid labels and
vector field. Changing the view's dimension returns it to the dataset default.
Figure export draws each view with the embedding it shows.

A vector field is bound to the embedding it was computed on. In a view that
shows another embedding, the overlay is unavailable and says which embeddings
it is bound to.

## What stays fixed

**Fixed** lists what else the dataset is read with:

| Row | Prepared dataset | File read in the browser | Served file |
|---|---|---|---|
| Neighbour graph | fixed by `prepare()` (not listed) | `obsp['connectivities']` when present | `--knn` |
| Gene names | the name set chosen by `prepare()` | the var index | `--gene-id-column` |
| Expression | fixed by `prepare()` (not listed) | `X` | `--expression` |
| Latent space | fixed by `prepare()` (not listed) | none | `--latent-key` |
| Vector fields | the bindings `prepare()` wrote | `<field>_umap_<d>d`, bound to `X_umap_<d>d` or `X_umap` | `--vector-fields` |

## Sessions

A session file records the data mapping in its dataset fingerprint, and each
view's own embedding with that view. Loading the session restores both:

- If the session names embeddings this dataset can show, they are restored
  before the views, whatever the panel shows now.
- If it names an embedding this dataset does not have, the session is refused
  and names it: "This session uses obsm['X_tsne'] columns [0, 1] in 2D, which
  this dataset does not have."
- If it was saved with a different neighbour graph, gene names, expression
  matrix, latent space or vector fields, it is refused and names the group.

See {doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

## When an embedding is refused

Apply, a per-view choice and a session all check an embedding's values before
it is drawn. Every value must be a finite number that fits in 32-bit floating
point, and the columns must not all be constant. A refusal names the first
offending cell:

| Message | Cause |
|---|---|
| `NaN at cell 17` | a missing value in row 17 |
| `Infinity at cell 3` / `-Infinity at cell 3` | an infinite value |
| `not representable as float32 at cell 5` | a value too large (or too small, but not zero) for 32-bit floats |
| `no coordinate variation` | every cell has the same coordinates |

The previous mapping stays in place, the row shows the reason, and the array is
listed disabled with that reason until the dataset is opened again. Because a
single cell has no coordinate variation, a one-cell file cannot be opened.

---

## Next steps

- {doc}`03_browser_file_picker_tutorial` — open an `.h5ad` or Zarr ZIP.
- {doc}`08_troubleshooting_data_loading` — when a file does not open.
- {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples` — multiview.
