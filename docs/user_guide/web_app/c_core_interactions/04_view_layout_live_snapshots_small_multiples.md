# Multiview: Live view and kept snapshots

**Audience:** everyone
**Time:** 10 minutes

Multiview makes controlled side-by-side comparisons of one loaded dataset.
Cellucid supports the live view plus at most **three** kept views—up to four
visible panels—and each view can show any embedding the dataset has turned on.

## A useful comparison

For a differentiation dataset, keep one view colored by cell type, then change
the selected view to a cell-cycle score or another measured covariate. Views of
one embedding can share one camera (**Link cameras**), so the same region stays
under the cursor in both panels. This helps answer "is this apparent branch
associated with biology or with a technical/continuous covariate?" without
relying on memory.

Views can also show different embeddings of the same cells: UMAP beside t-SNE,
two 3D embeddings side by side, a 3D view beside a 2D one. A highlight or a
selection made in one view marks the same cells in every view.

Multiview compares render state; it does not perform a statistical comparison.
Use {doc}`../h_analysis/index` for quantitative follow-up.

## Exact controls

Open **Compare views**. The **Views** group lists the open views, and the
controls under the list edit the view you select. The small `i` buttons say
what each control does.

| Control | Current behavior |
|---|---|
| Views list | One row per open view: its number, what it shows (`UMAP · 3D`), and under that the field it is coloured by. Select a row to edit that view. |
| `View 2 shows` | The embedding the selected view shows, chosen from the embeddings turned on under **Session › Data mapping**, grouped by dimension. With a single view it reads `Shows`. |
| `Keep view` | Adds a view: a copy of the selected one, with its own colours, filters and embedding. The new view is selected. |
| `Link cameras` | Views that show the same embedding share one camera and navigation mode. Available only when two open views show the same embedding; see [Linked cameras](#linked-cameras). |
| `Layout: Grid compare` | Shows the visible views together. |
| `Layout: Edit selected view` | Shows only the selected view at full canvas size for editing. |
| `Clear` | Removes every kept view and returns to the live view. |

`Layout` and `Clear` appear once a view has been kept. `Keep view` is available
only in `Render mode: Points`. Smoke disables it and forces the live,
single-view layout; see {doc}`03_render_modes_points_vs_volumetric_smoke`.

You can keep any number of views. In a grid, two or three views share one
row and four or more take a grid of ⌈√n⌉ columns. GPU memory, not a count, is
the limit: views that are not on screen give up their colours and visibility
on the GPU when memory runs short and are rebuilt, unchanged, when shown
again. Only a set of views shown at once that does not fit in the device's
memory is refused.

## Create a reliable comparison

1. In the live view, choose the embedding, navigation mode, coloring, category
   transparency, and filters you want to preserve.
2. Frame the cells of interest.
3. Click `Keep view`. View 2 is added and selected.
4. Under `View 2 shows`, choose the embedding it should show, or change its
   field or filters. Each change applies to view 2 only.
5. Select row 1 to go back to editing the live view.

After step 3 the layout is `Grid compare`.

## What a kept view preserves

A kept view starts as a copy of the selected view's:

- active categorical or continuous field and its rendered colors;
- category transparency and filter-derived visibility;
- continuous/outlier filter state and filter summary;
- embedding (see {doc}`../b_data_loading/13_data_mapping`);
- camera and navigation mode;
- centroid labels.

From then on these belong to that view: coloring, filtering or switching the
embedding of one view does not change the others. Highlight pages are a
separate feature: they belong to cells, not to views, and every view shows
them.

## Where the panels sit

`Grid compare` lays the panels out in the part of the window the sidebar does
not cover. With the sidebar open, the panels start at its right edge; every
panel has the same size, so the left panel is never partly hidden under the
sidebar. Closing the sidebar (the `☰`/`✕` tab or `H`), dragging its edge to
resize it, resizing the window, or moving the window to a screen with another
pixel density re-lays the panels out at once. While the sidebar slides open or
closed, the panels keep the area beside the open sidebar and take the full
width when the closing slide has finished, so no frame shows a covered panel.
On a phone the menu is a full-screen sheet that you close to reach the plot,
and the panels then use the whole screen.

Everything you do on the plot follows the same layout: a click or tap focuses
the panel under it, picking, hover, every selection tool, projectiles, the
panel title chips, and the selection bar all use the panel rectangles you see.

A single view (one panel, or `Edit selected view`) still draws across the whole
window, but the camera centre - the orbit or pan target, or the free-fly look
direction - sits at the centre of the uncovered part, so the data is centred
where you can see it.

## Select and edit the intended panel

Select a row of the Views list (click it, or focus it with `Tab` and press
`Enter`) to make that view active. Field selectors, filters, `View N shows`
and **Navigation** then describe the selected view; the label of the
embedding picker names the view it edits. Clicking inside a panel in
`Grid compare` selects it too. A drag or a wheel zoom that starts in another
panel moves that panel's camera at once; the sidebar switches to the panel
when the gesture ends (on release, or once the wheel has rested), so the
gesture never waits for the sidebar to redraw.

`Grid compare` is best for inspection. `Edit selected view` is best for
changing one panel, because it removes ambiguity about which view owns the
sidebar changes. Returning to the grid does not merge view state.

## Each view's camera and navigation

Every view has its own camera and navigation mode. A view showing a 3D
embedding orbits; a 1D or 2D view pans (**Planar**). Clicking or starting a
drag inside a panel focuses it: the drag moves that panel, and the
**Navigation** panel (with its **Applies to** row) and the `O`, `P` and `G`
keys edit its mode. A 3D view beside a 2D view keeps orbiting while the 2D
view pans. See {doc}`01_navigation_modes_orbit_planar_free_fly`.

## Linked cameras

**Link cameras** is on by default. While it is on, views that show the same
embedding share one camera and one navigation mode: turning, panning or zooming
one moves them all, and their rows are marked `linked` (hover the mark to see
which views it links). Views that show different embeddings always keep their
own cameras, even while it is on: one camera cannot frame a 3D and a 2D layout,
or two embeddings of different shape.

So the button is available only when two open views show the same embedding.
Otherwise it is disabled and shows as off, and its tooltip says why ("No two
views show the same embedding"). The setting itself is kept: as soon as two
views show the same embedding again, they are linked.

When you link cameras, each group of views of one embedding takes one camera:
the focused view's for its group, the first view's (the live view first) for
the others. When you unlink them, every view keeps what it shows.

## Read and operate the view rows

Each row contains:

- a number;
- a `⌖` mark on the view whose camera has focus;
- what the view shows, `UMAP · 3D`, and under it the field it is coloured by
  (hover for the exact embedding key);
- `linked` when the view shares its camera with views of the same embedding;
- its navigation mode, `Orbit`, `Pan` or `Fly`; press it to switch the mode
  (Orbit → Planar → Free-fly) for that view;
- `×` when that view can be removed.

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: The whole window in Grid compare with the sidebar open: two equal panels fill the area to the right of the sidebar and show the Pancreas 3D embedding from one shared camera — the left panel, the live view, coloured by cell_type, and the right panel, a kept view, recoloured by S_score on a Viridis scale — each under a title chip naming its field. The sidebar, scrolled to Compare Views in an earlier interface, lists badge 1 cell_type and the selected badge 2 S_score above the Layout select set to Grid compare, and Coloring & Filtering shows S_score in the Continuous obs select.
:width: 1440px

Grid compare with the sidebar open (earlier interface: the rows now also name
the embedding each view shows). The live view (cell_type) and a kept view
recoloured by S_score (selected) sit in two equal panels to the right of the
sidebar, framed by one linked camera, and each panel keeps its own coloring.
```

Removing `×` from a kept view deletes that view. When at least two kept views
exist, removing the live row hides the live panel and selects a kept view.
With only one kept view, removing live instead removes that kept view and
returns to the live view, so the viewer never reaches an invalid zero-view
state.

## Session and export behavior

A matching-dataset `.cellucid-session` saves and restores:

- `Grid compare` versus `Edit selected view`;
- active view and whether live is hidden;
- whether **Link cameras** is on, and every view's own camera and navigation
  mode;
- every view's embedding (as the embedding key and columns, so it is restored
  onto the same array), field and filters;
- whether each view shows centroid points and centroid labels.

A view's label and filter summary are recomputed from its field and filters.
Restoration recreates each kept view once rather than storing GPU buffers. The
dataset itself is never embedded in the session. See
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored` and
{doc}`../l_sessions_sharing/07_versioning_compatibility_and_dataset_identity`.

Figure Export can render the current multiview grid (**Export all views**);
each panel is drawn with the embedding it shows, its own projection and its own
centroid labels. Verify panel order and legends in
{doc}`../k_figure_export/02_export_ui_walkthrough`.

## Failure and recovery

### `Keep view` is disabled

Set **Visualization → Render mode:** to `Points`. With **Sample dataset** set
to `None` there are no cells to keep: load a dataset first (the button's
tooltip names the reason). If the button is enabled but the new view reports
`Keep view failed` (the device ran out of memory for the views shown at
once), preserve the current live view, show fewer views at once (Layout:
single) or remove an unneeded kept view, and retry.

### A change affected the wrong panel

1. Stop changing controls.
2. Select the intended row.
3. Choose `Edit selected view`.
4. Confirm its embedding (`View N shows`), field and filter summary before
   continuing.

### The wrong panel moves

Click (or, on a touch screen, start dragging) inside the intended panel first:
a camera press focuses the panel it lands on. A selection gesture (`Alt`, or a
finger with **Select mode** on) never moves focus. Views of the same embedding
move together while **Link cameras** is on; views of different embeddings
never do.

### Performance drops

Every visible panel adds rendering work, and vector-field trails allocate
per-view buffers. Views of one embedding share its position buffer and spatial
tree, so a second view of an embedding costs no extra point memory. Hide live
when at least two kept views exist, remove unused views, or return to a single
panel. See
{doc}`../n_benchmarking_performance/03_large_dataset_best_practices`.
