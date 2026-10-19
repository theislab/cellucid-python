# Multiview: Live view and kept snapshots

**Audience:** everyone
**Time:** 10 minutes

Multiview makes controlled side-by-side comparisons of one loaded dataset.
Cellucid supports the live view plus at most **three** kept snapshots—up to four
visible panels.

## A useful comparison

For a differentiation dataset, keep one view colored by cell type, then change
the selected view to a cell-cycle score or another measured covariate. With
cameras locked, the same region stays under the cursor in both panels. This
helps answer “is this apparent branch associated with biology or with a
technical/continuous covariate?” without relying on memory.

Multiview compares render state; it does not perform a statistical comparison.
Use {doc}`../h_analysis/index` for quantitative follow-up.

## Exact controls

Open **Compare Views → Multiview:**. The small `i` button states the important
scope rule: snapshots keep independent coloring and filters; select a panel and
use `Edit selected view` before changing only that panel.

| Control | Current behavior |
|---|---|
| `Keep view` | Freezes the selected view as a new panel and selects that new snapshot. |
| `Locked Cam` | Cameras are linked. This is the default once multiple views exist. |
| `Unlocked Cam` | Each view keeps an independent camera/navigation mode. |
| `Layout: Grid compare` | Shows the visible live/snapshot panels together. |
| `Layout: Edit selected view` | Shows only the selected panel at full canvas size for editing. |
| `Clear` | Removes every kept snapshot and returns to the live view. |

`Keep view` is available only in `Render mode: Points`. Smoke disables it and
forces the live, single-view layout; see
{doc}`03_render_modes_points_vs_volumetric_smoke`.

## Create a reliable comparison

1. In the live view, choose the dimension (and, where it has several, the
   embedding), navigation mode, coloring, category transparency, and filters
   you want to preserve.
2. Frame the cells of interest.
3. Click `Keep view`.
4. Leave `Locked Cam` on for the same-framing comparison.
5. Select the live or kept badge, choose `Edit selected view`, and change that
   panel's field, filters or **Embedding** (for example t-SNE in one panel and
   UMAP in the other).
6. Return to `Grid compare`.

After step 3, the new snapshot is active and the layout is `Grid compare`.
Repeating `Keep view` at the three-snapshot limit fails without replacing an
existing panel; remove a snapshot before trying again.

## What a kept snapshot preserves

A snapshot captures the selected view's:

- active categorical or continuous field and its rendered colors;
- category transparency and filter-derived visibility;
- continuous/outlier filter state and filter summary;
- dimension, the embedding it shows (the dataset default or its own, see
  {doc}`../b_data_loading/13_data_mapping`) and that embedding's positions;
- camera state and navigation mode;
- centroid render state associated with that view.

These are independent view contexts after capture. Later coloring or filtering
changes in one selected view do not propagate to the other panels. Highlight
pages are a separate feature and should not be treated as part of the
per-snapshot contract.

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

Click a view badge to make that view active. Field selectors, filters,
**Dimension:**, and **Navigation:** then describe the selected view.

`Grid compare` is best for inspection. `Edit selected view` is best for
changing one panel, because it removes ambiguity about which view owns the
sidebar changes. Returning to the grid does not merge view state.

## Camera behavior

Cellucid starts multiview with linked cameras:

- `Locked Cam`: dragging, zooming, and navigation mode are shared.
- `Unlocked Cam`: each panel owns its camera and navigation mode. Click inside
  a panel before navigating it.

Unlock cameras when comparing different framing or when deliberately combining,
for example, a 2D **Planar** view with a 3D **Orbit** view. Keep them locked when
the scientific comparison requires identical framing.

## Read and operate view badges

Each badge contains:

- a number pill;
- a `⌖` indicator for the currently focused camera when cameras are unlocked;
- a view label;
- a clickable `1D`, `2D`, or `3D` badge when another dataset dimension is
  available; a view showing its own embedding names it there, e.g.
  `2D · X_tsne`;
- a clickable `Orb`, `Pan`, or `Fly` badge when cameras are unlocked; it cycles
  **Orbit → Planar → Free-fly**;
- `×` when that view can be removed.

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: The whole window in Grid compare with the sidebar open: two equal panels fill the area to the right of the sidebar and show the Pancreas 3D embedding from the same locked camera — the left panel, the live view, coloured by cell_type, and the right panel, the kept snapshot, recoloured by S_score on a Viridis scale — each under a title chip naming its field. The sidebar, scrolled to Compare Views, lists badge 1 cell_type and the selected badge 2 S_score above the Layout select set to Grid compare, and Coloring & Filtering shows S_score in the Continuous obs select.
:width: 1440px

Grid compare with the sidebar open. The live view (badge 1, cell_type) and a
kept snapshot recoloured by S_score (badge 2, selected) sit in two equal panels
to the right of the sidebar, framed by one locked camera, and each panel keeps
its own coloring.
```

Removing `×` from a snapshot deletes that snapshot. When at least two snapshots
exist, removing the live badge hides the live panel and selects a snapshot.
With only one snapshot, removing live instead removes that snapshot and returns
to the live view, so the viewer never reaches an invalid zero-view state.

## Session and export behavior

A matching-dataset `.cellucid-session` saves and restores:

- `Grid compare` versus `Edit selected view`;
- active view and whether live is hidden;
- camera-lock state plus live/per-snapshot cameras;
- every snapshot's field, filters, dimension, embedding, label, and navigation
  state.

Restoration recreates the snapshot graph rather than storing GPU buffers. The
dataset itself is never embedded in the session. See
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored` and
{doc}`../l_sessions_sharing/07_versioning_compatibility_and_dataset_identity`.

Figure Export can render the current multiview grid; verify panel order and
legends in {doc}`../k_figure_export/02_export_ui_walkthrough`.

## Failure and recovery

### `Keep view` is disabled

Set **Visualization → Render mode:** to `Points`. With **Sample dataset** set
to `None` there are no cells to keep: load a dataset first (the button's
tooltip names the reason). If the button is enabled but snapshot creation
reports `Keep view failed`, preserve the current live view, remove an unneeded
snapshot if three already exist, and retry.

### A change affected the wrong panel

1. Stop changing controls.
2. Click the intended badge.
3. Choose `Edit selected view`.
4. Confirm its field, dimension, and filter summary before continuing.

### The wrong panel moves

With `Unlocked Cam`, click (or, on a touch screen, start dragging) inside the
intended panel first: a camera press focuses the panel it lands on. A selection
gesture (`Alt`, or a finger with **Select mode** on) never moves focus. Use
`Locked Cam` if all panels should move together.

### Performance drops

Every visible panel adds rendering work, and vector-field trails allocate
per-view buffers. Hide live when at least two snapshots exist, remove unused
snapshots, or return to a single panel. See
{doc}`../n_benchmarking_performance/03_large_dataset_best_practices`.
