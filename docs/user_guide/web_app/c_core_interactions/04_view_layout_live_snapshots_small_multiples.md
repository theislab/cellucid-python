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

1. In the live view, choose the dimension, navigation mode, coloring, category
   transparency, and filters you want to preserve.
2. Frame the cells of interest.
3. Click `Keep view`.
4. Leave `Locked Cam` on for the same-framing comparison.
5. Select the live or kept badge, choose `Edit selected view`, and change that
   panel's field or filters.
6. Return to `Grid compare`.

After step 3, the new snapshot is active and the layout is `Grid compare`.
Repeating `Keep view` at the three-snapshot limit fails without replacing an
existing panel; remove a snapshot before trying again.

## What a kept snapshot preserves

A snapshot captures the selected view's:

- active categorical or continuous field and its rendered colors;
- category transparency and filter-derived visibility;
- continuous/outlier filter state and filter summary;
- dimension and dimension-specific positions;
- camera state and navigation mode;
- centroid render state associated with that view.

These are independent view contexts after capture. Later coloring or filtering
changes in one selected view do not propagate to the other panels. Highlight
pages are a separate feature and should not be treated as part of the
per-snapshot contract.

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
  available;
- a clickable `Orb`, `Pan`, or `Fly` badge when cameras are unlocked; it cycles
  **Orbit → Planar → Free-fly**;
- `×` when that view can be removed.

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: Cellucid showing two side-by-side kept views of the same dataset.
:width: 1440px

Two panels compare one dataset while retaining independent coloring and filter
contexts.
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
- every snapshot's field, filters, dimension, label, and navigation state.

Restoration recreates the snapshot graph rather than storing GPU buffers. The
dataset itself is never embedded in the session. See
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored` and
{doc}`../l_sessions_sharing/07_versioning_compatibility_and_dataset_identity`.

Figure Export can render the current multiview grid; verify panel order and
legends in {doc}`../k_figure_export/02_export_ui_walkthrough`.

## Failure and recovery

### `Keep view` is disabled

Set **Visualization → Render mode:** to `Points`. If the button is enabled but
snapshot creation reports `Keep view failed`, preserve the current live view,
remove an unneeded snapshot if three already exist, and retry.

### A change affected the wrong panel

1. Stop changing controls.
2. Click the intended badge.
3. Choose `Edit selected view`.
4. Confirm its field, dimension, and filter summary before continuing.

### The wrong panel moves

With `Unlocked Cam`, click inside the intended panel first. Use `Locked Cam` if
all panels should move together.

### Performance drops

Every visible panel adds rendering work, and vector-field trails allocate
per-view buffers. Hide live when at least two snapshots exist, remove unused
snapshots, or return to a single panel. See
{doc}`../n_benchmarking_performance/03_large_dataset_best_practices`.
