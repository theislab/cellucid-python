# Troubleshooting (filtering)

This page is a **symptom → diagnosis → fix** guide for filtering issues.

It follows the troubleshooting template in:
- `cellucid/markdown/DOCUMENTATION_MASTER_GUIDE.md`

---

## Quick triage checklist (do this first)

1) Look at **Active filters**:
   - What does it say: `Showing X of Y points`?
   - Which filters are listed?
2) Confirm you’re in the **right view**:
   - In Live + Snapshots mode, click the snapshot you care about first.
3) If a filter exists, try **disable-before-delete**:
   - uncheck a filter to test impact, then decide to remove/reset.
4) If you used continuous filtering:
   - check whether `Live filtering` is `Off` (then you must click `FILTER`).
5) If you used gene filtering:
   - confirm the gene is still the active field (gene filters don’t apply once you switch to obs).
6) If you used outlier filtering:
   - confirm the outlier slider is visible and set below `100%`,
   - and confirm the outlier filter is not marked `(disabled)` in Active filters.

---

## “My plot is empty (0 visible cells)”

### Symptom

- The canvas looks empty.
- Active filters shows: `Showing 0 of N points`.

### Likely causes (ordered)

1) A continuous Min/Max filter range excludes everything (common).
2) You hid categories until none remain visible.
3) You are on a snapshot with a different filter stack than you expect.
4) Outlier filtering is active and set aggressively (rare but possible).
5) A gene-range filter is active (only while the gene is active).

### How to confirm

1) In **Active filters**, disable filters one-by-one.
2) Watch the count line after each toggle:
   - when it jumps above 0, you found (at least) one culprit.
3) If you’re in grid mode, click the live view and compare counts.

### Fix

Use the fastest safe resets:

- Continuous filter: click `RESET` in the continuous legend (or remove the filter row).
- Categorical filter: click `Show All` in the categorical legend (or remove the filter row).
- Outlier filter: set outlier slider to `100%` (or remove the filter row).
- Gene filter: widen the range, or switch away from gene field if you did not intend to filter by gene.

### Prevention

- Apply filters one at a time and watch the count line.
- For large datasets, turn `Live filtering` off and apply once via `FILTER`.
- In multiview, name snapshots after the filter state you intended (so you notice when you’re in the wrong view).

---

## “I removed a filter but cells don’t come back”

### Symptom

- You click the `×` next to a filter (or reset a slider), but the visible count stays lower than expected.

### Likely causes (ordered)

1) Another filter is still active (most common).
2) You removed a filter in one view, but you’re looking at a different snapshot/view.
3) You disabled a filter but didn’t remove/reset it (so the legend state still de-emphasizes categories).
4) A gene filter is active because the gene is still selected as the active field.
5) Outlier filtering is active on the current active field.

### How to confirm

1) In **Active filters**, verify the list is actually empty.
2) If it is empty but you still see fewer points, confirm the count line says `Showing all N points`.
3) Switch to the live view (if snapshots exist) and compare the Active filters list there.

### Fix

1) Remove or reset the remaining filters (Active filters is the checklist).
2) If you’re in multiview:
   - click the view you want to fix,
   - then remove/reset filters there.
3) If categories still look “off” after removing category filters:
   - use `Show All` in the categorical legend to fully restore.

### Prevention

- Treat Active filters as the single source of truth; don’t rely on memory of what you changed.
- Prefer disabling filters while debugging, then remove once you’re sure.

---

## “Outlier slider does nothing”

### Symptom

- You move the Outlier filter slider, but the visible count doesn’t change (or outliers remain).

### Likely causes (ordered)

1) The slider is still at `100%` (no outlier filtering).
2) The active field does not provide outlier quantiles (so filtering cannot apply).
3) The outlier filter is disabled in Active filters (shows `(disabled)`).
4) Many cells have missing outlier quantiles (e.g., categories too small), so fewer cells are affected.
5) Your dataset simply doesn’t have many “extreme” points for this field.

### How to confirm

1) Confirm the outlier slider is visible (it only appears when supported by the active field).
2) Set it to something obviously aggressive like `50%` as a diagnostic.
3) Look for an Active filters entry like: `Field: outlier ≤ 50%`.
4) If it exists, confirm it is not marked `(disabled)`.

### Fix

1) Switch to a categorical field that supports outlier filtering (often clusters).
2) Set outlier to `95–99%` for typical usage (avoid `50%` unless you truly mean it).
3) If the filter is disabled, re-enable it in Active filters.
4) If you created a derived/edited categorical field and the slider disappeared:
   - switch back to the original field, or re-export with the desired field upstream.

### Prevention

- Use outlier filtering as a “clean up clusters” tool, not a general QC replacement.
- Keep in mind it is active-field scoped (it stops applying if you switch to fields without outlier stats).

---

## “My continuous sliders change, but nothing happens until I click FILTER”

### Symptom

- You drag Min/Max sliders, but the plot doesn’t update.

### Likely causes (ordered)

1) `Live filtering` is set to `Off` (expected behavior).

### How to confirm

- In the continuous legend, check the `Live filtering` toggle state.

### Fix

- Click `FILTER`, or toggle `Live filtering` back to `On`.

### Prevention

- For very large datasets, keep Live filtering `Off` and apply changes with `FILTER` to stay responsive.

---

## “Filtering is extremely slow / UI lags”

### Symptom

- Sliders feel laggy, the UI stutters, or the browser tab becomes temporarily unresponsive while filtering.

### Likely causes (ordered)

1) Large `n_cells` combined with Live filtering ON (most common).
2) Many stacked filters are enabled.
3) You are in a heavy render mode (e.g., smoke) while changing filters.
4) Many snapshots/views exist and the app is updating multiple view overlays/badges.
5) Your machine/browser GPU is resource-constrained (especially on laptops).

### How to confirm

1) Turn Live filtering `Off` and apply once with `FILTER`. If performance improves, the bottleneck is repeated recomputation.
2) Disable filters one-by-one and see which one causes the biggest delay.
3) Switch to a simpler render mode (points) and retry.

### Fix

1) Use apply-once workflows:
   - Live filtering `Off` → adjust → `FILTER`.
2) Reduce the number of enabled filters:
   - remove “dead” filters that don’t change visible/available counts.
3) Reduce visualization load while filtering:
   - use points mode while tuning filters,
   - reduce the number of snapshots temporarily.
4) If you routinely need heavy gating:
   - pre-filter in Python before export.

### Prevention

- For large datasets, treat UI filtering as interactive refinement, not the primary data-cleaning step.
- Name snapshots after key filters to avoid repeatedly debugging “why is this view slow/empty?”
