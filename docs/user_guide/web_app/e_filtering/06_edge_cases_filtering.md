# Edge cases (filtering)

**Audience:** computational users + power users  
**Time:** 10–20 minutes  
**What you’ll learn:**
- The most common “surprising” filtering states and how to recognize them
- How missing values affect visibility (NaN and small categories)
- When some filters apply only conditionally (active gene / outlier-supporting fields)
- How multiview changes what “my filters” means

**Prerequisites:**
- A dataset loaded

---

## 1) All cells filtered out (0 visible points)

### What you see

- The canvas looks empty.
- Active filters shows: `Showing 0 of N points`.

### Why it happens

One or more filters reduced visibility to zero:
- a continuous Min/Max range excludes everything,
- categories were hidden until none remain,
- outlier threshold is very low (and outlier filtering is active),
- you are looking at a snapshot with a different filter stack than you expect.

### What to do

Use a deterministic “unwind” sequence:

1) In **Active filters**, disable filters one-by-one.
2) When you find the culprit, either:
   - relax it (widen range / show categories / raise outlier percent), or
   - remove it (×), or
   - reset it in the legend (`RESET`, `Show All`, outlier `100%`).

---

## 2) NaN / missing values (continuous filters)

### What you see

- Applying a continuous filter removes more cells than expected.
- Even with Min≈min and Max≈max, some cells never come back while the filter is active.

### Why it happens

Continuous range filters treat `NaN` (and `null`) as “outside range”.

So those cells become invisible whenever the filter is active.

### What to do

- If you want those cells to remain visible, fill or impute missing values upstream and re-export.
- If you want to exclude them, current behavior already does this.

---

## 3) “Some categories are disabled / I can’t click them”

### What you see

- In the categorical legend, some category rows are disabled.
- Tooltip: “No cells available in this category after other filters.”

### Why it happens

Cellucid computes:
- **available** counts = after other filters (ignoring this field’s category visibility),
- **visible** counts = after all filters (including this field).

If upstream filters remove all cells in a category, that category becomes unavailable in the current view.

### What to do

1) Relax upstream filters (continuous ranges, other category hides, gene/outlier filters if active).
2) Revisit the categorical legend.
3) Use `Show All` to make categories visible again.

---

## 4) “I disabled a category filter, but those categories look gray”

### What you see

- You disable a category filter in Active filters.
- Cells return, but some categories are colored gray/de-emphasized.

### Why it happens

Cellucid tracks:
- which categories are “visible” in the legend, and
- whether that visibility is currently enforced as a filter.

Disabling the category filter stops it from hiding cells, but the legend may still treat those categories as “hidden” for display purposes (so they appear gray).

### What to do

- Use `Show All` in the categorical legend (or remove the filter) to fully restore.

---

## 5) Outlier slider disappears after editing a categorical field

### What you see

- Outlier filtering was available on a categorical field.
- After merging/deleting categories (creating a derived field), the outlier slider is gone.

### Why it happens

Derived categorical fields typically do **not** carry latent outlier statistics, because those stats are tied to the original category definitions at export time.

### What to do

- Switch back to the original (non-derived) field if you need outlier filtering, or
- re-export data with the desired categorical field defined upstream.

---

## 6) “My gene filter (or outlier filter) stops applying when I switch fields”

### What you see

- You filter by gene expression or outliers.
- You switch to another field.
- The filtered-out cells reappear.

### Why it happens

Current behavior:
- gene filters apply only while that gene is the active field,
- outlier filtering applies only while the active field supports outlier stats (that’s when the control is visible).

### What to do

- Treat gene/outlier filtering as “active-field scoped” tools.
- If you need persistent gating while coloring by something else, export the gate as an obs field (e.g., a score/boolean) and filter it as obs.

---

## 7) Multiview: “filters look different across snapshots”

### What you see

- In Live + Snapshots mode, different views show different visible counts.
- Active filters changes when you click different snapshots.

### Why it happens

Each view/snapshot can have its own filter state.

### What to do

- Click the view you care about first, then inspect Active filters.
- If you want two views to match, apply the same filter stack to both (or recreate the snapshot from the filtered live view).

---

## 8) Quantization and slider resolution (export-time edge case)

### What you see

- Sliders feel “chunky” (small changes don’t affect visibility).
- Very fine thresholds are hard to set.

### Why it happens

Two effects can compound:

- UI sliders are `0–100` percent steps (1% increments).
- Exports can quantize continuous fields and/or outlier quantiles (8-bit/16-bit), reducing numeric precision.

### What to do

- Use UI filtering for coarse/interactive gating.
- Do fine-grained gating upstream in Python before export when precision matters.

---

## Next steps

- `07_troubleshooting_filtering` (symptom-based debugging)
