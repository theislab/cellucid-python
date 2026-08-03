# Color-by behavior

**Audience:** everyone (users interpreting colors) + computational users (exact semantics)  
**Time:** 20–30 minutes  
**What you’ll learn:**
- What “color by” actually does in Cellucid (and what it does *not* do)
- How **categorical** vs **continuous** fields map to colors
- What the neutral gray “None” color means (missing / invalid values)
- How legend interactions create **real filters** (visibility changes)
- How log scale, rescaling, and outlier filtering affect what you see

---

## First principle: coloring is data → colors, but the legend can also filter

Cellucid uses the active field to compute an RGBA color for every cell (point).

However, many legend interactions change **visibility**, not just color:
- unchecking a category hides those cells,
- narrowing a numeric range hides cells outside the range,
- enabling outlier filtering hides outliers (when supported).

This is why the UI groups it as **Coloring & Filtering**.

---

## Categorical color-by (discrete labels)

### What is colored, exactly?

For a categorical obs field, each cell has a category code:
- If the code is valid, the cell is colored with that category’s RGB color.
- If the code is missing/invalid, the cell uses a neutral gray “None” color.

Category colors come from a fixed categorical palette by default, and can be overridden per category.

A missing code is not a quiet variant of a category — it is outside the
category system entirely, and that has two visible consequences:

- Those cells are drawn in the neutral grey at **full opacity**, they get no
  legend row, and they are counted in no per-category count. They still count
  toward `Showing X of Y points`.
- **No categorical filter can hide them.** Unchecking every row, or clicking
  `Hide All`, leaves them on screen. This is the usual explanation for “I hid
  every category and grey points are still there”. See
  {doc}`../e_filtering/01_filtering_mental_model` for how that fits the rest of
  the filter stack.

### What does “hiding a category” do?

Unchecking a category in the legend is not a “dim”:
- the category becomes **fully transparent** (invisible),
- and the cells are treated as **filtered out** for downstream counts and summaries.

This means you can:
- filter by one categorical field (e.g., hide `doublets`)
- while coloring by another field (e.g., clusters)

### Category order (important for interpretation)

The legend list is sorted alphabetically (human-friendly), but:
- colors are attached to category **identity** (the category index), not the legend row order.

So “the first visible row” is not a special color; it’s just the first in alphabetical sort order.

```{figure} ../../../_static/screenshots/web_app/legend-categorical-clusters.png
:alt: A legend headed Click a swatch to pick a color, then Show All and Hide All buttons, then eight rows in alphabetical order — Alpha, Beta, Delta, Ductal, Epsilon, Ngn3 high EP, Ngn3 low EP, Pre-endocrine — each with a checkbox, a coloured swatch, a label, a cell count, a pencil icon and a bin icon.
:width: 428px

The `clusters` legend on the Pancreas sample. The rows run in alphabetical
order — `Alpha`, `Beta`, `Delta`, `Ductal`, … — which is *not* the order the
categories were defined in, and the palette is attached to the category, not to
the row position.
```

---

## Continuous color-by (numeric values)

Continuous fields (continuous obs and gene expression) map numbers → colors using:
- a **colormap** (palette), and
- a **color domain** (min/max used for color mapping).

### Neutral “None” gray (missing values)

Continuous values use the neutral gray “None” color when:
- the value is `NaN` *and* the field also holds at least one finite value, or
- log scale is enabled and the value is ≤ 0 (log is undefined there).

`NaN` is a drawable value. Cellucid reads it as “not measured here” and gives it
the neutral grey, which makes missing measurements obvious without forcing you
to filter them out.

An **infinity is not drawable, and it is not greyed**: it has no position on a
colour scale, and it would stretch the field’s range so far that every other
cell collapses onto one colour. A field or gene containing one is refused
instead — the active field does not change, the picture stays as it was, and you
get a counted message naming what is wrong.

:::{important}
**Two continuous payloads are refused rather than coloured.**

- **Any `+Infinity` or `-Infinity`.** The message reads
  `Field "X" contains 12 infinite values, so it has no colour scale.`, then the
  count of each kind out of the total cells, then the first few affected cell
  indices, then the one line of Python that repairs the values.
- **Every value `NaN`.** There is no range to scale, so the field is refused
  too, with `Field "X" has no value in any of 18,142,044 cells, so it has no
  colour scale.` This is usually a column that was never computed, or an empty
  expression layer.

For a gene these arrive as a multi-line notification titled
`Gene cannot be shown`. For a continuous obs field the same text arrives after
the prefix `Failed to load field:`. Both are walked through in
{doc}`05_troubleshooting_fields_legends`.
:::

:::{note}
In practice the log-scale case is the one you will meet. A **prepared export
cannot contain `NaN` or an infinity** in a continuous obs field or a gene: both
the Python and the R writer refuse to publish one, and both refuse with counts
rather than a bare complaint — *“Continuous obs field 'qc_score' cannot be
published: of 18,142,044 cells, 12 NaN, 3 +Inf. First affected cells: …”*.

So a `NaN` grey only appears when you open a raw `.h5ad` or Zarr archive
directly in the browser, or serve one with `serve_anndata` / `show_anndata`,
where the values are read as they are. In server mode a gene or continuous
column that is not entirely finite is refused at the request: the server answers
HTTP 422 and puts the same counted reason in the body, which the viewer quotes
back to you.
:::

### Filtering vs coloring (two separate concepts)

Continuous legends expose a **visible range** (min/max sliders).

- The *filter* decides which cells are visible at all.
- The *color domain* decides how the remaining values map onto the colormap.

If you filter but keep the color domain wide, visible points can look “washed out”.
That’s why Cellucid also provides **Rescale colorbar to slider range**.

### Log scale (coloring only)

The **Log color scale** toggle affects the *color mapping*, not the filter math:

- With log scale **off**: values map linearly across the color domain.
- With log scale **on**: values map by `log10(value)` across the log-domain.

Important consequences:
- Values `<= 0` cannot be placed on a log axis, so they render as the neutral “None” gray.
- Your **Min/Max filter sliders still operate in the original units** (not log units).

The legend prints this rule beside the toggle itself, as
`Zero/negative/NaN values use the None color`.

:::{important}
**The toggle is refused rather than approximated when it cannot work.**

- A field with **no positive value at all** rejects the change with
  `Log color scale requires at least one positive field value.` The toggle
  stays off. No substitute `0–1` range is invented for it.
- A field that *does* have positive values but whose **currently selected
  colour range excludes them** rejects the change with
  `Log color scale requires the selected color range to include a positive
  value.` Widen the range, or turn off **Rescale colorbar to slider range**,
  and try again.

So “I turned log scale on and everything went grey” can never be caused by a
field with no positives — that case cannot turn the toggle on in the first
place.
:::

---

## Rescale colorbar to slider range (contrast control)

The toggle **Rescale colorbar to slider range** controls the color domain:

- **Off**: colors are mapped using the field’s full data range  
  (useful for consistent comparisons across different filters / views).
- **On**: colors are mapped using your current slider Min/Max range  
  (useful for bringing out contrast inside a narrow subset).

Practical example:
1) Filter a QC score to a tight range (e.g., “good cells only”).
2) Turn **Rescale colorbar to slider range** on so the remaining variation is visible.

---

## Live filtering vs FILTER (performance + predictability)

Continuous legends support two ways to apply the numeric filter:

- **Live filtering: On** (default)  
  Moving sliders immediately updates visibility — and `FILTER` is greyed out
  and unclickable while it is on, with the tooltip
  `Turn off Live filtering to use Filter`. That is the intended state, not a
  fault.
- **Live filtering: Off**  
  Move sliders freely, then click **FILTER** to apply once.

When you might turn Live filtering off:
- very large datasets where recomputing visibility on every slider tick feels slow,
- remote desktop / low-power machines where UI responsiveness matters.

The **RESET** button is clickable in both states and restores the full numeric
range for that field, snapping both sliders back to `0` and `100`.

:::{warning}
**`Live filtering` is not sticky.** It returns to `On` every time the legend is
rebuilt — which includes changing the active field, changing the colormap,
switching between kept views, and loading a session — and it is not stored in a
`.cellucid-session`. On a large dataset, check the toggle again after any of
those before you start dragging.
:::

The screenshot of the disabled `FILTER` button, and the rest of the filter
inventory, are in
{doc}`../e_filtering/04_common_filter_types_document_every_filter_the_ui_exposes`.

---

## How range sliders work (what “0–100” means)

The Min/Max sliders are normalized to **0–100%** of the field’s numeric range:

- `0` means the field’s minimum value
- `100` means the field’s maximum value

The UI shows the corresponding numeric value next to each slider so you don’t have to reason in percent.

The two sliders cannot cross. Dragging `Min` past `Max` pushes `Max` along with
it (and vice versa), so pushing one all the way collapses the range to a single
value rather than being refused.

:::{note}
**The sliders move in whole percent steps.** The finest bound you can reach
from the UI is therefore one hundredth of the field's full data range, and
there is no box to type an exact number into. Read the printed numeric readout
beside each slider to see what you actually selected, and do finer gating
upstream — in Python or R, before export. Export-time quantization compounds
with this; see {doc}`../e_filtering/06_edge_cases_filtering`.
:::

---

## Gene expression specifics (var fields)

Gene expression fields behave like continuous fields with a few practical differences:

- Genes are selected via search and **loaded on demand** (first selection can be slower).
- Many genes are sparse: a large fraction of cells can be exactly zero.
  - With **log scale on**, zeros become “None” gray.
  - With log scale off, zeros map to the low end of the colormap.

Current behavior to know:
- **Gene filters apply only while that gene is the active field.**  
  If you switch to an obs field, the gene range filter is not applied to visibility.

If you want “filter by gene, color by clusters” workflows, treat this as a current limitation and use:
- categorical/continuous obs filters, or
- analysis outputs exported as obs fields (then filter as obs).

---

## Outlier filter (latent space)

Some fields provide a per-cell **latent outlier quantile** (a dataset-provided statistic).
When available, the UI shows:

- **Outlier filter (latent space)** slider in Display options

Semantics:
- `100%` means **no outlier filtering**.
- Lower values remove more outliers:
  - `95%` removes cells above the 95th percentile of the outlier score (top ~5% most outlier-like).

Important details:
- Outlier filtering only applies when the current active field supports it (that’s when the control is visible).
- Derived/user-defined categorical fields created by merges/deletions typically do **not** carry latent outlier stats, so the outlier control may disappear after you “edit” a field.

---

## Filter stacking (why points can disappear “unexpectedly”)

Cellucid supports a filter stack across **obs fields**:
- hidden categories (from categorical legends),
- numeric ranges (from continuous legends),
- and (when active and supported) outlier filtering.

Key workflow pattern:
1) Filter using one field (e.g., continuous QC metric).
2) Switch the active field to something else (e.g., clusters) to color by it.

The filters remain active (and are summarized in **Active filters**) until you reset them.

:::{important}
Clearing the active field (choosing **None**) makes points gray, but it does **not** clear existing obs filters.
:::

---

## Quick checks when “colors look wrong”

Before assuming something is broken:

1) Look at **Coloring & Filtering**:
   - Is the active field actually the one you intended (categorical vs continuous vs gene)?
2) Look at the **legend**:
   - Is log scale on?
   - Is rescale on (which changes contrast)?
3) Look at **Active filters**:
   - Are you hiding categories on another field?
   - Did you narrow a continuous range earlier?
4) Look for the neutral “None” gray:
   - It usually means missing values, invalid codes, or log-scale-incompatible values.
5) Look at the notification centre:
   - If the colours did not change **at all** after you picked a gene or field,
     the payload was refused rather than drawn. Read the message: it names the
     field, counts the offending values and gives the first affected cells.

For a full symptom catalog, see {doc}`05_troubleshooting_fields_legends`.

---

## Next steps

- For the exact UI controls (what each toggle/button does), go to {doc}`04_legend_behavior`.
- For selection, rename, delete, and restore semantics, go to {doc}`02_field_selector_ux`.
