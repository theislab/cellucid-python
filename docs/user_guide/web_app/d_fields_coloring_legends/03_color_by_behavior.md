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

This is why the UI groups it as **Coloring & filtering**.

---

## Categorical color-by (discrete labels)

### What is colored, exactly?

For a categorical obs field, each cell has a category code:
- If the code is valid, the cell is colored with that category’s RGB color.
- If the code is missing/invalid, the cell uses a neutral gray “None” color.

Category colors come from the file when it stores them, else from a fixed
categorical palette, and can be overridden per category. A file colours a field
the way Scanpy keeps a plot's colours: `uns['<field>_colors']`, one hex colour
(`#rrggbb`) per category, in category order, as `sc.pl.umap(adata,
color="leiden")` leaves it. Cellucid reads it as Scanpy does: a list shorter
than the categories is not used, and a longer one gives its first colours. A
list holding anything but hex colours is not used, nor is it for a text column
whose categories are found when it is loaded. Once categories are merged, the
merged field is coloured by Cellucid's palette.

A text or true/false column (no stored categories) lists its labels in
code-point order and gives each the colour its label names (a hash of the
label picks it; of two labels naming one colour, the one first in that order
keeps it and the other takes the next free one), so the legend does not depend
on which cells are open: opening more cells adds the labels they hold, and
every label shown keeps its place in that order and its colour.

A categorical with more categories than there are open cells, and too many
for one request of the connection to read (barcodes or clone ids kept as a
categorical: one category per cell, a million or more on a web host), lists
only the categories its open cells hold, in the file's order, each in the colour it has
with every category listed (the file's colour at its place, else the palette's
colour at its place), and the legend says how many of the file's categories
those are ("12,034 of its 1,203,456 categories are among them"). Only those
categories' names and colours are read. Opening more cells adds the ones they
hold, in the file's order, every one shown keeping its colour and whether it
is hidden; a session or starting view that names a category no open cell holds
keeps what it set for it until cells holding it open. Any other categorical
lists every category it stores, the unused ones too.

A missing code is not a quiet variant of a category — it is outside the
category system entirely, and that has two visible consequences:

- Those cells are drawn in the neutral grey at **full opacity**, they get no
  legend row, and they are counted in no per-category count. They still count
  toward `Showing X of Y cells`.
- **No categorical filter can hide them.** Unchecking every row, or clicking
  `Hide all`, leaves them on screen. This is the usual explanation for “I hid
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
:alt: The categorical legend for clusters: the Show all and Hide all buttons above eight rows, each holding a checkbox, a colour swatch, a label, a cell count, a rename pencil and a delete bin.
:width: 488px

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
to filter them out. The grey is drawn beneath the measured cells (a missing
category's grey too), so a missing value never covers a measured one.

The grey **follows the background**. It is a neutral grey one fixed step of
lightness (CIE L* 15.4) from the background, toward the middle: on white it is
the light grey Cellucid always drew (RGB 211), on the light grid a shade darker
(202), and on dark backgrounds a mid-dark grey (55 on the dark grid, 39 on
black). A light grey on a dark background outshone the expressing cells; at
this step most of every colour map stands out from the background more than
the grey does, on every theme. The same grey is used for a gene's zeros, for
missing values of every field kind, in kept views and multiview, in smoke
mode, in the cell info card's swatch, in the legend's dot beside **Gray**, and
in exported figures (made from the figure's own paper, in every export
strategy).

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
In practice the log-scale case is the one you will meet. A file written by
`prepare()` or `cellucid_prepare()` **cannot contain an infinity** in a
continuous obs field, or any non-finite value in gene expression: both
writers refuse one, with counts rather than a bare complaint. A `NaN` in a
continuous obs column is kept: it is stored as a missing value and drawn in
the neutral grey.

Every file is read as it is, however it is opened, so a `NaN` grey appears
wherever the file holds one. An infinity is refused when the gene or column is
read, with the same counted reason.
:::

### Filtering vs coloring (two separate concepts)

Continuous legends expose a **visible range** (min/max sliders).

- The *filter* decides which cells are visible at all.
- The *color domain* decides how the remaining values map onto the colormap.

If you filter but keep the color domain wide, visible points can look “washed out”.
That’s why Cellucid also provides **Rescale to range**.

### Log scale (coloring only)

The **Log color scale** toggle affects the *color mapping*, not the filter math:

- With log scale **off**: values map linearly across the color domain.
- With log scale **on**: values map by `log10(value)` across the log-domain.

Important consequences:
- Values `<= 0` cannot be placed on a log axis, so they render as the neutral “None” gray.
- Your **Min/Max filter sliders still operate in the original units** (not log units).

The toggle's (i) states this rule: zero, negative and missing values have no
logarithm and take the None color.

:::{important}
**The toggle is refused rather than approximated when it cannot work.**

- A field with **no positive value at all** rejects the change with
  `Log color scale requires at least one positive field value.` The toggle
  stays off. No substitute `0–1` range is invented for it.
- A field that *does* have positive values but whose **currently selected
  range holds none** (with **Rescale to range** on) is not refused: every cell
  in that range is zero or negative and draws in the None grey, and the range
  filter still applies. Widen the range, or turn off **Rescale to range**.

So “I turned log scale on and everything went grey” is never a field with no
positives — that case cannot turn the toggle on — but it can be a slider range
that holds no positive value.
:::

---

## Rescale to range (contrast control)

The toggle **Rescale to range** controls the color domain:

- **Off**: colors are mapped using the field’s full data range  
  (useful for consistent comparisons across different filters / views).
- **On**: colors are mapped using your current slider Min/Max range  
  (useful for bringing out contrast inside a narrow subset).

Practical example:
1) Filter a QC score to a tight range (e.g., “good cells only”).
2) Turn **Rescale to range** on so the remaining variation is visible.

---

## Live filtering vs Filter (performance + predictability)

Continuous legends support two ways to apply the numeric filter:

- **Live filtering: On** (default)  
  Moving sliders immediately updates visibility, so there is no **Filter**
  button: only **Reset** sits under the sliders.
- **Live filtering: Off**  
  **Filter** appears beside **Reset**. Move sliders freely, then click
  **Filter** to apply once.

When you might turn Live filtering off:
- very large datasets where recomputing visibility on every slider tick feels slow,
- remote desktop / low-power machines where UI responsiveness matters.

The **Reset** button is clickable in both states and restores the full numeric
range for that field, snapping both sliders back to `0` and `100`.

:::{warning}
**`Live filtering` is not sticky.** It returns to `On` every time the legend is
rebuilt — which includes changing the active field, changing the colormap,
switching between kept views, and loading a session — and it is not stored in a
`.cellucid-session`. On a large dataset, check the toggle again after any of
those before you start dragging.
:::

The screenshot of the disabled `Filter` button, and the rest of the filter
inventory, are in
{doc}`../e_filtering/04_common_filter_types_document_every_filter_the_ui_exposes`.

---

## How the range bounds work

Each bound — **Minimum** and **Maximum** — has a slider and a number box:

- **Type** a bound in the box (a decimal point or comma; `1e-3` too) and press
  Enter or leave the box: it applies exactly as typed. A value past the
  field's minimum or maximum is that end; text that is not a number is not
  applied.
- **Drag** the slider: it moves in thousandths of the field's range (spaced in
  log units when **Log color scale** is on, the zeros at its first step), and
  the bound it sets is the shortest decimal that tells one step from the next,
  which the box then shows. The track's ends are the field's minimum and
  maximum exactly.

The box always holds the exact bound applied, and only the bound you move or
type changes: the other keeps its value, also one restored from a session.
The two cannot cross. Setting `Minimum` past `Maximum` pushes `Maximum` along
with it (and vice versa), so pushing one all the way collapses the range to a
single value rather than being refused.

---

## Gene expression specifics (var fields)

Gene expression fields behave like continuous fields with a few practical differences:

- Genes are selected via search and **loaded on demand** (first selection can be slower).
- Many genes are sparse: a large fraction of cells can be exactly zero.
  The cells where the gene is 0 (the cells that do not express it, for a gene
  with no negative value) are drawn in the neutral gray of a missing value,
  not the colormap's lowest colour, so the colours show where the gene is
  expressed; the legend says so under the colour bar. The gray is the gene's
  zero, never the colour range's minimum: with **Rescale to range** following
  a filter of "≥ 1" on counts, the cells at 1 keep the lowest colour, and any
  cell below a moved-up minimum takes the lowest colour, as scanpy's `vmin`
  clamps it.
  The gray cells are drawn **beneath** the expressing ones: in 2D every
  expressing cell lies on top of the gray around it, and in 3D it shows
  through gray on the same surface (a gray cell clearly in front still hides
  it), at every level of detail, in kept views and in exported figures.
  With **log scale on**, zeros are gray as well (they have no logarithm).
  A gene with negative values (scaled expression) has no gray floor, and
  neither has an obs field: a zero there is a value.
- A gene read for a **sample** of the open cells (a large atlas,
  {ref}`genes-on-a-sample-of-the-cells`) draws three states, which the legend
  names in one line, **Faded: not read · gray: 0**. The cells it was not read
  for are faded almost to the background (a shade of it, never the gray
  of 0), so the map keeps its shape. The cells read are drawn on top and
  larger, sized from how dense the sample is on screen, so the pattern shows
  at once: 2,000 cells on a 2.3M-cell map draw as clear dots, while a sample
  already dense on screen draws at the normal size. Their zeros are gray
  beneath the expressing cells. This holds in 2D and 3D, at every level of
  detail (every read cell stays drawn), in kept views and multiview, when you
  pick a cell (a read cell is found at its drawn size), and in exported
  figures, whose note says the same. The fade and the gray follow the
  background: light or dark, and a figure's own paper; the legend shows a
  dot of each colour beside its words.

A gene's range filter keeps applying after you colour by another field, so
"filter by gene, colour by clusters" works; **Filters** lists it.

---

## Outlier filter

For a categorical field, each cell has a **latent outlier quantile**: its
distance to its category's centre in a latent space (an `obsm` array such as
`X_scVI` or `X_pca`, never the 2D/3D layout), as the share of its category at
that distance or closer. The **Legend** box shows:

- an **Outlier filter** slider, below **Show centroids**, with a line under it
  when there is something to say.

Semantics:
- `100%` means **no outlier filtering**.
- Lower values remove more outliers:
  - `95%` hides the farthest ~5% of every category.

Important details:
- An outlier filter keeps applying after you colour by another field, and
  **Filters** lists it (`cell_type: outlier ≤ 90%`).
- It needs a latent space, which Cellucid never picks for you: the file's
  (`prepare(latent_key=...)`), the launch's (`--latent-key`), or yours under
  **Data mapping › Latent space**. Without one there is no slider: the filter
  is its label, its (i) and one **Choose latent space** button
  ([Turning it on](../e_filtering/02_outlier_filtering_per_active_field.md#turning-it-on-choose-a-latent-space)).
- The quantiles are got on the slider's first move below `100%` for a field
  (picking the field reads nothing): read from the file when it stores them
  for that latent space, at once; otherwise computed in the browser, reading
  the latent space once, with a progress bar and **Cancel**, after which every
  other field filters at once. A session saved with a threshold below `100%`,
  or a latent-space change while a field's filter is in use, gets them at once
  ([Stored or computed](../e_filtering/02_outlier_filtering_per_active_field.md#stored-or-computed)).
- A categorical field created by merging or deleting categories has no outlier
  quantiles, so the outlier control is hidden after you “edit” a field.

---

## Filter stacking (why points can disappear “unexpectedly”)

Cellucid supports a filter stack across **obs fields**:
- hidden categories (from categorical legends),
- numeric ranges (from continuous legends),
- and (when active and supported) outlier filtering.

Key workflow pattern:
1) Filter using one field (e.g., continuous QC metric).
2) Switch the active field to something else (e.g., clusters) to color by it.

The filters remain active (and are summarized in **Filters**) until you reset them.

:::{important}
Clearing the active field (choosing **None**) makes points gray, but it does **not** clear existing obs filters.
:::

---

## Quick checks when “colors look wrong”

Before assuming something is broken:

1) Look at **Coloring & filtering**:
   - Is the active field actually the one you intended (categorical vs continuous vs gene)?
2) Look at the **legend**:
   - Is log scale on?
   - Is rescale on (which changes contrast)?
3) Look at **Filters**:
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
