# Legend behavior

**Audience:** everyone (basic legend use) + power users (category editing workflows)  
**Time:** 25–40 minutes  
**What you’ll learn:**
- Where the legend lives and when it appears
- Everything you can do in a **categorical legend** (hide/show, keep only some cells, recolor, rename, merge, delete-to-Unassigned)
- Everything you can do in a **continuous legend** (colormaps, log scale, rescale, filters)
- How counts and disabled rows behave when other filters are active
- What extra controls exist (centroids, latent outlier filter) and how they interact with views

---

## Where the legend lives (Legend)

The legend appears under:

- **Coloring & filtering → Legend**

It is hidden until you select an active field (categorical obs, continuous obs, or a gene).

:::{note}
The sidebar legend reflects the **currently active view** in multiview.
If you have multiple snapshots, click a panel to activate it before editing legend settings.
:::

---

## Categorical legend (checkboxes + swatches + category tools)

You see a categorical legend when the active field is **Categorical obs**.

### Anatomy of one category row

Each category row contains:
- a **checkbox** (show/hide the category)
- a **color swatch** (click to choose a color)
- a **category label**
- a **cell count** (how many cells are visible/available); click it to keep
  only some of the category's cells ([below](#keep-only-some-cells-of-a-category-the-count))
- optional tools: **rename** (pencil) and **delete** (trash)

```{figure} ../../../_static/screenshots/web_app/legend-categorical-clusters.png
:alt: The categorical legend for clusters: the Show all and Hide all buttons above eight rows, each holding a checkbox, a color swatch, a label, a cell count, a rename pencil and a delete bin.
:width: 488px

One row, read left to right: **checkbox · swatch · label · count · rename ·
delete**. The pencil renames the category; the bin merges it into
`Unassigned` and does **not** delete any cells. The rows run in the field's
own category order (see {doc}`03_color_by_behavior`).
```

### Show all / Hide all

Buttons at the top:
- **Show all**: checks all currently available categories
- **Hide all**: unchecks all currently available categories

“Available” matters when other filters are active:
- Categories with **0 available cells** are disabled (you can’t toggle them via **Show all** / **Hide all**).

### Category counts: “visible / available”

Counts are shown as:
- `1,234 cells` when visible equals available
- `visible / available cells` when they differ

Interpretation:
- **available** ≈ cells remaining after *other* filters (but before hiding this category)
- **visible** ≈ cells currently shown after applying category visibility

A category that keeps only some of its cells shows its count in the link
color (`200 / 1,200 cells`); see
[Keep only some cells of a category](#keep-only-some-cells-of-a-category-the-count).

If a category has `0 available` cells, its checkbox and color picker are
disabled, the row is grayed out, and hovering it shows
`No cells available in this category after other filters`.

On a subsample the counts are of the open cells, and a line above them says
so: *Counted on the 240 of 1,200 cells open.*
({doc}`../b_data_loading/21_subsampled_loading`).

```{figure} ../../../_static/screenshots/filtering/legend-category-unavailable.png
:alt: The clusters legend under a narrow S_score range filter: Alpha, Ductal, Ngn3 high EP and Ngn3 low EP stay ticked with small counts, while Beta, Delta, Epsilon and Pre-endocrine are grayed out with disabled checkboxes and read 0 cells.
:width: 488px

Two of the three count states at once: `0 cells` on four disabled rows, and
plain counts on the four rows that still have cells. The third,
`visible / available`, is printed only when the two numbers differ *and*
available is above zero.
```

### Change a category color (swatch)

Click a swatch to open a color picker and choose a color.
This updates:
- point colors in the viewer,
- centroid colors (if shown),
- and the saved session state (if you save a session bundle).

:::{tip}
Pick high-contrast colors for categories you plan to compare side-by-side in multiview.
:::

### Hide/show categories (checkboxes)

Unchecking a category:
- makes those cells invisible (fully transparent),
- and adds a categorical filter entry to **Filters**.

This is why hiding categories is a powerful workflow tool (not just cosmetic).

### Keep only some cells of a category (the count)

A large category can crowd out the others. Its count is a button: click it
to show only some of its cells.

1. Click the category's count (for example `1,200 cells`). A small box takes
   its place, reading `all`.
2. Type how many cells to keep, for example `200`, and press {kbd}`Enter`
   (clicking elsewhere applies it too; {kbd}`Escape` cancels).
3. The category now shows 200 of its cells, chosen at random. Its count reads
   `200 / 1,200 cells` in the link color, and hovering it says *Keeps a
   random 200 cells of Beta (seed 0); click to change, or type all*.

To show every cell of it again, click the count and type `all` (or any number
at least as large as its count). **Show all**, **Hide all**, and removing the
field's row under **Filters**, also show every cell again.

What to expect:

- **A seeded random sample.** The cells are drawn at random with a seed (0
  from the legend), so the same count always shows the same cells, whichever
  cells are open and in whatever order.
- **A larger count adds cells.** Raising 200 to 500 keeps the 200 shown and
  adds 300 more; lowering it keeps some of those shown.
- **It is a filter.** It hides cells as unticking a category does, and is
  listed under **Filters** with its field (`cell_type: keeping 200 of Beta`);
  unticking that row turns it off ({doc}`../e_filtering/03_filter_stack_ui_active_filters`).
- **It is undoable and saved.** Each count you type is one step of the
  history, so **Undo** takes it back
  ({doc}`../m_state_history/01_the_history`); a session file and the saved
  states this browser keeps hold it.
- **On a subsample**, the cells are drawn from the category's open cells.

From a notebook, `viewer.filter_category_cell_count(field={"source": "obs",
"key": "cell_type"}, category="Beta", count=200, seed=1)` does the same, with
any seed.

---

### Cells without a category (Missing)
When some cells have no category (an unlabeled cell, a classifier that gave
none), a **Missing** row under the categories counts them (`visible /
available`) and shows or hides them with its checkbox. They are drawn in the
neutral “None” gray. Hiding them is a filter like hiding a category: listed
under **Filters**, applied whatever is colored, saved in sessions. The row is
absent when every cell has a category.

### Categories with no name: (empty) and (blank)
A category whose name is empty is shown as *(empty)*, and one made of spaces
alone as *(blank)*, in a quieter style, so its row can be told apart and
found; a name with spaces at an end is shown in quotes (`" x "`), so it is
told from `x`. The same names appear wherever a category is named: the hover
card, the filters, menus, the Group categories dialog, figure legends,
analysis groups and what an agent reads (an agent can name the category so).
Only the display changes: the value in the data, in sessions and in exported
values stays what the file holds.

## Editing category labels (rename, merge, delete)

These operations are designed for “clean up labels for interpretation/figures”.
Cellucid tries to keep the workflow safe by preserving originals.

:::{important}
Editing behavior depends on whether the current field is a **source field** or a **user-defined derived field**:

- For source fields: destructive edits usually create a **new derived field** and move the original into **Deleted fields** (restorable).
- For already derived fields: repeated edits typically happen **in place** (faster, less clutter).

An edit in place changes the field in every view at once: a kept view that
shows the field is drawn again with its new categories, and keeps its own
colors and hidden categories for the categories it already had.
:::

:::{tip}
To put several categories together at once (a lineage from cell types),
use **Group** on the Categorical obs row instead of merging them one pair at
a time: it makes a new field of the groups in one step and leaves this one as
it is ({ref}`group-categories-into-a-new-field`).
:::

### A long category name

Each category keeps one line in the sidebar, so a long name is cut with "…".
Click (or tap) the name to see it whole; click it again, or anywhere else, to
close it. Hovering shows it too, above the row's hint.

### Rename a category

Two ways:
- click the pencil icon, or
- double-click the category label

Type the name and press {kbd}`Enter`, or click elsewhere; {kbd}`Escape` leaves
it as it was. A name the rules below refuse is said on the box, which stays
open until you change it.

Rules (the same for every label you give: a rename, a merged category's
label, a group of **Group**, a page's label in **Create categorical**):
- label cannot be empty
- one line, with no leading or trailing whitespace and no control or
  zero-width character (`U+200B`, `U+2060`, `U+FEFF`) —
  `"<label>" cannot be a category label: …`
- label must be unique in the field

A label a dataset brings is shown as it is, whatever it holds.

Each rename is one step of the history: **Undo** (`Ctrl+Z` / `⌘Z`) takes it
back.

### Merge categories (drag-and-drop)

To merge `A` into `B`:
1) drag the label for `A`
2) drop it onto the row for `B`
3) confirm in the dialog

Cellucid creates a merged bucket label like:
- `merged A + B` (made unique if needed)

The merged bucket color comes from the **dragged** category.

An agent (or a notebook, `merge_categories(..., label=...)`) merging `A` into
`B` can name the merged category: given `label` `B`, the merged category
keeps the name `B` (the color is still the dragged one's); without a label it
is named as the drag names it.

### Delete a category label (merge into “Unassigned”)

Click the trash icon to “delete” a label.

This does not delete cells. It changes the category assignment by merging that category into:
- `Unassigned` (a special catch-all bucket, made when the field has none; the
  same spelling **Create categorical** gives the cells on no page). An agent
  or a notebook can name another (`unassigned_label`). A new `Unassigned`
  takes the palette color of its place when it can be told apart from the
  categories' colors, else the next palette color that can: it never shares a
  kept category's color.

Colors the app picks (a new `Unassigned`, a category a relabel adds, a
group, the categories **Create categorical** adds beyond its pages, a new
highlight page) are told apart by one rule: two colors closer than 0.08 in
OKLab, a color space in which equal distances look equally different, are
taken for one, and the next palette color that is far enough from all the
others is used instead. Scanpy's default ten colors are all farther apart than
that. Colors you or a file chose are never changed by it; a file's categories
keep the color of their place, so a share of the cells and all of them show
one color.

This is useful for:
- removing junk labels like `doublet`/`unknown`/`low_quality`,
- iteratively collapsing labels while keeping a reversible history.

A merge and a delete are each one step of the history: **Undo** (`Ctrl+Z` /
`⌘Z`) takes either back, whether it edited a field of your own in place or
made a new field and moved the file's field to **Deleted fields**. The
confirmation says which of the two it does.

---

## Continuous legend (colorbar + palette + range sliders)

You see a continuous legend when the active field is:
- **Continuous obs**, or
- **Gene expression**

### Color scale (palette / colormap)

At the top, you’ll see:
- `Color scale (<ColormapName>)`
- a gradient bar you can click
- for a gene shown in transformed values, a line naming their units, as
  `Values: log1p(norm.) (Data mapping › Expression values)` (`norm.`,
  `log1p` or `log1p(norm.)`; none for values as stored or an obs field)

Click the gradient bar to open the **colormap menu**, then click a palette to apply it.

The numbers under the bar are the scale's ends, written as every measured
value in the app is: three significant figures, in scientific notation below
0.01 and from 10,000 (`0.0123`, `1230`, `1.23e4`), the notation chosen after
rounding (`9999.6` is `1e4`, never `10.0e3`). Exported figures write their
color bars and axis ticks by the same rule (ticks with as many figures as
tell them apart).

```{figure} ../../../_static/screenshots/web_app/legend-colormap-menu-open.png
:alt: The colormap list open under the S_score color bar: seven named gradients, Viridis, Plasma, Inferno, Magma, Cividis, Turbo and Twilight, with Viridis marked as the active one.
:width: 488px

The colormap menu, open. The active palette is marked, and the heading above
the bar always names it — `Color scale (Viridis)` — so you can read the current
palette without opening the menu.
```

Practical guidance:
- Prefer **Viridis** / **Cividis** for scientific readability (especially colorblind-friendly workflows).
- Use diverging palettes (e.g., **RdBu**) only when your values have a meaningful “center” (like 0).
  A diverging palette (Coolwarm, Blue-Red, Spectral, RdYlBu, RdBu, PiYG,
  PRGn, BrBG, Seismic) on a linear scale whose range spans 0 is **centered at
  0**: the scale runs from −m to m, m the larger of |minimum| and |maximum|,
  so 0 is drawn in its neutral middle (z-scores, scores, fold changes). The
  legend says *Centered at 0* under the bar, and exported figures note it
  under theirs. Over data of one sign, or with a sequential palette, the scale
  spans the data's own range.

(gray-for-zero)=
### Gray for zero

**Gray for zero** sits above **Log color scale**, with an (i) that explains
it. It is shown for a gene with no negative value (a scaled layer has none of
this, nor has an obs field: a zero there is a value).

- **Off** (the default): a gene's 0 is a value like any other and is drawn in
  the color scale's lowest color.
- **On**: the cells where the gene is 0 are drawn in the theme's neutral gray
  (made from the background: a step darker than a light one, a step lighter
  than a dark one), beneath the cells that express it, so the colors show
  only where the gene is expressed.
- On a **log scale** a 0 has no logarithm and is gray either way; the switch is
  then grayed out and says so on hover.

It is per field, like the log scale and Rescale to range: each gene keeps its
own setting, kept views and multiview keep the setting of the view it was set
in, and a `.cellucid-session` saves it (only when it is on). Everything that
draws the coloring follows it: 2D and 3D, smoke, kept and compared views,
the cell info card's swatch, and exported figures (PNG, SVG and hybrid). A
figure whose zeros are gray says **Gray: 0** under its color bar (with a gene
sample, **Faded: not read · gray: 0**).

### Log scale toggle

**Log color scale** changes the color mapping to log10 for positive values.

Notes:
- Zeros, negatives, and missing values render as neutral “None” gray.
- Log scale affects coloring only; filtering remains in original units. With
  **Rescale to range** on, a slider range that holds no positive
  value simply draws its cells gray; the range filter still applies.

### Rescale to range

**Rescale to range** controls whether the colormap spans:
- the full data range (Off), or
- your current slider range (On).

On is best for bringing out contrast in a filtered subset; Off is best for consistent cross-view comparisons.

---

## Filtering controls inside the continuous legend

Under **Visible range** you’ll see:
- **Live filtering** toggle
- **Minimum** and **Maximum**, each on one row: its label, a slider and a
  number box, the two rows aligned in columns. Type an exact bound, or drag
  (thousandths of the range; log-spaced on a log color scale). Only the bound
  you change moves. The box shows the bound in as many figures as fit it (at
  most seven characters: `5.65685` of 5.656854249492381, `1.23e11`), never cut
  off; hovering it shows the exact bound applied, and a bound you type is
  applied exactly as typed.
- **Reset**, and **Filter** beside it while Live filtering is off

Behavior summary:
- With Live filtering on, visibility updates as you drag sliders, and there is
  no **Filter** button: there would be nothing for it to do.
- With Live filtering off, **Filter** appears: you drag first, then click
  Filter to apply.
- Reset works in both states and restores the full numeric range for that field.

Filtering notes:
- `Live filtering` starts `On`. Turned off, it stays off for the page — for
  every continuous field and gene, every view and palette — until you turn it
  on again or reload the page. It is a way of working rather than something
  shown, so a `.cellucid-session` does not save it.
- Missing (`NaN`) values are excluded when a numeric filter is active.
- These filters appear in the **Filters** box and (for obs fields) stack with other obs filters.

The complete filter inventory is in
{doc}`../e_filtering/04_common_filter_types_document_every_filter_the_ui_exposes`.

---

## Extra controls in the Legend section

### Outlier filter

When the active field is a categorical field of the file, you’ll see, below
**Show centroids**:
- the **Outlier filter**: until a latent space is chosen, only its label,
  its (i) and one **Choose latent space** button; then the slider, with at
  most one short line under it saying what its first use reads (for example
  "First use reads X_scVI (12 MB).") or, while that read runs, a progress bar
  and **Cancel**

See {doc}`03_color_by_behavior` for semantics and edge cases.

### Centroids (categorical fields only)

When the active field is categorical, you’ll also see:
- **Show centroids**
- **Show labels**

These controls toggle additional overlays that summarize categories in the embedding.

```{figure} ../../../_static/screenshots/web_app/window-centroids-and-labels.png
:alt: The Pancreas embedding colored by clusters with Show centroids and Show labels both ticked, so each cluster carries a marker at its center with its category name printed on it.
:width: 1440px

**Show centroids** puts a marker at each category's center; **Show labels**
prints its name there. Together they turn a color legend you have to
cross-reference into a picture you can read directly — which is usually what
you want in a figure.
```

Where the centers come from:
- a dataset's own categorical fields get each category's center computed in
  the browser from its codes and the cells' positions, once per embedding the
  first time a view draws it, and cached on the field;
- a field you create — with `Create categorical`, by merging or deleting a
  category, by duplicating a field — gets each category's center as the plain
  mean of its cells' positions in the embedding the view draws. Centers are
  computed when the field is made, for every embedding already loaded; an
  embedding loaded later gets its centers when a view switches to it. Merging
  or deleting its categories recomputes them, and a saved session stores the
  field's codes, from which its centers come back exactly where they were.

Notes for multiview users:
- centroid visibility is tracked per view, so different snapshots can show/hide centroids independently.

Notes on counts:
- a category no cell belongs to has no center, so it gets no marker and no
  label: a field whose category list is far longer than the categories its
  cells use draws only the ones they use. A category that gains its first
  cell (a merge, a rebuilt consensus column) gets its centroid then;
- a centroid disappears when its category has no visible cells left, so
  centroids follow your filters rather than the raw data.

---

## Optional variant: legend behavior in Community Annotation voting mode

If community annotation voting is enabled for the current categorical field:
- category labels become clickable to open a voting modal,
- rename/merge/delete label tools are disabled (to keep votes stable),
- additional per-category summary lines can appear (consensus suggestions),
- a small `◉` button can appear to quickly create/remove a highlight for that category. It highlights the category's currently visible cells, so it is disabled while the category has none (for example while its checkbox is cleared).

If you don’t use Community Annotation, you can ignore this.

---

## Next steps

- For selection/rename/delete/restore behavior around fields, read {doc}`02_field_selector_ux`.
- For a symptom-driven debug guide, read {doc}`05_troubleshooting_fields_legends`.
