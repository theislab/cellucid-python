# Legend behavior

**Audience:** everyone (basic legend use) + power users (category editing workflows)  
**Time:** 25–40 minutes  
**What you’ll learn:**
- Where the legend lives and when it appears
- Everything you can do in a **categorical legend** (hide/show, recolor, rename, merge, delete-to-unassigned)
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
- a **cell count** (how many cells are visible/available)
- optional tools: **rename** (pencil) and **delete** (trash)

```{figure} ../../../_static/screenshots/web_app/legend-categorical-clusters.png
:alt: The categorical legend for clusters: the Show all and Hide all buttons above eight rows, each holding a checkbox, a colour swatch, a label, a cell count, a rename pencil and a delete bin.
:width: 488px

One row, read left to right: **checkbox · swatch · label · count · rename ·
delete**. The pencil renames the category; the bin merges it into
`unassigned` and does **not** delete any cells.
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

If a category has `0 available` cells, its checkbox and color picker are
disabled, the row is greyed out, and hovering it shows
`No cells available in this category after other filters`.

On a subsample the counts are of the open cells, and a line above them says
so: *Counted on the 240 of 1,200 cells open.*
({doc}`../b_data_loading/21_subsampled_loading`).

```{figure} ../../../_static/screenshots/filtering/legend-category-unavailable.png
:alt: The clusters legend under a narrow S_score range filter: Alpha, Ductal, Ngn3 high EP and Ngn3 low EP stay ticked with small counts, while Beta, Delta, Epsilon and Pre-endocrine are greyed out with disabled checkboxes and read 0 cells.
:width: 488px

All three count states at once: `0 cells` on four disabled rows, and small
plain counts on the four rows that still have cells. A row only prints
`visible / available` when the two numbers differ *and* available is above
zero.
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

---

## Editing category labels (rename, merge, delete)

These operations are designed for “clean up labels for interpretation/figures”.
Cellucid tries to keep the workflow safe by preserving originals.

:::{important}
Editing behavior depends on whether the current field is a **source field** or a **user-defined derived field**:

- For source fields: destructive edits usually create a **new derived field** and move the original into **Deleted fields** (restorable).
- For already derived fields: repeated edits typically happen **in place** (faster, less clutter).

An edit in place changes the field in every view at once: a kept view that
shows the field is drawn again with its new categories, and keeps its own
colours and hidden categories for the categories it already had.
:::

### Rename a category

Two ways:
- click the pencil icon, or
- double-click the category label

Rules:
- label cannot be empty
- max label length is enforced (very long labels are rejected)

Undo:
- rename the category back to its original label (Cellucid tracks originals).

### Merge categories (drag-and-drop)

To merge `A` into `B`:
1) drag the label for `A`
2) drop it onto the row for `B`
3) confirm in the dialog

Cellucid creates a merged bucket label like:
- `merged A + B` (made unique if needed)

The merged bucket color comes from the **dragged** category.

### Delete a category label (merge into “unassigned”)

Click the trash icon to “delete” a label.

This does not delete cells. It changes the category assignment by merging that category into:
- `unassigned` (a special catch-all bucket)

This is useful for:
- removing junk labels like `doublet`/`unknown`/`low_quality`,
- iteratively collapsing labels while keeping a reversible history.

Undo:
- restore the original field from **Deleted fields** (if a new derived field was created).

---

## Continuous legend (colorbar + palette + range sliders)

You see a continuous legend when the active field is:
- **Continuous obs**, or
- **Gene expression**

### Color scale (palette / colormap)

At the top, you’ll see:
- `Color scale (<ColormapName>)`
- a gradient bar you can click

Click the gradient bar to open the **colormap menu**, then click a palette to apply it.

```{figure} ../../../_static/screenshots/web_app/legend-colormap-menu-open.png
:alt: The colormap list open under the S_score colour bar: seven named gradients, Viridis, Plasma, Inferno, Magma, Cividis, Turbo and Twilight, with Viridis marked as the active one.
:width: 488px

The colormap menu, open. The active palette is marked, and the heading above
the bar always names it — `Color scale (Viridis)` — so you can read the current
palette without opening the menu.
```

Practical guidance:
- Prefer **Viridis** / **Cividis** for scientific readability (especially colorblind-friendly workflows).
- Use diverging palettes (e.g., **RdBu**) only when your values have a meaningful “center” (like 0).

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
- **Minimum** slider and **Maximum** slider
- **Filter** and **Reset** buttons

Behavior summary:
- With Live filtering on, visibility updates as you drag sliders — and **Filter
  is greyed out and unclickable**, with the tooltip
  `Turn off Live filtering to use Filter`. There is nothing for it to do.
- With Live filtering off, you drag first, then click Filter to apply.
- Reset works in both states and restores the full numeric range for that field.

Filtering notes:
- `Live filtering` starts `On`, and it goes back to `On` every time the legend
  is rebuilt — a different active field, a different colormap, a different kept
  view, or a loaded session. It is not saved in a `.cellucid-session`.
- Missing (`NaN`) values are excluded when a numeric filter is active.
- These filters appear in the **Filters** box and (for obs fields) stack with other obs filters.

The complete filter inventory, with a screenshot of the disabled `Filter`
button, is in
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
:alt: The Pancreas embedding coloured by clusters with Show centroids and Show labels both ticked, so each cluster carries a marker at its centre with its category name printed on it.
:width: 1440px

**Show centroids** puts a marker at each category's centre; **Show labels**
prints its name there. Together they turn a colour legend you have to
cross-reference into a picture you can read directly — which is usually what
you want in a figure.
```

Where the centres come from:
- a dataset's own categorical fields get each category's centre computed in
  the browser from its codes and the cells' positions, once per embedding the
  first time a view draws it, and cached on the field;
- a field you create — with `Create categorical`, by merging or deleting a
  category, by duplicating a field — gets each category's centre as the plain
  mean of its cells' positions in the embedding the view draws. Centres are
  computed when the field is made, for every embedding already loaded; an
  embedding loaded later gets its centres when a view switches to it. Merging
  or deleting its categories recomputes them, and a saved session stores the
  field's codes, from which its centres come back exactly where they were.

Notes for multiview users:
- centroid visibility is tracked per view, so different snapshots can show/hide centroids independently.

Notes on counts:
- a category no cell belongs to has no centre, so it gets no marker and no
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
