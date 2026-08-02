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

## Where the legend lives (Display options)

The legend appears under:

- **Coloring & Filtering → Display options**

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
:alt: A legend with the hint Click a swatch to pick a color at the top, Show All and Hide All buttons below it, and eight rows; each row has a blue ticked checkbox, then a coloured square swatch, then a category name such as Alpha or Ngn3 high EP, then a right-aligned cell count such as 481 cells, then a small pencil icon and a small bin icon.
:width: 428px

One row, read left to right: **checkbox · swatch · label · count · rename ·
delete**. The pencil renames the category; the bin merges it into
`unassigned` and does **not** delete any cells.
```

### Show All / Hide All

Buttons at the top:
- **Show All**: checks all currently available categories
- **Hide All**: unchecks all currently available categories

“Available” matters when other filters are active:
- Categories with **0 available cells** are disabled (you can’t toggle them via Show/Hide All).

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

```{figure} ../../../_static/screenshots/filtering/legend-category-unavailable.png
:alt: The same legend under an additional filter, where Beta, Delta, Epsilon and Pre-endocrine are drawn in grey with unticked disabled checkboxes and read 0 cells, while Alpha, Ductal, Ngn3 high EP and Ngn3 low EP remain black and ticked with counts of 1, 122, 10 and 16 cells.
:width: 428px

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
- and adds a categorical filter entry to **Active filters**.

This is why hiding categories is a powerful workflow tool (not just cosmetic).

---

## Editing category labels (rename, merge, delete)

These operations are designed for “clean up labels for interpretation/figures”.
Cellucid tries to keep the workflow safe by preserving originals.

:::{important}
Editing behavior depends on whether the current field is a **source field** or a **user-defined derived field**:

- For source fields: destructive edits usually create a **new derived field** and move the original into **Deleted Fields** (restorable).
- For already derived fields: repeated edits typically happen **in place** (faster, less clutter).
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
- restore the original field from **Deleted Fields** (if a new derived field was created).

---

## Continuous legend (colorbar + palette + range sliders)

You see a continuous legend when the active field is:
- **Continuous obs**, or
- **Gene Expression**

### Color scale (palette / colormap)

At the top, you’ll see:
- `Color scale (<ColormapName>)`
- a gradient bar you can click

Click the gradient bar to open the **colormap menu**, then click a palette to apply it.

```{figure} ../../../_static/screenshots/web_app/legend-colormap-menu-open.png
:alt: A legend headed COLOR SCALE (VIRIDIS) with its gradient bar expanded into a grid of small named gradient swatches — Viridis highlighted as active, then Plasma, Inferno, Magma, Cividis, Turbo, Twilight, Coolwarm and many more — above the numeric range and the Log color scale and Rescale colorbar toggles.
:width: 428px

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
- Log scale affects coloring only; filtering remains in original units.

### Rescale colorbar to slider range

**Rescale colorbar to slider range** controls whether the colormap spans:
- the full data range (Off), or
- your current slider range (On).

On is best for bringing out contrast in a filtered subset; Off is best for consistent cross-view comparisons.

---

## Filtering controls inside the continuous legend

Under **Filtering → Visible range** you’ll see:
- **Live filtering** toggle
- **Min** slider and **Max** slider
- **FILTER** and **RESET** buttons

Behavior summary:
- With Live filtering on, visibility updates as you drag sliders.
- With Live filtering off, you drag first, then click FILTER to apply.
- RESET restores the full numeric range for that field.

Filtering notes:
- Missing (`NaN`) values are excluded when a numeric filter is active.
- These filters appear in the **Active filters** box and (for obs fields) stack with other obs filters.

---

## Extra controls in Display options

### Outlier filter (latent space)

If the active field provides latent outlier quantiles, you’ll see:
- **Outlier filter (latent space)** slider

See {doc}`03_color_by_behavior` for semantics and edge cases.

### Centroids (categorical fields only)

When the active field is categorical, you’ll also see:
- **Show centroids**
- **Show labels**

These controls toggle additional overlays that summarize categories in the embedding.

```{figure} ../../../_static/screenshots/web_app/window-centroids-and-labels.png
:alt: The whole application window with the Pancreas embedding coloured by clusters, both Show centroids and Show labels ticked in the sidebar, and small markers with the category names Alpha, Beta, Delta, Ductal, Epsilon, Ngn3 high EP, Ngn3 low EP and Pre-endocrine drawn on the canvas at the centre of each coloured group.
:width: 1440px

**Show centroids** puts a marker at each category's centre; **Show labels**
prints its name there. Together they turn a colour legend you have to
cross-reference into a picture you can read directly — which is usually what
you want in a figure.
```

Notes for multiview users:
- centroid visibility is tracked per view, so different snapshots can show/hide centroids independently.

Notes on counts:
- a centroid disappears when its category has no visible cells left, so
  centroids follow your filters rather than the raw data.

---

## Optional variant: legend behavior in Community Annotation voting mode

If community annotation voting is enabled for the current categorical field:
- category labels become clickable to open a voting modal,
- rename/merge/delete label tools are disabled (to keep votes stable),
- additional per-category summary lines can appear (consensus suggestions),
- a small `◉` button can appear to quickly create/remove a highlight for that category.

If you don’t use Community Annotation, you can ignore this.

---

## Next steps

- For selection/rename/delete/restore behavior around fields, read {doc}`02_field_selector_ux`.
- For a symptom-driven debug guide, read {doc}`05_troubleshooting_fields_legends`.
