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

<!-- SCREENSHOT PLACEHOLDER
ID: categorical-legend-anatomy
Suggested filename: web_app/fields_legends/30_categorical-legend-anatomy.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 04_legend_behavior.md
Capture:
  - Select a categorical field with at least ~10 categories
  - Make sure some other filter is active so at least one row shows “visible / available”
  - Ensure rename (pencil) and delete (trash) icons are visible
Crop:
  - Include: the categorical legend section including Show All / Hide All
Annotations:
  - Call out: checkbox, swatch, label, count, rename icon, delete icon, Show All/Hide All
Alt text:
  - Categorical legend with category rows, checkboxes, color swatches, counts, and edit icons.
Caption:
  - Explain the meaning of each row control and that category visibility is a real filter.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for categorical legend anatomy.
:width: 100%

Categorical legends let you hide/show categories, adjust colors, and edit labels directly from the legend.
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

If a category has `0 available` cells, its checkbox and color picker are disabled and the row is greyed out.

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

<!-- SCREENSHOT PLACEHOLDER
ID: categorical-merge-confirm-dialog
Suggested filename: web_app/fields_legends/31_categorical-merge-confirm.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 04_legend_behavior.md
Capture:
  - Start a drag from one category label and drop onto another
  - Capture the confirm dialog (“Merge categories”)
Crop:
  - Include: the dialog + enough background legend rows to show what was merged
Annotations:
  - Call out: from-label, to-label, and the “this will create a new column vs edit in place” explanation text
Alt text:
  - Confirmation dialog for merging two categorical labels.
Caption:
  - Emphasize that merges are reversible via Deleted Fields (for source fields) and may edit in place for derived fields.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the merge categories confirmation dialog.
:width: 100%

Merging categories prompts for confirmation and explains whether the edit will create a new derived field or edit in place.
```

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

Practical guidance:
- Prefer **Viridis** / **Cividis** for scientific readability (especially colorblind-friendly workflows).
- Use diverging palettes (e.g., **RdBu**) only when your values have a meaningful “center” (like 0).

<!-- SCREENSHOT PLACEHOLDER
ID: continuous-legend-colormap-menu
Suggested filename: web_app/fields_legends/32_continuous-colormap-menu.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 04_legend_behavior.md
Capture:
  - Select a continuous field
  - Click the gradient bar so the colormap menu opens (show several options)
Crop:
  - Include: “Color scale (…)” header + gradient bar + the open colormap menu list
Annotations:
  - Call out: gradient bar click target, one palette option, and the active palette highlight
Alt text:
  - Continuous legend with the colormap menu open showing selectable palettes.
Caption:
  - Explain how to change palettes and why the palette choice matters for interpretation.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the continuous colormap menu.
:width: 100%

Click the gradient bar in a continuous legend to open the colormap menu and choose a palette.
```

### Log scale toggle

**Use log scale** changes the color mapping to log10 for positive values.

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

See `d_fields_coloring_legends/03_color_by_behavior` for semantics and edge cases.

### Centroids (categorical fields only)

When the active field is categorical, you’ll also see:
- **Show centroids**
- **Show labels**

These controls toggle additional overlays that summarize categories in the embedding.

Notes for multiview users:
- centroid visibility is tracked per view, so different snapshots can show/hide centroids independently.

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

- For selection/rename/delete/restore behavior around fields, read `d_fields_coloring_legends/02_field_selector_ux`.
- For a symptom-driven debug guide, read `d_fields_coloring_legends/05_troubleshooting_fields_legends`.
