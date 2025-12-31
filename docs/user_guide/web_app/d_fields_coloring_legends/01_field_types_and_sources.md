# Field types and sources

**Audience:** everyone (non-technical users + computational users)  
**Time:** 10–15 minutes  
**What you’ll learn:**
- What a “field” is in Cellucid
- The difference between **obs** (cell metadata) and **var** (genes / features)
- The difference between **categorical** and **continuous** fields
- What extra “field features” exist (centroids, latent outlier quantiles) and how they show up in the UI

---

## Mental model: a “field” is a per-cell array

In Cellucid, a **field** is simply “one value per cell”:

- A **categorical field** assigns each cell a label (e.g., `B cell`, `T cell`, `cluster_7`).
- A **continuous field** assigns each cell a number (e.g., `n_counts`, `pct_mito`, `pseudotime`).

The viewer then uses the active field to drive:
- how points are **colored** (color-by),
- what controls appear in the **legend** (checkboxes vs sliders),
- and often which points are **visible** (filters you create through the legend).

:::{tip}
If you’re not sure what something is: try selecting it and looking at the legend.

- If you see **checkboxes + colored swatches**, it’s categorical.
- If you see a **colorbar + min/max sliders**, it’s continuous.
:::

---

## “Obs” vs “Var” (what the words mean)

These names come from AnnData conventions, but you don’t need AnnData to use Cellucid.

### Obs (observations) = per-cell metadata

Obs fields are cell-level columns such as:
- cell type / cluster labels
- sample / batch / donor
- QC metrics (`n_counts`, `n_genes`, `pct_mito`)
- scores (`S_score`, `G2M_score`, marker module scores)
- pseudotime / diffusion pseudotime / latent dimensions (if exported as obs)

In the UI, obs fields appear in two dropdowns:
- **Coloring & Filtering → Categorical obs**
- **Coloring & Filtering → Continuous obs**

### Var (variables) = per-gene (feature) values per cell

Var fields are feature-level “things you can color by”, most commonly **genes**.

In the UI, var fields appear as:
- **Coloring & Filtering → Gene Expression** (search box + dropdown results)

:::{important}
Cellucid treats gene expression as **continuous**: selecting a gene always creates a continuous-style legend (colorbar + sliders).
:::

---

## Categorical vs continuous (what changes in the app)

### Categorical fields

Categorical fields have:
- a finite set of **category labels** (strings), and
- a per-cell **category code** (an integer index into the category list).

In the legend you can typically:
- show/hide categories (checkboxes),
- change category colors (swatches),
- rename categories,
- merge or delete categories (creating derived fields).

Categorical fields can also expose **centroids** (category centers) if present or computed.

<!-- SCREENSHOT PLACEHOLDER
ID: fields-types-categorical-legend-example
Suggested filename: web_app/fields_legends/02_categorical-legend-example.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 01_field_types_and_sources.md
Capture:
  - Select a categorical obs field with at least ~8 categories (e.g., clusters)
  - Ensure the legend shows: checkboxes, swatches, category names, and counts
  - Optional: enable “Show centroids” and “Show labels” so readers see what centroids mean
Crop:
  - Include: Display options box + a small portion of the point cloud
Annotations:
  - Call out: category checkbox, color swatch, category label, per-category cell count, (optional) centroid toggles
Alt text:
  - Categorical legend with per-category checkboxes and color swatches.
Caption:
  - State that categorical legends are for discrete labels and support category-level editing and visibility filtering.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a categorical legend example.
:width: 100%

Categorical fields produce a checkbox-based legend where you can hide/show categories and edit category colors and labels.
```

### Continuous fields

Continuous fields have:
- a per-cell numeric value (float-like),
- a color mapping (a **colormap**), and
- a visible numeric range filter (min/max).

In the legend you can typically:
- change the colormap (palette),
- optionally enable **log scale** for coloring,
- adjust the **visible range** with sliders (filtering),
- optionally rescale the colorbar to match your slider range.

<!-- SCREENSHOT PLACEHOLDER
ID: fields-types-continuous-legend-example
Suggested filename: web_app/fields_legends/03_continuous-legend-example.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 01_field_types_and_sources.md
Capture:
  - Select a continuous obs field with a wide dynamic range (or select a gene)
  - Ensure the legend shows: color scale header, gradient bar, log toggle, rescale toggle, min/max sliders
  - Optional: open the colormap menu (so readers see it exists)
Crop:
  - Include: Display options box + a small portion of the colored point cloud
Annotations:
  - Call out: gradient bar (colormap picker), colorbar min/max labels, log toggle, rescale toggle, range sliders
Alt text:
  - Continuous legend showing a colorbar, palette picker, and min/max range sliders.
Caption:
  - Explain that continuous legends control both coloring (palette/scale) and filtering (visible numeric range).
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a continuous legend example.
:width: 100%

Continuous fields produce a colorbar-based legend where you control the palette, scale, and visible numeric range.
```

---

## “Field sources” (where fields come from)

From a user perspective, a field can come from:

1) **Your dataset export** (the “source-of-truth” columns and genes you provided).
2) **A derived field created inside Cellucid** (user-defined fields), for example:
   - duplicating a field (“(copy)”),
   - deleting/merging categories (“(edited)” / “(merged)”),
   - building a categorical field from highlight pages (Category Builder).
3) **A renamed or restored field**, where the UI name differs from the original key used for loading.

This matters because derived fields are still “real fields”:
- they appear in selectors,
- can be filtered and saved in session bundles,
- and can be deleted/restored like anything else.

---

## Edge cases you should expect (and how Cellucid behaves)

### Missing values

- Continuous: `NaN`/missing values are shown as a neutral “None” gray, and they are **excluded** when you apply a numeric range filter.
- Categorical: invalid/missing category codes show as a neutral “None” gray; they typically do not appear as a named category in the legend.

### All-constant or all-missing fields

- If a continuous field has no finite values, Cellucid still renders a safe default range (so the UI doesn’t break).
- If all values are identical, the UI will still show a non-zero range so the colormap can render.

### Negative values + log scale

Log scale coloring requires positive values.
If you enable log scale on a field with zeros/negatives:
- those cells remain visible, but they use the neutral “None” gray color in the viewer.

---

## Next steps

- To learn how selection, rename, delete, and restore behave in detail, go to `d_fields_coloring_legends/02_field_selector_ux`.
- To understand exactly what the legend controls do (including filters), go to `d_fields_coloring_legends/04_legend_behavior`.
- If something looks wrong, jump to `d_fields_coloring_legends/05_troubleshooting_fields_legends`.
