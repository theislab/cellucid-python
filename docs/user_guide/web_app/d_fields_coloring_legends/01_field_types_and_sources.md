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

The var axis is whatever the export published — one name per gene, decided when
the dataset was prepared. It is not necessarily every gene the study measured,
and Cellucid neither knows nor can recover what was left out. For the built-in
samples this is documented: see {doc}`07_genes_in_the_built_in_samples`.

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

```{figure} ../../../_static/screenshots/web_app/legend-categorical-clusters.png
:alt: A legend with the hint Click a swatch to pick a color, Show All and Hide All buttons, and eight rows in alphabetical order each holding a checkbox, a coloured square swatch, a category name, a cell count, a pencil icon and a bin icon.
:width: 428px

A categorical legend. **Checkboxes and swatches** are the tell — if you see
these, the active field is categorical.
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

```{figure} ../../../_static/screenshots/web_app/legend-continuous-sscore.png
:alt: A legend headed COLOR SCALE (VIRIDIS) above a purple to yellow gradient bar labelled minus 0.37 and 1.14 at its ends, a Log color scale toggle reading Off, the note Zero negative and NaN values use the None color, a Rescale colorbar to slider range toggle reading On, the note On colors track sliders Off full data range, then a FILTERING heading, a VISIBLE RANGE subtitle, a Live filtering toggle reading On, Min and Max sliders with numeric readouts, and FILTER and RESET buttons.
:width: 428px

A continuous legend for the same dataset. **A colour bar with Min/Max sliders**
is the tell. The same legend appears for a gene — a gene *is* a continuous
field, just from the var axis instead of the obs axis.
```

:::{tip}
Two one-line notes the legend prints for you, which answer questions this page
would otherwise need paragraphs for:

- `Zero/negative/NaN values use the None color`
- `On: colors track sliders; Off: full data range`
:::

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
- Constant fields and constant genes are published, not dropped: the export
  format encodes them as `minValue == maxValue` with every code `0`, and the
  viewer decodes that back to the exact constant. A gene detected in no cell of
  a published subset is therefore still in the gene menu, at its exact value.

### Negative values + log scale

Log scale coloring requires positive values.
If you enable log scale on a field with zeros/negatives:
- those cells remain visible, but they use the neutral “None” gray color in the viewer.

---

## Next steps

- To learn how selection, rename, delete, and restore behave in detail, go to {doc}`02_field_selector_ux`.
- To understand exactly what the legend controls do (including filters), go to {doc}`04_legend_behavior`.
- If something looks wrong, jump to {doc}`05_troubleshooting_fields_legends`.
