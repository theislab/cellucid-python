# Outlier filtering (per active field)

**Audience:** everyone (wet lab + computational users)  
**Time:** 10–15 minutes  
**What you’ll learn:**
- What the outlier filter removes (and what it does *not* remove)
- How to turn it on: choosing the latent space it measures in
- Where its numbers come from: stored in the file, or computed once in the browser
- What the line under the slider tells you, and how to cancel a long computation

**Prerequisites:**
- A dataset loaded
- A categorical obs field selected (clusters, cell types)

---

## What this filter is (and is not)

For a categorical field, each cell gets a number in `[0, 1]`: how far it lies
from the centre of its category **in a latent space** — an `obsm` array such as
`X_scVI`, `X_scANVI` or `X_pca` — expressed as the share of its category that
lies at that distance or closer. The rule, for every category with at least
`centroid_min_points` cells (10 by default):

- the category's centre is the mean of its cells' latent coordinates;
- a cell's distance is its Euclidean distance to that centre;
- its quantile is the number of the category's cells at that distance or
  closer, divided by the category's size.

Interpretation:
- `0.50` means “around the median distance in this category”
- `0.95` means “in the farthest ~5% of this category”

The slider hides every cell whose quantile is above its value.

This filter is **not**:
- measured in the 2D/3D layout you see (UMAP, t-SNE): the latent space is
  where the biology is, the layout only draws it;
- a min/max range filter on a continuous field;
- a “global outlier” score across the entire dataset: it is relative to each
  category, so every category loses its own farthest share.

---

## Fast path (wet lab / non-technical)

Use this when your clusters/regions look “fuzzy” and you want to hide fringe points.

1) Pick a categorical field (clusters/cell types) in **Coloring & Filtering → Categorical obs**.
2) Look in the **Legend** box for **Outlier filter**. The line under it says
   what it can do for this field:
   - nothing under the slider: it works at once;
   - *“Choose the latent space in Data mapping”*: press **Choose latent
     space** and pick one (see [Turning it on](#turning-it-on-choose-a-latent-space));
   - *“The first use reads the latent space … once”*: the file stores no
     outlier distances for this space; moving the slider computes them (see
     [Stored or computed](#stored-or-computed)).
3) Drag the slider from `100%` down to `95%`.

```{figure} ../../../_static/screenshots/filtering/outlier-slider-95.png
:alt: A control labelled Outlier filter whose slider handle sits a short distance in from the right end, with the number 95 percent printed to its right and the mouse pointer resting on the handle.
:width: 436px

The control lives in the **Legend** box, below **Show centroids**. It starts
at `100%`, which means no outlier filtering at all.
```

4) Confirm in **Active filters** that you see a line like:
   - `<field>: outlier ≤ 95%`
5) If you removed too much:
   - increase the slider (e.g., `98%`), or reset to `100%`.

### What success looks like

```{figure} ../../../_static/screenshots/filtering/window-outlier-100-off.png
:alt: The whole application window coloured by clusters with the outlier slider at its right end reading 100 percent, loose single points scattered in the space around and between the coloured groups, and the sidebar reading Showing all 3,696 cells above No filters active.
:width: 1440px

**At `100%`** — no filtering. Note the loose points drifting between the
coloured groups, and the legend counts: Alpha 481, Ductal 916, Pre-endocrine
592.
```

```{figure} ../../../_static/screenshots/filtering/window-outlier-95-applied.png
:alt: The same window with the slider moved left and reading 95 percent, the mouse pointer on it, visibly fewer loose points between the coloured groups, every legend count reduced to Alpha 458, Ductal 872 and Pre-endocrine 564, and the sidebar reading Showing 3,518 of 3,696 cells above a row reading clusters colon outlier less than or equal to 95 percent.
:width: 1440px

**At `95%`** — the same camera, the same colouring, 178 fewer cells drawn.
Every category has lost roughly its own top 5%: Alpha 481 → 458, Ductal
916 → 872, Pre-endocrine 592 → 564. That per-category proportionality is the
signature of this filter, and it is what distinguishes it from a global
threshold.
```

So: the “halo” of far-away cells diminishes, the main structure becomes
clearer, and “Showing X of Y cells” drops by a few percent — not by half.

---

## Turning it on: choose a latent space

The filter needs a latent space, and Cellucid never picks one for you — a
2D/3D layout would give misleading distances, and a file may hold several
latent spaces. It is chosen in one of three ways:

- the file names it: `cellucid.prepare(..., latent_key="X_scVI")` in Python or
  `cellucid_prepare(..., latent_key = "X_scVI")` in R writes it as
  `uns['cellucid']['defaults']['latentKey']`;
- the launch names it: `cellucid serve data.h5ad --latent-key X_scVI`,
  `serve_anndata(..., latent_key="X_scVI")` or `show_anndata(..., latent_key="X_scVI")`;
- you choose it under **Session › Data mapping › Latent space** (“Outliers in”;
  see [Latent space](../b_data_loading/13_data_mapping.md#latent-space)).

Without one, the outlier control is shown **disabled** for every categorical
field, with what the filter measures, the file's arrays that look like a
latent space (more than three columns, e.g. `X_scVI or X_pca`) and a
**Choose latent space** button that opens that setting. The suggestions are
only suggestions.

---

## Stored or computed

The quantiles of every categorical field need the whole latent space: on a
2-million-cell atlas with a 30-dimensional scVI space that is 240 MB of values.
Cellucid therefore gets them in one of two ways, with the same bytes either
way:

- **Stored.** `prepare()` and `cellucid_prepare()` compute the quantiles of
  every categorical obs field in the latent space they are given and store
  them in the file (`uns['cellucid']['outlierQuantiles']`, one float32 array
  per field, about 1.8 MB per field for 562,000 cells). The browser reads only
  the field you filter, on the slider's first move: one small read, no line
  under the slider, and the filter works at once. The Data mapping list marks
  that latent space **outliers stored**.
- **Computed.** Any other file (an AnnData you wrote yourself, a CELLxGENE
  download), or another latent space than the stored one: before anything is
  read, the line under the slider says what the first use costs — *“The first
  use reads the latent space X_scVI (54 MB) once to compute the outlier
  distances; every other field then filters at once.”* Moving the slider below
  `100%` starts it: the latent space is read once, in the background, every
  field asked for at the same time is computed in the same pass, and the
  latent space is kept in memory (within Cellucid's memory budget), so every
  other field filters without reading it again. A latent space larger than
  that budget is read again for each field, and the line says so.

While the distances are read or computed, a bar under the slider and a
progress notification show how much of how much is read, how fast and how
long is left (also in the read log, `[cellucid:read]` in the browser
console). **Cancel** — under the slider or on the notification — stops the
read at once, keeps nothing, and puts the slider back at `100%`; moving it
again starts over. The slider stays usable while you wait: the filter applies
at the value it holds when the distances arrive.

Stored quantiles are used only when they still fit the file: written for the
latent space in use and the same `centroid_min_points`, and matching the
field's categories as the file holds them now. A file whose obs was edited
after it was prepared (categories renamed, cells subset) gets its quantiles
computed instead, and the console says why.

A latent space holding a missing value (NaN) or a value beyond float32 cannot
be measured in: computing stops at that cell, the line under the slider says
so, and the latent space is marked unusable in Data mapping.

---

## Practical path (computational users): exact semantics

### What the percent means

The slider is a **threshold** in `[0, 1]` displayed as a percent.

- `100%` means **no outlier filtering**.
- `p%` means: **hide cells whose outlier quantile is above `p/100`.**

Example:
- `95%` hides approximately the top 5% most outlier-like cells *within each category*, for fields where the outlier quantiles were computed that way.

### When a row appears in Active filters

The outlier row is written only when the threshold is **below `99.99%`**. A
slider parked at `100%` produces no row, which is correct — it is filtering
nothing — but it is also why “the slider does nothing” is such a common
report.

The row also carries **no `visible / available` cell count**, unlike category
and range rows. Its only numeric effect is on the `Showing X of Y cells` line.

### What happens with missing outlier quantiles

A cell whose category has fewer than `centroid_min_points` cells (default
`10`), or that has no category, gets no quantile, and is **never removed by
outlier filtering**. The web app holds it as exactly `-1`, so its visibility
loop can tell “none” from the valid `[0, 1]` interval without carrying a
non-finite value into a hot loop (the file stores NaN).

When *no* category of a field reaches `centroid_min_points`, the control is
shown disabled with *“No category of “field” has 10 cells, the fewest the
filter measures”*. Lower `centroid_min_points` (`prepare()`,
`cellucid_prepare()`, `serve_anndata()`, `show_anndata()`) to measure smaller
categories.

### When the outlier control appears

The outlier control appears for every **categorical obs field of the file**;
the line under it says whether it can filter and why not.

- It never appears for a **continuous** field or a **gene**: outlier quantiles
  are a per-category statistic.
- A field you create by merging or deleting categories has no outlier
  quantiles, so the control is hidden for it.
- Switching to another field stops applying the first field's outlier filter;
  its threshold is kept for when you come back.

### Sessions

A saved session keeps every field's outlier threshold. Restoring it gets the
quantiles the thresholds need — read from the file, or computed with the same
progress and Cancel. If they cannot be got (you cancel, no latent space is
chosen, the latent space cannot be read), the rest of the session is restored
and a notification names each outlier threshold that was not, with the
reason: a threshold is never dropped silently.

### Changing the latent space

Choosing another latent space reads nothing when you press **Apply**: every
threshold is kept (a threshold is a share of a category, so it means the same
in the new space), the old space's quantiles are dropped, and a field whose
outlier filter is in use gets its new quantiles at once, with progress and
Cancel. Choosing **None** turns outlier filters off until a latent space is
chosen again.

---

## Troubleshooting (quick)

### “I don’t see the outlier slider”

- The active field is continuous, a gene, or a field you created: outlier
  quantiles are per category of a field the file holds.

### “The slider is greyed out”

Read the line under it:
- *“Choose the latent space in Data mapping”* — press **Choose latent space**;
- *“No category … has 10 cells”* — the field's categories are too small to
  measure; lower `centroid_min_points`.

### “Moving the slider started a download”

The file stores no outlier distances for this latent space, so the first use
reads it once (the line under the slider said how much). Press **Cancel** to
stop; to avoid it, write the file with `prepare(latent_key=...)` or
`cellucid_prepare(latent_key = ...)`, which store them.

### “The slider does nothing”

Common causes:
- it’s still at `100%`,
- the distances are still being got (see the bar under the slider),
- your categories are small and many cells have no quantile (they are never filtered),
- the dataset has very tight clusters so few cells exceed `95%`.

### “It removes almost everything”

Common causes:
- threshold is too low (try `98–99%`),
- the latent space does not match the structure you are looking at: choose another under **Data mapping › Latent space**.

---

## Next steps

- {doc}`03_filter_stack_ui_active_filters` (how to disable/remove the outlier filter)
- {doc}`06_edge_cases_filtering` (what happens when all cells are filtered out)
