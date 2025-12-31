# Troubleshooting (fields & legends)

This page is a symptom-driven guide for problems with:
- field lists (obs/gene expression),
- coloring (categorical/continuous),
- legends (checkboxes/sliders),
- and “why did points disappear?” filter interactions.

If you’re new to the concepts, read these first:
- `d_fields_coloring_legends/01_field_types_and_sources`
- `d_fields_coloring_legends/03_color_by_behavior`
- `d_fields_coloring_legends/04_legend_behavior`

---

## First-response checklist (fast)

Before deep debugging, check these three things:

1) **Active field**: In **Coloring & Filtering**, is the correct selector active?
   - Categorical obs vs Continuous obs vs Gene Expression
2) **Legend toggles**:
   - Is **Use log scale** on (and your data has zeros/negatives)?
   - Is **Rescale colorbar to slider range** on (changing contrast)?
3) **Active filters**:
   - Do you have hidden categories or a narrowed numeric range from earlier?

Most “bugs” are one of these.

---

## Symptom: “The legend is missing / Display options is empty”

### Likely causes (ordered)

1) No active field is selected (all selectors are on **None**).
2) You’re editing a different view than you think (multiview active view mismatch).
3) The field failed to load (error notification, console error).

### How to confirm

- In **Coloring & Filtering**, check:
  - Categorical obs dropdown value
  - Continuous obs dropdown value
  - Gene Expression input (is a gene selected?)
- If you are in multiview:
  - click a view panel, then re-check the sidebar (legend follows active view)
- Open DevTools → Console and look for messages like:
  - `Failed to load field: ...`
  - `Failed to load gene: ...`

### Fix

1) Select a field (categorical/continuous) or a gene.
2) If selection fails, reload the dataset and try again.
3) If the error persists, jump to:
   - `b_data_loading/04_troubleshooting_data_loading`

### Prevention

- When working in multiview, use **Edit selected view** when making legend changes.

---

## Symptom: “Everything is gray”

### Likely causes (ordered)

1) No active field selected (Cellucid renders the “no-field” state as neutral gray).
2) Log scale enabled on data that is mostly zero/negative.
3) The field is all missing (`NaN`) or did not load correctly.

### How to confirm

- If all selectors show **None** → it’s (1).
- If **Use log scale** is on and you’re coloring by a sparse gene → it’s often (2).
- If you see a load failure notification → it’s (3).

### Fix

- For (1): pick a categorical/continuous field or gene.
- For (2): turn off **Use log scale**, or choose a gene/field with positive values.
- For (3): reload the dataset; if it persists, inspect export integrity.

### Prevention

- Treat “gray” as a meaning-bearing state:
  - gray can mean “no field selected” or “value is missing / log-incompatible”.

---

## Symptom: “The field list is empty / I can’t find my expected column”

### Likely causes (ordered)

1) The dataset export did not include that obs column.
2) The column exists but is classified differently than you expect:
   - categorical vs continuous split is decided during export.
3) The field was deleted (soft-deleted) and is sitting in **Deleted Fields**.
4) You’re looking at a different dataset (dataset identity mismatch).

### How to confirm

- Check whether the dropdown shows:
  - `(no categorical obs fields)` or `(no continuous obs fields)`
- Scroll to the bottom of **Coloring & Filtering** and look for **Deleted Fields**.
- If you recently edited/merged categories:
  - the original field may have been moved to Deleted Fields and replaced by a derived field (e.g., `clusters (merged)`).

### Fix

1) If it’s in **Deleted Fields** → click **Restore**.
2) If it’s missing from the export:
   - re-export with the column included (cellucid-python `prepare()` / your export pipeline).
3) If it landed in the “wrong” dropdown:
   - treat it as that type for now (you can still use it), and adjust export typing later.

### Prevention

- Establish a naming convention for derived fields (e.g., suffixes like `(cleaned)`), and avoid deleting original fields until you have a session saved.

---

## Symptom: “Gene Expression is missing / the gene search box isn’t there”

### Likely causes (ordered)

1) The dataset has no gene expression (var fields) available in the export.
2) The dataset loaded, but var metadata failed to load.

### How to confirm

- If the **Gene Expression** row is not visible at all, the dataset likely has no var fields.
- Check dataset info (if available) for gene count.
- DevTools → Network:
  - look for `var_manifest.json` (or equivalent) requests failing.

### Fix

- Use obs fields instead (categorical/continuous).
- If you expected genes, re-export with gene expression enabled.

### Prevention

- For very large gene matrices, prefer server-backed or lazy-loading workflows (server mode / Jupyter integration) so the app can load genes on demand efficiently.

---

## Symptom: “Gene search returns nothing / Enter selects the wrong gene”

### Likely causes (ordered)

1) You’re searching with a gene naming convention that isn’t in the dataset (symbols vs Ensembl IDs).
2) The search is substring-based; a short query matches many genes and Enter picks the top one.
3) The gene list is present but very large; your query is too broad.

### How to confirm

- Try a longer, more specific substring.
- If results show `...and N more`, you’re matching too many.
- Try searching by an alternate identifier style if you have it (e.g., `ENSG...` vs `MS4A1`).

### Fix

1) Type a more specific query and click the correct gene explicitly (don’t rely on Enter).
2) If you don’t know the identifier convention:
   - check the dataset’s gene IDs in your preprocessing notebook/export code.

### Prevention

- Standardize gene identifiers in your export and document it for collaborators (symbols vs Ensembl).

---

## Symptom: “Points disappeared / the dataset looks empty”

### Likely causes (ordered)

1) You applied filters earlier and forgot (hidden categories, continuous min/max range, outlier filtering).
2) You hid almost all categories (or clicked Hide All).
3) Continuous filter range is extremely narrow (or inverted).
4) Outlier filter is set below 100% and is hiding outliers.

### How to confirm

- Check **Active filters**:
  - if it lists anything other than “No filters active”, filters are applied.
- For categorical fields:
  - look for unchecked categories.
- For continuous fields:
  - look at Min/Max sliders and the numeric readouts.
- For outliers:
  - check **Outlier filter (latent space)** percent (if visible).

### Fix

1) For categorical filters:
   - click **Show All** in the categorical legend.
2) For continuous filters:
   - click **RESET** in the continuous legend filtering section.
3) For outlier filtering:
   - set the outlier slider back to `100%` (or switch to a field without outlier filtering).

### Prevention

- Watch **Active filters** as your “single source of truth” for why points are missing.

---

## Symptom: “A category is greyed out / checkbox is disabled”

### Likely causes (ordered)

1) There are **0 available cells** in that category after other filters.

### How to confirm

- Hover the greyed-out row; it typically indicates “No cells available in this category after other filters”.
- The count may read `0 cells`.

### Fix

- Clear the filter(s) that removed those cells:
  - categorical: Show All on the other field(s)
  - continuous: RESET range filters

### Prevention

- When filtering with multiple fields, remember that a category can become empty even if it exists in the dataset.

---

## Symptom: “Colors look washed out / everything is the same color”

### Likely causes (ordered)

1) You filtered to a narrow range but kept the color domain wide (rescale is off).
2) The field is nearly constant (very low dynamic range).
3) You’re coloring by a quantized/low-bit-depth export (banding/steps).

### How to confirm

- If you have a narrow Min/Max slider range:
  - check whether **Rescale colorbar to slider range** is Off.
- Toggle rescale On and see if contrast improves.

### Fix

1) Turn **Rescale colorbar to slider range** On.
2) Choose a more appropriate colormap (Viridis/Cividis are good defaults).
3) If values are quantized and you need more precision:
   - re-export with higher precision settings (cellucid-python export knobs).

### Prevention

- For publication-like figures, prefer rescale On when showing a restricted subset.
- For cross-view comparisons, prefer rescale Off for consistent mapping.

---

## Symptom: “Log scale makes everything gray”

### Likely causes (ordered)

1) The field has mostly zeros (common for sparse gene expression).
2) The field has negative values (some normalized/scaled outputs).
3) The field has no positive values at all.

### How to confirm

- Turn log scale off: if colors return, log incompatibility is the cause.

### Fix

- Use log scale only on fields with meaningful positive values.
- For sparse genes:
  - log scale can be useful to highlight the non-zero tail, but expect many zeros to remain gray.

### Prevention

- If you want log-like visualization, export log1p values (common in scRNA workflows) and keep log scale off in the viewer.

---

## Symptom: “Rename/Delete is disabled (field or category)”

### Likely causes (ordered)

1) Community Annotation voting mode is enabled for that categorical field (labels are locked).
2) A gene is loading and the UI is temporarily disabled (busy state).

### How to confirm

- If you see `🗳️` next to a field name in the categorical dropdown, voting mode is enabled.
- If you just clicked a gene and the app is loading, controls may be disabled briefly.

### Fix

- For voting mode:
  - disable voting mode (repo authors only) or use a duplicate field for local edits.
- For loading/busy:
  - wait for the load to complete; then retry.

### Prevention

- Duplicate a field before doing heavy category edits, especially when collaborating via annotation workflows.

---

## Symptom: “Merging/deleting a category created a new field / my original column vanished”

### Likely causes (ordered)

1) You edited a **source** categorical field; Cellucid created a derived field and soft-deleted the original for safety.

### How to confirm

- Look at the dropdown: you may now be on a field named like:
  - `<original> (edited)` or `<original> (merged)`
- Check **Deleted Fields**: the original is likely there and restorable.

### Fix

- If you didn’t mean to edit the source:
  - restore the original field from Deleted Fields and continue from that.
- If you did mean it:
  - keep working on the derived field (subsequent edits may happen in place).

### Prevention

- Treat derived-field suffixes as a history trail (great for reproducibility).

---

## Symptom: “I restored a field and its name changed”

### Likely causes (ordered)

1) The restored field name conflicted with an existing visible field, so Cellucid auto-renamed it (e.g., `(... restored)`).

### How to confirm

- Look for a notification like `Restored "X" as "X (restored)"`.

### Fix

- Rename the restored field to a clearer name (unique keys are required).

### Prevention

- Avoid creating multiple derived fields with identical names; use meaningful suffixes.

---

## Symptom: “I can’t restore a deleted field”

### Likely causes (ordered)

1) You clicked **Confirm** in Deleted Fields (purged restore capability).

### How to confirm

- The Deleted Fields row may no longer appear, or restore may do nothing.
- Session restores will also not bring it back if it was confirmed.

### Fix

- There is no in-app undo for confirmed deletion.
- Re-load/re-export the dataset if you need the original column again.

### Prevention

- Only use **Confirm deletion** when you are truly done with a field.

---

## Symptom: “Gene loading is slow / selecting a gene freezes the UI”

### Likely causes (ordered)

1) Large dataset + slow storage/network (remote export host).
2) Browser memory pressure (many genes loaded in one session).
3) You’re not using server-backed / lazy workflows for large matrices.

### How to confirm

- DevTools → Network: gene value requests are slow or failing.
- DevTools → Performance/Memory: memory grows as you load many genes.

### Fix

- Prefer server-backed or local server workflows for very large gene matrices.
- Avoid loading hundreds of genes in one session unless necessary.
- If the browser becomes unstable, reload and use a smaller set of genes.

### Prevention

- For interactive exploration, export the specific genes/scores you care about as obs fields when possible.

---

## Symptom: “Legend changes when I click different panels (multiview confusion)”

### Likely causes (ordered)

1) The active view changed; the sidebar always shows controls for the active view.

### How to confirm

- Click a different panel; observe the view badges (active view indicator) and the legend updating.

### Fix

- Use **Edit selected view** layout when you are configuring fields/legends.
- Then switch to **Grid compare** to compare.

### Prevention

- Keep cameras locked until you’re comfortable; it reduces “which panel am I editing?” confusion.

---

## Still stuck?

If you can reproduce the issue, capture:
- a screenshot of **Coloring & Filtering** (showing selectors + legend + active filters),
- the browser console output around the failure,
- and your dataset export metadata (obs/var manifests).

Then continue debugging from:
- `q_troubleshooting_index/index` (broader issues)
