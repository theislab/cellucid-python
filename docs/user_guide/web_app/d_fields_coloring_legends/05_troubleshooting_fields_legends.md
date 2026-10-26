# Troubleshooting (fields & legends)

This page is a symptom-driven guide for problems with:
- field lists (obs/gene expression),
- coloring (categorical/continuous),
- legends (checkboxes/sliders),
- and “why did points disappear?” filter interactions.

If you’re new to the concepts, read these first:
- {doc}`01_field_types_and_sources`
- {doc}`03_color_by_behavior`
- {doc}`04_legend_behavior`

---

## First-response checklist (fast)

Before deep debugging, check these three things:

1) **Active field**: In **Coloring & filtering**, is the correct selector active?
   - Categorical obs vs Continuous obs vs Gene expression
2) **Legend toggles**:
   - Is **Log color scale** on (and your data has zeros/negatives)?
   - Is **Rescale to range** on (changing contrast)?
3) **Filters**:
   - Do you have hidden categories or a narrowed numeric range from earlier?

Most “bugs” are one of these.

---

## Symptom: “The legend is missing / the Legend section is empty”

### Likely causes (ordered)

1) No active field is selected (all selectors are on **None**).
2) You’re editing a different view than you think (multiview active view mismatch).
3) The field failed to load (error notification, console error).

### How to confirm

- In **Coloring & filtering**, check:
  - Categorical obs dropdown value
  - Continuous obs dropdown value
  - Gene expression input (is a gene selected?)
- If you are in multiview:
  - click a view panel, then re-check the sidebar (legend follows active view)
- Read the notification centre before the console. A selection that failed says
  so there, and says why:
  - a gene raises a multi-line notification titled `Gene cannot be shown`,
    naming the gene, counting the offending values and giving the first affected
    cell indices;
  - a field raises the same text after the prefix `Failed to load field:`.
- The console carries the thrown error object for the same failure, which is
  what to copy into a bug report — not a different sentence to search for.

### Fix

1) Select a field (categorical/continuous) or a gene.
2) If the notification names non-finite values, this is not a transport problem:
   see {ref}`undrawable-continuous-payload` below.
3) If selection fails with nothing else to go on, reload the dataset and try again.
4) If the error persists, jump to:
   - {doc}`../b_data_loading/08_troubleshooting_data_loading`

### Prevention

- When working in multiview, use **Edit selected view** when making legend changes.

---

## Symptom: “Everything is gray”

### Likely causes (ordered)

1) No active field selected (Cellucid renders the “no-field” state as neutral gray).
2) **Log color scale** enabled on a sparse field, where most cells are zero and
   therefore drop to the “None” grey.
3) The gene or column you picked **cannot be drawn** — it contains an infinity,
   or it has no value at all — so it was refused and the active field never
   changed. What you are looking at is still the previous (or no-field) grey.
4) The field did not load correctly (transport failure).

### How to confirm

- If all selectors show **None** → it’s (1).
- If **Log color scale** is on and you’re coloring by a sparse gene → it’s often (2).
- If a notification titled `Gene cannot be shown` (or one beginning
  `Failed to load field:`) counts `+Infinity` / `-Infinity` / `NaN` values and
  lists the first affected cells → it’s (3). Note that the selector snaps back:
  the field was never activated.
- If you see a load failure that names a status code or a URL instead → it’s (4).

### Fix

- For (1): pick a categorical/continuous field or gene.
- For (2): turn off **Log color scale**, or choose a gene/field with positive values.
- For (3): repair the values at the source — reloading cannot change the
  outcome. See {ref}`undrawable-continuous-payload`.
- For (4): reload the dataset; if it persists, check the column in the AnnData file itself.

:::{note}
“The field has **no** positive values at all” is not a possible cause here:
in that case the toggle refuses to turn on
(`Log color scale requires at least one positive field value.`) and the colours
never change. Only a field with *some* positives can produce a mostly-grey
picture this way.
:::

### Prevention

- Treat “gray” as a meaning-bearing state:
  - gray can mean “no field selected” or “value is missing / log-incompatible”.

---

(undrawable-continuous-payload)=
## Symptom: “Gene cannot be shown” / a field that will not colour

### Likely causes (ordered)

1) The gene or continuous column contains `+Infinity` or `-Infinity`. An
   infinity has no position on a colour scale, so the payload is refused.
2) Every value in the payload is `NaN`. There is no range to scale, so it is
   refused too — usually a column that was never computed, or an empty
   expression layer served instead of the normalized one.
3) The values are outside what `float32` can carry (an overflow or an
   underflow). Only the writers report this kind; it is refused for the same
   reason.

### How to confirm

Read the notification. It answers each question in turn:

- **A gene** raises a multi-line notification titled `Gene cannot be shown`.
- **A continuous obs field** raises the same body after the prefix
  `Failed to load field:`.

For an infinity the body reads like:

> Gene "GATA3" contains 3 infinite values, so it has no colour scale.  
> Of 18,142,044 cells: 3 +Infinity, 12 NaN.  
> First affected cells: 4,118,237, 9,002,551, 17,441,908.  
> NaN is drawn as missing; an infinity is not a value on any scale.  
> Repair the expression matrix before serving or exporting it — …

For a payload with nothing in it:

> Field "leiden_score" has no value in any of 18,142,044 cells, so it has no colour scale.  
> Every value is NaN, which Cellucid reads as "not measured here".

A file served by `cellucid serve`, `serve_anndata` or `show_anndata` is read
by the viewer itself, so its gene or column is refused with these same
messages.

### Fix

Repair the values at the source. Reloading, re-picking the gene or clearing the
cache cannot change the outcome — nothing in the viewer can place an infinity on
a colour scale.

1) For an expression matrix:
   - sparse: `adata.X.data[~np.isfinite(adata.X.data)] = 0`
   - dense: `np.nan_to_num(adata.X, copy=False)`
2) For a continuous obs column:
   - `adata.obs["<key>"] = adata.obs["<key>"].replace([np.inf, -np.inf], np.nan)`
3) Then re-run `prepare()` (or save the file again and restart
   `cellucid serve` / `serve_anndata` / `show_anndata`) and reload the dataset.

If you meant the infinities as “missing”, replace them with `NaN` rather than
with `0`: Cellucid draws `NaN` as the neutral grey meaning “not measured here”,
whereas `0` is a real measurement and will be coloured as one.

### Prevention

- Check a column before you publish it: `np.isfinite(values).all()`.
- A file written by `prepare()` or `cellucid_prepare()` cannot carry this
  problem — both writers refuse an infinity in a continuous obs field, and any
  non-finite value in gene expression, with the same counts and the same
  first-affected indices. If you hit this at all, the file was written another
  way.

---

## Symptom: “The field list is empty / I can’t find my expected column”

### Likely causes (ordered)

1) The AnnData file does not hold that obs column.
2) The column exists but is classified differently than you expect:
   - categorical vs continuous follows the column's dtype in the file (`category` vs numeric).
3) The field was deleted (soft-deleted) and is sitting in **Deleted fields**.
4) You’re looking at a different dataset (dataset identity mismatch).

### How to confirm

- Check whether the dropdown shows:
  - `(no categorical obs fields)` or `(no continuous obs fields)`
- Scroll to the bottom of **Coloring & filtering** and look for **Deleted fields**.
- If you recently edited/merged categories:
  - the original field may have been moved to Deleted fields and replaced by a derived field (e.g., `clusters (merged)`).

### Fix

1) If it’s in **Deleted fields** → click **Restore**.
2) If it’s missing from the file:
   - check **Session › Data mapping**, which may have turned the column off;
   - otherwise add the column to `obs` and save the file again (or write it with
     `prepare(obs=...)`).
3) If it landed in the “wrong” dropdown:
   - treat it as that type for now (you can still use it), and change the
     column's dtype in `obs` later (for example `.astype("category")`).

### Prevention

- Establish a naming convention for derived fields (e.g., suffixes like `(cleaned)`), and avoid deleting original fields until you have a session saved.

---

## Symptom: “Gene expression is missing / the gene search box isn’t there”

### Likely causes (ordered)

1) The dataset file has no expression matrix (`X`, `raw/X` or a layer).
2) **Session › Data mapping** has no expression matrix chosen.

### How to confirm

- If the **Gene expression** row is not visible at all, the dataset likely has no var fields.
- Check the gene count in the **Dataset** box at the top of **Session**.
- Open **Session › Data mapping › Genes** and check which
  matrix is chosen.

### Fix

- Use obs fields instead (categorical/continuous).
- If you expected genes, write the file with its expression matrix
  (`prepare(gene_expression=...)`, or keep `X` when saving the AnnData).

### Prevention

- For very large gene matrices, store them column-major (CSC): each gene is then one contiguous read, with no gene index to build first.

---

## Symptom: “Gene search returns nothing / Enter selects the wrong gene”

### Likely causes (ordered)

1) The gene is published under a different name than the one you typed.
2) The gene is not in this dataset at all — the file holds only some of the
   source’s genes.
3) The search matches anywhere in a name; a short query matches many genes.
   The list puts the gene whose name is exactly what you typed first, then
   names that start with it, then the rest, and Enter takes the top one unless
   you first move the highlight with the arrow keys.
4) The gene list is present but very large; your query is too broad.

### How to confirm

- Try a longer, more specific substring. Matching ignores letter case, so
  `gata3` finds `GATA3`.
- If the list keeps growing as you scroll (every match is listed, 100 at a time), you’re matching too many.
- Open the field selector and read a few gene names directly. Search matches
  the gene names the dataset lists, and there is no alias mapping.
- Compare the count in the empty-state message (“This dataset publishes *N* gene
  names…”) against the gene count you expect from the source. A large gap means
  preparation published a subset.

### If you are using a built-in sample

The four human samples — Suo, Garcia, He, and Kanemaru — publish only the genes
whose Ensembl accession resolved to a symbol, which is between 45.1% and 70.2%
of the gene axis of the file each was built from. A gene that is genuinely in
the study can therefore be absent, under any spelling. The rule, the per-sample
counts, and what to do instead are in
{doc}`07_genes_in_the_built_in_samples`.

### Fix

1) Type a more specific query, then click the correct gene, or move to it with
   ArrowDown/ArrowUp and press Enter.
2) If the name you expected is not there, it is either published under another
   name or not published at all:
   - check which `var` column names the genes, and choose another under
     **Session › Data mapping › Genes** (the file's default is
     `uns['cellucid']['defaults']['geneIdColumn']`, written by
     `prepare(var_gene_id_column=...)`);
   - if no column contains your name, the gene is absent from the file and no
     viewer setting brings it back — write the dataset again with the genes
     you need.

### Prevention

- Name genes in your dataset with the vocabulary your readers will type — the
  symbols on your figures, not an accession they would have to look up — and
  say which vocabulary it is for collaborators.

---

## Symptom: “Points disappeared / the dataset looks empty”

### Likely causes (ordered)

1) You applied filters earlier and forgot (hidden categories, continuous min/max range, outlier filtering).
2) You hid almost all categories (or clicked **Hide all**).
3) Continuous filter range is extremely narrow (or inverted).
4) Outlier filter is set below 100% and is hiding outliers.

### How to confirm

- Check **Filters**:
  - if it lists anything other than “No filters active”, filters are applied.
- For categorical fields:
  - look for unchecked categories.
- For continuous fields:
  - look at Min/Max sliders and the numeric readouts.
- For outliers:
  - check the **Outlier filter** percent in the Legend box (categorical fields only).

### Fix

1) For categorical filters:
   - click **Show all** in the categorical legend.
2) For continuous filters:
   - click **Reset** in the continuous legend's **Visible range** section.
3) For outlier filtering:
   - set the outlier slider back to `100%` (or switch to a field without outlier filtering).

### Prevention

- Watch **Filters** as your “single source of truth” for why points are missing.

---

## Symptom: “A category is greyed out / checkbox is disabled”

### Likely causes (ordered)

1) There are **0 available cells** in that category after other filters.

### How to confirm

- Hover the greyed-out row; it typically indicates “No cells available in this category after other filters”.
- The count may read `0 cells`.

### Fix

- Clear the filter(s) that removed those cells:
  - categorical: Show all on the other field(s)
  - continuous: Reset range filters

### Prevention

- When filtering with multiple fields, remember that a category can become empty even if it exists in the dataset.

---

## Symptom: “Colors look washed out / everything is the same color”

### Likely causes (ordered)

1) You filtered to a narrow range but kept the color domain wide (rescale is off).
2) The field is nearly constant (very low dynamic range).
3) The stored values themselves take few distinct levels (integer counts,
   binned scores), which shows as banding/steps.

### How to confirm

- If you have a narrow Min/Max slider range:
  - check whether **Rescale to range** is Off.
- Toggle rescale On and see if contrast improves.

### Fix

1) Turn **Rescale to range** On.
2) Choose a more appropriate colormap (Viridis/Cividis are good defaults).
3) If the stored values are coarse and you need more precision:
   - store the column as floating point upstream and write the file again.

### Prevention

- For publication-like figures, prefer rescale On when showing a restricted subset.
- For cross-view comparisons, prefer rescale Off for consistent mapping.

---

## Symptom: “Log scale makes everything gray”

### Likely causes (ordered)

1) The field has mostly zeros (common for sparse gene expression).
2) The field has negative values (some normalized/scaled outputs).
3) Only a handful of cells are positive, so almost everything else is grey. A
   field with *no* positive value cannot produce this: the toggle refuses to
   turn on (`Log color scale requires at least one positive field value.`).

### How to confirm

- Turn log scale off: if colors return, log incompatibility is the cause.

### Fix

- Turn **Log color scale** on only for fields with meaningful positive values.
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
- Check **Deleted fields**: the original is likely there and restorable.

### Fix

- If you didn’t mean to edit the source:
  - restore the original field from Deleted fields and continue from that.
- If you did mean it:
  - keep working on the derived field (subsequent edits may happen in place).

### Prevention

- Treat derived-field suffixes as a history trail (great for reproducibility).

---

## Symptom: “Restore did nothing and I got an error about the name”

### Likely causes (ordered)

1) A **visible field already holds that name**. Restore is refused rather than
   renamed, so the field stays in **Deleted fields**.

### How to confirm

- The error reads `Cannot restore "X" while that visible field name exists`.
- The row is still listed under **Deleted fields**.

### Fix

1) Rename the *visible* field that is occupying the name (or delete it).
2) Restore again. The success message is a plain `Restored "X"`.

### Prevention

- Avoid creating a derived field whose name matches something you have deleted;
  use meaningful suffixes.

:::{note}
Cellucid never appends a `(restored)` suffix. If you are looking for one, it
does not exist — the restore either succeeds under the original name or is
refused.
:::

---

## Symptom: “I can’t restore a deleted field”

### Likely causes (ordered)

1) You clicked **Confirm** in Deleted fields (purged restore capability).

### How to confirm

- The Deleted fields row may no longer appear, or restore may do nothing.
- Session restores will also not bring it back if it was confirmed.

### Fix

- There is no in-app undo for confirmed deletion.
- Reload the dataset if you need the original column again.

### Prevention

- Only use **Confirm deletion** when you are truly done with a field.

---

## Symptom: “Gene loading is slow / selecting a gene freezes the UI”

### Likely causes (ordered)

1) Large dataset + slow storage/network (remote data host).
2) Browser memory pressure (many genes loaded in one session).
3) The matrix is row-major (CSR), so the first gene waits while the browser
   builds a gene index, once per matrix.

### How to confirm

- DevTools → Network: gene value requests are slow.
- DevTools → Performance/Memory: memory grows as you load many genes.

:::{note}
A gene that is *refused* is not a performance problem: a column that is not
entirely finite is refused with a notification naming the gene and counting the
offending values. See {ref}`undrawable-continuous-payload`. Slowness is the
separate problem this section is about.
:::

### Fix

- Store a very large gene matrix column-major (CSC), so no gene index is built first.
- Avoid loading hundreds of genes in one session unless necessary.
- If the browser becomes unstable, reload and use a smaller set of genes.

### Prevention

- For interactive exploration, export the specific genes/scores you care about as obs fields when possible.

---

## Symptom: “Legend changes when I click different panels (multiview confusion)”

### Likely causes (ordered)

1) The active view changed; the sidebar always shows controls for the active view.

### How to confirm

- Click a different panel; observe the selected row in **Compare views** and the legend updating.

### Fix

- Use **Edit selected view** layout when you are configuring fields/legends.
- Then switch to **Grid compare** to compare.

### Prevention

- Check the label of the embedding picker (`View 2 shows`): it names the view the sidebar is editing.

---

## Still stuck?

If you can reproduce the issue, capture:
- a screenshot of **Coloring & filtering** (showing selectors + legend + active filters),
- the browser console output around the failure,
- and the dataset's `obs` and `var` columns (`adata.obs.dtypes`, `adata.var.head()`).

Then continue debugging from:
- {doc}`../q_troubleshooting_index/index` (broader issues)
