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

1) **Active field**: In **Coloring & Filtering**, is the correct selector active?
   - Categorical obs vs Continuous obs vs Gene Expression
2) **Legend toggles**:
   - Is **Log color scale** on (and your data has zeros/negatives)?
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
- For (4): reload the dataset; if it persists, inspect export integrity.

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
   underflow). Only the writers and the AnnData server report this kind; it is
   refused for the same reason.

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

If the values live on an AnnData server rather than in an export, the server
refuses the request before it sends a byte, and the message shows that instead:

> Gene "GATA3" could not be loaded: the server answered 422.  
> The server said: `{"error": "non_finite_continuous_payload", …}`  
> The server examined this payload and will not publish it.  
> Repair the values at the source and reload. Nothing in the viewer can work around a value that has no position on a colour scale.

The JSON body carries `kind`, `name`, `counts` (`total`, `offending`, `nan`,
`positive_infinity`, `negative_infinity`, `float32_overflow`,
`float32_underflow`), `examples` and a `message`, so DevTools → Network names
the gene and the reason without any further digging.

### Fix

Repair the values at the source. Reloading, re-picking the gene or clearing the
cache cannot change the outcome — nothing in the viewer can place an infinity on
a colour scale.

1) For an expression matrix:
   - sparse: `adata.X.data[~np.isfinite(adata.X.data)] = 0`
   - dense: `np.nan_to_num(adata.X, copy=False)`
2) For a continuous obs column:
   - `adata.obs["<key>"] = adata.obs["<key>"].replace([np.inf, -np.inf], np.nan)`
3) Then re-run `prepare()` (or restart `serve_anndata` / `show_anndata`) and
   reload the dataset.

If you meant the infinities as “missing”, replace them with `NaN` rather than
with `0`: Cellucid draws `NaN` as the neutral grey meaning “not measured here”,
whereas `0` is a real measurement and will be coloured as one.

### Prevention

- Check a column before you publish it: `np.isfinite(values).all()`.
- Prepared exports cannot carry this problem — both the Python and the R writer
  refuse a non-finite continuous obs field or gene at preparation time, with the
  same counts and the same first-affected indices. If you hit this at all, you
  are reading a raw `.h5ad`/Zarr in the browser or serving one live.

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

1) The gene is published under a different name than the one you typed.
2) The gene is not in this export at all — preparation published only some of
   the source’s genes.
3) The search is substring-based; a short query matches many genes and Enter picks the top one.
4) The gene list is present but very large; your query is too broad.

### How to confirm

- Try a longer, more specific substring. Matching ignores letter case, so
  `gata3` finds `GATA3`.
- If results show `...and N more`, you’re matching too many.
- Open the field selector and read a few gene names directly. Every gene has
  exactly one name in an export, and search matches that one name — there is no
  second identifier to search instead, and no alias mapping.
- Compare the count in the empty-state message (“This dataset publishes *N* gene
  names…”) against the gene count you expect from the source. A large gap means
  preparation published a subset.

### If you are using a built-in sample

The four human samples — Suo, Garcia, He, and Kanemaru — publish only the genes
whose Ensembl accession resolved to a symbol, which is between 45.1% and 70.2%
of the gene axis of the file each was prepared from. A gene that is genuinely in
the study can therefore be absent, under any spelling. The rule, the per-sample
counts, and what to do instead are in
{doc}`07_genes_in_the_built_in_samples`.

### Fix

1) Type a more specific query and click the correct gene explicitly (don’t rely on Enter).
2) If the name you expected is not there, it is either published under another
   name or not published at all:
   - check which `var` column the export named its genes from
     (`prepare(var_gene_id_column=...)`, `show_anndata(gene_id_column=...)`);
   - if that column does not contain your name either, the gene is absent from
     the export and no viewer setting brings it back — re-export the dataset
     with the genes and names you need, or load the source through server mode
     or Jupyter.

### Prevention

- Name genes in your export with the vocabulary your readers will type — the
  symbols on your figures, not an accession they would have to look up — and
  say which vocabulary it is for collaborators.

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
- Check **Deleted Fields**: the original is likely there and restorable.

### Fix

- If you didn’t mean to edit the source:
  - restore the original field from Deleted Fields and continue from that.
- If you did mean it:
  - keep working on the derived field (subsequent edits may happen in place).

### Prevention

- Treat derived-field suffixes as a history trail (great for reproducibility).

---

## Symptom: “Restore did nothing and I got an error about the name”

### Likely causes (ordered)

1) A **visible field already holds that name**. Restore is refused rather than
   renamed, so the field stays in **Deleted Fields**.

### How to confirm

- The error reads `Cannot restore "X" while that visible field name exists`.
- The row is still listed under **Deleted Fields**.

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

- DevTools → Network: gene value requests are slow.
- DevTools → Performance/Memory: memory grows as you load many genes.

:::{note}
A gene request that *fails* against the Python AnnData server is not a
performance problem and no longer looks like one. The server validates the
column before it sends a byte, and a column that is not entirely finite is
answered **HTTP 422** with a JSON body naming the gene and counting the offending
values — so the Network tab tells you the cause directly. See
{ref}`undrawable-continuous-payload`. Slowness is the separate problem this
section is about.
:::

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
- {doc}`../q_troubleshooting_index/index` (broader issues)
