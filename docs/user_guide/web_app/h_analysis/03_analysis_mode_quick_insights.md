# Analysis mode: Quick (automatic composition + stats)

**Audience:** everyone (wet lab + computational)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- What Quick mode summarizes (and what it deliberately does not)
- How **Dynamic** vs **Manual** page selection works
- How to interpret the two Quick sections: **Composition** and **Statistics**
- Common failure modes (“nothing shows up”, “wrong page”, “missing fields”)

**Prerequisites:**
- A dataset loaded
- At least one highlight page (Quick can follow the active page automatically)

---

## What Quick mode is for (and what it is not)

Quick mode is the “sanity check / first glance” analysis:

- **Composition**: “What labels make up these cells?” (categorical obs)
- **Statistics**: “What are typical QC/score values here?” (continuous obs)

Quick mode is intentionally:
- **fast** (no configuration and no “Run” button — it recomputes as you switch pages),
- **cheap** (it reads only obs columns already in memory, so it never waits on gene expression),
- **obs-focused** (cell metadata), not gene-level hypothesis testing.

Quick mode is **not** for:
- gene–gene or field–gene relationships → use {doc}`05_analysis_mode_correlation_analysis`
- differential expression → use {doc}`06_analysis_mode_differential_expression_de`
- marker discovery across many groups → use {doc}`08_analysis_mode_genes_panel`

---

## Inputs and scope (what cells does it summarize?)

Quick mode summarizes **highlight pages**.

### Dynamic mode (default): “follow the active page”

By default, Quick uses **Dynamic** mode:
- it always summarizes the currently active highlight page,
- switching pages in Highlighted Cells updates Quick automatically.

### Manual mode: “combine pages intentionally”

Quick also supports manual selection:
- you choose one or more pages,
- Quick summarizes the **union** of those pages.

:::{note}
Quick mode summarizes across the selected pages *as one combined set*.

If you want side-by-side per-page comparisons with statistical tests, use {doc}`04_analysis_mode_detailed_analysis`.
:::

### Filtering note (visibility vs membership)

Quick uses page membership (stored cell indices), not canvas visibility.
If you expected filters to affect analysis, see {doc}`01_analysis_mental_model`.

---

## Outputs (what you should see)

Quick mode renders two sections:

### 1) Composition (categorical obs)

For each selected composition field:
- Quick reports the **top 5 categories by count** across the selected pages.
- It renders a stacked bar where segment widths are **percent of cells** in those top categories.
- It prints the top 3 categories as short labels for quick scanning.

Interpretation notes:
- The bar shows only the **top categories**; the long tail is not shown.
- Missing values may appear as `(missing)` when present.

### 2) Statistics (continuous obs)

For each selected statistics field, Quick reports:
- **Mean** — the average value
- **Median** — the middle value, half above and half below
- **Std** — how spread out the values are around the mean

These numbers are exact at any page size. Nothing here is sampled or estimated.

:::{note}
**How the three numbers are computed, and one thing to watch.**

Quick makes a single streaming pass over every cell in every selected page,
accumulating the mean and the sum of squared deviations, then allocates a
`Float32Array` holding *all* the finite values, sorts it, and reads the median
out of the middle (averaging the two middle values for an even count). Cost is
therefore O(n log n) in the number of selected cells, dominated by the sort — and
the result is the exact median, not an approximation, no matter how large the
selection is.

**Std is the population standard deviation** — it divides by *n*, not by *n − 1*.
That matters if you are checking the number against pandas, whose
`Series.std()` defaults to the *sample* standard deviation and divides by
*n − 1*. On a page of a few thousand cells the difference is invisible; on a page
of 20 cells Quick's value reads about 2.5% lower. `numpy.std(values)` reproduces
it exactly, because numpy defaults the other way. Detailed mode's **Std Dev**
column uses the same denominator as Quick, so the two panels agree with each
other.
:::

---

## Fast path (wet lab / non-technical)

Goal: “Tell me what this group of cells looks like in 2 minutes.”

1) Create or select a highlight page
   - Go to Highlighted Cells.
   - If needed, add a page and put a selection into it.
2) Open **Analysis → Quick**
   - You should immediately see a header like “\<PageName\>: X cells”.
3) In **Composition**, choose fields you care about
   - Common picks: `cell_type`, `cluster`, `sample`, `batch`.
4) In **Statistics**, choose QC fields
   - Common picks: `n_counts`, `n_genes`, `pct_mito`, or any score field.
5) Interpret
   - If composition looks wrong: you may be on the wrong page.
   - If statistics look extreme: you may have selected low-quality cells or a biased subset.

What success looks like:
- the page name matches what you intended,
- the cell count is plausible,
- composition and QC metrics match your expectations.

---

## Practical path (computational users)

### Quick is aggregated, not per-page

Quick does *not* show a per-page breakdown.
It computes summaries over the combined selected cells.

If you need:
- per-page distributions,
- plot customization,
- statistical tests,

use {doc}`04_analysis_mode_detailed_analysis`.

### “Effective n” can differ by field

- Composition counts exclude missing values.
- Statistics treat `NaN` as missing and drop it from the count.

So a page with 50,000 cells may effectively contribute fewer values for a field
with missingness. An `Inf` is a different matter: it is not a missing value and
Quick refuses to summarise a column containing one, because a mean or a spread
that silently swallowed an infinity would be meaningless. Fix the column
upstream rather than looking for a setting.

---

## Edge cases and pitfalls

- **No pages exist** → Quick has nothing to summarize. Create a highlight page first.
- **Empty page** (0 cells) → header shows 0 cells; sections are empty.
- **No categorical obs fields** → Composition has nothing to choose.
- **No continuous obs fields** → Statistics has nothing to choose.
- **Huge number of categories** → only top categories are shown; use Detailed + Bar Plot for full distribution.
- **Very small pages** (< 20–50 cells) → mean/std are unstable; interpret cautiously.

---

## Troubleshooting (Quick mode)

### Symptom: “Quick is blank / says there are no pages”

Likely causes:
- no highlight pages exist yet.

How to confirm:
- open Highlighted Cells and check if page tabs exist.

Fix:
- create a page and add cells to it (see {doc}`../f_highlighting_selection/index`).

### Symptom: “Quick shows the wrong group”

Likely causes:
- Quick is in **Dynamic** mode and a different page is active than you think.

How to confirm:
- open Quick’s **Page Selection** collapsible at the bottom and check the mode indicator.

Fix:
- switch to **Manual** mode and explicitly pick the page(s) to summarize.

### Symptom: “Composition/Statistics says ‘No … fields selected’”

Fix:
- use the “Choose … fields” dropdown and select fields.

---

## Interface reference

```{figure} ../../../_static/screenshots/analysis/quick-insights.png
:alt: The Quick panel headed "Page 1: 262 cells" with a Composition section showing one stacked bar per field with the top categories named beneath it, and a Statistics section with a Field / Mean / Median / Std table listing S_score and G2M_score.
:width: 488px

Quick opens with no configuration at all. **Composition** draws one stacked bar
per selected categorical field, segment width proportional to the share of cells
in that category, with the top three named underneath. **Statistics** gives mean,
median and standard deviation for each selected continuous field. The `⋯` button
on each section heading chooses which fields are shown; by default it is the
first two of each kind.
```

Quick is the only mode with **no expanded view** — there is no `⤢ Expand` button,
because there is no plot to enlarge and no table that would not already fit. If
you need a distribution rather than a summary, that is
{doc}`04_analysis_mode_detailed_analysis`.

Before any page has cells, the panel says so rather than showing an empty table:

```{figure} ../../../_static/screenshots/analysis/analysis-empty-no-pages.png
:alt: The Quick panel reading "Page 1: 0 cells" with a Composition section stating "No composition data available." and a Statistics section stating "No statistics data available."
:width: 488px

The no-pages state. `Page 1: 0 cells` is the diagnosis: the page exists but
nothing has been confirmed into it yet. Fix it in **Highlighting**, not here —
see {doc}`11_screenshots` step 1.
```

---

## Next steps

- {doc}`04_analysis_mode_detailed_analysis` (side-by-side comparisons, plots, statistical tests)
- {doc}`05_analysis_mode_correlation_analysis` (relationships between variables/genes)
- {doc}`10_troubleshooting_analysis` (analysis-wide failures: missing expression, slow DE, etc.)
