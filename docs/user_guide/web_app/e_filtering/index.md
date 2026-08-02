# Filtering (Visibility, Outliers, and Filter Stacks)

These pages explain how Cellucid’s filtering system controls **what is visible**, without deleting data.

They are written for mixed audiences:
- **Wet lab / non-technical**: click-by-click recipes and “what success looks like”.
- **Computational users**: exact semantics, performance scaling, and data edge cases (NaN/Inf, missing fields).
- **Power users**: how stacked filters interact with highlights, analysis, and session persistence.

:::{important}
Filtering is about **visibility**:

- A filtered-out cell still exists in the dataset.
- Other features may still reference it (highlights, analysis caches, saved sessions), depending on how the app scopes state.
:::

---

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```

---

## Fast path (choose what you’re trying to do)

| You want to… | Best filter | Where it lives | Notes |
|---|---|---|---|
| Hide specific clusters/samples/batches | Category visibility filter | Categorical legend checkboxes | Fast and intuitive; use `Show All` to undo |
| Gate by QC / score ranges | Numeric range filter | Continuous legend Min/Max sliders | Turn Live filtering off for large datasets |
| Hide “fringe” cells within clusters | Outlier filter (latent space) | Display options (when available) | `95%` is a common starting point |
| Show only gene-expressing cells | Gene expression range filter | Gene legend Min/Max sliders | Applies only while that gene is active |

---

## Where all of this lives

Everything on this page happens inside one sidebar accordion,
**Coloring & Filtering**. Open it once and keep it open.

```{figure} ../../../_static/screenshots/filtering/panel-clusters-unfiltered.png
:alt: The Coloring and Filtering panel with a Categorical obs select holding clusters, empty Continuous obs and Gene Expression controls, a dashed Display options box containing an Outlier filter (latent space) slider at 100 percent, Show centroids and Show labels checkboxes, Show All and Hide All buttons and eight category rows with colour swatches and cell counts, and below the box a heading reading ACTIVE FILTERS (SELECTED VIEW ONLY), the line Showing all 3,696 points, and a bordered area reading No filters active.
:width: 516px

The panel in its resting state on the Pancreas sample. The three selectors at
the top decide **what colours the points**. The dashed **Display options** box
and everything below it decide **which points are drawn**. The last two lines —
the count, and the filter list — are the only place the app tells you the truth
about what is hidden.
```

Once you actually filter something, that last block fills up:

```{figure} ../../../_static/screenshots/filtering/active-filters-three-rows.png
:alt: The Active filters block reading Showing 709 of 3,696 points above three rows, each with a ticked checkbox and a circled cross; the first reads S_score minus 0.15 to 1.14 followed by a cell count, the second reads clusters colon hiding Ductal, Epsilon and is cut off, the third reads clusters colon outlier less than or equal to 95 percent.
:width: 480px

Three of the four filter types active at once on the same 3,696 cells. Filters
combine with **AND**: a cell is drawn only if it survives every ticked row.
```

---

## Recommended reading order

1) {doc}`01_filtering_mental_model` (how filtering works conceptually)
2) {doc}`03_filter_stack_ui_active_filters` (how to see what’s active and undo safely)
3) {doc}`02_outlier_filtering_per_active_field` (the “fastest” common filter)
4) {doc}`04_common_filter_types_document_every_filter_the_ui_exposes` (reference: every filter you can create)
5) {doc}`05_performance_considerations` + {doc}`06_edge_cases_filtering` (avoid surprises)
6) {doc}`07_troubleshooting_filtering` (symptom → diagnosis → fix)
7) {doc}`08_screenshots` (every state above, photographed)

---

## Pages in this section

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`book;1.5em;sd-mr-1` Filtering mental model
:link: 01_filtering_mental_model
:link-type: doc

What filtering means (visibility, not deletion), how stacking works, and how filtering interacts with other features.
:::

:::{grid-item-card} {octicon}`sliders;1.5em;sd-mr-1` Outlier filtering
:link: 02_outlier_filtering_per_active_field
:link-type: doc

How to hide extreme values for the active field, what the slider means, and what happens when quantiles are missing.
:::

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` Filter stack UI
:link: 03_filter_stack_ui_active_filters
:link-type: doc

How to inspect, disable, and remove filters using the Active filters panel.
:::

:::{grid-item-card} {octicon}`list-ordered;1.5em;sd-mr-1` Filter types reference
:link: 04_common_filter_types_document_every_filter_the_ui_exposes
:link-type: doc

Catalog of every filter type exposed in the UI, with semantics, defaults, and edge cases.
:::

:::{grid-item-card} {octicon}`pulse;1.5em;sd-mr-1` Performance
:link: 05_performance_considerations
:link-type: doc

What gets slow, why it gets slow, and safe workflows for large datasets.
:::

:::{grid-item-card} {octicon}`info;1.5em;sd-mr-1` Edge cases
:link: 06_edge_cases_filtering
:link-type: doc

“Surprising” states like 0 visible cells, NaN/Inf handling, and missing fields referenced by filters.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting
:link: 07_troubleshooting_filtering
:link-type: doc

Symptom-based debugging for empty plots, filters that “do nothing”, and slow filtering.
:::

:::{grid-item-card} {octicon}`image;1.5em;sd-mr-1` Verified captures
:link: 08_screenshots
:link-type: doc

Every filtering state this section describes, photographed from the running app
on one reproducible sample.
:::

::::
