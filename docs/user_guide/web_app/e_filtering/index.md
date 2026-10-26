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
| Hide specific clusters/samples/batches | Category visibility filter | Categorical legend checkboxes | Fast and intuitive; use `Show all` to undo |
| Gate by QC / score ranges | Numeric range filter | Continuous legend Min/Max sliders | Turn Live filtering off for large datasets |
| Hide “fringe” cells within clusters | Outlier filter | Legend box, for a categorical field, once a latent space is chosen | `95%` is a common starting point |
| Show only gene-expressing cells | Gene expression range filter | Gene legend Min/Max sliders | Applies only while that gene is active |

---

## Where all of this lives

Everything on this page happens inside one sidebar accordion,
**Coloring & filtering**. Open it once and keep it open.

```{figure} ../../../_static/screenshots/filtering/panel-clusters-unfiltered.png
:alt: The Coloring & filtering panel with the clusters field selected: Color by with its three field selectors, Show cell info, the Legend group with Show centroids, Show labels and the Outlier filter at 100 percent, the Show all and Hide all buttons above the eight-category legend, and Filters reading Showing all 3,696 cells above No filters active.
:width: 532px

The panel in its resting state on the Pancreas sample. The three selectors at
the top decide **what colours the points**. The **Legend** group (with the
**Outlier filter**) and everything below it decide **which points are
drawn**. The last two lines —
the count, and the filter list — are the only place the app tells you the truth
about what is hidden.
```

Once you actually filter something, that last block fills up:

```{figure} ../../../_static/screenshots/filtering/active-filters-three-rows.png
:alt: The Filters block holding three rows of different kinds at once: an S_score range row ending in its cell count, cut off with an ellipsis, a clusters hiding row, and a clusters: outlier ≤ 95% row with no count.
:width: 496px

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

How to hide cells far from their category's centre in a latent space, choosing that latent space, and where the distances come from (stored in the file or computed once).
:::

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` Filter stack UI
:link: 03_filter_stack_ui_active_filters
:link-type: doc

How to inspect, disable, and remove filters using the Filters panel.
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
