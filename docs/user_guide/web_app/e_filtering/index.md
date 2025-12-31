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

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: filtering-panel-overview
Suggested filename: filtering/00_filtering-panel-overview.png
Where it appears: User Guide → Web App → Filtering → index.md
Capture:
  - UI location: left sidebar with “Coloring & Filtering”, “Display options”, and “Active filters” visible
  - State prerequisites: dataset loaded; at least one filter active (so Active filters is not empty)
  - Action to reach state: hide a category OR apply a continuous range filter, then open Active filters
Crop:
  - Include: field selectors + legend area + outlier slider (if visible) + Active filters
  - Include: enough canvas to see the filtered point cloud
Redact:
  - Remove: private dataset identifiers, local file paths, private org/repo names
Alt text:
  - Filtering controls in the left sidebar with Active filters showing at least one active filter.
Caption:
  - Point readers to Active filters as the source of truth for visibility.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the filtering panel overview.
:width: 100%

Filtering controls live in the left sidebar; Active filters is the source of truth for what currently affects visibility in the active view.
```

---

## Recommended reading order

1) `01_filtering_mental_model` (how filtering works conceptually)
2) `03_filter_stack_ui_active_filters` (how to see what’s active and undo safely)
3) `02_outlier_filtering_per_active_field` (the “fastest” common filter)
4) `04_common_filter_types_document_every_filter_the_ui_exposes` (reference: every filter you can create)
5) `05_performance_considerations` + `06_edge_cases_filtering` (avoid surprises)
6) `07_troubleshooting_filtering` (symptom → diagnosis → fix)
7) `08_screenshots` (capture checklist)

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

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Screenshots
:link: 08_screenshots
:link-type: doc

Screenshot capture checklist for the filtering section.
:::

::::
