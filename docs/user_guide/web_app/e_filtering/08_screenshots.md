# Screenshots

This page is a **screenshot capture checklist** for the Filtering section.

It exists so you (or a collaborator) can capture screenshots once, systematically, without hunting through every page.

---

## How placeholders work

All pages in this section can use:

- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

Recommended screenshot storage:

- `cellucid-python/docs/_static/screenshots/filtering/`

The “suggested filenames” below assume you place images under that folder.

Each screenshot placeholder block should include an HTML comment describing:
- what state to capture,
- what to crop/redact,
- suggested filename conventions,
- and suggested alt text + caption content.

General guidance lives in:
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

---

## Recommended screenshot set (filtering)

### Overview (1 screenshot)

1) **Filtering controls overview** (Coloring & Filtering + Active filters visible)

Page:
- `e_filtering/index`

Suggested filename:
- `filtering/00_filtering-panel-overview.png`

### Outlier filtering (2 screenshots)

2) **Outlier filter slider visible** (default state)

Page:
- `e_filtering/02_outlier_filtering_per_active_field`

Suggested filename:
- `filtering/01_outlier-filter-slider-default.png`

3) **Outlier filter effect** (slider moved; outliers gone)

Page:
- `e_filtering/02_outlier_filtering_per_active_field`

Suggested filename:
- `filtering/02_outlier-filter-slider-after.png`

### Filter stack (2–3 screenshots)

4) **Active filters empty state**

Page:
- `e_filtering/03_filter_stack_ui_active_filters`

Suggested filename:
- `filtering/10_active-filters-empty.png`

5) **Active filters with multiple stacked filters** (2–3 filters enabled)

Page:
- `e_filtering/03_filter_stack_ui_active_filters`

Suggested filename:
- `filtering/11_active-filters-multiple.png`

6) **Active filters with one disabled** (to show disable vs remove)

Page:
- `e_filtering/03_filter_stack_ui_active_filters`

Suggested filename:
- `filtering/12_active-filters-one-disabled.png`

### Failure modes (high value support screenshots)

7) **All cells filtered out** (plot empty, UI message visible)

Page:
- `e_filtering/06_edge_cases_filtering` and/or `e_filtering/07_troubleshooting_filtering`

Suggested filename:
- `filtering/20_all-cells-filtered-out.png`

8) **Outlier slider unavailable** (disabled/hidden with explanation, if applicable)

Page:
- `e_filtering/02_outlier_filtering_per_active_field`

Suggested filename:
- `filtering/21_outlier-slider-unavailable.png`

---

## Nice-to-have screenshots (high instructional value)

These aren’t strictly required, but they make the section dramatically easier for beginners:

9) **Continuous filtering: Live filtering OFF with FILTER enabled**

Page:
- `e_filtering/04_common_filter_types_document_every_filter_the_ui_exposes`

Suggested filename:
- `filtering/30_continuous-live-filtering-off.png`

10) **Categorical legend: categories disabled due to upstream filters**

Page:
- `e_filtering/06_edge_cases_filtering`

Suggested filename:
- `filtering/31_categories-disabled-no-available.png`

11) **Category filter disabled but categories de-emphasized (gray)**

Page:
- `e_filtering/06_edge_cases_filtering`

Suggested filename:
- `filtering/32_category-filter-disabled-gray.png`
