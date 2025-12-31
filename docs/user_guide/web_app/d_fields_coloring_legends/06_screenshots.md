# Screenshots

This page is a **screenshot capture checklist** for the Fields/Coloring/Legends section.

It exists so you (or a collaborator) can capture screenshots once, systematically, without hunting through every page.

---

## How placeholders work

All pages in this section currently use:

- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

Each placeholder is preceded by an HTML comment with:
- what state to capture,
- what to crop/redact,
- suggested filename conventions,
- suggested alt text + caption content,
- and callouts to annotate in the image editor.

General guidance lives in:
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

---

## Recommended screenshot set (fields, coloring, legends)

### Overview (1 screenshot)

1) **Coloring & Filtering panel overview** (selectors + Display options + Active filters)

Page:
- `d_fields_coloring_legends/index`

Suggested filename:
- `web_app/fields_legends/01_coloring-filtering-panel-overview.png`

### Field selection UX (2 screenshots)

2) **Three selector rows + action buttons** (Categorical obs / Continuous obs / Gene Expression)

Page:
- `d_fields_coloring_legends/02_field_selector_ux`

Suggested filename:
- `web_app/fields_legends/10_field-selector-rows-actions.png`

3) **Gene expression dropdown open** (search results visible)

Page:
- `d_fields_coloring_legends/02_field_selector_ux`

Suggested filename:
- `web_app/fields_legends/11_gene-search-dropdown.png`

### Legend types (2 screenshots)

4) **Categorical legend example** (checkboxes + swatches + counts)

Page:
- `d_fields_coloring_legends/01_field_types_and_sources`

Suggested filename:
- `web_app/fields_legends/02_categorical-legend-example.png`

5) **Continuous legend example** (colorbar + toggles + sliders)

Page:
- `d_fields_coloring_legends/01_field_types_and_sources`

Suggested filename:
- `web_app/fields_legends/03_continuous-legend-example.png`

### Legend interactions (3 screenshots)

6) **Categorical legend anatomy** (show/hide all + rename/delete icons)

Page:
- `d_fields_coloring_legends/04_legend_behavior`

Suggested filename:
- `web_app/fields_legends/30_categorical-legend-anatomy.png`

7) **Continuous colormap menu open**

Page:
- `d_fields_coloring_legends/04_legend_behavior`

Suggested filename:
- `web_app/fields_legends/32_continuous-colormap-menu.png`

8) **Live filtering OFF + FILTER enabled** (continuous filtering UI)

Page:
- `d_fields_coloring_legends/03_color_by_behavior`

Suggested filename:
- `web_app/fields_legends/21_continuous-filter-live-toggle.png`

### Category editing + restoration (2–3 screenshots)

9) **Merge categories confirm dialog**

Page:
- `d_fields_coloring_legends/04_legend_behavior`

Suggested filename:
- `web_app/fields_legends/31_categorical-merge-confirm.png`

10) **Deleted Fields panel** (Restore vs Confirm)

Page:
- `d_fields_coloring_legends/02_field_selector_ux`

Suggested filename:
- `web_app/fields_legends/12_deleted-fields-panel.png`

### Outlier filtering (1 screenshot, if your dataset supports it)

11) **Outlier filter slider visible and set below 100%**

Page:
- `d_fields_coloring_legends/03_color_by_behavior`

Suggested filename:
- `web_app/fields_legends/22_outlier-filter-slider.png`

---

## Highly recommended “failure mode” screenshots (support gold)

These save huge amounts of time in onboarding and debugging:

- “Everything is gray” state (no field selected)
- “Log scale makes everything gray” state on a sparse gene (log on, mostly zeros)
- “Category row disabled (0 available cells)” state (after applying another filter)

Suggested page to add these to:
- `d_fields_coloring_legends/05_troubleshooting_fields_legends`
