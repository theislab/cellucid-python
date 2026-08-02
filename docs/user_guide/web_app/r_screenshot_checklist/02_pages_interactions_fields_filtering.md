# Screenshot coverage — sections C, D, E

Scope: `user_guide/web_app/c_core_interactions/` (8 pages),
`user_guide/web_app/d_fields_coloring_legends/` (8 pages),
`user_guide/web_app/e_filtering/` (9 pages).

:::{note}
**This page was an assignment; the assignment was carried out.** It was written
as a pre-work audit and committed in `f304688` — the same commit that captured
the images it asked for — so from the moment it landed it described a state that
no longer existed. Every `has` record and every count in it was wrong on arrival.

It is now a **status record, measured against the tree on 2026-08-02**, not a
work order. The original per-page specification — dataset, interaction sequence,
crop target, cursor placement and filename for each requested capture — is
preserved in git: `git show f304688:docs/user_guide/web_app/r_screenshot_checklist/02_pages_interactions_fields_filtering.md`.
Read it there if you are producing one of the captures still outstanding below.
:::

## What changed

| Metric | At audit time (`f304688`) | Now |
| --- | --- | --- |
| `{figure}` blocks across C+D+E | 24 | **103** |
| Distinct images behind them | 10 | **46** |
| Pages with zero figures | 10 | **2** |
| Uses of `filtering/coloring-filtering-cell-type-panel.png` in C+D+E | 11 | **0** |
| Native capture scale | 1× | **2×**, emitted at the intrinsic width |

`coloring-filtering-cell-type-panel.png` — the single image that stood in for
eleven different topics — is no longer referenced anywhere in these three
sections. It survives on exactly one page,
`python_package/f_notebooks_tutorials/11_open_a_dataset_and_color_by_clusters.md:264`,
which is a legitimate single use.

## Current inventory

### `c_core_interactions/` — 26 blocks, 15 distinct images

| Page | Blocks | Images |
| --- | ---: | --- |
| `index.md` | 1 | `navigation-controls-orbit` |
| `01_navigation_modes_orbit_planar_free_fly.md` | 3 | `navigation-controls-orbit`, `-planar`, `-freefly` |
| `02_camera_controls_advanced.md` | 2 | `navigation-controls-freefly`, `camera-path-configured-panel` |
| `03_render_modes_points_vs_volumetric_smoke.md` | 1 | `render-mode-select` |
| `04_view_layout_live_snapshots_small_multiples.md` | 1 | `multiview-two-panels` |
| `05_dimension_switching_1d_2d_3d.md` | 3 | `window-dimension-1d`, `-2d`, `-3d` |
| `06_troubleshooting_core_interactions.md` | **0** | — |
| `07_screenshots.md` | 15 | the section's gallery |

### `d_fields_coloring_legends/` — 25 blocks, 11 distinct images

| Page | Blocks | Images |
| --- | ---: | --- |
| `index.md` | 3 | `field-selectors-three-routes`, `legend-categorical-clusters`, `legend-continuous-sscore` |
| `01_field_types_and_sources.md` | 2 | `legend-categorical-clusters`, `legend-continuous-sscore` |
| `02_field_selector_ux.md` | 3 | `field-selectors-three-routes`, `gene-search-dropdown-open`, `gene-search-no-match` |
| `03_color_by_behavior.md` | 1 | `legend-categorical-clusters` |
| `04_legend_behavior.md` | 4 | `legend-categorical-clusters`, `legend-colormap-menu-open`, `window-centroids-and-labels`, `filtering/legend-category-unavailable` |
| `05_troubleshooting_fields_legends.md` | **0** | — |
| `06_screenshots.md` | 10 | the section's gallery |
| `07_genes_in_the_built_in_samples.md` | 2 | `gene-search-dropdown-open`, `gene-search-no-match` |

### `e_filtering/` — 52 blocks, 20 distinct images

Every page in this section now carries at least one figure, and every image in
it is a purpose-built `filtering/` capture. `08_screenshots.md` is the gallery
(20 blocks); the teaching pages carry 1–7 each.

## Still outstanding

Confirmed against the tree, not inherited from the original audit.

1. **Two pages still have no figure at all.**
   `c_core_interactions/06_troubleshooting_core_interactions.md` and
   `d_fields_coloring_legends/05_troubleshooting_fields_legends.md`. Both are
   troubleshooting pages whose symptoms are visual, and `d_/05` still ends by
   telling the reader to capture a screenshot of a panel it never shows.
2. **No motion asset was ever produced.** The audit asked for recorded loops, or
   failing that before/after pairs, for orbit / planar / free-fly drag and for
   the multiview focus rule. Neither exists: there is no `.webm`, `.gif`, or
   `*-before`/`*-after` pair under `_static/screenshots/`. The three
   `navigation-controls-*` stills show the *control panels*, not the motion the
   page is about, so the "why did the wrong view move" gap is still open.
3. **Thirteen images in these sections have no provenance record**, so
   `node capture.mjs check` cannot re-verify them — 12 under `web_app/`
   (`app-overview-cell-type`, `camera-path-configured-panel`,
   `camera-path-transport-visible`, `camera-path-unconfigured`,
   `connectivity-edges-controls`, `connectivity-multiview`,
   `dark-theme-multiview`, `dimension-2d-planar-default`,
   `dimension-navigation-controls-2d`, `multiview-two-panels`,
   `startup-loaded-build`, `welcome-startup`) and one under `filtering/`
   (`coloring-filtering-cell-type-panel`). They predate the capture tool.
   `tests/test_docs_current_api_contract.py::test_screenshot_provenance_records_describe_the_published_images`
   holds the total across all topics at a number that may only decrease.
4. **KNN Connectivity has UI, working controls on all five samples and two
   screenshots — and still no prose page.** `index.html:857` labels the block
   `KNN Connectivity:`; `c_core_interactions/07_screenshots.md:95` is the only
   heading in the whole web-app guide that names the feature. The two teaching
   pages mention it only in passing (`03_render_modes…:30,50`).

### Resolved since the audit

- **Dimension switching and navigation mode** — `c_/05` now documents that a
  mode you chose yourself survives a dimension change, matching
  `viewer-contracts.js:27 getDefaultNavigationMode`.
- **`(restored)` suffix** — both `d_/02:221` and `d_/05:392` now state plainly
  that no such suffix exists, matching `overlay-public.js`.
- **`Use log scale`** — corrected to the app's actual label `Log color scale`
  (`assets/js/app/ui/modules/legend-renderer.js:200`) across `d_/03`, `d_/04`,
  `d_/05`, and — found in this pass, in a section the audit did not cover —
  `k_figure_export/04_quality_knobs_and_best_practices.md:161`. The string
  `Use log scale` no longer appears in the documentation or in the app.
- **`r_screenshot_checklist/` had no index** — it now has one, plus seven pages.
- **Community-annotation figures** — `_static/screenshots/community_annotation/`
  held one image at audit time and holds 16 now, so the hand-off the audit
  recorded is closed.

## Reference material (re-verified 2026-08-02)

### Multiview panels are not DOM elements

All panels render into the single `#glcanvas` through GL viewports
(`viewer.js:6349`, `gl.viewport(vx, vy, vw, vh)`). There is no per-panel DOM
node. To photograph one panel, crop the canvas region or set
`Layout: Edit selected view` so the panel fills the canvas. Never write a
capture scenario that crops to a per-panel selector.

### DOM crop targets

All forty confirmed present in `cellucid/index.html`:

`#visualization-section` · `#render-mode` · `#depth-controls` · `#renderer-controls` ·
`#smoke-controls` · `#points-controls` · `#connectivity-controls` ·
`#compare-views-section` · `#split-view-controls` · `#split-keep-view-btn` ·
`#camera-lock-btn` · `#split-view-badges-box` · `#split-view-badges-list` ·
`#view-layout-mode` · `#split-clear-btn` · `#dimension-controls` ·
`#dimension-select` · `#navigation-controls` · `#navigation-mode` ·
`#orbit-controls` · `#planar-controls` · `#freefly-controls` · `#reset-camera-btn` ·
`#coloring-filtering-section` · `#categorical-field-row` · `#continuous-field-row` ·
`#gene-expression-row` · `#gene-expression-dropdown` · `#display-options-container` ·
`#outlier-filter-container` · `#centroid-controls` · `#legend` ·
`#active-filters-container` · `#filter-count` · `#active-filters` ·
`#deleted-fields-section` · `#cinematic-camera-section` · `#shortcuts-section` ·
`#glcanvas` · `#sidebar`

### Dataset facts

Read from `cellucid-datasets/exports/*/dataset_identity.json`:

| Sample id | Cells | Genes | Dimensions | Default | Connectivity | Useful fields |
| --- | ---: | ---: | --- | ---: | --- | --- |
| `suo` | 561,947 | 5,103 | 1, 2, 3 | 3 | yes | `cell_type`, `organ`, `LVL0`–`LVL3`, `age`, `n_genes` |
| `garcia` | 219,731 | 5,754 | 1, 2, 3 | 3 | yes | same schema as suo |
| `he` | 71,650 | 5,152 | 1, 2, 3 | 3 | yes | same schema as suo |
| `kanemaru` | 131,636 | 3,691 | 1, 2, 3 | 3 | yes | same schema as suo |
| `pancreas` | 3,696 | 3,753 | 1, 2, 3 | 3 | yes | `clusters`, `clusters_coarse`, `cell_type`, `S_score`, `G2M_score` |

Because all five publish 1D, 2D and 3D, the "dataset provides only one embedding
dimension" state is **unreachable with any shipped sample**. Only the local
fixture `cellucid/tests/browser/fixtures/exports/current-ui-prepared` reproduces
it: 120 cells, **2D only**, one categorical field (`cell_type` =
alpha/beta/gamma, `centroid_outlier_quantile: 0.95`), one continuous field
(`score`), six genes (GAPDH, ACTB, CD3E, MS4A1, NKG7, LYZ). It cannot be used
for anything 3D, orbit, smoke, or realistic-looking.

### Label facts worth keeping

- The filter block is labelled `Active filters (selected view only):`
  (`index.html:1564`), not `Active filters`.
- `Zoom to cursor` appears twice with different wording — `Zoom to cursor
  (pinch-style)` in Compare Views (`index.html:1359`) and `Zoom to cursor` in the
  Camera Path block (`index.html:1716`). The duplicate navigation controls under
  Camera Path are still undocumented in `c_/02`.
- There is no `Original:` tooltip on the obs dropdown options; that tooltip
  exists only on gene-dropdown results
  (`field-selector-gene-expression.js:528`) and Deleted Fields rows
  (`field-selector-deleted-fields.js:432`).
