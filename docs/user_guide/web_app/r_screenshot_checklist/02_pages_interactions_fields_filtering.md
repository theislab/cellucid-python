# Screenshot coverage — sections C, D, E

:::{important}
Audited before Cellucid 0.9.5, which reads AnnData only. The capture
instructions use the current controls, fixtures and sample catalog; where a
record says what a page or an image showed at audit time (the **Prepared**
button, a prepared export), that is history. {doc}`index` maps the old names.
:::

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

### `c_core_interactions/` — 28 blocks, 16 distinct images

| Page | Blocks | Images |
| --- | ---: | --- |
| section index | 1 | `navigation-controls-orbit` |
| {doc}`/user_guide/web_app/c_core_interactions/01_navigation_modes_orbit_planar_free_fly` | 3 | `navigation-controls-orbit`, `-planar`, `-freefly` |
| {doc}`/user_guide/web_app/c_core_interactions/02_camera_controls_advanced` | 2 | `navigation-controls-freefly`, `camera-path-configured-panel` |
| {doc}`/user_guide/web_app/c_core_interactions/03_render_modes_points_vs_volumetric_smoke` | 2 | `render-mode-select`, `window-smoke-cloud` |
| {doc}`/user_guide/web_app/c_core_interactions/04_view_layout_live_snapshots_small_multiples` | 1 | `multiview-two-panels` |
| {doc}`/user_guide/web_app/c_core_interactions/05_dimension_switching_1d_2d_3d` | 3 | `window-dimension-1d`, `-2d`, `-3d` |
| {doc}`/user_guide/web_app/c_core_interactions/06_troubleshooting_core_interactions` | **0** | — |
| `07_screenshots.md` | 16 | the section's gallery |

### `d_fields_coloring_legends/` — 25 blocks, 11 distinct images

| Page | Blocks | Images |
| --- | ---: | --- |
| section index | 3 | `field-selectors-three-routes`, `legend-categorical-clusters`, `legend-continuous-sscore` |
| {doc}`/user_guide/web_app/d_fields_coloring_legends/01_field_types_and_sources` | 2 | `legend-categorical-clusters`, `legend-continuous-sscore` |
| {doc}`/user_guide/web_app/d_fields_coloring_legends/02_field_selector_ux` | 3 | `field-selectors-three-routes`, `gene-search-dropdown-open`, `gene-search-no-match` |
| {doc}`/user_guide/web_app/d_fields_coloring_legends/03_color_by_behavior` | 1 | `legend-categorical-clusters` |
| {doc}`/user_guide/web_app/d_fields_coloring_legends/04_legend_behavior` | 4 | `legend-categorical-clusters`, `legend-colormap-menu-open`, `window-centroids-and-labels`, `filtering/legend-category-unavailable` |
| {doc}`/user_guide/web_app/d_fields_coloring_legends/05_troubleshooting_fields_legends` | **0** | — |
| {doc}`/user_guide/web_app/d_fields_coloring_legends/06_screenshots` | 10 | the section's gallery |
| {doc}`/user_guide/web_app/d_fields_coloring_legends/07_genes_in_the_built_in_samples` | 2 | `gene-search-dropdown-open`, `gene-search-no-match` |

### `e_filtering/` — 52 blocks, 20 distinct images

Every page in this section now carries at least one figure, and every image in
it is a purpose-built `filtering/` capture. `08_screenshots.md` is the gallery
(20 blocks); the teaching pages carry 1–7 each.

## Still outstanding

Confirmed against the tree, not inherited from the original audit.

1. **Two pages still have no figure at all.**
   {doc}`/user_guide/web_app/c_core_interactions/06_troubleshooting_core_interactions` and
   {doc}`/user_guide/web_app/d_fields_coloring_legends/05_troubleshooting_fields_legends`. Both are
   troubleshooting pages whose symptoms are visual, and `d_/05` still ends by
   telling the reader to capture a screenshot of a panel it never shows.
2. **No motion asset was ever produced.** The audit asked for recorded loops, or
   failing that before/after pairs, for orbit / planar / free-fly drag and for
   the multiview focus rule. Neither exists: there is no `.webm`, `.gif`, or
   `*-before`/`*-after` pair under `_static/screenshots/`. The three
   `navigation-controls-*` stills show the *control panels*, not the motion the
   page is about, so the "why did the wrong view move" gap is still open.
3. **Ten images in these sections have no provenance record**, so
   `node capture.mjs check` cannot re-verify them — 9 under `web_app/`
   (`app-overview-cell-type`, `camera-path-configured-panel`,
   `camera-path-transport-visible`, `camera-path-unconfigured`,
   `connectivity-edges-controls`, `dimension-2d-planar-default`,
   `dimension-navigation-controls-2d`, `startup-loaded-build`,
   `welcome-startup`) and one under `filtering/`
   (`coloring-filtering-cell-type-panel`). They predate the capture tool.
   The three multiview frames that used to be on this list —
   `multiview-two-panels`, `connectivity-multiview` and
   `dark-theme-multiview` — went stale when Grid compare moved its panels out
   from under the sidebar, and were recaptured on `pancreas` by
   `buildMultiviewScenarios()` (`core-multiview-two-panels`,
   `core-connectivity-multiview`, `performance-dark-theme-multiview`), so they
   now carry provenance.
   `tests/test_docs_current_api_contract.py::test_screenshot_provenance_records_describe_the_published_images`
   holds the total across all topics at a number that may only decrease.
4. **KNN Connectivity has UI, working controls on all five samples and two
   screenshots — and still no prose page.** `index.html:857` labels the block
   `KNN Connectivity:`; `c_core_interactions/07_screenshots.md:95` is the only
   heading in the whole web-app guide that names the feature. The two teaching
   pages mention it only in passing (`03_render_modes…:30,50`).
5. **The outlier filter's captures show the earlier control, and its new
   states have none.** `filtering/outlier-slider-95`,
   `filtering/panel-clusters-unfiltered`, `filtering/window-outlier-100-off`
   and `filtering/window-outlier-95-applied` show a control labelled
   `Outlier filter (latent space)` above the centroid checkboxes in a box titled
   `Display options`; it is now `Outlier filter`, below **Show centroids**, in
   the box titled **Legend** (`index.html` `#display-options-container`,
   `#outlier-filter-container`). Recapture those four, then add one capture per
   state of `ui/modules/outlier-filter-control.js`, each cropped to
   `#outlier-filter-container` +8px at 1440×1000, 2×, for
   {doc}`/user_guide/web_app/e_filtering/02_outlier_filtering_per_active_field`:

   | Filename (`filtering/`) | State | How to reach it |
   | --- | --- | --- |
   | `outlier-no-latent` | no latent space chosen: only the label, its (i) and one **Choose latent space** button (no slider, no paragraph; the candidate arrays are the button's title) | scenario `filtering-outlier-no-latent`: `fixtures/pca-only.h5ad` opened with **Open… › File…**, colour by `cell_type`; or `pancreas`, colour by `clusters`, **Session › Data mapping › Latent space** set to **None**, Apply |
   | `outlier-cost-line` | quantiles to compute, nothing read yet: the slider and one short line, *“First use reads X_pca (… MB).”* | `suo`, colour by `cell_type`, Latent space set to `X_pca` (the file stores them for `scanvi` only), Apply; do not move the slider |
   | `outlier-progress-cancel` | computing: the progress bar, the line with bytes read, rate and time left, and **Cancel** | as `outlier-cost-line`, then move the slider to `95%` and capture while the bar runs (also capture the progress notification, `outlier-progress-notification`, cropped to the notification) |
   | `outlier-too-small` | no category reaches `centroid_min_points`: no slider, one line, *“No category has 1000 cells.”* (the full reason is its title) | `serve_anndata` on the `pancreas` store with `centroid_min_points=1000` (its largest cluster, Ductal, has 916 cells), colour by `clusters` |
   | `outlier-failed` | the latent space holds a NaN: no slider, the line *“This latent space cannot be used.”* with **Choose latent space**, and the latent space marked unusable in Data mapping | an `.h5ad` written with anndata itself (`prepare()` refuses non-finite embeddings) whose 10-column `obsm['X_bad']` holds one NaN, `cellucid serve … --latent-key X_bad`, move the slider |

   Each needs a scenario in `docs/_tooling/screenshots/scenarios.mjs` and a
   provenance record before it is referenced from a page.

### Resolved since the audit

- **Dimension switching and navigation mode** — `c_/01` and `c_/05` now
  document that each view owns its camera and mode, that a mode you chose is
  remembered per view and dimension, and that the default follows the
  dimension of the embedding the view shows
  (`viewer-contracts.js getDefaultNavigationMode`).
- **Multi-representation interface** — the Dimension dropdown and the
  `Compare views › Embedding` menu are gone: each view shows one embedding,
  chosen in `View N shows` (`#view-representation-select`), and keys 1, 2 and 3
  show a dimension's default, then cycle. Every frame captured before that
  (`window-dimension-1d/-2d/-3d`, `multiview-two-panels`,
  `connectivity-multiview`, `dimension-2d-planar-default`,
  `dimension-navigation-controls-2d`, `dark-theme-multiview`) shows the
  earlier interface; their captions say so until they are recaptured.
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

`Grid compare` lays the panels out in the part of the canvas the sidebar does
not cover, as equal cells starting at the sidebar's right edge (the whole width
when the sidebar is closed). The one DOM anchor per panel is its title chip in
`#view-title-layer`, centred on the panel it labels; `buildMultiviewScenarios()`
asserts the panel geometry through those chips before every multiview capture.

### DOM crop targets

All forty confirmed present in `cellucid/index.html`:

`#visualization-section` · `#render-mode` · `#depth-controls` · `#renderer-controls` ·
`#smoke-controls` · `#points-controls` · `#connectivity-controls` ·
`#compare-views-section` · `#split-view-controls` · `#split-keep-view-btn` ·
`#camera-lock-btn` · `#split-view-badges-box` · `#split-view-badges-list` ·
`#view-layout-mode` · `#split-clear-btn` · `#representation-controls` ·
`#view-representation-select` · `#navigation-controls` · `#navigation-mode` ·
`#orbit-controls` · `#planar-controls` · `#freefly-controls` · `#reset-camera-btn` ·
`#coloring-filtering-section` · `#categorical-field-row` · `#continuous-field-row` ·
`#gene-expression-row` · `#gene-expression-dropdown` · `#display-options-container` ·
`#outlier-filter-container` · `#centroid-controls` · `#legend` ·
`#active-filters-container` · `#filter-count` · `#active-filters` ·
`#deleted-fields-section` · `#cinematic-camera-section` · `#shortcuts-section` ·
`#glcanvas` · `#sidebar`

### Dataset facts

Read from each sample's identity in the sample catalog at audit time (prepared, before 0.9.5):

| Sample id | Cells | Genes | Dimensions | Default | Connectivity | Useful fields |
| --- | ---: | ---: | --- | ---: | --- | --- |
| `suo` | 561,947 | 5,103 | 1, 2, 3 | 3 | yes | `cell_type`, `organ`, `LVL0`–`LVL3`, `age`, `n_genes` |
| `garcia` | 219,731 | 5,754 | 1, 2, 3 | 3 | yes | same schema as suo |
| `he` | 71,650 | 5,152 | 1, 2, 3 | 3 | yes | same schema as suo |
| `kanemaru` | 131,636 | 3,691 | 1, 2, 3 | 3 | yes | same schema as suo |
| `pancreas` | 3,696 | 3,753 | 1, 2, 3 | 3 | yes | `clusters`, `clusters_coarse`, `cell_type`, `S_score`, `G2M_score` |

Because all five publish 1D, 2D and 3D, the "dataset provides only one embedding
dimension" state is **unreachable with any shipped sample**. Only the local
fixture `cellucid/tests/browser/fixtures/catalog/current-ui.h5ad` reproduces
it: 120 cells, **2D only**, one categorical field (`cell_type` =
alpha/beta/gamma; latent space `latent_space`, 3 columns, named as `latentKey`
with no stored outlier quantiles), one continuous field
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
