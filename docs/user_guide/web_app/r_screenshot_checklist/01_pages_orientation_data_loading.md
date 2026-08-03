# Screenshot coverage audit — sections A, B, G, Q

Scope: `user_guide/web_app/a_orientation/` (6 pages), `user_guide/web_app/b_data_loading/`
(13 pages), `user_guide/web_app/g_cross_highlighting/` (empty directory),
`user_guide/web_app/q_troubleshooting_index/` (1 page).

All paths below are relative to `cellucid-python/docs/` unless stated otherwise.
Read-only audit; nothing outside this file was modified.

---

## Summary

| Metric | Count |
| --- | ---: |
| Pages audited | 20 (6 A + 13 B + 1 Q; G contributes 0 pages) |
| Records below | 21 (20 pages + 1 for the empty G directory) |
| `{figure}` blocks in these pages | 18 |
| Distinct images they reference | 11 |
| Pages with at least one figure | 9 of 20 |
| Pages with zero figures | 11 of 20 |

| Verdict (primary) | Records |
| --- | ---: |
| OK | 1 |
| ADD | 6 |
| REPLACE | 4 |
| SEQUENCE | 4 |
| DIAGRAM | 2 |
| NONE | 4 |
| STALE-RISK | 0 as a sole verdict — see cross-cutting note 3; 3 pages carry an already-wrong build string alongside another verdict |

Records sum to 21. Each record is counted once under its **primary** verdict.
Four pages are hybrids: ``b/02_local_demo_tutorial`` and ``b/04_server_tutorial``
are counted under REPLACE but also need a SEQUENCE; ``a/03_quick_tour_60_seconds``
is counted under SEQUENCE but also needs one REPLACE.
{doc}`/user_guide/web_app/b_data_loading/09_screenshots` is counted under REPLACE because its contents
should be redistributed and the page deleted, not extended.

### Cross-cutting findings (apply to every page below)

1. **Every image in the entire docs tree is 1x device pixel ratio.** Measured
   with `sips`: `web_app/app-overview-cell-type.png` is 1440x900 and is
   displayed at `:width: 1440px`; `startup-loaded-build.png`,
   `welcome-startup.png` and both `jupyter/*.png` are 1440x1000 at
   `:width: 1440px`;
   `data_loading/data-loading-session-panel.png` is 246x560 at `:width: 246px`;
   the five `data_loading/h5ad-*` / `zarr-*` files are 1280x720 at
   `:width: 1280px`. Displayed 1:1 at native size, so on any HiDPI screen —
   which is most readers — they resample and blur. The user's first stated
   requirement is "high resolution". **Every recapture in this audit should use
   Playwright `deviceScaleFactor: 2` and keep the existing `:width:` value**,
   i.e. capture 2880x2000 and render at 1440px. This is a mechanical change that
   applies to all 10 images and needs no editorial decision.
2. **No image anywhere in the docs shows a mouse cursor**, and no image is
   annotated (no callouts, numbered badges, or arrows). There are also **no
   diagrams at all** in the docs — `grep` for `mermaid`, `graphviz`, and `.svg`
   across `user_guide/web_app/` returns only prose mentions of SVG *export*;
   the only SVG asset in `_static/` is the logo.
3. **Three captions hard-code a build string that is already wrong.** The
   current web build is `2026-07-27.24`
   (`cellucid/index.html:17`, `<meta name="cellucid-web-build-id" content="2026-07-27.24" />`).
   The docs claim `Build 2026-07-27.1` in {doc}`/user_guide/web_app/a_orientation/03_quick_tour_60_seconds` (line 31
   alt, line 36 caption), {doc}`/user_guide/web_app/a_orientation/04_ui_glossary_terminology` (line 48 alt, line 53
   caption), and {doc}`/user_guide/web_app/b_data_loading/10_standard_pancreas_dataset` (line 51 caption). A reader
   who follows the instruction to compare the footer will find a mismatch.
   Recommendation: drop the exact build number from alt text and captions and
   say "the footer prints the build identity" instead. Naming the build in prose
   guarantees the caption rots on the next deploy.
4. **Two images carry 8 of the 16 figure blocks in these sections and are reused
   docs-wide.** `web_app/app-overview-cell-type.png` appears on 8 pages total
   (4 of them mine); `data_loading/data-loading-session-panel.png` appears on
   10 pages total (4 of them mine). Reuse of a UI-map image is defensible; reuse
   of the *same* image as the "what success looks like" proof for three
   different loader paths is not.
5. **The gallery-page pattern is used in 10 web-app sections**
   (`b/09_screenshots.md`, `c/07`, `d/06`, `e/08`, `f/07`, `h/11`, `i/08`,
   `k/08`, `l/11`, `n/08`). See the analysis under
   {doc}`/user_guide/web_app/b_data_loading/09_screenshots` — the pattern should be retired.
6. **Two `_static/screenshots/` topic directories are empty**:
   `_static/screenshots/jupyter_hooks/` and `_static/screenshots/sessions_sharing/`.
   {doc}`/user_guide/web_app/l_sessions_sharing/11_screenshots` is therefore a gallery page holding
   zero images of its own — it reuses `data_loading/data-loading-session-panel.png`.
   Outside my sections, but it is the same defect and worth reporting.
7. **The style guide's own rule is being violated.**
   {doc}`/user_guide/python_package/h_developer_docs/15_docs_development_and_style_guide`
   (screenshot section, ~line 105) says: "confirm that the caption describes only
   what is visibly present." The `app-overview-cell-type.png` caption claims a
   legend that is not in the frame (see {doc}`/user_guide/web_app/a_orientation/03_quick_tour_60_seconds` below).

### DOM selectors verified in `cellucid/index.html` (for scripting captures)

Confirmed by reading the source, so `needs` entries below can name real targets:

- Sidebar accordions (`<details class="accordion-section">` → `<summary>`):
  `#session-section` "Session", `#visualization-section` "Visualization",
  `#compare-views-section` "Compare Views", `#coloring-filtering-section`
  "Coloring & Filtering", `#community-annotation-section` "Community Annotation",
  `#highlighted-cells-section` "Highlighting", `#page-analysis-section`
  "Analysis", `#cinematic-camera-section` "Camera Path", `#figure-export-section`
  "Figure Export", `#shortcuts-section` "Keyboard Shortcuts", `#benchmark-section`
  "Performance Benchmark".
- Session panel internals: `#dataset-info` (the DATASET/CELLS/GENES/OBS
  FIELDS/CONNECTIVITY block), `#dataset-select` ("Sample datasets:"),
  `#user-data-block` ("Local data:" with `#user-data-h5ad-btn` "H5AD",
  `#user-data-zarr-archive-btn` "Zarr ZIP", `#user-data-browse-btn` "Prepared"),
  `#remote-server-block` (`#remote-server-url`, `#remote-connect-btn` "Connect"),
  `#github-repo-block` (`#github-repo-url` placeholder
  `owner/repo/path/to/exports`, `#github-connect-btn` "Connect"),
  `#session-state-controls` (`#save-state-btn`, `#load-state-btn`).
- Other: `#welcome-modal` / `#welcome-demo-btn` ("Choose a dataset"), `#sidebar`,
  `#legend`, `#split-keep-view-btn` ("Keep view"),
  `.highlight-mode-btn[data-mode="lasso"]` ("Lasso"), `#navigation-controls`
  (inside Compare Views), `#web-build-version` (footer build span).
- Loadable fixtures under `cellucid/tests/browser/fixtures/`:
  `current-ui-smoke.h5ad` (120 cells, 6 genes), `current-ui-smoke.zarr.zip`,
  `exports/current-ui-prepared/`, `demo-custom-exports/` (three synthetic
  datasets: `synthetic-cell-types-2d`, `synthetic-development-3d`,
  `synthetic-trajectory-1d`), and `broken-exports/` (9 deliberately invalid
  exports: `corrupt-payload`, `mismatched-manifest`, `missing-payload`,
  `truncated-catalog`, `wrong-format-catalog`, and mixed variants).
- Catalog sample ids confirmed in `cellucid-datasets/exports/datasets.json`:
  `suo` (default), `garcia`, `he`, `kanemaru`, `pancreas`.

---

# A — `user_guide/web_app/a_orientation/`

## user_guide/web_app/a_orientation/index.md

- **path** — {doc}`/user_guide/web_app/a_orientation/index`
- **explains** — Landing card grid and recommended reading order for the four
  orientation pages, with a globbed `toctree`.
- **has** — none.
- **verdict** — NONE
- **needs** — n/a.
- **notes** — A pure navigation page; the `{grid-item-card}` octicons already
  give it visual structure. Adding a screenshot here would duplicate the quick
  tour one click away. One inconsistency worth noting: the reading order
  (lines 11–15) puts ``01_what_is_cellucid`` **last**, after
  ``02_system_requirements``, while the filename numbering puts it first. A
  wet-lab reader arriving at `01_...` by its number gets the conceptual page
  before the tour. Not a screenshot issue, but it undercuts the "wet-lab entry
  point" goal.

## user_guide/web_app/a_orientation/01_what_is_cellucid.md

- **explains** — What Cellucid is and is not, the dataset → embedding → field →
  view → artifact mental model, its relationship to `cellucid-python`, and the
  three persistence layers (in-memory, session bundle, `localStorage`).
- **has** — none. **There is a hole where a diagram is promised.** Line 68 reads
  "This diagram is the quickest way to understand 'what is a thing' in Cellucid
  and what changes when you click." Lines 69, 70, 71 are empty, and line 72
  resumes with "Key terms". The sentence has no referent. This is a broken
  reference, not an omission of taste.
- **verdict** — DIAGRAM
- **needs**
  1. `_static/diagrams/mental-model-object-graph.svg` — the missing artefact at
     line 69. A schematic, not a screenshot: boxes for **Dataset** →
     **Embedding (1D/2D/3D)** → **View** (with **Live view** and **Snapshot
     view** as two instances), with **Field** (categorical obs / continuous obs
     / gene expression), **Filter**, and **Highlight group** attaching to a
     View, and **Artifacts** (`.cellucid-session`, exported figure) leaving the
     View. Annotate each edge with scope — "per view" on field/filter/highlight,
     "global" on dataset/embedding — because the page's own "Scope note"
     paragraphs in {doc}`/user_guide/web_app/a_orientation/04_ui_glossary_terminology` depend on that distinction.
     A screenshot cannot express this; no arrangement of the UI shows the object
     graph. Must be authored as an SVG (or a MyST `mermaid` block if the theme
     is configured for it — it currently is not; no mermaid usage exists
     anywhere in the docs).
  2. Second, smaller diagram for the "State & persistence" section (line 109):
     three nested boxes labelled **In-memory state** (lost on refresh),
     **Session bundle** (`.cellucid-session`, explicit Save State),
     **localStorage preferences** (theme/background only), with a "refresh"
     arrow showing what survives. Suggested
     `_static/diagrams/persistence-layers.svg`.
- **notes** — The tab-set at lines 21–62 gives three depths of the same
  explanation; that is good for the wet-lab reader and needs no image. The page
  is the only one in the section with a genuine content hole rather than a
  coverage gap — fixing the empty diagram slot should rank above any recapture
  work in this section.

## user_guide/web_app/a_orientation/02_system_requirements.md

- **explains** — That WebGL2 is a hard requirement, which browsers/OS
  combinations are supported, what hardware matters (VRAM over CPU), and four
  long symptom-driven troubleshooting blocks (WebGL2 unsupported, blank canvas,
  context lost, general slowness).
- **has** — none.
- **verdict** — ADD
- **needs**
  1. **The WebGL2 failure overlay itself.** Load the app with WebGL2 forced off
     — `cellucid/tests/browser/webgl-harness.html` exists as a fixture for
     exactly this, and `tests/browser/startup-failure.spec.mjs` already
     exercises a startup-failure path, so the state is scriptable. Full window
     1440x900 (capture at 2x), cropped to the error overlay plus enough
     surrounding chrome to prove it is Cellucid. This is the single highest-value
     image on the page: the prose quotes the exact string `WebGL2 is required
     but not supported in this browser.` (line 18) and a reader who is *seeing*
     that string needs to confirm they are on the right page in one glance.
     Cursor: not useful. Filename `_static/screenshots/web_app/webgl2-required-overlay.png`.
  2. **The "WebGL context lost. Reload required" overlay** (section at line 177).
     Same harness; the app shows a distinct overlay for this. Crop to the overlay.
     Filename `_static/screenshots/web_app/webgl-context-lost-overlay.png`.
     Without it, a reader cannot tell symptom 1 from symptom 3, which is the
     whole navigational purpose of the section.
  3. **Blank canvas that is actually a filter problem** (section at line 149,
     cause 2: "the dataset is loaded but all points are invisible"). Load `suo`,
     open **Coloring & Filtering** (`#coloring-filtering-section`), drag the
     outlier threshold until the visible count reads `0`. Crop to the whole
     window so both the empty canvas and the `0 visible` readout are in frame —
     the pairing *is* the lesson. Cursor on the threshold slider handle.
     Filename `_static/screenshots/web_app/blank-canvas-zero-visible.png`.
- **notes** — The self-check at line 71 tells the reader to evaluate
  `document.createElement("canvas").getContext("webgl2") !== null` "in any
  supported browser console". A wet-lab reader does not know what a browser
  console is, how to open it, or where to type. This is the clearest place in my
  four sections where the docs fail the stated "even a wet lab person should
  understand and use the app" bar. Either add a browser-devtools screenshot
  (Chrome DevTools Console with the expression typed and `true` returned —
  cursor not needed, but the typed line and the `true` result must both be
  legible), or replace the instruction with a link to an external WebGL report
  page. I recommend the screenshot: `_static/screenshots/web_app/devtools-webgl2-check.png`,
  cropped to the DevTools Console panel only, at 2x.

## user_guide/web_app/a_orientation/03_quick_tour_60_seconds.md

- **explains** — A six-step click-through: dismiss the welcome overlay and
  confirm/replace the default dataset, move the camera, color by a field, keep a
  snapshot, make a lasso highlight, then export a figure or save a session.
- **has** — three figures, two of them shared:
  - line 15 `web_app/welcome-startup.png`, `:width: 1440px` — native 1440x1000.
    Used **only here**. Appropriate.
  - line 30 `web_app/startup-loaded-build.png`, `:width: 1440px` — native
    1440x1000. **Also used on** `04_ui_glossary_terminology.md:47`.
  - line 45 `web_app/app-overview-cell-type.png`, `:width: 1440px` — native
    1440x900. **Reused on 8 pages docs-wide**:
    `a_orientation/03`, `a_orientation/04:22`, `b_data_loading/02:196`,
    `b_data_loading/03:81`, {doc}`/user_guide/web_app/d_fields_coloring_legends/06_screenshots`,
    {doc}`/user_guide/python_package/a_landing_pages/01_what_is_cellucid_python`,
    {doc}`/user_guide/python_package/f_notebooks_tutorials/21_prepare_exports_with_quantization_and_compression`,
    {doc}`/user_guide/python_package/g_api_reference_coverage/api/export`.
- **verdict** — SEQUENCE (plus one REPLACE)
- **needs**
  1. **REPLACE `app-overview-cell-type.png` at line 45.** I opened the image.
     The heading above it says "One screenshot 'map' (recommended)" and the
     prose says "A single **annotated** screenshot makes the rest of the docs
     much easier to follow" (line 43) — **the image is not annotated**. Worse,
     its caption on every one of the 8 pages reads "the sidebar controls the
     active view while the categorical legend maps directly to the colored
     points", and **there is no legend in the frame**: the sidebar is scrolled to
     Session + the top of Visualization, and `#legend` is far below the fold.
     The filename says "cell-type" but no coloring control is visible to
     corroborate that either. Replacement: load `suo`, scroll `#sidebar` so that
     **Coloring & Filtering** and `#legend` are both visible, full window
     1440x900 at 2x, then composite numbered callout badges (1–6) over the
     regions the tour names in order — 1 Session accordion, 2 Compare Views →
     Navigation, 3 Coloring & Filtering, 4 Compare Views → Multiview / Keep
     view, 5 Highlighting, 6 Figure Export + Save State — with a matching
     numbered key in the caption. Cursor: no (callouts replace it). Filename
     `_static/screenshots/web_app/ui-map-annotated.png`. This one image is the
     wet-lab reader's whole orientation; it is currently the weakest asset in the
     section.
  2. **SEQUENCE: one capture per numbered step**, because the page is a
     click-by-click tour and currently illustrates none of the six steps. All on
     `suo` unless noted, full window 1440x900 at 2x, cropped as stated:
     - Step 1 (line 56): after clicking `#welcome-demo-btn` "Choose a dataset",
       crop to `#session-section`, cursor on `#dataset-select`.
       `_static/screenshots/web_app/tour-01-choose-dataset.png`.
     - Step 2 (line 76): open **Compare Views**, crop to `#navigation-controls`
       showing the mode select set to **Orbit** with the orbit sub-controls
       expanded, cursor on the Navigation select.
       `_static/screenshots/web_app/tour-02-navigation.png`.
     - Step 3 (line 92): open **Coloring & Filtering**, pick the `cell_type`
       categorical field; crop to a vertical strip containing
       `#coloring-filtering-section` **and** `#legend` populated with categories,
       cursor on the field selector. `_static/screenshots/web_app/tour-03-color-by.png`.
       (Do not reuse `filtering/coloring-filtering-cell-type-panel.png` — that
       one is already spread across 6 pages in sections D and E.)
     - Step 4 (line 105): after clicking `#split-keep-view-btn`, full window
       showing two panels and the new view badge, cursor on "Keep view".
       `_static/screenshots/web_app/tour-04-keep-view.png`. Note
       `web_app/multiview-two-panels.png` already exists (1440x900) and is used
       in section C and N — check it before capturing a near-duplicate; if it
       shows the badge clearly, reuse it here and skip the capture.
     - Step 5 (line 113): `.highlight-mode-btn[data-mode="lasso"]` active, a
       lasso path drawn mid-gesture over a cluster, crop to the canvas plus
       `#highlighted-cells-section`, cursor at the lasso path head — **this is
       the single best cursor candidate in the section**, because a lasso is a
       drag and a static frame without a pointer does not read as a gesture.
       `_static/screenshots/web_app/tour-05-lasso-highlight.png`.
     - Step 6 (line 123): crop to `#session-state-controls` with cursor on
       `#save-state-btn`. `_static/screenshots/web_app/tour-06-save-state.png`.
- **notes** — Stale build string, see cross-cutting note 3: line 31 alt and
  line 36 caption both assert `Build 2026-07-27.1`; the app ships
  `2026-07-27.24`. Line 130 tells the reader to "scroll to the sidebar footer and
  include its exact **Build** value" — that instruction deserves a small crop of
  the footer (`#web-build-version`) so a beginner knows what they are looking
  for; suggest `_static/screenshots/web_app/sidebar-footer-build.png`, cropped
  to the footer line only, ~246px wide at 2x, with the build value legible but
  **not** repeated in the caption.

## user_guide/web_app/a_orientation/04_ui_glossary_terminology.md

- **explains** — Definitions of the app's vocabulary grouped by area: startup and
  build identity, the compact `i` help buttons, data model, fields/coloring/
  legends, views and snapshots, navigation and camera, highlighting, analysis and
  export, and the global/per-view/session/localStorage scope words.
- **has** — two figures, both shared:
  - line 22 `web_app/app-overview-cell-type.png` `:width: 1440px` — same
    8-page-reused image described above, under the heading "UI map (one
    screenshot you'll reuse everywhere)".
  - line 47 `web_app/startup-loaded-build.png` `:width: 1440px` — also on
    `03_quick_tour_60_seconds.md:30`.
- **verdict** — REPLACE
- **needs**
  1. **Swap line 22 to the annotated `ui-map-annotated.png` produced for
     {doc}`/user_guide/web_app/a_orientation/03_quick_tour_60_seconds`.** Two pages call the same image "the UI map"
     and "one screenshot you'll reuse everywhere"; the callout-annotated version
     is the one that earns that description. No new capture — reuse is correct
     *here*, because a glossary genuinely wants the same reference frame as the
     tour.
  2. **ADD one crop per terminology cluster that a beginner cannot picture from
     text.** Prioritised, all `suo`, 2x:
     - **View badge** (line 115: "the clickable pill/row representing a view
       (with indicators like `3D`, `Orb/Pan/Fly`)"). This is pure jargon without
       a picture. Keep two views, crop tightly to the badge row so both badges
       and their per-view indicators are legible; cursor on the inactive badge
       to show it is clickable.
       `_static/screenshots/web_app/glossary-view-badges.png`.
     - **Cameras locked vs unlocked** (lines 116–117). Two crops of the same
       control in its two states, or one crop with the toggle plus two badges
       showing differing nav indicators.
       `_static/screenshots/web_app/glossary-cameras-locked.png` /
       `-unlocked.png`.
     - **Categorical legend vs continuous legend** (lines 91–93). Two narrow
       `#legend` crops side by side — categories with swatches, and a numeric
       colour scale. `_static/screenshots/web_app/glossary-legend-categorical.png`
       and `-continuous.png`, ~246px display width each.
     - **The `i` help button open state** (lines 58–65). The page describes the
       dialog's open/close behaviour in detail but never shows one. Click
       `#user-data-info-btn`; crop to `#user-data-block` with
       `#user-data-info-tooltip` visible; cursor on the `i` button.
       `_static/screenshots/web_app/glossary-info-tooltip-open.png`. There is a
       browser test (`tests/browser/info-popovers.spec.mjs`) that already drives
       this state.
- **notes** — Stale build string in line 48 alt and line 53 caption (see
  cross-cutting note 3). The glossary is the right place to *define* "Build" but
  the wrong place to *quote* a build number.

## user_guide/web_app/a_orientation/05_which_workflow_is_for_me_decision_tree.md

- **explains** — How to pick between prepared-folder loading, `cellucid-python`,
  and Community Annotation, via a plain-language tab-set, an ASCII decision tree
  for computational users, and a "when to use what" table with gotchas.
- **has** — none. **A second suspicious gap**: lines 15 and 16 are empty between
  the last bullet of "The short answer" (line 13) and the `---` at line 17,
  matching the shape of the confirmed missing-diagram hole in
  {doc}`/user_guide/web_app/a_orientation/01_what_is_cellucid`. Whether or not something was removed there, the
  section is the natural home for the visual form of the decision tree.
- **verdict** — DIAGRAM
- **needs**
  1. `_static/diagrams/workflow-decision-tree.svg` — replace the ASCII tree at
     lines 42–51 with a real flowchart, and place it at line 15 so both the
     wet-lab and computational tabs can point at one artefact. Entry question
     "What do you have?" branching to: **prepared `export/` folder** → Local data
     → Prepared; **`.h5ad` / `.zarr`** → server mode or Jupyter; **AnnData in a
     notebook** → `show_anndata()`; **nothing yet** → ask a colleague to
     `prepare()`. Then a second tier on sharing: **one-to-one / reproducible** →
     `.cellucid-session`; **public, no server** → GitHub exports;
     **many-to-many labelling** → Community Annotation. Each leaf should carry
     the page reference it currently links to in prose. The ASCII block is
     invisible to the wet-lab tab (it lives only in the "Computational / Power
     User" tab) and unreadable on a phone.
- **notes** — No screenshot belongs on this page; it is a routing page and any
  UI capture would date faster than the routing logic. The table at line 66 says
  "Compare multiple hypotheses → Multiview snapshots ... Common gotcha: Smoke
  mode is disabled when snapshots exist" — this is the only statement of that
  constraint in section A, and it is also asserted in
  `03_quick_tour_60_seconds.md:148` ("You are likely in smoke mode (multiview is
  points-only)"). Consistent, so no defect, but worth a diagram annotation.

---

# B — `user_guide/web_app/b_data_loading/`

Section-level note on the five loader paths and the user's explicit request for
terminal and Jupyter captures:

| Loader path | Page | Terminal shot needed? | Jupyter shot needed? | Covered today |
| --- | --- | --- | --- | --- |
| Public sample / local demo | ``02_local_demo_tutorial`` | no | no | reused generic image only |
| Browser file picker | ``03_browser_file_picker_tutorial`` | no | no | reused generic images only |
| Server / CLI | ``04_server_tutorial`` | **yes — 3 more** | no | 1 shot, and it shows a different command than the step it sits under |
| Jupyter | ``05_jupyter_tutorial`` | one (`ssh -L` tunnel) | **yes — 4** | **zero figures on the page** |
| GitHub custom repo | ``11_custom_dataset_repository`` | one (`cellucid serve ./exports`) | no | **zero figures** |

## user_guide/web_app/b_data_loading/index.md

- **explains** — Routes the reader to the right loader page via a "you have… →
  best first choice" table, card grids for getting started / viewing methods /
  concepts, and a method-vs-server comparison table.
- **has** — line 41 `data_loading/data-loading-session-panel.png`, `:width: 246px`
  (native 246x560). **Reused on 10 pages docs-wide**: this page,
  `01_loading_options_overview.md:254`, `03_browser_file_picker_tutorial.md:49`,
  `09_screenshots.md:13`, `l_sessions_sharing/{index,01_session_mental_model,
  03_save_restore_ux,11_screenshots}.md`,
  {doc}`/user_guide/python_package/b_concepts_mental_models/03_state_persistence_and_scope`,
  {doc}`/user_guide/python_package/c_data_preparation_api/01_prepare_export_overview`,
  {doc}`/user_guide/python_package/g_api_reference_coverage/api/sessions`.
- **verdict** — OK
- **needs** — Recapture at 2x only (cross-cutting note 1). The image is the right
  image for this page: an index that lists five loader paths should show the one
  panel that contains all five. I opened it and confirmed it shows DATASET/
  DESCRIPTION/SOURCE/URL/CELLS/GENES/OBS FIELDS/CONNECTIVITY, SAMPLE DATASETS,
  LOCAL DATA (H5AD / ZARR ZIP / PREPARED), REMOTE SERVER, GITHUB DATA, SESSION
  STATE — every path in the table above.
- **notes** — This is the one place the 10-way reuse of this image is *correct*.
  Consider adding coloured region labels ("① sample ② local ③ remote ④ GitHub
  ⑤ session") keyed to the Fast Path table rows, which would make it a genuine
  index graphic rather than a bare panel photo.

## user_guide/web_app/b_data_loading/01_loading_options_overview.md

- **explains** — The canonical 14-row loading-options matrix, what the three
  input formats mean, what lazy loading is and which modes actually have it,
  where vector fields come from per format, plus copy/paste install/export/serve/
  notebook snippets and a compact troubleshooting index.
- **has** — line 254 `data_loading/data-loading-session-panel.png` `:width: 246px`
  — the 10-page reused image, placed under "Interface reference" at line 250.
- **verdict** — ADD
- **needs**
  1. **A rendering of the 14-option matrix as a decision graphic** is *not* what
     I recommend — the table at line 131 is genuinely better as a table. Instead:
  2. **One "lazy vs not lazy" comparison capture** for the section at line 64,
     which is the page's central concept and is currently pure prose. Load
     `pancreas` from the sample catalog, open DevTools Network, filter to the
     dataset origin, and search a gene in **Coloring & Filtering**; crop to the
     Network panel showing the single on-demand gene request arriving *after*
     the initial load. Then the contrasting shot: select
     `cellucid/tests/browser/fixtures/current-ui-smoke.h5ad` through
     `#user-data-h5ad-btn` and capture the Network panel showing the whole file
     fetched up front. Two images, side by side, make the "❌* Browser `.h5ad`
     loading is **not truly lazy**" footnote at line 148 concrete. Cursor: no.
     `_static/screenshots/data_loading/lazy-gene-request.png` and
     `_static/screenshots/data_loading/h5ad-eager-whole-file.png`.
  3. Keep the Session panel figure but move it **above** the matrix rather than
     leaving it at line 254 after the code blocks — options #1–#5 in the matrix
     are all clicks in that panel, and the reader meets the matrix 120 lines
     before they see what it refers to.
- **notes** — The page title says "All 14 Options" and the matrix has exactly 14
  rows — verified. The 512 MiB browser limit cited in the sibling page is real:
  `cellucid/assets/js/data/h5ad.js:163` `MAX_H5AD_BROWSER_FILE_BYTES = 512 * 1024 * 1024`
  and `data-source.js:237` `MAX_PREPARED_BROWSER_BYTES`. This page does not state
  the limit anywhere in the matrix, though it is the deciding factor for rows
  #4 and #5; consider adding it to the footnote at line 148.

## user_guide/web_app/b_data_loading/02_local_demo_tutorial.md

- **explains** — Export a dataset once with `prepare()`, generate
  `datasets.json` with `generate_datasets_manifest()`, validate locally, then
  publish the `exports/` root to a public GitHub repo and share a
  `?github=owner/repo/exports&dataset=id` link; plus the optional local-demo
  `exportsBaseUrl` route.
- **has** — line 196 `web_app/app-overview-cell-type.png` `:width: 1440px` — the
  8-page reused Suo 3D overview, placed under "Step 3 — Validate Locally (Before
  You Upload)".
- **verdict** — REPLACE + SEQUENCE
- **needs**
  1. **REPLACE line 196.** The step is "validate *your PBMC export* by opening it
     with the **Prepared** picker"; the image shows the bundled Suo sample in a
     3D orbit view with no picker interaction and no legend. It proves nothing
     about the step. Replacement: click `#user-data-browse-btn` **Prepared**,
     select `cellucid/tests/browser/fixtures/exports/current-ui-prepared/`, and
     capture the settled result — crop to `#session-section` so the
     `#dataset-info` rows (DATASET, CELLS, GENES, OBS FIELDS) read out the
     just-validated export, with the success notification toast in frame if it
     is still visible. Cursor on `#user-data-browse-btn`.
     `_static/screenshots/data_loading/validate-prepared-export.png`.
  2. **SEQUENCE for Step 4 — Publish to GitHub** (line 203), which is where a
     wet-lab reader is most likely to fail and where there is currently no
     visual aid at all:
     - **The exports root on disk**, so the reader can compare their folder to
       the required layout in the code block at line 58. A file-manager or
       `tree` terminal capture of `exports/` with `datasets.json` at the root and
       two dataset folders beside it. Terminal is better — it shows the exact
       command. `_static/screenshots/data_loading/exports-tree-terminal.png`,
       cropped to the terminal, command `tree -L 2 exports/` visible with output.
     - **The GitHub connect form filled in**, cropped to `#github-repo-block`
       with `#github-repo-url` containing
       `theislab/cellucid-demo-custom-datasets/exports` and cursor on
       `#github-connect-btn` ("Connect"). This is the exact field described at
       line 234 ("What repo path to enter in Cellucid") and it has no picture.
       `_static/screenshots/data_loading/github-connect-form.png`, ~246px.
     - **The validated dataset list after connecting**, cropped to
       `#dataset-select` expanded showing the three synthetic entries — proof of
       the "show the validated dataset list" claim at line 244.
       `_static/screenshots/data_loading/github-connected-dataset-list.png`.
     - **The failure state**, cropped to `#github-repo-block` after connecting to
       a path with no `datasets.json`, showing the exact error surfaced by the
       app. The troubleshooting block at line 295 quotes the message
       `datasets.json not found` — a reader matching their screen to the docs
       needs to see it rendered.
       `_static/screenshots/data_loading/github-datasets-json-not-found.png`.
- **notes** — The privacy warning at line 39 ("you are effectively publishing
  your dataset") is the most consequential paragraph in the section and is pure
  text in the middle of a long page. Consider an admonition-styled callout image
  or, better, promote it into a `:::{danger}` directive — it is currently plain
  prose under an `##` heading and reads as optional.

## user_guide/web_app/b_data_loading/03_browser_file_picker_tutorial.md

- **explains** — The three browser-only loading controls (Prepared folder, H5AD
  file, Zarr ZIP archive), when each applies, the hard browser-memory limit on
  direct H5AD, minimum `obsm` requirements, and a long symptom list.
- **has** — two figures, both shared:
  - line 49 `data_loading/data-loading-session-panel.png` `:width: 246px` —
    the 10-page reused panel.
  - line 81 `web_app/app-overview-cell-type.png` `:width: 1440px` — the 8-page
    reused Suo overview, under "Option #3 … What success looks like".
- **verdict** — SEQUENCE
- **needs** — This page has three distinct controls and shows a picture of
  none of them being used. It is the wet-lab reader's most likely first
  self-service path and it needs a per-option sequence:
  1. **Prepared (Option #3, line 56).** Two frames.
     (a) Cursor on `#user-data-browse-btn`, cropped to `#user-data-block`, with
     the button in hover state — the "What you click" instruction at line 68 is
     one sentence with no picture.
     (b) The OS directory chooser open on
     `cellucid/tests/browser/fixtures/exports/current-ui-prepared/`, cropped to
     the dialog. The page's own troubleshooting says "The file picker won't let
     me select a folder" and "The **Prepared** control accepts directories only"
     (line 168) — showing a *directory* chooser rather than a *file* chooser is
     the whole point and cannot be conveyed in prose.
     `_static/screenshots/data_loading/picker-prepared-button.png`,
     `_static/screenshots/data_loading/picker-prepared-directory-dialog.png`.
  2. **H5AD (Option #4, line 88).** Cursor on `#user-data-h5ad-btn`; then the
     loaded result cropped to `#session-section` showing SOURCE `H5AD file`.
     `_static/screenshots/data_loading/picker-h5ad-loaded.png`. Note the existing
     `data_loading/h5ad-current-loaded.png` already shows this state but is
     stranded on the gallery page {doc}`/user_guide/web_app/b_data_loading/09_screenshots` — **move it here** rather
     than capturing a duplicate (see that page's entry).
  3. **Zarr ZIP (Option #5, line 138).** Same treatment with
     `#user-data-zarr-archive-btn` and
     `cellucid/tests/browser/fixtures/current-ui-smoke.zarr.zip`; the existing
     `data_loading/zarr-zip-loaded.png` should move here.
  4. **The oversize-H5AD refusal.** The page states the 512 MiB ceiling at
     line 41 of the sibling page and warns "large `.h5ad` files can freeze the
     tab" (line 101). The app has an explicit refusal message —
     `cellucid/assets/js/data/h5ad.js:192`: "H5AD direct browser files must have
     a positive safe size no larger than 512 MiB; use the Cellucid server or
     prepared format". Capture that error as rendered, cropped to the
     notification. `_static/screenshots/data_loading/h5ad-too-large-refusal.png`.
     This turns a scary-sounding warning into a recognisable, recoverable state.
- **notes** — The "Recommended alternative for larger `.h5ad`" block at line 104
  gives a terminal command and the resulting URL `http://127.0.0.1:8765/?anndata=true`.
  That URL is correct for direct AnnData and matches
  `04_server_tutorial.md:59`. Verified consistent.

## user_guide/web_app/b_data_loading/04_server_tutorial.md

- **explains** — Options #6–#11: run `cellucid serve` over a prepared folder,
  `.h5ad`, or `.zarr`; the localhost security model; every CLI flag; the Python
  `serve()` / `serve_anndata()` equivalents; and the SSH-tunnel workflow for
  remote machines.
- **has** — no terminal figure. The page's four `cellucid serve` transcripts are
  verbatim ```` ```text ```` blocks; only `server/open-served-dataset.png` (the
  browser result) and `server/serve-cli-help.png` are images.
- **verdict** — SEQUENCE
- **needs**
  1. ~~Recapture the prepared-export startup in a plain terminal.~~
     **WITHDRAWN — terminal captures are forbidden.** The audit recommendations
     below originally asked for four more terminal images. They must not be
     taken. `cellucid serve` prints the **absolute** path of the directory it
     resolved (`server/_server.py`, `print_detail("Path", str(self.data_dir))`)
     and AnnData logs the absolute file path beside it, so an image of these
     commands contains a host path by construction — a neutral prompt does not
     help, and `strings` cannot audit rendered pixels. Two such images shipped
     with a home directory and a pyenv path in them, and two more with a scratch
     path; all four were removed and transcribed. See "Never photograph a
     terminal. Transcribe it." in {doc}`00_capture_tooling_and_conventions`.
  2. ~~Move the prepared-catalog capture to Option #6.~~ **DONE as text.** The
     transcript now sits under "Fast Path (CLI)" with `?source=remote`, and the
     direct-AnnData transcript under Option #7/#8 with `?anndata=true`, so the
     two URLs can no longer be confused for one another.
  3. ~~Add the SSH-tunnel terminal pair.~~ **WITHDRAWN, and the worst of the
     four:** `ssh -L 8765:localhost:8765 user@remote-host` puts a username and a
     hostname in the command itself. Document the tunnel as two `text` blocks
     labelled "on the server" and "on your laptop".
  4. **ADD the health-endpoint check** referenced at line 336
     (`http://127.0.0.1:8765/_cellucid/health`). A browser tab showing the small
     JSON response, cropped to the tab plus response. The troubleshooting step
     says "you should get a small JSON response" — showing it removes the
     ambiguity about what "small JSON" looks like.
     `_static/screenshots/server/health-endpoint-json.png`.
- **notes** — `cellucid serve --help` is offered by the page and its output then
  re-described by hand below it. `server/serve-cli-help.png` now carries the
  real output. It is the one terminal image that survives, because `--help`
  prints only the parser's own text — no resolved path, no prompt beyond a
  directory name — so it cannot leak. It is also the one that would most benefit
  from becoming a ```` ```text ```` block: 879 KB of PNG for text a reader
  cannot copy, search, or diff.

## user_guide/web_app/b_data_loading/05_jupyter_tutorial.md

- **explains** — Options #12–#14: `show()` for a prepared directory and
  `show_anndata()` for AnnData/`.h5ad`/`.zarr` inside a notebook; how the
  localhost server + same-origin iframe work; the exact `show_anndata` keyword
  signature; vector fields in notebooks; `viewer.highlight_cells` /
  `set_color_by` / `set_visibility` control; the `@viewer.on_*` hooks; pulling
  `viewer.state` and session bundles back into AnnData; `debug_connection()`;
  cleanup; and remote/HPC tunnelling.
- **has** — **none.** At 26 KB this is the largest page in the section and the
  only loader page with zero figures. It defers to
  ``../../python_package/f_notebooks_tutorials/05_jupyter_embedding_hooks_sessions_gallery``
  (tip at line 14, next-steps at line 786), which owns both existing Jupyter
  images: `jupyter/pancreas-notebook-embed.png` and
  `jupyter/pancreas-debug-connection.png` (1440x1000 each). Neither is
  referenced from the web-app section at all.
- **verdict** — SEQUENCE — **the single biggest gap in my four sections.**
- **needs** — The user asked specifically for Jupyter captures on the loader
  paths. This page is that path and has nothing. Four captures, all in
  JupyterLab at 1440x1000, 2x, using the `pancreas` dataset so the output is
  recognisable single-cell data rather than a toy fixture:
  1. **The minimal embed.** A notebook cell containing exactly the "Minimal:
     show an in-memory `AnnData`" snippet from line 66
     (`viewer = show_anndata(adata, height=600, dataset_name=…, dataset_id=…)`),
     with the rendered Cellucid iframe below it. Crop to the cell **and** its
     output together — the input/output pairing is the lesson; an image of the
     viewer alone teaches nothing about Jupyter. Cursor: no.
     `_static/screenshots/jupyter/show-anndata-minimal.png`.
     `jupyter/pancreas-notebook-embed.png` already covers roughly this state and
     could be **reused here immediately** at zero capture cost — I opened it and
     it shows the file browser, the notebook, the embedded viewer, and cells `[3]`
     and `[4]` with `viewer.wait_for_ready` / `set_color_by` / `highlight_cells`
     and their printed `ViewerState`. Reusing it is strictly better than the
     current zero. If recapturing, crop out the left file-browser panel — it
     shows unrelated tutorial filenames and wastes a third of the frame.
  2. **A hook firing.** Cell with `@viewer.on_selection` (line 450) registered,
     then a lasso drawn in the embedded viewer, then the printed
     `Selected: <n>` output visible below. Crop to the notebook column showing
     decorator cell, viewer, and output. Cursor at the lasso path in the iframe —
     this is the second strong cursor case in my sections, because the image has
     to convey "I dragged here and Python printed that".
     `_static/screenshots/jupyter/on-selection-hook-fired.png`.
  3. **`viewer.debug_connection()` output** (line 519). The page describes what
     the report checks in a 4-line paragraph and never shows it.
     `jupyter/pancreas-debug-connection.png` exists and is unreferenced from
     this section — **reuse it here.** Verify first that its content matches the
     described checks (`/_cellucid/health`, `/_cellucid/info`,
     `/_cellucid/datasets`, ping/pong, event counts, frontend debug snapshot).
  4. ~~The remote/HPC tunnel (line 578), as a notebook + laptop-terminal
     composite.~~ **WITHDRAWN — the terminal half cannot be photographed.**
     `ssh -N -L 8765:127.0.0.1:8765 <user>@<remote-host>` carries a username and
     a hostname, and a capture would show the real ones. Show the notebook cell
     printing `viewer.viewer_url` as an image if it helps, and the `ssh`
     invocation as a ```` ```text ```` block beside it.
- **notes** — The page's own tip (line 14) sends the reader to the Python-package
  section for "real captures", which means the web-app data-loading section
  outsources its only visual evidence to a different guide. That is a structural
  problem, not just a missing image: a reader in `b_data_loading/` following the
  five loader paths hits four pages with pictures and one that says "the
  pictures are over there". At minimum, cross-reference the two existing
  `jupyter/*.png` files inline here; they cost nothing and already exist.

## user_guide/web_app/b_data_loading/06_dataset_identity_why_it_matters.md

- **explains** — That the dataset id is the app's primary key; where it comes
  from per loading method; what must stay stable (id, cell ordering, field names)
  versus what may change (name, description, timestamps, compression); how to
  choose and verify one; and the session/annotation failures caused by changing
  it.
- **has** — none.
- **verdict** — ADD
- **needs**
  1. **The mismatch failure as the app renders it** — the page's two
     troubleshooting symptoms ("My session won't restore / says dataset
     mismatch", line 175; "Community annotation shows 'no votes' / 'wrong
     dataset'", line 191) are exactly the states a reader arrives here with,
     and neither is shown. Load `pancreas`, then `#load-state-btn` with a
     `.cellucid-session` saved against a different dataset id; crop to the
     resulting refusal message. Cursor on `#load-state-btn`.
     `_static/screenshots/data_loading/session-dataset-mismatch.png`.
     This is the one capture that turns an abstract page into a diagnostic one.
  2. **Where the id is actually visible in the UI.** Crop to `#dataset-info`
     with the DATASET row legible, ~246px, so a reader can answer "what is my
     dataset id right now" without opening JSON.
     `_static/screenshots/data_loading/dataset-info-identity-row.png`.
- **notes** — The page is otherwise correctly text-heavy; identity is a concept,
  not a screen. Do not over-illustrate it. The `dataset_id` character rules at
  lines 53–58 (1–180 ASCII, letter/digit first, no trailing `.`, no reserved
  Windows device name) match the rules restated in
  `02_local_demo_tutorial.md:141` and `05_jupyter_tutorial.md:260` — verified
  consistent across all three.

## user_guide/web_app/b_data_loading/07_folder_file_format_expectations_high_level_link_to_spec.md

- **explains** — The four things Cellucid can load (prepared folder, GitHub
  exports root, `.h5ad`, Zarr v2 store), the required and typical files inside a
  prepared export, the `datasets.json` catalog rules, the minimum `obsm`
  embedding keys, and the vector-field naming/`vectors/` layout with its identity
  schema.
- **has** — none.
- **verdict** — NONE
- **needs** — n/a.
- **notes** — Correctly figure-free. The page's content is directory trees and
  JSON schemas, both of which are already rendered as code blocks that copy,
  search, and diff — a screenshot of a file tree would be strictly worse. The one
  thing that *could* be visual is the relationship between
  `dataset_identity.json["vector_fields"]` and the integer-indexed files under
  `vectors/` (lines 156–190), because the "`<index>` is the vector field's
  position, not its name" rule at line 168 is the most-missed detail here; but a
  small annotated JSON-to-filename mapping diagram would serve better than a
  screenshot, and it is low priority. Everything on this page cross-checks
  against the spec at
  {doc}`/user_guide/python_package/c_data_preparation_api/09_output_format_specification_exports_directory`
  by reference rather than restatement, which is the right structure.

## user_guide/web_app/b_data_loading/08_troubleshooting_data_loading.md

- **explains** — A triage flow for load failures: identify the loading path,
  run a 2-minute checklist, then eight symptom → likely-causes → how-to-confirm →
  fix blocks, plus when to abandon a workflow and what to include in a help
  request.
- **has** — none.
- **verdict** — ADD
- **needs** — Troubleshooting pages usually earn NONE, but this one is the
  exception: **four of its eight symptoms are defined by an on-screen message
  the docs quote but never show.** A reader matching their screen against the
  docs needs the rendered form. All crops to the notification/error region only,
  ~600px display width, 2x, no cursor:
  1. `datasets.json not found` (line 215) — reproduce by connecting
     `#github-repo-url` to a repo path with no catalog.
     `_static/screenshots/data_loading/error-datasets-json-not-found.png`.
     (Shared with ``02_local_demo_tutorial`` item 2d — capture once, use twice.)
  2. "No embedding / no UMAP" (line 110) — reproduce with an H5AD lacking every
     `X_umap_*d` key. `_static/screenshots/data_loading/error-no-embedding.png`.
  3. "Dataset loads, but I see no points (blank canvas)" (line 131) — shared
     with ``a_orientation/02_system_requirements`` item 3.
  4. **The broken-export refusals.** `cellucid/tests/browser/fixtures/broken-exports/`
     ships nine deliberately invalid exports (`corrupt-payload`,
     `mismatched-manifest`, `missing-payload`, `truncated-catalog`,
     `wrong-format-catalog`, `truncated-payload`, and three `mixed-*` variants).
     These are already scripted for the test suite, so capturing their rendered
     error states is nearly free. Pick the three most likely in the wild —
     `missing-payload`, `mismatched-manifest`, `truncated-catalog` — and show
     each refusal. `_static/screenshots/data_loading/error-broken-export-<name>.png`.
     This directly serves "make the docs complete": a reader currently has no way
     to tell a corrupt export from a network failure.
- **notes** — The DevTools instructions recur here (line 113 "DevTools → Network:
  look for 404/CORS/timeout errors (don't guess)") with the same wet-lab
  accessibility problem flagged under ``02_system_requirements``. One annotated
  DevTools Network capture, reused from both pages, would fix it in both places.

## user_guide/web_app/b_data_loading/09_screenshots.md

- **explains** — Nothing procedural. It is a gallery of six CI acceptance
  captures produced by the synthetic-fixture browser test, grouped as "Loading
  and session controls", "Direct H5AD load", "Direct Zarr ZIP load", and
  "Browser-engine checks", plus a pointer to where the real server and Jupyter
  evidence lives.
- **has** — six figures, five of them used **nowhere else**:
  - line 13 `data_loading/data-loading-session-panel.png` `:width: 246px` (the
    10-page reused panel)
  - line 23 `data_loading/h5ad-current-loaded.png` `:width: 1280px` — only here
  - line 31 `data_loading/h5ad-current-visualization.png` `:width: 1280px` — only here
  - line 40 `data_loading/zarr-zip-loaded.png` `:width: 1280px` — only here
  - line 51 `data_loading/h5ad-firefox.png` `:width: 1280px` — only here
  - line 58 `data_loading/h5ad-webkit.png` `:width: 1280px` — only here
- **verdict** — REPLACE (dissolve the page; redistribute its contents)
- **needs** — **The gallery-page pattern should be retired, here and in the nine
  sibling sections that copy it** (`c/07`, `d/06`, `e/08`, `f/07`, `h/11`,
  `i/08`, `k/08`, `l/11`, `n/08`). Reasons, specific to this page:
  - **It strands evidence away from the instruction.** Five images that would
    each answer a "what does success look like?" question on
    {doc}`/user_guide/web_app/b_data_loading/03_browser_file_picker_tutorial` sit on a page nobody reaches while
    troubleshooting. Meanwhile `03_...` illustrates its H5AD and Zarr options
    with a reused Suo overview that shows neither.
  - **The dataset is a 120-cell, 6-gene synthetic fixture.** I opened
    `h5ad-current-loaded.png`: it is roughly a hundred blue/orange/purple dots
    scattered on a grid, with `#dataset-info` reading `CELLS 120 / GENES 6 /
    OBS FIELDS 2 / CONNECTIVITY None`. That is a correct acceptance artefact and
    a poor teaching image — a wet-lab reader learning what a loaded dataset
    looks like should see real single-cell structure, not a test fixture.
  - **Four of the six are near-duplicates.** `h5ad-current-loaded`,
    `h5ad-current-visualization`, `h5ad-firefox`, and `h5ad-webkit` are the same
    fixture in the same state; the last two differ only by browser engine, which
    is invisible in the crop. Three of the four earn nothing for a reader.
  - Redistribution plan: move `h5ad-current-loaded.png` to
    {doc}`/user_guide/web_app/b_data_loading/03_browser_file_picker_tutorial` Option #4; move `zarr-zip-loaded.png`
    to Option #5; **delete** `h5ad-current-visualization.png`,
    `h5ad-firefox.png`, `h5ad-webkit.png` from the docs (the browser matrix is a
    CI claim, and CI already asserts it — a reader cannot verify engine parity
    from a picture, and the repo's rule against carrying superseded material argues against keeping
    three images that duplicate one state); then delete this page and drop its
    card from `index.md:163`. Recapture the two survivors against `pancreas`
    rather than the synthetic fixture if the acceptance test can be pointed at a
    real export; if it cannot, keep the fixture but say so plainly in the
    caption.
- **notes** — Lines 67–75 already concede the pattern's weakness: "The synthetic
  picker captures above isolate browser-file behavior. Two data-backed workflows
  use the standard Pancreas dataset" — and then link out to the two pages that
  *do* embed their evidence inline. That is the argument for inline placement,
  written on the gallery page itself.

## user_guide/web_app/b_data_loading/10_standard_pancreas_dataset.md

- **explains** — What the `pancreas` sample is (3,696 cells, 3,753 genes, the
  scVelo vignette dataset), how to load it from the sample picker, which
  categorical/continuous fields, genes, edges, embeddings and the `velocity_umap`
  vector field it exposes, what is source-owned versus derived, and its pinned
  build contract.
- **has** — line 47 `vector_field_velocity/pancreas-velocity-3d.png`
  `:width: 1440px` (native 1440x1000). **Also used on**
  {doc}`/user_guide/web_app/i_vector_field_velocity/08_screenshots`.
- **verdict** — ADD
- **needs**
  1. **The load step has no picture** (lines 10–20), and it is four numbered
     clicks. Crop to `#session-section` with `#dataset-select` open and
     **Pancreatic endocrinogenesis (scVelo)** highlighted in the list, cursor on
     that option. `_static/screenshots/data_loading/pancreas-sample-select.png`.
  2. **The verification step has no picture** (lines 27–45). The page makes a
     precise, checkable claim: the compact Session statistics round to **4K**
     while **Coloring & Filtering** reports the exact **Showing all 3,696
     points**. That contrast between rounded and exact is a real source of reader
     confusion and is perfectly suited to one capture. Crop a vertical strip
     containing `#dataset-info` (CELLS 4K, GENES 4K) **and** the
     `#coloring-filtering-section` count line (Showing all 3,696 points) in the
     same frame. `_static/screenshots/data_loading/pancreas-counts-rounded-vs-exact.png`.
     No cursor.
  3. The existing velocity figure is appropriate and stays — it is the one
     figure in my sections that clearly illustrates the sentence above it.
- **notes** — Stale build string: line 51 caption says "The standard Pancreas
  sample in Build 2026-07-27.1"; the app ships `2026-07-27.24`. Cross-cutting
  note 3. Every numeric claim I could check against
  `cellucid-datasets` holds: `pancreas` is present in
  `cellucid-datasets/exports/datasets.json` with `suo` as `default`, matching
  line 7 ("The catalog id is exactly `pancreas`; Suo remains the catalog
  default").

## user_guide/web_app/b_data_loading/11_custom_dataset_repository.md

- **explains** — How to publish a public prepared-dataset catalog on GitHub: try
  the reference `theislab/cellucid-demo-custom-datasets/exports` repo, the exact
  `datasets.json` contract and its strict rules, building the dataset
  directories from Python or R, a four-step validation sequence, branch
  addressing, optional GitHub Pages, size discipline, and a fast-diagnosis table.
- **has** — **none.** This is the most operationally detailed page in the section
  and it is entirely unillustrated.
- **verdict** — SEQUENCE
- **needs** — The page opens with a one-minute hands-on ("Try the reference
  repository first", line 26) that is five clicks long and has no visual
  anchor. Capture that flow end to end; it doubles as the GitHub-loader
  documentation the whole section lacks. All at 1440x900, 2x:
  1. **The empty GitHub form.** Crop to `#github-repo-block` showing the
     placeholder `owner/repo/path/to/exports` and the **Load** button
     (`#github-connect-btn`). `_static/screenshots/data_loading/github-form-empty.png`,
     ~246px.
  2. **The form filled with the exact reference value.**
     `#github-repo-url` containing `theislab/cellucid-demo-custom-datasets/exports`
     (the exact string at line 36), cursor on `#github-connect-btn`.
     `_static/screenshots/data_loading/github-form-reference-repo.png`.
     (Same capture as ``02_local_demo_tutorial`` item 2b — take it once.)
  3. **The connected catalog.** `#dataset-select` expanded with all three
     entries — `synthetic-cell-types-2d`, `synthetic-development-3d`,
     `synthetic-trajectory-1d` — visible, matching the table at line 48.
     `_static/screenshots/data_loading/github-three-dataset-catalog.png`.
  4. **One loaded synthetic dataset per geometry**, to substantiate line 54
     ("Planar in 1D and 2D, Orbit in 3D") — a claim about navigation defaults
     that the reader is asked to take on faith. `synthetic-development-3d` in
     Orbit with its 3D `velocity_umap` overlay on is the most persuasive single
     frame. `_static/screenshots/data_loading/github-synthetic-development-3d.png`.
  5. **Terminal: the local catalog validation** from "2. Test the complete
     catalog through the local server" (line 226). Plain terminal running
     `cellucid serve ./exports --no-browser` with the printed
     `Viewer URL: …/?source=remote` visible. This is the second terminal capture
     the loader paths need and it belongs here, not only on the server page.
     `_static/screenshots/server/serve-exports-catalog-validation.png`.
  6. **The raw-catalog check** from step 3 (line 236). Browser tab showing
     `https://raw.githubusercontent.com/theislab/cellucid-demo-custom-datasets/main/exports/datasets.json`
     returning the JSON. The page warns the response "must return the JSON file
     itself, not an HTML sign-in page, a Git LFS pointer, or a 404" (line 251) —
     showing the good response makes the three bad ones recognisable by
     contrast. `_static/screenshots/data_loading/raw-catalog-json-response.png`.
- **notes** — The reference repository, its three dataset ids, and their exact
  contents are corroborated by the fixture copy at
  `cellucid/tests/browser/fixtures/demo-custom-exports/`, which contains
  `datasets.json` plus `synthetic-cell-types-2d/`, `synthetic-development-3d/`,
  and `synthetic-trajectory-1d/` — matching the table at line 48. The
  `synthetic-cell-types-2d` fixture has `connectivity/`, `obs/` (4 fields),
  `var/` (8 genes), and `points_2d.bin.gz`, consistent with "72 cells, 8 genes,
  four `obs` fields, a 2D embedding". Good — the fixtures make every capture in
  this list scriptable offline, without depending on GitHub availability at
  capture time.

## user_guide/web_app/b_data_loading/12_sample_dataset_provenance.md

- **explains** — Which of the five built-in samples can be rebuilt from a pinned
  recipe (only `pancreas`) and which four are recorded artefacts; exactly what
  the Pancreas build contract pins; and for the other four, what is recorded,
  what is not established, and what is not recoverable.
- **has** — none.
- **verdict** — ADD (minimal)
- **needs**
  1. **The sample picker with all five entries**, so the catalog table at line 18
     has a visual counterpart and a reader can match a row to what they see.
     Crop to `#dataset-select` expanded showing `suo`, `garcia`, `he`,
     `kanemaru`, `pancreas` with their cell counts. Cursor: no.
     `_static/screenshots/data_loading/sample-catalog-five-entries.png`.
     Verified against `cellucid-datasets/exports/datasets.json`, which lists
     exactly those five with `"default": "suo"`.
  2. Nothing else. The page's substance is a reproducibility argument about
     digests, environments, and what is *not* knowable; there is no screen that
     shows any of it, and inventing one would misrepresent the page.
- **notes** — This is the most carefully-hedged page in my four sections and I
  found no claim I could contradict. The "Not established" and "Not recoverable"
  lists (lines 98–125) are the kind of writing the rest of the docs should
  imitate. One small opportunity: line 34 says the four human samples "publish
  only the genes whose Ensembl accession resolved to a symbol, out of files that
  had already been reduced to 8,192 features upstream" and points at
  ``d_fields_coloring_legends/07_genes_in_the_built_in_samples``; a reader hunting
  a missing gene would be helped by a capture of the gene-search empty state
  (there is a browser test `tests/browser/gene-search-empty-state.spec.mjs` that
  produces it), but that image belongs on the section-D page, not here.

---

# G — `user_guide/web_app/g_cross_highlighting/`

## user_guide/web_app/g_cross_highlighting/ (EMPTY DIRECTORY)

- **path** — `user_guide/web_app/g_cross_highlighting/` — contains no files at
  all (`ls -la` shows only `.` and `..`).
- **explains** — Nothing. There is no page.
- **has** — n/a. There is also no `_static/screenshots/g_cross_highlighting/`
  directory.
- **verdict** — NONE — and specifically **do not commission any screenshot for
  this section, now or as part of this objective.** The feature it was written
  for does not exist in the app. Justification below.
- **needs** — Nothing to capture. There is no UI to point a camera at. The
  correct action is a two-line cleanup, not a documentation or screenshot task:
  1. `rmdir cellucid-python/docs/user_guide/web_app/g_cross_highlighting/` and
     `rmdir cellucid-python/docs/user_guide/web_app/r_screenshot_checklist/` —
     both are empty on-disk residue. Git cannot represent an empty directory, so
     they survived the commit that deleted their contents and were never cleaned
     up locally.
  2. Optionally remove `cellucid/assets/js/app/ui/modules/cross-highlight/`
     (contains only a zero-byte `.gitkeep`) and the orphan
     `// Cross-highlighting` comment at
     `cellucid/assets/js/app/analysis/ui/index.js:143`, which has no code under
     it. Both are reservations for unbuilt work and are exactly the kind of
     placeholder the repo's clean-break rule argues against keeping. This
     is a web-app change and is out of scope for a docs objective — flagging it,
     not proposing it here.
- **notes** —

  **What the directory used to hold, and why it is empty.** All seven files were
  deleted in one commit in `cellucid-python`:
  `617dbb5` ("Update", 2026-07-26), which removed
  `01_what_cross_highlighting_is_user_story.md`, `02_data_requirements.md`,
  `03_ux_design.md`, `04_performance_correctness_notes.md`,
  `05_troubleshooting_cross_highlighting.md`,
  `06_reference_implementation_notes.md`, and section index — plus, in the same
  commit, both files of `r_screenshot_checklist/`. The same commit also removed
  the `toctree` entries and the `{grid-item-card}` for cross-highlighting from
  {doc}`/user_guide/web_app/index`. This was a deliberate, complete de-registration,
  not an accident.

  **Why it was deleted: it documented a feature that was never built.** The
  deleted section index was titled "Cross-highlighting (analysis ↔ embedding;
  planned)" and opened with a `:::{warning}` saying the feature is "**planned**
  and **under development**" and that "some or all of the described UI may be
  missing or non-functional". The deleted `01_...user_story.md` defined it as:
  pick a subset from an analysis plot (a bar/bin/distribution) and see those same
  cells light up in the embedding viewer, with hover-preview, click-to-select,
  and a "Save as Page" action.

  **The feature does not exist in the source.** Confirmed directly, not inferred:
  - `cellucid/assets/js/app/ui/modules/cross-highlight/` contains exactly one
    file, a **zero-byte `.gitkeep`**. `ls -la` output: `total 0`.
  - `cellucid/assets/js/app/analysis/ui/index.js:143` is a bare
    `// Cross-highlighting` comment with **no export or code beneath it** — the
    next non-blank lines are the `DEFAULT EXPORT` banner comment.
  - Grep across `cellucid/assets/js` for `crossHighlight`, `cross_highlight`,
    `crossFilter`, `linkedSelection`, `propagateHighlight`, `sharedHighlight`
    returns **zero** hits. `syncHighlight` returns three hits, all
    `syncHighlightBufferForLod` in `rendering/viewer.js` and
    `rendering/highlight-renderer.js` — a GPU level-of-detail buffer refresh,
    unrelated to linking views. `brush` returns one hit, a comment about
    "Brushed brass/bronze bezel" styling in `orbit-anchor.js`.
  - The Analysis accordion (`#page-analysis-section`) and the highlight accordion
    (`#highlighted-cells-section`) are separate and unwired; there is no
    plot→embedding interaction to photograph.

  **It is not a Sphinx build hazard.** {doc}`/user_guide/web_app/index` lines
  151–170 use an **explicit** `toctree` list, not `:glob:`, and it runs
  ``f_highlighting_selection/index`` → ``h_analysis/index`` with no `g` entry (it
  also skips `m` and `r`). Sibling `:glob:` patterns are `[0-9]*` resolved
  relative to their own directory, so none can reach in. A directory with zero
  source files produces no document, so there is nothing for
  "document isn't included in any toctree" to fire on. The strict CI build
  (`cellucid-python/.github/workflows/docs-check.yml`,
  `sphinx-build -W --keep-going -b html docs docs/_build/html`) passes today.
  I also confirmed `docs/_static/screenshots/g_cross_highlighting/` does not
  exist and is referenced by nothing.

  **No dangling references anywhere.** A case-insensitive grep of the whole docs
  tree for `cross_highlight`, `cross-highlight`, `cross highlight`, `g_cross`,
  and `crosshighlight` produced exactly one hit, and it is unrelated:
  `user_guide/web_app/h_analysis/07_analysis_mode_gene_signature.md:24` —
  "which can be compared across highlight pages (groups)". A workspace-wide grep
  for `g_cross_highlighting|r_screenshot_checklist` returns zero hits. I
  separately confirmed that no `{doc}` or `{ref}` target on
  {doc}`/user_guide/web_app/q_troubleshooting_index/index` points into this directory.

  **The adjacent topic people will assume is missing is, in fact, covered.**
  "Do highlights propagate across views/panels?" is a different question from
  cross-highlighting, and it is documented well in
  {doc}`/user_guide/web_app/f_highlighting_selection/04_selection_synchronization`
  ("Selection synchronization (views, pages, filters, Python)"), whose
  "Sync between views and snapshots" section states that confirmed groups are
  **global** — highlighted across all panels — but still subject to each panel's
  own filter visibility, while an in-progress lasso gesture is tied to the panel
  it starts in. {doc}`/user_guide/web_app/f_highlighting_selection/01_highlight_mental_model` covers
  highlight groups, pages, and the "what you see is per view" rule.
  **That page is in section F, not mine** — but it is the page that would carry
  any screenshot a reader might come to `g_` looking for, and it is worth telling
  whoever audits section F that `web_app/multiview-two-panels.png` (1440x900,
  already in the tree) and `highlighting_selection/highlighting-selected-page.png`
  already cover adjacent ground for it.

  **What the app actually calls this area** (verified against
  `cellucid/index.html` lines 1540–1573, correcting one label): the accordion is
  `#highlighted-cells-section` and its visible `<summary>` reads **"Highlighting"**
  (the element id says "highlighted-cells", but the id is not the label — do not
  caption a screenshot "Highlighted Cells"). Inside it: `.highlight-mode-buttons`
  (`role="group"`, `aria-label="Highlight mode selection"`) holding four
  `.highlight-mode-btn` entries labelled **"Annotation based"**
  (`data-mode="annotation"`, the default, `aria-pressed="true"`), **"KNN drag"**,
  **"Proximity drag"**, and **"Lasso"**; `#highlight-mode-description` ("Alt+click
  a cell to highlight its group. Alt+drag to select a range.");
  `#highlight-pages-tabs` (`role="tablist"`, `aria-label="Highlight pages"`) with
  `#add-highlight-page` ("+"); `#highlight-count` (default text "No cells
  highlighted"); and `#clear-all-highlights` ("Clear"). Keyboard `X` clears
  highlights. The real model is **global highlight state rendered per view** —
  one `state.highlightArray` — not cross-highlighting.

  **Bottom line for this objective:** removing speculative docs for an unbuilt
  feature was the right call, and there is no content or screenshot gap to
  backfill until the feature ships. The only defect here is two empty
  directories on disk. Do not let the alphabetical gap in the section letters
  (`g`, `m`, `r` are all absent) be read as missing coverage.

---

# Q — `user_guide/web_app/q_troubleshooting_index/`

## user_guide/web_app/q_troubleshooting_index/index.md

- **explains** — The cross-section troubleshooting triage map: a 2-minute
  checklist, a symptom-picker table routing to seven anchored subsections
  (install/environment, data loading, rendering/GPU, selection/highlighting,
  analysis, export, community annotation), a "before you assume it's a bug" list
  of four common state mismatches, and a bug-report checklist.
- **has** — none. It is the only page in its directory.
- **verdict** — NONE
- **needs** — n/a. Justification below.
- **notes** —
  - **This page correctly needs no screenshots, and adding them would make it
    worse.** Its entire function is routing: every one of its ~45 links sends
    the reader to a deep-dive page that owns the symptom. A screenshot here
    would either duplicate an image from the destination page (adding a second
    place to keep in sync for zero navigational value) or illustrate one symptom
    among seven and imply a priority the triage table deliberately avoids.
  - **The four "before you assume it's a bug" mismatches (lines 51–63) are the
    exception worth noting, but the fix belongs elsewhere.** Filters-vs-
    visibility, active-view confusion, membership-vs-visibility, and dataset
    identity are all states that *are* worth showing — and each already has a
    destination page that should show it. Specifically, "Dataset identity
    mismatch" routes to ``b_data_loading/06_dataset_identity_why_it_matters``,
    where I have recommended exactly that capture. Keep the images at the
    destinations; keep this page textual.
  - **STALE-RISK if illustrated.** This page enumerates other sections'
    anchors. Any screenshot placed here would go stale whenever the destination
    section changes its UI, and — because the page is a hub — the staleness
    would be maximally visible. That is a second, independent reason to leave it
    figure-free.
  - Link integrity: I spot-checked the `{doc}` targets and the `{ref}` anchors
    (`install-env`, `data-loading`, `rendering`, `selection-highlighting`,
    `analysis`, `export`, `community-annotation`, `always-check`, `bug-report`)
    — every `{ref}` used in the symptom-picker table at lines 26–32 is defined
    later in the same file. **No target anywhere on this page points into
    `g_cross_highlighting/`**; selection/highlight problems route to
    ``f_highlighting_selection/06_troubleshooting_highlighting`` (line 168). So
    the empty G directory does not orphan any link from this page.
  - One prose accuracy point: line 74 says "Use a current stable Chrome, Edge,
    Firefox, or Safari release", matching the support table in
    `a_orientation/02_system_requirements.md:33`. Consistent.
