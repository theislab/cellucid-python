# Screenshot coverage audit — sections I, J, K

Scope: `cellucid-python/docs/user_guide/web_app/{i_vector_field_velocity, j_community_annotation, k_figure_export}`.
Every `.md` in all three sections was read in full. Docs were **not** modified; this file is the only output.

All paths below are relative to `cellucid-python/docs/` unless absolute.

---

## Summary

| Section | Pages | `{figure}` directives | **Unique** images | Pages with ≥1 figure |
|---|---:|---:|---:|---:|
| `user_guide/web_app/i_vector_field_velocity/` | 9 | 3 | 3 | **1 of 9** |
| `user_guide/web_app/j_community_annotation/` | 4 | 4 | **1** | 4 of 4 (all the same file) |
| `user_guide/web_app/k_figure_export/` | 10 | 11 | 6 | 6 of 10 |
| **Total** | **23** | **18** | **10** | **11 of 23** |

### Verdict counts

| Verdict | Count | Pages |
|---|---:|---|
| `OK` | 2 | `k/index.md`, `k/04` |
| `REPLACE` | 5 | `i/08`, `k/02`, `k/03`, `k/05`, `k/08` |
| `ADD` | 8 | ``i/index``, `i/01`, `i/03`, `i/04`, `i/05`, `i/06`, `k/01`, `k/06` |
| `SEQUENCE` | 7 | `i/02`, `i/07`, ``j/index``, `j/01`, `j/02`, `j/03`, `k/07` |
| `NONE` | 1 | `k/09` |
| `DIAGRAM` | 0 | — (see note under `i/01`) |
| `STALE-RISK` | (overlay on 3) | `k/05` + `k/08` (verified stale), `k/03` (verified stale) |

`OK` means *content is correct*; it does **not** exempt the image from the global 2× re-shoot below.

### Three global defects (apply to all 10 images, tracked once, not repeated per page)

1. **Every image is 1× device-pixel-ratio.** `overlay-controls.png` is 1440×1000 shown at `:width: 1440px`; `download.png` is 224×246 at `:width: 224px`; `panel-overview.png` is 246×256 at `:width: 246px`. On any retina/HiDPI display and in print these are soft. The user's stated objective is "high resolution screenshots … scaled down depending on the size". **Fix: capture at `deviceScaleFactor: 2` and keep the existing `:width:` values**, so a 224 px-wide panel ships as a 448 px PNG rendered at 224 px. This is a mechanical re-shoot of all 10 files and should be done in the same pass as any content change.
2. **No mouse pointer anywhere in any of the 18 figure references.** The user explicitly asked for a visible cursor "whenever relevant". Every "click X" instruction in these three sections is currently unillustrated.
3. **Sizes are inconsistent for the same class of thing.** Sidebar-section crops are `224px` in `k/`, `246px` in `j/` and ``k/index``. Pick one (`246px` for a whole accordion section including its header chevron, `224px` for a sub-block inside it) and apply it. `community_annotation/disconnected-panel.png` is 246×123 — a 123-pixel-tall sliver.

### Two verified staleness hits (not guesses — checked against source)

- **`figure_export/labels-annotations.png` is stale.** `cellucid/assets/js/app/ui/modules/figure-export/figure-export-ui.js:541` creates `#figure-export-reference-grid` with label `' Reference grid'`, inserted in `annotationChecks` between *Legend* and *3D orientation*. Lines 580–591 add a `chromeDisclosure` hint reading `` `Never exported: ${NON_EXPORTED_VIEWER_CHROME.join(', ')}.` ``. **Neither appears in the screenshot**, which shows only `AXES / LEGEND / 3D ORIENTATION / DEPTH SORT`. Meanwhile `k/02` documents "**Reference grid** is enabled by default", `k/06` has a whole *Reference grid* subsection, and `k/01` says the chrome list "is disclosed in the Annotations block". The prose is current; the picture predates it.
- **`figure_export/download.png` is stale.** The screenshot's hint text ends at "…Hybrid embeds the shader-rendered point pass." The current source (`figure-export-ui.js:761–765`) continues: "Vector points carry the viewer's atmospheric fog but are flat discs: only Hybrid and PNG reproduce the 3D sphere lighting." That second sentence is exactly the claim `k/06` §*Atmospheric fog and 3D lighting in vector exports* makes in prose.

### Two empty documentation directories found in passing

`user_guide/web_app/r_screenshot_checklist/` and `user_guide/web_app/g_cross_highlighting/` contain **zero files**. `r_screenshot_checklist` is directly relevant here — the checklist that would govern this work does not exist.

---

## Capture mechanics that the Playwright tool needs (shared by all three sections)

Verified from the existing browser specs, so the screenshot tool can reuse a proven path.

**Real datasets.** `cellucid-datasets/exports/datasets.json` — `pancreas` is the **only bundled sample carrying `vector_fields`** (`suo`, `garcia`, `he`, `kanemaru` have none, so the entire velocity block is `display:none` on them). Pancreas = *Pancreatic endocrinogenesis (scVelo)*, 3696 cells, `default_field: "velocity_umap"`, dims `[1,2,3]`, `default_dimension: 3`. **The "3,696 cells / independent 1D, 2D, 3D" claim in `i/08` is verified correct.**

**Deterministic fixture.** `/?exportsBaseUrl=${ENCODED_EXPORTS_BASE_URL}&dataset=current-ui-prepared` → "Current UI prepared fixture", 120 points. `ENCODED_EXPORTS_BASE_URL` is from `cellucid/tests/browser/helpers/origins.mjs`; dismiss the welcome modal with `dismissWelcome(page)` from `cellucid/tests/browser/helpers/welcome.mjs` (asserts `#welcome-modal`, presses Escape). **This fixture is right for contract tests and wrong for documentation** — see `k/02`.

**Velocity is a Figure-Export blocker.** `figure-export/utils/overlay-fidelity.js` (`VELOCITY_OVERLAY_CONTROL_ID = 'velocity-overlay-enabled'`) blocks export while the overlay is on. **All velocity imagery must be canvas/viewport capture of `#glcanvas`, never Figure Export output.**

**Animation determinism.** The overlay is a 5-pass GPU particle simulation with curl-noise turbulence and a persistent trail framebuffer. Two captures of the *same* settings differ. For any before/after pair the tool must: same dataset, same camera (`Reset camera`), toggle the overlay off→on to clear trail buffers, change exactly one slider, then settle a **fixed** number of frames before capturing both halves. Without that, a "density 5K vs 200K" pair also differs in particle placement and the reader cannot attribute the difference to the parameter.

---

# Section I — `user_guide/web_app/i_vector_field_velocity/`

## Verified control inventory (drives the per-parameter guidance below)

Source of truth is `cellucid/index.html` lines 854–1118 (markup) + `cellucid/assets/js/app/ui/modules/velocity-overlay-controls.js`.

Container structure — crop targets:

| Crop target | Selector | Note |
|---|---|---|
| Whole feature block | `#velocity-overlay-controls` | `display:none` when the dataset has no vector fields |
| Header row + info button | `#velocity-overlay-controls > div.d-flex.items-center.gap-1` | holds `Vector Field Overlay:` and `#vector-field-help-btn` |
| Info popover | `#vector-field-help-tooltip` | `role=dialog`, `hidden` by default |
| Settings box | `#velocity-overlay-settings` | `display:none` until `#velocity-overlay-enabled` is checked |
| Advanced disclosure | `#velocity-overlay-settings > details` (summary text `Advanced Visual Settings`) | **collapsed by default**, no `open` attribute |
| Advanced body | `#velocity-overlay-settings > details > div.mt-1` | |
| Status line | `#velocity-overlay-info` | |
| Canvas | `#glcanvas` | |

**There is no wrapper element per advanced sub-group.** `Particle Rendering`, `Trail Settings`, `HDR & Bloom`, `Color Grading`, `Cinematic Effects` are bare `div.legend-help.font-medium` headers followed by sibling `.control-block`s. Per-group crops must be **geometric** (from the header div's top to the last slider before the next header), not selector-based.

Controls (label verbatim / selector / min–max–step / default). All 24 sliders confirmed present and all documented defaults confirmed against the markup:

| Group | Label | Input id | min | max | step | default | Badge |
|---|---|---|---|---|---|---|---|
| core | `Vector field:` | `#velocity-field` | — | — | — | dataset default | — |
| core | `Particle density:` | `#velocity-density` | 1 | 500 | 1 | 15 | `15K` |
| core | `Flow speed:` | `#velocity-speed` | 5 | 500 | 1 | 300 | `3.0×` |
| core | `Trail length:` | `#velocity-lifetime` | 10 | 1500 | 10 | 800 | `8.0s` |
| core | `Particle size:` | `#velocity-size` | 0.5 | 30 | 0.5 | 1 | `1` |
| core | `Opacity:` | `#velocity-opacity` | 0 | 100 | 1 | 60 | `60%` |
| core | `Color scheme:` | `#velocity-colormap` | — | — | — | `viridis` | — |
| core | `Sync with LOD` | `#velocity-sync-lod` | — | — | — | **checked** | — |
| Particle Rendering | `Intensity:` | `#velocity-intensity` | 0.05 | 1.5 | 0.05 | 0.25 | |
| Particle Rendering | `Glow amount:` | `#velocity-glow` | 0 | 1 | 0.05 | 0.3 | |
| Particle Rendering | `Comet stretch:` | `#velocity-comet-stretch` | 0 | 2 | 0.05 | 0.6 | |
| Particle Rendering | `Core sharpness:` | `#velocity-core-sharpness` | 0 | 1 | 0.05 | 0.7 | |
| Trail Settings | `Trail persistence:` | **`#velocity-trail-fade`** | 0.9 | 0.995 | 0.005 | 0.925 | |
| Trail Settings | `Chromatic fade:` | `#velocity-chromatic-fade` | 0 | 1 | 0.05 | 0 | |
| Trail Settings | `Turbulence:` | `#velocity-turbulence` | 0 | 1 | 0.05 | 0.3 | |
| HDR & Bloom | `Exposure:` | `#velocity-exposure` | 0.1 | 2 | 0.05 | 0.5 | |
| HDR & Bloom | `Bloom strength:` | `#velocity-bloom-strength` | 0 | 0.5 | 0.01 | 0.08 | |
| HDR & Bloom | `Bloom threshold:` | `#velocity-bloom-threshold` | 0.1 | 1 | 0.05 | 0.75 | |
| HDR & Bloom | `Anamorphic ratio:` | `#velocity-anamorphic` | 1 | 3 | 0.1 | 1.2 | |
| Color Grading | `Saturation:` | `#velocity-saturation` | 0.5 | 2 | 0.05 | 1.15 | |
| Color Grading | `Contrast:` | `#velocity-contrast` | 0.5 | 2 | 0.05 | 1.05 | |
| Color Grading | `Highlights:` | `#velocity-highlights` | 0.5 | 1.5 | 0.05 | 0.85 | |
| Color Grading | `Shadows:` | `#velocity-shadows` | 0.5 | 1.5 | 0.05 | 1.05 | |
| Cinematic Effects | `Vignette:` | `#velocity-vignette` | 0 | 1 | 0.05 | 0 | |
| Cinematic Effects | `Film grain:` | `#velocity-film-grain` | 0 | 0.1 | 0.002 | 0 | |
| Cinematic Effects | `Chromatic aberr.:` | `#velocity-chromatic-aberration` | 0 | 1 | 0.05 | 0 | |

**Gotcha for the tool:** the label reads `Trail persistence:` but the id is `#velocity-trail-fade`. The JS internal validation label is also `Trail fade`. Docs and screenshot both correctly say *Trail persistence* — no staleness, but a scripting trap.

### Honest triage: still vs pair vs recording

The overlay is animated, so ask per parameter *"does a frozen frame carry the information?"*

- **A single still is enough** (the parameter changes what a frame looks like): `Particle density`, `Trail length`, `Particle size`, `Opacity`, `Color scheme`, `Intensity`, `Glow amount`, `Comet stretch`, `Core sharpness`, `Trail persistence`, `Chromatic fade`, `Exposure`, `Bloom strength`, `Bloom threshold`, `Anamorphic ratio`, `Saturation`, `Contrast`, `Highlights`, `Shadows`, `Vignette`, `Film grain`, `Chromatic aberr.` → each of these wants a **low/high PAIR**.
- **A still cannot carry it** — `Flow speed`. It is a rate. At a fixed lifetime, changing 0.05× → 5.0× changes how far a particle travels per second, which a frozen frame renders as roughly similar streaks. **Do not fake a `Flow speed` before/after still.** Either ship a short recorded sequence (2×3 s loops side by side) or state in prose that the effect is temporal and show only the slider. `i/03`'s claim "Reduce `Flow speed:` before increasing trails" is a motion claim and cannot be evidenced by a PNG.
- **A still cannot carry it either** — `Turbulence` is marginal: at 0 vs 1 a frozen frame *does* differ (coherent streamlines vs wandering paths), so a pair is legitimate, but it is much weaker than the motion difference. Ship the pair, caption it honestly.
- **Needs a paired capture plus a number** — `Sync with LOD`. Its visible effect is *fewer particles when zoomed out* (LOD factor 0/0.25/0.5/1.0). A pair must be zoomed **out** to a coarse LOD with the box checked vs unchecked, and is only convincing next to an FPS readout.

---

## {doc}`/user_guide/web_app/i_vector_field_velocity/index`

**explains** — Section landing page: what the GPU particle overlay is, the Pancreas starting question, and a card grid linking the eight sub-pages.

**has** — `none`.

**verdict** — **ADD** (1 hero image).

**needs**
- Dataset: `pancreas` (the only bundled sample with velocity).
- Interaction: load `?dataset=pancreas` → `dismissWelcome` → open `Visualization` accordion → `Compare Views → Dimension:` = `2D`, `Mode:` = `Planar` → check `#velocity-overlay-enabled` → wait for `#velocity-overlay-info` to read `Available for 1D, 2D, 3D.` and for the `Velocity overlay ready` toast → settle 120 frames.
- Crop: full app viewport, 1440×900 logical, but **crop the right edge to where the data ends** — do not ship 300 px of empty grid (the current section images all do).
- Cursor: no.
- Filename: `vector_field_velocity/pancreas-2d-overview.png` at `:width: 1440px`.

**notes** — The "Pages in this section" card grid has **no card for ``06_edge_cases``**, although the "Recommended reading order" list above it does reference it. Every other numbered page has a card.

---

## {doc}`/user_guide/web_app/i_vector_field_velocity/01_what_vector_fields_are_user_facing`

**explains** — The mental model: one displacement vector per cell in the same coordinate space as the embedding, dimension-specific, rendered as particles seeded from visible cells.

**has** — `none`.

**verdict** — **ADD** (2 images). *Considered DIAGRAM and rejected*: the two ideas that would want a drawing (row-order correspondence; "vectors are dimension-specific") are both better served by real UI captures, because the failure the page warns about is invisible in a diagram and visible in the app.

**needs**
1. *"What good looks like vs what noise looks like"* — same Pancreas 2D view, two captures: (a) `Turbulence:` 0.00, `Bloom strength:` 0.00, `Trail length:` 3.0s — coherent flow along the endocrine continuum; (b) `Turbulence:` 1.00 — the same field rendered as procedural wander. Crop `#glcanvas` only, square, ~800×800 logical each. No cursor. Filenames `vector_field_velocity/flow-coherent.png` / `flow-turbulent.png`, `:width: 400px` each so they sit side by side. **This is the single highest-value velocity image in the section**, because §"The flow looks random" and `i/04`'s repeated "turbulence can suggest structure that isn't in the data" warning are currently pure assertion.
2. *The availability hint* — see `i/02` step 2 below; reference the same file rather than duplicating.

**notes** — §"How to tell whether your dataset contains vector fields" says "If the dataset has **no vector fields**, the **Vector Field Overlay** block will not appear." **Verified correct** — `syncAvailability()` sets `#velocity-overlay-controls.style.display = 'none'`; it is hidden, not greyed. The three "Implementation references" at the foot of the page all resolve to real files.

---

## {doc}`/user_guide/web_app/i_vector_field_velocity/02_enabling_overlay_and_selecting_field`

**explains** — Where the toggle lives, the three deliberate availability states, what dimension switching does, and dataset/session reset behaviour.

**has** — `none`. (Points at ``08_screenshots`` in prose only.)

**verdict** — **SEQUENCE** (4 images). This page is a click-by-click enabling procedure with three named UI states and zero pictures; it is the page a wet-lab reader hits first and gives up on.

**needs** — one ordered strip, all on `pancreas`, all cropped to `#velocity-overlay-controls`, `:width: 246px`:

1. **Block present, overlay off.** State before: dataset loaded, `Visualization` open. Action: none. State after: `#velocity-overlay-enabled` unchecked, `#velocity-overlay-settings` hidden. Cursor: **yes, hovering `#velocity-overlay-enabled`** — this is the "step 6: Check Show overlay" instruction. → `vector_field_velocity/enable-01-off.png`
2. **Disabled for the current dimension.** How to reach it without a bespoke dataset: use the `_test` export (`cellucid-datasets/exports/_test/`, field `T_fwd_umap`) or, cleanest, reuse the spec trick from `cellucid/tests/browser/velocity-multiview-lifecycle.spec.mjs` — `page.route()` rewrites `dataset_identity.json` to advertise `available_dimensions: [2]` only, then set the view to 3D. State after: checkbox **disabled**, `#velocity-overlay-info` reads `Vector fields available for 2D. Switch embedding dimension to enable.` Crop `#velocity-overlay-controls`. Cursor: no. → `vector_field_velocity/enable-02-wrong-dimension.png`
3. **Loading.** Route-delay `vectors/0_2d.bin.gz` by ~2 s, check the box, capture while `#velocity-overlay-info` reads `Loading vector field…`. Cursor: no. → `vector_field_velocity/enable-03-loading.png`
4. **Enabled.** `#velocity-overlay-settings` visible, `#velocity-overlay-info` = `Available for 1D, 2D, 3D.`, `Vector field:` = `velocity_umap`. Cursor: no. → `vector_field_velocity/enable-04-ready.png`

Plus one small extra: the **info popover**. The page says "A small `i` button beside the heading explains that the animation uses dimension-specific per-cell velocity or drift vectors" — that popover has never been shown. Un-`hidden` `#vector-field-help-tooltip`, crop the header row + tooltip, cursor on `#vector-field-help-btn`. → `vector_field_velocity/overlay-info-popover.png`, `:width: 246px`.

**notes** — The four literal status strings quoted in this page (`Vector fields available for …`, `Available for …. Select a vector field …`, `Loading vector field…`, `Failed to load vector field.`) all match the source exactly, including the U+2026 ellipsis. No staleness.

---

## {doc}`/user_guide/web_app/i_vector_field_velocity/03_core_parameters_document_exact_ui_labels`

**explains** — The 8 core controls with exact labels, ranges and initial values, plus a "scientific reading" and a "performance-first" preset.

**has** — `none`.

**verdict** — **ADD** (1 panel crop + 5 parameter pairs). Per-parameter, not one blanket entry:

**needs**

*Panel crop* — `#velocity-overlay-settings` **with the `details` collapsed**, so the reader sees exactly the 8 rows this page's table describes and nothing else. The current `overlay-controls.png` is a whole-app shot at 1440 px; this page needs the sidebar block alone. Cursor: no. → `vector_field_velocity/core-controls.png`, `:width: 224px`.

*Parameter pairs* — all Pancreas 2D Planar, `#glcanvas` cropped to the data bounding box, identical camera, overlay toggled off→on between captures, 120 settle frames, everything else at default. Each pair `:width: 400px` side by side.

| Parameter | Low capture | High capture | Filename stem | Cursor |
|---|---|---|---|---|
| `Particle density:` | `#velocity-density` = 5 (`5K`) | = 200 (`200K`) | `param-density-{5k,200k}` | on the slider thumb in the **low** image only |
| `Trail length:` | `#velocity-lifetime` = 20 (`0.2s`) | = 1500 (`15.0s`) | `param-trail-length-{short,long}` | no |
| `Particle size:` | `#velocity-size` = 0.5 | = 8 | `param-size-{0p5,8}` | no |
| `Opacity:` | `#velocity-opacity` = 20 | = 100 | `param-opacity-{20,100}` | no |
| `Color scheme:` | `#velocity-colormap` = `viridis` | = `coolwarm` | `param-colormap-{viridis,coolwarm}` | no |

*Do not produce* a `Flow speed:` pair — see the triage note above. Instead crop the `Flow speed:` `.control-block` alone (slider + `3.0×` badge) so the label and range are documented, and let the prose own the temporal claim.

*`Sync with LOD`* — one pair, both zoomed out far enough to hit a coarse LOD, checked vs unchecked, with the performance HUD visible if one can be enabled. Filename `param-sync-lod-{on,off}`. Mark clearly in the caption that the difference is particle count at coarse LOD, not at full zoom.

**notes** — Every range and initial value in this page's table matches the markup exactly (verified row by row against `index.html` 892–947). §"Dense haze" references `Trail persistence:` — an advanced control — from a core-controls page; that is fine but means the reader needs `04` open too.

---

## {doc}`/user_guide/web_app/i_vector_field_velocity/04_advanced_parameters_document_every_setting`

**explains** — All 18 advanced settings across five groups, with ranges, defaults, performance impact and a "don't mislead" warning.

**has** — `none`.

**verdict** — **ADD** (5 group crops + 8 parameter pairs). This page documents 18 controls and shows zero.

**needs**

*Group crops* — five geometric crops of the expanded `#velocity-overlay-settings > details > div.mt-1`, one per header. There is **no per-group container element**, so the tool must crop from each `div.legend-help.font-medium` header's top to the bottom of the last `.control-block` before the next header. Cursor: no. `:width: 224px`.
`vector_field_velocity/adv-particle-rendering.png`, `adv-trail-settings.png`, `adv-hdr-bloom.png`, `adv-color-grading.png`, `adv-cinematic-effects.png`.

*Parameter pairs* — same capture discipline as `03`. Prioritised, because 18 pairs is 36 images:

| Priority | Parameter | Low | High | Why it earns an image |
|---|---|---|---|---|
| **must** | `Trail persistence:` | `#velocity-trail-fade` = 0.900 | = 0.995 | The page calls it "the 'looks like fog' knob" and spends a whole subsection distinguishing it from `Trail length:`. Those two are the most confused pair in the section and **the two must be shown together in one 2×2 contact sheet** — `Trail length` short/long × `Trail persistence` low/high. Filename `adv-trail-length-vs-persistence.png`, `:width: 700px`. |
| **must** | `Turbulence:` | 0.00 | 1.00 | The scientific-honesty claim of the whole section. May reuse the `i/01` pair. |
| **must** | `Bloom strength:` | 0.00 | 0.50 | Named as the first thing to zero for both performance and honesty, in four separate pages. |
| **should** | `Comet stretch:` | 0.00 | 2.00 | The only control that makes *direction* readable in a still. |
| **should** | `Intensity:` | 0.05 | 1.50 | First remedy for "too faint" in `05` and `07`. |
| **should** | `Glow amount:` | 0.00 | 1.00 | |
| **nice** | `Anamorphic ratio:` | 1.0 | 3.0 | Visually dramatic, cheap to capture. |
| **nice** | Colour grading | 2×2 sheet: `Saturation` 0.5/2.0 × `Contrast` 0.5/2.0 | | Four subtle controls do not each deserve a pair; one contact sheet does. |

*Skip entirely* (state so rather than shipping weak images): `Core sharpness`, `Chromatic fade`, `Exposure`, `Bloom threshold`, `Highlights`, `Shadows`, `Vignette`, `Film grain`, `Chromatic aberr.` — the group crop shows the control, the table gives the range, and a low/high pair at documentation scale would be near-indistinguishable. Document that decision in the section's screenshot page instead of silently omitting.

**notes** — All 18 ranges/defaults match the markup. One engine detail the page could state and does not: particle colour is `texture(u_colormapTex, vec2(v_speedNorm, 0.5))` — colour is **only** normalized speed, there is no direction encoding and no alternative colour mode. The page says this correctly in `03` ("it does not encode direction") but a reader on `04` looking at *Color Grading* might infer otherwise.

---

## {doc}`/user_guide/web_app/i_vector_field_velocity/05_performance_and_quality`

**explains** — What actually scales (particles vs full-screen trail/bloom passes), three named presets, a triage workflow, and "it looks wrong but is working" pitfalls.

**has** — `none`.

**verdict** — **ADD** (1 triptych + 2 pairs).

**needs**
1. **Preset triptych.** Three captures of the same Pancreas 2D view at *Laptop-safe* (`density 5K, trail 2.0s, opacity 50%, bloom 0.00, persistence 0.900, turbulence 0.00`), *Desktop balanced* (defaults), *Presentation/cinematic* (`density 200K, trail 15.0s, comet stretch 1.0, persistence 0.99, bloom 0.30`). Crop `#glcanvas` to the data bounds. No cursor. → `vector_field_velocity/preset-{laptop,balanced,cinematic}.png`, `:width: 460px` each. This is the page's whole argument and it is currently three bullet lists.
2. **"Too faint" / "too bright" pair.** §"Quality pitfalls" names both failure looks. `opacity 10% + size 0.5 + intensity 0.05` vs `density 300K + trail 15s + persistence 0.99 + bloom 0.50 + opacity 100%`. → `vector_field_velocity/pitfall-{too-faint,too-bright}.png`, `:width: 400px`.
3. **`Sync with LOD` pair** — may be shared with `03`.

**notes** — §4 tells the reader to enable **Visualization → Renderer settings → Level-of-Detail (LOD)**. That control lives outside `#velocity-overlay-controls` and is **not shown anywhere in this section**; a `Renderer settings` crop with the LOD checkbox and a cursor on it would close the loop. Treat as an optional 4th image, `:width: 224px`.

---

## {doc}`/user_guide/web_app/i_vector_field_velocity/06_edge_cases`

**explains** — 13 named edge cases across data / scale / UI-state / environment, each as *what you see → why → confirm → do*.

**has** — `none`.

**verdict** — **ADD**, but **selective**. Most entries are correctly text-only; three are visual and currently unillustrated.

**needs**
1. *Edge case UI/state 3 — "Render mode = Volumetric smoke cloud"*: "**What you see:** you can't find the overlay controls." Capture `Visualization` with `Render mode:` = `Volumetric smoke cloud (alpha)`, showing that `#velocity-overlay-controls` is genuinely absent, next to the `Points` state where it is present. Cursor on the `Render mode:` select in the smoke capture. → `vector_field_velocity/render-mode-{smoke,points}.png`, `:width: 224px` each.
2. *Edge case UI/state 4 — "All cells filtered out"*: reachable exactly as `velocity-multiview-lifecycle.spec.mjs` does it — `#legend` → `Hide All` button, assert `#filter-count` = `Showing 0 of N points`. Capture the empty canvas plus the `Showing 0 of 3,696 points` counter in one crop. → `vector_field_velocity/all-cells-filtered.png`, `:width: 900px`.
3. *Data edge case 4 — "All-zero or extremely small vectors"*: route-inject a zero-filled vector binary (the spec's `createVelocityBytes` pattern with `vx = vy = 0`) and capture the near-static overlay. This is the one data edge case with a distinct visual signature. → `vector_field_velocity/zero-magnitude-field.png`.

**NONE (justified)** for edge cases 3 (wrong basis), 5 (NaN/Inf), 6 (huge magnitudes), 7 (row-order mismatch) and all three *Environment* cases: these are diagnosed in Python or the console, and the failure mode of 7 is explicitly "renders but is misleading" — a screenshot of a misleading render, presented as documentation, is worse than no screenshot.

**notes** — Edge case 3 says "the current detection rules are UMAP-oriented (e.g. `velocity_umap_2d`)". Confirmed against `cellucid-python/src/cellucid/vector_fields.py` naming `<field>_umap_<dim>d`. No staleness found on this page.

---

## {doc}`/user_guide/web_app/i_vector_field_velocity/07_troubleshooting_velocity_overlay`

**explains** — Symptom → cause → confirm → fix → prevention for nine overlay failures, plus how to add vector fields in Python and a table of three literal error strings.

**has** — `none`.

**verdict** — **SEQUENCE** (3 error-state captures). This is the longest page in the section (13 kB) and the one a stuck user reaches; every symptom is described in words.

**needs** — three `#velocity-overlay-controls` crops, `:width: 246px`, no cursor, each reached by route interception:

1. *"Overlay toggle instantly turns off"* — `page.route()` the vector binary to a 500. Capture with `#velocity-overlay-info` = `Failed to load vector field.` **and** the error toast (`.notification-error`) visible in the same frame — the page says "You may also see a toast/notification containing an error message" and that pairing is exactly what makes the picture worth having. Crop must include both the sidebar block and the toast, so this one is a viewport crop, `:width: 900px`. → `vector_field_velocity/error-failed-to-load.png`
2. *"No fields appear in dropdown"* — the `h5ad-ui-smoke.spec.mjs` path: after switching to a no-vector sample, `#velocity-field` has value `''` and **0 options**. Capture the open, empty select. → `vector_field_velocity/error-empty-field-select.png`
3. *"Velocity rendering unavailable."* — `#velocity-overlay-info` in that state (the spec drives it by failing the transform-feedback allocation). → `vector_field_velocity/error-rendering-unavailable.png`

**notes**
- The three "Common error messages" quoted at the bottom (`No 2D/3D vector field found for "…" in obsm`, `Vector field "…" has N dimensions, expected D`, `vectors length X !== expected Y`) were **not** re-verified against source in this pass — flag for the doc-accuracy sweep, not for screenshots.
- The Python snippets reference `cellucid.add_transition_drift_to_obsm(...)` and `prepare(...)`. The user's objective mentions *"Sometimes you may need jupyter or terminal screenshots when showing the data loader paths."* This page is the strongest candidate in section I: **a Jupyter capture of `adata.obsm.keys()` showing `velocity_umap_2d` next to `X_umap_2d`** would make §"How do I add vector fields to my dataset?" concrete. Suggested: `vector_field_velocity/jupyter-obsm-keys.png`, notebook cell + output only, `:width: 700px`. Note that `cellucid.prepare(...)` is used in the snippets while the agent-verified Python entry point is `prepare_data(...)` — worth confirming which alias the package exports before capturing a terminal/notebook image that would bake the wrong name in.

---

## {doc}`/user_guide/web_app/i_vector_field_velocity/08_screenshots`

**explains** — The section's "verified captures" page: three full-app images (2D core controls, 2D advanced expanded, 3D orbit).

**has** — 3 figures, all unique to this page:
- L10 `vector_field_velocity/overlay-controls.png`, `:width: 1440px` (file is 1440×1000)
- L21 `vector_field_velocity/advanced-settings.png`, `:width: 1440px` (1440×1000)
- L32 `vector_field_velocity/pancreas-velocity-3d.png`, `:width: 1440px` (1440×1000)

No reuse elsewhere in the docs.

**verdict** — **REPLACE** (all three).

The *content* is good and, unusually, fully verifiable: I read all three at full size and every documented default is legible and correct (15K / 3.0× / 8.0s / 1.0 / 60% / Viridis / Sync with LOD on; and in the advanced shot 0.25 / 0.30 / 0.60 / 0.70 / 0.925 / 0.00 / 0.30 / 0.50 / 0.08 / 0.75 / 1.2 / 1.15 / 1.05 / 0.85 / 1.05 / 0.00 / 0.000 / 0.00). The captions are accurate. **The framing is the problem:**

1. **Clipped at the top.** `overlay-controls.png` begins with `VECTOR FIELD OVERLAY:` flush against row 0 — the label is shaved. `advanced-settings.png` begins mid-`ADVANCED VISUAL SETTINGS`. The sidebar was mid-scroll when the shot was taken.
2. **Stray chrome.** All three carry a floating round **✕ button** at roughly (309, 27) over the canvas, unexplained and unmentioned in any caption.
3. **Bleed.** A ~25 px vertical sliver of a second panel column is visible at x≈0–25 in all three.
4. **Wasted frame.** In the 2D captures the data ends around x≈1150 of 1440; a quarter of every image is blank grid.
5. **1×.**

**needs** — three re-shoots, `pancreas`, `dismissWelcome`, sidebar scrolled so the target block starts ~16 px below the crop top, the ✕ dismissed, `deviceScaleFactor: 2`:
- `overlay-controls.png` — 2D Planar, advanced **collapsed**, crop from the top of `#velocity-overlay-controls` through the right edge of the data. Cursor: no. Keep `:width: 1440px`.
- `advanced-settings.png` — same view, `details` open, sidebar scrolled so all five group headers are in frame. Keep `:width: 1440px`.
- `pancreas-velocity-3d.png` — 3D Orbit. Note the orbit-anchor dot is visible mid-canvas; either turn `Show orbit anchor` off or mention it. Keep `:width: 1440px`.

**notes**
- The page's provenance claim — "3,696 cells with independent 1D, 2D, and 3D embeddings and dimension-matched `velocity_umap` vectors" — is **verified** against `cellucid-datasets/exports/pancreas/dataset_identity.json` (`available_dimensions [1,2,3]`, `default_dimension 3`, three `vectors/0_{1,2,3}d.bin.gz`).
- The 3D caption says the view is "colored by `clusters`". `CATEGORICAL OBS: clusters` is indeed visible in the crop, but the particle layer at 60 % opacity almost completely occludes the point cloud, so the *claim* is not visible in the *image*. Either lower `Opacity:` to ~35 % for that capture or drop the colour claim from the caption.
- This page is where the "skipped parameter pairs" decision from `i/04` should be recorded, and where the honest statement about `Flow speed:` being unrepresentable in a still belongs.

---

# Section J — `user_guide/web_app/j_community_annotation/`

## {doc}`/user_guide/web_app/j_community_annotation/index`

**explains** — Feature overview: offline-first, scope-based (dataset + repo + branch + user), role split, consensus maths, local vs GitHub storage, and the annotator/author quickstart tabs.

**has** — 1 figure, L40 `community_annotation/disconnected-panel.png`, `:width: 246px`. **Reused on all four pages.**

**verdict** — **SEQUENCE** (see the lifecycle section below).

**needs** — steps 1, 6, 9, 11 of the lifecycle: the disconnected panel it already has, the **connected status panel**, a **Pending** voting modal and a **Consensus** voting modal. The consensus-maths section (`voters` / `netVotes` / `confidence`) is the page's core content and is completely unillustrated even though the modal renders the numbers literally.

**notes** — §"Shareable Links" documents `?annotations=owner/repo` and `?annotations=owner/repo@branch`; the **Copy share link** icon action that produces them is documented in `03` and never shown.

---

## {doc}`/user_guide/web_app/j_community_annotation/01_annotator_guide`

**explains** — The annotator's whole round: connect, Pull, find the 🗳️ column, open the voting modal, read the consensus line, vote, comment, propose a suggestion, Publish, wrap up; then a large troubleshooting catalogue.

**has** — 1 figure, L82 `community_annotation/disconnected-panel.png`, `:width: 246px`. Same file as the other three pages.

**verdict** — **SEQUENCE**.

**needs** — lifecycle steps 2–5 (wizard), 8 (🗳️ badge), 9–12 (modal: Pending → vote → Consensus / Disputed), 13 (new suggestion + CAP). This is a 16 kB click-by-click procedure with one picture of an empty panel.

**notes** — §7 documents `▲` / `▼` and "clicking the same button again removes your vote". The buttons are `.community-annotation-vote-btn.vote-up` / `.vote-down` — trivially croppable and currently invisible to the reader. §8 documents Enter-submits / Shift+Enter-newline for comments; the comment bar is `.community-annotation-comment-bar`.

---

## {doc}`/user_guide/web_app/j_community_annotation/02_author_guide`

**explains** — The author's setup and operations: dataset-id stability, repo layout, `annotations/config.json` schema, GitHub App install, self-hosted Worker contract, enabling columns, tuning thresholds, merging duplicates, derived columns, consensus export, and a large troubleshooting catalogue. 38 kB — the largest page in all three sections.

**has** — 1 figure, L61 `community_annotation/disconnected-panel.png`, `:width: 246px`. Same file again.

**verdict** — **SEQUENCE**.

**needs** — lifecycle steps 6, 7, 14, 15, 16 (connected status panel, MANAGE ANNOTATION with the threshold slider, moderation merge, derived consensus column, consensus snapshot download). **A 38 kB author guide illustrated only by the state before the author has done anything is the worst single ratio in the audit.**

**notes**
- **There is a visible gap at lines 214–216** where §"Confirm the dataset id in the UI (recommended)" ends with "the Community Annotation status panel displays the dataset id" followed by two blank lines and a horizontal rule. That reads like a removed figure. The status-panel crop (lifecycle step 6) belongs exactly there.
- Repo-layout claims **verified** against `cellucid-annotation`: `annotations/{config.json, config.schema.json, schema.json}`, `annotations/users/` (contains only `.gitkeep` — **no example user file**), `annotations/moderation/{merges.json, merges.schema.json}`, `.github/workflows/validate.yml`, `scripts/validate_user_files.py`. All present. `merges.json` is `{"version":1,"updatedAt":"1970-01-01T00:00:00Z","merges":[]}`.
- The template `config.json` uses `datasetId: "example-dataset-id"` with `fieldsToAnnotate: ["cell_type","batch"]`, not the `leiden`/`my_atlas_v1` used in this page's examples. Not wrong, but a screenshot of the template file next to the guide's example would remove a real point of confusion for a non-technical author copying the template.
- §6 "self-host the GitHub OAuth + API proxy" documents an eight-route Worker contract and a health JSON. That is a **terminal-screenshot** opportunity per the user's objective: `curl https://<worker>/ | jq` showing `{"status":"ok","service":"Cellucid GitHub Auth","contractVersion":1,…}`. Suggested `community_annotation/worker-health-curl.png`, `:width: 700px`. Redact the origin.

---

## {doc}`/user_guide/web_app/j_community_annotation/03_ui_reference`

**explains** — Explicitly a "**button-by-button** reference" for every Community Annotation surface: the sidebar accordion, the 4-step GitHub sync wizard, the identity modal, MANAGE ANNOTATION, DERIVED CONSENSUS COLUMN, CONSENSUS SNAPSHOT + LOCAL CACHE, and the voting modal — plus a troubleshooting catalogue.

**has** — 1 figure, L34 `community_annotation/disconnected-panel.png`, `:width: 246px`. Same file, fourth use.

**verdict** — **SEQUENCE** — and this is the page where the single-image reuse is *least* defensible. A page whose stated purpose is button-by-button reference has one picture, of the state in which almost none of those buttons exist.

**needs** — essentially the whole lifecycle, but specifically one crop per documented accordion: the wizard's four steps, `Profile (optional)` + the "Your identity" modal, `MANAGE ANNOTATION`, `DERIVED CONSENSUS COLUMN`, `CONSENSUS SNAPSHOT + LOCAL CACHE`, and an annotated voting modal with the suggestion-card anatomy (label / `net …` / `Ontology:` / `Markers:` / `Evidence:` / `▲ n` / `▼ n`) called out.

**notes**
- **There is a visible gap at lines 326–328**, right after the comments bullets and before the `---` that ends §"Voting Modal (Per Category)". Two blank lines where a figure would sit. Same pattern as `02` line 214. Both look like removed images.
- **Prose defect found while gathering selectors: §"Entry button states" lists only two of three states.** The page says the connected entry button reads either `Connect GitHub…` (not signed in) or `GitHub sync…` (signed in + connected). `controls-panel.js` emits **three** labels: `Connect GitHub…`, **`Choose repo…`** (authenticated but no repository chosen), and `GitHub sync…`. The middle state is undocumented — and it is precisely the state a user lands in after signing in but before finishing the wizard, i.e. a likely place to get stuck. Worth a sentence and, ideally, its own crop.
- Everything else on this page that I checked against source is verbatim-accurate: the four wizard step titles, `Add repo` / `Reload` / `Filter repositories…` / `Connect` / `Switch repo` / `Pull latest` / `Publish`, the `Auto pull` toggle with 10/15/60-minute options (`600_000` / `900_000` / `3_600_000` ms), the `Threshold` slider and `Min annotators` input with `Apply` / `Reset`, `Build derived column`, `Consensus snapshot (cellucid-consensus.json)` + `Download`, `Clear session` / `Clear downloads`, and every field label in the "Your identity" modal. The five repository paths validated on connect are also correct.

---

## Judgement on the shared `disconnected-panel.png`

The image itself is not *wrong*: 246×123, correct content — `COMMUNITY ANNOTATION` header, the `i` info button, `No annotation repo connected.`, and the `CONNECT REPO` button. It matches the prose in `03` §"Entry button states" exactly.

It is **wrong to use it four times**. The consequence is that the docs depict this feature *exclusively in its empty state*: there is no published image anywhere in the Cellucid documentation of a connected repository, a signed-in user, a suggestion, a vote, a consensus line, a merge, or a consensus snapshot. A reader who has never used the feature cannot tell from the docs what it looks like when it is working — which is precisely the "wet lab person should understand and use the app" bar in the objective.

Keep the file. Use it **once**, on {doc}`/user_guide/web_app/j_community_annotation/03_ui_reference` §"Entry button states", where it is the literal subject. Replace its three other uses with the state each page is actually about (``index`` → connected status panel; `01` → the wizard; `02` → MANAGE ANNOTATION). Also re-shoot at 2× and give it ~8 px of breathing room — 123 px tall is a sliver.

---

# LIFECYCLE: community annotation

**The headline deliverable.** Ordered image sequence covering the full feature, with the mechanism for reaching each state without a real GitHub repository.

## The mock, verified

Four spec files exist, all in `cellucid/tests/browser/`:

| Spec | Lines | Reaches |
|---|---:|---|
| `community-annotation-startup-lifecycle.spec.mjs` | 4549 (34 tests) | full connection lifecycle **including Publish** |
| `community-annotation-cap-voting.spec.mjs` | 807 (2 tests) | voting modal, suggestions, comments, CAP search, moderation merge |
| `community-annotation-export.spec.mjs` | 69 (1 test) | consensus snapshot download |
| `community-annotation-oauth-worker-lifecycle.spec.mjs` | 308 (2 tests) | Worker source called directly in Node — **no UI, not useful here** |

`cellucid/tests/browser/fixtures/community-annotation-harness.html` is an 11-line bare page (`<main id="harness">`) for module-level tests — **not** useful for screenshots. Screenshots must run against the real app.

There are **two** working strategies. Pick per image.

### :warning: Hard precondition for Strategy A — the local-dev host gate

`window._simulate_repo_connected` and `window._author_mode` are installed by `access-store.js` → `installDevConsoleOverrides()`, and **both `readDevOverrideRole()` and `isSimulateRepoConnectedEnabled()` first call `isLocalDevHost()`** (`cellucid/assets/js/app/utils/local-dev.js`). That returns true only for `file:`, `localhost`, `*.localhost`, `127.x.x.x`, `::1`, `0.0.0.0`, `*.local`, `10.*`, `192.168.*`, `172.16–31.*`.

**On any other origin the flags are silently inert.** The Playwright server binds `127.0.0.1:4173` (`cellucid/scripts/browser-test-ports.mjs`), so it works there — but a screenshot job pointed at a staging URL or `www.cellucid.com` will produce disconnected panels with no error. This is the single most likely way for the capture run to fail confusingly.

### Strategy A — dev-simulate (fastest; no network at all)

```js
await page.addInitScript(() => {
  window._simulate_repo_connected = true;
  window._author_mode = true;            // omit/false for annotator-role captures
  sessionStorage.setItem(
    'cellucid:github-app-auth:session',
    JSON.stringify({ token: 'doc-token', user: { id: 42, login: 'docs-user' } })
  );
});
```
`isAnnotationRepoConnected()` returns true immediately and `getEffectiveRole()` returns `'author'`. **Nothing is fetched** — no `/api/repos`, no `/auth/user`; the token string is never validated. `_author_mode` is what makes **MANAGE ANNOTATION** appear.

Limitation: the status chips read from the simulated state, so a shot that must show a real `owner/repo@branch` and `@login` needs Strategy B.

### Strategy B — synthetic Worker via `page.route` (needed for the wizard and Publish)

```js
await page.addInitScript(({ workerOrigin, user }) => {
  window.__CELLUCID_GITHUB_WORKER_ORIGIN__ = workerOrigin;     // 'https://worker.example'
  sessionStorage.setItem('cellucid:github-app-auth:session',
    JSON.stringify({ token: 'synthetic-token', user }));       // { id: 4242, login: 'tester' }
}, { workerOrigin: 'https://worker.example', user: { id: 4242, login: 'tester' } });
```
The origin override is accepted **only on a local-dev host** (`github-auth.js` throws `GITHUB_WORKER_ORIGIN_UNTRUSTED` elsewhere) — the same gate as Strategy A.

Then route `https://worker.example/**`. Exact fixture shapes, all verified from the spec:

| Route | Response |
|---|---|
| `GET /` (capability handshake, runs **first**) | `{"status":"ok","service":"Cellucid GitHub Auth","contractVersion":1,"endpoints":["/auth/login","/auth/callback","/auth/user","/auth/installations","/auth/installation-repos","/cap/lookup-cells","/cap/search-datasets","/api/repos/*"]}` — **must match `WORKER_ENDPOINTS` element-by-element in order** or sign-in aborts before redirect |
| `GET /auth/user` | `{"id":4242,"login":"tester"}` |
| `POST /auth/installations` | `{"installations":[{"id":71,"account":{"login":"team"}}]}` |
| `POST /auth/installation-repos` (body `{installation_id}`) | `{"repositories":[{"id":1,"full_name":"team/labels","private":false}]}` |
| `GET /api/repos/<owner>/<repo>` | `{full_name, default_branch:"main", private:false, allow_forking:true, permissions:{pull:true,triage:false,push:true,maintain:<author>,admin:<author>}}` — `maintain/admin true` ⇒ role `author`; all false + `allow_forking:false` ⇒ Publish disabled |
| `GET /api/repos/*/*/git/trees/*` | `{"tree":[{"type":"tree","path":"annotations/users","sha":"d…"},{"type":"blob","path":"annotations/users/.gitkeep","sha":"8b137891791fe96927ad78e64b0aad7bded08bdc","size":1}],"truncated":false}` |
| `GET /api/repos/*/*/contents/<path>` | `{type:"file", encoding:"base64", content:<base64 JSON>, sha:<char×40>}` for `annotations/schema.json`, `config.schema.json`, `moderation/merges.schema.json`, `config.json`, `users/ghid_4242.json`, `moderation/merges.json` |
| `PUT /api/repos/*/*/contents/*` (**Publish**) | `200`, headers `x-cellucid-operation-id: <echo>`, `x-cellucid-operation-outcome: applied`, `access-control-expose-headers: X-Cellucid-Operation-Id, X-Cellucid-Operation-Outcome`; body `{"content":{"sha":"4444…"}}` |

The three `*.schema.json` documents only need their identity pair — the app checks `$id` only:
each carries the draft 2020-12 `$schema` plus an `$id` on the frozen bare-apex
contract host, one per document (`user-v1`, `config-v1`, `merges-v1`). The exact
`$id` strings are pinned in ``docs/contributing/web_app.md``; do not retype them here,
because `tests/test_canonical_url_contract.py` pins bare-apex URLs to that page.

Deep link used by the spec:
`/?exportsBaseUrl=<encoded>&dataset=<id>&annotations=team%2Flabels%40main&acceptance=<tag>`

**(B2) Load the app on a real dataset**, then `dismissWelcome(page)` (`#welcome-modal` → Escape).
For documentation use `?dataset=pancreas` (3696 real cells, real cluster labels) — **not** `current-ui-prepared`.

**(C) Drive session state directly from the page context.**
```js
const { getCommunityAnnotationSession } =
  await import('/assets/js/app/community-annotations/session.js');
const s = getCommunityAnnotationSession();
s.setProfile({ username:'ghid_42', githubUserId:42, login:'…', displayName:'…',
               title:'', orcid:'', linkedin:'' });
s.setFieldCategories('clusters', [...]);
s.setFieldAnnotated('clusters', true);
s.setAnnotatableConsensusSettings('clusters', { minAnnotators: 3, threshold: 0.5 });
s.rebuildMergedViewFromUserFiles([ /* fabricated ghid_* user documents */ ]);
```
`rebuildMergedViewFromUserFiles(remoteUserDocs, options)` (`session.js:3283`) is **the key to the whole lifecycle**: it is the exact function a real Pull feeds, so handing it fabricated documents produces a genuine multi-user merged view — multiple voters, real consensus arithmetic, real card ordering — with no network at all. Other useful methods on the same object: `vote(fieldKey, catIdx, suggestionId, direction)`, `addSuggestion`, `addComment`, `computeConsensus(fieldKey, catIdx, {minAnnotators, threshold})` (`:2854`), `addModerationMerge({fieldKey, catIdx, fromSuggestionId, intoSuggestionId, note})` (`:2659`), `buildConsensusDocument()` (`:2954`), `setFieldClosed`.

**The fabricated user document schema** is fixed by `cellucid-annotation/annotations/schema.json` (`user-v1`). Required: `version` (const 1), `username`, `githubUserId` (int ≥1), `updatedAt` (ISO Z), `suggestions`, `votes`. Shape:
```json
{ "version": 1, "username": "ghid_101", "githubUserId": 101,
  "login": "a-reviewer", "displayName": "A. Reviewer", "title": "Postdoc",
  "updatedAt": "2026-01-15T10:00:00Z",
  "suggestions": { "clusters:Beta": [
      { "id": "sug_beta_01", "label": "Beta cell", "ontologyId": "CL:0000169",
        "markers": [{"gene":"INS","logFC":4.2,"pval":1e-40}],
        "evidence": "High INS, IAPP; low GCG.",
        "proposedBy": "ghid_101", "proposedAt": "2026-01-15T10:00:00Z" } ] },
  "votes": { "sug_beta_01": "up" },
  "comments": { "sug_beta_01": [
      { "id":"c1", "text":"Agrees with the reference atlas.",
        "authorUsername":"ghid_101", "createdAt":"2026-01-15T10:05:00Z" } ] } }
```
Bucket keys are `<fieldKey>:<categoryLabel>`; `proposedBy`/`authorUsername` must match `^ghid_[1-9][0-9]*$`; votes are `"up"` | `"down"`. **Every displayed name in a screenshot is chosen by you** — use plausible, clearly-fictional identities and never a real collaborator's GitHub login.

**(D) Route interception, only where a network hop is unavoidable.**
- CAP search: `page.route('**/cap/lookup-cells', …)` → `{ "contractVersion": 1, "results": [...], "omittedInvalidCount": 0 }` (handle the `OPTIONS` preflight with a 204 + CORS headers, as the spec does).
- ORCID lookup in the identity modal goes **directly to `https://pub.orcid.org`** from the browser — route-mock or leave the field untouched.
- The Worker routes (`/auth/installations`, `/auth/installation-repos`, `/api/repos/*`) are only needed for wizard steps 3–4 and for a real Publish.

**Voting-modal entry point** (used by the spec — bypasses the legend entirely):
```js
const { openCommunityAnnotationVotingModal } =
  await import('/assets/js/app/ui/modules/community-annotation-voting-modal.js');
openCommunityAnnotationVotingModal({ defaultCatIdx, defaultFieldKey, state });
```
`state` is a two-method stub: `async ensureFieldLoaded(){}` and `getFields(){ return [{_isDeleted:false, categories:[…], key, kind:'category', loaded:true}] }`. Readiness sentinel: `.community-annotation-voting-detail`.

**Consensus line selector and exact text** (`community-annotation-voting-modal.js:2956–2966`):
`div.community-annotation-consensus.status-consensus|status-disputed|status-pending`, rendering
`Consensus: "Beta cell" (100% • 4 voters)` / `Disputed: "…" (50% • 4 voters)` / `Pending (1 voters)`.

### Crop-target crib sheet (all verified against source)

| Shot | Selector |
|---|---|
| Sidebar accordion (any state) | `#community-annotation-section` (`<details>`; summary text **`Community Annotation`**) |
| Rendered panel body | `#community-annotation-controls` |
| Status chips | `.community-annotation-status-list` (rows `.community-annotation-status-row`, chips `.community-annotation-status-chip--ok/--warn/--danger`, `-key`/`-val`, actions `.community-annotation-status-action`) |
| Wizard (any step) | `page.getByRole('dialog', { name: 'GitHub sync', exact: true })` → `.community-annotation-modal` |
| Stepper only | `.community-annotation-stepper` (4 × `.community-annotation-step[data-step]`, state `--active`/`--done`/`--locked`) |
| Repo chooser grid | `[aria-label="Selectable annotation repositories"]`; cards `.community-annotation-repo-card[data-repo-full-name]` |
| Sync buttons / auto-pull | `.community-annotation-sync-actions`, `.community-annotation-auto-pull-row` |
| Wizard status line | `.community-annotation-wizard-status[role=status]` |
| MANAGE ANNOTATION | `.analysis-accordion-item.open` whose `.analysis-accordion-title` is `MANAGE ANNOTATION` |
| Derived column | same, title `DERIVED CONSENSUS COLUMN` |
| Snapshot + cache | same, title `CONSENSUS SNAPSHOT + LOCAL CACHE` |
| Voting modal | `.community-annotation-modal-overlay`, title **`Community voting`** |
| Voting detail only | `.community-annotation-voting-detail` |
| One suggestion card | `.community-annotation-suggestion-card` (> `-top`, `-label`, `-ontology`, `-net`, `-meta`, `-evidence`, `-actions`) |
| Vote buttons | `.community-annotation-vote-btn.vote-up` / `.vote-down` (text `▲ 7` / `▼ 2`, `aria-pressed`, `.is-mine`, `.is-delegated`) |
| New-suggestion form | `.community-annotation-dashed-box` containing `.community-annotation-new-title` = `New suggestion` |
| CAP results | `.cap-search-results` (items `.cap-search-item`) |
| Comments / merged / evidence | `.community-annotation-secondary-overlay` (titles `Comments · <label>`, `Merged suggestions`, `Evidence`) |
| Identity modal | `.community-annotation-modal.community-annotation-modal--narrow`, title **`Your identity`** |
| Confirm dialog | `.confirm-dialog-overlay` |

**Do not target** `.community-annotation-panel`, `.community-annotation-panel-title`, `.community-annotation-badge*`, `.community-annotation-cluster-*`, `.community-annotation-category-*`, `.legend-vote-btn`, `.community-annotation-repo-row` — these exist in `_community-annotations.css` but **no JS emits them**. There is likewise no `#community-annotation-panel`.

**Exact wizard step titles** (for caption accuracy): `Sign in with GitHub`, `Install the GitHub App`, `Select an annotation repository`, `Sync (pull / publish)`.

**Exact status strings worth capturing**: `Pulled latest annotations.`, `Publish complete.`, `Downloaded cellucid-consensus.json`, `No annotation repo connected.`

---

## The ordered sequence

Dataset for every step: `pancreas`, annotatable field `clusters`, target category `Beta` (a real, biologically meaningful cluster). `deviceScaleFactor: 2` throughout.

### 1. Disconnected — the entry point
- **Before:** app loaded, `#community-annotation-section` collapsed.
- **Action:** click the accordion header.
- **After:** panel reads `No annotation repo connected.` with a `CONNECT REPO` button.
- **Crop:** `#community-annotation-section` (header + body), plus ~8 px padding.
- **Cursor:** on `CONNECT REPO`.
- **File:** `community_annotation/lifecycle-01-disconnected.png`, `:width: 246px`.
- **Mock-reachable:** **YES** — no mock at all. This is the existing image; re-shoot at 2× with the cursor.

### 2. Wizard step 1 — Sign in with GitHub
- **Before:** step 1 of 4 active, not signed in.
- **Action:** click `Connect repo`.
- **After:** `Continue with GitHub` button visible; stepper shows 1/4.
- **Crop:** `.community-annotation-modal` (whole wizard incl. `.community-annotation-stepper` and `.community-annotation-wizard-nav`).
- **Cursor:** on `Continue with GitHub`.
- **File:** `community_annotation/lifecycle-02-signin.png`, `:width: 560px`.
- **Mock-reachable:** **YES** — open the wizard with **no** sessionStorage token. Nothing is fetched until the button is pressed.

### 3. Wizard step 2 — Install the GitHub App
- **Before:** signed in, installations listed.
- **Action:** advance to step 2 (`.community-annotation-step[data-step="2"]` active; title `Install the GitHub App`).
- **After:** `Add repo` and `Reload` (`.community-annotation-reload-btn`) buttons, plus `.community-annotation-repo-grid[aria-label="Available annotation repositories"]` and its pager (`Previous` / `Next`, status `Showing 1–N of M repositories. Page 1 of P.`).
- **Crop:** `.community-annotation-modal`.
- **Cursor:** on `Add repo`.
- **File:** `community_annotation/lifecycle-03-install-app.png`, `:width: 560px`.
- **Mock-reachable:** **YES, fully.** Strategy B + `POST /auth/installations` → `{"installations":[{"id":71,"account":{"login":"team"}}]}`. Shape verified from `community-annotation-startup-lifecycle.spec.mjs`.

### 4. Wizard step 3 — Select an annotation repository
- **Before:** step 3 (title `Select an annotation repository`), repo list populated.
- **Action:** click a repo card.
- **After:** card gets `.community-annotation-repo-card--selected` + `aria-pressed="true"`; the wizard nav button label changes to `Connect` (or `Switch repo` when replacing, then `Connecting…` / `Switching…`).
- **Crop:** `.community-annotation-modal`, ensuring `[aria-label="Selectable annotation repositories"]` with 2–3 cards, the `Filter repositories…` input, and the `Public`/`Private` + `Selected`/`Connected` `.community-annotation-repo-meta-pill`s are all in frame.
- **Cursor:** on the selected card.
- **File:** `community_annotation/lifecycle-04-select-repo.png`, `:width: 560px`.
- **Mock-reachable:** **YES, fully.** Strategy B + `POST /auth/installation-repos` → `{"repositories":[{"id":1,"full_name":"team/labels","private":false}]}` + `GET /api/repos/team/labels` → the repoInfo fixture. Use obviously fictional names; `your-lab/pancreas-annotations` reads better in docs than the spec's `team/labels`.

### 5. Wizard step 4 — Sync (Pull / Publish)
- **Before:** repo connected.
- **Action:** advance to step 4.
- **After:** `Pull latest` and `Publish`, the `Auto pull` toggle and its 10/15/60-minute select, and the wizard status line.
- **Crop:** `.community-annotation-modal` incl. `.community-annotation-sync-actions` and `.community-annotation-auto-pull-row`.
- **Cursor:** on `Pull latest`.
- **File:** `community_annotation/lifecycle-05-sync.png`, `:width: 560px`.
- **Mock-reachable:** **YES.** Strategy A renders the panel on its own (this is what `community-annotation-export.spec.mjs` relies on). Use Strategy B only if the caption needs a real `Pulled latest annotations.` status line, which `:2076` demonstrates.

### 6. Connected — the status panel
- **Before:** connected + signed in.
- **Action:** close the wizard.
- **After:** `.community-annotation-status-panel` shows `Dataset`, `GitHub` (login), `Repo` (`owner/repo@branch`), the `Copy share link` icon, and the entry button now reading `GitHub sync…`.
- **Crop:** `#community-annotation-controls`.
- **Cursor:** on the `Copy share link` icon.
- **File:** `community_annotation/lifecycle-06-connected-status.png`, `:width: 246px`.
- **Mock-reachable:** **YES** — levers (A)+(C) only.
- **This is the image {doc}`/user_guide/web_app/j_community_annotation/02_author_guide` line 214 is missing** and the natural replacement for the disconnected panel on {doc}`/user_guide/web_app/j_community_annotation/index`.

### 7. Author controls — MANAGE ANNOTATION
- **Before:** connected as author.
- **Action:** expand `MANAGE ANNOTATION`, select `clusters` in `Categorical obs:`, click `Add`.
- **After:** `Close`/`Reopen` available and the per-column **Annotatable consensus settings** appear: `Threshold` slider, `Min annotators` input, `Apply`, `Reset`.
- **Crop:** the `MANAGE ANNOTATION` accordion body.
- **Cursor:** on the `Threshold` slider thumb.
- **File:** `community_annotation/lifecycle-07-manage-annotation.png`, `:width: 246px`.
- **Mock-reachable:** **YES** — `window._author_mode = true` + `s.setFieldAnnotated('clusters', true)`.

### 8. The 🗳️ badge in the field dropdown
- **Before:** `clusters` annotatable.
- **Action:** open the categorical field select.
- **After:** `clusters` carries 🗳️; a closed column would carry 🗳️🏁.
- **Crop:** the open dropdown list.
- **Cursor:** on the 🗳️ entry.
- **File:** `community_annotation/lifecycle-08-votable-badge.png`, `:width: 300px`.
- **Mock-reachable:** **YES**. Consider a second frame with `s.setFieldClosed('clusters', true)` to show 🗳️🏁 — `01` and ``index`` both cite it as a top troubleshooting cause and it has never been shown.

### 9. Voting modal — Pending
- **Before:** `minAnnotators: 3`, one fabricated voter.
- **Action:** click the `Beta` category in the legend.
- **After:** modal open; `.community-annotation-consensus.status-pending` reads `Pending (1 voters)`; one suggestion card.
- **Crop:** `.community-annotation-voting-detail`.
- **Cursor:** no.
- **File:** `community_annotation/lifecycle-09-pending.png`, `:width: 560px`.
- **Mock-reachable:** **YES** — `rebuildMergedViewFromUserFiles([oneUserDoc])`.

### 10. Casting a vote
- **Before:** step 9 state.
- **Action:** hover `▲` on the leading card.
- **After:** unchanged (capture the hover, not the result) — this documents `01` §7.
- **Crop:** a single `.community-annotation-suggestion-card`, showing label, `net …`, `Ontology:`, `Markers:`, `Evidence:`, `▲ n`, `▼ n`.
- **Cursor:** **yes — on `.community-annotation-vote-btn.vote-up`.** This is the most cursor-dependent image in the whole audit.
- **File:** `community_annotation/lifecycle-10-vote.png`, `:width: 460px`.
- **Mock-reachable:** **YES** — `s.vote('clusters', catIdx, id, 'up')` for the after-state if a second frame is wanted.

### 11. Consensus reached
- **Before:** four fabricated voters all upvoting `Beta cell`.
- **Action:** none.
- **After:** `.community-annotation-consensus.status-consensus` reads `Consensus: "Beta cell" (100% • 4 voters)`.
- **Crop:** `.community-annotation-voting-detail`, whole modal.
- **Cursor:** no.
- **File:** `community_annotation/lifecycle-11-consensus.png`, `:width: 560px`.
- **Mock-reachable:** **YES** — four user documents, each with `"votes": {"sug_beta_01": "up"}`, through `rebuildMergedViewFromUserFiles`.

### 12. Disputed
- **Before:** four voters, two upvoting `Beta cell`, two upvoting `Beta (immature)`.
- **After:** `.status-disputed` reads `Disputed: … (50% • 4 voters)` with two tied cards.
- **Crop:** `.community-annotation-voting-detail`.
- **Cursor:** no.
- **File:** `community_annotation/lifecycle-12-disputed.png`, `:width: 560px`.
- **Mock-reachable:** **YES**. This is the exact row 5 of `02`'s worked-examples table ("4 users split: 2 upvote A, 2 upvote B → Disputed") and pairs with it directly.

### 13. Propose a suggestion (+ CAP search)
- **Before:** modal open.
- **Action:** fill `Label`, `Ontology id`, `Marker genes`, `Evidence`; then press `Search CAP`.
- **After:** `.cap-search-results` populated.
- **Crop:** two images — (a) the filled `.community-annotation-new` form, cursor on `Add`; (b) the CAP result list, cursor on a `.cap-search-item`.
- **Files:** `community_annotation/lifecycle-13a-new-suggestion.png`, `lifecycle-13b-cap-search.png`, `:width: 460px`.
- **Mock-reachable:** **YES** — `page.route('**/cap/lookup-cells', …)` exactly as the spec does, including the `OPTIONS` 204. Note `01` §"About CAP search" documents the privacy surface; the caption should not imply CAP was really queried.

### 14. Moderation — merge duplicates
- **Before:** two duplicate cards (`Beta cell`, `Beta-cell`), author role.
- **Action:** drag one card onto the other; confirm; add a note.
- **After:** a `Merged bundle …` row (`.community-annotation-bundle-row`) with a `View merged` button; de-duplicated bundle totals.
- **Crop:** the merged card + bundle row; a second frame of the `View merged` secondary modal (`.community-annotation-secondary-modal`).
- **Cursor:** on `View merged`.
- **Files:** `community_annotation/lifecycle-14a-merge.png`, `lifecycle-14b-view-merged.png`, `:width: 460px`.
- **Mock-reachable:** **YES** — `s.addModerationMerge({fieldKey, catIdx, fromSuggestionId, intoSuggestionId, note})`. Entirely local; `merges.json` is only written on Publish.

### 15. Derived consensus column
- **Before:** consensus computed.
- **Action:** open `DERIVED CONSENSUS COLUMN`, pick `clusters`, type `community_cell_type`, set threshold/min annotators, click `Build derived column`.
- **After:** the new column exists; the embedding recolours to consensus labels plus `Disputed` / `Pending`.
- **Crop:** two images — the accordion body (`:width: 246px`), and the recoloured `#glcanvas` + legend (`:width: 900px`). The second is the payoff image for the entire feature: *this is what community annotation gets you.*
- **Cursor:** on `Build derived column` in the first.
- **Files:** `community_annotation/lifecycle-15a-derived-column-controls.png`, `lifecycle-15b-derived-column-view.png`.
- **Mock-reachable:** **YES** — purely local.

### 16. Export the consensus snapshot
- **Before:** consensus computed.
- **Action:** open `CONSENSUS SNAPSHOT + LOCAL CACHE`, click `Download`.
- **After:** `cellucid-consensus.json` downloaded.
- **Crop:** the accordion body showing `Consensus snapshot (cellucid-consensus.json)` + `Download`, and both `Clear session` / `Clear downloads`.
- **Cursor:** on `Download`.
- **File:** `community_annotation/lifecycle-16-consensus-snapshot.png`, `:width: 246px`.
- **Mock-reachable:** **YES** — `buildConsensusDocument()` is local. `community-annotation-export.spec.mjs` already exercises this. Optionally pair with a **terminal screenshot** of `jq '.consensus' cellucid-consensus.json` showing the `status/label/confidence/voters/netVotes/suggestionId` shape that `02` §12 documents in prose — a strong candidate for the objective's "terminal screenshots" ask.

### 17. Publish — direct push
- **Before:** local changes pending, `permissions.push: true`.
- **Action:** `GitHub sync…` → step 4 → `Publish`.
- **After:** `.community-annotation-wizard-status` reads **`Publish complete.`** The client issued `PUT /api/repos/<owner>/<repo>/contents/annotations/users/ghid_<id>.json` with an `x-cellucid-operation-id` UUIDv4.
- **Crop:** `.community-annotation-modal` with the status line in frame.
- **Cursor:** on `Publish`.
- **File:** `community_annotation/lifecycle-17-publish.png`, `:width: 560px`.
- **Mock-reachable:** **YES, fully.** Strategy B `handleMutation` returns `200` + `{"content":{"sha":"4444…"}}` with `x-cellucid-operation-outcome: applied`. Exercised by `community-annotation-startup-lifecycle.spec.mjs:3694` *"an exact-byte no-op publication seeds the SHA baseline for the next edit"* and `:4031` for the author's `config.json` write.
- **Also worth capturing:** the **Publish-disabled** state, which `01` and `03` both troubleshoot at length. Serve repoInfo with `permissions` all false **and** `allow_forking: false` → `selectAnnotationPublicationMode()` returns `null` and Publish is unavailable. `:2076` *"read-only collaborator keeps Pull and the repository while Publish stays explicitly unavailable"* is exactly that state. → `community_annotation/lifecycle-17b-publish-unavailable.png`.
- **Still NOT reachable:** a genuine **fork + Pull Request** confirmation. The route is selected by `selectAnnotationPublicationMode()` from repository metadata, but no spec drives a fork/PR mutation end to end, so the fork-creation and PR-creation responses are unspecified. **Do not fabricate a "Pull request opened" screenshot** — document the PR path in prose only.

### 18. Profile / "Your identity" modal
- **Action:** `Profile (optional)` → `Edit`.
- **After:** modal titled **`Your identity`** with `Name:`, `Affiliation / role:`, `LinkedIn:`, `ORCID:`, `Save` / `Cancel` — all four label strings confirmed verbatim against `identity-profile-modal.js`, matching `03` exactly.
- **Crop:** `.community-annotation-modal.community-annotation-modal--narrow`; **include `.community-annotation-external-lookup-disclosure`**, whose text is `Typing 3 or more characters in Name or ORCID searches the public ORCID registry. Requests omit credentials and referrer information.` — `03` documents this disclosure at length and it has never been shown.
- **Cursor:** on `Save`.
- **File:** `community_annotation/lifecycle-18-identity.png`, `:width: 460px`.
- **Mock-reachable:** **YES for the empty/filled form** (Strategy A is enough). The ORCID **auto-suggest listbox** (`[role=listbox][aria-label="ORCID search suggestions"]`) requires `page.route('https://pub.orcid.org/**', …)` because `orcid-client.js` calls ORCID **directly from the browser**, not through the Worker; `/v3.0/{id}/person` and `/v3.0/expanded-search/` response shapes are parsed by `parseOrcidPersonName` / `parseOrcidExpandedSearch` but were not exercised by any browser spec. Ship the plain form; add the suggest dropdown only after verifying that contract.

## Lifecycle summary

**18 numbered steps (20 images). All 18 are mock-reachable without a real GitHub repository.**

- **13 need no network at all** — Strategy A only: 1, 2, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, and the form half of 18.
- **5 need `page.route` against the synthetic Worker** — Strategy B, with every response shape verified from `community-annotation-startup-lifecycle.spec.mjs`: 3, 4, 5, 17, and 17b.
- **1 needs a non-Worker route mock that *is* verified** — 13b (CAP), shape taken from `community-annotation-cap-voting.spec.mjs`.
- **2 sub-states remain genuinely unreachable and must not be faked**: the **fork + Pull Request** confirmation (step 17, no spec drives it; shapes unspecified) and the **ORCID auto-suggest listbox** (step 18, direct browser→`pub.orcid.org` call, no browser spec covers it).

Both exclusions are narrow: the *feature* is fully documentable in pictures, and only two optional sub-states are not.

---

# Section K — `user_guide/web_app/k_figure_export/`

## {doc}`/user_guide/web_app/k_figure_export/index`

**explains** — What figure export captures and does not, plot size vs DPI arithmetic, filenames and privacy, and the reading order.

**has** — 1 figure, L109 `figure_export/panel-overview.png`, `:width: 246px` (file 246×256). **Also used on `08_screenshots.md` L5.**

**verdict** — **OK** (content correct; global 2× re-shoot applies).

I checked it at full size: it shows the `FIGURE EXPORT` header, the intro sentence, the four collapsed sections (`FRAMING` / `LABELS & ANNOTATIONS` / `STYLE` / `DOWNLOAD`) with their summary hints, and the `EXPORT` button. That is exactly what the caption claims and exactly what the source builds (`figure-export-ui.js:3253–3271`, titles `'Framing'`, `'Labels & Annotations'`, `'Style'`, `'Download'`). Good crop, no bleed.

**needs** — re-shoot only, at `deviceScaleFactor: 2`, cropped to `#figure-export-section`. Optionally add a cursor on `EXPORT`.

**notes** — The DPI arithmetic (`output_pixels ≈ plot_pixels × DPI/96`) is stated three times across the section (``index``, `02`, `04`) with the same worked example. Consistent — no defect, but a single small annotated diagram would serve all three; that is the one place in section K where a **DIAGRAM** would beat a screenshot.

---

## {doc}`/user_guide/web_app/k_figure_export/01_figure_export_goals_wysiwyg_and_reproducibility`

**explains** — WYSIWYG vs reproducible, the three levels of reproducibility, the capture checklist, what is not reproduced, and a packaging recipe.

**has** — `none`.

**verdict** — **ADD** (2 images).

**needs**
1. **Viewer ↔ exported figure side by side.** The page's central promise is "the exported figure should match the state you are currently looking at". Capture the Pancreas 2D viewer, export a PNG of that view, and place them side by side. Dataset `pancreas`, colour by `clusters`, one filter active so the `(hidden)` legend marking is visible. → `figure_export/wysiwyg-viewer-vs-export.png`, `:width: 900px`. No cursor.
2. **An `Export blocked` dialog.** The page introduces blocking as a design principle ("Where a layer cannot be reproduced, Cellucid says so instead of quietly dropping it"). Easiest reproduction: enable the velocity overlay on `pancreas` and press `Export` → the `overlay-fidelity.js` blocker fires with *"The velocity vector field is enabled in the viewer, but figure export… Turn the velocity overlay off before exporting."* Crop `.fidelity-warning-dialog` (component at `figure-export/components/fidelity-warning-dialog.js`). Cursor: no. → `figure_export/export-blocked-velocity.png`, `:width: 560px`. **This image serves four pages** (`01`, `04`, `06`, `07`) — capture once, reference from all.

**notes** — §"WYSIWYG" claims interaction chrome "is disclosed in the Annotations block". That disclosure is real (`figure-export-ui.js:582–586`, `Never exported: …`) but is **absent from the one screenshot of that block** — see `k/05`.

---

## {doc}`/user_guide/web_app/k_figure_export/02_export_ui_walkthrough`

**explains** — The click-by-click export: where the panel is, preview vs export, size/DPI, titles/legends/axes/reference grid, selection emphasis, download options, the Frame-export crop workflow, and multi-panel legend rules.

**has** — 1 figure, L395 `figure_export/preview.png`, `:width: 1440px` (file 1440×900). Also on `08_screenshots.md` L15.

**verdict** — **REPLACE**, plus **SEQUENCE** for the crop workflow.

**Why REPLACE — I read the file at full size and it is the worst image in the audit:**
1. **It is the 120-point `current-ui-prepared` fixture.** Roughly 95 % of the 1440×900 frame is empty white grid with ~120 scattered single-pixel dots. A wet-lab reader learns nothing about what a figure export looks like.
2. **The sidebar is clipped mid-section.** `FRAMING` runs off the bottom edge; the last legible line is `Preview sample: 120 points • shader-accurate points`. `LABELS & ANNOTATIONS`, `STYLE`, `DOWNLOAD` and the `EXPORT` button are all out of frame.
3. **Stray chrome.** The same unexplained floating ✕ button at ~(309, 27) as the velocity captures.
4. **The caption does not describe the image.** Both captions say Preview "exposes the current framing and visual state before any file is written". The large area in the image is the **live viewer canvas**, not the export preview; the actual preview is the ~175 px thumbnail buried inside `FRAMING`. As shipped, the figure invites the reader to mistake the viewer for the preview.
5. 1×.

**needs**
- *Replacement for the page's "Interface reference"*: `pancreas`, coloured by `clusters`, 2D Planar, `Figure Export` panel open with `FRAMING` expanded and `SHOW PREVIEW` checked, sidebar scrolled so the preview thumbnail is fully visible, viewport cropped so the canvas shows real data rather than empty grid. Cursor: no. → keep the filename `figure_export/preview.png` (so `08` picks it up too), `:width: 1440px`, captured at 2×.
- *The Frame-export sequence* — §"How to use Frame export (step-by-step)" is six numbered steps and one of the most-asked features; it has no images:
  1. `SHOW PREVIEW` on, `FRAME EXPORT` off. Crop the `FRAMING` section. Cursor on `FRAME EXPORT` (`#figure-export-crop-enabled`). → `figure_export/crop-01-enable.png`
  2. `FRAME EXPORT` on, drag handle mid-drag on `#figure-export-plot-frame`. Crop the preview area. Cursor on a corner handle. → `figure_export/crop-02-drag.png`
  3. After `Confirm` — the framed preview. Crop the preview area. Cursor on `Confirm`. → `figure_export/crop-03-confirmed.png`
  Each `:width: 300px`; `#figure-export-crop-lock` (`Lock aspect`) should be visible in all three.
- *Selection emphasis* — §"Selection emphasis" documents the `n = 1,234 selected` badge and the α input. Select a cluster, enable `#figure-export-selection-emphasis`, export, and show the exported figure with its badge. → `figure_export/selection-emphasis-export.png`, `:width: 700px`.
- *Multi-panel legends* — §"Legends in a multi-panel export" is the densest passage in section K (shared legend vs per-panel legend; the 96 px / 44 px legibility floor). Two exported grids: (a) panels agreeing → one shared legend beside the grid; (b) panels with different hidden categories → a legend inside each cell. → `figure_export/multipanel-legend-{shared,per-panel}.png`, `:width: 700px` each. `cellucid/tests/browser/figure-export-multiview-parity.spec.mjs` is the reference for reaching grid-compare state.

**notes** — Everything documented on this page that I checked against source is accurate: the five size presets, `Auto text`, the `Legend: Right/Bottom` select, `Reference grid`, `3D orientation`, `Depth sort`, `Emphasize selection` + α, the five download modes, the three SVG strategies, and the `Keep` target (`#figure-export-target-group`, default `100000`, min 1000, max 5000000, step 1000).

---

## {doc}`/user_guide/web_app/k_figure_export/03_export_formats_and_renderers`

**explains** — SVG vs PNG, the three SVG strategies, why SVGs get huge, browser contract, and Illustrator/Inkscape notes.

**has** — 1 figure, L192 `figure_export/download.png`, `:width: 224px` (file 224×246). Also on `08_screenshots.md` L47.

**verdict** — **REPLACE** — **STALE-RISK, confirmed stale.**

The image shows the hint ending at "…Hybrid embeds the shader-rendered point pass." The current source (`figure-export-ui.js:761–765`) appends: *"Vector points carry the viewer's atmospheric fog but are flat discs: only Hybrid and PNG reproduce the 3D sphere lighting."* That sentence is the exact claim `k/06` §*Atmospheric fog and 3D lighting in vector exports* makes. The screenshot predates it.

**needs**
- *Re-shoot* `download.png`: `DOWNLOAD` section expanded, `SVG + PNG` / `300` / `Full vector`, full hint text visible, `EXPORT ALL VIEWS (SPLIT-VIEW)` in frame. Crop the `DOWNLOAD` disclosure. Cursor: no. `:width: 224px`, 2×.
- *Add a second `DOWNLOAD` state*: set the strategy select to `Optimized vector` so `#figure-export-target-group` (`Keep` + the numeric input, default `100000`) is revealed. **That control is documented on this page's sibling `02` and in `07` §"tune the 'Keep' point count" and has never been shown** — it is `display:none` under any other strategy, so the existing screenshot could not contain it. → `figure_export/download-optimized-keep.png`, `:width: 224px`. Cursor on the `Keep` input.
- *Add a strategy comparison*: the same dense Pancreas region exported three ways — Full vector / Optimized vector (Keep 5,000) / Hybrid — cropped identically and shown with their file sizes in the caption. The decision table at the top of this page is the section's main navigational aid and is entirely numeric. → `figure_export/strategy-{full,optimized,hybrid}.png`, `:width: 300px` each.

**notes** — The caption on this page ("Download controls with output format, DPI, and split-view choice") and the different caption on `08` for the same file ("with SVG plus PNG, 300 DPI, Full vector strategy, and split-view choice") are both accurate to the image. Fine, but two captions for one file is a maintenance trap: when the image is re-shot, both must be checked.

---

## {doc}`/user_guide/web_app/k_figure_export/04_quality_knobs_and_best_practices`

**explains** — Vector vs raster decision, sizing recipes, large point clouds without misleading downsampling, colour/contrast/accessibility, and safe vs dangerous post-processing.

**has** — 1 figure, L220 `figure_export/style.png`, `:width: 224px` (file 224×450). Also on `08_screenshots.md` L40.

**verdict** — **OK** (content correct; global 2× re-shoot applies) — **plus ADD**.

The image is a clean crop of the expanded `STYLE` section: `BACKGROUND` = `Match viewer`, `FONT` = `Arial`, `TEXT SIZES` with `AUTO TEXT` checked and Base 15 / Legend 15 / Ticks 15 / Axis 17 / Title 20 / Centroids 11 greyed out. Matches the prose exactly.

**needs**
1. **The colourblind simulation.** §"Colorblindness preview (high value, low effort)" calls this out as the cheapest quality win in the section and shows nothing. Enable `SHOW PREVIEW`, set the select (`.figure-export-colorblind-select`, aria-label `Colorblind preview`) to `Normal` then `Deuteranopia`, and capture the preview thumbnail in both states on a `clusters`-coloured Pancreas view where two categories genuinely collide. → `figure_export/colorblind-{normal,deuteranopia}.png`, `:width: 300px` each. Cursor on the select in the first. **This is the highest-value addition in section K** — it is a *quality* feature the reader cannot evaluate from words.
2. **The `(hidden)` legend entry.** Both this page and `06` insist that hiding a category keeps its legend entry with a hollow grey swatch. Export a figure with two categories hidden and crop the legend. → `figure_export/legend-hidden-entries.png`, `:width: 300px`.
3. **The `Log10 color scale` colorbar.** Same argument: documented twice in prose, never shown. Colour by a continuous field with `Use log scale` on, export, crop the colorbar. → `figure_export/legend-log-colorbar.png`, `:width: 300px`.

**notes** — There is a **visible gap at lines 226–228** (two blank lines between the figure block and the `---`), matching the pattern in `j/02` and `j/03`. Possibly a removed second figure.

---

## {doc}`/user_guide/web_app/k_figure_export/05_metadata_and_provenance`

**explains** — What provenance means, the per-panel metadata model, PNG `iTXt` fields, the `Description` string format, the structured `Comment` JSON, SVG RDF/Dublin Core, filename rules, and a methods-text template.

**has** — 1 figure, L262 `figure_export/labels-annotations.png`, `:width: 224px` (file 224×397). Also on `08_screenshots.md` L32.

**verdict** — **REPLACE** — **STALE-RISK, confirmed stale** (the strongest staleness finding in the audit).

The image shows `TITLE:` = `Current UI prepared fixture`, then `ANNOTATIONS:` with exactly four checkboxes — `AXES`, `LEGEND`, `3D ORIENTATION`, `DEPTH SORT` — then axis-label fields and `Legend: Right`, then `EMPHASIZE SELECTION`.

Current source builds `annotationChecks` as **five** entries in this order (`figure-export-ui.js:548–553`): `Axes`, `Legend`, **`Reference grid`**, `3D orientation`, `Depth sort`; and appends a `chromeDisclosure` hint reading `` `Never exported: ${NON_EXPORTED_VIEWER_CHROME.join(', ')}.` `` (lines 580–591). **Neither is in the image.** Meanwhile `k/02` documents "**Reference grid** is enabled by default", `k/06` devotes a subsection to it, and `k/01` says the chrome list "is disclosed in the Annotations block". The prose moved; the picture did not.

The image also has a second, softer problem shared with the whole section: the `TITLE:` field displays `Current UI prepared fixture`, leaking the 120-point test fixture's name into user documentation.

**needs**
- *Re-shoot* on `pancreas`: `LABELS & ANNOTATIONS` expanded, title auto-filled from the real dataset, **all five** checkboxes visible, the `Never exported: …` hint line in frame, `Legend: Right`, `EMPHASIZE SELECTION` visible. Crop the disclosure. Cursor: no. Keep the filename, `:width: 224px`, 2×.
- *Add a metadata inspection image.* The page's whole subject is embedded provenance and there is no picture of any. A **terminal screenshot** is the right medium and matches the user's stated objective: `exiftool -iTXt:all figure.png` or `python -c "import json,png…"` showing `Software`, `Creation Time`, `Dataset`, `Color Field`, `Description`, and the pretty-printed `Comment` JSON with its `views[]` array. → `figure_export/metadata-inspect-terminal.png`, `:width: 900px`. **Redact local paths in the capture** — the page itself warns that `Source File` may be a local path.

**notes** — The `Description` grammar quoted on this page (`Views: A. Live (Field: …; Filters: …) | B. … • Source: …`) and the `Comment` JSON keys are consistent with `09`'s description of `utils/figure-provenance.js` and `unanimousFieldKey()`. No contradiction found between these two pages.

---

## {doc}`/user_guide/web_app/k_figure_export/06_edge_cases`

**explains** — 0 visible points, tiny groups, NaN/Inf, missing fields, SVG size explosion, legend overflow, many panels, reference grid, 3D shader points, atmospheric fog vs 3D lighting, depth ordering, 3D axes, smoke/connectivity/velocity blockers, viewer chrome, dark backgrounds, WebGL2, fonts, downloads, privacy.

**has** — `none`.

**verdict** — **ADD** (4 images), selectively.

**needs**
1. **The three blockers.** §"Volumetric smoke cloud render mode not exported", §"Connectivity overlay not exported", §"Velocity vector field not exported" each quote an exact blocker string. Capture the `.fidelity-warning-dialog` once per blocker (velocity is the cheapest — see `k/01`; smoke via `Render mode:` = `Volumetric smoke cloud (alpha)`; connectivity via the KNN connectivity overlay). → `figure_export/blocked-{velocity,smoke,connectivity}.png`, `:width: 560px` each. Reuse the velocity one from `01`.
2. **Reference grid on vs off.** §"Reference grid" explains that the viewer's `Background` control decides both the clear colour and whether the grid box is drawn, and that the checkbox is disabled for `White`/`Black`. Two exported figures — `Grid (light)` with `Reference grid` on, and the same view with it off — plus, ideally, a third crop showing the checkbox **disabled** under `Background: White`. → `figure_export/reference-grid-{on,off,disabled}.png`, `:width: 300px` each.
3. **Dark-background ink.** §"Dark-background figures" claims titles/ticks/axis labels/frame/badge take their colour from the figure background. Export the same view under `Background: Match viewer` from `Grid (dark)` and from `Grid (light)`, side by side. → `figure_export/dark-vs-light-ink.png`, `:width: 700px`.
4. **Category explosion.** §"Category explosion (legend overflow)" is about an unreadable figure; showing one is legitimate and instructive. Export a `clusters`-style field with many categories at `Legend: Right` vs `Legend: Bottom`. → `figure_export/legend-overflow-{right,bottom}.png`, `:width: 400px` each.

**NONE (justified)** for NaN/Inf, font substitution, browser download restrictions and the privacy subsection: all four are either invisible, environment-specific, or already covered by the terminal capture proposed on `05`.

**notes** — §"Exporting with 0 visible points" and §"Missing/renamed/deleted fields" both describe outcomes that are trivially screenshottable but low value — a picture of an empty plot teaches nothing the sentence does not. Deliberately excluded.

---

## {doc}`/user_guide/web_app/k_figure_export/07_troubleshooting_figure_export`

**explains** — Quick triage, the `Export blocked` catalogue, and eleven symptom → cause → confirm → fix entries.

**has** — `none`.

**verdict** — **SEQUENCE** (the triage checklist) + reuse.

**needs**
1. **The `Export blocked` dialog** — reuse `figure_export/export-blocked-velocity.png` from `k/01`. §"If you see 'Export blocked'" lists four blockers and shows none; one real dialog anchors all four.
2. **The triage checklist as a strip.** §"Quick triage checklist" step 3 says "Export a small **PNG** at a small plot size (600×450, 150 DPI)". Three crops of the `FRAMING` + `DOWNLOAD` sections configured exactly that way, so a stuck user can pattern-match their own panel against the picture. Cursor on the `W` input (`aria-label` `Export width (px)`) in the first. → `figure_export/triage-{size,dpi,export}.png`, `:width: 224px` each.
3. **`[FigureExport]` console errors.** Four separate entries tell the reader to "check the browser console" and "look for `[FigureExport]` errors". A devtools-console screenshot showing one real `[FigureExport]` line is the fastest way to teach a wet-lab user what that instruction even means. → `figure_export/console-figureexport-error.png`, `:width: 700px`.

**notes** — §"Symptom: 'Export button does nothing'" and §"'Export failed'" both reference a `#figure-export-btn`-driven flow; the button is visible in `panel-overview.png`, so a cross-reference to that figure would already help without new capture.

---

## {doc}`/user_guide/web_app/k_figure_export/08_screenshots`

**explains** — The section's "verified captures" page: six figures covering panel structure, preview/framing, and labels/style/download.

**has** — 6 figures, **all reused from other pages**, none unique to this page:
| Line | File | `:width:` | Also on |
|---:|---|---|---|
| 5 | `figure_export/panel-overview.png` | `246px` | {doc}`/user_guide/web_app/k_figure_export/index` L109 |
| 15 | `figure_export/preview.png` | `1440px` | `02` L395 |
| 22 | `figure_export/framing.png` | `224px` | **nowhere else — unique to this page** |
| 32 | `figure_export/labels-annotations.png` | `224px` | `05` L262 |
| 40 | `figure_export/style.png` | `224px` | `04` L220 |
| 47 | `figure_export/download.png` | `224px` | `03` L192 |

**verdict** — **REPLACE** (inherits every problem catalogued above).

Specifically: it carries both confirmed-stale images (`labels-annotations.png`, `download.png`), the badly-framed `preview.png`, and all six at 1×. Its own captions must be re-checked after each re-shoot — `download.png`'s caption here is more specific than the one on `03` and will need to stay true.

**needs**
- All six re-shot per the per-page entries above, on `pancreas` rather than the 120-point fixture.
- `framing.png` — the only file unique to this page. Content is correct (`PLOT SIZE: Screen (half)`, `W 1200 / H 887`, `SHOW PREVIEW` checked, colourblind select, `REFRESH`, `FRAME EXPORT` unchecked, the preview thumbnail, and `Preview sample: 120 points • shader-accurate points`). But **the thumbnail shows the near-empty 120-point fixture and the status line says so in words** — on `pancreas` it would read a real point count and show a real embedding. Re-shoot. Cursor on `FRAME EXPORT`.
- Its caption "Cellucid Figure Export preview for the deterministic 120-cell fixture" is *honest* about the fixture — which is exactly the tell that the wrong dataset was used for user-facing documentation.

**notes** — Compare with `i/08`, which is the same idea done better: three unique images, none reused. `k/08` is a pure index of images owned by other pages. Consider whether the page should survive re-shooting at all, or become a short "how these captures were produced" page (dataset, viewport, DPR, app version) — which is what `r_screenshot_checklist/` was presumably meant to hold.

---

## {doc}`/user_guide/web_app/k_figure_export/09_reference_implementation_notes`

**explains** — Developer reference: the three-layer architecture, module paths, the payload/metadata builder, known limitations with their exact constants, and the minimum bug-repro package.

**has** — `none`.

**verdict** — **NONE (justified).**

The page is a file-path map and a constants list for developers reading source (`LAYOUT_CONSTANTS`, `PANEL_LEGEND_MAX_SHARE` 0.5, `PANEL_LEGEND_GAP` 8 px, `PANEL_LEGEND_MIN_WIDTH` 96 px, `PANEL_LEGEND_MIN_HEIGHT` 44 px). There is no UI state on it. A screenshot would be decoration.

One partial exception: §"What to capture when filing a bug" asks for "screenshot of fidelity warning dialog text if it appeared" — satisfied by cross-referencing the `Export blocked` image from `k/01` rather than adding a new one.

**notes** — All twenty-odd module paths listed on the page resolve to real files under `cellucid/assets/js/app/ui/modules/figure-export/`, including `utils/overlay-fidelity.js`, `utils/render-mode.js`, `utils/viewer-background.js`, `components/reference-grid.js`, `utils/figure-provenance.js`, `utils/panel-label.js`, `utils/selection-badge.js`, `utils/fog.js`, `utils/figure-ink.js` and `utils/zip-archive.js`. No broken references found on this page.

---

# Consolidated new-image count

| Section | Re-shoots of existing | New images | Total captures |
|---|---:|---:|---:|
| I — velocity | 3 | ~34 (incl. 17 pair halves + 5 group crops) | ~37 |
| J — community annotation | 1 | ~24 (18 lifecycle steps; several have 2 frames) | ~25 |
| K — figure export | 6 | ~30 | ~36 |
| **Total** | **10** | **~88** | **~98** |

# The three biggest gaps

1. **Community annotation is documented exclusively in its empty state.** 82 kB of prose across four pages, one 246×123 image of a disconnected panel, reused four times. No published image anywhere shows a connected repo, a suggestion, a vote, a consensus line, a merge, or a derived consensus column. **There is no technical reason for this**: all 18 lifecycle steps are mock-reachable, 13 of them with no network whatsoever (`window._simulate_repo_connected` + a `sessionStorage` token + `session.rebuildMergedViewFromUserFiles()`), and the remaining five with `page.route` fixtures whose exact shapes are already written down in `community-annotation-startup-lifecycle.spec.mjs`. The one trap is the `isLocalDevHost()` gate — the capture run must be served from `127.0.0.1`/`localhost` or the dev flags are silently inert.
2. **Velocity documents 26 controls and illustrates none of them individually.** Eight of nine pages have zero figures. Every "reduce turbulence for honest reading", "bloom washes out the data", "trail length is not trail persistence" claim in the section is unevidenced — and those are scientific-integrity claims, not cosmetic ones.
3. **Two figure-export screenshots are demonstrably stale, and the section's flagship image shows a 120-point test fixture.** `labels-annotations.png` is missing the `Reference grid` checkbox and the `Never exported: …` disclosure that three pages of prose describe; `download.png` is missing the atmospheric-fog sentence that `06` describes. `preview.png` shows a near-empty synthetic canvas with a clipped sidebar and a stray ✕ button, under a caption describing something else in the frame.
