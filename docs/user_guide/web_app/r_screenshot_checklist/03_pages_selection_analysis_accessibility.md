# Screenshot coverage audit — sections F, H, O

Scope: `user_guide/web_app/f_highlighting_selection/`,
`user_guide/web_app/h_analysis/`,
`user_guide/web_app/o_accessibility_privacy_security/`.
All paths below are relative to `cellucid-python/docs/`.
Read-only pass; the only file written is this one.

## Summary

| Metric | Value |
| --- | --- |
| Pages audited | **23** (F=8, H=12, O=3) |
| `{figure}` directives present | **17** (F=1, H=15, O=1) |
| Unique image files behind them | **9** |
| Pages with **zero** figures | **12 of 23** (F: 7/8 · H: 3/12 · O: 2/3) |
| Images captured at 1× DPR (blurry on HiDPI) | **9 of 9** — every one |
| Images demonstrably stale vs current app/datasets | **4** |

### Verdicts

| Verdict | Count | Pages |
| --- | ---: | --- |
| OK | **0** | — |
| REPLACE | **10** | `f/07`, ``h/index``, `h/02`, `h/03`, `h/04`, `h/05`, `h/07`, `h/08`, `h/09`, `o/01` |
| ADD | **8** | ``f/index``, `f/03`, `f/04`, `f/05`, `f/06`, `h/06`, `h/10`, `o/02` |
| SEQUENCE | **2** | `f/02`, `h/11` |
| DIAGRAM | **2** | `f/01`, `h/01` |
| NONE (justified) | **1** | ``o/index`` |
| *(secondary flag)* STALE-RISK | **9** | `f/02`, `f/03`, `f/05`, `f/06`, `h/06`, `h/07`, `h/08`, `h/09`, `o/01` |

### The five findings that dominate this audit

1. **Nothing is high-resolution.** Every image is captured at deviceScaleFactor 1 and
   displayed at its own natural pixel width (`246px` image rendered `:width: 246px`).
   On any HiDPI display all nine are visibly soft. The fix is mechanical and applies
   to all of them: capture at DPR 2, keep `:width:` at the CSS size.
2. **Two images encode facts the app no longer produces.**
   `analysis/marker-genes-expanded.png` and `analysis/gene-signature.png` show Ensembl
   IDs (`ENSG00000158488.16`, `ENSG00000017483.16`). No published sample dataset
   contains a single `ENSG` name any more — `suo`, `garcia`, `he`, `kanemaru` publish
   HGNC symbols and `pancreas` publishes MGI symbols (verified against
   `cellucid-datasets/exports/*/var_manifest.json`; ENSG count = 0 in all five).
   Both shots predate the gene-symbol switch.
3. **{doc}`/user_guide/web_app/h_analysis/06_analysis_mode_differential_expression_de` has no figure at all** — confirmed,
   the file contains zero `{figure}` directives. DE is the flagship statistical feature,
   it is the only explicit-run mode with a two-phase progress tracker, a volcano plot,
   a live threshold panel and an FDR-denominator summary table, and none of that is shown.
4. **Two existing captures are botched.** `analysis/detailed-categorical.png` has a
   stray Plotly modebar tooltip reading "Zoom" parked over the plot, and the bar chart
   underneath shows exactly one bar. `analysis/correlation-results.png` shows a
   *deselected, zero-cell* `Page 1 (0)` with the correlation computed on
   `Rest of Page 1` only, reporting `r = -0.005`, `R² = 0.000***` over a discrete
   `age` axis that renders as vertical stripes — a null result presented as the
   canonical example of correlation.
5. **The lifecycles are entirely missing.** Section F has one still of a *finished*
   selection and nothing of the act of selecting. Section H has six stills of six
   modes and nothing of the path selection → configure → run → read → export.

---

# Per-page records

## {doc}`/user_guide/web_app/f_highlighting_selection/index`

**explains** — Routes the reader through the section and states the section's
governing distinction: selection is the temporary candidate set, highlighting is the
persisted group inside a named page.

**has** — `none`.

**verdict** — **ADD**

**needs** — One orientation figure, because the page names three sub-areas
("Highlight mode", "Highlight pages", "Highlighted groups") the reader cannot yet
locate. Load the **pancreas** sample (`?dataset=pancreas`), colour by `cell_type`,
open the **Highlighting** accordion, confirm one small lasso so the group list is
non-empty, and crop to `#highlighted-cells-section .accordion-content` at DPR 2. No
cursor — this is a map, not an action. Three call-out labels overlaid (1 mode box,
2 page tabs, 3 group list) matching the numbered list in {doc}`/user_guide/web_app/f_highlighting_selection/03_highlight_ui`.
`_static/screenshots/highlighting_selection/panel-anatomy.png`, `:width: 246px`.

**notes** — The fast-path list is accurate against the source. No broken references.

---

## {doc}`/user_guide/web_app/f_highlighting_selection/01_highlight_mental_model`

**explains** — Establishes that filters control visibility while highlights control
membership, and that LOD/downsampling narrows what is *drawn* but never what a gesture
can *select* (except the seed `Alt+click`).

**has** — `none`.

**verdict** — **DIAGRAM** (primary), plus one two-shot proof pair.

**needs** —
- The `text`-fenced "mental picture" block at lines 41–55 should become a real
  diagram. It is a layering diagram (dataset → per-view visibility masks → global
  highlight pages → per-view renderer), it is the single most-linked concept in the
  section, and ASCII art in a `text` fence is unreadable at mobile widths. Mermaid
  `flowchart LR` is enough; the docs already build with `myst_parser`, so verify
  `sphinxcontrib-mermaid` is in `docs/conf.py` before committing to that syntax.
  No dataset, no capture.
- Then **two real screenshots** to make the one-sentence claim ("Filters control
  visibility. Highlights control membership.") observable. **suo** sample
  (561,947 cells — LOD is actually active there, which the page talks about),
  colour by `organ`, confirm an annotation selection on `Liver`. Shot A: crop to
  `#highlighted-groups` with the count reading `N cells highlighted`. Shot B:
  identical crop after adding a filter that hides Liver, count now reading
  `0 of N highlighted cells visible`. Cursor: none — the point is the text delta.
  `highlighting_selection/membership-vs-visibility-a.png` /
  `-b.png`, `:width: 246px` each, presented as a two-item `{grid}`.

**notes** — Line 108 says "Cellucid currently renders highlighted cells using a single
highlight style (halo/ring)" — verified true; groups carry no colour, pages do.
Line 104 gives the auto label as `Lasso (12,345 cells)`. The stored label is exactly
that (`lasso-selection.js:122`), but `highlight-summary-ui.js:151` re-appends the count,
so the list actually renders `Lasso (12,345 cells) (12,345)`. Any new screenshot of the
group list will show the doubled count and contradict this line — fix the doc text or
the app before capturing.

---

## {doc}`/user_guide/web_app/f_highlighting_selection/02_selection_tools_document_each_tool`

**explains** — Documents all four selection tools (annotation, lasso, proximity, KNN),
the `Alt` / `Shift+Alt` / `Ctrl+Alt` set-operation contract, and the 2D-vs-3D and
multi-view pitfalls of each.

**has** — `none`. This is the longest tool reference in the section and it is entirely
text.

**verdict** — **SEQUENCE** — see **LIFECYCLE A** below, which is anchored on this page.

**needs** — Steps 2–13 of Lifecycle A land here. Minimum viable subset if the full
sequence is too much: one *mid-gesture* shot per tool (four images), because the
gesture is the thing prose cannot convey — the dashed lasso polygon, the proximity
radius circle, the KNN degree growth, the annotation range label.

**notes — STALE-RISK** — Line 256 claims the UI indicates *"neighbor graph not
available"*. That string does not exist anywhere in `cellucid/assets/js`. The real
strings are `Neighbor graph ready (N edges)` (success toast, `main.js:3420`),
`Failed to load neighbor graph` (failure toast), and the mode help line *"Neighbor
graph status appears in notifications."* (`mode-copy.js:13`). Do not caption a
screenshot with the invented string.

Second, and more important: the page's "Practical fix — load connectivity/edges for
your dataset (often via a 'Show edges' toggle)" understates the real mechanism. In
`highlight-renderer.js` (~line 3300) the **first `Alt+click` in KNN mode itself
triggers the edge load and returns without selecting**. The user's first KNN gesture
always appears to do nothing, then a toast arrives. That behaviour is undocumented and
is exactly what a screenshot should capture (Lifecycle A step 9).

---

## {doc}`/user_guide/web_app/f_highlighting_selection/03_highlight_ui`

**explains** — The UI map: what each button in the Highlighting accordion does, how
pages are created/renamed/recoloured/deleted/combined, how groups are
enabled/removed/cleared, and how to read the "visible vs total" count.

**has** — `none`. This is the page a wet-lab reader lands on to find a control, and it
contains no picture of any control.

**verdict** — **ADD** (5 figures), with a **STALE-RISK** flag.

**needs** — All on **pancreas**, DPR 2, sidebar width 246 CSS px:
1. *Step controls exist* — after finishing one lasso, crop to `.highlight-mode-box`
   (contains the four `.highlight-mode-btn`s, `#highlight-mode-description` showing
   `Step 1: 412 cells`, and `#lasso-step-controls` with Confirm / ↩ / ↪ / Cancel).
   Cursor composited on `#lasso-confirm-btn`.
   `highlighting_selection/step-controls.png`, `:width: 246px`.
2. *Page tabs, mid-rename* — create a second page, double-click the name so
   `.highlight-page-tab-input` is live with `Beta cells` typed. Crop
   `#highlight-pages-container`. Cursor on the tab name.
   `highlighting_selection/page-rename.png`, `:width: 246px`.
3. *Recolour* — the `input[type=color]` swatch open on a page tab. Crop
   `#highlight-pages-container` plus the native picker.
   `highlighting_selection/page-recolor.png`, `:width: 246px`.
4. *Combine menu* — click `.highlight-page-tab-combine` (title "Combine with another
   page"); the floating `role="menu"` shows `∩ Intersection with <name>` and
   `∪ Union with <name>` as separate `.page-combine-option` buttons. Crop the tab strip
   plus the menu. Cursor on the `∩` option.
   `highlighting_selection/page-combine-menu.png`, `:width: 246px`.
5. *Create Categorical* — expand `#cat-builder-accordion-item` with two pages dropped
   into `#cat-builder-dropzone` and the `#conflict-section` visible (overlap strategy
   radios). Crop `#category-builder-container`. Cursor on the
   `overlap-strategy=intersections` radio.
   `highlighting_selection/create-categorical.png`, `:width: 246px`.

**notes — STALE-RISK, two defects:**
- Lines 197–204: *"Some builds include a 'Category builder' UI … If you don't see it,
  it may be disabled or hidden in your build."* This is false and self-contradicted by
  the section's own existing screenshot, which shows the panel. The container
  `#category-builder-container` is hard-coded in `index.html` and
  `field-selector.js:1345-1359` constructs it unconditionally (it *throws* if already
  initialised — there is no disabled path). Its on-screen label is **"Create
  Categorical"**, not "Category builder". Rewrite as a documented feature.
- Lines 108–115: combining pages is described as drag-and-drop only. There is also an
  explicit `.highlight-page-tab-combine` button on every tab. The button is the
  discoverable and keyboard-reachable path and is undocumented.

---

## {doc}`/user_guide/web_app/f_highlighting_selection/04_selection_synchronization`

**explains** — The scope table: what is global (pages, groups, candidate set) vs
per-view (visibility), what switching pages does and does not do, and the Python/Jupyter
event direction (selection fires on **confirm**, not on drag).

**has** — `none`.

**verdict** — **ADD** (2 figures).

**needs** —
1. *The same group, two panels* — **suo**, colour by `organ`, confirm one annotation
   group on `Liver`. Enter grid multi-view with two panels; apply a filter in panel B
   that excludes `Liver`. Full window 1440×1000 at DPR 2, cropped to the canvas region
   only (exclude the sidebar) so both panels are side by side; panel A shows haloed
   cells, panel B shows none. Overlay two thin labels "live view" / "snapshot ·
   organ ≠ Liver". No cursor.
   `highlighting_selection/multiview-highlight-scope.png`, `:width: 720px`.
2. *Confirm is the event boundary* — a **Jupyter** screenshot, which the brief
   explicitly wants for the data-loader/bridge paths. JupyterLab notebook with
   `@viewer.on_selection` registered in cell 1, the embedded viewer in cell 2, and cell
   3's output showing the callback fired **once** with a cell count after Confirm was
   clicked (not once per drag frame). Crop to the notebook column including the viewer
   iframe. Cursor on the viewer's Confirm button.
   `jupyter/on-selection-fires-on-confirm.png`, `:width: 860px`.
   (Reuses the existing `_static/screenshots/jupyter/` topic; that folder already holds
   `pancreas-notebook-embed.png` and `pancreas-debug-connection.png`.)

**notes** — Scope table verified correct against `selection-state.js` and
`mode-ui.js`. `viewer.debug_connection()` reference is live.

---

## {doc}`/user_guide/web_app/f_highlighting_selection/05_edge_cases_highlighting`

**explains** — The "weird but expected" checklist: zero visible cells, field changes,
cells in multiple groups, page switching mid-selection, million-cell groups, low-contrast
highlights, 3D lasso projection, multi-view, missing KNN graph, dataset changes.

**has** — `none`.

**verdict** — **ADD** (3 figures, chosen for the cases where the *visual* is the whole
point).

**needs** —
1. *3D projection ambiguity, before/after* — the single most-explained and
   least-shown idea in the section. **pancreas**, 3D mode, colour by `cell_type`.
   Shot A: lasso drawn around what looks like one cluster, `Step 1: 890 cells`
   readout visible in an inset. Shot B: camera rotated ~90°, same confirmed set now
   visibly spanning two separated lobes. Crop the canvas only, 720×540 each,
   presented as an A/B pair. Cursor on the lasso head in A only.
   `highlighting_selection/lasso-3d-projection-a.png` / `-b.png`, `:width: 360px` each.
2. *Zero visible cells* — crop to `#highlighted-groups` with the count reading
   `0 of 55,096 highlighted cells visible` and the group checkbox still ticked. Cursor
   on the checkbox to make it obvious the group is enabled.
   `highlighting_selection/zero-visible-highlighted.png`, `:width: 246px`.
3. *KNN with no neighbour graph* — this **cannot be reproduced on any sample dataset**:
   `suo`, `pancreas`, `garcia`, `he`, `kanemaru` all ship a `connectivity/` directory.
   Use the local fixture `cellucid/tests/browser/fixtures/exports/current-ui-prepared`
   (120 cells, no `connectivity/`). Select KNN drag, `Alt+click` a cell, capture the
   resulting notification toast plus `#highlight-mode-description`. Crop to the
   sidebar column plus the toast.
   `highlighting_selection/knn-no-graph.png`, `:width: 400px`.

**notes — STALE-RISK** — Line 172 repeats the non-existent *"Neighbor graph not
available"* string. Capture the real toast text first, then write the caption from it.

---

## {doc}`/user_guide/web_app/f_highlighting_selection/06_troubleshooting_highlighting`

**explains** — Symptom → cause → confirm → fix → prevention, and — uniquely — a table
of the five literal sentences the app writes under the Highlight mode buttons.

**has** — `none`.

**verdict** — **ADD** (1 composite figure — the highest value-per-pixel image in the
whole audit).

**needs** — The "What the app tells you" table (lines 54–60) quotes five exact strings
that all render into the same element, `#highlight-mode-description`. Capture that one
element five times and stack them into a single labelled composite:

| Row | How to produce it | Dataset |
| --- | --- | --- |
| `That click did not land on a cell…` | Proximity drag, zoom far out with LOD on, `Alt+click` empty space, no candidate set | suo |
| `That gesture ended outside the view…` | Lasso, `Alt+drag`, release the mouse over the sidebar | pancreas |
| `That cell has no value on the active field…` | Annotation based on a field with a missing-code sentinel, `Alt+click` a sentinel cell | he |
| `That gesture did not change the selection.` | Lasso the identical polygon twice with `Alt` | pancreas |
| `Annotation selection needs an active field…` | Annotation based with colour-by set to None | pancreas |

Crop each to `#highlight-mode-description` (plus the `Step N:` line above it where
present, since the doc says the readout distinguishes the two cases). Assemble into one
vertical strip with the trigger written beside each row. No cursor. DPR 2.
`highlighting_selection/selection-notices.png`, `:width: 640px`.

**notes — STALE-RISK** — Line 291 quotes *"Neighbor graph not available"* again (third
occurrence in the section). Line 479 tells bug reporters to "use ``07_screenshots`` for
capture specs", but `07_screenshots.md` currently contains one figure and zero capture
specs — either write specs there or drop the promise.

---

## {doc}`/user_guide/web_app/f_highlighting_selection/07_screenshots`

**explains** — Nine lines. One figure of a confirmed highlight page, captioned as the
group that later analysis and session workflows consume.

**has** — `highlighting_selection/highlighting-selected-page.png`, **line 3**,
`:width: 246px`. Native size 246×363 → **1× DPR, soft on HiDPI**.
**Reused in two other repositories' pages**:
`user_guide/python_package/e_jupyter_hooks/02_quickstart_minimal_roundtrip.md:108` and
`user_guide/python_package/f_notebooks_tutorials/32_session_persistence_and_restoring_analysis_artifacts.md:247`.
Any replacement must be checked against those two captions as well.

**verdict** — **REPLACE**, and expand this page into the home of **Lifecycle A**.

**needs** —
- Recapture the existing shot at DPR 2 from the **suo** sample (the current one is
  suo — 55,096 cells, `Annotation` group — and is still factually correct). Same crop:
  `#highlighted-cells-section .accordion-content`. Keep the filename so the two
  Python-package pages keep working.
- Then host the ordered sequence. Right now the page's title promises "Verified
  highlighting capture" (singular) and the troubleshooting page points bug reporters
  here for "capture specs" that do not exist.

**notes** — The current image already disproves {doc}`/user_guide/web_app/f_highlighting_selection/03_highlight_ui`'s "some builds
include a Category builder": **Create Categorical** is visible at the bottom of it.
The group label in it reads `Annotation (55,096 cells) (5…` — the doubled count
described under `01`.

---

## {doc}`/user_guide/web_app/h_analysis/index`

**explains** — Section router: a "what are you trying to answer" table mapping six
questions onto the six analysis modes, plus the recommended reading order.

**has** — `analysis/analysis-panel-tabs.png`, **line 52**, `:width: 246px`.
Native 246×262 → 1× DPR. **Reused twice more**: `02_analysis_ui_overview.md:165` and
`11_screenshots.md:5` — the same 262-pixel-tall image appears three times in one section.

**verdict** — **REPLACE**

**needs** — Recapture at DPR 2, but fix the real problem first: at the 246 px sidebar
width every mode description truncates (`Full control over opt…`, `Explore variabl…`,
`Fi…`, `Compute si…`, `Discover mar…`). The one figure whose job is to name the six
modes cannot show what four of them do. Capture the **floating/undocked** analysis
panel instead — click `.analysis-accordion-copy-btn`, which produces a
`.accordion-section.analysis-window-panel` that can be widened — or widen the sidebar
before capture. Dataset: **pancreas**. Crop to `.analysis-accordion`. No cursor.
`analysis/analysis-panel-modes.png`, `:width: 420px`. Retire
`analysis-panel-tabs.png` and update all three call sites in the same change.

**notes** — Mode names in the caption match `comparison-module.js` exactly (Quick,
Detailed, Correlation, Differential Expression, Gene Signature, Marker Genes). The
`ANALYSIS_TYPES` registry names them differently again (`Quick Insights`,
`Correlation Explorer`, `Gene Signature Score`) but those strings never reach the DOM.

---

## {doc}`/user_guide/web_app/h_analysis/01_analysis_mental_model`

**explains** — The three nouns (dataset / page / variable), membership vs visibility,
`Rest of <page>` complements, the caching layers, and the statistical scope limits.

**has** — `none`.

**verdict** — **DIAGRAM** (primary) + 1 supporting capture.

**needs** —
- A diagram for the three nouns and the scope boundary: dataset (all indices) → pages
  (union of enabled groups) → variables (obs categorical / obs continuous / gene) →
  the mode that consumes them, with the filter/visibility layer drawn *outside* the
  analysis boundary. This is the page's entire thesis and it currently has no visual.
  Mermaid, no capture.
- One capture for `Rest of <page>`, which is pure vocabulary until seen. **pancreas**,
  one page named `Beta`. Open **Detailed**, click the **Compare pages:** selector; the
  derived chip `Rest of Beta` is rendered alongside the base page with its own cell
  count (`page-derivation-utils.js:60` formats it as exactly `Rest of ${base}`). Crop
  to `.page-selector` / the `Compare pages:` block. Cursor on the `Rest of Beta` chip.
  `analysis/rest-of-page-selector.png`, `:width: 246px`.

**notes** — Statistical claims spot-checked and accurate (chi-squared with expected ≥ 5,
automatic Fisher for sparse 2×2, N/A for larger sparse tables, Welch + Mann–Whitney with
exact/asymptotic switch, tie-corrected Kruskal–Wallis). No broken references.

---

## {doc}`/user_guide/web_app/h_analysis/02_analysis_ui_overview`

**explains** — Where Analysis lives, accordion behaviour (one mode open at a time),
auto-run vs explicit-run modes, the Compare pages selector, preview vs expanded modal,
and what "Copy" (floating windows) does.

**has** — `analysis/analysis-panel-tabs.png`, **line 165**, `:width: 246px`. Third
appearance of the same file (see ``h/index``).

**verdict** — **REPLACE** (the reused tabs shot) **+ ADD** (3 figures for the three
behaviours the page describes and nothing shows).

**needs** —
1. *Auto-run vs Run* — an A/B pair cropped to `#analysis-panel-correlation` (result
   appears with no button) and `#analysis-panel-differential` (a
   `Run Differential Expression` button sits between the form and any result). Cursor
   on the Run button in B. **pancreas**.
   `analysis/autorun-vs-run-a.png` / `-b.png`, `:width: 246px` each.
2. *Preview → modal* — the `⤢ Expand` button in a result's action row, cursor on it,
   cropped to the bottom of `#analysis-panel-detailed`. The page states the button is
   the *only* way in and the preview is not clickable — worth showing.
   `analysis/expand-button.png`, `:width: 246px`.
3. *Copy → floating window* — after clicking `.analysis-accordion-copy-btn` on
   Correlation, a `.accordion-section.analysis-window-panel` floats over the canvas
   with its own header, copy and close buttons. Full window 1440×1000 at DPR 2, cropped
   to the region containing both the sidebar accordion and the floating window so the
   relationship is legible. Cursor on the floating window's title bar.
   `analysis/floating-analysis-window.png`, `:width: 860px`.

**notes** — "Settings are copied. Results are not copied." verified against
`analysis-window-manager.js`. The safety cap on open windows is real. Correlation's
`Color by:` and `Correlation method:` labels match `correlation-analysis-ui.js`.

---

## {doc}`/user_guide/web_app/h_analysis/03_analysis_mode_quick_insights`

**explains** — Quick mode: Dynamic (follow active page) vs Manual page selection, the
Composition section (top-5 categories as a stacked bar) and the Statistics section
(mean/median/std per continuous field).

**has** — `analysis/quick-insights.png`, **line 182**, `:width: 246px`.
Native 246×546 → 1× DPR. Also at `11_screenshots.md:15`.

**verdict** — **REPLACE**

**needs** — The content of the current shot is correct (suo, `Page 1: 55,096 cells`,
composition on `sample_ID` and `organ`, statistics on `age` and `n_genes`, the
`▶ Page Selection (Dynamic)` collapsible). Two problems to fix on recapture:
- 1× DPR.
- The composition bars render **greyscale**, which makes the "segment widths are
  percent of cells" claim hard to read and looks broken beside the coloured page chips.
  Verify against `quick-insights-ui.js` whether that is the intended palette before
  captioning it.

Recapture on **suo**, `?dataset=suo`, colour by `organ`, one annotation-selected page
renamed `Liver`, crop `#analysis-panel-simple`, DPR 2, no cursor.
Then **ADD** a second figure for the Dynamic/Manual distinction, which is the page's
main failure mode: the `Page Selection` collapsible expanded with the **Manual** radio
chosen and two pages ticked. Cursor on the Manual control.
`analysis/quick-page-selection-manual.png`, `:width: 246px`.

**notes** — No stale claims found. The troubleshooting entries name real UI ("open
Quick's **Page Selection** collapsible at the bottom and check the mode indicator").

---

## {doc}`/user_guide/web_app/h_analysis/04_analysis_mode_detailed_analysis`

**explains** — Detailed mode: one variable across many pages, plot-type choice per
variable kind, the summary table, and the statistical annotations (chi-squared/Fisher;
Welch + Mann–Whitney for 2 pages; ANOVA + Kruskal–Wallis for 3+).

**has** — `analysis/detailed-categorical.png` (`:width: 488px`, native 488×1014) and
`analysis/detailed-expanded.png` (`:width: 1440px`, native 1440×1105), the
`analysis-detailed-categorical-preview` and `analysis-detailed-expanded` scenarios, both
on **pancreas** at DPR 2.

**verdict** — done. Nothing outstanding.

**needs** —
- ~~A **Plotly modebar tooltip reading "Zoom" is frozen over the plot**.~~ **Done** —
  the pointer sits on the `Plot Type:` select, and the modebar is not hovered.
- ~~The bar plot shows **exactly one bar**.~~ **Done** — `cell_type` on **pancreas**,
  eight bars, `Page 1` and `Rest of Page 1` both selected and both in the legend.
- ~~**ADD** the expanded view.~~ **Done** — `.analysis-modal-content`, cursor on the
  `Export: CSV` button.

Two things the recapture corrected in the instruction above. **Percentages** is left
off: the figure is captioned as counts and the summary table beside it prints counts, so
turning it on would put a percentage axis under a counts caption. And "8 categories, all
legible" is not achievable and is not what the page now claims — a categorical axis is
thinned to the labels its box can show, which is one name in the sidebar and five in the
expanded view; `cellucid c177e77ce` made that a clean drop instead of a smear.

**notes** — The test descriptions match `statistical-tests.js`. "Export … per-cell
values … `detailed-analysis-data.csv` … columns `page`, `cell_index`, `<variable>`"
matches {doc}`/user_guide/web_app/h_analysis/09_exporting_analysis_results`; both should be re-verified after any
recapture that shows the export row.

---

## {doc}`/user_guide/web_app/h_analysis/05_analysis_mode_correlation_analysis`

**explains** — Correlation mode: per-page Pearson/Spearman on paired finite values,
what is reported (`r`, `r²`, `p`, `n`, slope, intercept), and that plotting downsamples
to ~50k points while statistics use the full pairing.

**has** — `analysis/correlation-results.png`, **line 191**, `:width: 224px`.
Native 224×604 → 1× DPR. Also at `11_screenshots.md:35`.

**verdict** — **REPLACE** — the current image demonstrates a degenerate case.

**needs** — Three defects in the current file:
- **`Page 1 (0)` is greyed out and deselected** — the user's own page is empty, so the
  only thing correlated is `Rest of Page 1`. The caption ("plots two continuous
  observation fields for the selected cells") describes a state that is not on screen.
- The result reads **`r = -0.005`, `R² = 0.000***`** — a null correlation carrying
  three significance stars because `n = 46,829`. This is precisely the
  "significance with tiny effect size" trap the *Detailed* page warns about, presented
  here as the canonical example.
- X is `age`, which in `suo` is a small integer range (4–17). The scatter renders as
  **vertical stripes**, which is not what a correlation plot should teach.

Recapture on **pancreas**, whose continuous obs are `S_score` and `G2M_score` — the
only genuinely continuous, genuinely structured pair among all five samples. Set
`X Axis Variable: Continuous obs → S_score`, `Y: Continuous obs → G2M_score`,
`Compare pages:` one non-empty page, `Color by: cell_type`,
`Correlation method: Pearson (linear)`. Crop `#analysis-panel-correlation`. Cursor on
the `Correlation method:` select. DPR 2.
`analysis/correlation-results.png` (same name), `:width: 224px`.

Optionally **ADD** the Spearman counterpart of the same pair as a B-shot to make the
"try Spearman if you suspect non-linearity" advice actionable.

**notes** — Statistics list verified. The `p-value` derivation (two-sided Student *t*
on `n-2` df) matches the implementation.

---

## {doc}`/user_guide/web_app/h_analysis/06_analysis_mode_differential_expression_de`

**explains** — DE between two highlight pages: Wilcoxon (default) or Welch, the
`log2FC = log2((meanA + 0.01)/(meanB + 0.01))` formula, Benjamini–Hochberg over
*tested* genes only, the volcano plot, the threshold controls, and the summary rows
including the FDR denominator.

**has** — **`none`.** Confirmed: zero `{figure}` directives in a 287-line page. **This
is the single largest gap in the section**, on the flagship statistical feature.

**verdict** — **ADD** (5 figures), with a **STALE-RISK** flag.

**needs** — All on **pancreas** (3,696 cells / 3,753 genes — DE completes in seconds,
and `Ins1`/`Ins2`/`Iapp`/`Gcg` give a result a wet-lab reader recognises). Two named
pages: `Beta` (annotation-selected on `cell_type = Beta`) and `Alpha`.
1. *The configured form* — `#analysis-panel-differential` before running: the
   `Select pages to compare:` label with `Page A:` and `Page B:` `<select>`s. Open
   `Page B:` so the two `<optgroup>`s are visible — **Pages** and **Wildcards**
   (the latter holds `Rest of Beta`). The page's prose already claims those exact
   group labels; this is the shot that proves it. Cursor on the open `Page B:` select.
   `analysis/de-page-selection.png`, `:width: 246px`.
2. *Method + performance* — same panel scrolled to `Statistical method:` (Wilcoxon /
   t-test, each with its description line) and `Performance Settings` **expanded**
   (batch size, memory budget, network parallelism, compute parallelism). The page
   devotes a whole section to these four controls and never shows them. Cursor on the
   Performance Settings disclosure.
   `analysis/de-performance-settings.png`, `:width: 246px`.
3. *Running* — captured mid-run: the run button reading `Running...` and the progress
   tracker on phase 1 of 2, `Loading & Computing` (phase names are literal, from
   `de-analysis-ui.js:343`). Crop `#analysis-panel-differential`. Cursor on the
   disabled run button. Pair it with a phase-2 `Multiple Testing Correction` variant if
   the run is long enough to catch — on pancreas it may not be, in which case use
   **suo** with `Rest of <page>` and accept the longer runtime.
   `analysis/de-running.png`, `:width: 246px`.
4. *The volcano, expanded* — `.analysis-modal-content` after `⤢ Expand`:
   `.analysis-modal-plot` with the volcano (up/down/not-significant colouring and
   labelled top-N genes), the `PLOT OPTIONS` column with the p-value threshold, the
   |log2FC| slider and the **Use FDR-corrected p-values** checkbox, the `Export:`
   row (PNG / SVG / CSV), the `SUMMARY STATISTICS` table and the
   `Top Differentially Expressed Genes` table with its Top 5/10/20/100/All select.
   This is the flagship image of section H. **No cursor** — the |log2FC| slider asked
   for here is a *read-only state* in this frame, and the frame that does move it
   (`volcano-threshold-b.png`) crops to `.analysis-modal-plot`, which the slider is
   outside of. A pointer drawn on the slider would either be outside the crop or point
   at a control the caption does not say was touched.
   `analysis/de-volcano-expanded.png`.
5. *The FDR denominator* — a tight crop of just the `.de-modal-stats` table, at a
   readable width, showing the five rows the page documents in prose:
   `Genes tested (FDR denominator)`, `Not tested (< 10 cells with a value)`,
   `Significant (…)`, `Upregulated`, `Downregulated`. Choose a page pair small enough
   that *Not tested* is non-zero, so the page's central argument ("quote genes tested,
   not panel size") is visible rather than asserted. No cursor.
   `analysis/de-fdr-denominator.png`, `:width: 420px`.

**notes — STALE-RISK, two:**
- Line 237: *Symptom: "Need at least 1 page for comparison"*. **That string does not
  exist in the source.** The two real strings are
  `Select two non-empty cell groups to run differential expression`
  (`de-analysis-ui.js:228`) and
  `Differential expression requires two non-empty cell groups.`
  (`selectors.js:456`, shown when fewer than two non-empty options exist). Fix the
  symptom heading before anyone screenshots it.
- Line 171: the summary row is quoted as `Significant (FDR < 0.05)`. The rendered
  label is dynamic and uses `≤` plus the fold-change threshold:
  `Significant (FDR ≤ 0.05, |log₂FC| ≥ 1)`. A screenshot of figure 5 will contradict
  the prose as written.
- Also worth documenting: `createPageComparisonSelector` throws rather than degrades if
  a zero-cell page is pre-selected, and swaps to a `Comparing: A vs B` badge display
  when exactly two base pages exist and wildcards are off. Neither branch is described.

---

## {doc}`/user_guide/web_app/h_analysis/07_analysis_mode_gene_signature`

**explains** — Gene signature scoring: comma-separated gene input, mean/sum/median
scoring, z-score and min-max normalisation computed across all selected pages combined,
and the `gene_signature_scores.csv` export (page, score — no cell index).

**has** — `analysis/gene-signature.png`, **line 176**, `:width: 224px`.
Native 224×623 → 1× DPR. Also at `11_screenshots.md:45`.

**verdict** — **REPLACE**, flagged **STALE-RISK**.

**needs** — The current image is **stale against the shipped datasets**. Its
`SIGNATURE GENES` box contains
`ENSG00000017483.16, ENSG00000019169.11, ENSG00000021826.18, ENSG00000023171.20`.
No sample dataset publishes Ensembl IDs any more (0 `ENSG` names across all five
`var_manifest.json` files). The page's own prose uses `CD3E, CD4, IL7R` and `MS4A1`.
A reader who copies what is on screen gets zero matches — and the page's headline
warning is precisely "matching is exact, no alias mapping".

It also shows a **single violin** (only `Page 1` selected; `Rest of…` is unticked), so
"compare across pages" is not demonstrated.

Recapture on **pancreas**: `Signature Genes: Ins1, Ins2, Iapp` (all three present in
`pancreas/var_manifest.json`), two pages `Beta` and `Rest of Beta` both selected,
`Scoring Method: Mean expression`, `Normalization: None`,
`Visualization: Violin plot` — which yields two clearly separated violins and a
biologically obvious result. Crop `#analysis-panel-signature`. Cursor on the
`Normalization:` select. DPR 2.
`analysis/gene-signature.png` (same name), `:width: 224px`.

**notes** — Line 165 suggests `MS4A1` as a probe gene. `MS4A1` is **not** in
`pancreas` (mouse) and should be `Ms4a1` there; it is present in `suo`/`garcia`. If the
new screenshot uses pancreas, align the prose example with it. The cross-reference to
``d_fields_coloring_legends/07_genes_in_the_built_in_samples`` is the right place to
resolve this and should be checked in the same pass.

---

## {doc}`/user_guide/web_app/h_analysis/08_analysis_mode_genes_panel`

**explains** — Marker Genes: one-vs-rest discovery across every category of a
categorical obs field, the three modes (Ranked Genes / Clustered / Custom Genes),
per-group BH denominators, the "genes tested / not tested" line, and the two CSV shapes.

**has** — `analysis/marker-genes-ranked.png` (`:width: 488px`, native 488×1584) and
`analysis/marker-genes-clustered.png` (`:width: 1160px`, native 1160×746), the
`analysis-marker-genes-ranked` and `analysis-marker-genes-clustered` scenarios in
`docs/_tooling/screenshots/scenarios.mjs`, both on **pancreas** at DPR 2.

**verdict** — **ADD** 1. The **REPLACE** is done: see needs below.

**needs** —
- ~~**Caption/image mismatch.**~~ **Done**, and the instruction that replaced it was
  itself wrong. It asked for "the ranked marker list visible" inside
  `#analysis-panel-genesPanel`; there is no such list in the sidebar. All three modes
  draw the *same heatmap* there — `_showResult` in
  `assets/js/app/analysis/ui/analysis-types/genes-panel-ui.js` imports the gene heatmap
  and renders it whatever the mode is — and the ranked per-group table exists only in
  the expanded view. `08:141-148` had this right all along. The figure now shows the
  Ranked Genes form with the heatmap it really draws, and the caption says so.
- ~~The 224 px preview heatmap is **illegible**.~~ **Fixed in the app**, not in the
  capture: `cellucid c177e77ce` thins a categorical axis to the labels its box can show
  and drops the rest whole, so the gene axis is now a handful of separated names rather
  than a smear. The colour-bar tick labels below the `Z-score` bar still overprint — see
  the standing defect noted at the end of this record.
- ~~Group-by is `organ` on **suo**.~~ **Done** — both figures are `cell_type` on
  **pancreas**, 8 groups.

Widths: the `:width: 246px` / `:width: 720px` above were guesses made before the size
rule existed. The rule now fixes them — 2 × the crop's CSS width, capped at 1440 — and
the tool prints the value.

Still **ADD**: a tight crop of the *genes tested* line above the ranked table, on a
dataset where the "not tested" half is non-zero, so the page's central rule is
demonstrated rather than asserted:
`analysis/marker-genes-tested-line.png`.

**notes** — `Use cached results`, `CLUSTERING`, `PERFORMANCE SETTINGS`,
`DISCOVER MARKERS` all verified as real controls. The two "no markers" messages quoted
at lines 165–170 are exact.

**standing defects, reported not fixed here** — Two things the recapture made visible
and neither is a capture problem:

- ~~The group dropdown above the ranked table lists **`category-code:0`,
  `category-code:1`, …** rather than the group names.~~ **Fixed in the app**, not in
  the capture: the picker, the heatmap axis and the ranked CSV's `group` column now
  all render the category name through the one exported rule, `encodeGroupName` in
  `assets/js/app/analysis/genes-panel/expression-matrix-builder.js`; the handle
  survives only as the `<option>`'s `value`, because it still keys `markers.groups`,
  the matrix columns, the hover lookup and the saved `modalSelectedGroupId`.
  **Still ADD:** `analysis/marker-genes-expanded.png` predates the fix and shows
  `category-code:0` in the picker beside `Ductal` on the heatmap row — it must be
  re-shot from the `analysis-marker-genes-expanded` scenario. It is embedded on
  {doc}`/user_guide/web_app/h_analysis/09_exporting_analysis_results`, which is the only page still using it.
- The horizontal colour bar under the sidebar heatmap draws its tick labels on top of
  one another. `tests/browser/analysis-plot-label-collisions.spec.mjs:144` measures
  `.xaxislayer-above .xtick text` only, so colour-bar ticks are outside its acceptance.

---

## {doc}`/user_guide/web_app/h_analysis/09_exporting_analysis_results`

**explains** — What each mode exports, the exact CSV schemas and filenames, that
exports live in the expanded modal only, and that most CSVs carry no provenance.

**has** — `analysis/marker-genes-expanded.png`, `:width: 1440px`, native 1440×1105, a
860 CSS px crop of `.analysis-modal-content` captured at DPR 2 by the
`analysis-marker-genes-expanded` scenario. This is now its only use in the section —
`11_screenshots.md` no longer embeds it.

**verdict** — **RE-SHOOT.** The image now predates the group-picker fix recorded under
{doc}`/user_guide/web_app/h_analysis/08_analysis_mode_genes_panel` above: it shows `category-code:0` in the picker where
the app draws the category name. Everything else on the page is settled.

**needs** —
- ~~**The caption is factually wrong** ("PDF, PNG, and CSV export actions").~~ **Done**
  on both pages; the toolbar in `components/export.js` emits exactly **PNG**, **SVG**
  and **CSV**, and both captions now say so.
- ~~The image is **stale**: `ENSG…` gene names.~~ **Done** — recaptured on **pancreas**,
  whose var manifest publishes MGI symbols (`Aplp1`, `Ins2`, `Gcg`).
- ~~The crop cuts the `PLOT OPTIONS` column mid-control.~~ **Unchanged and correct as
  is.** That column scrolls; a crop of `.analysis-modal-content` frames the whole modal,
  and the column is genuinely taller than it. Sizing the modal until the column fits
  would photograph a window the app never opens at.

~~Then **ADD** the two shots that close the export loop.~~ **Both withdrawn.**
- The browser download shelf renders a real local path and cannot be sanitised.
- The terminal capture is forbidden outright: see *Never photograph a terminal.
  Transcribe it.* in ``00_capture_tooling_and_conventions``. The column list
  (`gene, meanA, meanB, log2FoldChange, pValue, adjustedPValue`) belongs in a verbatim
  ```` ```text ```` block if it is ever shown, not in an image. This record was missed
  by the sweep that withdrew the other four commissioned terminal captures.

**notes** — All five CSV schemas were spot-checked against the serialisers and are
accurate. The `⤢ Expand` claim ("the button is the only way in — the preview itself is
not clickable") matches `createExpandButton`. The Marker Genes half of that check has
since been re-verified against `_buildModalCSVExport`, which dispatches on
`result.metadata.mode` and refuses an unrecognised one.

---

## {doc}`/user_guide/web_app/h_analysis/10_troubleshooting_analysis`

**explains** — Symptom → cause → confirm → fix for the eight most common analysis
failures (empty analysis, slow DE, wrong volcano, NaN correlation, windows not
restoring, tiny marker groups, empty signature, filters not affecting results).

**has** — `none`.

**verdict** — **ADD** (2 figures — the two symptoms whose *appearance* is the
diagnostic).

**needs** —
1. *No pages* — the state a first-time user actually hits. **pancreas** freshly
   loaded, no highlight page created, **Analysis → Quick** open. Crop
   `#analysis-panel-simple` showing the empty/no-pages state. Pair it with the
   Highlighting panel above showing the empty group list, so the causal link is on one
   image. No cursor.
   `analysis/analysis-empty-no-pages.png`, `:width: 246px`.
2. *Volcano hidden by thresholds* — an A/B pair of `.analysis-modal-plot` from the same
   DE result: A at `|log2FC| ≥ 3` with almost no points coloured, B at `|log2FC| ≥ 0.5`
   with the same run showing a populated volcano. This is the page's most-reported
   symptom and it is purely visual. Cursor on the threshold slider in A.
   `analysis/volcano-threshold-a.png` / `-b.png`, `:width: 420px` each.

**notes** — The Marker Genes error quoted at line 266,
`Group "X" has only N cells. Minimum required: 10.`, **is accurate** — it is the
template at `analysis/genes-panel/constants.js:172`,
`Group "{name}" has only {n} cells. Minimum required: {min}.`. This is the one quoted
error string in section H that survives verification; the DE symptom heading in `06`
does not. Safe to use as a screenshot caption.

---

## {doc}`/user_guide/web_app/h_analysis/11_screenshots`

**explains** — A gallery of "verified analysis captures", one per mode.

**has** — Seven figures, all duplicates of images already embedded on their own mode
pages: `analysis-panel-tabs.png` (**line 5**, 246px), `quick-insights.png`
(**line 15**, 246px), `detailed-categorical.png` (**line 25**, 224px),
`correlation-results.png` (**line 35**, 224px), `gene-signature.png` (**line 45**,
224px), `marker-genes.png` (**line 55**, 224px), `marker-genes-expanded.png`
(**line 62**, 860px). **Zero unique images.** All seven are 1× DPR.

**verdict** — **SEQUENCE** — rebuild this page as the home of **LIFECYCLE B**.

**needs** — In its current form the page is pure duplication: every image and nearly
every caption is copied from `03`–`09`, so it doubles the maintenance surface without
adding information, and it carries the wrong "PDF" caption a second time. Replace the
gallery with the ordered lifecycle below. Keep a short "verified capture conditions"
preamble naming the dataset, viewport and DPR used, which is what the troubleshooting
pages already promise readers will find here.

**notes** — Line 66 repeats the false "PDF, PNG, and CSV export actions" caption.
Line 59 ("Marker Genes ranks candidate genes for the selected group") repeats the
Ranked-vs-Clustered mismatch from `08`.

---

## {doc}`/user_guide/web_app/o_accessibility_privacy_security/index`

**explains** — Two-page router with a fast-path table (sensitive data / sharing
sessions / colourblind figures) and two grid cards.

**has** — `none`.

**verdict** — **NONE** — justified.

**needs** — Nothing. The page contains no UI claim, no procedure and no state; it is a
table of three rows and two `{grid-item-card}`s pointing at the two real pages. A
screenshot here would decorate, not explain. The one thing it *does* assert visually —
the two `{octicon}` glyphs — is theme-rendered and needs no capture.

**notes** — All five cross-references resolve. No stale prose.

---

## {doc}`/user_guide/web_app/o_accessibility_privacy_security/01_accessibility`

**explains** — Colour/contrast defaults for a WebGL point cloud, the colourblind
simulation in Figure Export, the `i` help-dialog keyboard contract, the global
shortcut table, what keyboard cannot do, and motion-sensitivity guidance.

**has** — `web_app/keyboard-shortcuts.png`, `:width: 516px`. Native 516×1154, a
258 CSS px crop of `#shortcuts-section` captured at DPR 2 by the
`accessibility-panel-keyboard-shortcuts` scenario. Not reused elsewhere in the
docs.

**verdict** — **ADD** 3. The **REPLACE** is done: see needs 1 below.

### Is accessibility screenshot-shaped? Partly — here is the honest split.

**Warrants an image (visual state the reader must recognise):**
- The **focus ring**. It is a single global rule —
  `:focus-visible { outline: 2px solid var(--color-border-focus); outline-offset: 2px }`
  in `assets/css/base/_accessibility.css:15` — so one image teaches it everywhere. The
  page currently never shows what a focused control looks like.
- The **`i` help dialog**, whose open/closed/focus-return behaviour is described in
  eleven lines of prose (lines 137–152) and is trivially demonstrable in one shot.
- The **colourblind simulation**, whose entire value proposition is a before/after of
  the same plot. Prose cannot do this at all.
- The **shortcuts overlay** itself — but properly cropped (see below).

**Does not warrant an image (and should not get one):**
- **Screen-reader behaviour.** `aria-controls`, `aria-haspopup="dialog"`,
  `aria-expanded`, `role="dialog"`, accessible names — none of this is visible in a
  screenshot. A VoiceOver/NVDA rotor capture would show the assistive tool's UI, not
  Cellucid's, and would go stale with the AT's own releases. Keep this as prose, or as
  a transcript of announced strings if evidence is wanted.
- **"What keyboard cannot do"** (lasso, fine orbit, fast slider drags). An absence has
  no picture.
- **Motion sensitivity.** A still frame cannot convey motion; if evidence is wanted it
  is a short muted video or an animated comparison, not a PNG.
- **Colourmap recommendations** (Viridis/Cividis vs rainbow). Better served by the
  colormap documentation in `d_fields_coloring_legends` than by duplicating swatches
  here.

**needs** —
1. ~~**REPLACE the shortcuts figure.**~~ **Done.** The old 1440×900 frame gave roughly
   **19%** of its area to the Keyboard Shortcuts panel it was captioned as showing; the
   rest was the expanded Performance Benchmark panel, a canvas rendering **synthetic
   benchmark data** (`POINTS 4K · FPS 120`, `Data Pattern: Model Surface`) and a stray
   `×` close button. It is now the `accessibility-panel-keyboard-shortcuts` scenario in
   `docs/_tooling/screenshots/scenarios.mjs`: **pancreas**, the six sections that open on
   first load collapsed, `#shortcuts-section` opened and cropped whole so its summary
   stays in frame, DPR 2. The emitted width is the size rule's, not a chosen number —
   258 CSS px × 2 = `:width: 516px`.
2. **ADD focus ring** — `Tab` to the `#clear-all-highlights` button (a high-contrast
   target on a light panel) and capture with the ring rendered. Crop to
   `.highlighted-groups-header` so the ring is large in frame. **No** composited
   cursor — a mouse pointer beside a focus ring teaches the wrong thing; the whole
   point is that no pointer is involved. Add a `⇥ Tab` glyph overlay instead.
   `accessibility/focus-ring.png`, `:width: 246px`.
3. **ADD the `i` dialog** — an `.info-btn` focused and its `.info-tooltip`
   (`role="dialog"`) open beneath it, with the adjacent control label still visible,
   which is the page's stated design intent ("without hiding the nearby action label").
   Use `#cinematic-navigation-info-btn`. Crop to the enclosing `.control-block`.
   `accessibility/info-dialog-open.png`, `:width: 300px`.
4. **ADD the colourblind A/B** — the payload figure of this page. **suo**, colour by
   `organ` (9 categories, a genuinely hard case), Figure Export → `Show preview` ticked
   → the `Colorblind preview` select (`.figure-export-colorblind-select`, real options
   `Normal / Deuteranopia / Protanopia / Tritanopia`). Shot A `Normal`, shot B
   `Deuteranopia`, both cropped to the preview canvas at identical bounds so they can
   be flipped between. Cursor on the select in B.
   `accessibility/colorblind-normal.png` / `accessibility/colorblind-deuteranopia.png`,
   `:width: 360px` each.
   Note the existing `figure_export/preview.png` shows the right control but is another
   full-window 1440×900 shot of a **120-point fixture** — unusable here.

Create `_static/screenshots/accessibility/` for 2–4; there is no such topic folder yet.

**notes** — The old image baked in an app string that contradicted the selection docs:
its Highlighting block read `Shift + click — Add to selection`, while the renderer
requires `altKey` on every selection path (`highlight-renderer.js` lines ~3139, ~3200,
~3299, ~3352) and the real binding is **`Shift+Alt`+click**. That was fixed in
`index.html` — and pinned by `tests/keyboard-shortcut-panel-contract.test.mjs`, which
now ties every row of the panel to the handler that serves it — before the recapture, so
the current image ships the corrected panel and its five groups.

---

## {doc}`/user_guide/web_app/o_accessibility_privacy_security/02_privacy_model`

**explains** — What leaves the machine per loading workflow, the three storage
lifetimes (memory / HTTP cache / `localStorage` + `sessionStorage` + IndexedDB), how to
clear traces, what session bundles and figure exports embed, and the DNS-rebinding
`Host` validation with its `421 Misdirected Request` behaviour.

**has** — `none`.

**verdict** — **ADD** (2 figures — both DevTools/terminal, both load-bearing).

**needs** — Most of this page is correctly *not* screenshot-shaped: a workflow matrix,
a threat model and an IRB checklist are tables, and "no network traffic" has no
picture. Two exceptions where a capture converts an assertion into evidence:
1. **DevTools → Application → Storage**, which the page walks through step by step
   (lines 149–160) and never shows. Capture with Cellucid loaded and Community
   Annotation used once, so `localStorage` shows `cellucid`-prefixed keys,
   `sessionStorage` shows the token slot, and IndexedDB shows the annotation cache —
   i.e. exactly the three lifetimes described. **Redact the token value**, and use a
   throwaway profile so no real repo name or username is in frame. Crop to the DevTools
   Application pane. Cursor on the `Clear site data` button.
   `privacy/devtools-storage.png`, `:width: 860px`.
2. **A terminal proving the `Host` check**, since the page makes a strong security
   claim (`421` for any unexpected `Host`) and the reader has no way to verify it.
   `curl -sS -H 'Host: evil.example' http://127.0.0.1:PORT/ -i | head -5` beside a
   passing `curl -sS -H 'Host: localhost:PORT' ...`, showing `421 Misdirected Request`
   and `200` respectively. Light theme, no username/hostname in the prompt, port
   visible. `privacy/host-check-421.png`, `:width: 720px`.

**Explicitly NOT worth an image on this page:** the workflow matrix, the
`localStorage` key list (better as a code block), the "no network for local folder"
claim (a DevTools Network panel of a *quiet* tab is an unconvincing negative — the
existing prose about app assets still loading covers it better), and anything showing
a real session bundle or figure metadata, since the page's own advice is to redact
before sharing screenshots.

**notes** — Lines 41–43 already warn that screenshots leak identifiers. Every capture
specified above must be taken in a clean browser profile with a synthetic repo name.
The `--allowed-host` / `allowed_hosts=[...]` guidance and the
`jupyter-server-proxy` consequence match the Python-side documentation reference.

---

# LIFECYCLE A: Selection — from nothing selected to a named page used elsewhere

**Anchor pages:** {doc}`/user_guide/web_app/f_highlighting_selection/02_selection_tools_document_each_tool`
(steps 2–13) and `07_screenshots.md` (the sequence as a whole).

**Real tools** (`STEP_CONTROL_TOOLS`, and toolbelt order in `index.html:1546-1549`):
`annotation` (Annotation based) · `knn` (KNN drag) · `proximity` (Proximity drag) ·
`lasso` (Lasso). No others exist.

**Real states** (`lasso-selection.js:344-390`, `mode-ui.js`, `step-controls.js`):

| State | What is on screen |
| --- | --- |
| **S0 idle** | `#highlight-mode-description` = `HIGHLIGHT_MODE_COPY[mode]`; no `#<tool>-step-controls` element exists |
| **S1 in gesture** | preview highlight repainted every pointer-move; lasso polygon / proximity radius circle / KNN growth / annotation range label overlaid on canvas |
| **S2 step recorded** | `Step N: X cells` + a `.lasso-mode-tag` of `intersected` \| `+added` \| `−removed` \| `current selection`; step controls present; Undo enabled from step 1 onward; Confirm enabled iff X > 0 |
| **S2′ empty** | `No selection`; step controls present; **Confirm disabled** |
| **S3 abandoned** | one `<small>` notice appended under the readout; step counter unchanged, no undo spent |
| **S4 confirmed** | step controls removed, focus returns to the pressed `[data-mode]` button; new `.highlight-item` in `#highlighted-groups-list`; `#highlight-count` updated |
| **S5 page level** | tab strip rename input / colour input / `∩`/`∪` combine menu / Create Categorical |

**Standing capture conditions** for every step: viewport 1440×1000, DPR 2, light theme,
sidebar at its default 246 CSS px. Dataset **pancreas** (`?dataset=pancreas`, 3,696
cells, 8 `cell_type` categories, ships `connectivity/`) unless a step names another —
small enough that individual points are hittable without LOD confusion, and its cell
types are recognisable to a wet-lab reader. Steps needing LOD or scale use **suo**;
the "no neighbour graph" step needs the local fixture.

| # | Before | Action | After | Crop target | Cursor on | Filename (`_static/screenshots/highlighting_selection/`) |
| --: | --- | --- | --- | --- | --- | --- |
| 1 | Dataset loaded, nothing selected | Open the **Highlighting** accordion; click `[data-mode="lasso"]` | S0 — lasso help text, no step controls, empty group list | `#highlighted-cells-section .accordion-content` | `[data-mode="lasso"]` | `life-01-idle-lasso-selected.png` |
| 2 | S0 | Hold `Alt`, drag a closed loop around the Beta cluster — **capture mid-drag** | S1 — dashed polygon on canvas, enclosed points preview-haloed | canvas region only, 720×540 | polygon head (the moving vertex) | `life-02-lasso-mid-drag.png` |
| 3 | S1 | Release over the plot | S2 — `Step 1: 1,062 cells`, step controls appear, Undo enabled, Confirm enabled | `.highlight-mode-box` | `#lasso-confirm-btn` | `life-03-step-1-controls.png` |
| 4 | S2 (step 1) | Switch to 3D, rotate ~90°, `Alt`+lasso the same region | S2 — `Step 2: 640 cells` + `intersected` tag | `.highlight-mode-box` | none | `life-04-step-2-intersect.png` |
| 5 | S2 (step 2) | `Shift+Alt`+lasso a neighbouring lobe | S2 — `Step 3: 812 cells` + `+added` tag | `.highlight-mode-box` | none | `life-05-step-3-union.png` |
| 6 | S2 (step 3) | `Ctrl/Cmd+Alt`+lasso a corner | S2 — `Step 4: 705 cells` + `−removed` tag | `.highlight-mode-box` | none | `life-06-step-4-subtract.png` |
| 7 | S2 (step 4) | Click **↩ Undo** | S2 back at step 3; **↪ Redo** now enabled | `#lasso-step-controls` (tight) | `#lasso-undo-btn` | `life-07-undo-redo.png` |
| 8 | S2 | Click **Confirm** | S4 — group `Lasso (3 views) (812 cells) (812)` in the list, `812 cells highlighted`, step controls gone | `#highlighted-groups` | `#lasso-confirm-btn` (pre-click frame) | `life-08-confirmed-group.png` |
| 9 | S4, no prior KNN use | Click `[data-mode="knn"]`, `Alt+click` a cell — **the first click loads edges and selects nothing** | toast `Neighbor graph ready (N edges)` | sidebar column + toast | the clicked cell | `life-09-knn-graph-loads.png` |
| 10 | edges loaded | `Alt+click` a seed and drag upward | S1→S2 — degree grows in discrete steps, `Step 1: N cells` | canvas 720×540 + inset `.highlight-mode-box` | seed cell | `life-10-knn-degree-drag.png` |
| 11 | S0 | Click `[data-mode="proximity"]`, `Alt+click` a cell, drag outward | S1 — radius circle overlay, points inside preview-haloed | canvas 720×540 | circle edge (the drag point) | `life-11-proximity-radius.png` |
| 12 | S0, colour-by `cell_type` | Click `[data-mode="annotation"]`, `Alt+click` one Beta cell | S2 — the whole `Beta` category selected (`Alt` = replace on categorical) | canvas 720×540 + inset readout | the clicked cell | `life-12-annotation-categorical.png` |
| 13 | S0, colour-by `S_score` (continuous) | `Alt+drag` vertically from a cell | S1 — **range preview label** near the pointer, band preview-haloed | canvas 720×540 | the range label | `life-13-annotation-range.png` |
| 14 | various | Reproduce the five refusal messages (see the `06` record for triggers) | S3 ×5 | `#highlight-mode-description` ×5, stacked | none | `selection-notices.png` |
| 15 | S4, one page | Double-click the tab name, type `Beta cells`, `Enter`; click the colour swatch, pick a colour | page renamed and recoloured | `#highlight-pages-container` | tab name (rename frame) | `life-15-page-rename-recolor.png` |
| 16 | Two pages `Beta cells`, `Alpha cells` | Click `.highlight-page-tab-combine` on `Beta cells` | combine menu: `∩ Intersection with Alpha cells` / `∪ Union with Alpha cells` | tab strip + floating menu | the `∩` option | `life-16-page-combine.png` |
| 17 | Named pages exist | Open **Analysis → Quick** | header reads `Beta cells: 1,062 cells` — the named page is now the analysis unit | `#analysis-panel-simple` (top third) | none | `life-17-page-drives-analysis.png` |

**Handover:** step 17 is deliberately the same frame that opens Lifecycle B. The two
sequences should be captured in one session so the page names, colours and cell counts
agree across sections F and H.

**Cursor policy for this sequence.** Composite a cursor only where the pointer *is* the
information: the moving vertex of a lasso (2), the drag point that sets a radius (11),
the seed cell (10, 12), the range label anchor (13), and the button about to be clicked
(1, 3, 7, 8, 15, 16). Never on a pure state readout (4, 5, 6, 14, 17) — a cursor there
implies an action that is not happening.

---

# LIFECYCLE B: Analysis — from a selection to a configured analysis to a result to an export

**Anchor page:** {doc}`/user_guide/web_app/h_analysis/11_screenshots` (rebuilt), with individual steps
cross-referenced from `03`–`09`.

**Real modes** (`comparison-module.js:200-253` — these are the DOM `data-mode` values,
which is what a capture script must target):

| Label in UI | `data-mode` | Panel id | Runs |
| --- | --- | --- | --- |
| Quick | `simple` | `#analysis-panel-simple` | auto |
| Detailed | `detailed` | `#analysis-panel-detailed` | auto |
| Correlation | `correlation` | `#analysis-panel-correlation` | auto |
| Differential Expression | `differential` | `#analysis-panel-differential` | **`Run Differential Expression`** |
| Gene Signature | `signature` | `#analysis-panel-signature` | auto |
| Marker Genes | `genesPanel` | `#analysis-panel-genesPanel` | **`Discover Markers`** |

**Real result surfaces:**
- sidebar preview inside the mode's `.analysis-accordion-content`
- `⤢ Expand` → `.analysis-modal` › `.analysis-modal-content` with four regions:
  `.analysis-modal-plot`, `.analysis-modal-options` (holding the `Export:` row —
  **PNG · SVG · CSV**, no PDF — and `PLOT OPTIONS`),
  `.analysis-modal-stats-panel` (`SUMMARY STATISTICS`),
  `.analysis-modal-annotations-panel` (`STATISTICAL ANALYSIS`)
- `.analysis-accordion-copy-btn` → a floating `.accordion-section.analysis-window-panel`
- DE progress phases, literally: `Loading & Computing` → `Multiple Testing Correction`

**Standing capture conditions:** viewport 1440×1000, DPR 2, light theme. Dataset
**pancreas** for every step (3,696 cells / 3,753 genes → DE and marker discovery finish
in seconds; `Ins1`, `Ins2`, `Iapp`, `Gcg` are present and recognisable). Two pages
prepared by Lifecycle A: **`Beta cells`** and **`Alpha cells`**.

| # | Before | Action | After | Crop target | Cursor on | Filename (`_static/screenshots/analysis/`) |
| --: | --- | --- | --- | --- | --- | --- |
| 1 | Lifecycle A finished | — | Two named, coloured pages with non-zero counts | `#highlighted-cells-section .accordion-content` | none | `life-01-two-pages-ready.png` |
| 2 | Analysis section closed | Open **Analysis** | All six modes listed, all collapsed, each with a copy button and chevron | `.analysis-accordion` (widened so descriptions do not truncate) | none | `analysis-panel-modes.png` |
| 3 | — | Click **Quick** | Auto-run: `Beta cells: 1,062 cells`, Composition + Statistics, `▶ Page Selection (Dynamic)` | `#analysis-panel-simple` | `#analysis-header-simple` | `life-03-quick-autorun.png` |
| 4 | — | Open **Detailed**; set `Categorical obs → cell_type`; select both pages; `Bar Plot` + Percentages | Grouped bars, both pages, no modebar tooltip in frame | `#analysis-panel-detailed` | the `Plot Type:` select | `detailed-categorical.png` |
| 5 | step 4 | Click `⤢ Expand` | Modal: plot + `PLOT OPTIONS` + `SUMMARY STATISTICS` + `STATISTICAL ANALYSIS` (chi-squared / Cramér's V card) + `Export:` row | `.analysis-modal-content` | `⤢ Expand` (pre-click frame) | `detailed-expanded.png` |
| 6 | — | Open **Correlation**; `X = S_score`, `Y = G2M_score`, `Color by = cell_type`, Pearson | Structured scatter with a real `r`, `n`, trend line — not a null result | `#analysis-panel-correlation` | `Correlation method:` select | `correlation-results.png` |
| 7 | — | Open **Differential Expression**; open the `Page B:` select | Two `<optgroup>`s visible: **Pages** and **Wildcards** (`Rest of Beta cells`) | `#analysis-panel-differential` (top) | the open `Page B:` select | `de-page-selection.png` |
| 8 | step 7 | Set `Statistical method: Wilcoxon`; expand `Performance Settings` | Batch size / memory budget / network parallelism / compute parallelism all visible | `#analysis-panel-differential` (mid) | the Performance Settings disclosure | `de-performance-settings.png` |
| 9 | configured | Click `Run Differential Expression`; capture mid-run | Button reads `Running...`; progress on phase 1, `Loading & Computing` | `#analysis-panel-differential` | the disabled run button | `de-running.png` |
| 10 | phase 1 | Wait for phase advance | Progress reads `Multiple Testing Correction` (use **suo** + `Rest of` if pancreas is too fast to catch) | progress tracker (tight) | none | `de-running-phase2.png` |
| 11 | run complete | — | Sidebar result: the volcano preview and the `⤢ Expand` button, and nothing else — `Top Differentially Expressed Genes` is drawn only in the modal (`de-analysis-ui.js:28`) | `.analysis-preview-container` ∪ `.analysis-expand-btn` | the `⤢ Expand` button | `de-sidebar-result.png` |
| 12 | step 11 | Click `⤢ Expand` | Full modal: volcano with labelled top-N genes, threshold controls, FDR checkbox, `Export: PNG SVG CSV`, summary table, ranked table with Top 5/10/20/100/All | `.analysis-modal-content` | the `\|log2FC\|` threshold slider | `de-volcano-expanded.png` |
| 13 | step 12 | Drag `\|log2FC\|` from 0.5 to 3 | A/B pair: populated volcano → nearly empty volcano | `.analysis-modal-plot` ×2, identical bounds | the slider (B frame) | `volcano-threshold-a.png` / `-b.png` |
| 14 | step 12 | — | Tight crop of `.de-modal-stats`: `Genes tested (FDR denominator)` / `Not tested (< 10 cells with a value)` / `Significant (…)` / `Upregulated` / `Downregulated`, with *Not tested* non-zero | `.de-modal-stats` | none | `de-fdr-denominator.png` |
| 15 | — | Open **Gene Signature**; paste `Ins1, Ins2, Iapp`; both pages; Mean expression; None; Violin | Two clearly separated violins + gene chips confirming all three resolved | `#analysis-panel-signature` | `Normalization:` select | `gene-signature.png` |
| 16 | — | Open **Marker Genes**; `Group By: cell_type`; `Mode: Ranked Genes`; click `Discover Markers` | The gene heatmap every mode draws in the sidebar, with the `Discover Markers` form above it — **not** a ranked list; that table exists only in the expanded view | `#analysis-panel-genesPanel` | none | `marker-genes-ranked.png` |
| 17 | step 16 | Switch to `Mode: Clustered`, `⤢ Expand` | Heatmap with gene and group dendrograms, its gene axis printing many more names than the sidebar box has room for | `.analysis-modal-plot` | none | `marker-genes-clustered.png` |
| ~~18~~ | ~~any expanded result~~ | ~~Click `CSV`~~ | **WITHDRAWN** — the browser download shelf renders a real local path and cannot be sanitised | — | — | ~~`export-download.png`~~ |
| ~~19~~ | ~~file downloaded~~ | ~~`head -3` in a terminal~~ | **WITHDRAWN** — terminal output is transcribed as a verbatim `text` block, never photographed | — | — | ~~`export-csv-terminal.png`~~ |

**Cursor policy for this sequence.** Composite a cursor on the control that produces
the next state (2→3, 4, 6, 7, 8, 9, 12, 13B, 15) and nowhere else. Steps 10, 11, 14,
16 and 17 are read-only states; a cursor in them is noise. Step 12 is the exception the
capture forced: the `|log2FC|` slider is in the options column, outside the
`.analysis-modal-plot` crop of step 13, so a cursor drawn on it would not be in either
frame — both threshold figures ship without one, and their alt text says so.

**Why pancreas and not suo for this sequence.** The old marker-genes shot recorded a
17.8 s run over 561,947 cells and 2,185 genes grouped by `organ` into 9 groups, and its
heatmap was unreadable at sidebar width. `pancreas` has 3,696 cells, 3,753 genes and
8 `cell_type` groups; every mode completes fast enough to capture reliably and
`Ins1`/`Ins2`/`Gcg` make the biology self-evident to the wet-lab reader the objective
names. Keep **suo** for the two steps that need scale (10, and any LOD demonstration in
Lifecycle A). Note that "every axis label fits" was never true of a 224 px box and is
not what fixed the smear: the axis is thinned to the labels the box can show, so
**pancreas** in the sidebar prints one cell-type name out of eight, and five out of
eight in the expanded view.

---

# Cross-cutting recommendations

1. **Fix the capture pipeline before capturing anything.** Every image in these three
   sections is DPR 1 rendered at 1:1. Standardise on DPR 2 with `:width:` set to the
   CSS width, and record the viewport, DPR, theme, dataset and app build in the two
   `*_screenshots.md` pages. That single change improves all 17 existing figure
   references.
2. **Fix four prose defects before they are photographed**, because a new screenshot
   will make each one visibly wrong:
   - the non-existent DE symptom string `Need at least 1 page for comparison`
     (`h/06:237`),
   - the non-existent `neighbor graph not available` string (`f/02:256`, `f/05:172`,
     `f/06:291`),
   - `PDF, PNG, and CSV export actions` — there is no PDF export (`h/09:239`,
     `h/11:66`),
   - `Some builds include a "Category builder" UI` — it is unconditional and named
     **Create Categorical** (`f/03:197-204`).
3. **Retire {doc}`/user_guide/web_app/h_analysis/11_screenshots` as a duplicate gallery.** All seven of its
   figures are copies of images already on the mode pages, so every caption defect is
   duplicated there. Rebuild it as the Lifecycle B sequence.
4. **Pick pancreas as the standard docs dataset for procedural captures.** It is the
   only sample with genuinely continuous obs (`S_score`, `G2M_score`), it is small
   enough that every mode completes in seconds, its `cell_type` cardinality (8) fits
   every legend and axis, and it ships both `connectivity/` and `vectors/`.
5. **Create `_static/screenshots/accessibility/` and `_static/screenshots/privacy/`.**
   Section O currently borrows one image from the `web_app` topic and has no home of
   its own.
