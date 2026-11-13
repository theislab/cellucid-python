# Decision record: analysis

This page answers: **"Why are the analyses built the way they are: one
function on a worker or on the page, genes read as whole columns outside the
field state, results kept in the history as they were exported, and every
change of a panel one action?"** It records the decisions with their dates,
the measurements that justify them, what was tried and taken out, and what
is still open. How the subsystem works is {doc}`11_analysis_architecture`,
and this page does not repeat it. Update this page in the same change as the
code it describes.

## At a glance

**Audience**
- Developers changing `app/analysis/`: a mode, a plot, a compute operation,
  the gene reads behind an analysis, or what a result holds.

**Time**
- 5 minutes today; the tables grow as history moves here from the code.

**Prerequisites**
- {doc}`11_analysis_architecture`.

:::{note}
This record was started on 2026-11-12 from the reference page and the code,
and the same day the history that stood in the comments of `app/analysis/`,
its tests and its browser specs moved here. A date written "by 2026-11-12" is
the day the entry was first written into these docs, where the day of the
decision itself was not recorded.
:::

---

## What it does now

Paths are under `cellucid/assets/js/app/analysis/` unless they name another
place.

1. **Analysis is created by the application, not by the interface's
   coordinator** (`comparison-module.js` `createComparisonModule`, from
   `app/main.js`), and has its own interface under `ui/`.
2. **An analysis runs on highlight pages, the visible cells and variables**
   (a categorical or continuous `obs` field, a gene). Pages are state, so
   they are an analysis's inputs.
3. **One data layer reads for every analysis** (`data/data-layer.js`): page
   data cached by page version, and whole gene columns read from the source
   for a computation, never installed in the field state.
4. **A compute operation is one function that runs in a worker or on the
   page** (`compute/`), so where it runs never changes a number.
5. **Every analysis change is an action** (`app/actions/definitions/analysis.js`,
   provided by `actions.js`): open, close, pages, run, one `set` a mode made
   from the mode's options schema, and the window actions. A person's change
   in a panel or a window is recorded by observation, one step each.
6. **What the analyses hold is the `analysis` state domain**: each window's
   and the docked panel's mode, settings and result. A result is kept as the
   analysis exported it and put back as it was, never computed again by a
   restore; where a window stands is the `layout` domain's.
7. **A share that leaves results out** writes a record of each, and opening
   that file runs each analysis again and says where the new result was
   computed on something else (`app/session/recompute-results.js`).
8. **The options of a mode are one schema** (`options/`): the mode's `set`
   action takes its parameters from it, and the defaults tables of the user
   guide are held to it by a test.

## Decisions and why

Paths are under `cellucid/assets/js/app/analysis/`.

### The panel, runs and actions

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | The panel shows one analysis at a time. An action that draws on an analysis of the panel shows it when the panel shows another, as a person opens an analysis to use it, and is refused only when it still cannot draw (the section or the sidebar closed); the shown one is put back (`actions.js`) | An agent that set up two analyses was refused "not shown" when it ran, set the options and chose the pages of the first after the second was shown; one that had run differential expression and then markers was refused "not shown" twice for the first one's options, whether or not it passed `show: true` | Refusing an action for the panel's showing another analysis |
| by 2026-11-12 | An action's Stop reaches a run before its job exists: stopped while the run waits to start, it throws why and nothing runs; stopped once it runs, the job is cancelled as Cancel does. A run that ended without a result is reported as stopped (by an agent's Stop, Cancel or another run), never as a run that changed nothing | A run stopped before its job existed used to run to its end; an agent was told "it was so already" | A Stop that only reaches a started job |
| by 2026-11-12 | `runNow` waits for the groups' read and refuses only with a reason a read cannot change; a failed gene read is reported with its own cause; reading an analysis (`analysis.<mode>.get`) makes no mode's analysis and so changes the `analysis` domain with no step | An agent or a notebook running right after setting the options was refused with "Reading the groups..." when the read had not ended (WebKit, a CI run, `analysis.run` of the markers); a failed read was reported as "invalidated by an analysis request ownership change" whatever the cause; undo of a turn left the domain changed after a get | Refusing a run for a read in flight |
| by 2026-11-12 | A change a person makes in a panel is a step of its own even while an action runs (the dispatcher's `recordObserved`) | A person's change is a step whatever else is under way | Folding it into the action |
| by 2026-11-12 | An analysis option whose choices are highlight pages (Differential expression's Group A and Group B) takes the page's name in the actions: the set, the get and a person's recorded change speak names | An agent read `pageA: "restof__page_1"`, was refused "Group A has no choice" for the page's name, then guessed "page_2" | Page ids in the actions |
| by 2026-11-12 | A set of groups has one form: the schema keeps a caller's order and puts ids the data does not list after them, so the picker reading it back changes nothing (`options/option-schema.js`). Gene names do not go in the "Exclude genes" families list (`means` points to the pattern field) | An agent's list of groups in its own order was rewritten by the next form refresh into the field's order and the history recorded "Changed outside the history"; an agent sent gene names to the families list ("Exclude genes has no choice "GSC"") | One set in two orders |
| by 2026-11-12 | Changing the fields of the dataset (created, renamed, deleted, categories edited) makes every analysis read their categories again, in the panel and in every window. Highlight and filter changes reach every built analysis: a result kept in another mode is marked stale, and a filter change marks stale an analysis scoped to "Visible cells" | An agent's new field was "no choice" of Differential expression until another change; only the active mode heard highlight changes, so a result kept in another mode was never marked stale and its groups were cached against old highlights, and nothing heard filter changes | Only the active mode hearing a change |
| by 2026-11-12 | One registry of analysis modes (`AnalysisUIManager`, filled by `ComparisonModule._registerAnalysisTypes()`) | A second hand-written `ANALYSIS_TYPES` table in the `analysis-types` barrel was never read, so it drifted unnoticed (`Quick Insights` and `Correlation Explorer` where the app shows `Quick` and `Correlation`; Differential Expression as `minPages: 2, maxPages: 2` where the registration accepts one page) | A second table of modes |
| by 2026-11-12 | The host markup owns the container id of an analysis (the target of the accordion header's `aria-controls`); the manager never overwrites it. With no dataset open an analysis says so, never "Requires gene expression data." | Overwriting it left every header pointing at an element that does not exist; a missing dataset read as a missing kind of data | Overwriting the host's id |
| by 2026-11-12 | The page selection shared by all modes is validated against the selectable domain with the derived "Rest of X" pages (`selectablePageIdSet`), in the manager, in `BaseAnalysisUI` and in `page-derivation-utils.js`: one domain, one function | Validating against the base pages alone made choosing "Rest of Page 1" in Detailed analysis throw `page "restof__page_1" was not found`, and Differential Expression and Marker Genes (which build no page selector) threw the same when a rest-of page selected in another mode was handed to them on a mode switch; `_requireCustomPageColors` read the widget's current view and refused a color for a rest-of page out of one saved session | Validating against the base pages |
| by 2026-11-12 | Detailed analysis builds its variable selector once and re-attaches it across rebuilds, so the control that raises a change is never destroyed by it | The selector reports `('', '')` when a type is chosen while its variable is still unchosen; the panel read that as "selection cleared", wiped `dataSource.type` and rebuilt the selector from the empty config, so the type just chosen came back as "None" | Rebuilding the control that raised the change |

### Results, tests and what a result says

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | Composition tests: a type is tested over the pairs holding any of its cells; the t-tests drop a pair with none on either side, as the signed-rank test does (a zero difference); below `MIN_PAIRS_WITH_CELLS` (3) such pairs a type is not tested | In a pair with none its two transformed proportions are the units' sizes alone, so the difference has the same sign for every donor whose one organ gave more cells than the other, which a paired t reads as a change of a type that is not there: on Suo, liver against spleen (7 donors with both), myelocytes (cells in one donor's spleen only) came out adjusted p 0.03 from seven pairs. With two pairs a paired t has one degree of freedom and the signed-rank test cannot reach p < 0.5 | Testing a type over pairs with none of its cells |
| by 2026-11-12 | Composition tests, paired: the mean proportions and their ratio are over the paired samples only, and a type with cells only in units left out of the pairs is not tested | With every sample, a donor with 2 liver cells and none in the marrow turned the log2 ratio of LYVE1+ macrophages toward the liver while the paired test (6 of 6 donors) found more in the marrow; one Langerhans cell in a spleen-only donor came out "significant", adjusted p 0.04, from seven pairs without it | A mean over units the test never used |
| by 2026-11-12 | The composition table file carries the statistic, the samples each condition's mean is over and, paired, the pairs the test took, beside the rows of the agent's table; without a condition the agent's rows carry numbers (cells, pooled proportion, per-sample shares) | The file once gave every sample as the test's (`n_samples_Liver 14` beside a paired t over 7 donors) and left out the per-sample proportions; without a condition the rows held only the names, and the per-sample shares a donor-level argument needs were never given | A file that says more samples than the test used |
| by 2026-11-12 | A 2 x 2 table takes Yates's correction, as SciPy `chi2_contingency` and R `chisq.test` do and Composition does (8.1, p 0.004426525857919832 for the 15/5 against 5/15 case; Cramer's V stays on the uncorrected statistic). A Mann-Whitney U above 2^24 is never held in a Float32 | An odd integer above 2^24 (16,777,559 for a fixture) cannot be held by Float32: stored in one and shipped to the engine, it was thrown away | An uncorrected 2 x 2 test; Float32 for a rank statistic |
| by 2026-11-12 | A result's "about" block (differential expression, markers, gene signature, dot plot) states what the rows alone cannot: the two groups and their cells, the test, the correction, the matrix and its values, the cutoffs passed; each tested marker group's cells; the tests between pages for a signature; what a mean is of. In the differential expression rows `top` counts each direction, and the wording says so. The result of "As stored" differential expression says that a ratio of means of log-scale values is not an "X-fold" change of expression | An agent read differential expression rows alone and had to guess the group sizes and that the values were raw counts; it weighed "a marker of 23 cells" without knowing it was 23; it had signature means alone and could not say whether a difference between pages held; it read raw-count means as expression levels; it asked for 40 each way and read top 3,753 first; it wrote "the means differ about 16-fold on this log scale" | Rows with no statement of what they are of |
| by 2026-11-12 | A differential expression read on a gene sample says each group's cells among the sample and that a group is skipped for too few sampled cells | A walk on a CELLxGENE atlas (2026-10-31; genes read for 6,320 of 584,944 cells, comparing 108 against 94 cells) said "Every cell: computed on all cells of its groups" above "Genes read for 6,320 of 584,944 cells", and "15 skipped" of 50 cell types with no reason given | A result that says all cells for a sample |
| 2026-11-10 | Where p-values tie (many underflow to 0 on hundreds of thousands of cells) the stronger test statistic comes first, before the larger fold change, in the volcano's labels, its CSV and the result's own table; the score reader accepts an infinite statistic (a pseudobulk t with no spread) | On Suo (NK, 30,283 cells, against the other 531,664) 486 genes had p = 0; a volcano labelled HBE1, KLF1, RHAG (strongly down, log2FC about -10) while the result table's top genes were KLRD1, PRF1, KLRF1, IL2RB, GZMH (score 544 to 450, log2FC about 5); a pseudobulk DE had failed with "score for CD3E must be a finite number or unavailable" | Ranking ties by fold change |
| 2026-11-10 | Quick insights accepts a continuous field as Float32 or Float64 (the reader gives either) | Suo's `n_genes`, a float64 integer count, failed with "requires a Float32Array values payload" | Refusing a Float64 field |
| by 2026-11-12 | Analysis reads categorical codes at the width the field stores them | Fields were copied into a Uint16Array, which turned a uint32 code into code mod 65,536 (a cell of category 65,536 counted as category 0, no error) | A narrower copy of a code |
| by 2026-11-12 | An analysis's page cache identity is the set of the page's cell indices, and the data layer keeps the page's index set beside its digest from one traversal (`_pageIndexSnapshots`) | Versions were built from the first five, last five and (above 100 cells) a tenth-step middle sample, and the cache key reduced to `pageId@count`, so replacing almost every cell of a page while keeping its cardinality gave a byte-identical key and Correlation, Detailed, Quick Insights, Gene Signature and the comparison module served results for the previous cell set; the digest and the values were built by two traversals with an `await` between them, so a membership change in that window was cached under a version it no longer matched | A key built from a sample of the cells |
| by 2026-11-12 | A categorical analysis's plot data that is degenerate gives an honest empty state, and a marker result with no gene passing the filters shows one empty state in the panel and in the expanded view alike, with how many genes were tested in how many groups and Relax filters | When no trace remains the render pipeline handed Plotly the normal axis layout, and a single-cell page drew bare 0 to 6 axes around nothing; the expanded view said only "No genes pass the filters." with no way on | A plot of nothing; two empty states |

### The expanded view and the controls

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | The expand affordance is drawn once by the base form for every mode, as a real `.analysis-expand-btn`; it is disabled ("Opens once the result is drawn") until the result is published | DE, Marker Genes and Gene Signature each appended an identical `.analysis-actions` row, while Correlation asked only for the preview container's `clickable` shortcut, a bare click listener on a `<div>` with `tabindex="-1"` and no role, so its expanded view was unreachable by keyboard and unannounced; the preview is drawn first and the result published once its drawing has committed, and a click in between (DE of 562k cells on Suo, 2026-11-10) said "Expanded analysis view requires a completed analysis result" | A click listener on a div; an Expand that can be pressed too early |
| by 2026-11-12 | An open expanded view is drawn again when a new result is published, in every form-based mode, through one result owner that refreshes it | The overlay was bound to whichever result was current when it was opened (each form-based mode assigned `_lastResult` inline in five places), so changing the inputs left the two surfaces showing different results: on Gene Signature (Chromium) the sidebar preview moved to mean 1.448292 while the overlay still showed mean 3.966388 under a title naming the previous gene set; Composition drew its preview from a copy of the result, so its Expand never matched | A result assigned inline by each mode |
| by 2026-11-12 | The expanded view's card opens with a transform easing from scale(0.96) to scale(1), and nothing is measured off its client rectangle while it eases; a plot is measured from the host's layout size. An expanded view keeps the size it was given, with a `max-width` of the viewport less the margins, and a window resize moves a view that is out of the window back inside, once a frame | A client rectangle is the transformed box, so the old opening code that anchored the card from `getBoundingClientRect()` in a 200 ms timer and placed the panel-resizer corner in a 250 ms timer was off whenever the transition had not finished, by a different amount on every open (screenshots of one view differed run to run); a host measured by `getBoundingClientRect()` reports the visual box of the opening and a transform change never notifies a ResizeObserver, so the plot rendered 4% small in the overlay; a plot type pinning `layout.height` overflowed the 180 px sidebar preview and swallowed clicks on Expand; the cap was `none` (to lift the stylesheet's 860 px), so a view sized on a large window ran past a smaller window's edges | Measuring a transformed box; a cap of none |
| by 2026-11-12 | The expanded view is a `<dialog>` opened through `wireAppDialog`; the multi-select dropdown is a native popover and no modal | The view said `aria-modal` with a hand-made Tab trap and the page behind stayed reachable; the dropdown also stated `aria-modal` and held Tab in a hand-made trap, and its popup, appended to body, never received focus and Escape stranded it | A hand-made modal |
| by 2026-11-12 | The analysis controls are real controls: page chips and collapsible sections are buttons with a name and a state, selects have an accessible name, a copied window has no duplicate ids, lazily built analysis controls are in the session's UI-control inventory, and the form is not rebuilt on an accordion reopen | Found by driving the running app: `aria-controls` pointing at overwritten ids; chips as click-listener divs (no role, name, selected state); sections as divs (hidden from keyboard users); a copied window with duplicate ids (the copy's label toggled the sidebar's checkbox); a session unrestorable unless the same panel was open; the rebuilt form discarding choices | Div controls |
| by 2026-11-12 | The gene search list is a combobox: Enter, Escape and the arrow keys; a typed gene is active over None | The list was driven by `mousedown` alone (no keys, no combobox semantics), so picking a gene was impossible without a pointer; typing a gene left None active, so Enter cleared the choice instead of taking the gene typed | A list only a pointer can drive |
| 2026-11-10 | A page selector exposes the color every page it offers is drawn in (its swatch), so a plot colors pages as their swatches show, whatever order they were picked in | Detailed drew pages CD4 T (red swatch) and Page 1 (blue swatch), picked in that order, as blue and red: colored by pick order | Coloring pages by pick order |

### Plots

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | A box plot that draws points hands Plotly precomputed statistics and only the drawn points instead of every value (`plots/plotly-box-statistics.js`); box and violin points are drawn as one image per box by a worker (`plots/point-layer.js`), not as an SVG path per point; a histogram's values are grouped into the bins Plotly would assign (at most three of each bin's own values with weights, `plots/plotly-histogram-bins.js`) and a violin's density is computed here operation for operation as Plotly does (`plots/plotly-violin-density.js`), so the analysis workers sum it in parallel off the main thread | Plotly's own passes cost seconds on a large page (Measurements); the spec's reference is the raw-value trace the module built before statistics were precomputed | Plotly computing statistics of every value |
| by 2026-11-12 | The legend of a plot is left out of the automatic placement when the rows would squeeze the plot below `LEGEND_MIN_PLOT_HEIGHT` (100 px) | Entries that wrapped to more rows than the fixed top margin held were drawn over the plot's axis labels | Legend rows over the axis labels |
| by 2026-11-12 | Four label collisions were fixed, each reported against the committed documentation screenshots: the volcano's x-axis title `log2 Fold Change` against the `Down (n)` legend entry; categorical tick labels in the 224-pixel sidebar preview; the dot plot's size legend title against its color bar, its dots against each other in the expanded view and its gene names where the panel is narrow; the correlation plot's R2 boxes | See the rows below and Measurements; the correlation plot's R2 boxes were stacked 15% of the plot's height apart, less than one box is tall; the dot plot's legends were in fractions of a plot area whose height the preview does not fix | Offsets that hold at one size only |
| by 2026-11-12 | A volcano plot lays out the band below the x axis top-down in pixels (tick labels, the axis title at a pinned standoff, the legend below it), the margin being exactly the total, and `legend.y` is converted from pixels with the plot height the layout produces; a plot that already places its legend below the plot area keeps its own offset (`placesLegendBelowPlot`) | The legend was at `y: -0.14`, a fraction of the plot area's height, and the axis title wherever Plotly put it inside a fixed 90 px margin; two offsets measured from different origins into one band collide as soon as the box changes size, and `log2 Fold Change` was drawn across the `Down (n)` legend entry in the sidebar preview and in the expanded overlay alike | A constant that only holds at one size |
| 2026-11-10 | Volcano y ticks over the +infinity band are round values (1/2/5 x 10^n), none within 90% of the boundary | Quarters of the top read 80.2, 160, 241, 321 and 321 and overlapped +infinity | Quarter ticks |
| by 2026-11-12 | A categorical axis is anchored at both ends: a thinned axis still names its first and last category, so the capacity floor is two. The fit belongs in `PlotFactory.render`, where the box is known | Eight categories (`cell_type` and `clusters` of the published Pancreas dataset, one short of the rotation threshold) draw upright and `Pre-endocrine` needs a slot of about 89 px of the 149 px the sidebar leaves (a capacity of one), so eight bars stood under a single tick reading `Ductal`; a pair at the ends of an eight-band axis sits 7/8 of the box apart, 130 px, against the 41 px `Ductal` and `Epsilon` need; eighteen cell-type names overlapped in ten pairs at 224 x 180 and none at 720 x 380 (Measurements) | A fixed number of labels |
| 2026-11-10 | A dot plot makes the room of a name one whole pixel more than the box the engine draws its label in (`nameSlot`), and the box is the SVG label's as Plotly draws it, not the canvas's font box or `getBBox()` (`plots/text-measure.js`) | Each engine sizes the box its own way and Plotly puts each name on a pixel, so names 10.8 px apart are 10 px apart one row and 11 the next; a box of 10.8 overlapped its neighbor in Firefox (Measurements) | A canvas font box; `getBBox()` |
| by 2026-11-12 | A density curve is cached per page-values array and bandwidth; its support is a page's own Float32 array with unit weights past the distinct-value bound, not a Float64 copy with an array of ones | The plot is redrawn for every option change and once more for CSV export, and each redraw re-ran the kernel sum; the analysis workers each receive the support, so the copy was four times the bytes to send to each of them | Re-running the kernel on a redraw |

### Reading, sorting, caching and cost

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | The analysis cost model is fitted per set of genes and a cost is known only for the genes it was measured with; an unrecorded work gets no estimate | Scaled per gene, the time of a run of 6 genes (the live-agents fixture: 120 cells, 546 ms) told an agent's permission card that a 1-second marker run of 3,691 genes (Kanemaru) would take 88 hours (the live-agents run estimated about 102 h) | Extrapolating a cost per gene |
| by 2026-11-12 | The result cache is bounded by bytes only, with no entry count and no age | A count cap evicted small results while memory was free; an age let go of results nothing had made stale | An entry cap, a time to live |
| by 2026-11-12 | Order statistics, ranks and sorted copies come from a radix pass over the IEEE-754 key, not from the native numeric sort or a comparator index sort (below `NATIVE_SORT_MAX_LENGTH`, 131,072, the native sort stays); a numeric typed array is checked for non-numbers by constructor (the two BigInt views), not by scanning; one gene's record is found through the key index | The comparator form runs a callback per comparison; the scan ran twice per group on every gene of a differential-expression run; describing every gene of an 8,192-gene panel to read one cost a record per gene per read | Comparator index sorts for ranks |
| by 2026-11-12 | Quick insights' observation fields on first open and a gene signature's genes start every unit at once and the fetch scheduler paces the requests | Four fields (six requests with the categorical outlier masks) were read one at a time, and a 60-gene signature loaded in fixed batches of ten with a macrotask barrier between them (Measurements) | Serial rounds and fixed batches |
| by 2026-11-12 | A gene stream's read over some cells starts once half the gate is free (or every task left fits), not as each task leaves the gate; the genes of a task are sized from the bytes measured as read (a sparse gene's entries are tens of times fewer than a dense column), not held at the densest form (`shared/gene-stream-runner.js`, `shared/run-tuner.js`) | A read over some cells costs its cells whatever its genes: started as each task left the gate, a 10M-cell differential expression of two large groups made a read of every task (54 s against 4 s for the tree it replaced), and held at the densest form it went one gene a task, 30,000 tasks and page messages where the tree sent 128 (Measurements) | A read per task; one gene a task |
| by 2026-11-12 | A page of one group is that group's own cell list, and Differential expression's group resolution is kept | Copying a selection into an array cost 85 ms of a change at 5,000,000 cells, and Run's state is asked on every option, highlight and page change (Measurements) | Copying a page's membership per read |
| by 2026-11-12 | Chunked gene tables stop when their expanded view or UI goes away; `createSlicedYield` yields one macrotask | Each copied loop was guarded by a render revision nothing advanced on close and kept appending rows to a detached table; `setTimeout(0)` is clamped to 4 ms, which over a long chain of yields turns into seconds of idle time (a platform fact) | A loop that outlives its table |

## Measurements that justify it

Where a row names no machine, the comment it came from named none: the
numbers are orders of magnitude on one machine, not a benchmark.

| What | Where and when | Measured |
|---|---|---|
| Density plot, `PARALLEL_MIN_KERNEL_TERMS` = 3,000,000 terms (about 40 ms of work, within one main-thread task) | Machine not named | A 25,000-value page drew in 66 ms on the main thread and 105 ms through the workers; a 100,000-value continuous page held the main thread for 135 ms at the earlier threshold of 20 million terms. Summed on the main thread, a million-cell continuous field takes 1.4 s (a hundred evaluations per distinct value) |
| Box plot, `WORKER_MIN_VALUES` = 100,000 | Machine not named | The linear-time sort for the box statistics takes 30 ms of main-thread time for a million values and ten times that for ten million |
| Plotly's box `calc` on one million raw values (a point object per value, binned, sorted through a comparator) | Machine not named | About 0.7 s of main-thread work before anything is drawn; a zero-inflated gene whose quartiles are all zero makes every expressing cell an outlier: 170,000 SVG elements and about a second for a one-million-cell page |
| Plotly's histogram `calc` on raw values | Machine not named | About 0.2 s for a million cells, two seconds for ten million |
| Plotly's violin `calc` | Machine not named | About 400 million kernel terms for a million distinct values and a 2.4 s freeze |
| Box and violin points as SVG paths | A page of half a million cells with all points shown | Half a million DOM nodes and 5 s of main-thread work (one 2.5 s task), and each hover move 0.7 s |
| Observing the whole document subtree for Plotly's notes | One box plot of a million-cell page | A subtree query for every inserted element: 170,000 queries |
| Page-data CSV export | Half a million cells | A 15 MB text; joining one string per row on the main thread held it for about 70 ms and built the text twice |
| Counting category labels | Ten million labels | Walking them through a `Map` takes about 0.2 s of main-thread time each time |
| `DISPATCH_CLONE_BYTES_PER_TASK` = 96 MiB | The density plot of a ten-million-cell page handing each of eight workers its support and weights | `postMessage` structured clone runs at about 9 GB/s; dispatching every slice at once is about 1 GB in one task (114 ms); a task at the budget copies in about 10 ms |
| Gene stream runner: a read over some cells started as each task left the gate | A 10M-cell differential expression of two groups (1.9M cells) | A read of every task, 15,000 reads, 54 s against 4 s for the tree it replaced, which read a task a worker at a time |
| Run tuner: genes per task sized from the densest form | The same 10M-cell run | One gene a task, 30,000 tasks and page messages where the tree it replaced sent 128 |
| Re-deriving a page's membership (a `Set` of every index, a comparator sort, a digest) | A 250,000-cell page, and a fifty-gene signature | 26 ms for one page and 1.3 s of a 1.8 s signature |
| Copying a selection into an array for a count; `Array.from` of a page membership | A 5,000,000-cell selection; 10 million cells | 85 ms of the change that made it; 0.3 s and 40 to 80 MB per read. Differential expression's group resolution (a pass over every cell) took 8 ms a pass at 2M cells; the kept answer costs 0.01 ms |
| Ranks by a comparator index sort; a `Map` keyed by heap numbers | One million distinct values; a 500,000-cell gene | About a second for the sort; roughly 45 ns per lookup, so the lookup, not the statistic, was the cost of a rank test. The native sort beats the radix pass below about 200,000 values on a first call (`NATIVE_SORT_MAX_LENGTH` is 131,072) |
| Summary and tests of continuous pages (`WORKER_MIN_VALUES` = 200,000) | The interface's thread | Linear-time order statistics and tests still take tens of milliseconds per million values |
| A run of 6 genes extrapolated per gene | To 3,691 genes (a marker run) | Told a permission card that a 1-second run would take 88 hours |
| Quick insights' four observation fields on first open | Chromium at 150 ms emulated latency | Six strictly serial round trips spanning 2.07 s, never more than one request in flight; a 60-gene signature in fixed batches of ten was six sequential waves however fast the link |
| Categorical tick labels through the real renderer | Eighteen cell-type names | Ten overlapping pairs of glyphs at 224 x 180 and none at 720 x 380 |
| A 7 px gene name turned to run down | Chromium and Firefox, 2026-11-10 | 9 px wide in Chromium and 10.8 px in Firefox, where the names overlapped; the label's `getBBox()` reports 9.1 px in Firefox against the 10.8 px the box takes |
| Composition tests (the p-values quoted under Decisions) | Suo, liver against spleen, seven pairs | Adjusted p 0.03 (myelocytes) and 0.04 (one Langerhans cell) from seven pairs, before a type was tested over its own pairs only |

## Tried and reverted

Nothing is recorded here: every behavior replaced by a fix is in Decisions
with the reason.

## Known limits and open items

Nothing has been recorded here yet. The edge cases a contributor meets
(variable selection, visibility against highlight membership, a dataset read
again, the per-gene streaming path) are on {doc}`11_analysis_architecture`.

## Where the code and tests live

| Concern | Module(s) | Tests |
|---|---|---|
| The orchestrator and its interface | `cellucid/assets/js/app/analysis/comparison-module.js`, `cellucid/assets/js/app/analysis/ui/` | the `analysis-*` contract files of `cellucid/tests/`, the `analysis-*` specs of `cellucid/tests/browser/` |
| Reading for an analysis | `cellucid/assets/js/app/analysis/data/` | `cellucid/tests/analysis-subsample-differential.test.mjs`, `cellucid/tests/browser/analysis-no-long-tasks.spec.mjs` |
| Compute operations and statistics | `cellucid/assets/js/app/analysis/compute/`, `cellucid/assets/js/app/analysis/stats/` | `cellucid/tests/analysis-scanpy-rank-genes-parity.test.mjs`, `cellucid/tests/browser/analysis-density-workers-identity.spec.mjs` |
| Plots | `cellucid/assets/js/app/analysis/plots/`, `cellucid/assets/js/app/analysis/core/` | the `analysis-*` contract files of `cellucid/tests/` |
| Options and their schema | `cellucid/assets/js/app/analysis/options/` | `cellucid/tests/analysis-options-schema-contract.test.mjs` |
| Actions, results in the history | `cellucid/assets/js/app/analysis/actions.js`, `cellucid/assets/js/app/state/domains/analysis.js`, `cellucid/assets/js/app/analysis/shared/analysis-result-session.js`, `cellucid/assets/js/app/session/recompute-results.js` | `cellucid/tests/analysis-action-stop-contract.test.mjs`, `cellucid/tests/round-trip-analyses-contract.test.mjs`, `cellucid/tests/browser/analysis-cancel.spec.mjs` |
| Pacing of reads and tasks | `cellucid/assets/js/app/analysis/shared/run-tuner.js`, `cellucid/assets/js/app/analysis/shared/progressive-cells.js` | `cellucid/tests/browser/analysis-no-long-tasks.spec.mjs` |

## Related pages

- {doc}`11_analysis_architecture`: the reference.
- {doc}`17_extension_point_add_analysis_mode`.
- {doc}`23_decisions_partial_reading`: the gene reads an analysis makes on a
  share of the cells.
- {doc}`33_decisions_actions_domains_and_history`: results as values of the
  history, and stale results.
- {doc}`../h_analysis/index`: what a user sees.
