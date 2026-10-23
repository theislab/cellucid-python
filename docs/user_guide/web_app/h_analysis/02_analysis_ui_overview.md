# Analysis UI overview

**Audience:** everyone (wet lab + computational + power users)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- Where Analysis lives and what it depends on (highlight pages)
- How the mode accordion behaves (open/close; one mode at a time)
- Which modes auto-run vs require a “Run” button
- How “Copy” creates floating analysis windows (what is copied vs recomputed)
- Where exports live (CSV) and where plot options live (modal/expanded view)

**Prerequisites:**
- A dataset loaded in the web app
- At least one highlight page (recommended for nearly all modes)

---

## Where Analysis lives (UI map)

In the left sidebar, Analysis is a dedicated accordion section:

- **Sidebar → Analysis** (directly below the Highlighted Cells / highlight pages UI).

Analysis is intentionally placed right after highlighting because:
- you define groups with highlight pages, and
- you analyze/compare those pages here.

If Analysis looks empty, the first thing to check is whether you have highlight pages:
see {doc}`../f_highlighting_selection/index`.

---

## The Analysis mode accordion (how it behaves)

Inside the Analysis section you’ll see a mode accordion with these items:

- **Quick** (automatic composition + stats summaries)
- **Detailed** (choose a variable + plot type + options)
- **Correlation** (choose X and Y variables; scatter plot + stats)
- **Differential Expression** (compare two pages; volcano plot + gene table)
- **Gene Signature** (paste gene list; score across selected pages)
- **Marker Genes** (discover one-vs-rest markers for a categorical obs field)

```{figure} ../../../_static/screenshots/analysis/analysis-panel-modes.png
:alt: Six stacked rows reading Quick (Automatic insights), Detailed (Full control over options), Correlation (Explore variable relationships), Differential Expression (Find DE genes between groups), Gene Signature (Compute signature scores) and Marker Genes (Discover markers across groups). Each row carries a small copy-window icon and a chevron.
:width: 808px

All six modes, collapsed. Each row shows the mode name, a one-line description, a
**copy** button that undocks the mode into a floating window, and a chevron that
opens it in place. In the sidebar’s real width the descriptions truncate; this
capture widens the sidebar so all six read in full.
```

Accordion behavior:
- Only one mode can be open at a time.
- Clicking an already-open mode closes it.
- Closing the mode hides results in the sidebar; reopening triggers a refresh if inputs changed.

---

## Auto-run vs “Run” button (important)

Cellucid uses two interaction styles:

### Auto-run modes (no “Run” button)

These recompute automatically when inputs become valid:
- **Quick**
- **Detailed**
- **Correlation**
- **Gene Signature**

Typical trigger events:
- selecting pages,
- selecting variables,
- changing plot options,
- adding/removing cells from pages.

### Explicit-run modes (you click a button)

These are potentially expensive and require an explicit action:
- **Differential Expression** → `Run Differential Expression`
- **Marker Genes** → `Discover Markers`

Under the button they say what a run will read before you press it, and while
it runs they say what it is doing, with **Cancel** (see
[Running, cancelling and running again](#running-cancelling-and-running-again)).

The difference is visible in the form itself: an auto-run mode ends at its last
control, and an explicit-run mode ends at a button.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/autorun-vs-run-a.png
:alt: The Detailed panel showing a Variable select, a Compare pages list, a Plot Type select and a bar plot underneath, with no run button anywhere between them.
:width: 488px

**Auto-run** (Detailed). The result sits directly under the controls. Change a
control and it redraws.
```
:::

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/autorun-vs-run-b.png
:alt: The Differential Expression panel showing Page A and Page B selects, a Statistical method select, a collapsed Performance Settings row, and a full-width Run Differential Expression button at the bottom with a pointer on it. No result is shown.
:width: 488px

**Explicit run** (Differential Expression). Nothing is computed until the button
is pressed, and changing a control afterwards cancels the run in flight.
```
:::
::::

---

## Running, cancelling and running again

Every analysis that reads genes or computes for more than a moment runs as one
**job** you can stop: Differential Expression, Marker Genes (and the heatmap
genes it reads afterwards), Correlation, Gene Signature, Quick Insights and
Detailed (a gene's values, a box, violin, histogram or density plot).

**Before it runs.** An explicit-run mode says what a run will read under its
button, from what is chosen now: Differential Expression
`Tests 3,753 genes in 1,234 vs 5,678 cells.`, Marker Genes
`Reads 3,753 genes in 3,696 cells, testing each category of celltype against
the rest.`

**While it runs.** A box under the form (in Quick Insights and Detailed, where
the result appears) shows, once the job has run for a quarter of a second:

- the mode's name and, when it is known, what the run covers
  (`3,753 genes in Alpha (1,234 cells) vs Beta (5,678 cells)`);
- what it is doing: the step (`Testing genes`, `Reading genes`,
  `Discovering markers`, `Building expression matrix`, `Clustering`,
  `Reading <variable>`, `Drawing the plot`), how many genes are done of how
  many, how much it has read and how fast, and how long is left;
- a progress bar, and **Cancel**.

A job that takes longer than 0.7 s also shows one progress notification with
the same line and its own cancel button (×). It is the run's only
notification: the gene reads it makes are not announced again by the general
read notifications -- neither is a gene a view started reading before the run
needed it too, which the run takes as its own read from then on -- and it
closes with the outcome (`Differential expression: 412 significant genes`).

**Cancel.** Cancelling stops the job where it is: its gene reads stop in the
file reader (a remote file stops downloading), and its computations stop in
the analysis workers (a worker still busy a quarter of a second after the
cancel is replaced by a fresh one, so nothing keeps running for a cancelled
job). Nothing of a cancelled run is kept: no partial result is shown, cached
or exported. A read that something else also waits for, the gene the
embedding is coloured by, or the same gene another analysis window reads,
goes on for it.

The same happens, without asking, when you close the analysis window, switch
to another dataset, start the same analysis again, or change a control of a
running explicit-run mode.

**After a cancel** the box reads `<mode> cancelled` and
`Cancelled after 12 s. 1,240 of 3,753 genes were done. Nothing was kept.`,
with **Run again**. A cancel is never reported as an error. **Run again** runs
the form as it stands.

**When a run fails** the box reads `<mode> failed` with the reason, again with
**Run again**, and the notification says the same once.

---

## Page selection UI (“Compare pages”)

Most modes include a page selector labeled **Compare pages**.

What you can do there:
- click tabs to select which pages are included,
- see per-page cell counts,
- customize page colors (used consistently in plots),
- select all/deselect all quickly.

Some modes also show derived pages:
- **Rest of \<page\>** (the complement of a page across the full dataset).

```{figure} ../../../_static/screenshots/analysis/rest-of-page-selector.png
:alt: A "Compare pages:" block with the help line "Click pages to select. Click color to customize." and two chips: "Page 1" with a colour swatch and the count 262, and "Rest of Page 1" with its own colour swatch and the count 3.4K. A pointer rests on the second chip.
:width: 448px

`Rest of Page 1` is the derived complement — every cell **not** in Page 1. You do
not create it; it appears next to each base page with its own count and colour,
which is what makes a two-group comparison possible from a single selection.
```

:::{tip}
When a plot “doesn’t match what you expect”, the most common cause is selecting the wrong pages.
Always confirm which pages are selected in the mode you’re currently using.
:::

---

## Variable selection UI (obs vs genes)

### Detailed mode

Detailed uses a two-step variable selector:
1) choose the **type**: Categorical obs / Continuous obs / Gene expression
2) choose the **specific variable** (field key or gene)

If gene expression is not available in the dataset, the selector will show an explicit “No gene expression data available” message.

The gene list is the same windowed list as the sidebar's gene search (see
{doc}`../d_fields_coloring_legends/02_field_selector_ux`): it lists **None**
and every matching gene, opens with the chosen gene highlighted, and only
builds the rows near what is shown. `↓` / `↑` move the highlight and wrap at
either end, **Home** and **End** jump to the first and the last entry,
`Enter` chooses the highlighted entry and `Esc` closes the list.

A gene that stays under the pointer, or highlighted with the arrow keys, in the
gene list for 150 ms starts loading in the background, so choosing it usually
needs no wait. Moving on cancels that load, and choosing the gene uses it.

### Correlation mode

Correlation uses two selectors:
- **X axis variable**: (continuous obs or gene)
- **Y axis variable**: (continuous obs or gene)

It also includes:
- **Color by** (optional categorical obs; otherwise points are colored by page)
- **Correlation method**: Pearson (linear) or Spearman (rank)

---

## Sidebar preview vs expanded view

Every result-producing mode renders **twice**, into two very different surfaces.
Knowing which one you are looking at explains most “where is that number?”
questions.

### The sidebar preview

The preview lives inside the mode’s own accordion panel, roughly 224 pixels wide.
It contains **the plot and nothing else** — no tables, no thresholds, no export.
Its job is to confirm that the run produced something and to give you a shape to
recognise.

Plots trade away decoration rather than legibility at that width. A category
axis prints only the names that fit and drops the rest whole, so the labels you
can see are real categories rather than truncated ones; the volcano keeps its
axis titles and drops its legend, because an unlabelled axis cannot be read at
any size while a red/blue/grey volcano still reads without a key. Both come back
as the box grows, so expand when you need the names.

### The expanded view

Below every preview is an **⤢ Expand** button.

```{figure} ../../../_static/screenshots/analysis/expand-button.png
:alt: A right-aligned pill button reading "⤢ Expand" beneath the bottom edge of a preview plot, with a pointer resting on it.
:width: 472px

The **⤢ Expand** button is the only way into the full view. The preview above it
is not clickable — it is a plot, not a control.
```

The button is a real `<button>`: it takes focus with `Tab`, is announced as
“Expand”, and opens on `Enter` or `Space`. Its tooltip reads *“Open in full view
with statistics and export options”*.

The expanded view is a modal over the whole window, laid out in four regions
separated by drag handles you can move:

| Region | Holds |
| --- | --- |
| **Plot** (top left) | the same plot, with room to be legible |
| **Options** (right) | the **Export:** row — PNG · SVG · CSV — and **Plot options** |
| **Summary statistics** (bottom left) | the numeric table for this mode |
| **Statistical analysis** (bottom right) | tests, interpretation, ranked tables |

```{figure} ../../../_static/screenshots/analysis/detailed-expanded.png
:alt: A modal divided into four regions: a grouped bar plot top-left, an Export row and PLOT OPTIONS column on the right, a SUMMARY STATISTICS table bottom-left, and a STATISTICAL ANALYSIS panel bottom-right showing a chi-squared card that reads N/A.
:width: 1440px

The four-region expanded layout, here for **Detailed**. Every mode with an Expand
button uses this same layout, so learning it once is enough.
```

Close it with the `×` in its header or with `Escape`.

:::{important}
**Tables, thresholds and exports exist only in the expanded view.** Nothing is
lost by working in the sidebar, but nothing is gained either — as soon as you
need a number, expand.
:::

Two consequences worth remembering:

- Changing a **plot option** in the expanded view re-draws; it does **not**
  re-run the analysis. Statistics computed by the run are fixed. This is why
  moving a differential-expression threshold is instant.
- Changing a **form control** back in the sidebar does re-run, and for an
  auto-run mode that happens the moment the control changes.

---

## “Copy” (floating analysis windows)

Each analysis accordion item has a small **Copy** control.

Copy creates a floating analysis window so you can:
- keep one result visible while exploring another mode, or
- run the same mode with different settings side-by-side.

```{figure} ../../../_static/screenshots/analysis/floating-analysis-window.png
:alt: The full application window with the sidebar on the left and a small floating panel headed "Correlation" placed over the point cloud, carrying its own copy and close buttons and its own copy of the correlation form.
:width: 1440px

**Copy** undocks a mode into a floating window that sits over the viewer. It has
its own header, its own copy and close buttons, and its own independent copy of
the form — so you can keep one configuration on screen while changing another.
```

Important behavior:
- **Settings are copied. Results are not copied.**  
  The copied window recomputes from current data/pages.
- Floating windows have their own UI instance, so they can be configured independently.
- There is a safety cap on the number of open analysis windows (to prevent memory explosions).

---

## Sessions (what restores)

If you save/restore a session:
- the docked panel reopens the mode that was open, with every mode's settings,
- floating analysis windows reopen with their geometry and settings,
- the shown results of Differential expression, Marker genes, Gene signature and
  Correlation are saved and restored as they were, without recomputing; Quick and
  Detailed recompute from their settings,
- a saved result that names a gene, field or category the open dataset lacks is
  refused rather than shown against different data,
- gene-expression caches are not saved; genes are read again when needed.

See {doc}`../l_sessions_sharing/index` for session semantics.

---

## Next steps

- {doc}`11_screenshots` (one analysis end to end, eight screenshotted steps)
- {doc}`03_analysis_mode_quick_insights` (Quick: composition + stats at a glance)
- {doc}`04_analysis_mode_detailed_analysis` (Detailed: variable + plot type + statistical tests)
- {doc}`09_exporting_analysis_results` (what exports exist and what they contain)
