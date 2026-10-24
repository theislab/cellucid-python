# Analysis UI overview

**Audience:** everyone (wet lab + computational + power users)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- Where Analysis lives and what it depends on (highlight pages)
- How the mode accordion behaves (open/close; one mode at a time)
- How an analysis form is laid out (essentials, More options, Run)
- Which modes auto-run and which wait for **Run**
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

Inside the Analysis section you’ll see the analyses in three groups:

- **Genes**
  - **Marker genes** (the markers of every category of a field)
  - **Differential expression** (compare two groups of cells; volcano plot + gene table)
  - **Dot plot** (mean expression and % expressing of genes per group)
  - **Gene signature** (score gene sets across selected pages)
- **Cells**
  - **Composition** (cell-type proportions across samples and conditions)
- **Explore**
  - **Quick** (automatic composition + stats summaries)
  - **Detailed** (choose a variable + plot type + options)
  - **Correlation** (choose X and Y variables; scatter plot + stats)

```{figure} ../../../_static/screenshots/analysis/analysis-panel-modes.png
:alt: The Analysis panel in three groups: Genes (Marker genes, Differential expression, Dot plot, Gene signature), Cells (Composition) and Explore (Quick, Detailed, Correlation). Each row carries an (i), a copy-window icon and a chevron.
:width: 808px

Every analysis, collapsed. Each row shows its name, an (i) that says what it
does, a **copy** button that opens it in a floating window, and a chevron that
opens it in place.
```

Accordion behavior:
- Only one mode can be open at a time.
- Clicking an already-open mode closes it.
- Closing the mode hides results in the sidebar; reopening triggers a refresh if inputs changed.

---

## Auto-run vs “Run” button (important)

**Quick** and **Detailed** recompute as soon as their inputs are valid: selecting
pages or a variable, changing plot options, or editing a page’s cells.

Every other analysis — **Marker genes**, **Differential expression**, **Dot
plot**, **Gene signature**, **Composition** and **Correlation** — is a form that
ends in **Run**, and nothing is computed until you press it.

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

**Run** (Differential expression). Nothing is computed until the button is
pressed.
```
:::
::::

---

## The analysis form

Every analysis with **Run** lays its form out the same way, so it stays short
however many options it has:

1. **Essentials** — the two or three choices that define the analysis (Group by,
   Group A and Group B, the gene sets), each with its **(i)**. A mode that compares
   pages shows its page selector here too.
2. **More options** — one collapsed disclosure holding a collapsed group per topic
   (Cells, Test, Filters, Display, …). Each header counts the options changed
   from their defaults (`2 changed`), and **Reset to defaults** puts them back. An
   option that does not apply (the reference group while comparing against the
   rest) is hidden.
3. **Run** — with a muted line above it when the analysis cannot run, saying why
   (`Choose two pages`, `Add genes to a set`). Run is disabled until then.

Every default is stated in a **Defaults** table on the mode’s page. Which groups
you left open is remembered with the session.

**After a run.** Changing an option that would change the numbers leaves the
result in place, dimmed, with a **Stale · Run** chip: press it (or Run) to run
again. Options that only change what is shown — the filters and labels of a
volcano, the ranking and number of markers, a signature’s plot and rescaling —
redraw the result at once, without computing again. Changing pages or page
contents marks a result stale the same way.

---

## Running, cancelling and running again

Every analysis that reads genes or computes for more than a moment runs as one
**job** you can stop: Differential expression, Marker genes, Dot plot,
Composition, Correlation, Gene signature, Quick Insights and Detailed (a gene's
values, a box, violin, histogram or density plot).

**While it runs.** Run reads `Running…`, and a box under the form (in Quick
Insights and Detailed, where the result appears) shows, once the job has run for
a quarter of a second:

- the mode's name and, when it is known, what the run covers
  (`3,753 genes in Alpha (1,234 cells) vs Beta (5,678 cells)`);
- what it is doing: the step (`Preparing groups`, `Testing genes`,
  `Correcting for multiple testing`, `Reading genes`, `Scoring cells`,
  `Computing the correlation`, `Reading <variable>`, `Drawing the result`), how
  many genes are done of how many, how much it has read and how fast, and how
  long is left;
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

The same happens, without asking, when you close the analysis window or switch
to another dataset. Changing an option while a run goes does not stop it: the
result arrives and is marked **Stale · Run** if the form no longer matches it.

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

Correlation’s **X axis** and **Y axis** use the same picker (continuous obs or
gene). **Method** (Pearson or Spearman) and **Colour by** are under More options.

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
separated by splitters you can move:

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

**Moving and resizing.** Drag the header to move the view, with a mouse, a pen
or a finger. Drag any edge or corner to resize it; the bottom-right corner shows
a grip. Double-click the header to fill the window, and again to go back. The
least size is what the plot, the options and the statistics need side by side.
From the keyboard, `Tab` reaches the splitters and the bottom-right corner; the
arrow keys move them 16 px (64 px with `Shift`), and `Home` / `End` take the
corner to the least and the largest size.

The size and the splits you give an analysis’s expanded view are remembered
for that analysis: its next expanded view opens the same way. They are saved
in the session, and forgotten when another dataset opens.

On a screen 768 px wide or less the expanded view fills the screen, in one
column; the grip between the plot and the statistics moves that split.

:::{important}
**Tables, thresholds and exports exist only in the expanded view.** Nothing is
lost by working in the sidebar, but nothing is gained either — as soon as you
need a number, expand.
:::

Two consequences worth remembering:

- Changing a **plot option** in the expanded view re-draws; it does **not**
  re-run the analysis. Statistics computed by the run are fixed. This is why
  moving a differential-expression threshold is instant.
- Changing a **form control** back in the sidebar re-runs an auto-run mode the
  moment the control changes; in a mode with **Run** it marks the result
  **Stale · Run**, unless it only changes what is shown.

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

### Moving, resizing and maximizing a window

- **Move**: drag the header (mouse, pen or finger), or press `Alt` + an arrow
  key on the header or an edge (16 px; 64 px with `Shift`).
- **Resize**: drag any of the four edges or four corners. The edges are thin
  until the pointer reaches them, then show an accent line; the bottom-right
  corner always shows a grip. A size readout (`640 × 480`) shows while you
  drag. Near another window’s edge, the sidebar, the middle of the free area,
  the screen edges or the height of the content, the edge snaps and a thin
  guide line shows where.
- **Fit**: double-click the top or bottom edge to fit the height to the
  content; double-click a side edge to return to the width the window opened
  at.
- **Maximize**: the □ button in the header, or a double-click on the header,
  fills the free area beside the sidebar; the same again restores the size.
- **Keyboard**: once focus is inside a window, `Tab` also reaches its edges and
  corners. The arrow keys move the focused edge (16 px; 64 px with `Shift`),
  `Home` and `End` go to the least and the largest size, and a screen reader
  reads the size on each key.

A window is at least as wide as the sidebar’s least width and as tall as its
header and first row; the only largest size is the screen. A new window of an
analysis opens at the size you last gave one; the sizes are saved in the
session and forgotten when another dataset opens.

### On a phone

On a screen 720 px wide or less, analysis windows are bottom sheets: full
width, at 40%, 75% or all of the screen. Drag the grip above the header to
change the height (it settles on the nearest of the three), tap the grip to
step through them, or use the arrow keys on it. A tap on the header collapses
the sheet to its header. One sheet shows at a time; with several open, a tab
strip above the sheet names them and brings one to the front. Turning the phone
or widening the window gives every window back its own place and size.

---

## Sessions (what restores)

If you save/restore a session:
- the docked panel reopens the mode that was open, with every mode's settings
  and which of its More options groups were open,
- floating analysis windows reopen with their place, size, open or collapsed
  state, stacking order and settings, and new windows and expanded views keep
  the sizes remembered for each analysis,
- the shown results of Differential expression, Marker genes, Gene signature,
  Correlation, Dot plot and Composition are saved and restored as they were, without recomputing; Quick and
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
