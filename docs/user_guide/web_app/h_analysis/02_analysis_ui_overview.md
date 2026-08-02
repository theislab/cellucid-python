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

They also show progress feedback (phases / spinner / notifications).

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

## Page selection UI (“Compare pages”)

Most modes include a page selector labeled **Compare pages:**.

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

### Correlation mode

Correlation uses two selectors:
- **X Axis Variable:** (continuous obs or gene)
- **Y Axis Variable:** (continuous obs or gene)

It also includes:
- **Color by:** (optional categorical obs; otherwise points are colored by page)
- **Correlation method:** Pearson (linear) or Spearman (rank)

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
| **Options** (right) | the **Export:** row — PNG · SVG · CSV — and **PLOT OPTIONS** |
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
- floating analysis windows can be reopened with their geometry and settings,
- heavy analysis caches (notably bulk gene caches) may be restored lazily to speed up gene-heavy modes,
- results may still recompute depending on what the mode needs and what changed.

See {doc}`../l_sessions_sharing/index` for session semantics.

---

## Next steps

- {doc}`11_screenshots` (one analysis end to end, eight screenshotted steps)
- {doc}`03_analysis_mode_quick_insights` (Quick: composition + stats at a glance)
- {doc}`04_analysis_mode_detailed_analysis` (Detailed: variable + plot type + statistical tests)
- {doc}`09_exporting_analysis_results` (what exports exist and what they contain)
