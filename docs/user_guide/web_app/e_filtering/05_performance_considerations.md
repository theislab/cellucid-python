# Performance considerations

**Audience:** computational users + power users  
**Time:** 10–15 minutes  
**What you’ll learn:**
- Why filtering can lag on large datasets (what scales with `n_cells`)
- Which UI controls are “hot paths” (and how to use them safely)
- Practical workflows that stay fast for large exports

**Prerequisites:**
- A dataset loaded (ideally large enough to notice performance differences)

---

## The core cost model

Every time you change a filter, Cellucid recomputes a per-cell visibility mask
and then uploads the result to the GPU. Those are two different costs and they
scale differently.

**Recomputing the mask** is roughly **O(n_cells × n_enabled_filters)**. A cell
that is already hidden short-circuits — the mask is tested before any filter
predicate — so the constant is smaller than the formula suggests, but the
sweep is still over every cell.

**Uploading the mask** costs what actually changed, not what the dataset
contains. The upload compares row by row and sends only runs of changed rows,
up to thirty-two separate regions, or one bounding range once there are more
than that.
Hiding a handful of cells in a ten-million-cell dataset moves kilobytes, not
megabytes. A filter scattered across the whole dataset still touches every row
and still costs a full upload — that is the honest limit of a row-granular
texture, not an oversight.

This is why “small” UI actions (like moving a slider) can feel expensive on
million-cell datasets: it is the mask sweep, not the upload.

---

## Filtering makes drawing cheaper

:::{important}
**A heavily filtered view now renders faster than an unfiltered one.** With
ninety percent of the cells hidden, drawing is close to six times faster than
drawing all of them.
:::

This has not always been true, and the reason it was not is worth knowing if
you are comparing against an older recollection or an older screenshot.

Hidden points used to be shrunk to zero size in the vertex shader. That is not
enough: the driver clamps point size to a minimum of one pixel, so every hidden
point was still assembled, rasterised and shaded before being thrown away in
the fragment stage. Filtering therefore made rendering *slower*. Hidden points
are now rejected in clip space instead, before rasterisation, so they cost
essentially nothing.

The same change fixed a **correctness** bug that matters more than the speed:

:::{warning}
In the **Ultra-light (square points)** shader quality, the fragment shader has
no discard, so a hidden point still wrote depth and could **occlude visible
points behind it**. Filtering could silently remove cells you had asked to see.

If you have an older screenshot, exported figure, or memory of a filtered view
taken in Ultra-light quality, it may show **fewer** cells than the same state
shows now. The other two quality modes were unaffected and are byte-identical
before and after.
:::

---

## Which interactions are most expensive

### Continuous filtering with Live filtering ON

With `Live filtering = On`, every slider movement can trigger a recomputation.

On large datasets, “scrubbing” sliders back and forth is the most common cause of lag.

### Outlier slider dragging

Outlier filtering also recomputes visibility.

The UI throttles the expensive part to at most one recomputation every 50 ms
(about twenty per second) while you drag, and applies the exact final value
once you let go — so the readout stays smooth even when the recomputation
cannot. On large datasets it can still feel slow.

### Many stacked filters

Each additional enabled filter adds work per cell.

If you stack:
- multiple continuous ranges,
- multiple categorical hides,
- plus outlier filtering,

you can easily multiply the amount of per-cell work.

### Heavy render modes and modules reacting to visibility changes

Filtering affects more than points:
- while Smoke is active, committed visibility changes coalesce and rebuild
  density after the updated counts and controls have had one paint opportunity
- while Smoke is inactive, visibility changes only mark density dirty; the next
  Smoke entry builds it from the then-visible cells
- zero visible cells clear the current smoke volume
- edges/graphs may update which edges are visible
- legends recompute category counts (visible/available)

So even if the filter computation is fast, “everything that reacts to it” can compound the cost.

---

## Best practices (fast and safe)

### 1) Use “apply once” workflows

For continuous fields:

1) Turn **Live filtering** to `Off`
2) Drag Min/Max to the rough range you want
3) Click `FILTER` once

This avoids recomputing on every slider step.

:::{warning}
**`Live filtering` does not stay off.** The toggle is rebuilt with the legend
and returns to `On` whenever the legend re-renders — changing the active field,
changing the colormap, switching between kept views, and loading a session all
do that. It is not saved in a `.cellucid-session` either. On an
eighteen-million-cell dataset, switching palette from Viridis to Cividis is
enough to put you back on live recomputation with no other signal, so glance at
the toggle before every drag.
:::

```{figure} ../../../_static/screenshots/filtering/type-numeric-range.png
:alt: The continuous legend Filtering block with the Live filtering toggle reading Off, the hint below it reading Drag sliders then click Filter, Min and Max sliders with numeric readouts, and a black clickable FILTER button beside RESET with the mouse pointer on FILTER.
:width: 428px

The state you want on a large dataset: `Live filtering` **Off**, so `FILTER`
becomes clickable and one recomputation replaces one per slider step. While
`Live filtering` is On the button is greyed out entirely, which is the quickest
way to tell which mode you are in.
```

### 2) Change one thing at a time

When you’re exploring:
- adjust one filter,
- check the count line (“Showing X of Y points”),
- then add the next filter.

This makes it obvious which filter causes lag and which filter actually changes the visible set.

### 3) Prefer categorical filtering for coarse gating

If your goal is “remove one sample/batch/cluster”, categorical filters are usually fast and easy to reason about.

### 4) Keep multiview lean while tuning filters

If you are in Live + Snapshots:
- tune filters in the live view first (or with fewer snapshots),
- then create snapshots after you reach a stable filter state.

### 5) Use Active filters counts to detect “dead” filters

If a filter row shows `visible / available` with the same numbers, that filter is currently a no-op *given the other filters*.

Remove or disable it to reduce work.

---

## When you need to go faster than UI filtering

If you routinely need heavy gating on very large datasets:

- do coarse filtering in Python before export (e.g., remove extreme QC failures),
- export fewer fields (or precompute derived fields you actually need),
- and treat UI filtering as the *interactive refinement layer*, not the primary data cleaning step.

---

## Next steps

- {doc}`06_edge_cases_filtering` (what to expect when filters eliminate most cells)
- {doc}`07_troubleshooting_filtering` (performance-specific symptoms)
