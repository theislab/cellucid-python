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

Every time you change a filter, Cellucid recomputes a per-cell visibility mask.

High-level scaling:
- **Visibility recomputation is roughly O(n_cells × n_enabled_filters)**.

This is why “small” UI actions (like moving a slider) can feel expensive on million-cell datasets.

---

## Which interactions are most expensive

### Continuous filtering with Live filtering ON

With `Live filtering = On`, every slider movement can trigger a recomputation.

On large datasets, “scrubbing” sliders back and forth is the most common cause of lag.

### Outlier slider dragging

Outlier filtering also recomputes visibility.

The UI throttles updates while dragging, but on large datasets it can still feel slow.

### Many stacked filters

Each additional enabled filter adds work per cell.

If you stack:
- multiple continuous ranges,
- multiple categorical hides,
- plus outlier filtering,

you can easily multiply the amount of per-cell work.

### Heavy render modes and modules reacting to visibility changes

Filtering affects more than points:
- smoke/volumetric modes may rebuild density from visible points
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

- `06_edge_cases_filtering` (what to expect when filters eliminate most cells)
- `07_troubleshooting_filtering` (performance-specific symptoms)
