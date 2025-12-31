# What cross-highlighting is (user story)

:::{warning}
Cross-highlighting is **planned** and **under development**.
This page describes the intended user experience and the mental model so you can recognize when it becomes available (and file high-quality bugs while it’s being built).
:::

## At a glance

**Audience**
- Wet lab / beginner: learn “click a plot → see those cells in the embedding”.
- Computational users: learn what is actually being mapped (cell indices) and why some plot types cannot map reliably.
- Developers: learn what the UI is trying to achieve so the implementation matches user expectations.

**Time**
- 5–15 minutes

**Prerequisites**
- Familiarity with:
  - {doc}`../h_analysis/index` (analysis plots)
  - {doc}`../f_highlighting_selection/index` (highlights/pages, and “preview highlight”)

**What you’ll learn**
- The core user story (why this feature exists)
- What happens on hover vs click (planned)
- How cross-highlighting relates to selection tools and highlight pages

---

## What this feature is (one sentence)

Cross-highlighting lets you **pick a subset of cells from an analysis plot** (a bar/bin/distribution) and immediately **see where those same cells are in the embedding viewer**, with an option to save that set as a reusable highlight page/group.

---

## The user story (wet lab / non-technical)

You ran an analysis (e.g., “compare conditions”, “marker genes”, “QC distributions”) and you see a plot.

Now you want to answer:
- “Where are the cells that make this bar high?”
- “Are these outliers all from one cluster?”
- “If this gene is high in this group, do those cells sit in one region of UMAP?”

**Planned workflow**
1) In the analysis window, hover a bar/bin → preview-highlight those cells in the embedding.
2) Click the bar/bin → select that set (still non-destructive).
3) A small notification appears: “N cells selected” with actions:
   - **Save as Page** (turn it into a highlight page or highlight group)
   - **Clear**
4) Switch back to the embedding view and continue exploring, filtering, or running more analysis.

:::{note}
Cross-highlighting is not a replacement for lasso/KNN selection. It’s a faster way to “pull cells out of a plot”.
:::

---

## Practical path (computational users)

### What is being mapped (the contract)

Everything depends on a simple contract:

- The app stores and moves around sets of **cell indices**: integers in `0..(n_cells-1)` for the loaded dataset order.
- An analysis plot element must be able to map “this bar/bin/point” → “these cell indices”.

If that mapping is missing or ambiguous, cross-highlighting must be disabled (or it will highlight the wrong cells).

### What it does *not* do

Cross-highlighting does not:
- change your filters (visibility),
- permanently modify highlight pages/groups automatically,
- change the dataset ordering,
- “know biology” (it is a UI mapping feature, not an analysis method).

### Filters can hide the highlight

Cross-highlighting can “work” but look like it didn’t if the target view is filtering those cells out.
You may see a notification like “3,000 cells selected” but 0 visible highlighted points.

Common fix: temporarily disable filters in the target view, or switch to a view where those cells are visible.

### Scatterplots are special (sampling breaks mapping)

Many analysis scatterplots are sampled for speed (e.g., max ~500 points). If the plot only contains a sampled subset, `pointIndex` no longer matches the original dataset index.

Planned behavior (initial release):
- scatterplot cross-highlighting is **disabled** rather than wrong.

---

## Deep path (experts / developers)

### How this fits with highlights and pages

The intended design is:
- Hover: uses **preview highlight** (transient).
- Click: sets a **cross-highlight selection** and shows a toast.
- “Save as Page”: converts the selection into the existing highlight system (page/group) so it can persist and be used for analysis and sessions.

This is why the plan emphasizes routing through the same DataState highlight pipeline used by manual selection tools.

### Where the plan lives

Start here:
- `cellucid/markdown/CROSS-HIGHLIGHTING-FIX-PLAN.md` (root causes + phases)
- `cellucid/markdown/CROSS-HIGHLIGHTING-IMPLEMENTATION-CHECKLIST.md` (step-by-step checklist)
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md` (viewId, snapshots, dimensionality pitfalls)

---

## Next steps / related pages

- {doc}`02_data_requirements`
- {doc}`03_ux_design`
- {doc}`../f_highlighting_selection/01_highlight_mental_model`
- {doc}`../h_analysis/index`
