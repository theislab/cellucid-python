# Data requirements

:::{warning}
Cross-highlighting is **planned** and **under development**.

This page documents the *data contract* required for cross-highlighting to be correct.
If the contract isn’t met, cross-highlighting must be disabled (better “no highlight” than “wrong highlight”).
:::

## At a glance

**Audience**
- Wet lab / beginner: understand “why cross-highlighting might be unavailable”.
- Computational users: know what must be present in exports/analysis results for mapping to work.
- Developers: know what metadata must be carried through analysis plotting code.

**Time**
- 10–20 minutes

**Prerequisites**
- Familiarity with the “cell index” concept:
  - {doc}`../f_highlighting_selection/01_highlight_mental_model`

**What you’ll learn**
- The mapping contract: plot element → list of cell indices
- Which plot types can satisfy the contract (and which cannot)
- Common ways the contract breaks (sampling, aggregation, dataset mismatch)

---

## The core requirement: plot elements must map to cell indices

Cross-highlighting requires that **every selectable plot element can be translated into a concrete set of cells**.

In practical terms, for an analysis plot:
- you need access to the per-cell values used to construct the plot, **and**
- a `cellIndices` array that maps those values back to dataset indices.

If an analysis result is *only aggregated* (counts per category, means per group, etc.) and the per-cell mapping is not retained, cross-highlighting cannot recover it later.

---

## What “cell indices” means (and why it matters)

Cellucid uses **dataset order** as the canonical cell identity:
- Cell `i` is “whatever row `i` is” in the loaded dataset.

This is powerful (fast and simple) but it means:
- cross-highlighting is only correct if analysis and viewer are referring to the *same* dataset ordering.

:::{important}
If you re-export/reorder your dataset, saved highlights and cross-highlight selections can become meaningless.
This is why sessions and highlights are dataset-fingerprint-bound in the design.
:::

---

## Plot type requirements (planned support)

| Plot type | Can map to cells? | What must exist |
|---|---:|---|
| Barplot | ✅ Yes | Per-cell category values + cellIndices |
| Histogram | ✅ Yes | Per-cell numeric values + cellIndices |
| Violin/Box | ✅ Yes | `cellIndices` for the distribution group |
| Pieplot | ✅ Yes | Per-cell category values + cellIndices |
| Scatterplot | ⚠️ Not initially | Requires a sampling map (see below) |
| Heatmap (aggregated) | ❌ Usually no | Heatmaps generally represent aggregate summaries, not per-cell points |

### Scatterplots: sampling breaks index mapping

Many analysis scatterplots sample to a small number of points for speed.
When this happens, `point.pointIndex` refers to a **sampled index**, not the original per-cell index.

Planned behavior (initial release):
- scatterplot cross-highlighting is disabled until a sampling map is carried through.

Reference:
- `cellucid/markdown/CROSS-HIGHLIGHTING-IMPLEMENTATION-CHECKLIST.md` (scatterplot warning)
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md` (sampling map pattern)

---

## Multi-view / multi-embedding requirements (view targeting)

Cell indices are **dimension-agnostic** (index 123 is the same cell in 2D and 3D).
But cross-highlighting still needs to decide:
- which *view* should receive the highlight (live vs snapshot, embedding A vs embedding B).

Planned requirement:
- cross-highlighting events carry a `viewId` (or follow a consistent default like “active view”).

Why this matters:
- each view can have different filters/visibility,
- snapshot views can have different positions/buffers,
- per-view isolation prevents “highlights appearing in the wrong panel”.

Reference:
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md` (render vs renderWithSnapshot + viewId propagation)

---

## If you are preparing data (cellucid-python context)

Cross-highlighting is a web app feature, but its correctness depends on **stable cell identity**.
If you use `cellucid-python` to export datasets:

- keep the dataset row order stable across exports if you expect highlights/sessions to be comparable,
- avoid silently dropping/reordering cells between “analysis in Python” and “view in browser” workflows,
- treat subsetting as creating a *new dataset identity* (new fingerprint).

If you plan to send selections/highlights between Python and the UI, treat “cell indices” as the shared ID space.

---

## Minimal checklist (for developers wiring new analysis plots)

When you add a plot that should support cross-highlighting, ensure:
- The plot’s backing `pageData` retains `cellIndices`.
- Any per-cell `values` used to compute bins/categories are retained (or computable from retained data).
- If the plot samples points, it also retains a **sampling map** from plotted points → original indices.
- The plot registers itself with the cross-highlight manager after render (planned pattern).

---

## Next steps / related pages

- {doc}`03_ux_design` (how these requirements show up as UX)
- {doc}`05_troubleshooting_cross_highlighting` (symptoms when the contract is broken)
