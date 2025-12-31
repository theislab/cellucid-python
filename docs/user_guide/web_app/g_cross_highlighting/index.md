# Cross-highlighting (analysis ↔ embedding; planned)

:::{warning}
Cross-highlighting is **planned** and **under development**.

This section documents the intended user experience and the implementation plan.
Depending on the app build you are using, some or all of the described UI may be missing or non-functional.
:::

Cross-highlighting is the feature that lets you **jump between “summaries” and “where those cells live”**:

- In the **analysis window**, you see plots like barplots, histograms, box/violin plots, etc.
- In the **embedding viewer**, you see the same cells as points.
- Cross-highlighting connects them: click/hover a plot element → highlight the corresponding cells in the embedding.

This feature is especially useful when you want to answer questions like:
- “Which cells contribute to this bar/bin?”
- “Where are the outliers that dominate the tail?”
- “If this group has a marker, are they spatially localized in UMAP/3D?”

## At a glance

**Audience**
- Wet lab / beginner: “click a plot → see cells → optionally save as a selection”.
- Computational users: understand the *mapping contract* (cell indices), limitations (e.g., sampling), and how filters affect what you see.
- Developers: understand the wiring plan and the current gaps.

**Time**
- 10–25 minutes

**Prerequisites**
- A dataset loaded in the web app
- Familiarity with:
  - {doc}`../h_analysis/index` (where the plots come from)
  - {doc}`../f_highlighting_selection/index` (what “highlight pages/groups” mean)

**What you’ll learn**
- What cross-highlighting *is* (and what it is not)
- What data must exist for plot → cell mapping
- What UX is planned (hover preview, click selection, “Save as Page”)
- Known limitations and how to debug mismatches

---

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```

---

## Fast path (pick your goal)

| You want to… | Start here |
|---|---|
| Understand the idea and when it’s useful | {doc}`01_what_cross_highlighting_is_user_story` |
| Know why cross-highlighting can’t “find the right cells” | {doc}`02_data_requirements` |
| See intended UI behavior (what you should expect to happen) | {doc}`03_ux_design` |
| Avoid performance pitfalls (large datasets, many views) | {doc}`04_performance_correctness_notes` |
| Diagnose “nothing happens” / “wrong cells highlighted” | {doc}`05_troubleshooting_cross_highlighting` |
| Contribute / wire it up (developer notes) | {doc}`06_reference_implementation_notes` |

---

## Status and roadmap (high level)

The current cross-highlighting code exists in the web app, but **the integration is not complete**, so the feature is considered **planned / under development**.

The implementation plan is tracked in:
- `cellucid/markdown/CROSS-HIGHLIGHTING-FIX-PLAN.md`
- `cellucid/markdown/CROSS-HIGHLIGHTING-IMPLEMENTATION-CHECKLIST.md`
- `cellucid/markdown/CROSS-HIGHLIGHTING-NEW-MODULES.md`
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md`

In one page, the plan is:
1) **Core wiring (critical)**: add missing viewer APIs and connect the manager to the viewer lifecycle.
2) **State integration (high)**: route highlights through the DataState highlight/preview system (so behavior matches the rest of the app).
3) **Hover + plot coverage (medium)**: make hover previews consistent across plot types (with safe fallbacks).
4) **Multi-view correctness (medium)**: ensure highlights target the correct view/snapshot and respect per-view visibility.

If you are using Cellucid via `cellucid-python` (Jupyter/CLI), note that cross-highlighting itself is a **web app feature**.
`cellucid-python` can eventually drive/receive highlight events via hooks, but the cross-highlighting wiring described here lives in the web app repo (`cellucid/`).
