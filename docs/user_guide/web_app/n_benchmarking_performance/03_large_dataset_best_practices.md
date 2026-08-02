# Large dataset best practices

**Audience:** everyone (especially laptop users and teams sharing large exports)  
**Time:** 20–45 minutes  
**What you’ll learn:**
- A “large dataset survival checklist” that keeps Cellucid responsive
- Which loading workflows are safe for large `.h5ad`/`.zarr`/exports
- How to avoid the most common performance cliffs (views, pixels, overlays, category explosion)
- What to ask your computational collaborator to change upstream when the UI can’t save you

**Prerequisites:**
- Ideally: a dataset large enough to notice performance differences
- Optional: access to server mode (recommended for large `.h5ad`/`.zarr`)

---

## First: what do we mean by “large”?

“Large” is not a single number—it depends on:
- your GPU (VRAM and driver quality),
- your window size and device pixel ratio,
- whether you use overlays/smoke,
- and how many views you keep open.

Practical (very rough) intuition:
- **< 100k cells**: usually fine on most modern machines.
- **100k–500k**: you need good habits (avoid recompute loops, keep the window modest).
- **500k–2M+**: drawing is usually still fine; what bites first is the *CPU* cost
  of each filter change and the memory footprint of large highlight groups.
- **Beyond that**: the built-in benchmark offers presets up to 20 million points
  precisely because the renderer is not the first thing to give way. Run it on
  your own machine before assuming a dataset is too big — the answer is a
  property of your GPU, not of the cell count.

If you hit “WebGL context lost”, jump to {doc}`07_troubleshooting_performance` and {doc}`../a_orientation/02_system_requirements`.

:::{tip}
The single most useful habit on a large dataset is to **filter early**. Cells
that filters have hidden are rejected before they are drawn, so a working view
that shows the population you care about costs less to navigate than the full
atlas — the reverse of what most people expect.
:::

---

## Fast path: large-dataset survival checklist (copy/paste mental defaults)

If you only do one thing from this page, do this list in order.

### 1) Use the right loading workflow

- For large `.h5ad` and very large `.zarr`: **Server Mode** is the default recommendation.  
  See {doc}`../b_data_loading/04_server_tutorial`.
- Loading a large `.h5ad` *directly in the browser* is usually a trap because it is not truly lazy.  
  See {doc}`../b_data_loading/03_browser_file_picker_tutorial` → “Important performance limitation”.

### 2) Start in the lowest-risk render configuration

- Prefer **Points** render mode while exploring (smoke is cinematic and heavy).
- Keep the browser window modest (pixels matter; retina is expensive).
- Close other GPU-heavy tabs (video calls + WebGL + Cellucid = sadness).

### 3) Tune in one view, then choose a layout deliberately

Cellucid gives you the live view plus at most three snapshots.

Workflow:
1) Tune fields/filters in the live view — it is unambiguously the cheapest state
   to iterate in.
2) Only then add snapshots for the comparison you actually want.

Do not assume “fewer views = faster”. Two and three views are laid out as a row
and keep full-size points; four are laid out as a 2×2 grid with half-size points
and often cost about the same as one view. See {doc}`06_edge_cases_performance`.

### 4) Filter like a grown-up (avoid recompute loops)

On large datasets:
- turn **Live filtering** off for continuous fields,
- move the sliders,
- click **FILTER** once.

Why: repeated recomputation is the most common CPU bottleneck.  
See {doc}`../e_filtering/05_performance_considerations`.

### 5) Treat overlays as “optional modules”

Vector field overlays can be amazing—but they are GPU-expensive.

For large datasets, start with a laptop-safe preset (low particle density, no bloom).  
See {doc}`../i_vector_field_velocity/05_performance_and_quality`.

### 6) Be intentional about categories (avoid “category explosion”)

If a categorical field has thousands–tens of thousands of unique values:
- the legend becomes unusable,
- color assignment and counts can become a performance cliff,
- and the UI becomes cognitively overwhelming.

If you see this, treat it as a **data preparation problem**, not a “tweak a knob” problem.

### 7) Separate “exploration” from “publication/export”

Exporting publication figures can be expensive (large output sizes, huge SVGs).

Best practice:
- explore with performance-safe settings,
- then switch to an export-focused workflow once your view is final.

### 8) Save sessions intentionally

On large datasets, sessions can become large if you:
- keep many views,
- create huge highlight groups,
- or accumulate many artifacts.

When collaborating, prefer:
- fewer views,
- fewer huge highlights,
- and documented “what this session represents”.

See {doc}`../l_sessions_sharing/index`.

---

## Wet-lab friendly workflow (what to do, what to ask for)

This is written for the situation:
“I just want to explore the data, but I’m not the one exporting it.”

### What you can do in the UI (no code)

1) Use Server Mode if you have a `.h5ad` or `.zarr` (ask your collaborator to run it if needed).
2) Start by coloring with a simple categorical field (clusters, sample, batch).
3) Make one snapshot only after you find an interesting view.
4) If filters feel slow: turn Live filtering off and use FILTER.
5) If the app becomes choppy: clear snapshots and disable overlays.

### What to ask your computational collaborator for (copy/paste request)

Ask for one of these deliverables:

- **A server endpoint** you can open in the web app (recommended for very large data), or
- **A pre-exported folder** (fast, reliable, no Python needed after export).

Ask them to also provide (if possible):

- a “human-friendly” cluster/sample label (no 20,000-unique-ID fields as your primary categorical),
- a short list of “recommended fields” (QC + biology),
- a “lite export” if the full dataset is too large for your laptop.

If you can only ask for one improvement: ask for a **lite export** that is guaranteed to run smoothly on a normal laptop.

---

## Computational workflow (how to make exports that stay fast)

This section is intentionally practical and non-dogmatic: it’s about what actually helps users.

### Principle 1 — Make two exports when the dataset is huge

For million-cell datasets, a single “one size fits all” export often fails.

Create two bundles:

1) **Lite** (interactive-first)  
   - fewer fields (only what people actually use),
   - simplified categories (collapsed labels),
   - conservative quality defaults.

2) **Full** (everything, slower)  
   - complete field coverage,
   - full-resolution expression,
   - intended for power users and final exports.

This is often the fastest way to reduce “it’s slow” support load while still preserving access to the full dataset.

### Principle 2 — Treat category explosion as a correctness issue

Even if the app could render it, a categorical field with tens of thousands of unique values is rarely usable.

Best practice:
- collapse rare categories into “Other”,
- use stable, human-readable labels (avoid raw UUIDs),
- keep the “primary” categorical fields to a manageable size.

### Principle 3 — Do coarse filtering upstream

UI filtering is great for interactive refinement, but it is not a substitute for data cleaning.

If users consistently apply the same heavy gate:
- remove obvious QC failures upstream,
- export fewer “junk” fields,
- and keep the UI workflows focused on discovery.

### Principle 4 — Document the intended workflow

Many performance problems are really “workflow mismatch” problems.

If your team expects:
- “use server mode, then open in browser”,
- “always start in points mode”,
- “keep 2–4 views”,

write it down and link it in your lab’s internal handoff doc.

---

## Performance patterns for very large data (how to stay productive)

### Pattern A — Coarse-to-fine exploration

1) Start with a coarse categorical field (clusters, sample).
2) Use a *small number* of filters to isolate a region.
3) Only then inspect gene expression or run heavier analysis.

Why this works:
- early steps are mostly GPU drawing + light UI,
- you avoid repeatedly recomputing on “the entire dataset”.

### Pattern B — “One view while tuning, snapshots only for comparison”

Use the live view to tune settings, then snapshot for the final state you want to compare.

Why this works:
- avoids multiplying expensive updates across many views.

### Pattern C — Export artifacts instead of pushing the UI further

If a task is inherently heavy (huge DE, massive export, complex overlays):
- do it once,
- export the result (tables/figures),
- and iterate offline or with smaller subsets.

This is often more reproducible than doing many incremental heavy operations in the UI.

---

## Common “performance cliffs” to avoid

These are so common that it’s worth recognizing them immediately.

- **The two-view row**: adding one snapshot doubles the pixels shaded, while
  adding two more (reaching the 2×2 grid) usually does not.
- **High-DPI + large window**: smooth on an external monitor, choppy on a laptop retina display.
- **Smoke mode on a laptop**: looks great, then crashes (context lost) when grid density is too high.
- **Overlay across views**: overlays scale with pixels × views in every layout.
- **Huge highlight groups**: highlighting costs GPU memory per highlighted cell
  per view, and it does not go away when you look elsewhere.
- **Category explosion**: the UI becomes slow and unusable even if rendering is fine.

For a deeper list, see {doc}`06_edge_cases_performance`.

---

## Mini-troubleshooting (fast)

If large datasets are unusable, triage in this order:

1) Confirm WebGL2 + hardware acceleration: {doc}`../a_orientation/02_system_requirements`.
2) Reduce pixels: make the window smaller and retry.
3) Clear snapshots and retry.
4) Disable smoke mode and overlays and retry.
5) Change filtering workflow: Live filtering off → FILTER once.
6) Clear large highlight groups and retry.
7) If still slow: switch loading method (server mode) and retry.
8) If still slow: treat it as a dataset/export issue; go to {doc}`09_reporting_performance_bugs`.

---

## Interface reference

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: Cellucid showing two side-by-side views of the same dataset.
:width: 1440px

Two kept views support side-by-side comparison while retaining one dataset
identity. Note that both panes are full canvas height: this row layout keeps
points at full size, which is why it is the costliest arrangement per pane.
```

---

## Next steps

- {doc}`04_benchmarking_methodology_and_metrics` (measure the impact of a change reliably)
- {doc}`06_edge_cases_performance` (catalog of performance cliffs)
- {doc}`07_troubleshooting_performance` (symptom → diagnosis → fix)
