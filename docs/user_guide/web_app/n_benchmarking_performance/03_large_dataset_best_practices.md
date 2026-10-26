# Large dataset best practices

**Audience:** everyone (especially laptop users and teams sharing large datasets)  
**Time:** 20–45 minutes  
**What you’ll learn:**
- A “large dataset survival checklist” that keeps Cellucid responsive
- Which loading workflows are safe for large `.h5ad`/`.zarr`/`.zarr.zip` files
- How to avoid the most common performance cliffs (views, pixels, overlays, category explosion)
- What to ask your computational collaborator to change upstream when the UI can’t save you

**Prerequisites:**
- Ideally: a dataset large enough to notice performance differences
- Optional: access to server mode (for a large `.h5ad` or Zarr store on a remote machine)

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

- A large `.h5ad` or Zarr store opens *directly in the browser*, read from disk as each view needs it (an R object or a MatrixMarket folder is read once from the start instead). A cell-major (CSR) matrix is indexed by gene first; a file written by `prepare()` or `cellucid_prepare()` is gene-major and skips that.  
  See {doc}`../b_data_loading/03_browser_file_picker_tutorial` → “How a direct `.h5ad` is read”.
- For data on another machine: **Server Mode** over an SSH tunnel, so only the bytes a view needs cross the network. A served file is read by the same browser reader, so serving a local file does not make it faster.  
  See {doc}`../b_data_loading/04_server_tutorial`.

### 2) Start in the lowest-risk render configuration

- Prefer **Points** render mode while exploring (smoke is cinematic and heavy).
- Keep the browser window modest (pixels matter; retina is expensive).
- Close other GPU-heavy tabs (video calls + WebGL + Cellucid = sadness).
- Two of the biggest levers have already been pulled for you, and are worth
  checking rather than assuming: `Antialiasing` starts off
  at or above five million cells in 3D, or 200,000 cells in a 1D or 2D view, and
  the point size starts small, both chosen from the cell count. Turning antialiasing off saved 19–34% of the frame in 3D
  and about 90% in a dense 2D view, in the measurements on
  {doc}`02_performance_considerations_what_gets_slow_and_why`; it applies to the
  next frame either way.
- `Level of detail (LOD)` starts off. Turn it on whenever the frame is slow while
  zoomed out — it works at every dataset size. `Auto` draws fewer points the
  further back you pull, down to a floor of `min(2,000,000, cells ÷ 8)` points:
  never more than an 8× reduction, and never below two million on a large atlas.
  If you want fewer than that on weak hardware, set `Force LOD level` yourself:
  the forced level is deliberately not held to that floor. The exact rule, and why the reduction lands where it does, is in
  {doc}`02_performance_considerations_what_gets_slow_and_why`.

### 3) Tune in one view, then choose a layout deliberately

Cellucid keeps the live view plus any number of snapshots; every view shown is
drawn every frame.

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
- click **Filter** once.

Why: repeated recomputation is the most common CPU bottleneck.  
See {doc}`../e_filtering/05_performance_considerations`.

### 5) Treat overlays as “optional modules”

Vector field overlays can be amazing—but they are GPU-expensive.

For large datasets, start with a laptop-safe preset (low particle density, no bloom).  
See {doc}`../i_vector_field_velocity/05_performance_and_quality`.

### 6) Be intentional about categories (avoid “category explosion”)

If a categorical field has thousands–tens of thousands of unique values:
- the legend still scrolls (only the rows in view are drawn),
- but colours stop being distinguishable and counts take longer,
- and the field is cognitively overwhelming as a colouring.

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

1) Open the file directly (**Open…** or drop it); use Server Mode when it stays on another machine (ask your collaborator to run it if needed).
2) Start by coloring with a simple categorical field (clusters, sample, batch).
3) Make one snapshot only after you find an interesting view.
4) If filters feel slow: turn Live filtering off and use Filter.
5) If the app becomes choppy: clear snapshots and disable overlays.

### What to ask your computational collaborator for (copy/paste request)

Ask for one of these deliverables:

- **A server endpoint** you can open in the web app (when the data stays on their machine), or
- **A file written by `prepare()`** (gene-major, fast to read, no Python needed to open it).

Ask them to also provide (if possible):

- a “human-friendly” cluster/sample label (no 20,000-unique-ID fields as your primary categorical),
- a short list of “recommended fields” (QC + biology),
- a “lite” subset if the full dataset is too large for your laptop.

If you can only ask for one improvement: ask for a **lite version** of the dataset that is guaranteed to run smoothly on a normal laptop.

---

## Computational workflow (how to write datasets that stay fast)

This section is intentionally practical and non-dogmatic: it’s about what actually helps users.

### Principle 1 — Write two versions when the dataset is huge

For million-cell datasets, a single “one size fits all” file often fails.

Write two files with `prepare()`, each with its own `dataset_id`:

1) **Lite** (interactive-first)  
   - fewer fields (only what people actually use),
   - simplified categories (collapsed labels),
   - conservative quality defaults.

2) **Full** (everything, slower)  
   - complete field coverage,
   - full-resolution expression,
   - intended for power users and final figures.

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
- write fewer “junk” fields (`prepare(obs_keys=...)` lists only those that matter),
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

- **Multiview without Level of detail**: every view draws every point again,
  so on a large dataset each kept view adds about one view's frame time (four
  views of 10,000,000 cells: 30 → 147 ms a frame). With Level of detail on, the
  same four views took 13.5 ms.
- **High-DPI + large window**: smooth on an external monitor, choppy on a laptop retina display.
- **Smoke mode on a laptop**: looks great, then crashes (context lost) when render resolution and ray quality push one frame past the GPU watchdog.
- **Overlay across views**: overlays scale with pixels × views in every layout.
- **Huge highlight groups**: highlighting costs GPU memory per highlighted cell
  per view, and it does not go away when you look elsewhere.
- **Category explosion**: the field loads and scrolls, but colours stop being distinguishable.

For a deeper list, see {doc}`06_edge_cases_performance`.

---

## Mini-troubleshooting (fast)

If large datasets are unusable, triage in this order:

1) Confirm WebGL2 + hardware acceleration: {doc}`../a_orientation/02_system_requirements`.
2) Reduce pixels: make the window smaller and retry.
3) Clear snapshots and retry.
4) Disable smoke mode and overlays and retry.
5) Change filtering workflow: Live filtering off → Filter once.
6) Clear large highlight groups and retry.
7) If still slow: move the data to a local disk (or serve it from the machine that holds it) and retry.
8) If still slow: treat it as a dataset issue; go to {doc}`09_reporting_performance_bugs`.

---

## Interface reference

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: The whole window in Grid compare with the sidebar open: two equal panels to the right of the sidebar show the Pancreas 3D embedding from the same linked camera — view 1 coloured by cell_type, view 2 recoloured by S_score on a Viridis scale — each under a title chip naming its embedding and field. The sidebar shows Compare views with the two view rows, each marked linked and Orbit, View 2 shows UMAP · 3D, Keep view, Link cameras, Layout Grid compare and Clear, and Coloring & filtering with S_score in the Continuous obs select.
:width: 1440px

Two kept views support side-by-side comparison while retaining one dataset
identity. Both panels are the full canvas height and split the width to the
right of the sidebar equally: this row layout keeps points at full size, which
is why it is the costliest arrangement per panel.
```

---

## Next steps

- {doc}`04_benchmarking_methodology_and_metrics` (measure the impact of a change reliably)
- {doc}`06_edge_cases_performance` (catalog of performance cliffs)
- {doc}`07_troubleshooting_performance` (symptom → diagnosis → fix)
