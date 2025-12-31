# Performance mental model

**Audience:** everyone (wet lab → computational → developer)  
**Time:** 10–25 minutes  
**What you’ll learn:**
- The three performance bottlenecks in Cellucid (GPU, CPU, I/O) and how they “feel”
- The multipliers that turn “fine” into “unusable” (cells, views, categories, pixels)
- A reliable triage workflow: identify → confirm → change one knob → measure
- What’s “normal slow” vs “something is wrong”

**Prerequisites:**
- None (recommended: have any dataset loaded so you can try the checks)

---

## The one-sentence model

Cellucid stays fast when it can keep **(1) per-cell computation**, **(2) GPU rendering**, and **(3) data loading** within your machine’s limits; “performance problems” happen when one of these becomes the bottleneck.

If you can identify which of the three is limiting you, the fix is usually obvious.

---

## Three explanations (pick your depth)

::::{tab-set}

:::{tab-item} Wet‑lab / Non‑technical

Think of Cellucid like a microscope with three “power supplies”:

- **Graphics (GPU)**: draws the points smoothly when you move/zoom.
- **Computation (CPU)**: recalculates which cells are visible when you filter or run analysis.
- **Loading (I/O)**: reads data from your disk/server (especially when you switch genes or load a big dataset).

When Cellucid is slow, it’s usually because **one** of these is maxed out. Your goal is to do a quick check to figure out which one.

:::

:::{tab-item} Computational / Power user

Cellucid performance is dominated by:

- **O(n_cells)** work (filter masks, selections/highlights, some aggregations),
- **GPU fill + draw cost** (points, post-processing, overlays, smoke raymarch),
- **I/O and caching** (remote data, expression fetches, “first-time” loads).

The practical trick: many actions you’d consider “UI-only” trigger large buffer updates or recomputation; *avoid repeated recompute loops* and *avoid multiplying work by opening many views*.

:::

:::{tab-item} Developer / Deep

At a high level:

- The viewer maintains large arrays (positions, colors, visibility, selection, sometimes per-view buffers).
- Many UI actions update these arrays and then **upload to GPU** or trigger dependent recompute.
- Some modules add full-screen passes or per-view framebuffers, making performance scale with **pixel count** and **view count**.

So “performance” is usually about:

1) reducing hot-path allocations/recompute frequency, or  
2) reducing GPU work per frame, or  
3) reducing bytes/requests on the critical path.

:::

::::

---

## The 3 bottlenecks (what they feel like)

### 1) GPU-bound (rendering)

**What it feels like:**
- orbit/pan/zoom is choppy (“low FPS”),
- but the rest of the UI may still respond (buttons click, menus open).

**Common triggers in Cellucid:**
- many views/snapshots visible at once (grid view multiplies work),
- high visual quality (heavy shaders, large point sizes, volumetric smoke),
- GPU-heavy overlays (vector field / velocity),
- large browser window on a retina/high-DPI screen (more pixels).

**Fast confirmation tests (no tools required):**
1) Make the browser window smaller. If it immediately gets smoother, you were pixel/GPU-bound.
2) Clear snapshots (go back to fewer views). If it immediately gets smoother, you were view/GPU-bound.
3) Disable GPU-heavy modes (smoke mode, vector overlay). If it immediately gets smoother, you were GPU-bound.

**Typical fixes:**
- reduce the number of views,
- reduce visual quality knobs (see {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke` and {doc}`../i_vector_field_velocity/05_performance_and_quality`),
- keep the window smaller while exploring.

---

### 2) CPU-bound (computation)

**What it feels like:**
- the whole page “hitches” or freezes briefly after an action,
- sliders feel laggy, typing/search feels delayed,
- the fan spins up even if you’re not navigating in 3D.

**Common triggers in Cellucid:**
- repeated filtering recomputation while scrubbing sliders (especially with Live filtering),
- many enabled filters (stacking multiplies per-cell checks),
- heavy analysis operations on large groups,
- extremely large category accounting (very high category counts can make counting/legends expensive).

**Fast confirmation tests:**
1) If orbit/pan/zoom is smooth but filtering/analysis actions stutter, you’re likely CPU-bound.
2) Turn Live filtering off (continuous fields) and apply once with `FILTER`. If it becomes smooth, you were recompute-bound.

**Typical fixes:**
- avoid “recompute loops” (turn Live filtering off; make changes then apply once),
- reduce the number of enabled filters,
- do coarse gating upstream in Python for very large datasets (treat UI filtering as refinement).

See {doc}`../e_filtering/05_performance_considerations` for the filtering-specific cost model.

---

### 3) I/O-bound (loading from disk/network)

**What it feels like:**
- loading a dataset takes a long time,
- switching to a gene pauses while data loads,
- things are fast *after* loading, but slow at the moment data is requested.

**Common triggers in Cellucid:**
- loading large `.h5ad` directly in the browser (not truly lazy; can freeze/crash),
- remote server latency or throttled networks,
- first-time loads of gene expression chunks (cold cache),
- loading from slow disks / network drives / remote desktop mounts.

**Fast confirmation tests:**
1) If the slowdown happens only when “loading something new” (dataset/gene), suspect I/O.
2) If the same action is faster the second time, you’re seeing caching effects (cold vs warm).

**Typical fixes:**
- use **Server Mode** for large `.h5ad` / huge datasets (see {doc}`../b_data_loading/04_server_tutorial`),
- prefer exports that support lazy loading over “browser-only big file” workflows,
- benchmark with controlled “cold vs warm cache” runs (see {doc}`04_benchmarking_methodology_and_metrics`).

---

## The multipliers (what turns “fine” into “slow”)

Most costs in Cellucid are not “mysterious”—they scale with a small set of multipliers.

| Multiplier | Why it matters | Typical symptoms |
|---|---|---|
| `n_cells` | Many operations are per-cell (visibility masks, colors, selections) | filters/updates take longer; memory pressure grows |
| `n_views` (snapshots) | Each view adds GPU work, and sometimes per-view state | grid view stutters; context lost appears sooner |
| Window pixel count (`width × height × dpr²`) | Overlays and post-processing scale with pixels | smoother when window is smaller |
| `n_enabled_filters` | Filtering often scales like `O(n_cells × n_filters)` | slider scrubbing becomes the main lag source |
| Category count | Legends/counts and color mapping can degrade with huge category explosions | legends slow; UI becomes cluttered/unusable |
| “First-time” loads | cold cache → extra I/O | first gene switch slow; second is faster |

:::{tip}
If you’re doing performance work, treat **views** and **pixels** as first-class multipliers.

Reducing views (clear snapshots) and reducing pixels (smaller window) are the two fastest “no-risk” ways to confirm a GPU bottleneck.
:::

---

## A triage workflow that works (every time)

Use this workflow whenever something is slow.

### Step 0 — Name the action

Be specific. “It’s slow” is not a reproducible problem statement.

Good examples:
- “Dragging Min/Max slider for a QC field lags.”
- “Orbiting the 3D view is choppy in grid view.”
- “Switching genes takes 10–20 seconds the first time.”

### Step 1 — Identify the likely bottleneck (GPU vs CPU vs I/O)

Do one fast confirmation test:

- **GPU suspicion:** shrink the window; disable overlays; clear snapshots.
- **CPU suspicion:** turn Live filtering off; do one filter apply; avoid repeated actions.
- **I/O suspicion:** repeat the same action (cold vs warm); try a local server instead of remote.

### Step 2 — Change one knob that targets that bottleneck

Avoid changing 5 settings at once—you won’t know what helped.

Examples:
- GPU: reduce views, reduce overlay density, avoid smoke.
- CPU: apply filters once, reduce enabled filters, reduce analysis scope.
- I/O: use server mode, move data to local disk, reduce data requested.

### Step 3 — Measure (even roughly)

You don’t need perfect metrics. A stopwatch + “feels smoother” is enough to pick a direction.

If you are comparing versions/settings, use `04_benchmarking_methodology_and_metrics` for a repeatable approach.

### Step 4 — Lock in a “safe preset”

Once you find a fast, stable configuration:

- keep that as your “baseline preset” (especially on laptops),
- then increase quality/complexity gradually (one knob at a time).

---

## What’s normal vs suspicious

### Normal (expected)

- First-time gene loads are slower than the second time (cache warming).
- Filtering is slower when `n_cells` is huge, especially with Live filtering on.
- Grid view is slower than single view (multiple views draw in parallel).
- Smoke mode and vector overlays can be dramatically slower than points mode (they do more GPU work).

### Suspicious (worth troubleshooting)

- The app becomes permanently slow even after clearing snapshots and disabling overlays.
- You repeatedly hit “WebGL context lost” on moderate datasets (could be a driver/browser issue).
- The tab’s memory usage grows continuously during normal interaction (possible leak or unbounded caching).
- The slowdown happens on “small” datasets in the same way it happens on huge ones (could be an environment problem).

If any of these match, start with {doc}`../a_orientation/02_system_requirements` and then go to {doc}`07_troubleshooting_performance`.

---

## If you remember only five rules

1) **Views multiply work.** Keep snapshots lean while exploring.  
2) **Avoid recompute loops.** Turn Live filtering off and apply once on large data.  
3) **Pixels matter.** Smaller windows (or lower-DPI monitors) can dramatically improve GPU-bound workflows.  
4) **Use the right loading mode.** Big `.h5ad` in the browser is a trap; server mode is your friend.  
5) **Report performance with context.** Dataset size + hardware + steps + a number beats “it’s slow”.

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: performance-three-bottlenecks-diagram
Suggested filename: benchmarking_performance/01_three-bottlenecks-diagram.png
Where it appears: User Guide → Web App → Benchmarking and Performance → 01_performance_mental_model.md
Capture:
  - Asset type: diagram (SVG preferred) or annotated screenshot
  - Content: GPU-bound vs CPU-bound vs I/O-bound, with “what it feels like” examples
  - Goal: give non-technical users a quick decision tree
Crop:
  - Include only the diagram (no browser chrome)
Redact:
  - Use generic labels (no private dataset names)
Alt text:
  - Diagram explaining GPU-, CPU-, and I/O-bound performance bottlenecks in Cellucid.
Caption:
  - Cellucid performance triage: decide whether you are GPU-bound, CPU-bound, or I/O-bound, then change the one knob that targets that bottleneck.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder diagram for performance bottlenecks.
:width: 100%

Cellucid performance triage: decide whether you are GPU-bound, CPU-bound, or I/O-bound, then change the one knob that targets that bottleneck.
```

---

## Next steps

- {doc}`02_performance_considerations_what_gets_slow_and_why` (the concrete cost model + “which knob helps which problem”)
- {doc}`03_large_dataset_best_practices` (step-by-step workflow for huge datasets)
- {doc}`07_troubleshooting_performance` (symptom → diagnosis → fix)
- {doc}`09_reporting_performance_bugs` (copy/paste bug report template)
