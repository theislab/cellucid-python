# Benchmarking methodology and metrics

**Audience:** computational users + developers (still usable for “I just need numbers”)  
**Time:** 20–45 minutes  
**What you’ll learn:**
- What to measure (and what each metric actually means)
- How to run a benchmark that is reproducible and comparable across runs
- How to separate cold-cache vs warm-cache behavior
- A simple benchmark worksheet you can reuse for regressions and bug reports

**Prerequisites:**
- A dataset you can load in Cellucid
- A current stable Chrome, Edge, Firefox, or Safari release; use the same
  browser version for before/after comparisons

---

## Why benchmark at all?

Benchmarks turn:

- “it feels slower”

into:

- “load time increased from ~6s to ~11s on the same dataset/hardware”

That is the difference between a performance issue that can be fixed and one that becomes a long back-and-forth.

---

## The core rule: benchmark a scenario, not “the app”

Cellucid is interactive. “Performance” depends on what you’re doing.

Always define a benchmark scenario in one sentence, for example:
- “Cold-load dataset X into points mode and reach a stable first render.”
- “Switch from a categorical field to gene expression (first gene load).”
- “Apply a continuous filter range on a 1M-cell dataset with Live filtering off.”
- “Enable vector field overlay at 15K particles in a 2×2 grid.”

If you can’t name the scenario, you can’t compare results.

---

## What to measure (recommended metrics)

You do not need all of these every time. Pick the ones that match your scenario.

### 1) Load metrics

**Time to first render (TTFR)**  
Time from “user initiates load” to “points are visible and the UI is responsive”.

**Time to interactive (TTI)**  
Time until common interactions (orbit/zoom, color change) respond without long stalls.

Why both matter:
- TTFR can be fast but the app can still be “doing work” and feel unresponsive for a bit.

---

### 2) Interaction metrics (responsiveness)

**Frame rate (FPS)** (GPU-bound indicator)  
Lower FPS usually means the GPU is overloaded (or pixels/views are too high).

**Input latency / hitching** (CPU-bound indicator)  
Long pauses after actions (filter apply, selection updates) usually mean CPU-bound work.

Practical note:
- Humans notice hitching more than they notice a small FPS change.

---

### 3) Memory metrics (stability)

**Tab memory footprint**  
High memory can cause:
- crashes,
- “WebGL context lost” (GPU memory pressure),
- slowdowns due to garbage collection (CPU stalls).

You typically care about:
- whether memory stabilizes after a while (good),
- or grows unbounded during normal use (suspicious).

**GPU memory, roughly.** Before the renderer has published a measured figure,
Cellucid's own reports estimate GPU memory at **28 bytes per point** — about
27 MiB at one million points and 267 MiB at ten million. Read it as an order of
magnitude for whether a dataset will fit, not as a measurement: the report marks
it `GPU Memory (estimate)` while it is one, and switches to `GPU Memory` once
the renderer reports the real value. Highlights, snapshot views and overlays all
add on top of it.

---

### 4) Network / I/O metrics (remote workflows)

For remote servers or lazy-loading gene expression:
- request count,
- total bytes transferred,
- and latency distributions can dominate perceived performance.

---

## A reproducible benchmark workflow (step-by-step)

This is the workflow that produces comparable results without fancy tooling.

### Step 1 — Record your benchmark context (always)

At minimum record:
- dataset identifier and size (cells, genes, fields; approximate is fine),
- loading mode (browser-export, server mode, remote server),
- browser + version,
- machine + GPU (even “M1 MacBook Air” is better than nothing),
- window size + whether retina/high-DPI (pixels matter),
- the exact Cellucid settings that affect the scenario (views, render mode, overlays),
- `Point size (log):`, whether `Antialiasing (smooth point edges)` is ticked, and
  whether `Level-of-Detail (LOD)` is on and forced to a level.

:::{important}
Record those last three even though you did not set them. The point size and the
antialiasing state are both chosen from the dataset's cell count when it opens,
so **two runs on two datasets can differ in them without anyone touching a
control** — and they are the two rendering settings measured to move frame time
in {doc}`02_performance_considerations_what_gets_slow_and_why`. A comparison that
does not pin them is not a comparison of the thing you changed.
:::

If you are reporting a performance bug, use the copy/paste template in {doc}`09_reporting_performance_bugs`.

---

### Step 2 — Control the big confounders

Before you measure:

1) Close other heavy tabs/apps (especially video calls and other WebGL content).
2) Keep the window size fixed for all runs (don’t benchmark full-screen once and small-window the next time).
3) Keep the **layout** fixed, not just the view count. One view, a two-view row
   and a 2×2 grid are three different workloads because point size follows pane
   height; they are not points on one curve.
4) Keep the filter state fixed. Hidden cells are rejected before rasterisation,
   so a run on a filtered view and a run on the full dataset measure different
   scenes.
5) Avoid changing multiple quality knobs between runs — and set `Point size
   (log):` and `Antialiasing (smooth point edges)` explicitly rather than
   inheriting whatever each dataset opened with.

Thermals matter:
- On laptops, performance can degrade after a few minutes due to throttling.
- If you’re doing a long benchmark session, allow cooldown breaks or use a fan/plugged-in power.

---

### Step 3 — Decide: cold-cache vs warm-cache

Performance can differ dramatically between:

- **cold cache**: first load, nothing cached
- **warm cache**: repeated actions where data is already in memory/cache

You usually want to measure both.

How to approximate a cold run:
- reload the page, and (optionally) disable cache in devtools for the duration of the run.

How to approximate a warm run:
- repeat the same scenario without reloading (e.g., switch gene A → gene B → gene A again).

Always label results as cold vs warm.

---

### Step 4 — Run multiple trials and summarize

Do not trust a single run.

Recommended:
- 5 runs for “quick check”
- 10 runs for “regression verification”

Summarize with:
- median (robust to outliers),
- and a rough range (min–max) or interquartile range if you want to be more rigorous.

---

## How to measure without special tools (practical)

You can get useful numbers with only your browser.

### Option A — Stopwatch (works for TTFR/TTI)

This is good enough for many “is it slower?” questions.

Define start/end points:
- Start: click “Load” / select dataset / press enter.
- End: points visible + you can orbit smoothly for 2–3 seconds.

Write down the time. Repeat 5×.

---

### Option B — Browser performance tooling

Chrome and Edge expose a tab-level task manager (`Shift+Esc`), a Performance
Monitor, and the Network panel:

**Chrome/Edge Task Manager (tab-level CPU/GPU/memory)**
Open with `Shift+Esc`.

Use it to confirm:
- CPU pegged during filtering/analysis (CPU-bound),
- GPU memory pressure (if the GPU memory column is enabled),
- whether the tab memory grows over time.

**DevTools → Performance Monitor (FPS + CPU + memory signals)**  
Open DevTools → Command menu → “Show Performance Monitor” (or More tools → Performance monitor).

Use it to watch:
- FPS while navigating (GPU-bound),
- JS heap growth (memory),
- CPU usage spikes during interactions.

**DevTools → Network (remote / lazy-loading)**  
Use the Network tab to see:
- request waterfalls,
- bytes transferred,
- long-latency outliers.

Firefox and Safari expose the same core evidence through their developer
tools: record a performance/timeline trace for CPU and frame activity, and use
the Network panel for request timing and transferred bytes. Record the browser
name and version with every result; never compare numbers from different
engines as if they were the same benchmark environment.

---

### Option C — Firefox equivalents

Firefox has similar tools:
- `about:performance` for tab impact,
- Performance profiling tools,
- Network monitor for request waterfalls.

The exact UI differs, but the logic is the same: separate GPU-bound vs CPU-bound vs I/O-bound.

---

## Benchmark worksheet template (copy/paste)

Use this table format in a GitHub issue, lab notebook, or internal doc.

### Context

- Dataset: `<name or id>`
- Size: `<n_cells> cells`, `<n_genes> genes`, `<n_fields> fields`
- Layout: `<1 view / 2-view row / 3-view row / 2x2 grid>`
- Loading mode: `<export folder / server mode / remote>`
- Browser: `<Chrome 123>`
- Machine/GPU: `<e.g., M2 Pro / RTX 3070 / Intel Iris>`
- Window: `<width × height>`, `<devicePixelRatio>`
- Render mode: `<Points/Smoke>`
- Point size: `<slider position>` → `<rendered size>`
- Antialiasing: `<on / off>`
- LOD: `<off / Auto / forced level N>`
- Overlays: `<none / vector field overlay settings>`
- Filters: `<none / N enabled>`, visible cells: `<count>`

### Results

| Scenario | Cold (median) | Cold (range) | Warm (median) | Warm (range) | Notes |
|---|---:|---:|---:|---:|---|
| TTFR (load → first render) | | | | | |
| Switch gene (first time) | | | | | |
| Apply filter (FILTER once) | | | | | |
| Navigate in the fixed layout (FPS) | | | | | |

---

## Interpretation guide (what the numbers usually mean)

This section is intentionally heuristic—use it to choose next steps.

### TTFR got worse (but FPS is fine)

Often I/O-bound:
- more bytes loaded,
- slower server/disk,
- more metadata to parse,
- or a regression in loading pipeline.

Next step:
- capture a Network waterfall (if remote),
- or compare cold vs warm to see caching effects.

### FPS got worse (especially when moving camera)

Often GPU-bound:
- more pixels,
- a layout change (a two- or three-view row is the expensive shape),
- fewer cells filtered out than in the run you are comparing against,
- heavier shaders/overlays,
- or a driver/browser regression.

Next step:
- reduce pixels first; see if it recovers immediately. Then check that the
  layout and the filter state match the baseline run.

### Actions “hitch” (stalls after clicking/filtering)

Often CPU-bound:
- visibility recomputation loops,
- large group operations,
- garbage collection due to memory churn.

Next step:
- avoid repeated slider scrubbing; apply once,
- reduce enabled filters; benchmark again.

---

## Next steps

- {doc}`05_benchmark_tools` (in-app synthetic benchmark and situation report)
- {doc}`07_troubleshooting_performance` (turn symptoms into diagnoses)
- {doc}`09_reporting_performance_bugs` (what to include when filing a performance issue)
