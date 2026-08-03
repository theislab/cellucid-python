# Benchmark tools

**Audience:** developers + performance-minded users  
**Time:** 10–25 minutes  
**What you’ll learn:**
- How to run the built-in **synthetic rendering benchmark**
- How to interpret the metrics (FPS, frame time, GPU memory, LOD)
- How to use “Copy Situation Report” for bug reports and regressions
- What this benchmark does *not* measure (and how to benchmark those separately)

**Prerequisites:**
- Cellucid open in a current stable desktop release of Chrome, Edge, Firefox,
  or Safari with WebGL2 enabled

---

## Where the benchmark lives (UI location)

Open:

- Left sidebar → **Performance Benchmark** accordion

It includes:
- six quick presets: 100K, 500K, 1M, 5M, 10M, 20M points,
- a custom point count input accepting any whole number from 1 to 50,000,000,
- eight synthetic **data patterns**,
- a **Load Synthetic Data** button,
- a metric grid that stays hidden until the first run finishes,
- a **Copy Situation Report** button,
- and an **Analyze Performance** tool that suggests what to change.

```{figure} ../../../_static/screenshots/benchmarking_performance/benchmark-panel-controls.png
:alt: The Performance Benchmark panel before any run, showing the six point-count presets from 100K to 20M, the custom point count field, the data pattern dropdown, and the Load Synthetic Data, Copy Situation Report and Analyze Performance buttons.
:width: 516px

The panel before its first run. Everything above **Load Synthetic Data** is
input; everything below it appears only once a run has produced numbers.
```

:::{important}
The synthetic benchmark **loads synthetic data into the viewer**.

Use it when:
- you want to test rendering limits without needing a real dataset, or
- you want to compare hardware/browser builds.

Do **not** use it to benchmark:
- real dataset loading,
- gene expression fetching,
- analysis runtime,
- or figure export.

For those, use {doc}`04_benchmarking_methodology_and_metrics`.
:::

---

## What this benchmark measures (and why it’s useful)

This benchmark is primarily a **rendering stress test**:

- how many points your GPU can render smoothly,
- how well LOD/frustum culling helps,
- and how GPU memory scales with point count.

It is especially useful for answering questions like:
- “Can this laptop handle a 1M-cell dataset interactively?”
- “Did a renderer change reduce FPS on the same machine?”
- “Is the bottleneck likely GPU vs CPU?”

---

## Quick start (3 minutes)

1) Open **Performance Benchmark** in the left sidebar.
2) Click a preset like **500K** or **1M**.
3) Click **Load Synthetic Data**.
4) Orbit/pan/zoom for ~10 seconds and watch:
   - **FPS**
   - **Frame Time**
   - **GPU Memory**
5) If it feels slow, click **Analyze Performance** and read the “What to do” list.

If you want a copy/paste summary for a bug report, click **Copy Situation Report**.

---

## Recommended workflow (reproducible)

If you’re comparing two builds or changes, make the benchmark comparable:

1) Fix window size (pixels matter).
2) Fix the layout (single view only; clear snapshots). A two-view row and a 2×2
   grid are different workloads, not “more of the same one”.
3) Start in points mode (smoke mode is a different workload).
4) Use the same preset count and the same pattern.
5) Let FPS stabilize for ~5–10 seconds, then record the metric.
6) Wait for **Timing Details** to appear before trusting p95 — it needs at least
   30 sampled frames.

If you want statistics, repeat the run 3–5 times and report the median.

Because generation is seeded, the same pattern and count give you the same
points every time, so a difference between two runs is a difference in the
renderer or the machine — never in the data.

---

## Data patterns (how to choose)

The pattern dropdown controls the spatial distribution of synthetic points. The
eight options, in the order they appear:

| Pattern | What it produces | Use it to test |
|---|---|---|
| **Model Surface** *(default)* | Points sampled over a loaded 3D model surface | A dense, structured shell — the pattern the acceptance samples used |
| **Atlas-like** | Many groups spread across the space | “Many groups across space” intuition |
| **Batch Effects** | Offset duplicate structure | Overlapping-population intuition |
| **Gaussian Clusters** | Compact blobs | The closest analogue to a clustered embedding |
| **Uniform Random** | Points everywhere, no structure | The harshest case for LOD and frustum culling |
| **Octopus** | Radiating arms | Elongated, trajectory-like structure |
| **3D Spirals** | Interleaved spiral arms | Fine structure that LOD can destroy visibly |
| **Flat UMAP** | A dense 2D manifold | The common 2D-embedding shape |

If you are testing *rendering limits*, Uniform Random is often the harshest,
because nothing can be culled. If you are testing *realistic single-cell
scenes*, Gaussian Clusters or Atlas-like are more representative.

Generation is seeded, so the same pattern and count produce the same points on
every run, and it happens on a worker thread — the interface stays responsive
while a large set is built.

---

## What the panel reports (every field, and its unit)

After a run, the grid shows six headline tiles, then a timing breakdown, then
the generation time.

| Field | Unit / format | What it means |
|---|---|---|
| **Points** | count | How many synthetic points were published |
| **FPS** | whole frames per second | Derived from recent frame times |
| **Frame Time** | `X.XX ms`, sometimes with a `(render Y.YY ms)` suffix | Wall time for a frame |
| **GPU Memory** | `N.N MB`, or `~N.N MB (estimate)` | Reported where the browser exposes it; otherwise estimated from the point count |
| **LOD Level** | `Full`, or `Level N` | Which level of detail the renderer settled on |
| **Visible Points** | count | How many points survived LOD and culling |
| **Timing Details → min / p95 / max** | `X.XX ms` each | Frame-time distribution. Appears only once at least 30 frames have been sampled |
| **Generated in … ms** | whole ms | How long the synthetic set took to build |

:::{important}
The six headline tiles are **CPU-side wall-clock measurements** taken between
frames. They are not GPU timer readings, so they include anything else the tab
was doing.

Real GPU timing exists, but only inside **Analyze Performance**, which uses the
browser’s disjoint-timer-query extension where the driver provides it and falls
back to CPU estimates where it does not. If you need to attribute time to the
GPU specifically, that is the tool — not the headline FPS.
:::

### Reading them

**FPS and frame time**

- **60 FPS** is the “smooth” target on most displays.
- 16.7 ms/frame ≈ 60 FPS; 33 ms/frame ≈ 30 FPS.
- If FPS drops as you increase point count: you are hitting a GPU limit.
- If FPS is unstable with long spikes (high p95/max relative to min): you are
  seeing jank, which is usually a CPU or allocation problem rather than a
  drawing one.

**GPU memory**

- It is the main predictor of “WebGL context lost” on large datasets.
- If it climbs close to your GPU’s limit, expect context loss.
- If memory is fine but FPS is low, you are pixel- or compute-bound, not
  memory-bound.
- If the value is prefixed with `~`, the browser did not expose a real figure
  and it is an estimate derived from the point count. Do not compare an
  estimate against a measured value.

**LOD level and visible points**

- LOD exists to reduce the number of rendered points when zoomed out.
- If LOD reduces **Visible Points** and FPS improves, the renderer is adapting
  successfully.
- If **Visible Points** stays high and FPS stays low, try LOD in Visualization —
  and lower `Force LOD level:` too when you need to go past `Auto`'s floor of
  `min(2,000,000, cells ÷ 8)` points,
  a smaller window, a smaller point size, or turning antialiasing off. Lowering
  `Shader quality:` is *not* on that list — it did not change frame time on the
  hardware this was measured on, and the reason is in
  {doc}`02_performance_considerations_what_gets_slow_and_why`.

Related docs:
- {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke` (LOD, culling, shader quality)

---

## “Copy Situation Report” (use this for bug reports)

The **Copy Situation Report** button is meant to produce a compact block of text you can paste into:
- a GitHub issue,
- a performance regression note,
- or a lab internal thread.

Best practice:
- include the report *and* your scenario description (“what I was doing when it was slow”).

If you’re filing an issue, follow {doc}`09_reporting_performance_bugs`.

---

## “Analyze Performance” (why is it slow?)

The **Analyze Performance** button runs a short measurement and produces:
- a **Problems Found** list naming one primary bottleneck and any contributing
  ones,
- a prioritized **What to do** checklist,
- and a collapsed **Show detailed stats** block with the per-phase overheads.

The bottleneck names it can print are `GPU`, `CPU`, `Fragment/Fill Rate`,
`Vertex Processing`, `LOD Selection`, `Frustum Culling`, `Frame Stuttering`,
`JS/CPU Pressure`, `Garbage Collection` and `Main Thread Blocked`. Use it as a
triage assistant:

- `Fragment/Fill Rate` — you are covering too many pixels. Reduce the point size
  first, then the window size, then the view count.
- `Vertex Processing` — you are drawing too many points. Enable LOD, or filter.
- `LOD Selection` / `Frustum Culling` — the *selection* work, not the drawing,
  is showing up. Avoid forcing a high LOD level.
- `CPU`, `JS/CPU Pressure`, `Main Thread Blocked` — the frame is not
  GPU-limited at all; look at filtering, analysis and category accounting.

:::{note}
Earlier builds printed a `Shader Complexity` bottleneck and a
`Shader overhead:` figure, and suggested switching to `Ultra-light`. **They are
gone**, and the detailed stats now report `Point size response:` in that row
instead — the multiplier by which frame time grows as the point size is swept,
computed from a sweep the analysis already performs.

The shader figure was removed rather than corrected because the instrument
cannot support it. It timed the three qualities in a fixed order, one run each,
for about sixteen frames apiece, and compared the spread against fixed
thresholds. The effect itself is not there — the three qualities are within
noise at every point size, because nearly every rasterised fragment dies at the
depth test before the fragment shader runs. Adding counterbalancing and a
measured noise floor was tried and still produced a false positive on a busy
machine, so no threshold on that instrument is trustworthy.

**Point size** is the rendering setting that responds, and **antialiasing** is
the largest single lever left on the plain path. Both are measured in
{doc}`02_performance_considerations_what_gets_slow_and_why`. Running the
analysis no longer changes your shader quality as a side effect.
:::

---

## Safety notes and footguns

- Synthetic benchmark data is not your real dataset. Reload the page to restore your original data.
- Don’t start with 20M points on a laptop. Ramp up: 100K → 500K → 1M → 5M.
- If you hit “WebGL context lost”, reduce point count and/or reduce GPU load (smaller window, no snapshots, no smoke).
- The custom point count is capped at 50,000,000. That is a deliberate refusal,
  not a hardware limit: the position and colour arrays for that many points
  already run to hundreds of megabytes before anything reaches the GPU. A value
  outside 1–50,000,000 is rejected with a **Point count out of range** notice
  rather than attempted.
- The controls in this panel start out disabled and become usable once the app
  has finished wiring them. If a preset looks greyed out for a moment during a
  slow start, that is the intended behaviour — see
  {doc}`../q_troubleshooting_index/index`.

---

## Next steps

- {doc}`04_benchmarking_methodology_and_metrics` (benchmark scenarios beyond synthetic rendering)
- {doc}`06_edge_cases_performance` (performance cliffs and limits)
- {doc}`09_reporting_performance_bugs` (how to file a performance issue with the right context)
