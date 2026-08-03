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

The practical trick: the *CPU* side of an action usually costs a pass over every
cell even when the action changes almost nothing, while the *GPU* side only
moves what actually changed. So the thing to avoid is a repeated recompute loop
(slider scrubbing with Live filtering on), not the size of the resulting change.

:::

:::{tab-item} Developer / Deep

At a high level:

- The viewer maintains large arrays (positions, colors, visibility, selection, sometimes per-view buffers).
- A settled frame uploads nothing. Buffer and texture writes are gated behind
  dirty flags, and the frame that follows a camera move reuses its published
  index buffer unless the set of admitted spatial nodes actually changed. The
  steady-state cost is uniforms plus draw calls.
- A visibility change uploads only the rows of the alpha texture whose bytes
  moved, so the transfer is proportional to how scattered the change is rather
  than to the dataset size. Deciding *which* cells changed is still a pass over
  every cell on the CPU.
- Some modules add full-screen passes or per-view framebuffers, making performance scale with **pixel count**.

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
- a two- or three-view layout (see the note below — it is the *shape* of the
  layout that matters, not just the count),
- large point sizes, or volumetric smoke (its shaders really are heavy; the
  points-mode ones are not — see
  {doc}`02_performance_considerations_what_gets_slow_and_why`),
- GPU-heavy overlays (vector field / velocity),
- large browser window on a retina/high-DPI screen (more pixels).

**Fast confirmation tests (no tools required):**
1) Make the browser window smaller. If it immediately gets smoother, you were pixel/GPU-bound.
2) Clear snapshots (go back to a single view). If it immediately gets smoother, you were view/GPU-bound.
3) Disable GPU-heavy modes (smoke mode, vector overlay). If it immediately gets smoother, you were GPU-bound.

**Typical fixes:**
- keep the window smaller while exploring,
- reduce `Point size (log):` once dots overlap on screen, and check whether
  `Antialiasing (smooth point edges)` is on — those are the two rendering
  settings measured to move frame time
  ({doc}`02_performance_considerations_what_gets_slow_and_why`),
- go back to one view, or go up to four (again, see the note),
- reduce visual quality knobs (see {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke` and {doc}`../i_vector_field_velocity/05_performance_and_quality`).

:::{note}
**More views is not automatically slower.** Cellucid allows at most four views —
the live view plus three snapshots — and it lays them out as one row for two or
three views, and as a 2×2 grid for four.

Point size scales with the height of the pane a point is drawn in, so in the 2×2
grid every point is drawn at half the diameter, covering roughly a quarter of
the pixels. Four panes × a quarter of the area each means the four-view grid
shades about the same total number of pixels as one full-size view. Two and
three views keep the full pane height, so they really do cost about twice and
three times a single view.

The practical consequence is counter-intuitive but real: on a fill-rate-limited
scene, **four views can be no more expensive than one, while two or three are
the expensive layouts.** This only holds while size attenuation is on (it is by
default); with attenuation at zero, point size no longer follows the pane and
every extra view is straightforwardly more work.
:::

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
| `n_cells` | Many operations are per-cell (visibility recomputation, colors, selections) | filter applies take longer; memory pressure grows |
| Point size | A dot's area grows with the square of its size, so point size sets how many pixels the renderer has to cover | frame time rises steeply once points start overlapping on screen |
| Layout shape (1, 2, 3 or 4 views) | Each view submits its own draw, but point size follows pane height — so the *row* layouts (2, 3) cost most and the 2×2 grid costs about as much as one view | a two-view comparison stutters where a four-view grid does not |
| Window pixel count (`width × height × dpr²`) | Overlays and post-processing scale with pixels | smoother when window is smaller |
| `n_enabled_filters` | Each filter change re-derives visibility for every cell, testing each enabled filter | slider scrubbing becomes the main lag source |
| Category count | Legends/counts and color mapping can degrade with huge category explosions | legends slow; UI becomes cluttered/unusable |
| “First-time” loads | cold cache → extra I/O | first gene switch slow; second is faster |

:::{important}
**Two of these were already set for you, from the cell count, when the dataset
opened.** The point size, and whether `Antialiasing (smooth point edges)` is on
— on below 5,000,000 cells, off at or above. Both are live settings you can
change at any time, and both are the settings measured in
{doc}`02_performance_considerations_what_gets_slow_and_why`. Check what they are
before assuming they sit at some fixed default.
:::

:::{tip}
If you’re doing performance work, treat **pixels** as the first-class multiplier
and **layout shape** as the second.

Reducing pixels (smaller window) is the fastest “no-risk” way to confirm a GPU
bottleneck. Changing the view count is a useful second test, but read it
carefully: going from two views to one should help, and going from three to
four may also help.
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

If you are comparing versions/settings, use {doc}`04_benchmarking_methodology_and_metrics` for a repeatable approach.

### Step 4 — Lock in a “safe preset”

Once you find a fast, stable configuration:

- keep that as your “baseline preset” (especially on laptops),
- then increase quality/complexity gradually (one knob at a time).

---

## What’s normal vs suspicious

### Normal (expected)

- First-time gene loads are slower than the second time (cache warming).
- Filtering is slower when `n_cells` is huge, especially with Live filtering on:
  the cost of *deciding* what is visible is per-cell and does not shrink when
  the filter hides most of the data.
- A two- or three-view row is slower than a single view; a 2×2 grid usually is
  not (see the note above).
- Smoke mode and vector overlays can be dramatically slower than points mode (they do more GPU work).
- Navigating a heavily filtered view is *faster* than navigating the same
  dataset unfiltered. Hidden cells are rejected before rasterisation, so they
  cost nothing to draw.

### Suspicious (worth troubleshooting)

- The app becomes permanently slow even after clearing snapshots and disabling overlays.
- You repeatedly hit “WebGL context lost” on moderate datasets (could be a driver/browser issue).
- The tab’s memory usage grows continuously during normal interaction (possible leak or unbounded caching).
- The slowdown happens on “small” datasets in the same way it happens on huge ones (could be an environment problem).

If any of these match, start with {doc}`../a_orientation/02_system_requirements` and then go to {doc}`07_troubleshooting_performance`.

---

## If you remember only five rules

1) **Pixels matter most.** Smaller windows (or lower-DPI monitors) can dramatically improve GPU-bound workflows — and so does a smaller point size, for the same reason.  
2) **Avoid recompute loops.** Turn Live filtering off and apply once on large data.  
3) **Measure the layout, don’t assume it.** One view and a 2×2 grid cost about the same; two- and three-view rows are the expensive ones.  
4) **Use the right loading mode.** Big `.h5ad` in the browser is a trap; server mode is your friend.  
5) **Report performance with context.** Dataset size + hardware + steps + a number beats “it’s slow”.

---

## Next steps

- {doc}`02_performance_considerations_what_gets_slow_and_why` (the concrete cost model + “which knob helps which problem”)
- {doc}`03_large_dataset_best_practices` (step-by-step workflow for huge datasets)
- {doc}`07_troubleshooting_performance` (symptom → diagnosis → fix)
- {doc}`09_reporting_performance_bugs` (copy/paste bug report template)
