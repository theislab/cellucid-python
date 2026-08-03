# Edge cases (performance)

**Audience:** everyone (support-oriented; useful before filing a bug)  
**Time:** 15–30 minutes  
**What you’ll get:**
- The specific configurations where Cellucid goes from smooth to unusable
- A confirm-then-fix recipe for each, in the same shape every time

This page catalogs **high-impact performance edge cases**—the “performance cliffs” where Cellucid can go from smooth to unusable with a small change.

It is written in a “support-playbook” style:

- **What you see**
- **Why it happens**
- **How to confirm**
- **Fix**
- **Prevention**

If you’re here because something is already slow, you can also jump straight to {doc}`07_troubleshooting_performance`.

---

## Edge case: the two-view row (the layout that actually costs)

### What you see
- A single view feels fine, but adding one snapshot makes navigation choppy.
- Adding a *fourth* view does not make it noticeably worse — and may improve it.

### Why it happens
This one surprises people, so it is worth stating precisely.

Cellucid renders at most four views: the live view plus up to three snapshots.
Two or three views are laid out as **one row**; four are laid out as a **2×2
grid**.

Point size scales with the height of the pane a point is drawn in. In a one-row
layout every pane is still full canvas height, so the points stay full size and
you pay for each pane in full — two views shade about twice the pixels of one.
In the 2×2 grid each pane is half height, so every point is drawn at half the
diameter and covers about a quarter of the area. Four panes at a quarter of the
area each add up to roughly the fill cost of a single view.

So the cost curve is not monotonic: **1 ≈ 4 < 2 < 3.**

Overlay buffers and post-processing passes are allocated per view and do not get
this benefit, so a heavy vector-field overlay still scales with the view count
in every layout.

### How to confirm
1) Go back to a single view. If FPS recovers, the layout was the cost.
2) If you were on two or three views, try adding one more to reach four. If FPS
   recovers there too, you were fill-rate bound and the smaller point size fixed
   it.

### Fix
- Tune filters/fields in the live view first.
- Choose the layout for the comparison you actually want, then check both
  directions if it is slow — one view and four views are the cheap ones.
- If overlays are on, reduce their density before blaming the layout.

### Prevention
- Treat “Keep view” as a deliberate comparison action, not a running history of
  everything you tried. There are only three snapshot slots.

Related docs:
- {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`
- {doc}`03_large_dataset_best_practices` (workflow patterns)

---

## Edge case: high-DPI / retina screens (pixels are expensive)

### What you see
- The app is smooth on an external monitor but choppy on a laptop screen.
- Performance improves dramatically when you shrink the browser window.

### Why it happens
GPU work scales with pixels. On retina screens, the browser can render at 2× (or more) device pixel ratio, which can mean ~4× pixels to shade.

This hurts:
- post-processing,
- overlays,
- smoke mode,
- and any full-screen passes.

### How to confirm
1) Shrink the window and retry.
2) If performance improves immediately, you were pixel-bound.

### Fix
- Use a smaller window while exploring.
- Reduce overlay quality (density/bloom) if using overlays.
- Avoid smoke mode at high settings on laptops.

### Prevention
- For workshops: recommend an external monitor or “don’t full-screen” guidance for older GPUs.

---

## Edge case: two settings changed themselves when you opened another dataset

### What you see
- The same machine, the same window, the same layout — and a dataset that draws
  visibly differently, or measures differently, from the one you had open before.
- Points look chunkier or finer than you expected, or edges look harder.
- A benchmark comparison between two datasets does not reproduce.

### Why it happens
Two rendering settings are derived from the cell count each time a dataset is
published, and datasets are switched **in place**:

- **Point size** is chosen so the total drawn area stays roughly constant, so
  the dot diameter falls with the square root of the cell count.
- **`Antialiasing (smooth point edges)`** is chosen automatically — on below
  5,000,000 cells, off at or above.

Neither is a bug and neither is a reset: they are the right opening answers for
two datasets that are not the same size. But nothing announces the point-size
change, so it reads as the app behaving differently for no reason.

### How to confirm
1) Read the line under the antialiasing checkbox. While the choice is still
   automatic it names the count, as in `Chosen automatically: 18,142,044 cells.`
   Once you have clicked the checkbox yourself the line goes quiet and the box
   holds your answer on every dataset from then on.
2) Compare `Point size (log):` against what you remember. The opening size is
   applied on every dataset publication — but *before* a saved session or a
   published preset is replayed, so a size you stored still wins.

### Fix
- Set both explicitly before comparing anything: they are live settings and
  antialiasing applies to the next frame.

### Prevention
- Record both in any benchmark or bug report — see
  {doc}`04_benchmarking_methodology_and_metrics`.

Related docs:
- {doc}`02_performance_considerations_what_gets_slow_and_why` (what each one costs, measured)

---

## Edge case: integrated GPU or low VRAM (context loss risk)

### What you see
- “WebGL context lost. Reload required…”
- The plot area becomes blank after a heavy action.
- The app repeatedly crashes on large datasets.

### Why it happens
Large datasets and GPU-heavy modes allocate big GPU buffers. If the browser/GPU runs out of VRAM, WebGL contexts can be lost (this is common behavior; the browser is protecting system stability).

Highlight overlays are a quiet contributor: a highlighted cell costs GPU memory
per view it is highlighted in, so a million-cell highlight group across a
four-view layout is a large standing allocation that never appears in the
dataset size.

### How to confirm
- It often happens right after:
  - enabling smoke mode with high grid density,
  - enabling overlays at high density,
  - loading a very large dataset,
  - selecting a very large highlight group,
  - or adding snapshots.

### Fix (safe → aggressive)
1) Reload the page.
2) Reduce GPU load:
   - points mode (avoid smoke),
   - fewer views,
   - disable overlays or lower density,
   - smaller window.
3) If it persists on moderate data: confirm hardware acceleration and WebGL2,
   then update the browser, operating system, and GPU driver when applicable.

### Prevention
- Keep a “laptop-safe preset” (points mode + fewer views + no heavy overlays).

Related docs:
- {doc}`../a_orientation/02_system_requirements` (context lost troubleshooting)

---

## Edge case: smoke mode quality cliffs (grid density + ray quality)

### What you see
- Smoke mode is extremely slow or immediately causes context loss.
- Smoke mode looks blank or “too faint”.

### Why it happens
Smoke mode is a ray-marched volumetric rendering:
- **Grid density** increases density-build time and GPU memory;
  the current exact choices run from 32³ through 128³ and start at 128³.
- **Visible-cell count** increases density-build validation and submission.
- **Render resolution** increases the pixels rendered each frame.
- **Ray quality** increases ray-march steps per pixel.

### How to confirm
- Lower render resolution first, then ray quality, when frame rate is poor.
- Lower grid density when entering Smoke or rebuilding after a filter is slow,
  or when GPU memory pressure is the problem. Grid changes are applied after
  the approximately 300 ms slider debounce and then rebuild the volume.

### Fix
- Lower render resolution and ray quality in **Volumetric smoke:**.
- Lower grid density for build latency or memory pressure.
- Lower noise detail if memory pressure remains.
- Return to points mode when doing scientific work; use smoke selectively.

### Prevention
- Treat smoke mode as “presentation”: start low and scale up slowly.

Related docs:
- {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke`

---

## Edge case: vector field overlay across views (pixels × views × post-processing)

### What you see
- Overlay is fine in one view but unusable once snapshots are added.
- Overlay becomes slow on large screens or with bloom enabled.

### Why it happens
The overlay often allocates per-view buffers and runs post-processing passes.
Costs scale with:
- pixel count,
- number of views,
- particle density,
- and bloom/blur passes.

### How to confirm
1) Disable bloom (or set bloom strength to 0.00) and retry.
2) Reduce particle density and retry.
3) Return to a single view and retry. Unlike point drawing, overlay cost really
   does scale with the number of views in every layout — the smaller points of a
   2×2 grid do not shrink a particle system or its post-processing passes.

### Fix
- Use a conservative preset (low density, short trails, no bloom) during exploration.
- Enable overlays only when you need them.

### Prevention
- Separate a “scientific preset” (fast, conservative) from a “cinematic preset”.

Related docs:
- {doc}`../i_vector_field_velocity/05_performance_and_quality`

---

## Edge case: category explosion (the UI becomes unusable)

### What you see
- A categorical field has thousands–tens of thousands of categories.
- Legends become enormous; interactions and counts lag; colors are meaningless.

### Why it happens
Even if rendering points is fine, the UI cost of:
- listing categories,
- computing counts,
- assigning colors,
- and presenting a usable legend

can become a CPU and UX cliff.

### How to confirm
- Switch to a simpler categorical field; if the UI becomes responsive, the category explosion is the issue.

### Fix
- Use a different field for interactive work.
- Collapse or group categories upstream (data preparation).
- Treat unique IDs (barcodes, UUIDs) as *identifiers*, not visualization categories.

### Prevention
- Export “human-facing” categorical fields with manageable category counts.

Related docs:
- {doc}`../d_fields_coloring_legends/index`

---

## Edge case: repeated slider scrubbing (death by recompute)

### What you see
- Dragging a Min/Max slider causes stutters or a delayed response.
- CPU usage spikes during slider movement.

### Why it happens
With Live filtering enabled, each slider movement triggers a full re-derivation
of visibility across every cell, testing every enabled filter.

The important detail: that cost depends on the dataset size, **not** on how much
the filter changed. Nudging a slider so that eight more cells are hidden costs
the same pass as hiding half the dataset. Only the resulting upload to the GPU
is proportional to what moved, and the upload is not the part you feel.

This is especially punishing when:
- `n_cells` is large, and/or
- many filters are enabled.

### How to confirm
- Turn Live filtering off and apply once; if it becomes smooth, you found it.

### Fix
- Live filtering off → stage slider changes → click `FILTER` once.
- Reduce enabled filters (disable or remove no-op filters).

### Prevention
- Teach the “apply once” pattern for large datasets.

Related docs:
- {doc}`../e_filtering/05_performance_considerations`

---

## Edge case: huge highlight groups and session bloat

### What you see
- Selecting/highlighting huge groups becomes slow.
- Saving/restoring sessions takes a long time.
- Session files become very large.

### Why it happens
Huge per-cell sets (highlights/selections) can:
- require large arrays/bitsets,
- require merge/priority resolution for rendering,
- and bloat session artifacts if stored.

### How to confirm
- Clear highlights and retry; if responsiveness returns, highlight scale is the culprit.

### Fix
- Keep highlight groups smaller when possible.
- Use highlights for targeted comparisons, not as a full export of membership for millions of cells.
- If you need full membership tables, export results and work offline.

### Prevention
- Define conventions for highlights in team workflows (e.g., “no >500k-cell highlight groups in shared sessions”).

Related docs:
- {doc}`../f_highlighting_selection/05_edge_cases_highlighting`
- {doc}`../l_sessions_sharing/index`

---

## Edge case: large `.h5ad` loaded in-browser (memory cliff)

### What you see
- The tab freezes during load.
- The browser crashes or becomes unresponsive.

### Why it happens
Browser `.h5ad` loading is not truly lazy; large files can exceed browser memory budgets.

### How to confirm
- Try the same file via Server Mode; if it works there, the in-browser load path was the cliff.

### Fix
- Use Server Mode for large `.h5ad` / `.zarr`.

### Prevention
- Treat in-browser `.h5ad` load as a “small quick preview” workflow only.

Related docs:
- {doc}`../b_data_loading/03_browser_file_picker_tutorial`
- {doc}`../b_data_loading/04_server_tutorial`

---

## Edge case: thermal throttling (laptops get slower over time)

### What you see
- The first minute is smooth, then it gets slower.
- Fans spin up, the laptop gets hot, and performance degrades.

### Why it happens
GPU and CPU clocks can drop under sustained load to protect the device (thermal throttling).

### How to confirm
- Wait 2–5 minutes; if performance steadily degrades without changing settings, thermals are likely involved.

### Fix
- Lower GPU load (fewer views, smaller window, disable overlays/smoke).
- Plug in power (some laptops downclock on battery).
- Allow cool-down time.

### Prevention
- Use conservative presets for long sessions.

---

## Edge case: extension / corporate policy interference

### What you see
- Performance is unexpectedly bad on “normal” datasets.
- Some interactions lag only in certain browsers or managed environments.

### Why it happens
Extensions and policies can:
- disable hardware acceleration,
- inject scripts,
- block resources,
- or interfere with file APIs.

### How to confirm
- Use a clean/private profile with extensions disabled.
- Confirm WebGL2 from the console with
  `document.createElement("canvas").getContext("webgl2") !== null`.

### Fix
- Disable interfering extensions for the site.
- In corporate environments, ask IT to allow hardware acceleration and required origins.

### Prevention
- For enterprise deployments: document the required WebGL2, hardware
  acceleration, download, storage, and origin policies.

Related docs:
- {doc}`../a_orientation/02_system_requirements`

---

## Next steps

- {doc}`07_troubleshooting_performance` (symptom → diagnosis → fix)
- {doc}`09_reporting_performance_bugs` (how to report a performance cliff with the right context)
