# Edge cases (performance)

This page catalogs **high-impact performance edge cases**—the “performance cliffs” where Cellucid can go from smooth to unusable with a small change.

It is written in a “support-playbook” style:

- **What you see**
- **Why it happens**
- **How to confirm**
- **Fix**
- **Prevention**

If you’re here because something is already slow, you can also jump straight to {doc}`07_troubleshooting_performance`.

---

## Edge case: too many views (snapshots multiply everything)

### What you see
- Single view feels fine, but grid view becomes choppy or unresponsive.
- Adding “just one more snapshot” suddenly tanks FPS.

### Why it happens
Each view is effectively another render workload.

Many costs scale with:
- **`n_views`** (more viewports to render),
- and **pixels × views** for overlay buffers (vector trails, post-processing).

### How to confirm
1) Clear snapshots (return to a single view).
2) If FPS recovers immediately, you found the multiplier.

### Fix
- Tune filters/fields in the live view first.
- Keep snapshots to a small number (e.g., 2–4) for comparisons.

### Prevention
- Treat “Keep view” as a deliberate comparison action, not a running history of everything you tried.

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

## Edge case: integrated GPU or low VRAM (context loss risk)

### What you see
- “WebGL context lost. Reload required…”
- The plot area becomes blank after a heavy action.
- The app repeatedly crashes on large datasets.

### Why it happens
Large datasets and GPU-heavy modes allocate big GPU buffers. If the browser/GPU runs out of VRAM, WebGL contexts can be lost (this is common behavior; the browser is protecting system stability).

### How to confirm
- It often happens right after:
  - enabling smoke mode with high grid density,
  - enabling overlays at high density,
  - loading a very large dataset,
  - or opening many views.

### Fix (safe → aggressive)
1) Reload the page.
2) Reduce GPU load:
   - points mode (avoid smoke),
   - fewer views,
   - disable overlays or lower density,
   - smaller window.
3) If it persists on moderate data: try a different browser and update GPU drivers.

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
- **Grid density** increases 3D volume resolution (major GPU memory + compute).
- **Ray quality** increases raymarch steps (major GPU compute).

### How to confirm
- Reduce grid density first; if FPS recovers immediately, you found the bottleneck.

### Fix
- Lower grid density and ray quality.
- Lower render resolution (if exposed).
- Return to points mode when doing scientific work; use smoke selectively.

### Prevention
- Treat smoke mode as “presentation”: start low and scale up slowly.

Related docs:
- {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke`

---

## Edge case: vector field overlay + many views (pixels × views × post-processing)

### What you see
- Overlay is fine in one view but unusable in grid view.
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
3) Reduce number of views and retry.

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
With Live filtering enabled, each slider movement can trigger an expensive recomputation.

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
- Try a private window (no extensions) or a different browser profile.
- Check `chrome://gpu` (hardware acceleration).

### Fix
- Disable interfering extensions for the site.
- In corporate environments, ask IT to allow hardware acceleration and required origins.

### Prevention
- For enterprise deployments: document supported browsers and required policies.

Related docs:
- {doc}`../a_orientation/02_system_requirements`

---

## Next steps

- {doc}`07_troubleshooting_performance` (symptom → diagnosis → fix)
- {doc}`09_reporting_performance_bugs` (how to report a performance cliff with the right context)
