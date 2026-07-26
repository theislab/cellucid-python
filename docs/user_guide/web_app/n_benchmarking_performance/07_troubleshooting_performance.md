# Troubleshooting (performance)

This page is organized by **symptom → diagnosis → fix**.

If you want the conceptual model first, read {doc}`01_performance_mental_model`.

If you want a catalog of “performance cliffs”, read {doc}`06_edge_cases_performance`.

---

## Troubleshooting template (use this structure)

When adding a new entry, use:

### Symptom
What the user sees (include on-screen text if possible).

### Likely causes (ordered)
3–7 plausible causes, each testable.

### How to confirm
Concrete checks (UI actions, view count, toggles, browser tools).

### Fix
Step-by-step actions (safe fixes first).

### Prevention
What to do earlier to avoid this next time.

---

## Symptom: “Navigation is choppy / low FPS when I move the camera”

### Likely causes (ordered)
1) GPU-bound rendering (too many pixels, too many views, heavy shaders).
2) GPU-heavy modes enabled (smoke, vector overlays, heavy post-processing).
3) Hardware acceleration is disabled or GPU drivers are unstable.

### How to confirm
- Make the window smaller. If it becomes smoother immediately, you’re pixel/GPU-bound.
- Clear snapshots (go to single view). If it becomes smoother immediately, you’re view/GPU-bound.
- Disable overlays/smoke. If it becomes smoother immediately, you’re GPU-bound.
- In the browser console, confirm
  `document.createElement("canvas").getContext("webgl2") !== null`; then use
  the browser's graphics diagnostics to confirm hardware acceleration.

### Fix
1) Reduce **views**: clear snapshots; keep 1–3 views.
2) Reduce **pixels**: keep the window smaller while exploring.
3) Use **Points** mode (avoid smoke except intentionally).
4) If using vector overlays: reduce density and disable bloom first.  
   See {doc}`../i_vector_field_velocity/05_performance_and_quality`.
5) If you still see low FPS on moderate datasets, update the browser,
   operating system, and GPU driver when applicable.

### Prevention
- Keep a “laptop-safe preset” (points mode, few views, no heavy overlays).
- Treat smoke mode as cinematic/presentation, not a default analysis mode.

---

## Symptom: “Filtering is slow / sliders lag / the UI stutters when I filter”

### Likely causes (ordered)
1) CPU-bound visibility recomputation (especially with Live filtering on).
2) Too many enabled filters (work scales with `n_cells × n_filters`).
3) Multiview multiplier (filters/updates interact with many views).

### How to confirm
- Turn Live filtering off and click `FILTER` once; if it becomes smooth, you were stuck in a recompute loop.
- Disable all but one filter; if performance returns, filter stacking was the multiplier.
- Clear snapshots and retry; if performance improves, view count was amplifying the cost.

### Fix
1) For continuous fields: **Live filtering off → adjust sliders → click FILTER once**.  
   See {doc}`../d_fields_coloring_legends/03_color_by_behavior`.
2) Disable/remove no-op filters (filters that don’t change counts given the other filters).
3) Keep multiview lean while tuning filters; snapshot after you stabilize.
4) For very large datasets: do coarse gating upstream in Python; treat UI filtering as refinement.

### Prevention
- Teach “apply once” workflows for large datasets.
- Avoid using the UI as the primary data-cleaning layer on million-cell data.

Related docs:
- {doc}`../e_filtering/05_performance_considerations`
- {doc}`../e_filtering/07_troubleshooting_filtering`

---

## Symptom: “Loading a dataset freezes / takes forever / tab crashes”

### Likely causes (ordered)
1) You are loading a large `.h5ad` directly in the browser (memory cliff; not truly lazy).
2) Disk/network I/O is slow (network drive, remote mount, throttled server).
3) The dataset is too large for the current machine/GPU (VRAM pressure during initialization).
4) Browser extensions or corporate policies interfere with file APIs or GPU acceleration.

### How to confirm
- If it’s a `.h5ad` browser load: try the same file via Server Mode; if it works, the in-browser load path was the issue.
- If it’s remote: open Network tab and look for long requests / timeouts.
- Check Task Manager / Activity Monitor: if memory usage spikes rapidly, you’re hitting a memory cliff.

### Fix
1) Prefer **Server Mode** for large `.h5ad` / huge `.zarr`.  
   See {doc}`../b_data_loading/04_server_tutorial`.
2) Move data to a fast local disk (avoid network mounts when possible).
3) Reduce dataset complexity (fewer fields, fewer categories) or use a lite export.
4) Try a clean browser profile / private window to rule out extensions.

### Prevention
- Treat in-browser `.h5ad` load as “small quick preview only”.
- For large teams, standardize on a supported loading workflow (server mode or exports).

Related docs:
- {doc}`../b_data_loading/03_browser_file_picker_tutorial`
- {doc}`../b_data_loading/01_loading_options_overview`

---

## Symptom: “Switching genes is slow / gene coloring pauses”

### Likely causes (ordered)
1) Gene expression is loaded on demand (cold cache).
2) The data is being fetched over a remote server/network.
3) The gene/feature set is extremely large and search/metadata is heavy.

### How to confirm
- Switch to the same gene again; if it’s faster the second time, you’re seeing cache warming.
- If remote: watch Network requests while switching genes; long latency suggests I/O-bound.

### Fix
1) Expect the first gene load to be slower; benchmark cold vs warm explicitly.
2) For remote workflows: move the server closer (LAN) or use local server mode when possible.
3) If gene switching is consistently slow: reduce other multipliers (views, overlays) while exploring expression.

### Prevention
- For large datasets: choose a workflow that supports efficient lazy-loading (server mode).

Related docs:
- {doc}`../b_data_loading/04_server_tutorial`
- {doc}`04_benchmarking_methodology_and_metrics` (cold vs warm measurement)

---

## Symptom: “WebGL context lost” or “the canvas goes blank”

### Likely causes (ordered)
1) GPU out of memory (too many points/views, smoke grid too high, overlay too heavy).
2) GPU driver crash/reset (unstable drivers, integrated GPUs under pressure).
3) Corporate/managed environment disables/blocks hardware acceleration.

### How to confirm
- Did it happen immediately after enabling smoke, increasing quality, or opening many views? That points to VRAM pressure.
- Check `chrome://gpu` / `about:support` to confirm hardware acceleration and WebGL2.

### Fix (safe → aggressive)
1) Reload the page (context loss requires reload).
2) After reload, reduce GPU load:
   - points mode (avoid smoke),
   - fewer views,
   - disable overlays or reduce density/bloom,
   - smaller window.
3) If it still happens on moderate datasets:
   - confirm hardware acceleration and WebGL2 are enabled,
   - update the browser, operating system, and GPU driver when applicable.

### Prevention
- Keep a conservative default preset for laptops and workshops.

Related docs:
- {doc}`../a_orientation/02_system_requirements` (context lost troubleshooting)

---

## Symptom: “It’s fine in single view but slow in grid view”

### Likely causes (ordered)
1) You are simply multiplying GPU work by `n_views`.
2) Overlays/post-processing scale with pixels × views.

### How to confirm
- Clear snapshots; if performance recovers, the view multiplier is the culprit.

### Fix
- Keep grid view to the minimum number of comparisons you need.
- Disable overlays while in grid view, or lower overlay density and bloom.

### Prevention
- Adopt the workflow: tune in live view → snapshot for comparison.

---

## Symptom: “Performance gets worse over time”

### Likely causes (ordered)
1) Thermal throttling (laptops).
2) Memory growth (large caches, large sessions, possible leak).
3) Many accumulated views/highlights/analysis artifacts in one session.

### How to confirm
- If the laptop gets hot and performance degrades gradually: thermals are likely.
- Watch memory in Chrome Task Manager; if it steadily rises without stabilizing, investigate memory behavior.
- Clear snapshots/highlights and see if responsiveness returns.

### Fix
1) Reduce sustained GPU load (fewer views, smaller window, no smoke/overlays).
2) Reload the page to reset memory state.
3) Avoid accumulating lots of snapshots and massive highlights in one long-running session.

### Prevention
- For long sessions, keep conservative presets and take periodic “reset” breaks (reload + restore a saved session).

---

## Symptom: “It’s slow on one machine/browser but fine on another”

### Likely causes (ordered)
1) Different GPU class (integrated vs discrete), different VRAM.
2) Hardware acceleration disabled on the slow machine.
3) Driver differences (outdated or unstable GPU drivers).
4) Browser extension/policy interference.

### How to confirm
- Compare the browser-neutral WebGL2 console check and the browser's graphics
  diagnostics on both machines.
- Run the synthetic benchmark with the same preset and compare FPS.
  See {doc}`05_benchmark_tools`.

### Fix
- Enable hardware acceleration in the current stable Chrome, Edge, Firefox, or
  Safari release.
- Update drivers/OS on the slow machine.
- Use a lighter export/workflow for weaker machines.

### Prevention
- For teams/workshops, publish a minimum supported environment checklist.

---

## If you still can’t resolve it

At this point, treat it as a reportable issue.

Go to {doc}`09_reporting_performance_bugs` and include:
- your scenario,
- dataset size and loading method,
- hardware/browser,
- and at least one number (TTFR, FPS, or filter apply time).
