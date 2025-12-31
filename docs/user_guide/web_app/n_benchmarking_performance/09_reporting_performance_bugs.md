# Reporting performance bugs

This page is a **copy/paste playbook** for reporting performance problems in a way that is actionable.

Goal: minimize back-and-forth by including the context that matters:

- what you were doing,
- what hardware/browser you were on,
- how big the dataset is,
- and at least one number.

---

## Before you file (2-minute sanity checklist)

Performance problems are often workflow/environment issues. Before filing:

1) Confirm WebGL2 + hardware acceleration: {doc}`../a_orientation/02_system_requirements`.
2) Reduce views (clear snapshots) and retry.
3) Reduce pixels (smaller window) and retry.
4) Disable smoke mode and overlays and retry.
5) If filtering is slow: Live filtering off → `FILTER` once.
6) If loading is slow: switch to Server Mode for large `.h5ad`/`.zarr`.

If none of these change the behavior, it is likely worth reporting.

---

## Choose the right kind of report

### A) “Rendering is slow” (FPS / navigation)

Include:
- the synthetic benchmark result (if available), and/or
- your measured FPS and view/window settings.

See {doc}`05_benchmark_tools_if_exposed` and {doc}`04_benchmarking_methodology_and_metrics`.

### B) “Filtering/analysis is slow” (CPU hitches)

Include:
- dataset size (cells),
- number of enabled filters,
- whether Live filtering was on,
- and the time for a single “apply once” operation.

### C) “Loading is slow / freezes” (I/O + memory)

Include:
- loading mode (browser `.h5ad` vs server mode vs export folder),
- whether the slowdown is cold-cache only,
- and any network/storage constraints (remote server, network drive).

---

## Copy/paste template (recommended)

Paste this into a GitHub issue, email, or internal tracker.

### Summary

One sentence:
- “When I ___, Cellucid ___.”

Example:
- “When I apply a QC filter on a 1.2M-cell dataset, the UI freezes for ~8–12 seconds.”

### Scenario (what exactly were you doing?)

- What action triggers the problem?
- Is it consistent (every time) or intermittent?
- Is it worse in grid view (many snapshots) than in single view?
- Does window size change it?

### Expected vs observed

- Expected: `<what you expected>`
- Observed: `<what actually happens>`

### Dataset

If you can share details:
- Dataset name/ID: `<id>`
- Approx size: `<n_cells> cells`, `<n_genes> genes`, `<n_fields> fields`
- Notable complexity:
  - category counts for key categorical fields (if relevant),
  - number of views/snapshots in the session,
  - overlays enabled (which ones, settings).

If the dataset is private:
- provide the sizes and complexity numbers anyway,
- avoid sharing names/identifiers,
- and consider providing a small synthetic/minimal reproduction (see below).

### Loading method

- `<export folder / server mode / remote server / browser file picker>`
- If server/remote:
  - server location: `<local / LAN / remote>`
  - approximate latency/bandwidth if known

### Environment

- OS: `<macOS / Windows / Linux> <version>`
- Browser: `<Chrome/Edge/Firefox/Safari> <version>`
- GPU: `<e.g., M2 / RTX 3060 / Intel Iris>` (approx is fine)
- WebGL2 hardware acceleration: `<yes/no/unknown>`
- Window size: `<width × height>`, device pixel ratio if known

### Steps to reproduce (numbered)

1) …
2) …
3) …

### Measurements (at least one)

Pick what matches the issue:

- TTFR (load → first render): `<median>`, `<range>`, `<cold/warm>`
- FPS while navigating: `<approx>`
- Filter apply time (Live filtering off, click FILTER once): `<seconds>`
- “Analyze Performance” result from benchmark panel (if used): paste output

If the app exposes it:
- Paste **Copy Situation Report** from the **Performance Benchmark** panel.

### Attachments (optional but high value)

- Screenshot(s) of the relevant UI state (views count, toggles, warnings)
- Console log excerpt (errors/warnings)
- Network HAR file (if remote loading is involved)

### Workarounds you tried

Bullet list:
- “Cleared snapshots”
- “Smaller window”
- “Disabled overlays”
- “Server mode”
- “Different browser”

---

## How to get the most useful “machine + GPU” info (quick)

### Chrome / Edge

- WebGL status: open `chrome://gpu` and look for WebGL2 hardware acceleration.
- Per-tab metrics: open Chrome Task Manager with `Shift+Esc` (CPU/memory; sometimes GPU memory).

### Firefox

- WebGL status: open `about:support` and look for WebGL2 info.

---

## How to create a minimal reproduction (when data is private)

The goal is to preserve the performance characteristics without sharing sensitive data.

Options:

1) **Use the synthetic benchmark** to show baseline rendering performance (point count, FPS, GPU memory).
2) **Downsample + anonymize** a small subset that still triggers the issue (e.g., 200k cells instead of 2M).
3) **Share structure, not values**:
   - “Field X has 25k categories”
   - “We have 12 snapshots open”
   - “Overlay density was 50K”

Even without sharing data, these numbers let developers reason about the bottleneck.

---

## Redaction and privacy checklist

Before posting logs/screenshots:
- remove dataset names/IDs if private,
- remove local file paths (often contain usernames),
- remove private org/repo names (GitHub integration),
- remove tokens/URLs with secrets.

If you’re unsure, err on the side of redacting.

---

## Next steps

- {doc}`04_benchmarking_methodology_and_metrics` (how to produce comparable numbers)
- {doc}`05_benchmark_tools_if_exposed` (synthetic benchmark + “Copy Situation Report”)
- {doc}`07_troubleshooting_performance` (if you want to try more fixes first)
