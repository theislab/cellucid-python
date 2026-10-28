# System requirements

**Audience:** everyone (especially if anything is blank/slow)  
**Time:** 5–15 minutes (depending on troubleshooting)  
**What you’ll learn:**
- What Cellucid *requires* (WebGL2) vs what is “nice to have”
- How to diagnose browser/GPU issues quickly
- How to interpret common failure modes (blank canvas, context lost)

---

## Hard requirement: WebGL2

Cellucid’s viewer is **WebGL2‑only**. If your browser/device cannot create a WebGL2 context, the app cannot run.

Common on-screen message (exact):

- `WebGL2 is required but not supported in this browser.`

If you hit this, your fastest fix is usually:

1) update to the current stable release of your browser,
2) enable browser hardware acceleration,
3) update GPU drivers / OS,
4) disable strict “graphics blocking” policies (common in managed/corporate environments).

---

## Supported desktop environments

Cellucid supports the current stable releases of:

| Browser | macOS | Windows | Linux |
|---|---:|---:|---:|
| Chrome | Yes | Yes | Yes |
| Edge | Yes | Yes | Yes |
| Firefox | Yes | Yes | Yes |
| Safari | Yes | — | — |

All supported browsers must provide WebGL2 and the native
`DecompressionStream('gzip')` API. Browser-native file pickers differ visually,
but opening a file or a folder (**Open… › File…** or **Folder…**, or a drop)
reads it the same way on every supported platform.

### Hardware (what matters)

Cellucid’s limiting resource is **memory**, in two places — the browser tab's
JavaScript memory and the GPU — not CPU speed:

- **System RAM** (the tab's memory): about 18–19 bytes per cell for a 2D
  embedding and 23–24 for 3D once loaded, plus about 4 bytes per cell for each
  field or gene you colour by.
- **GPU memory (VRAM)**: about 13 bytes per cell in 2D and 17 in 3D. On
  computers with unified memory (Apple silicon) both come out of the same RAM.
- **At the largest sizes, count the whole footprint.** With levels of detail
  on, the browser's own overheads come on top: measured in Chromium on Apple
  silicon at 100 million cells, the tab and its GPU process held 54 bytes a
  cell in 3D once drawn (5.4 GB), 78 while the spatial order of its levels
  of detail was taken in, and 50 after (5.0 GB); 53 in 2D. Cells Cellucid
  wrote in priority order (or a subsample) draw their levels of detail at
  once, before that order exists
  ({doc}`../n_benchmarking_performance/05_benchmark_tools` has every
  measurement). Opening a dataset larger than the device holds, Cellucid
  opens a subsample sized to fit (77 bytes a cell against half the device's
  memory: about 104 million cells on a 16 GB machine).
- **GPU driver quality**: outdated drivers can cause rendering glitches or failures.

Rule of thumb:

- If you can smoothly play modern 3D web content, you can usually run Cellucid demo datasets.
- For “millions of points + lots of interaction”, a discrete GPU is strongly recommended.

**Multisampling is optional, not required.** Antialiasing is drawn by the app
rather than requested from the browser, so a device that cannot multisample still
runs Cellucid — it just draws harder point edges. When that is the case, the
`Antialiasing` (i) tip reads (exact):

- `This browser is not providing antialiasing for this view.`

There is nothing to fix in the app if you see it; the setting has no effect on
that device.

### Input devices (what is expected)

- Mouse or trackpad recommended.
- Keyboard shortcuts are available (see the in-app “Keyboard shortcuts” accordion).
- Touch devices work for navigation and selection (one finger orbits — pans in
  planar mode — two fingers pinch to zoom and drag to pan), but are much more
  likely to run out of GPU memory.
- **Volumetric smoke** needs float blending (`EXT_float_blend`), which Safari
  on iPhone and iPad and many Android browsers do not offer. There, choosing
  smoke shows “Volumetric smoke isn't available on this device: …” and the view
  stays in Points; everything else works.
- On a phone-sized screen (720 px wide or less) the sidebar opens over the whole
  screen. Tap **✕** in its top-right corner to close it and see the plot; tap
  **☰** in the top-left corner to bring it back.
- A camera frames the same data on every screen shape: the 45° field of view
  applies to the plot's shorter side, so on a phone held upright a sample's
  preset or the default view shows across what a landscape screen shows down,
  and nothing is cut off at the sides. A landscape screen is unchanged.
- On a touch screen the small icon buttons grow to fingertip size, and the
  controls that otherwise wait for a mouse to hover over them (the analysis
  rows' copy buttons, the highlight page tabs' buttons, plot toolbars) are always
  shown. An undocked panel moves by dragging its header with a finger.
- Selecting cells (highlighting) is an {kbd}`Alt` + click or
  {kbd}`Alt` + drag gesture with a keyboard. A touch screen has no
  {kbd}`Alt`, so turn on **Select mode** instead: its switch sits in the
  **Highlighting** panel and, on a touch screen, on the selection bar at the
  bottom of the plot. While it is on, a one-finger tap or drag runs the active
  highlight tool exactly as {kbd}`Alt` + click or {kbd}`Alt` + drag does, and
  two fingers still zoom and pan. Its **Intersect / Add / Subtract** buttons
  stand in for {kbd}`Shift` and {kbd}`Ctrl`, and the bar also confirms,
  cancels and clears selections, so the menu can stay closed. Turn the switch
  off (or press {kbd}`Esc` with nothing in progress) to go back to
  navigating. See {ref}`selection-bar`.


::::{grid} 1 2 2 2
:gutter: 2

:::{grid-item}
```{figure} ../../../_static/screenshots/web_app/phone-sidebar-open.png
:alt: Cellucid on a phone-sized touch screen with the Pancreas sample open: the sidebar fills the screen, with its close button (✕) in a strip at the top right above the Session section.
:width: 780px

On a phone the sidebar is the whole screen; **✕** at the top right closes it.
```
:::

:::{grid-item}
```{figure} ../../../_static/screenshots/web_app/phone-plot.png
:alt: The same phone screen after tapping ✕: the whole Pancreas sample in view on the plot, the menu button (☰) at the top left, and the selection bar with its Select mode switch at the bottom.
:width: 780px

Closed, the plot fills the screen with the whole sample in view: **☰** brings
the sidebar back, and the selection bar at the bottom holds **Select mode**.
```
:::

::::

---

## Quick self-check (2 minutes)

Do these in order:

1) **Confirm WebGL2 works**
   - In any supported browser console, evaluate
     `document.createElement("canvas").getContext("webgl2") !== null`.
   - The result must be `true`; use the browser's graphics diagnostics to
     confirm hardware acceleration.

2) **Try the app in a clean state**
   - Disable “aggressive” ad blockers for the site.
   - Try a private window to rule out extension interference.

3) **If the app loads but is very slow**
   - Reduce smoke quality (if in smoke mode) or switch back to points mode.
   - Close other GPU-heavy tabs (3D, video calls, other WebGL apps).

---

## Corporate / managed environments (common pitfalls)

In locked-down environments, “the browser works” does not always imply “WebGL2 works”.

Common blockers:

- Group policies disabling hardware acceleration
- Strict Content Security Policy (CSP) that blocks required resources
- Network proxy rewriting responses (can corrupt binary/session downloads)
- Disabled local file access (affects some local loading workflows)

If you are in this situation, you’ll usually need help from IT:

- allow hardware acceleration,
- allow required origins,
- or use a supported workflow (e.g., a local `cellucid serve` handing the file to the browser).

---

## Notebook / iframe embedding notes (pointer lock + fullscreen)

If Cellucid is embedded in an iframe (a notebook widget's output, or your own
`<iframe>` embed):

- **Pointer lock** (“Capture pointer”) may be blocked unless every enclosing
  frame allows it.
- **Fullscreen** may be blocked unless every enclosing frame allows it. The
  notebook widget's own frame allows fullscreen, but some notebook front ends
  (VS Code, Colab) render outputs inside frames of their own that may not.

Symptom: clicking “Capture pointer” does nothing, or fullscreen immediately exits.

Fix: use the standalone web app (not embedded) — for a notebook,
`cellucid.AnnDataViewer(...)` opens the same viewer in a browser tab — or
adjust the iframe's permissions.

---

## Troubleshooting (massive)

This section is intentionally long. Use the symptom that matches what you see.

### Symptom: “WebGL2 is required…” (app fails immediately)

**Likely causes (ordered):**

1) Browser does not support WebGL2 (old version, unusual browser).
2) Hardware acceleration disabled (browser setting or corporate policy).
3) GPU/driver blacklisted by the browser (stability protection).
4) Remote desktop / VM environment without GPU passthrough.

**How to confirm**

- Evaluate the browser-neutral WebGL2 console check above.
- Use the browser's graphics diagnostics to confirm hardware acceleration.

**Fix**

1) Update to the current stable release of Chrome, Edge, Firefox, or Safari.
2) Enable hardware acceleration in browser settings.
3) Update browser + GPU drivers.
4) If in a VM/remote desktop, switch to a machine with a real GPU or enable GPU passthrough.

**Prevention**

- For workshops, ask attendees to test WebGL2 the day before.

---

### Symptom: Blank canvas (UI loads, but the plot area is empty)

**Likely causes (ordered):**

1) WebGL context creation succeeded, but rendering failed due to GPU/driver instability.
2) The dataset is loaded but **all points are invisible** (filters/outlier threshold).
3) You are in smoke mode with parameters that effectively make the cloud invisible (very low density).

**How to confirm**

- Open the browser devtools console and look for WebGL errors.
- Check the app’s “filter count”/visibility indicators (if present).
- Switch back to **Points** render mode.
- Click **Reset camera**.

**Fix**

1) Refresh the page once (sometimes the GPU recovers).
2) Use **Reset camera** and clear filters.
3) If it persists, confirm hardware acceleration and WebGL2 are enabled, then
   update the browser, operating system, and GPU driver when applicable.

**Prevention**

- Avoid enabling smoke mode at high render resolution and ray quality on laptops.

---

### Symptom: “WebGL context lost. Reload required…”

Cellucid shows an overlay and requires a reload when the GPU context is lost.

**Likely causes (ordered):**

1) GPU out of memory (too many points, too many views, or high smoke settings).
2) Driver crash/reset (especially on older integrated GPUs).
3) Another GPU-heavy app/tab causing pressure (video calls, 3D apps).

**How to confirm**

- It often happens right after a heavy action:
  - switching to smoke at a high render resolution,
  - increasing ray quality,
  - opening many views/snapshots,
  - or loading a very large dataset.

**Fix (safe → aggressive)**

1) Reload the page.
2) After reload, reduce GPU load:
   - avoid smoke mode, or lower its render resolution, ray quality, and noise detail,
   - keep fewer views,
   - reduce point size,
   - close other GPU-heavy tabs.
3) If the problem repeats, reduce the dataset/rendering load or use a machine
   with sufficient GPU memory.

**Prevention**

- Treat smoke mode as “cinematic”: great visuals, but it can be GPU-expensive.
- For large datasets, start in points mode and scale up quality slowly.

---

### Symptom: Everything is slow / laggy (but not crashing)

**Likely causes:**

- Dataset is near your GPU’s limit.
- Points are large enough to cover a lot of the screen.
- You are in grid mode with many views.

**Fix**

- Reduce the number of views (clear snapshots) — note that a two- or three-view
  row costs more than the 2×2 four-view grid.
- Lower `Point size` and avoid very large point sizes.
- Turn off `Antialiasing` under **Performance** — it
  applies to the next frame, and was measured at roughly a fifth of the frame
  time at point size 0.75 on one Apple M1 Pro, at the cost of visibly harder
  point edges. Look at the box before you reach for it: until you click it, the
  choice is automatic from the cell count and the view's dimension — on below
  5,000,000 cells in 3D and below 200,000 in a 1D or 2D view, off at or above —
  so on a large dataset it is already clear and there is nothing left to turn off.
- If using smoke: reduce render resolution and ray quality, then grid density; or set animation speed to 0 (`Still`) so a still view is not redrawn.

Lowering `Shader quality` is a reasonable-sounding move that does **not** help:
the three qualities measured the same at every point size tested. It changes how
points look, not how fast they draw. See
{doc}`../n_benchmarking_performance/02_performance_considerations_what_gets_slow_and_why`.

---

## Next steps

- If the app loads: go to {doc}`03_quick_tour_60_seconds`.
- If navigation feels confusing: go to {doc}`../c_core_interactions/01_navigation_modes_orbit_planar_free_fly`.
