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
but the local H5AD file, Zarr ZIP, and prepared-export workflows use the same
validation contract on every supported platform.

### Hardware (what matters)

Cellucid’s limiting resource is usually the **GPU** (not CPU):

- **GPU memory (VRAM)**: more VRAM = fewer “context lost” crashes on large datasets.
- **GPU driver quality**: outdated drivers can cause rendering glitches or failures.
- **System RAM**: helps for very large datasets and smoke mode (but does not replace VRAM).

Rule of thumb:

- If you can smoothly play modern 3D web content, you can usually run Cellucid demo datasets.
- For “millions of points + lots of interaction”, a discrete GPU is strongly recommended.

### Input devices (what is expected)

- Mouse or trackpad recommended.
- Keyboard shortcuts are available (see the in-app “Keyboard Shortcuts” accordion).
- Touch devices can work for basic orbit/planar navigation, but are much more likely to run out of GPU memory.

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
- or use a supported workflow (e.g., a local server hosting the export folder).

---

## Notebook / iframe embedding notes (pointer lock + fullscreen)

If Cellucid is embedded in an iframe (e.g., in Jupyter):

- **Pointer lock** (“Capture pointer”) may be blocked unless the iframe allows it.
- **Fullscreen** may be blocked unless the iframe allows it.

Symptom: clicking “Capture pointer” does nothing, or fullscreen immediately exits.

Fix: use the standalone web app (not embedded) or adjust iframe permissions.

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
- Click **Reset Camera**.

**Fix**

1) Refresh the page once (sometimes the GPU recovers).
2) Use **Reset Camera** and clear filters.
3) If it persists, confirm hardware acceleration and WebGL2 are enabled, then
   update the browser, operating system, and GPU driver when applicable.

**Prevention**

- Avoid enabling smoke mode at very high grid sizes on laptops.

---

### Symptom: “WebGL context lost. Reload required…”

Cellucid shows an overlay and requires a reload when the GPU context is lost.

**Likely causes (ordered):**

1) GPU out of memory (too many points, too many views, or high smoke settings).
2) Driver crash/reset (especially on older integrated GPUs).
3) Another GPU-heavy app/tab causing pressure (video calls, 3D apps).

**How to confirm**

- It often happens right after a heavy action:
  - switching to smoke with a high grid,
  - increasing ray quality,
  - opening many views/snapshots,
  - or loading a very large dataset.

**Fix (safe → aggressive)**

1) Reload the page.
2) After reload, reduce GPU load:
   - avoid smoke mode or lower grid/quality,
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
- Lower `Point size (log):` and avoid very large point sizes.
- Turn off `Antialiasing (smooth point edges)` under **Image quality** and
  reload — measured at roughly a fifth of the frame time at the default point
  size on one Apple M1 Pro, at the cost of visibly harder point edges.
- If using smoke: reduce grid density, ray quality, and render resolution.

Lowering `Shader quality:` is a reasonable-sounding move that does **not** help:
the three qualities measured the same at every point size tested. It changes how
points look, not how fast they draw. See
{doc}`../n_benchmarking_performance/02_performance_considerations_what_gets_slow_and_why`.

---

## Next steps

- If the app loads: go to {doc}`03_quick_tour_60_seconds`.
- If navigation feels confusing: go to {doc}`../c_core_interactions/01_navigation_modes_orbit_planar_free_fly`.
