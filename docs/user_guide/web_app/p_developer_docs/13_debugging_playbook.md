# Debugging playbook

This page is the “if it’s broken, do this” playbook for the Cellucid web app.

It is intentionally **procedural**: follow it top-to-bottom when you hit a bug, and you will usually end up with a minimal reproduction and the exact clues needed to fix it.

## At a glance

**Audience**
- Wet lab / non-technical: use this to capture actionable bug reports.
- Computational users: use this to diagnose data/loading/field issues quickly.
- Developers: use this to debug state/rendering/persistence regressions.

**Time**
- 5–10 minutes for most issues to get to a good bug report
- 20–60 minutes for complex performance/memory issues

---

## Step 1: Record the environment (always)

Before you click anything else, write down:

1) **Where you ran Cellucid**
   - Hosted: `https://www.cellucid.com`
   - Local app: `http://localhost:8000` (your own `cellucid/` checkout)
   - Embedded viewer in Jupyter (iframe)

2) **Browser + OS**
   - Browser name + version
   - OS name + version

3) **Dataset source**
   - local-demo / local-user / remote server / GitHub / Jupyter

4) **Dataset identity**
   - dataset id (from the UI dataset info block or from `dataset_identity.json`)

Why this matters:
- the same bug can be caused by CORS in one environment and by caching in another.

---

## Step 2: Turn on debug logging (recommended)

### App-wide debug (recommended default)

```js
localStorage.setItem('CELLUCID_DEBUG', 'true');
location.reload();
```

This enables:
- the logger in `cellucid/assets/js/utils/debug.js`
- many `[Main]`, `[UI]`, and source-selection logs

### Analysis debug (only when debugging analysis)

Either:
- add `?debug=1` to the URL, or
- set `localStorage.setItem('debug', '1')` and reload

This enables debug utilities in:
- `cellucid/assets/js/app/analysis/shared/debug-utils.js`

---

## Step 3: Identify the subsystem (pick one)

Use the symptom list to choose your debugging path:

| Symptom | Subsystem | Go to |
|---|---|---|
| Blank page / module load error | Boot / deployment | “Boot issues” |
| Dataset does not load / loads forever | Data loading | “Data loading issues” |
| Field selector empty / wrong colors | Fields + legends | “Fields/legend issues” |
| Filters don’t behave / count wrong | Filters + visibility | “Filtering issues” |
| Selection/highlights wrong | Highlights + tools | “Highlight issues” |
| Dimension switch breaks | Dimension manager + viewer positions | “Dimension issues” |
| WebGL context lost / slow | Rendering/GPU | “Rendering issues” |
| Session save/load wrong | Session serializer | “Session issues” |
| Analysis panel broken | Analysis module | “Analysis issues” |
| Figure export broken | Figure export | “Export issues” |
| Community annotation broken | GitHub auth/sync | “Community annotation issues” |

---

## Boot issues (page won’t start)

### Symptom: “Module script failed to load”

Likely causes (ordered):
1) You opened via `file://` instead of HTTP.
2) Server sends wrong MIME type for `.js`.
3) Deployed under a subpath and module URLs don’t resolve.

How to confirm:
- DevTools → Console: look for the module load error.
- DevTools → Network: click the failing `.js` request and check headers.

Fix:
- Run `python -m http.server` from `cellucid/` (see {doc}`02_local_development_setup`).
- Fix server MIME types (see {doc}`03_build_run_and_deployment`).

### Symptom: “WebGL2 is required but not supported”

Likely causes:
- WebGL disabled by browser or corporate policy
- remote desktop environment without GPU support

Fix:
- Try another browser/machine; confirm WebGL2 is enabled.

---

## Data loading issues

### Symptom: “Dataset loads forever”

Likely causes (ordered):
1) CORS/mixed-content blocked a fetch.
2) Wrong path/manifest missing and the UI is waiting on a fetch that never completes.
3) Very large file in browser file picker mode (no lazy loading).

How to confirm:
- Network tab: which request is pending/blocked?
- Console: enable `CELLUCID_DEBUG` and look for the last “loading…” log.

Fix:
- For large h5ad/zarr, use remote server mode (cellucid-python).
- For GitHub datasets, verify raw URLs are reachable.
- For CORS, run the local app and connect from the same origin when possible.

### Symptom: “Unexpected token `<`” during load

Likely cause:
- the server returned an HTML document for a JSON or binary request.

Fix:
- Open the failed request URL in a new tab; verify it’s real JSON/binary and returns correct status codes.

### Read the failure message before theorising

Four distinct load failures used to report themselves identically. They no
longer do, and the wording identifies the cause:

| Message contains | Cause | Error code |
|---|---|---|
| “was cancelled” | An abort — re-thrown untouched, never relabelled | (the original `AbortError`) |
| “is larger than the metadata size ceiling” | A bounded reader refused the payload | `VALIDATION_ERROR` |
| “response body transfer failed” | The connection dropped mid-body | `NETWORK_ERROR` |
| “must contain valid JSON” | Genuinely malformed content | `INVALID_FORMAT` |

The connection and catalogue UI classifies these by **error code and HTTP
status**, never by matching message text, so a message you rephrase will keep
its user-facing sentence — but a code you get wrong will silently move the
failure into the wrong bucket. Set the code deliberately.

Deep dive:
- {doc}`09_data_loading_pipeline_and_caching`

---

## Fields/legend issues

### Symptom: “Field selector is empty”

Likely causes (ordered):
1) `obs_manifest.json` missing/unreadable.
2) `obs_manifest.json` contains no usable fields.
3) You are in a snapshot view with a different view context.

How to confirm:
- Network tab: did `obs_manifest.json` return 200?
- Console: check `window._cellucidState.pointCount` and whether a field loader exists.

Fix:
- Fix exports to include `obs_manifest.json` (or use `prepare()`).
- Ensure the correct dataset base URL is in use.

### Symptom: “Colors look wrong / legend doesn’t match”

Likely causes:
- field change didn’t trigger a color recompute
- colormap override or category rename/delete registry is involved

How to confirm:
- Subscribe to `field:changed` events and see what fires when you change a field.
- Inspect `window._cellucidState.getActiveField()` in console.

Fix:
- Ensure state owns the recompute and calls `viewer.updateColors(...)`.

---

## Filtering issues

### Symptom: “Filter count changes but points don’t”

Likely causes (ordered):
1) transparency/visibility updated in state but not pushed to the viewer;
2) it *was* pushed, but the alpha bytes did not change, so the upload was
   correctly skipped — the count you are reading is derived from something else;
3) the change is below the visible-alpha threshold, so it rounds to the same
   byte.

How to confirm:
- Put a breakpoint in `viewer.updateTransparency` and look at its **return
  value**: it reports whether the published alpha generation actually moved.
  `false` means the bytes were identical, and no upload was needed.
- Check `window._cellucidState.categoryTransparency` length and values.

Fix:
- Ensure the filter manager calls `viewer.updateTransparency(...)` (state→viewer
  sync is state-owned).
- If the return value is `false` but you expected a change, the bug is upstream
  in how the alpha array was derived, not in the upload.

---

## Highlight issues

### Symptom: “Highlighted cells disappear after filtering”

Clarify expected behavior:
- Highlights are stored independently of visibility; filtering hides points but does not remove them from highlight groups.

If behavior differs:
- Confirm whether you are looking at “visible highlight count” vs “total highlight membership”.

How to confirm:
- Check highlight pages and group membership in `window._cellucidState.getHighlightPages()`.

---

## Dimension issues

### Symptom: “Dimension switch changes UI but not the embedding”

Likely causes (ordered):
1) Embedding missing for that dimension.
2) Fetch failed (network/CORS).
3) Positions buffer length mismatch caused viewer to reject it.

How to confirm:
- Network tab: embedding fetches.
- Console: look for `[Viewer] updatePositions: position count mismatch`.

Fix:
- Verify dataset identity + embeddings metadata.
- Validate positions array length before calling `viewer.updatePositions`.

---

## Rendering issues (GPU / performance)

### Symptom: “WebGL context lost”

Likely causes:
- GPU memory pressure, especially with smoke or overlays on large datasets,
  or a very large highlight group (highlight geometry costs bytes per selected
  cell **per view**).

How to confirm:
- DevTools console logs mention context loss.

Fix:
- Reduce dataset size and reproduce.
- Disable smoke/velocity and re-test to isolate.
- Clear highlights and re-test.

Deep dive:
- {doc}`07_rendering_pipeline_webgl_and_performance_notes`
- {doc}`../n_benchmarking_performance/index`

### Symptom: “Frames got slower and I cannot see why”

Start from the baseline: **an idle frame uploads nothing.** If your change made
frames slower, the most likely cause is that you put work on the draw path
without a gate.

How to confirm:
1) Open the **Performance Benchmark** panel once, so the harness module is
   published on `window._cellucidBenchmarkHarness`.
2) Build a harness against the live viewer and canvas. It shadows the GL context
   and counts, per frame: buffer allocations and bytes, sub-updates and bytes,
   texture uploads and bytes, draw calls, sync stalls, and `getError` /
   `finish` / `flush` / `readPixels` calls.
3) Hold the camera still. Every upload counter should be **zero**. A non-zero
   count while idle is the bug.
4) Then move the camera. Element-buffer uploads should appear only when the
   admitted node set changes, not on every frame.

Fix:
- Gate the new work behind a dirty flag, a cache comparison, or the allocation
  watermark, following the three gates described in
  {doc}`07_rendering_pipeline_webgl_and_performance_notes`.

:::{note}
The headline FPS and frame-time tiles in the benchmark panel are CPU wall-clock
measurements between frames, so they include everything else the tab is doing.
Only **Analyze Performance** uses real GPU timer queries, and only where the
driver provides the extension. Do not attribute a frame-time change to the GPU
on the strength of the headline tiles alone.
:::

### Symptom: “A frame rate measurement is not reproducible”

Check, in this order, that between the two runs you did not change: the window
size, the **layout** (one view / two- or three-view row / 2×2 grid — these are
three different workloads, not one curve), the filter state (hidden cells are
rejected before rasterisation and cost nothing), or the machine's background
load. See {doc}`../n_benchmarking_performance/04_benchmarking_methodology_and_metrics`.

---

## Session issues (save/load)

### Symptom: “Session loads but key state is missing”

Likely causes:
- the feature is intentionally excluded
- the one awaited restore is still running lazy chunks
- exact dataset identity or chunk validation rejected and rolled back the operation

How to confirm:
- Require terminal **Session fully restored**.
- Check the exact Console error and progress failure; there is no skip/partial mode.

Fix:
- Confirm you loaded the same dataset id/fingerprint.
- If feature should persist, add it to the closed current profile with
  completeness, rollback, and documentation tests.

Deep dive:
- {doc}`10_sessions_persistence_and_serialization`

---

## Analysis issues

### Symptom: “Analysis panel is blank”

Likely causes:
- analysis init threw
- Plotly load blocked

How to confirm:
- Check `window._comparisonModule` exists.
- Console error stack trace.

Deep dive:
- {doc}`11_analysis_architecture`

---

## Export issues

### Symptom: “SVG export freezes”

Likely cause:
- Full vector export on too many points.

Fix:
- Use optimized-vector or hybrid.

Deep dive:
- {doc}`12_figure_export_architecture`

---

## Community annotation issues

### Symptom: “Sign-in succeeds but repos list is empty”

Likely causes:
- GitHub App not installed on the org/repo
- wrong user signed in (multiple accounts)

How to confirm:
- In the UI: which GitHub user is shown?
- Confirm app installation on the target repo owner.

Fix:
- Install the GitHub App on the repo owner.
- Reconnect and Pull.

### Symptom: “Editing in two tabs causes disconnect”

This is intentional:
- a scope lock prevents silent data loss from concurrent edits.

How to confirm:
- You’ll see a notification that the lock was lost.

Fix:
- Close extra tabs for the same dataset/repo/user scope.

Developer references:
- `cellucid/assets/js/app/community-annotations/REPO_SETUP.md`
- `cellucid/assets/js/app/community-annotations/scope-lock.js`

---

Next: {doc}`14_testing_ci_and_release_process` (how to validate changes before sharing).
