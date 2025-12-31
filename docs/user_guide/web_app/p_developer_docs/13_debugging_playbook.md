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
- server returned HTML for a JSON/binary request (SPA fallback or error page).

Fix:
- Open the failed request URL in a new tab; verify it’s real JSON/binary and returns correct status codes.

Deep dive:
- {doc}`09_data_loading_pipeline_and_caching`

---

## Fields/legend issues

### Symptom: “Field selector is empty”

Likely causes (ordered):
1) `obs_manifest.json` missing/unreadable.
2) Dataset loaded with “empty obs” fallback (dev-phase behavior).
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

Likely causes:
- transparency/visibility updated in state but not pushed to viewer.

How to confirm:
- Put a breakpoint in `viewer.updateTransparency`.
- Check `window._cellucidState.categoryTransparency` length and values.

Fix:
- Ensure filter manager calls `viewer.updateTransparency(...)` (state→viewer sync is state-owned).

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
- GPU memory pressure, especially with smoke or overlays on large datasets.

How to confirm:
- DevTools console logs mention context loss.

Fix:
- Reduce dataset size and reproduce.
- Disable smoke/velocity and re-test to isolate.

Deep dive:
- {doc}`07_rendering_pipeline_webgl_and_performance_notes`
- {doc}`../n_benchmarking_performance/index`

---

## Session issues (save/load)

### Symptom: “Session loads but key state is missing”

Likely causes:
- the feature is intentionally excluded
- heavy artifacts are still loading (lazy chunks)
- dataset mismatch caused dataset-dependent chunks to skip

How to confirm:
- Watch notifications during lazy restore.
- Check console for “dataset mismatch” skips.

Fix:
- Confirm you loaded the same dataset id/fingerprint.
- If feature should persist, add a session contributor and document it.

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
