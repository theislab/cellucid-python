# Troubleshooting (velocity overlay)

**Audience:** everyone using the overlay  
**Time:** find your symptom in 1–2 minutes; fixes take 5–15  
**What you’ll get:** the overlay working, or a precise statement of why it cannot.

This page is a troubleshooting-first companion to the vector field overlay docs.
For behaviour that looks wrong but is correct, go to {doc}`06_edge_cases` instead.

Each entry follows the same structure:
**Symptom → Likely causes (ordered) → How to confirm → Fix → Prevention**.

---

## Quick triage (60 seconds)

If you’re stuck, do this first:

1) Confirm **Visualization → Render mode:** is `Points`.  
2) Confirm you are on a dimension that actually has vector fields (2D vs 3D).  
3) Enable `Show overlay`, then temporarily set:
   - `Particle density:` 5K
   - `Trail length:` 3.0s
   - `Opacity:` 60%
   - `Bloom strength:` 0.00 (Advanced Visual Settings → HDR & Bloom)
4) If it’s still broken, look at which symptom below matches what you see.

---

## Symptom: “I can’t find the Vector Field Overlay controls”

### Likely causes (ordered)

1) You are in `Render mode: Volumetric smoke cloud (alpha)` (overlay controls are currently in the Points-only controls area).  
2) The dataset has **no vector fields**, so the overlay block is hidden.  
3) You are not looking under **Visualization** in the sidebar (UI location confusion).  

### How to confirm

- In the left sidebar, open **Visualization** and look for:
  - `Render mode:` and then
  - a block titled `Vector Field Overlay:`

### Fix

1) Set **Visualization → Render mode:** to `Points`.
2) If the overlay block still does not appear, your dataset likely has no vector fields:
   - load a dataset known to contain vector fields, or
   - add/export vector fields (see “How do I add vector fields?” below).

### Prevention

- When preparing datasets, include vector fields explicitly (preferred via `cellucid.prepare(..., vector_fields=...)`).

---

## Symptom: The overlay block exists, but `Show overlay` is disabled

### Likely causes (ordered)

1) Vector fields exist, but **not for the current dimension** (e.g., only 2D vectors exist and you are in 3D).  
2) The dataset only contains vector fields for a different basis/dimension than you expect.  

### How to confirm

- The UI will show a hint like:
  - `Vector fields available for 2D, 3D. Switch embedding dimension to enable.`

### Fix

1) Switch the current view to a supported dimension (2D or 3D).
2) If you need the missing dimension, export/generate that dimension’s vector field and reload the dataset.

### Prevention

- Export both `*_umap_2d` and `*_umap_3d` vector fields when you expect users to switch dimensions.

---

## Symptom: “Overlay toggle instantly turns off” (checkbox unchecks itself)

### Likely causes (ordered)

1) The vector field failed to load (file missing, wrong shape, CORS/permission issue).  
2) You enabled the overlay in a dimension where no field exists (it auto-disables to avoid nonsense).  
3) You switched datasets (the UI resets the overlay state on dataset changes).  

### How to confirm

- Watch the info line under the overlay settings:
  - it may briefly show `Loading vector field…`, then `Failed to load vector field.`
- You may also see a toast/notification containing an error message.

### Fix

1) Confirm you are in a supported dimension for this dataset.
2) Try a different field in `Vector field:` (if available).
3) If you are loading via server/public hosting, confirm the vector files exist and are reachable.
4) If you are loading an AnnData, confirm `adata.obsm` contains the expected key and shape (see “How do I add vector fields?”).

### Prevention

- Keep vector field naming and shapes consistent (dimension-specific, row-aligned).
- Prefer pre-exported folders created by `cellucid.prepare(...)` for reproducibility.

---

## Symptom: “No fields appear in dropdown” (empty `Vector field:`)

### Likely causes (ordered)

1) There are no vector fields for the current dimension (the dropdown is dimension-filtered).  
2) Your dataset metadata says vector fields exist, but the actual arrays/files are missing.  
3) The vector fields are present but do not match naming conventions, so they aren’t detected.  

### How to confirm

- Switch dimensions and see if entries appear in 2D or 3D.
- If you’re loading an AnnData, inspect `adata.obsm.keys()` in Python and look for:
  - `velocity_umap_2d`, `velocity_umap_3d`, `T_fwd_umap_2d`, etc.
- If you’re loading an export folder, check for:
  - `dataset_identity.json` containing a `vector_fields` block
  - a `vectors/` directory with files like `velocity_umap_2d.bin` (or `.bin.gz`)

### Fix

1) Switch to a dimension that has vector fields.
2) If the dataset should have vector fields but doesn’t:
   - regenerate/export them and reload the dataset.

### Prevention

- Use required dimension-suffixed keys: `<field>_umap_<dim>d` (e.g.,
  `velocity_umap_2d`).
- Validate shapes before export: `(n_cells, dim)` with matching row order.

---

## Symptom: Overlay is enabled, but nothing is visible on the canvas

### Likely causes (ordered)

1) `Opacity:` is 0% or very low, or `Particle size:` is too small.  
2) All cells are filtered out / invisible, so there is nothing to spawn particles from.  
3) Vectors are all-zero. Every particle is then discarded before it is drawn, so
   nothing appears at any setting — see {doc}`06_edge_cases`.  
4) Vectors are extremely small but non-zero, so motion is imperceptible.  
5) You are GPU-saturated and frames are not updating smoothly.  

### How to confirm

- Set `Opacity:` to 60% and `Particle size:` to ~2.0.
- Temporarily relax filters so you have visible points.
- Reduce `Particle density:` to 5K (to reduce GPU load).
- If possible in Python: check vector magnitudes (e.g., `np.linalg.norm(vectors, axis=1)`).

### Fix

1) Adjust `Opacity:` and `Particle size:`.
2) Ensure you have visible cells (undo filters that hide everything).
3) Increase `Flow speed:` slightly (2×–4×) to see motion.
4) If steps 1–3 change nothing at all, stop tuning: an all-zero field cannot be
   revealed by any control. Recompute the vectors or verify you are exporting
   the intended array.

### Prevention

- Include a “known-good” vector field in your demo dataset so users can sanity-check the visualization pipeline.

---

## Symptom: Overlay is too faint / too bright

### Likely causes (ordered)

- Too faint:
  1) low `Opacity:` or tiny `Particle size:`
  2) low `Intensity:` / low `Exposure:` (advanced)
  3) dark background + low contrast
- Too bright:
  1) high `Particle density:` + long `Trail length:`
  2) high `Trail persistence:` (advanced)
  3) `Bloom strength:` too high (advanced)

### How to confirm

- Switch to conservative “debug” settings:
  - `Particle density:` 5K
  - `Trail length:` 3.0s
  - `Opacity:` 60%
  - `Bloom strength:` 0.00

### Fix

- If too faint:
  1) raise `Opacity:` and `Particle size:`
  2) then raise `Intensity:` slightly (advanced)
  3) only then consider `Exposure:` changes
- If too bright:
  1) reduce `Particle density:`
  2) reduce `Trail length:`
  3) reduce `Bloom strength:` (or set to 0.00)
  4) reduce `Opacity:`

### Prevention

- Treat “cinematic” tuning as a separate preset; keep a conservative preset for scientific interpretation.

---

## Symptom: Performance drops drastically when I enable the overlay

### Likely causes (ordered)

1) Particle count is too high for your GPU (`Particle density:`).  
2) Trails and bloom are dominating (large window + high DPI + bloom).  
3) You have many snapshot views open (each view has its own trail buffers).  
4) LOD is disabled, or it is enabled and the camera is close enough that `Auto`
   is drawing full detail — in which case `Sync with LOD` has nothing to reduce
   until you pull back.  

### How to confirm

- Disable the overlay and see if FPS immediately recovers.
- Switch to single view (or clear snapshots) and compare.
- Shrink the browser window and compare.

### Fix

1) Set `Particle density:` to 5K–10K.
2) Set `Trail length:` to 2–5s.
3) Set `Bloom strength:` to 0.00 (Advanced → HDR & Bloom). This is not only a
   visual change: it also frees the bloom render targets and skips two
   full-screen passes.
4) Reduce number of visible views.
5) Enable renderer LOD (**Visualization → Renderer settings → Level-of-Detail
   (LOD)**), keep `Sync with LOD` enabled, and pull the camera back — `Auto`
   coarsens at any dataset size. Lower `Force LOD level:` by a few steps when you
   need to go past `Auto`'s floor of `min(2,000,000, cells ÷ 8)` points.

### Prevention

- For large datasets, document a “default performance-safe preset” and encourage users to start there.
- See {doc}`05_performance_and_quality` for what the LOD floor is and when the
  forced slider is the lever to reach past it.

---

## Symptom: “The flow vanishes when I zoom out”

### Likely causes (ordered)

1) `Sync with LOD` is enabled and the renderer has dropped six or more levels
   below full detail. At that point the overlay's particle multiplier is `0` and
   the overlay is disposed until you come back up.  
2) You are only seeing the intermediate steps of the same ladder — the flow
   thins rather than disappearing (`0.5` at one or two levels below full detail,
   `0.25` at three to five).  

### How to confirm

- Zoom back in. If the flow returns without you touching any control, this is
  the ladder, not a failure.
- Uncheck `Sync with LOD` and repeat the zoom. If the flow now survives, the
  ladder was the cause.

### Fix

- Uncheck `Sync with LOD` if you need the flow at every zoom level, and accept
  the frame cost, or
- raise `Force LOD level:` so the renderer stays nearer full detail.

### Prevention

- Treat `Sync with LOD` as a performance setting, not a display setting, and
  read {doc}`05_performance_and_quality` before using it in a figure workflow.

---

## Symptom: `Preparing velocity overlay...` takes a long time

### Likely causes (ordered)

1) Your dataset is huge and Cellucid is building internal spawn tables (especially after filter changes).  
2) You are changing filters faster than the build completes, so it restarts.  
3) The browser is too busy to schedule the background work promptly — the build
   runs on idle time, not on the render loop.  

### How to confirm

- Stop touching the filters for a few seconds and watch whether the notification settles.
- Check whether the main thread is busy with another long task (a large export, a DE run).

### Fix

1) Stop rapidly changing filters; wait a moment for the overlay to finish preparing.
2) Reduce GPU load (lower density, shorter trails, bloom=0).
3) If needed, toggle the overlay off and on once after you settle on filters.

### Prevention

- On very large datasets, keep density modest and settle your filters before
  enabling the overlay.

:::{note}
`Velocity overlay ready (no visible cells)` is **not** this symptom, and not a
failure. It is the ordinary outcome when every cell is filtered out: the spawn
table built correctly and came back empty. Relax the filters and the flow
returns without touching the overlay. See {doc}`06_edge_cases`.
:::

---

## Symptom: the notification reads `Velocity overlay unavailable`

### Likely causes (ordered)

1) The spawn-table build threw. This is the terminal failure state of that
   build, and it is the *only* thing this message means. The overlay also
   disables itself.  
2) The build allocates two GPU textures (one visibility, one spawn table) at the
   end; under memory pressure that allocation is where it fails.  

### How to confirm

- Look for the second notification, `Velocity overlay preparation failed: <message>`.
  That one carries the actual error text.
- Open the browser console; the thrown error is logged there too.
- Check whether the same view had just been enlarged, duplicated into more
  snapshots, or pushed to a high `Particle density:`.

### Fix

1) Reduce GPU load: lower `Particle density:`, set `Bloom strength:` to `0.00`,
   close extra snapshot views, shrink the window.
2) Toggle `Show overlay` off and on to retry the build.
3) If the console reports WebGL context loss, reload the page.

### Prevention

- Keep the overlay's settings modest while you also have many views open; the
  trail buffers are per view.

---

## Symptom: “It works in 2D but not in 3D” (or vice versa)

### Likely causes (ordered)

1) You only exported/provided the vector field for one dimension.  
2) The 3D vector field exists but has the wrong shape (e.g., `(n_cells, 2)` stored under a 3D key).  

### How to confirm

- Switch dimension and see whether the dropdown changes.
- Check your export folder or AnnData `obsm` keys for `*_2d` and `*_3d` variants.

### Fix

- Provide the missing dimension’s vector field (and reload).

### Prevention

- Always export the vector field for every dimension you expect users to use (at least 2D and/or 3D).

---

## Symptom: The flow looks “wrong” (directionality contradicts expectation)

### Likely causes (ordered)

1) Row order mismatch between embedding and vector field.  
2) Vectors were computed in a different basis than the embedding you’re viewing.  
3) Aggressive tuning (high turbulence, extreme bloom/contrast) is hiding structure.  
4) The underlying velocity/drift model is noisy/uncertain.  

### How to confirm

- Set `Turbulence:` low (0–0.2) and disable bloom (Bloom strength = 0).
- Switch to a different field (if you have one) and compare.
- In Python, verify:
  - same cell order
  - same embedding basis
  - reasonable magnitude distribution

### Fix

1) Verify data alignment (cell order, basis) and regenerate vectors if needed.
2) Use conservative visualization settings.
3) Cross-check with an independent velocity visualization (e.g., scVelo stream plot) to validate interpretation.

### Prevention

- Bake a “sanity check” step into your data prep pipeline (e.g., plot a few arrow glyphs in Python to confirm directionality before exporting).

---

## How do I add vector fields to my dataset?

This is the most common root cause when the overlay is missing entirely, and it
is a data-preparation question rather than a UI one. The full recipe — both the
`AnnData.obsm` route and the `cellucid.prepare(..., vector_fields=...)` route,
with the exact shape and dtype requirements — is
{doc}`../../python_package/c_data_preparation_api/08_vector_fields_velocity_displacement`.

Two rules from it decide whether the field is detected at all, so they are worth
repeating here:

- **The key carries the dimension.** Cellucid looks for `obsm` keys of the form
  `<something>_umap_<dim>d` — `velocity_umap_2d`, `T_fwd_umap_3d`. A key without
  the suffix is not a vector field, and a key whose array width disagrees with
  its suffix is rejected at load time with an error naming both numbers.
- **The matching embedding must exist.** A 2D field needs `X_umap_2d`, a 3D field
  needs `X_umap_3d`. Export every dimension you expect people to switch to.

A finished export folder carries the evidence: `dataset_identity.json` gains a
`vector_fields` block, and a `vectors/` directory appears holding one `.bin`
(or `.bin.gz`) per field and dimension.

---

## Common error messages (exact strings)

These are the messages the app actually emits, verbatim. Match on the wording,
not on a paraphrase — each one has exactly one cause.

| Message | Raised by | The one condition |
|---|---|---|
| `Dataset metadata does not advertise an exact 2D vector field key for "<id>"` (or `3D`) | AnnData adapter | The overlay asked for a dimension this field never registered a key for. |
| `Metadata advertises '<obsmKey>' for 2D vector field "<id>", but it is missing from obsm` | AnnData adapter | The key was registered at load time but is no longer present in `obsm`. |
| `Vector field '<obsmKey>' has N dimensions, expected D` | AnnData adapter | A 2D array is stored under a 3D key, or the reverse. Note the **single** quotes around the `obsm` key. |
| `Vector field '<obsmKey>' has N columns, but its Dd suffix requires exactly D` | AnnData adapter, at dataset load | The same shape mismatch, caught during discovery — so the field never appears in `Vector field:` at all. |
| `VectorFieldManager.loadField: vectors length X does not equal expected length Y.` | Vector field manager | The payload length is not `n_cells × components`: usually a stale or mismatched `vectors/*.bin` in an export folder. |
| `VectorFieldManager.loadField: field "<id>" contains a non-finite value at cell N, component C.` | Vector field manager | A `NaN` or `Inf` in the field. The message names the exact cell and component; see {doc}`06_edge_cases`. |
| `VectorFieldManager.loadField: scaling field "<id>" overflowed at cell N, component C.` | Vector field manager | The values are finite but so large that scaling them leaves `float32` range. |

Whichever one fires, the visible consequence is the same: `Show overlay` turns
itself back off and the info line reads `Failed to load vector field.`

---

## Debug bundle to collect (for bug reports)

When reporting an overlay bug, include:

- Dataset loading path: file picker / server / Jupyter / GitHub-hosted exports
- Browser + OS + GPU (e.g., Chrome 121, macOS 14.2, Apple M2)
- Dataset size: `n_cells`, and which dimension you were viewing (2D vs 3D)
- Selected `Vector field:` label and a screenshot of the overlay settings panel
- Any toast/notification error text (copy/paste)
- Browser console logs around the failure (copy/paste)
