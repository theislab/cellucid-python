# Troubleshooting (velocity overlay)

This page is a troubleshooting-first companion to the vector field overlay docs.

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

1) You are in `Render mode: Volumetric smoke cloud` (overlay controls are currently in the Points-only controls area).  
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

- Use explicit keys: `<field>_umap_<dim>d` (e.g., `velocity_umap_2d`).
- Validate shapes before export: `(n_cells, dim)` with matching row order.

---

## Symptom: Overlay is enabled, but nothing is visible on the canvas

### Likely causes (ordered)

1) `Opacity:` is 0% or very low, or `Particle size:` is too small.  
2) All cells are filtered out / invisible, so there is nothing to spawn particles from.  
3) Vectors are all-zero (or extremely small), so motion is imperceptible.  
4) You are GPU-saturated and frames are not updating smoothly.  

### How to confirm

- Set `Opacity:` to 60% and `Particle size:` to ~2.0.
- Temporarily relax filters so you have visible points.
- Reduce `Particle density:` to 5K (to reduce GPU load).
- If possible in Python: check vector magnitudes (e.g., `np.linalg.norm(vectors, axis=1)`).

### Fix

1) Adjust `Opacity:` and `Particle size:`.
2) Ensure you have visible cells (undo filters that hide everything).
3) Increase `Flow speed:` slightly (2×–4×) to see motion.
4) If the data might be all-zero, recompute vectors or verify you’re exporting the intended array.

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
4) LOD is disabled, so `Sync with LOD` can’t reduce workload.  

### How to confirm

- Disable the overlay and see if FPS immediately recovers.
- Switch to single view (or clear snapshots) and compare.
- Shrink the browser window and compare.

### Fix

1) Set `Particle density:` to 5K–10K.
2) Set `Trail length:` to 2–5s.
3) Set `Bloom strength:` to 0.00 (Advanced → HDR & Bloom).
4) Reduce number of visible views.
5) Enable renderer LOD (**Visualization → Renderer settings → Level-of-Detail (LOD)**) and keep `Sync with LOD` enabled.

### Prevention

- For large datasets, document a “default performance-safe preset” and encourage users to start there.

---

## Symptom: “Velocity overlay unavailable” / “Preparing velocity overlay…” takes a long time

### Likely causes (ordered)

1) Your dataset is huge and Cellucid is building internal spawn tables (especially after filter changes).  
2) All cells are invisible (spawn table has no candidates).  
3) The browser is too busy to schedule the background work promptly.  

### How to confirm

- If you have 0 visible points, the overlay cannot spawn particles.
- If you are changing filters rapidly, the overlay may rebuild repeatedly.

### Fix

1) Ensure some cells are visible.
2) Stop rapidly changing filters; wait a moment for the overlay to finish preparing.
3) Reduce GPU load (lower density, shorter trails, bloom=0).
4) If needed, toggle the overlay off and on once after you settle on filters.

### Prevention

- On very large datasets, prefer LOD + `Sync with LOD`, and keep density modest.

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

## How do I add vector fields to my dataset? (common fixes)

This is the most common “root cause” when the overlay is missing.

### Option A: Add vectors to `AnnData.obsm` (UMAP-based naming)

Cellucid detects UMAP vector fields using keys like:

- `velocity_umap_2d` or `velocity_umap_3d`
- `T_fwd_umap_2d` / `T_bwd_umap_2d` (CellRank drift-style)

Example (2D):

```python
# adata.obsm["X_umap_2d"] should be (n_cells, 2)
# velocity_umap_2d must also be (n_cells, 2) in the same row order.
adata.obsm["velocity_umap_2d"] = velocity_umap_2d.astype("float32")
```

CellRank transition matrices → drift vectors helper:

```python
import cellucid

cellucid.add_transition_drift_to_obsm(
    adata,
    T_fwd,          # (n_cells, n_cells)
    basis="umap",
    field_prefix="T_fwd",
    dim=2,
)
```

Then load the dataset via Jupyter/server/browser and enable the overlay.

### Option B: Export vectors with `cellucid.prepare(..., vector_fields=...)`

If you are generating a pre-exported folder:

```python
from cellucid import prepare

prepare(
    latent_space=latent,
    obs=obs,
    var=var,
    gene_expression=X,
    X_umap_2d=X_umap_2d,
    vector_fields={
        "velocity_umap_2d": velocity_umap_2d,
        # Add more fields if you like:
        # "T_fwd_umap_2d": drift_umap_2d,
    },
    out_dir="exports/my_dataset",
    compression=6,
)
```

The resulting export folder will contain:

- `dataset_identity.json` with a `vector_fields` block
- `vectors/*.bin` (or `*.bin.gz`) files

---

## Common error messages (what they usually mean)

If you get a toast/notification error, these strings are helpful:

- `No 2D/3D vector field found for "…" in obsm`  
  The key exists in metadata, but no matching `obsm` array is found (or the wrong dimension).

- `Vector field "…" has N dimensions, expected D`  
  You provided a 2D array under a 3D key (or vice versa).

- `vectors length X !== expected Y`  
  The loaded binary vector file has the wrong length for `n_cells * components` (export mismatch / wrong file).

---

## Debug bundle to collect (for bug reports)

When reporting an overlay bug, include:

- Dataset loading path: file picker / server / Jupyter / GitHub-hosted exports
- Browser + OS + GPU (e.g., Chrome 121, macOS 14.2, Apple M2)
- Dataset size: `n_cells`, and which dimension you were viewing (2D vs 3D)
- Selected `Vector field:` label and a screenshot of the overlay settings panel
- Any toast/notification error text (copy/paste)
- Browser console logs around the failure (copy/paste)
