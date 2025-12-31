# Edge cases

**Audience:** computational users and power users  
**Time:** 5–15 minutes  
**Goal:** prevent surprises by explicitly listing known edge cases and what Cellucid does.

---

## Data edge cases

Each entry below follows: **what you see → why → how to confirm → what to do**.

### 1) No vector fields at all

- **What you see:** the **Vector Field Overlay** block is missing from the sidebar.
- **Why:** your dataset contains no vector field metadata / no detectable vector arrays.
- **Confirm:** no “Vector Field Overlay:” section appears under **Visualization**.
- **What to do:** add/export vector fields (usually in Python) and reload the dataset.

### 2) Vector fields exist, but not for the current dimension

- **What you see:** the overlay block exists, but `Show overlay` is disabled and you see a hint like “Vector fields available for 2D, 3D. Switch embedding dimension to enable.”
- **Why:** vector fields are dimension-specific; the dataset might only have 2D vectors but you are viewing 3D.
- **Confirm:** switch the view dimension and watch the toggle enable/disable.
- **What to do:** switch to a supported dimension, or export the missing dimension’s vector field.

### 3) You have a vector field, but it’s in a different basis than Cellucid expects

- **What you see:** no fields appear (same as “no vector fields”), even though you computed velocity/drift.
- **Why:** the current detection rules are UMAP-oriented (e.g., `velocity_umap_2d`). If your vectors are in a different basis (e.g., t-SNE), Cellucid won’t detect them automatically.
- **Confirm:** check your AnnData `obsm` keys (or export metadata) for `*_umap*` field ids.
- **What to do:** store/export vectors using Cellucid’s naming convention for the basis currently supported by your workflow (UMAP).

### 4) All-zero or extremely small vectors

- **What you see:** overlay appears enabled, but the flow looks static (or almost static).
- **Why:** magnitudes are ~0 (either biologically, or due to a preprocessing mistake such as scaling down too aggressively).
- **Confirm:** in Python, inspect `np.linalg.norm(vectors, axis=1)`; in the UI, try increasing `Flow speed:` and `Particle size:`—if nothing changes, data may be near-zero.
- **What to do:** verify how the vector field was computed and whether it matches the embedding; ensure you’re not accidentally using an all-zero layer.

### 5) NaN/Inf in vectors

- **What you see:** flow has “dead zones”, discontinuities, or looks weaker than expected.
- **Why:** non-finite values are sanitized to 0 for stability.
- **Confirm:** check for `np.isfinite(vectors).all()` in Python.
- **What to do:** clean the vector field (replace invalid values) or recompute upstream.

### 6) Magnitudes are huge (units/normalization mismatch)

- **What you see:** particles shoot off rapidly; trails look like streaks unrelated to local structure.
- **Why:** the vector field was computed in a coordinate system that doesn’t match the embedding scale (or magnitudes were not normalized).
- **Confirm:** reduce `Flow speed:` to ~0.5×. If it still looks explosive, the raw vectors are likely too large.
- **What to do:** verify that vectors are computed in the same embedding coordinate system as the points (and that row order matches).

### 7) Row order mismatch (the “hardest” failure mode)

- **What you see:** overlay renders, but directionality looks nonsensical or contradicts known biology.
- **Why:** the vector field rows do not correspond to the same cells as the embedding rows.
- **Confirm:** compare a few cells by index between your embedding and vector field source; ensure no reindexing/shuffling happened between them.
- **What to do:** regenerate/export the vector field with consistent ordering. There is no UI-side fix.

---

## Scale edge cases

### 1) Huge datasets (hundreds of thousands to millions of cells)

- **What you see:** enabling the overlay is slow, or it becomes the main bottleneck.
- **Why:** even though vectors load lazily, the overlay still needs to simulate and render many particles every frame, and build spawn tables from visible cells.
- **Confirm:** performance improves immediately when you disable `Show overlay`.
- **What to do:** use the laptop-safe preset (`Particle density:` 5K–10K, short trails, bloom=0), and consider enabling LOD + `Sync with LOD`.

### 2) Many views/snapshots visible

- **What you see:** overlay performance is fine in single view but collapses when comparing many views.
- **Why:** each view needs its own trail buffers (and potentially bloom buffers). GPU memory and full-screen passes multiply with view count.
- **Confirm:** switch to a single view; if performance returns, multiview is the multiplier.
- **What to do:** reduce number of views or disable bloom, then increase density only if you have headroom.

### 3) Large window + high-DPI (retina) screens

- **What you see:** the overlay is expensive even at moderate particle counts.
- **Why:** trail buffers scale with the number of pixels rendered.
- **Confirm:** shrink the browser window; if FPS improves, resolution is the bottleneck.
- **What to do:** use a smaller window for exploration; save “big window cinematic” for final screenshots with conservative settings.

---

## UI/state edge cases

### 1) Switching datasets while the overlay is enabled

- **What you see:** after loading a new dataset, the overlay turns off.
- **Why:** Cellucid resets the overlay on dataset changes to avoid carrying state across datasets.
- **Confirm:** the `Show overlay` checkbox is unchecked after a dataset swap.
- **What to do:** re-enable the overlay for the new dataset.

### 2) Switching dimensions while the overlay is enabled

- **What you see:** the overlay may briefly show “Loading vector field…” or auto-disable.
- **Why:** the overlay must load the vector field for the new dimension (if available).
- **Confirm:** the dropdown list changes between 2D and 3D (fields are filtered by dimension).
- **What to do:** wait for the reload; if it disables, switch back to a dimension that has a field.

### 3) Render mode = “Volumetric smoke cloud”

- **What you see:** you can’t find the overlay controls.
- **Why:** the overlay controls currently live in the Points rendering controls area.
- **Confirm:** **Visualization → Render mode:** is not `Points`.
- **What to do:** switch to `Points` render mode.

### 4) “All cells filtered out”

- **What you see:** overlay is enabled but nothing draws (and you may see a “Velocity overlay unavailable” toast).
- **Why:** particles spawn from visible cells only; if no cells are visible, there is nothing to spawn from.
- **Confirm:** your filtered/visible cell count is 0 (or the plot is empty).
- **What to do:** relax filters (or disable overlay until cells are visible again).

---

## Environment edge cases

### 1) WebGL/GPU issues (context loss)

- **What you see:** the canvas turns blank, or the overlay never renders again after a heavy operation.
- **Why:** GPU memory pressure can cause “WebGL context lost” events.
- **Confirm:** browser console shows WebGL context loss messages.
- **What to do:** reload the page; reduce overlay load (density, trails, bloom, fewer views).

### 2) Corporate/locked-down browsers

- **What you see:** some loading operations fail, or the app behaves inconsistently across machines.
- **Why:** browser policies, extensions, or GPU driver constraints can interfere.
- **Confirm:** try the same dataset on a different browser/machine.
- **What to do:** use a supported browser configuration and keep GPU drivers up to date; if needed, use server mode and test in a clean profile.
