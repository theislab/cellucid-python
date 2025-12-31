# Render modes (Points vs Volumetric “Smoke”)

**Audience:** everyone (smoke mode is optional, but easy to misuse)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- What each render mode is for (and what it is not)
- Why smoke mode can be slow/blank (and how to fix it)
- The key quality/performance knobs (with recommended defaults)

---

## Where render mode lives

Render mode is controlled in the left sidebar:

- **Visualization** accordion → **Render mode** dropdown
  - **Points**
  - **Volumetric smoke cloud**

---

## Points mode (default)

Points mode renders each cell as a point sprite with GPU-accelerated shading.

### When to use points mode

- you need precise, cell-level interpretation (“these are discrete cells”),
- you want to use **multiview snapshots** (Keep view),
- you are doing most scientific work (coloring, filtering, selecting, comparing).

### Controls that matter (points mode)

In the **Visualization** accordion you’ll see:

- **Depth perception**
  - Point size (log)
  - 3D lighting
  - Atmospheric fog
  - Perspective size scaling
- **Shader quality**
  - Full (lighting + fog)
  - Light (circular, no lighting)
  - Ultra-light (square points)
- **Renderer settings**
  - Level-of-Detail (LOD)
  - Force LOD level
  - Frustum culling

These are documented more deeply in the performance section, but the key idea is:

- lower quality = less GPU work = smoother interaction.

---

## Smoke mode (volumetric density rendering)

Smoke mode renders the dataset as a **ray-marched volumetric cloud**. It is intentionally cinematic and can be helpful for “density intuition”, but it is not the default scientific mode.

### When smoke mode is useful

- presentations and demos,
- giving intuition about density/structure at a glance,
- a “global shape” view of a very dense embedding.

### When smoke mode is *not* a good idea

- when you need to interpret individual cells,
- when you want multiview comparisons (snapshots),
- on machines with limited GPU memory.

:::{important}
Smoke mode is disabled when snapshots exist.

If you click “Volumetric smoke cloud” while you have kept views, the UI will force you back to points mode. Clear snapshots first if you really want smoke.
:::

---

## Smoke controls (what each knob does)

Smoke controls appear under **Visualization → Volumetric smoke** only when render mode is smoke.

### Quality & performance (start here)

- **Grid density**: resolution of the 3D density grid (major GPU cost driver).
- **Ray quality**: number of ray-marching steps (major GPU cost driver).
- **Render resolution**: scales the render target resolution (big speed lever; lower is faster).
- **Noise detail**: resolution of the noise texture (affects fine structure; also impacts cost).

### Shape & motion

- **Cloud density**: overall density/opacity of the volume.
- **Fine detail**: amount of high-frequency detail (scaled adaptively with resolution).
- **Turbulence**: warping strength of the noise field.
- **Animation speed**: motion speed of the smoke.
- **Edge softness**: softness of density boundaries.

### Lighting

- **Light absorption**: how much light is absorbed through the volume.
- **Light scattering**: how “glowy” the volume becomes.
- **Direct lighting**: intensity of direct light contribution.

:::{tip}
For “it’s too slow” issues, the biggest wins are usually:

1) lower Grid density
2) lower Ray quality
3) lower Render resolution
:::

---

## Deep path (exact behavior that matters for experts)

Smoke sliders are not independent:

- Changing **Grid density** and **Noise detail** changes adaptive scaling factors.
- Several parameters (density/detail/absorption/scatter/steps) are automatically scaled based on resolution so the cloud stays visually comparable across settings.

Notable discrete settings (implementation detail):

- Grid density maps to a discrete set of cubic sizes (e.g., 32³ … 1024³).
- Noise detail maps to discrete noise resolutions (e.g., 32³ … 256³).

Smoke rebuild behavior:

- Smoke volume is built on first switch to smoke mode.
- It is rebuilt after certain changes using a short debounce (~300 ms), so rapid slider changes don’t rebuild on every frame.

---

## Screenshots (placeholders)

<!-- SCREENSHOT PLACEHOLDER
ID: render-mode-points-vs-smoke
Where it appears: Render modes → Points vs Smoke
Capture:
  - Same dataset, same camera framing
  - One screenshot in Points mode, one in Smoke mode (before/after pair)
Crop:
  - Include: enough sidebar to show the selected Render mode
Annotations:
  - Optional: label which is points vs smoke
Alt text:
  - Side-by-side comparison of points rendering and volumetric smoke rendering.
Caption:
  - Points mode shows individual cells; smoke mode shows a volumetric density cloud for cinematic/global-shape views.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for points vs smoke comparison.
:width: 100%

Points mode shows individual cells; smoke mode shows a volumetric density cloud for cinematic/global-shape views.
```

---

## Edge cases

- **Smoke looks blank**: density may be too low, or the dataset has few visible points (filters/outliers). Increase cloud density and/or clear filters.
- **Smoke is extremely slow**: grid density and ray quality are likely too high for your GPU.
- **Switching render mode “does nothing”**: if snapshots exist, smoke is blocked; clear snapshots first.

---

## Troubleshooting (render modes)

For a full catalog, see `c_core_interactions/06_troubleshooting_core_interactions`.

### Symptom: Smoke mode is blank

**Likely causes:**

- Cloud density is too low.
- All/most points are filtered out (nothing to render).
- GPU precision/driver issue (rare; try points mode to confirm).

**Fix:**

1) Switch to points mode and confirm points exist.
2) Clear filters/outlier thresholds.
3) Switch back to smoke and increase **Cloud density**.

### Symptom: Smoke mode crashes / context lost

Go to `a_orientation/02_system_requirements` → “WebGL context lost”.
