# Performance and quality

**Audience:** everyone (especially laptop users)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- What makes the overlay slow (what actually scales)
- How to tune quality vs performance intentionally
- How to distinguish “my GPU is overloaded” from “my data is wrong”

---

## Performance mental model (what actually scales)

The Vector Field Overlay is GPU-intensive for two independent reasons:

1) **Particles** (simulation + drawing)  
2) **Trails + post-processing** (full-screen passes, bloom, compositing)

Understanding which part is limiting you makes performance tuning much faster.

### 1) Particle cost scales with particle count

The biggest “knob” is:

- `Particle density:` (15K means 15,000 particles)

More particles:

- look smoother and “more stream-like”
- cost more GPU every frame (simulation + draw calls)

### 2) Trail/post-processing cost scales with pixels (resolution)

The overlay draws trails into a framebuffer roughly the size of your current view.

This means performance is strongly affected by:

- your browser window size
- device pixel ratio (retina screens push more pixels)
- how many views/snapshots are visible at once (each view needs its own trail buffers)

### 3) Bloom and other “cinematic” effects can dominate GPU time

Bloom involves multiple blur passes. Even with the default “subtle bloom”, it can become expensive when:

- the window is large
- you have many views open
- particle density and trail persistence are high (more bright pixels → more bloom work)

If you need speed, set `Bloom strength:` to `0.00` first.

### 4) `Sync with LOD` only helps if LOD is enabled

`Sync with LOD` reduces particle work when zoomed out, but it relies on the renderer’s Level-of-Detail system.

To benefit:

1) Enable **Visualization → Renderer settings → Level-of-Detail (LOD)**, and keep
2) **Vector Field Overlay → Sync with LOD** enabled.

---

## Recommended presets (copy/paste mental defaults)

These are intentionally conservative. Adjust one knob at a time.

### Laptop-safe preset

- `Particle density:` 5K–10K
- `Flow speed:` 2.0×–3.0×
- `Trail length:` 2.0s–5.0s
- `Opacity:` 40%–60%
- Advanced:
  - `Bloom strength:` 0.00–0.05 (consider 0.00)
  - `Trail persistence:` 0.900–0.930
  - `Turbulence:` 0.00–0.20 (lower is more “honest”)

### Desktop balanced preset

- `Particle density:` 15K–50K
- `Flow speed:` 2.0×–4.0×
- `Trail length:` 5.0s–10.0s
- `Opacity:` 50%–70%
- Advanced:
  - `Bloom strength:` 0.05–0.10
  - `Trail persistence:` 0.920–0.950
  - `Turbulence:` 0.10–0.30

### Presentation / cinematic preset (use with scientific caution)

- `Particle density:` 50K–200K
- `Trail length:` 8.0s–15.0s
- Advanced:
  - `Comet stretch:` 0.8–1.2
  - `Trail persistence:` 0.950–0.990
  - modest bloom (avoid going “white hot”)

:::{important}
If you are preparing publication figures, treat cinematic settings as a separate “look” preset and keep a conservative scientific preset you can return to.
:::

---

## Performance triage workflow (fast, reliable)

When the overlay is slow, do this in order (each step is a large win):

1) Reduce `Particle density:` to **5K**.
2) Reduce `Trail length:` to **2.0s–3.0s**.
3) Expand **Advanced Visual Settings → HDR & Bloom** and set `Bloom strength:` to **0.00**.
4) If you have many snapshot views: switch to fewer views (or clear snapshots) and retry.
5) Shrink the browser window (fewer pixels = faster) and retry.
6) Enable renderer LOD (**Visualization → Renderer settings → Level-of-Detail (LOD)**) and keep `Sync with LOD` enabled.
7) If you previously changed many filters quickly: disable and re-enable the overlay (it rebuilds internal spawn tables).

---

## Quality pitfalls (it looks wrong but is “working”)

### Too faint / “I can barely see it”

Common causes:

- `Opacity:` too low (or 0%)
- `Particle size:` too small
- `Intensity:` too low, or `Exposure:` too low (advanced)
- dark background + low contrast

Fixes:

1) Raise `Opacity:` to 60%.
2) Raise `Particle size:` to ~2.0.
3) If needed: raise `Intensity:` slightly before raising bloom/exposure.

### Too bright / “it looks like glowing fog”

Common causes:

- high `Particle density:` + long `Trail length:`
- high `Trail persistence:` (advanced)
- `Bloom strength:` too high (advanced)
- too-high `Opacity:`

Fixes (fastest first):

1) Reduce `Particle density:`.
2) Reduce `Trail length:`.
3) Reduce `Bloom strength:` or set it to 0.00.
4) Reduce `Opacity:`.

### Jittery / flickery / hard to interpret

Common causes:

- low FPS (GPU saturated)
- `Flow speed:` too high (particles leap between frames)
- very short trails (not enough temporal integration)

Fixes:

1) Reduce `Particle density:` until FPS recovers.
2) Reduce `Flow speed:` to ~2.0×.
3) Increase `Trail length:` slightly (but keep density moderate).

### Trails smear when you move the camera

Some smear is inevitable because trails integrate over time. Cellucid also compensates by fading trails faster during large camera moves, which can make trails look “shorter” while you’re navigating.

If smearing is still too strong:

- reduce `Trail length:` and/or reduce `Trail persistence:`

---

## Debugging checklist (fast)

Use this when you want to quickly decide whether you have:
“wrong settings”, “missing data”, or “GPU/browser issues”.

1) Confirm **Visualization → Render mode:** is `Points`.
2) Confirm the overlay section exists (**Vector Field Overlay:** is visible).
3) Confirm `Show overlay` is enabled (checkbox stays on).
4) Confirm you are on a dimension that has a vector field (if the toggle is disabled, switch dimension).
5) Set conservative settings:
   - `Particle density:` 5K
   - `Trail length:` 3.0s
   - `Opacity:` 60%
   - `Bloom strength:` 0.00
6) If it is still slow: reduce window size and reduce number of views.
7) If you see errors/toasts: go to `07_troubleshooting_velocity_overlay` and match the symptom.
8) If the entire canvas goes blank or you see “WebGL context lost”: reload the page and consider lowering GPU load (density, trails, bloom, fewer views).
