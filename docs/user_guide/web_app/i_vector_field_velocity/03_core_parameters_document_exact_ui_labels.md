# Core parameters (document exact UI labels)

**Audience:** everyone (with a “fast path” for beginners)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- The meaning of each core slider/toggle (using the **exact UI labels**)
- Which knobs change *appearance* vs *performance*
- Safe “first settings” for laptops vs desktops

---

## Core controls (exact UI labels, ranges, defaults)

These controls appear under **Visualization → Vector Field Overlay** after you enable `Show overlay`.

| Setting (UI label) | Range (UI) | Default | What it changes (user-facing) | Performance impact | Common footguns |
|---|---:|---:|---|---|---|
| `Particle density:` | 1–500 (shown as `K`) | 15K | How many particles are simulated/drawn at once. Higher = denser “stream” look. | High (scales ~linearly with particle count) | High density can look like fog (misleading) and can tank FPS. |
| `Flow speed:` | 0.05×–5.0× | 3.0× | How fast particles move per frame along the field. It is a **visual multiplier**, not biological time. | Medium | Too high can look chaotic (particles “teleport” across structure). Too low can look “stuck”. |
| `Trail length:` | 0.1s–15.0s | 8.0s | How long a particle lives before it respawns. Longer lifetimes usually yield longer visible trajectories. | Medium–High | Very long trails can overwhelm the point cloud and hide structure. |
| `Particle size:` | 0.5–30 | 1.0 | The size of each particle (screen-space point size, scaled by device pixel ratio). | Low–Medium | Too small can look like “nothing is happening”; too big can hide points. |
| `Opacity:` | 0%–100% | 60% | Overall transparency of the overlay. | Low | `0%` makes the overlay invisible. High opacity can wash out points. |
| `Color scheme:` | discrete list | `Viridis` | How particle color is mapped from **vector magnitude** (speed) to color. | Low | Users sometimes assume color encodes direction; by default it encodes magnitude. |
| `Sync with LOD` | on/off | on | Reduces particle work when the renderer is zoomed out / using higher LOD (and uses LOD sampling for spawn candidates). | Usually helpful | Turning this off can be very expensive on large datasets; turning it on can reduce fine detail when zoomed far out. |

### Notes on units (important)

- `Particle density:` is displayed in **thousands**.  
  Example: `15K` means **15,000 particles**. The maximum is `500K`.

- `Flow speed:` is a multiplier (×). It is intentionally unitless.

- `Trail length:` is shown in seconds (s), but it is a visualization lifetime, not biological time.

---

## Recommended starting presets

### Laptop-safe preset (start here if anything stutters)

- `Particle density:` 5K–10K
- `Flow speed:` 2.0×–3.0×
- `Trail length:` 2.0s–5.0s
- `Particle size:` 1.0–2.0
- `Opacity:` 40%–60%
- `Sync with LOD`: on

### Desktop balanced preset

- `Particle density:` 15K–50K
- `Flow speed:` 2.0×–4.0×
- `Trail length:` 5.0s–10.0s
- `Particle size:` 1.0–2.0
- `Opacity:` 50%–70%
- `Sync with LOD`: on

:::{important}
If you are making publication figures, prefer **clarity over drama**:
shorter trails, lower density, and minimal post-processing.
:::

---

## Interaction notes (how controls combine)

### Density × trail length = “visual accumulation”

If the overlay looks like a glowing haze, the most common cause is:

- too many particles *and* too-long trails

Fix by reducing `Particle density:` first (largest win), then reducing `Trail length:`.

### Opacity × bloom (advanced) can wash out the view

If you later enable/adjust bloom in advanced settings, remember:

- `Opacity:` makes everything brighter
- bloom amplifies the brightest parts even more

Use a conservative opacity (40%–60%) when bloom is non-zero.

---

## Mini troubleshooting (core controls)

### “I enabled the overlay but see nothing”

Try, in this order:

1) Confirm `Opacity:` is not 0%.
2) Increase `Particle size:` to ~2.0.
3) Reduce `Flow speed:` to ~2.0× (very high speed can look like flicker/noise).
4) Confirm you are in a dimension that actually has a vector field (see `02_enabling_overlay_and_selecting_field`).

### “It’s a bright haze / I can’t see points anymore”

1) Reduce `Particle density:` (biggest impact).
2) Reduce `Trail length:`.
3) Reduce `Opacity:`.

### “It’s very slow / stuttering”

1) Reduce `Particle density:` to 5K.
2) Reduce `Trail length:` to 2–3s.
3) Keep `Sync with LOD` on.
4) If you have many snapshot views, switch to a single view (multiview multiplies GPU work).
