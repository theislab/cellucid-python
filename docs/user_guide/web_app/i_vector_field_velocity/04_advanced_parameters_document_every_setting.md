# Advanced parameters (document every setting)

**Audience:** power users and figure-makers  
**Time:** 15–30 minutes  
**What you’ll learn:**
- What each advanced setting does (and which ones are purely cosmetic)
- How to produce readable overlays for talks/papers (without washing out the data)
- How to avoid “looks cool but misleading” parameter combinations

---

## Where these controls live

All advanced settings are under:

- **Visualization → Vector Field Overlay → Advanced Visual Settings** (expand the dropdown triangle)

They only show up after you enable `Show overlay`.

:::{important}
Advanced settings are powerful, but they can also make the overlay scientifically misleading (e.g., adding turbulence or heavy post-processing can suggest structure that isn’t in the data). Use them intentionally.
:::

---

## Particle Rendering (exact UI labels)

These controls change how individual particles look.

| Setting (UI label) | Range (UI) | Default | What it does | Performance impact | Guidance |
|---|---:|---:|---|---|---|
| `Intensity:` | 0.05–1.5 | 0.25 | Overall brightness/visibility of particles and trails. | Low–Medium | Increase if the overlay is too faint; decrease if it washes out the point cloud. |
| `Glow amount:` | 0–1 | 0.30 | Adds a halo around particles (“neon” look). | Medium | Keep subtle for papers; glow can hide fine structure. |
| `Comet stretch:` | 0–2 | 0.60 | Elongates particles in the direction of motion (makes direction easier to read). | Low | Good for presentations; moderate values are usually safe. |
| `Core sharpness:` | 0–1 | 0.70 | Controls how crisp the particle core is (sharp vs soft). | Low | Lower sharpness can feel “smokier”; higher can look noisy at high density. |

---

## Trail Settings (exact UI labels)

These controls change how motion accumulates over time.

| Setting (UI label) | Range (UI) | Default | What it does | Performance impact | Guidance |
|---|---:|---:|---|---|---|
| `Trail persistence:` | 0.900–0.995 | 0.925 | How strongly the previous frame remains visible (higher = longer afterimage). | Medium | This is the “looks like fog” knob. Increase carefully. |
| `Chromatic fade:` | 0–1 | 0.00 | Makes trails fade with different rates per color channel (rainbow trails). | Low | Leave at 0 for scientific figures; use only for artistic demos. |
| `Turbulence:` | 0–1 | 0.30 | Adds organic noise to particle motion. | Low–Medium | For faithful representation, set low (0–0.2). High turbulence can obscure true directionality. |

### “Trail length” vs “Trail persistence” (don’t confuse these)

- `Trail length:` (core control) changes **how long a particle lives** before respawning.
- `Trail persistence:` (advanced) changes **how long pixels remain visible** as an afterimage.

If you want “short, readable” trails:

- reduce `Trail length:` and/or reduce `Trail persistence:`

If you want “long, cinematic” trails:

- increase `Trail length:` and increase `Trail persistence:` (but expect higher GPU cost).

---

## HDR & Bloom (exact UI labels)

These controls add film-like glow and tone mapping. They can dramatically change the look.

| Setting (UI label) | Range (UI) | Default | What it does | Performance impact | Guidance |
|---|---:|---:|---|---|---|
| `Exposure:` | 0.10–2.00 | 0.50 | Global brightness of the overlay’s HDR composite. | Medium | If everything is too dark, raise exposure slightly before raising opacity. |
| `Bloom strength:` | 0.00–0.50 | 0.08 | How strong the glow spill is around bright particles/trails. | High | For scientific use, keep low (0–0.10) or set to 0 to effectively disable bloom. |
| `Bloom threshold:` | 0.10–1.00 | 0.75 | How bright something must be before bloom applies. | Medium | Increase threshold to keep bloom only on the brightest pixels. |
| `Anamorphic ratio:` | 1.0–3.0 | 1.2 | Makes bloom stretched horizontally (“cinematic streaks”). | High | Keep close to 1.0 for scientific figures. |

---

## Color Grading (exact UI labels)

These controls change the final color/contrast of the overlay layer.

| Setting (UI label) | Range (UI) | Default | What it does | Performance impact | Guidance |
|---|---:|---:|---|---|---|
| `Saturation:` | 0.50–2.00 | 1.15 | Makes colors more/less vivid. | Low | Avoid extreme saturation in scientific figures (can imply stronger magnitude differences than real). |
| `Contrast:` | 0.50–2.00 | 1.05 | Expands or compresses intensity range. | Low | Too much contrast makes trails “clip” and hide subtle structure. |
| `Highlights:` | 0.50–1.50 | 0.85 | Scales the bright end of the tone mapping. | Low | Lower highlights if bloom/opacity cause washout. |
| `Shadows:` | 0.50–1.50 | 1.05 | Scales the dark end of the tone mapping. | Low | Increase slightly if trails disappear into a dark background. |

---

## Cinematic Effects (exact UI labels)

These are almost entirely aesthetic.

| Setting (UI label) | Range (UI) | Default | What it does | Performance impact | Guidance |
|---|---:|---:|---|---|---|
| `Vignette:` | 0–1 | 0.00 | Darkens image edges. | Low | Usually disable for scientific figures (can hide peripheral structure). |
| `Film grain:` | 0.000–0.100 | 0.000 | Adds noise texture over the overlay. | Low | Disable for scientific figures; it can look like measurement noise. |
| `Chromatic aberr.:` | 0–1 | 0.00 | Separates color channels near edges (lens effect). | Low | Disable for scientific figures; can reduce readability and create false color edges. |

---

## Suggested workflows

### If you want a “clean scientific overlay”

1) Keep `Turbulence:` low (0–0.2).
2) Keep `Bloom strength:` low (0–0.08) and `Anamorphic ratio:` close to 1.0.
3) Avoid chromatic effects (`Chromatic fade:`, `Chromatic aberr.:`, `Film grain:`).
4) Use shorter trails and moderate density (see `03_core_parameters_document_exact_ui_labels`).

### If you want a “presentation / cinematic overlay”

1) Increase `Comet stretch:` (direction readability).
2) Increase `Trail persistence:` a bit (but monitor washout).
3) Use modest bloom (avoid extreme).
4) Keep a separate “scientific export” preset so your visuals stay honest when you switch contexts.
