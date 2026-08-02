# Core particle-flow controls

**Audience:** everyone
**Time:** 10 minutes

These controls appear under **Visualization → Vector Field Overlay:** after
`Show overlay` is enabled.

```{figure} ../../../_static/screenshots/vector_field_velocity/set-core-overlay-controls.png
:alt: The Vector Field Overlay controls with Show overlay checked, the Vector field selector showing velocity_umap, and the Particle density, Flow speed, Trail length, Particle size, Opacity, Color scheme and Sync with LOD controls at their initial values, above a collapsed Advanced Visual Settings group.
:width: 472px

Every core control, at the values the table below records as initial: `15K`,
`3.0×`, `8.0s`, size `1`, `60%`, `Viridis`, `Sync with LOD` on. The numeric
readout beside each slider is the *transformed* value, not the slider position.
```

## Exact labels, ranges, and initial values

| UI label | Exact range/options | Initial value | What changes |
|---|---|---:|---|
| `Vector field:` | Fields declared for the selected dimension | Dataset default or explicit selection | The per-cell vector array. Changing it loads a different field. |
| `Particle density:` | 1K–500K, 1K steps | 15K | Simulated/drawn particle count. This is the primary particle-cost control. |
| `Flow speed:` | 0.05×–5.00×, 0.01× steps | 3.0× | Visual advection multiplier; it is not biological time. |
| `Trail length:` | 0.1s–15.0s, 0.1s steps | 8.0s | Particle lifetime before respawn; it is not elapsed biological time. |
| `Particle size:` | 0.5–30, 0.5 steps | 1.0 | Screen-space particle size. Large particles increase pixel fill and can hide cells. |
| `Opacity:` | 0%–100%, 1% steps | 60% | Composite opacity. `0%` is fully invisible. |
| `Color scheme:` | `Viridis`, `Plasma`, `Turbo`, `Cividis`, `Magma`, `Coolwarm` | `Viridis` | Maps normalized vector magnitude to particle color; it does not encode direction. |
| `Sync with LOD` | off/on | on | At coarser LODs, reduces particle count and uses the LOD-visible cells as spawn candidates. |

The controls update the renderer while the overlay is active. They do not
rewrite `adata.obsm`, prepared vector binaries, or the dataset identity.

:::{warning}
`Sync with LOD` is a fidelity control, not only a performance one. It scales
particle count down in steps as the level of detail coarsens, and past a
threshold it **disposes the overlay entirely** — the flow vanishes with the
camera zoomed out and no message. If the animation disappears while you are
zooming out, that is this control, not a missing field. Turn it off, or zoom
back in.

It only has an effect when **Visualization → Renderer settings →
Level-of-Detail (LOD)** is itself enabled; with LOD off, this checkbox does
nothing.
:::

:::{note}
At most 65,536 distinct cells can ever seed particles, chosen deterministically
from the currently visible cells. On a dataset larger than that, the flow is
drawn from a sample of your cells, not from all of them. The sample is stable
for a given visible set, so the picture is reproducible — but a sparse-looking
region in a very large dataset may be a sampling artefact rather than a
low-magnitude one.
:::

## Scientific reading preset

Use this when the question is “does the supplied direction agree with this
annotated continuum?”:

- `Particle density:` 10K–20K
- `Flow speed:` 1.5×–3.0×
- `Trail length:` 2.0s–5.0s
- `Particle size:` 1.0–2.0
- `Opacity:` 40%–60%
- `Color scheme:` `Viridis` or `Cividis`
- `Sync with LOD`: on

Then expand **Advanced Visual Settings** and reduce `Turbulence:` toward 0.
This makes the supplied vectors easier to distinguish from procedural motion.
Use {doc}`04_advanced_parameters_document_every_setting` for every advanced
control.

## Performance-first preset

Begin here on an integrated GPU, high-DPI display, or multiview layout:

- `Particle density:` 5K
- `Trail length:` 2.0s
- `Particle size:` 1.0
- `Opacity:` 50%
- `Sync with LOD`: on

Raise density only after interaction is stable. Density changes simulation and
draw work approximately with particle count; speed, lifetime, opacity, and
palette are mainly interpretive/visual controls. Trail buffers and bloom also
scale with viewport pixels, so reducing particle density alone may not cure a
high-resolution bottleneck. See {doc}`05_performance_and_quality`.

## Read combined controls correctly

### Dense haze

High `Particle density:`, long `Trail length:`, high `Trail persistence:`, and
high opacity/bloom can merge separate directions into a bright sheet. Reduce
density first, then lifetime, persistence, and opacity.

### Flicker or apparently discontinuous motion

Reduce `Flow speed:` before increasing trails. A large visual multiplier can
move particles too far between frames even when the source vectors are valid.

### No visible flow

Check in this order:

1. `Opacity:` is above 0%.
2. The selected view has visible cells.
3. `Particle size:` is large enough for the display.
4. The status does not read `Failed to load vector field.`
5. Vector magnitudes are not all zero.

See {doc}`07_troubleshooting_velocity_overlay` for data, network, dimension,
and WebGL recovery.

## Multiview and session scope

Core overlay settings are global across visible views; each view has its own
dimension, visibility-derived spawn table, and trail buffers. Do not expect
`Keep view` to freeze a separate particle preset.

For a matching dataset, `.cellucid-session` saves the current field, toggle,
core controls, and advanced controls. The vector arrays are loaded again from
the dataset. See {doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

The exact initial controls are shown on the real Pancreas dataset in
{doc}`08_screenshots`.
