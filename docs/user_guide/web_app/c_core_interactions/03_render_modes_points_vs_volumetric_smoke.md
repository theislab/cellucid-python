# Render modes: Points or Volumetric smoke

**Audience:** everyone
**Time:** 10 minutes

Cellucid starts in **Points**. Choose **Volumetric smoke cloud** only when a
density-shaped overview is more useful than individual-cell inspection.

## Choose the representation before tuning it

| Question | Choose **Points** | Choose **Volumetric smoke cloud** |
|---|---|---|
| Do individual cells, rare populations, exact selections, or counts matter? | Yes | No |
| Do you need kept views for a side-by-side comparison? | Yes | No; smoke is single-view only |
| Is the goal to read a dense manifold's broad occupied shape? | Possible | Often clearer |
| Is the result intended as quantitative density evidence? | Use cells and an explicit analysis | No; smoke appearance depends on rendering controls |

A useful biological example is a densely sampled developmental atlas. Smoke can
make the occupied shape of a broad differentiation manifold easier to see before
you return to Points to inspect progenitors and terminal populations. It does
**not** establish that intermediate cells exist, estimate a probability density,
or turn an embedding into a physical tissue volume.

## Open the relevant controls

Use **Visualization → Render mode:**:

- `Points` is the default.
- `Volumetric smoke cloud` reveals **Volumetric smoke:** and hides the
  points-only depth, renderer, connectivity, and vector-field controls.

The small `i` button beside **Volumetric smoke:** explains the mode without
keeping extra prose open in the sidebar.

## Points mode

Points renders one sprite per displayed cell and is the working mode for
coloring, filtering, highlighting, connectivity, vector fields, and multiview.
Its relevant controls are:

- **Depth perception:** `Point size (log):`, `3D lighting:`,
  `Atmospheric fog:`, and `Perspective size scaling:`
- **Shader quality:** `Full (lighting + fog)`, `Light (circular, no lighting)`,
  or `Ultra-light (square points)`
- **Renderer settings:** `Level-of-Detail (LOD)`, `Force LOD level:`, and
  `Frustum culling`

For large datasets, begin here and use
{doc}`../n_benchmarking_performance/03_large_dataset_best_practices`.

## Smoke controls: exact current inventory

Every smoke slider has a 0–100 control position. The value printed beside it is
the meaningful render value.

### Quality & performance

| UI label | Initial position/readout | Meaning |
|---|---|---|
| `Grid density:` | 80 / `128³` | Resolution of the density volume. The exact available cube sizes are 32³, 48³, 64³, 96³, and 128³. The 128³ ceiling is valid even on WebGL2 devices that provide only the required minimum 2D texture size; this is the largest memory lever. |
| `Ray quality:` | 75 / `High` | Ray-marching work. The displayed bands are `Fast`, `Balanced`, `High`, and `Ultra`. |
| `Render resolution:` | 15 / `0.51x` | Off-screen render scale from 0.25x through 2.00x. Lowering it is usually the fastest recovery from poor frame rate. |
| `Noise detail:` | 58 / `128³` | Noise-texture resolution. The exact choices are 32³, 48³, 64³, 96³, 128³, 192³, and 256³. |

### Shape & motion

| UI label | Initial position | Meaning |
|---|---:|---|
| `Cloud density:` | 56 | Volume opacity. Its displayed value is recalculated when grid or noise resolution changes. |
| `Fine detail:` | 60 | High-frequency structure. Its displayed value is also resolution-adaptive. |
| `Turbulence:` | 10 (`10%`) | Procedural warping; it changes appearance, not the data. |
| `Animation speed:` | 40 (`1.00`) | Procedural smoke motion. Set it to 0 when motion distracts from structure. |
| `Edge softness:` | 0 (`0.2`) | Softness of the density boundary. |

### Lighting

| UI label | Initial position | Meaning |
|---|---:|---|
| `Light absorption:` | 65 | Attenuation through the cloud; the displayed value is resolution-adaptive. |
| `Light scattering:` | 0 (`0.0`) | Added scattered glow. |
| `Direct lighting:` | 3 (`0.04x`) | Direct-light contribution. |

`Cloud density:`, `Fine detail:`, `Light absorption:`, `Light scattering:`, and
ray work are intentionally adjusted with grid/noise resolution. Record the
printed readouts—not only slider positions—when reproducing an appearance.

## What is rebuilt, and when

Smoke is lazy: no density volume is built during ordinary Points startup. The
first switch to smoke builds a GPU volume from the **currently visible** cells.
Filtering, changing the active field, or changing view-scoped visibility marks
that volume dirty. While smoke is active, a completed visibility change queues
one rebuild after the newly published counts and controls have had a paint
opportunity; other changes from the same publication are coalesced. Only direct
`Grid density:` slider movement uses the approximately 300 ms debounce, so
dragging the slider does not rebuild on every input event.

The build stays on the GPU: Cellucid scans the existing position and visibility
arrays without making a full-size position copy, streams visible positions
through one bounded 3 MiB staging batch, splats them into a density atlas,
reduces its maximum there, normalizes to 8-bit density, and publishes the native
3D texture as one ordered transaction. No pixel data or normalization crosses
into JavaScript. The previous volume remains published until the replacement
commands, validation, and notification publication succeed. If a replacement
fails, Cellucid keeps that volume, reports the exact failure, and returns the
visible mode controls to Points so stale smoke is never presented as current.

This has three practical consequences:

1. A smoke cloud after filtering describes only the visible subset.
2. If no cells are visible, Cellucid clears the smoke volume instead of showing
   stale density.
3. Returning to Points is the reliable way to verify which cells produced a
   surprising shape.

## Multiview and session behavior

Smoke and kept views are mutually exclusive:

- **Keep view** is disabled in smoke.
- If kept views already exist, selecting smoke is rolled back and Cellucid
  reports `Volumetric smoke requires a single view. Clear snapshots first.`
- Smoke forces the live, single-view layout.

A `.cellucid-session` saves `Render mode:` and the smoke slider/select values.
The session does not contain the dataset or a frozen GPU volume; loading it with
the matching dataset rebuilds smoke from that dataset's restored visible cells.
See {doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

## Failure and recovery

### The cloud is blank

1. Return to `Points`.
2. Confirm that points are visible and **Showing … points** is not zero.
3. Relax categorical, continuous, and outlier filters.
4. Select `Volumetric smoke cloud` again.
5. Increase `Cloud density:` only after visible cells are confirmed.

### Interaction is slow or the WebGL context is lost

Recover in this order:

1. Lower `Render resolution:`.
2. Lower `Ray quality:`.
3. Lower `Grid density:`.
4. Lower `Noise detail:`.
5. Return to `Points`; after a context loss, reload the page before retrying.

Changing lighting first rarely addresses the dominant cost. For platform and
GPU checks, see {doc}`../a_orientation/02_system_requirements`; for general
interaction failures, see
{doc}`06_troubleshooting_core_interactions`.

## Related workflows

- {doc}`04_view_layout_live_snapshots_small_multiples` — points-only multiview
- {doc}`../i_vector_field_velocity/index` — points-only directional overlay
- {doc}`../e_filtering/01_filtering_mental_model` — the visible subset used to
  build smoke
- {doc}`../k_figure_export/index` — explicit, reproducible figure output
