# Render modes: Points or Volumetric smoke

**Audience:** everyone
**Time:** 10 minutes

Cellucid starts in **Points**. Choose **Volumetric smoke cloud** only when a
density-shaped overview is more useful than individual-cell inspection. The
dropdown labels it `(beta)` and the sidebar tags it `Beta`: it is slower than
points, and figure export cannot reproduce it.

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

Use **Visualization → Render mode**:

- `Points` is the default.
- `Volumetric smoke cloud (beta)` reveals **Smoke quality**, **Smoke shape**
  and **Smoke lighting** and hides the points-only shader, point size, depth,
  level-of-detail, connectivity, and vector-field controls. `Antialiasing` (in
  **Performance**) stays, in both modes, because both renderers draw into the
  same application-owned multisampled scene target rather than into a buffer of
  their own.

A pill reading `Beta` sits beside the `Render mode` label and is shown only
while smoke is the selected mode. It names the maturity of the selected option,
not of the control, so `Points` never wears it, and its hover text reads
`Volumetric smoke cloud is in beta: it is slower than points, and figure export cannot reproduce it.`
Smoke is the only thing in Cellucid carrying a maturity tag; see
{doc}`../a_orientation/04_ui_glossary_terminology`.

The small `i` button beside **Smoke quality** explains the mode without
keeping extra prose open in the sidebar.

```{figure} ../../../_static/screenshots/web_app/render-mode-select.png
:alt: A control labelled RENDER MODE holding a dropdown that reads Points, ringed in blue with the mouse pointer pressing it, above a dashed box labelled DEPTH PERCEPTION whose first slider, POINT SIZE (LOG), reads 2.87.
:width: 480px

`Render mode` is the first control in **Visualization**, and it is the switch
that decides which of the controls beneath it exist at all. The depth block
under it, shown here (now titled **Depth cues**), is one of the blocks that
disappears in smoke mode.
```

## Points mode

Points renders one sprite per displayed cell and is the working mode for
coloring, filtering, highlighting, connectivity, vector fields, and multiview.
Its relevant controls are:

- `Point size`, then **Depth cues:** `3D lighting`,
  `Atmospheric fog`, and `Perspective size scaling`. `Point size` is a
  logarithmic slider, not a percentage, and a dataset opens at the size its cell
  count implies rather than at a fixed default — which is why the same slider
  position looks different on a 3,696-cell trajectory and on an atlas.
- **Shader quality** `Full (lighting + fog)`, `Sculpted (soft 3D light)`,
  `Light (circular, no lighting)`, or `Ultra-light (square points)` — a
  **look** choice: in 3D views `Full`, `Light` and `Ultra-light` measured the
  same frame time at every point size tested and `Sculpted` 3–5% above `Full`,
  and in 2D views all four measured within a few percent of one another
  (Apple M1 Pro). In 1D and
  2D views the round qualities draw a smooth one-pixel edge, which
  antialiasing cannot give a point sprite. Note that `Ultra-light` draws
  squares rather than round dots, so it also changes far more pixels when
  antialiasing is off.

  `Sculpted` shades every cell as a soft-lit, glossy ball: a wrapped key
  light, a sky-and-ground ambient, a broad gloss in place of a hard white
  highlight (no glitter on dense clouds, and none at all on dots too small to
  show it), and darkening toward each ball's edge, so touching cells stay
  separate. The lighting is computed in linear light and normalised so that
  the **centre of every ball is exactly its legend colour**. Its points are
  opaque: fog fades their colour toward the fog colour instead of thinning
  them, and with antialiasing on in a 3D view their edges are smoothed per
  sample against the depth of what is behind, so a soft rim never shows the
  background through a nearer cell. Highlighted cells get the same 3D ring as
  `Full`, held still instead of pulsing, so a selection costs nothing while the
  view is idle. `Sculpted` is not the default; `Full` is.
- **Performance:** `Level of detail (LOD)` draws a coarser subset of the
  cells when you are zoomed out; `Frustum culling` skips cells the camera is not
  pointing at. Both are checkboxes and **both start off**. Ticking
  `Level of detail (LOD)` reveals a `Force LOD level` slider directly beneath
  it — that row is hidden entirely while the checkbox is clear, which is the
  usual reason people cannot find it. The slider sits at its lowest position,
  reading `Auto`, until you move it. `Auto` follows the camera at every dataset
  size, down to a floor of `min(2,000,000, cells ÷ 8)` points; an explicitly
  forced level is not held to that floor and can go coarser.
- `Antialiasing`, also under **Performance** — a live setting, like
  every other: a change reaches the next frame. Unlike the two checkboxes above
  it, it is not hidden in smoke mode, because it governs the whole scene target,
  which both renderers share. Until you click
  the checkbox yourself it follows the dataset, and the first line of its (i)
  says so — `Chosen automatically: 3,696 cells in 3D.` on the Pancreas sample, or
  `Chosen automatically from the size of the dataset.` before a dataset is
  open. Clicking it ends automatic selection.

Every measured number behind those four bullets — the point-size curve and its
reachable range, when shader quality is and is not a speed lever, the two-million-point
LOD budget, and the dataset size at which antialiasing switches itself off —
lives in one place:
{doc}`../n_benchmarking_performance/02_performance_considerations_what_gets_slow_and_why`.
For why a `.cellucid-session` neither records nor changes the antialiasing
setting, see
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

For large datasets, begin here and use
{doc}`../n_benchmarking_performance/03_large_dataset_best_practices`.

## What the smoke cloud shows

Smoke draws the **currently visible** cells as a lit, semi-transparent cloud,
coloured by the active field: each region of the cloud takes the mean display
colour of the cells there, so a categorical field tints each cluster in its
legend colour and a continuous field blends its colormap. Density sets how
opaque a region is; the colour says which cells produced it.

```{figure} ../../../_static/screenshots/web_app/window-smoke-cloud.png
:alt: The whole window in Volumetric smoke cloud mode: the Pancreas cells drawn as a soft, lit cloud tinted with the colours of their clusters inside the grid box, and the sidebar showing Render mode set to Volumetric smoke cloud (beta) with its Beta tag, then the Volumetric smoke controls with Grid density 128³, Ray quality High, Render resolution 0.50x, Noise detail 128³, Cloud density 1.00 and Animation speed Still.
:width: 1440px

The Pancreas sample coloured by `clusters` and switched to smoke, with every
control at its initial value except `Animation speed`, which is at 0
(`Still`). The cloud keeps each cluster's legend colour; where two clusters
overlap, their colours mix.
```

A 2D or 1D embedding, or a 3D one that happens to be flat, has no thickness
along its missing axes. Smoke spreads such a cloud into a soft slab (a soft
tube for 1D) of fixed thickness rather than drawing a sheet a ray would step
straight over, so a 2D embedding renders as a thin, evenly lit layer.

## Smoke controls: exact current inventory

Every smoke slider has a 0–100 control position. The value printed beside it is
the render value, and each slider maps its position through one fixed scale: no
slider rescales another, so a readout means the same thing at every grid and
noise resolution, and a recorded readout reproduces an appearance.

### Smoke quality

| UI label | Initial position/readout | Meaning |
|---|---|---|
| `Grid density` | 80 / `128³` | Resolution of the density volume: 32³, 48³, 64³, 96³, or 128³, in five equal bands of the slider. 128³ is the maximum on every device: its float splat atlas is 1536×1408, inside WebGL2's guaranteed 2048 `MAX_TEXTURE_SIZE`, and the published volume is 8 MiB. Changing it rebuilds the volume. Because rays take a fixed number of samples per voxel, a finer grid also means proportionally more samples per ray. |
| `Ray quality` | 75 / `High` | Samples per voxel along each ray, 0.5 to 2 on a logarithmic scale (1.41 at the initial position). Frame cost is about proportional to it; empty space costs the same at every setting. The readout names the band: `Fast` (below 33), `Balanced` (33–65), `High` (66–89), `Ultra` (90 and above). |
| `Render resolution` | 33 / `0.50x` | Fraction of the canvas the ray march runs at, 0.25x to 1.00x in steps of 0.05; the result is upsampled into the scene. Whatever the setting, one frame marches at most 2560×1440 rays, so a 4K or high-DPI canvas at 1.00x is rendered at that size and upsampled. Lowering it is usually the fastest recovery from a poor frame rate. |
| `Noise detail` | 80 / `128³` | Resolution of the two procedural noise volumes: 32³, 48³, 64³, 96³, or 128³, in five equal bands. At 128³ they take about 18 MiB of GPU memory together. Changing it regenerates the noise, not the density volume. |

### Smoke shape

| UI label | Initial position/readout | Meaning |
|---|---|---|
| `Cloud density` | 65 / `1.00` | Opacity, logarithmic from 0.05 to 5, so equal slider movements make equally noticeable changes, from thin haze to solid cloud. |
| `Fine detail` | 40 / `2.0` | Strength and frequency of the high-frequency noise, 0 to 5. |
| `Turbulence` | 10 / `10%` | Curl-noise warping, 0% to 100%; it changes appearance, not the data. |
| `Animation speed` | 40 / `1.00` | Procedural smoke motion, 0 to 2.5. At 0 the readout is `Still`: the smoke stops moving, and once the view has settled it is not redrawn at all, so the GPU stays idle until the camera or a setting changes. Use it when motion distracts from structure. |
| `Edge softness` | 15 / `0.3` | Noise erosion of thin regions, 0 to 2; 0 turns it off. |

### Smoke lighting

| UI label | Initial position/readout | Meaning |
|---|---|---|
| `Light absorption` | 50 / `1.0` | Shadow depth, 0 (no shadows) to 2. |
| `Light scattering` | 50 / `1.0` | Forward-scattering glow towards the light, 0 to 2. |
| `Direct lighting` | 50 / `1.00x` | Key-light intensity, 0x to 2x. |

Only `Grid density` rebuilds the volume and `Noise detail` regenerates the
noise; every other control applies on the next frame.

## What is rebuilt, and when

Smoke is lazy: no density volume is built during ordinary Points startup. The
first switch to smoke builds a GPU volume from the **currently visible** cells
and their current colours. The volume is a function of three things — which
cells are visible, where they are, and what colour they are drawn in — so each
of these marks it dirty:

- a visibility change: filtering, outlier thresholds, view-scoped visibility;
- a colour change: switching the active field, a colormap, a category colour,
  or the log scale;
- a position change: switching dimension (3D, 2D, 1D), or a new dataset.

While smoke is active, a dirty volume is rebuilt once the newly published
counts and controls have had a paint opportunity; changes from the same
publication are coalesced. While Points is active, the volume only stays dirty
and the next switch to smoke rebuilds it. Only direct `Grid density` slider
movement uses the approximately 300 ms debounce, so dragging the slider does
not rebuild on every input event.

The build stays on the GPU. Cellucid scans the existing position, colour, and
visibility arrays without making a full-size copy, streams visible cells
through one bounded staging batch (positions plus 8-bit colours, at most
262,144 cells), and splats each batch twice: density into a 32-bit float atlas
and colour into a four-channel float atlas. It reduces the density maximum on
the GPU and normalizes both into one 8-bit RGBA volume — red is density, the
other three channels the mean colour premultiplied by that density — and
publishes the native 3D texture as one ordered transaction. From that volume it
then derives a half-resolution lighting volume (the depth of cloud towards the
light, and the sky occlusion) and a small occupancy volume that lets each ray
leap over empty space. No pixel data or normalization crosses into JavaScript.
The previous volume remains published until the replacement commands,
validation, and notification publication succeed. If a replacement fails,
Cellucid keeps that volume, reports the exact failure, and returns the visible
mode controls to Points so stale smoke is never presented as current.

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

A `.cellucid-session` saves `Render mode` and the smoke slider/select values.
It stores slider positions, and because every scale is fixed, the same
positions restore the same readouts. The session does not contain the dataset
or a frozen GPU volume; loading it with the matching dataset rebuilds smoke from
that dataset's restored visible cells and colours. The built-in samples' saved
views carry the initial values in the tables above. A session saved in smoke and
opened on a device that cannot run smoke restores everything else and opens in
`Points`, with the message below. See
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

## Failure and recovery

### Smoke is not available on this device

Smoke builds its cloud by adding every cell into high-precision (float32)
buffers on the GPU, which needs float blending. Safari on iPhone and iPad, and
many Android browsers, do not offer it. Cellucid checks this before building
anything: choosing smoke there shows a notification titled **Volumetric smoke
unavailable** — “Volumetric smoke isn't available on this device: its graphics
driver can't blend high-precision colour, which the smoke cloud needs to add up
the cells. Points rendering works as usual.” — and `Render mode` stays on
`Points`. Other missing capabilities get the same headline with their own
reason. Use a desktop browser for smoke.

### The cloud is blank

1. Return to `Points`.
2. Confirm that points are visible and **Showing … points** is not zero.
3. Relax categorical, continuous, and outlier filters.
4. Select `Volumetric smoke cloud (beta)` again.
5. Increase `Cloud density` only after visible cells are confirmed.

### Interaction is slow or the WebGL context is lost

A smoke frame costs about rays × samples per ray: rays are set by the window
and `Render resolution` (never more than 2560×1440 of them), and samples per
ray by `Ray quality` and, through the voxel size, `Grid density`. A frame that
runs past the GPU's watchdog — about two seconds on Windows and macOS — loses
the whole WebGL context. Memory is rarely the cause: the density grid is capped
at 128³, an 8 MiB volume (under 50 MiB of temporary atlases while it is built),
and the noise volumes at about 18 MiB.

Recover in this order:

1. Lower `Render resolution`.
2. Lower `Ray quality`.
3. Lower `Noise detail`.
4. Lower `Grid density` if entering smoke, or rebuilding it after a filter or
   recolour, is the slow part; it also cuts the samples every ray takes.
5. Return to `Points`; after a context loss, reload the page before retrying.

If nothing needs to move, set `Animation speed` to 0 (`Still`): a settled
still view draws no frames at all. Changing lighting first rarely addresses
the dominant cost. For platform and GPU checks, see
{doc}`../a_orientation/02_system_requirements`; for general interaction
failures, see {doc}`06_troubleshooting_core_interactions`.

## Related workflows

- {doc}`04_view_layout_live_snapshots_small_multiples` — points-only multiview
- {doc}`../i_vector_field_velocity/index` — points-only directional overlay
- {doc}`../e_filtering/01_filtering_mental_model` — the visible subset used to
  build smoke
- {doc}`../k_figure_export/index` — explicit, reproducible figure output
