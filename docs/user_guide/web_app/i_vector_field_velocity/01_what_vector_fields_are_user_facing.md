# What vector fields are (user-facing)

**Audience:** everyone (wet lab / beginner → expert)  
**Time:** 5–15 minutes  
**What you’ll learn:**
- What the vector field (“velocity”) overlay is, in plain language
- What data is required for it to exist (and why it sometimes doesn’t show up)
- How 1D/2D/3D “dimension switching” changes what vector fields are available
- How to interpret the overlay safely (what it is *and isn’t*)

---

## What this feature is (one sentence)

Cellucid’s **Vector Field Overlay** draws an animated **particle flow** on top of your embedding using **one vector per cell** (e.g., RNA velocity or transition drift vectors), so you can *see directionality* in a way that’s easy to understand.

## The same view, overlay off and on

```{figure} ../../../_static/screenshots/vector_field_velocity/01-embedding-without-overlay.png
:alt: The Pancreas 2D planar embedding coloured by cell_type, with the Vector Field Overlay block ringed in the sidebar and its Show overlay box unchecked.
:width: 1440px

**Off.** This is the ordinary point cloud. The **Vector Field Overlay** block is
present in the sidebar because this dataset ships a vector field; `Show overlay`
is unchecked, and no vector data has been downloaded yet.
```

```{figure} ../../../_static/screenshots/vector_field_velocity/02-embedding-with-overlay.png
:alt: The same Pancreas 2D planar embedding with Show overlay checked: bright particle trails stream across the embedding, the points are dimmed underneath, and the expanded overlay controls show velocity_umap at 15K particles and 3.0x flow speed.
:width: 1440px

**On.** The same cells, same camera. Checking `Show overlay` fetches the vector
payload, allocates the particle system, and starts the animation. The bright
streaks are particle *trails*, not cells and not arrows.
```

---

## Fast path (wet lab / non-technical)

### What you’re looking at

Think of each cell as having a tiny “push” in a certain direction on the embedding. The overlay shows many small particles “floating” along those directions, so you can see:

- where displayed per-cell directions are locally coherent
- where the field points along or across an annotated continuum
- where directionality is weak or ambiguous (flow looks noisy, mixed, or faint)

### What it’s good for (real questions)

- “Do cells appear to move from cluster A toward cluster B?”
- “Is there a consistent direction along this continuum?”
- “Do two conditions have different directional trends?”

For example, in the standard Pancreas sample you can color by `clusters`, show
`velocity_umap` in 2D **Planar**, and then inspect its independent 3D field in
**Orbit**. Treat agreement as a useful visual check, not as proof of endocrine
lineage or kinetic time.

### What success looks like

When it’s working, you see:

- a stable particle flow that follows the supplied vectors
- flow that changes appropriately when you switch to a different vector field (if multiple exist)
- flow that respects filtering (particles are seeded only from currently visible cells)

:::{important}
The animation is primarily **qualitative**. It is not a substitute for proper velocity uncertainty checks, stream plots, or model diagnostics.
:::

---

## Practical path (computational users)

### What a “vector field” means in Cellucid

In Cellucid, a vector field is:

- a **per-cell displacement vector** (one vector per cell),
- expressed in the **same embedding coordinate system** as your points,
- with dimensionality **matching the view** (1D/2D/3D).

In other words, if your embedding is `(n_cells, 2)`, then your vector field must be `(n_cells, 2)` and each row corresponds to the same cell as the embedding row.

### What the overlay does with that data

Cellucid turns your per-cell vectors into a visualization by:

- spawning particles from **visible cells** (it respects filters/visibility),
- moving each particle using its associated cell's dimension-matched vector,
- drawing motion as trails (so you can see direction and persistence),
- coloring particles using a **colormap** (the default mapping is based on **vector magnitude**).

The renderer does not reconstruct a continuous biological trajectory between
cells, infer transition probabilities, or follow a real cell through time.
`Flow speed:`, `Trail length:`, turbulence, and post-processing change the
animation, not the stored vector values.

### A trail is not a streamline

This is the single most important thing to understand before interpreting the
picture.

A particle is spawned from one cell and **keeps that cell for its whole life**.
Every frame it re-reads *that same cell's* vector and moves along it. It never
looks up the field at the position it has reached. Each trail is therefore a
straight ray in the direction of one cell's vector, faded in and out over its
lifetime — not a curve traced through the field.

Concretely:

- **A trail is not an integral curve.** Following a streak from one cluster into
  another follows one cell's vector extrapolated forwards. It is not a path the
  field implies, and it is not a trajectory.
- **Curvature you see is not in the data.** What bends the streaks is
  `Turbulence:`, a divergence-free noise term added to the motion. It ships
  **on**, at `0.30` of its `0–1` range. Set it to `0` before reading the field's
  geometry.
- **Coherence between neighbouring streaks is real.** Neighbouring cells whose
  vectors agree produce parallel streaks; disagreement produces a visible mess.
  That contrast — coherent versus incoherent — is the honest signal here.

### What the colour and the size mean

Particle colour is the magnitude of the spawn cell's vector, divided by the
**largest magnitude in this field for this dimension**, clamped to `0..1`, and
looked up in `Color scheme:`. Particle size carries the same quantity, and the
per-frame step length is proportional to the raw magnitude times `Flow speed:`.

Two consequences:

- **The ramp is relative to this exact payload.** The same colour in two
  datasets, or in the same dataset at two different embedding dimensions, does
  not mean the same magnitude. There is deliberately no numeric legend or colour
  bar for it: the number would have no interpretable unit.
- **Very slow cells are under-represented.** Particles below a floor magnitude
  are discarded outright, and slow particles are recycled sooner than fast ones.
  A low-magnitude region therefore looks *sparser* as well as dimmer. Sparse
  does not mean "fewer cells".

### What the overlay never claims

| It does not show | Why |
|---|---|
| A rate, a speed, or a time | The stored numbers are normalised embedding units per arbitrary step. `Flow speed:` is a display gain you choose. |
| Any motion out of the embedding plane | The vectors on disk are already a projection. The component of the high-dimensional velocity that the embedding cannot express is gone before Cellucid sees the data. |
| A trajectory, lineage, or transition probability | Nothing in the renderer integrates the field or follows a cell. |
| Magnitudes comparable across dimensions or datasets | See above; and for the built-in Pancreas sample the 1D payload is on a different scale from the 2D and 3D payloads. |
| Uncertainty | There is no confidence channel in the format and none in the drawing. |

### Dimensionality is not optional (1D vs 2D vs 3D)

Vector fields are **dimension-specific**:

- a 1D vector field can be shown only in a **1D** embedding view
- a 2D vector field can be shown only in a **2D** embedding view
- a 3D vector field can be shown only in a **3D** embedding view

If your dataset only contains a 2D field and you switch the view to 3D, the overlay is expected to become unavailable (or auto-disable).

### How to tell whether your dataset contains vector fields

In the web app UI:

- If the dataset has **no vector fields**, the **Vector Field Overlay** block will not appear.
- If the dataset has vector fields but not for the **current dimension**, the block appears, but **Show overlay** is disabled and the UI tells you which dimensions are supported (e.g., “Vector fields available for 2D, 3D. Switch embedding dimension to enable.”).

---

## Deep path (experts / developers)

### Contract: coordinate space and row order

The overlay assumes:

- **Row order matches**: the `i`-th row in the vector field corresponds to the `i`-th cell in the embedding (and the `i`-th cell in the dataset as loaded).
- **Same coordinate space**: vectors are interpreted in the same embedding coordinate space as the points for that dimension.

If either assumption is violated, the overlay may still render but it will be misleading (often the hardest failure mode to notice).

### Normalization/scaling (why magnitudes are “relative”)

Embeddings are normalized for stable WebGL rendering. Vector fields are scaled to match that normalized render space, which means:

- direction is meaningful (if your vectors are correct)
- **absolute magnitudes are not directly comparable** to your original units unless you account for normalization

This is why the overlay exposes a **Flow speed** multiplier: it lets you tune visual motion without claiming biological time units.

The scaling is deliberately boring, and knowing that is what lets you trust the
relative magnitudes:

- The exporter validates the array (finite, right shape, right row count) and
  applies **one isotropic scale factor** — the same one it applies to the
  embedding itself. No filtering, no smoothing, no unit-normalisation, no
  per-cell rescaling.
- The app applies the same kind of scale again at load, which is exactly `1.0`
  for a prepared export and the real normalisation for a direct
  `h5ad://`/`zarr://` source. It then computes one global maximum magnitude,
  used only for colour and size.

So **relative magnitudes are preserved end to end** and **absolute units are
lost**. What arrives on disk is what gets drawn, up to a single multiplier.

:::{warning}
That is a statement about Cellucid, not about your vectors. Whether the numbers
on disk mean anything is decided entirely upstream, by whatever produced the
projection. For the built-in Pancreas sample the answer is documented in
{doc}`../b_data_loading/10_standard_pancreas_dataset`, and it is not
"a rate": the projection was computed with scVelo's `retain_scale=False`, so
the arrow length reflects the directional coherence of the local transition
kernel rather than a speed. Its 2D and 3D payloads were additionally passed
through scVelo's quiver autoscaler and its 1D payload was not, so the three
dimensions of that one sample are not on a common scale.
:::

### Where the overlay data can come from (loading paths)

Vector fields can arrive via multiple data loading workflows:

- **Prepared exports**: vector field binaries live under `vectors/` and the dataset’s `dataset_identity.json` contains a `vector_fields` metadata block.
- **AnnData in the browser / server / Jupyter**: vector fields are detected from
  exact dimension-suffixed `adata.obsm` keys such as `velocity_umap_2d`.

The preparation contracts are documented in
{doc}`../../python_package/c_data_preparation_api/08_vector_fields_velocity_displacement`;
loading-path checks are in
{doc}`../b_data_loading/01_loading_options_overview`.

---

## Mini troubleshooting (common misconceptions)

### “Why don’t I see arrows?”

Cellucid’s vector overlay is a **particle flow visualization**; it does not
render static arrow glyphs. Generate a separate stream or arrow plot with a
tool such as scVelo when that representation is required.

### “The flow looks random”

Common causes:

- vectors are noisy/low magnitude (biologically ambiguous or underfit model)
- vectors were computed on a different basis/order than the embedding you’re viewing
- the overlay settings (trail length, turbulence, bloom) are too aggressive for your data

Start with conservative settings (shorter trails, low turbulence) and verify row order/basis alignment.

---

## Next steps

- Continue to {doc}`02_enabling_overlay_and_selecting_field` to learn where the
  controls live and how dimension switching affects availability.
- Compare the real Pancreas 2D and 3D captures in {doc}`08_screenshots`.

---

## Implementation references

- `cellucid/assets/js/data/vector-field-manager.js`
- `cellucid/assets/js/app/ui/modules/velocity-overlay-controls.js`
- `cellucid/assets/js/rendering/overlays/velocity/velocity-overlay.js`
