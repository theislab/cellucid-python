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

---

## Fast path (wet lab / non-technical)

### What you’re looking at

Think of each cell as having a tiny “push” in a certain direction on the embedding. The overlay shows many small particles “floating” along those directions, so you can see:

- where the population tends to “flow”
- where trajectories appear to start/end
- where directionality is weak or ambiguous (flow looks noisy or faint)

### What it’s good for (real questions)

- “Do cells appear to move from cluster A toward cluster B?”
- “Is there a consistent direction along this continuum?”
- “Do two conditions have different directional trends?”

### What success looks like

When it’s working, you see:

- a particle flow that roughly follows your biological expectations (e.g., along a trajectory)
- flow that changes appropriately when you switch to a different vector field (if multiple exist)
- flow that respects filtering (if you hide cells, the overlay should mostly disappear from those regions)

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
- moving particles through space using the local per-cell displacement,
- drawing motion as trails (so you can see direction and persistence),
- coloring particles using a **colormap** (the default mapping is based on **vector magnitude**).

### Dimensionality is not optional (1D vs 2D vs 3D)

Vector fields are **dimension-specific**:

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

### Where the overlay data can come from (loading paths)

Vector fields can arrive via multiple data loading workflows:

- **Prepared exports**: vector field binaries live under `vectors/` and the dataset’s `dataset_identity.json` contains a `vector_fields` metadata block.
- **AnnData in the browser / server / Jupyter**: vector fields are detected from `adata.obsm` keys using Cellucid conventions (UMAP-based naming).

---

## Mini troubleshooting (common misconceptions)

### “Why don’t I see arrows?”

Cellucid’s overlay is not a static arrow plot by default. It’s a **particle flow visualization**. If you want arrow glyphs, you would typically generate them outside Cellucid (e.g., scVelo stream/arrow plots) or export vector fields and use a dedicated arrow rendering (future feature).

### “The flow looks random”

Common causes:

- vectors are noisy/low magnitude (biologically ambiguous or underfit model)
- vectors were computed on a different basis/order than the embedding you’re viewing
- the overlay settings (trail length, turbulence, bloom) are too aggressive for your data

Start with conservative settings (shorter trails, low turbulence) and verify row order/basis alignment.

---

## Next steps

- Continue to `02_enabling_overlay_and_selecting_field` to learn where the controls live and how dimension switching affects availability.

---

## Related internal references (implementation/spec)

- `cellucid/markdown/VECTOR_FIELD_OVERLAY_CONVENTIONS.md`
- `cellucid/markdown/VELOCITY_OVERLAY_PLAN.md`
- `cellucid/markdown/VELOCITY_OVERLAY_IMPLEMENTATION_PLAN.md`
