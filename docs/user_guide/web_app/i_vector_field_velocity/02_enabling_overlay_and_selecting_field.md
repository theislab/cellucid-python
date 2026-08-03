# Enable the overlay and select a field

**Audience:** everyone  
**Time:** 5–10 minutes

The vector-field overlay is **off by default** and lazy: Cellucid does not load
vector values or allocate the particle overlay until you enable it.

## Fast path

1. Load the dataset and open **Visualization**.
2. Set `Render mode:` to `Points`.
3. Select the required 1D, 2D, or 3D view under **Compare Views → Dimension:**.
4. Find **Vector Field Overlay:**.
5. If `Vector field:` is already visible, choose a field. This is required when
   the dataset has multiple fields and does not declare a default.
6. Check `Show overlay`.
7. Wait for `Loading vector field…`, then for particle preparation to complete.

The settings panel opens after the overlay is enabled. A small `i` button beside
the heading explains that the animation uses dimension-specific per-cell
velocity or drift vectors.

## Availability states

The UI has three deliberate states.

### No **Vector Field Overlay:** block

The loaded dataset declares no vector fields in any dimension. There is no
message for this state, and none is needed: the whole block is simply absent,
and `Show connectivity edges` is followed straight by **Reset Camera**.

```{figure} ../../../_static/screenshots/vector_field_velocity/read-dataset-without-vector-field.png
:alt: The bottom of the Visualization panel with the Kanemaru developing heart sample loaded: Show connectivity edges is followed straight by Reset Camera, with no Vector Field Overlay block between them.
:width: 466px

**Kanemaru** — no vector field. Of the five built-in samples only the Pancreas
sample ships one, so this is what four of the five look like.
```

```{figure} ../../../_static/screenshots/vector_field_velocity/read-dataset-with-vector-field.png
:alt: The same stretch of the Visualization panel with the Pancreas sample loaded: a ringed Vector Field Overlay block with an unchecked Show overlay box sits between Show connectivity edges and Reset Camera.
:width: 466px

**Pancreas** — the same stretch of the same panel, with the block present. This
is the only difference; if you are looking for the overlay and this block is not
there, the dataset does not have one.
```

Add a field using the exact preparation contract in
{doc}`../../python_package/c_data_preparation_api/08_vector_fields_velocity_displacement`,
then reload the dataset.

### Block shown, `Show overlay` disabled

Read the message beneath the controls:

- `Vector fields available for 1D, 2D, 3D. Switch embedding dimension to enable.`
  means no field matches the selected dimension. The dimension list reflects
  the dataset and can contain any subset of 1D, 2D, and 3D.
- `Available for 1D, 2D, 3D. Select a vector field before enabling the overlay.`
  means matching fields exist, but none is selected.

Only fields declared for the selected dimension appear in `Vector field:`.
Cellucid does not silently choose the first field. A declared default is selected
only when that same field supports the current dimension.

### Overlay enabled

The status becomes `Available for …`. `Vector field:` and all core/advanced
controls are visible. Changing the field shows `Loading vector field…`; on
success the overlay uses the newly selected field.

If loading fails, the error notification contains the underlying message, the
checkbox is turned off, and the status reads `Failed to load vector field.`.
Follow {doc}`07_troubleshooting_velocity_overlay` rather than repeatedly
toggling it.

## Dimension changes

Vector data is dimension-specific:

- If the selected field also exists in the new dimension, Cellucid loads that
  array and keeps the overlay enabled.
- If that field is unavailable but an applicable declared default exists, the
  dimension-filtered selector can use that default.
- If there is no valid selected/default field for the new dimension, the
  overlay turns off and asks for an explicit selection.
- If the dataset has no field at all for the new dimension, the overlay turns
  off and reports which dimensions are available.

For a large remote or AnnData-served dataset, the first switch to a dimension
can require a separate vector request. Wait for the current request to finish
before changing dimension or field again.

## Filters and multiview

Particles are seeded from the selected view's visible cells. Filtering all
cells produces a ready overlay with no visible spawn cells; restore visibility
before judging the field.

The overlay is a shared render feature, not a snapshot-specific annotation.
Multiview maintains separate per-view visibility and trail buffers, but the
field and visual-control settings are global. Additional visible panels increase
GPU work. See
{doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`.

The overlay controls are points-only. Selecting
`Volumetric smoke cloud (alpha)` hides them and uses the smoke renderer instead.

## Dataset and session behavior

Loading or replacing a dataset always:

- clears the selected vector field;
- turns `Show overlay` off;
- hides the settings panel until the new dataset publishes its field inventory.

For the same loaded dataset, a `.cellucid-session` captures `Show overlay`,
`Vector field:`, every core/advanced slider, `Color scheme:`, and
`Sync with LOD`. Restore still validates the current field options; the session
does not contain the vector arrays or dataset. See
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

## Verified interface

{doc}`08_screenshots` shows the real Pancreas `velocity_umap` field in 2D
**Planar**, the expanded advanced controls, and the independent 3D field in
**Orbit**.
