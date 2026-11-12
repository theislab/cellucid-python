# Enable the overlay and select a field

**Audience:** everyone  
**Time:** 5–10 minutes

The vector-field overlay is **off by default** and lazy: Cellucid does not load
vector values or allocate the particle overlay until you enable it.

## Fast path

1. Load the dataset and open **Visualization**.
2. Set `Render mode` to `Points`.
3. Show an embedding the field is bound to (the embedding the vectors were
   computed on): choose it under **Compare views → View N shows** (with a
   single view, **Shows**), or press `1`, `2` or `3` for its dimension.
4. Find **Vector field**.
5. If `Vector field` is already visible, choose a field. This is required when
   the dataset has multiple fields and does not declare a default.
6. Check `Show overlay`.
7. Wait for `Loading vector field…`, then for particle preparation to complete.

The settings panel opens after the overlay is enabled. A small `i` button beside
the heading explains that the animation uses per-cell velocity or drift
vectors, each bound to one embedding.

## Availability states

The UI has three deliberate states.

### No **Vector field** block

The loaded dataset declares no vector fields for any embedding. There is no
message for this state, and none is needed: the whole block is simply absent
from **Visualization**.

```{figure} ../../../_static/screenshots/vector_field_velocity/read-dataset-without-vector-field.png
:alt: The Connectivity block of the Visualization panel with the Kanemaru heart sample loaded, and no Vector field block under it, because the dataset has none.
:width: 496px

**Kanemaru** — no vector field. Of the five built-in samples only the Pancreas
sample ships one, so this is what four of the five look like.
```

```{figure} ../../../_static/screenshots/vector_field_velocity/read-dataset-with-vector-field.png
:alt: The same part of the Visualization panel with the Pancreas sample loaded: the Connectivity block followed by a Vector field block whose Show overlay box is unchecked.
:width: 496px

**Pancreas** — the same stretch of the same panel, with the block present. This
is the only difference; if you are looking for the overlay and this block is not
there, the dataset does not have one.
```

Add a field using the exact preparation contract in
{doc}`../../python_package/c_data_preparation_api/08_vector_fields_velocity_displacement`,
then reload the dataset.

### Block shown, `Show overlay` disabled

Read the message beneath the controls:

- `Vector fields are bound to UMAP · 2D, UMAP · 3D. Show one of them in this view to enable.`
  means no field is bound to the embedding the selected view shows. The list
  names every embedding some field is bound to as **View N shows** names it
  (two embeddings named alike there add their `obsm` key: `UMAP · 2D (X_umap_raw)`),
  so it tells you which embedding to show — another dimension, or another
  embedding of the same dimension in **View N shows** under **Compare views**.
  The dimmed `Show overlay` says the same when you point at it.
- `Available on UMAP · 2D, UMAP · 3D. Select a vector field before enabling the overlay.`
  means several fields are bound to the embedding the view shows and none is
  selected.

Only fields bound to the embedding the view shows appear in `Vector field`.
When exactly one is, it is selected: there is nothing to choose. Among several,
Cellucid does not guess: the dataset's declared default is selected when it is
bound to that embedding, and otherwise you pick one.

### Overlay enabled

The status becomes `Available on …`. `Vector field` and all core/advanced
controls are visible. Changing the field shows `Loading vector field…`; on
success the overlay uses the newly selected field. While a field loads, `Show
overlay` and `Vector field` wait (dimmed, "the vector field is loading"); a
newer choice arriving then (a session, an undo, an agent) is the one shown.

If loading fails, one error notification reads `The vector field could not be
shown: …` with the reason (a problem inside Cellucid says so, its own words
under **Details**), the checkbox is turned off, and the status reads `The vector field could not be shown.`.
Follow {doc}`07_troubleshooting_velocity_overlay` rather than repeatedly
toggling it.

## Dimension and embedding changes

A vector field is bound to the embedding it was computed on, so what matters is
the embedding the view shows, not only its dimension:

- If the selected field is also bound to the newly shown embedding, Cellucid
  loads that array and keeps the overlay enabled.
- If that field is unavailable but an applicable declared default exists, the
  selector can use that default.
- If there is no valid selected/default field for the new embedding, the
  overlay turns off and asks for an explicit selection.
- If no field at all is bound to the new embedding (for example a view switched
  to a t-SNE, or a dataset whose 2D slot was changed under **Session › Data
  mapping**), the overlay turns off and names the embeddings that have one.

For a large remote or AnnData-served dataset, the first switch to an embedding
can require a separate vector request. Wait for the current request to finish
before changing dimension, embedding or field again.

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
`Volumetric smoke cloud (beta)` hides them and uses the smoke renderer instead.

## Dataset and session behavior

Loading or replacing a dataset always:

- clears the selected vector field;
- turns `Show overlay` off;
- hides the settings panel until the new dataset publishes its field inventory.

For the same loaded dataset, a `.cellucid-session` captures `Show overlay`,
`Vector field`, every core/advanced slider, `Color scheme`, and
`Sync with LOD`. Restore still validates the current field options; the session
does not contain the vector arrays or dataset. See
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

## Figures

A figure exported while the overlay is on (PNG, SVG, the export preview)
draws the same field as a still picture: arrows on a grid, the mean of nearby
cells' vectors, colored by speed with the overlay's `Color scheme`, with a key
naming the field. See
{doc}`../k_figure_export/06_edge_cases` (The vector field: arrows on a grid).

## Verified interface

{doc}`08_screenshots` shows the real Pancreas `velocity` field in 2D
**Planar**, the expanded advanced controls, and the independent 3D field in
**Orbit**.
