# Verified vector-field captures

These captures use the real **Pancreatic endocrinogenesis (scVelo)** sample:
3,696 cells with independent 1D, 2D, and 3D UMAP embeddings, a 2D PCA view,
and one `velocity` field bound to each of the three UMAPs. The images document the current interface; they do not
claim that the particle trails validate the underlying velocity model.

## Overlay controls

```{figure} ../../../_static/screenshots/vector_field_velocity/overlay-controls.png
:alt: The Pancreas sample in 2D with the velocity overlay running over the cells and the Vector field controls open in the sidebar.
:width: 1440px

The 2D **Planar** capture shows `velocity` with the initial 15K density,
3.0× flow speed, 8.0s lifetime, size 1.0, 60% opacity, `Viridis`, and
`Sync with LOD` enabled.
```

## Advanced visual settings

```{figure} ../../../_static/screenshots/vector_field_velocity/advanced-settings.png
:alt: The 2D velocity overlay with Advanced visual settings open in the sidebar: Particle rendering, Trail settings, HDR and bloom, Color grading and Cinematic effects.
:width: 1440px

The same real 2D field with **Advanced visual settings** expanded. Particle,
trail, HDR/bloom, color-grading, and cinematic controls change the rendering,
not the stored per-cell vectors.
```

## Standard Pancreas velocity

```{figure} ../../../_static/screenshots/vector_field_velocity/pancreas-velocity-3d.png
:alt: The Pancreas sample in its 3D view coloured by clusters with the velocity overlay running, and the Vector field controls open in the sidebar.
:width: 1440px

The independent 3D `velocity` payload in its declared **3D** / **Orbit**
default. Compare it with the 2D capture to inspect projection dependence; do not
interpret either animated trail as a reconstructed single-cell lineage.
```

For the sample's exact biological provenance and data contract, see
{doc}`../b_data_loading/10_standard_pancreas_dataset`. For control semantics,
continue with {doc}`03_core_parameters_document_exact_ui_labels` and
{doc}`04_advanced_parameters_document_every_setting`.
