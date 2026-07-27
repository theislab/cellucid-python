# Verified vector-field captures

These captures use a synthetic six-population 2D fixture produced by the current
Python exporter. The per-cell `velocity_umap_2d` vectors define a smooth
counter-clockwise flow; no biological interpretation is implied.

## Overlay controls

```{figure} ../../../_static/screenshots/vector_field_velocity/overlay-controls.png
:alt: A synthetic 2D dataset with the Vector Field Overlay enabled and animated particle-flow controls visible.
:width: 1440px

The enabled overlay exposes field, density, speed, trail, particle size,
opacity, palette, and LOD synchronization while the flow remains visible.
```

## Advanced visual settings

```{figure} ../../../_static/screenshots/vector_field_velocity/advanced-settings.png
:alt: A vector-field particle flow with Advanced Visual Settings expanded.
:width: 1440px

Advanced Visual Settings expose particle, trail, bloom, and color-grading
controls without changing the underlying vector values.
```

## Standard Pancreas velocity

```{figure} ../../../_static/screenshots/vector_field_velocity/pancreas-velocity-3d.png
:alt: Cellucid showing the Pancreatic endocrinogenesis (scVelo) sample in its 3D Orbit view, colored by clusters, with the velocity_umap overlay enabled and vector controls visible.
:width: 1440px

The public Pancreas sample exercises a real dimension-matched
`velocity_umap` payload. This Build 2026-07-26.4 capture shows the enabled
overlay together with its field controls and the dataset's declared **3D** /
**Orbit** defaults.
```
