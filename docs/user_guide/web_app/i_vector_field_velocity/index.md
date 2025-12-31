# Vector Field / Velocity Overlay (GPU Particle Overlay)

This section documents Cellucid’s **GPU particle overlay** for visualizing **per-cell vector fields** (e.g., RNA velocity, drift/displacement vectors) on top of an embedding.

**Audience:** wet lab / non-technical, computational users, power users  
**Time:** 10–30 minutes (depending on depth)  
**Prerequisites:** a dataset with at least one vector field for the current dimension (1D/2D/3D)

:::{important}
The overlay is primarily a **qualitative** visualization of directionality. Always cross-check interpretation with your underlying velocity/drift methodology and uncertainty diagnostics.
:::

---

## Recommended reading order

1) `01_what_vector_fields_are_user_facing` (what you’re looking at)
2) `02_enabling_overlay_and_selecting_field` (how to turn it on reliably)
3) `03_core_parameters_document_exact_ui_labels` (the knobs most users touch)
4) `05_performance_and_quality` (if it looks wrong or runs slow)
5) `06_edge_cases` and `07_troubleshooting_velocity_overlay` (when something breaks)
6) `04_advanced_parameters_document_every_setting` (polish / cinematic tuning)

---

## Pages in this section

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`info;1.5em;sd-mr-1` What Vector Fields Are
:link: 01_what_vector_fields_are_user_facing
:link-type: doc

Mental model and what data is required for the overlay to exist.
:::

:::{grid-item-card} {octicon}`check;1.5em;sd-mr-1` Enable + Select Field
:link: 02_enabling_overlay_and_selecting_field
:link-type: doc

Where the toggle lives, what disables it, and how the dropdown behaves across 1D/2D/3D.
:::

:::{grid-item-card} {octicon}`sliders;1.5em;sd-mr-1` Core Parameters
:link: 03_core_parameters_document_exact_ui_labels
:link-type: doc

Density, speed, lifetime, size, opacity, colormap, and LOD syncing (with exact UI labels).
:::

:::{grid-item-card} {octicon}`tools;1.5em;sd-mr-1` Advanced Parameters
:link: 04_advanced_parameters_document_every_setting
:link-type: doc

Trail, bloom/HDR, color grading, and cinematic effects—every setting documented.
:::

:::{grid-item-card} {octicon}`pulse;1.5em;sd-mr-1` Performance + Quality
:link: 05_performance_and_quality
:link-type: doc

What scales with particles/trails, recommended settings by hardware, and “why it looks weird” checks.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting
:link: 07_troubleshooting_velocity_overlay
:link-type: doc

Symptom → diagnosis → fix for the most common overlay failures.
:::

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Screenshot Checklist
:link: 08_screenshots
:link-type: doc

Production checklist + capture specs for screenshots used in this section.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
