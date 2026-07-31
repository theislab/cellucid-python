# Dimension switching (1D / 2D / 3D)

**Audience:** everyone (dimension switching is foundational)  
**Time:** 8–15 minutes  
**What you’ll learn:**
- What “dimension” means in Cellucid (it’s an embedding choice)
- How to switch dimensions (dropdown, badges, shortcuts)
- What changes (positions, overlays) vs what stays (most UI state)

---

## What “dimension” means

In Cellucid, “dimension” refers to which embedding dimensionality is active for a view:

- **1D**: a 1D coordinate per cell (often a timeline or pseudotime-like axis)
- **2D**: a 2D embedding (UMAP/tSNE/etc)
- **3D**: a 3D embedding (UMAP3D/etc)

Datasets can provide any subset of these. If only one dimension exists, the UI hides dimension controls.

---

## Where the dimension selector lives

Dimension controls live under:

- **Compare Views** accordion → **Dimension** → `1D / 2D / 3D` dropdown

The dropdown appears only when multiple dimensions are available for the current dataset.

---

## Three ways to switch dimensions

### 1) Dropdown (most explicit)

Use the **Dimension** dropdown to select `1D`, `2D`, or `3D` for the currently active/focused view.

### 2) View badge (fast, per view)

In multiview, each view badge may show a dimension badge (e.g. `3D`).

- Clicking the dimension badge cycles through available dimensions.

This is the fastest way to set different views to different dimensions.

### 3) Keyboard shortcuts (fastest)

- `1` → switch to 1D (if available)
- `2` → switch to 2D (if available)
- `3` → switch to 3D (if available)

These shortcuts are ignored while typing in inputs/selects.

---

## Per-view behavior (important)

Dimension is **per view**, not strictly global:

- the live view can be in 3D,
- a snapshot can be in 2D,
- another snapshot can be in 1D.

What the dropdown targets depends on your multiview context:

- In **Edit selected view**, it targets the selected view.
- In **Grid compare**, click inside the panel you want to change first (focus rule), then switch.

---

## What changes when you switch dimensions

### Changes

- **Point positions** change to the coordinates of the new embedding.
- Any **dimension-specific overlays** update (e.g., vector fields/velocity overlays).
- Some derived rendering structures (centroids, spatial indices) may rebuild; large datasets can take a moment.

### Usually does not change

- Your dataset choice (obviously).
- Most UI selections (fields, filters) unless a feature is inherently dimension-specific.
- Your camera mode setting (orbit/planar/free-fly), although the *feel* may be better if you switch modes to match the dimension.

:::{tip}
If switching to a new dimension leaves you “lost in space” (empty-looking view), click **Reset Camera**.
:::

---

## Interface reference

```{figure} ../../../_static/screenshots/web_app/dimension-navigation-controls-2d.png
:alt: Cellucid Compare Views controls showing 2D dimension and Planar navigation selected.
:width: 246px

For a 2D embedding, Cellucid selects 2D with Planar navigation; both settings remain explicit if the user chooses another valid configuration.
```

```{figure} ../../../_static/screenshots/web_app/dimension-2d-planar-default.png
:alt: Cellucid displaying a two-dimensional embedding with Planar navigation.
:width: 1440px

A two-dimensional embedding opens with Planar navigation and a flat, screen-aligned exploration surface.
```

---

## Edge cases

- **You don’t see the dimension dropdown**: the dataset likely provides only one embedding dimension.
- **Switching to 2D/1D does nothing**: that dimension may not exist for this dataset.
- **Switch fails with an error**: the embedding data may be missing/corrupt; check the console and the data loading docs.

---

## Troubleshooting (dimension switching)

For a full catalog, see {doc}`06_troubleshooting_core_interactions`.

### Symptom: “I can’t find 2D / the dropdown is missing”

**Likely cause:** dataset has only one available dimension.

**Fix:** confirm the export includes the missing embedding (data prep issue) or use a dataset that provides it.

### Symptom: “After switching, the plot looks empty”

**Likely causes:**

- camera is framed poorly for the new embedding,
- filters removed all visible points.

**Fix:**

1) Click **Reset Camera**.
2) Clear filters/outlier thresholds.
