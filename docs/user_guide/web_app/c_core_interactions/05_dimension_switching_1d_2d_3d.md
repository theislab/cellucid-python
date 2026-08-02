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
- **The navigation mode**, unless you have chosen one yourself: 3D selects
  **Orbit**, 1D and 2D select **Planar** — see below.
- Any **dimension-specific overlays** update (e.g., vector fields/velocity overlays).
- Some derived rendering structures (centroids, spatial indices) may rebuild; large datasets can take a moment.

### Usually does not change

- Your dataset choice (obviously).
- Most UI selections (fields, filters) unless a feature is inherently dimension-specific.
- **A navigation mode you chose yourself.** Once the **Mode** dropdown differs
  from the mode the current dimension implies, it is yours and dimension changes
  leave it alone.

:::{tip}
If switching to a new dimension leaves you “lost in space” (empty-looking view), click **Reset Camera**.
:::

---

## What switching actually looks like

The three frames below are one session on the Pancreas sample (which publishes
all three embeddings), changed only with the **Dimension** dropdown.

```{figure} ../../../_static/screenshots/web_app/window-dimension-3d.png
:alt: The whole application window with a three-dimensional point cloud inside a light grid box, the sidebar showing a DIMENSION dropdown reading 3D and a MODE dropdown reading Orbit with orbit options beneath it.
:width: 1440px

**3D**, the Pancreas default. `Dimension: 3D`, `Mode: Orbit`.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-2d.png
:alt: The same window after choosing 2D; the cells now lie on a single plane seen face-on inside the grid box, the DIMENSION dropdown reads 2D, and the MODE dropdown has followed it to Planar with the planar options — KEYBOARD SPEED, ZOOM TO CURSOR, INVERT AXES — beneath it, with the mouse pointer pressing the dimension dropdown.
:width: 1440px

**2D.** The coordinates changed — the cells are now coplanar — and the camera
followed. `Mode` moved to `Planar` on its own, so the embedding is seen
face-on rather than at an angle.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-1d.png
:alt: The same window after choosing 1D; the cells are strung along a single straight line of colour segments crossing the grid box, the DIMENSION dropdown reads 1D, and the MODE dropdown reads Planar with the planar options beneath it.
:width: 1440px

**1D.** Every cell keeps its colour but has only one coordinate, so the whole
dataset becomes a coloured line — an ordering, not a map. `Mode` is `Planar`
here too: only a 3D embedding gets Orbit.
```

:::{important}
**Switching dimension switches navigation mode with it.** The app carries a rule
that 3D uses **Orbit** and 1D/2D use **Planar**, and it re-applies that rule on
every dimension change, as the frames above show.

The rule stops applying to a view as soon as the **Mode** dropdown differs from
what the current dimension implies — pick **Free-fly**, or **Orbit** while you
are in 2D, and that choice is yours to keep across later dimension changes. Put
**Mode** back to the dimension's own default and the rule takes over again.
:::

---

## Edge cases

- **You don’t see the dimension dropdown**: the dataset provides only one
  embedding dimension. None of the five built-in samples do this — all five
  publish 1D, 2D and 3D — so if you meet it, you are on your own export.
- **Switching to 2D/1D does nothing**: that dimension may not exist for this dataset.
- **Switch fails with an error**: the embedding data may be missing/corrupt; check the console and the data loading docs.
- **The picture changed but the camera feels wrong**: you have chosen a
  navigation mode at some point in this session, so the dimension rule no longer
  touches this view. Set **Mode** back to the dimension's default — **Planar**
  for 1D/2D, **Orbit** for 3D — and it resumes following the dimension.

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
