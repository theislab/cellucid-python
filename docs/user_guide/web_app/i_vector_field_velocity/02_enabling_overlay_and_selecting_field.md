# Enabling overlay and selecting field

**Audience:** everyone  
**Time:** 5–15 minutes  
**What you’ll learn:**
- Where the Vector Field Overlay controls live in the UI (exact labels)
- Why you sometimes cannot enable the overlay (and what to do)
- How the **Vector field:** dropdown is populated (dimension-specific)
- What happens when you switch 1D/2D/3D while the overlay is on

---

## Prerequisites (before you look for the toggle)

1) You must have a dataset loaded.
2) The dataset must contain at least one vector field for **some** dimension (1D/2D/3D).
3) You must be in **Render mode: Points** (the overlay controls are currently placed in the Points controls area).

If any of these is false, the UI may hide the overlay controls entirely or disable the toggle.

---

## Fast path (just make it work)

1) Open the left sidebar → **Visualization**.
2) Set **Render mode:** `Points`.
3) Scroll to **Vector Field Overlay:**.
4) If you see **Show overlay** enabled, check it.
5) In **Vector field:** select the field you want (if there is more than one).
6) Wait for the “Loading vector field…” state to finish, then look for animated particle flow on the canvas.

If the checkbox disables itself or the overlay doesn’t appear, jump to `07_troubleshooting_velocity_overlay`.

---

## Where the controls live (UI map, exact labels)

The overlay controls are in the left sidebar under:

- **Visualization**
  - **Render mode:** (must be `Points`)
  - **Vector Field Overlay:**
    - `Show overlay` (checkbox)
    - (appears only when enabled) **Vector field:** (dropdown)
    - (appears only when enabled) sliders + advanced settings

:::{note}
If you don’t see **Vector Field Overlay:** at all, that usually means the dataset has *no* vector fields (this is normal for many datasets).
:::

---

## When the toggle is disabled (and what that means)

There are two common “disabled” states:

### A) The entire overlay block is missing

**Meaning:** your dataset provides **no vector fields**.  
**Fix:** export/add vector fields in Python (or load a dataset that includes them).

### B) The overlay block is visible, but `Show overlay` is disabled

**Meaning:** vector fields exist, but **not for the current embedding dimension**.

In this case, Cellucid will show a hint like:

- `Vector fields available for 2D, 3D. Switch embedding dimension to enable.`

**Fix:** switch the active view to a supported dimension (e.g., 2D or 3D), then try again.

---

## Field dropdown behavior (what appears in “Vector field:”)

### Only fields available for the current dimension are listed

The **Vector field:** dropdown is filtered to the active view’s current dimension:

- in a 2D view → only fields with 2D vectors appear
- in a 3D view → only fields with 3D vectors appear

This prevents “selecting a field that cannot possibly render” for the current view.

### Label vs internal id

Each dropdown entry has:

- a user-facing **label** (e.g., `Velocity (UMAP)`)
- an internal **id** (often something like `velocity_umap`)

You will usually see only the label in the dropdown. The id matters when you prepare data in Python (see `01_what_vector_fields_are_user_facing` for the naming contract).

---

## What happens when you switch 1D/2D/3D while the overlay is on

Cellucid tries to keep the overlay enabled across dimension changes:

- If a vector field exists for the *new* dimension, Cellucid loads it and keeps the overlay on.
- If no vector field exists for the new dimension, Cellucid **auto-disables** the overlay (to avoid rendering nonsense).

This is expected behavior. Vector fields are dimension-specific.

:::{note}
For large datasets, dimension switching may trigger a brief “Loading vector field…” state again because the app may need to fetch a different vector array.
:::

---

## State and persistence (be explicit)

### Overlay enable state

- The `Show overlay` checkbox is a **global render setting** for the current dataset.
- When it is enabled, the overlay is drawn in the live view and in snapshot views (where applicable).

### Dataset changes

When you load a **new dataset**, the UI intentionally resets the overlay:

- the checkbox is turned off
- the settings panel collapses

This prevents accidentally carrying a velocity overlay from dataset A into dataset B.

### Snapshots (“Keep view”) and multiview

- Overlay settings are not “per-snapshot annotations”; they are a render layer.
- If you are comparing many views at once, enabling the overlay can increase GPU load (see `05_performance_and_quality`).

---

## Screenshots to capture (spec only)

See `08_screenshots` for the full checklist, but the minimum set for *this page* is:

- Overlay controls visible in **Visualization → Render mode: Points**
- A screenshot showing the **disabled** toggle with the “switch dimension” hint
- A “before/after” pair: overlay off vs on (same camera)
