# UI glossary (terminology)

**Audience:** everyone (this is the shared vocabulary)  
**Time:** 10–20 minutes (skim now, return later)  
**What you’ll get:** exact meanings of UI terms so the rest of the docs feel predictable

---

## How to use this page

- If you’re new: read the **bold terms only**.
- If you’re debugging something subtle: read the “Scope / persistence” notes under each section.

:::{tip}
This glossary is written to match the app UI language, even when a “technical” term might be more precise. The goal is: when the docs say “view” or “snapshot”, you know exactly what screen element/state is meant.
:::

---

## UI map (one screenshot you’ll reuse everywhere)

<!-- SCREENSHOT PLACEHOLDER
ID: ui-glossary-ui-map
Where it appears: UI glossary → UI map
Capture:
  - Any dataset loaded
  - Left sidebar open
  - Show: plot canvas + sidebar accordions (Session, Coloring & Filtering, Compare Views, Navigation, Highlighting, Analysis, Figure Export)
Crop:
  - Include: full sidebar + ~60% of canvas
  - Exclude: browser chrome and private dataset identifiers
Annotations:
  - 6–12 callouts naming the core UI regions (sidebar, canvas, view badges, legend)
Alt text:
  - Cellucid interface showing the canvas and main sidebar sections.
Caption:
  - Cellucid UI overview: the left sidebar controls state; the canvas renders the active view(s).
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a Cellucid UI overview.
:width: 100%

Cellucid UI overview: the left sidebar controls state; the canvas renders the active view(s).
```

---

## Data model terms

- **Dataset**: the loaded data bundle (points + metadata + optional gene expression + optional connectivities/vector fields).
- **Cell / point**: one row/observation rendered as a point. (The UI often says “points” because the viewer is generic.)
- **Gene / feature**: one variable/feature; used for gene expression coloring and some analyses.
- **Embedding**: coordinates for points in 1D/2D/3D (UMAP/tSNE/PCA/etc).
- **Dimension (1D / 2D / 3D)**: which embedding dimensionality you are currently viewing for a given view.
  - **4D**: present in some UI lists as “reserved”; not currently supported.
- **Connectivity**: a neighbor graph (edges). When enabled, edges are drawn between visible points only.
- **Vector field / velocity overlay**: an optional overlay that visualizes per-cell vectors (dimension-specific).

---

## Fields, coloring, and legends

- **Field**: something you can color by.
  - **Categorical obs**: discrete labels (clusters, batch, sample).
  - **Continuous obs**: numeric values (QC metrics, scores).
  - **Gene expression / var**: expression of a selected gene/feature.
- **Active field**: the currently selected field that drives coloring (per view).
- **Legend**: the UI element that explains the current coloring.
  - Categorical legend shows categories and their colors.
  - Continuous legend shows a numeric color scale.
- **Color-by**: the act of selecting a field to drive coloring.

Scope note:

- In multiview, the “active field” is typically **per view** (each view can show a different coloring).

---

## Views, snapshots, and “small multiples”

These terms are the most important for understanding Cellucid.

- **View**: one panel that renders the dataset with its own state (camera, dimension, coloring, filters, highlights).
- **Live view**: your default working view (often labeled “All cells” in the view badges).
- **Snapshot view** (also “kept view”): a view created by clicking **Keep view**; used to compare multiple states side-by-side.
- **Multiview**: the overall feature of having multiple views/snapshots.
- **View layout**
  - **Grid compare**: shows all views in a grid.
  - **Edit selected view**: shows only the active view so you can edit it precisely.
- **Active view**: the view you are editing/configuring (the one selected in the view badges).
- **Focused view**: in grid mode, the view under your last click; this determines which view the camera controls target.
- **View badge**: the clickable pill/row representing a view (with indicators like `3D`, `Orb/Pan/Fly`).
- **Cameras locked**: all views share one camera (navigation is synchronized).
- **Cameras unlocked**: each view has its own camera and navigation mode; badges show per-view navigation indicators.

Scope note:

- When cameras are unlocked, Cellucid stores camera state **per view** and switches it when you focus a different view.

---

## Navigation and camera terms

- **Orbit**: rotate around an anchor (best for 3D).
- **Planar**: pan/zoom as if looking at a flat map (best for 2D).
- **Free-fly**: first-person navigation (immersive; supports pointer lock).
- **Pointer lock / Capture pointer**: a browser feature that hides the cursor and reports raw mouse movement (needed for FPS-style looking).
- **Orbit anchor**: a visual compass/anchor indicator shown in orbit mode (optional).
- **Reset Camera**: returns the camera to a default framing (and can reset related UI controls depending on how you trigger it).

---

## Highlighting and selection terms

- **Highlight mode**: the active interaction tool for selecting points (e.g., lasso).
- **Highlight group**: a named set of highlighted cells (used for comparison/analysis).
- **Highlight page**: a collection of highlight groups (useful when you want multiple “sets” of groups).

Scope note:

- Highlights can be global or per-view depending on the feature; always check the page you’re using (highlighting is documented under `f_highlighting_selection/index`).

---

## Analysis and export terms

- **Analysis**: computations driven by the current dataset state (often highlights/groups).
- **Figure Export**: export an image/vector figure of the current view(s) suitable for papers.
- **Save State / Load State**: save/restore a `.cellucid-session` bundle that captures the application state.

---

## Scope & persistence vocabulary (used throughout the docs)

- **Global**: one value shared across the entire app/dataset (not specific to a view).
- **Per view**: stored separately for each view/snapshot.
- **Session bundle**: a saved `.cellucid-session` file (explicit user action).
- **Local storage**: persistent browser storage used for small preferences (e.g., theme/background).
- **Session storage**: storage that is cleared when the tab closes (used for sensitive tokens like GitHub OAuth in community annotation).

---

## Next steps

- Choose a workflow: `a_orientation/05_which_workflow_is_for_me_decision_tree`
- Learn navigation + multiview: `c_core_interactions/index`
