# View layout: Live + Snapshots (small multiples)

**Audience:** everyone (multiview is a core Cellucid superpower)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- What a “view” is (live vs snapshot)
- How “Keep view” works and what it preserves
- Grid compare vs edit selected view
- Camera locking/unlocking and per-view indicators (badges)

---

## Mental model: “views” are separate contexts

Cellucid supports multiple simultaneous **views** of the same dataset:

- **Live view**: your working view.
- **Snapshot views**: created by clicking **Keep view**; used for side-by-side comparison.

Each view can carry its own state (at minimum: camera + dimension; often also coloring/filter context depending on feature).

---

## Where multiview controls live

Multiview is controlled from:

- **Compare Views** accordion → **Multiview**

Key UI elements:

- **Keep view**: create a snapshot panel.
- **Locked Cam** toggle: link/unlink cameras across views.
- **View badges**: clickable list of views with indicators.
- **Layout** dropdown:
  - Grid compare
  - Edit selected view
- **Clear**: remove all kept views.

---

## Create a snapshot (Keep view)

### Fast path

1) Get the live view into a state you want to preserve (camera + coloring + filters).
2) Click **Keep view**.
3) A new view badge/panel appears.

### What “Keep view” preserves (practical)

At minimum, a snapshot preserves:

- the current camera framing and navigation mode,
- the current dimension (1D/2D/3D),
- the current coloring and transparency buffers (so the view looks the same immediately),
- metadata such as a short “filters summary” text (used for context).

:::{important}
Snapshots are points-mode only.

- If you are in smoke mode, **Keep view** will not create snapshots.
- If you already have snapshots, switching to smoke mode is blocked.
:::

---

## Layout: Grid compare vs Edit selected view

### Grid compare

- Shows all views at once (small multiples).
- Best for comparing clusters, thresholds, and alternative colorings.
- Click inside a panel to focus it for navigation (focus rules matter).

### Edit selected view

- Shows only the currently selected view, at full size.
- Best for making precise changes (selecting fields, adjusting filters).

Practical workflow:

1) Use **Edit selected view** to set up each view.
2) Switch back to **Grid compare** to compare them.

---

## Camera locking (linked vs independent)

The **Locked Cam** toggle controls whether views share a single camera.

### Cameras locked (linked)

- One camera controls all views.
- Navigation mode is shared.
- Great for “same framing, different coloring”.

### Cameras unlocked (independent)

- Each view stores its own camera and navigation mode.
- View badges show extra indicators (navigation mode per view).
- Great for “different framing per view” or mixing 2D planar + 3D orbit across panels.

:::{tip}
If grid compare feels confusing, keep cameras locked until you’re comfortable. Unlocking is powerful, but it makes “which panel am I controlling?” more subtle.
:::

---

## View badges (what the indicators mean)

Each view badge is a compact status line.

Common elements:

- **Number pill**: the view index.
- **⌖ camera indicator**: which view is currently the “active camera” when cameras are unlocked.
- **Label**: typically derived from the active field or your snapshot label.
- **Dimension badge**: e.g. `3D` (click to cycle dimensions if multiple are available).
- **Navigation badge** *(cameras unlocked only)*: `Orb`, `Pan`, or `Fly` (click to cycle).
- **× remove**: removes the view (or hides live view in some cases).

<!-- SCREENSHOT PLACEHOLDER
ID: multiview-badges-indicators
Where it appears: View layout → View badges
Capture:
  - Create at least 2 snapshots
  - Switch to Grid compare
  - Unlock cameras so nav badges (Orb/Pan/Fly) appear
  - Ensure at least one view shows a dimension badge (e.g., 3D)
Crop:
  - Include: the view badges list and enough canvas to show multiple panels
Annotations:
  - Call out: number pill, ⌖, dimension badge, nav badge, remove ×
Alt text:
  - View badge list showing per-view dimension and navigation indicators.
Caption:
  - View badges summarize each panel and expose per-view dimension/navigation controls when cameras are unlocked.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the view badges and indicators.
:width: 100%

View badges summarize each panel and expose per-view dimension/navigation controls when cameras are unlocked.
```

---

## Remove views: × vs Clear

- **× on a snapshot view**: removes that single snapshot.
- **Clear**: removes *all* snapshot views.

Removing a view deletes its per-view context (camera/dimension and any view-scoped state) so it cannot be recovered unless you saved a session.

---

## Edge cases

- **Layout dropdown disabled**: layout editing is points-mode only; smoke forces a single view.
- **Live view “disappears”**: when you remove the live view badge while multiple snapshots exist, the live view can be hidden to maximize compare space.
- **Badges don’t show dimension or nav indicators**:
  - dimension badge appears only if the dataset provides multiple dimensions,
  - nav badge appears only when cameras are unlocked.

---

## Troubleshooting (multiview)

For a full catalog, see `c_core_interactions/06_troubleshooting_core_interactions`.

### Symptom: “Keep view does nothing”

**Likely causes:**

- you are in smoke mode,
- render mode is not points,
- an internal error prevented snapshot creation (check console).

**Fix:**

1) Switch render mode to **Points**.
2) Click **Keep view** again.

### Symptom: “The wrong panel moves when I drag/scroll”

Click inside the panel you want to control first (focus rules in grid compare).
