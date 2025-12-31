# Troubleshooting (core interactions)

This page is a **symptom → diagnosis → fix** catalog for navigation, camera, render modes, multiview, and dimension switching.

If you are brand-new to the UI, it may be faster to skim `c_core_interactions/01_navigation_modes_orbit_planar_free_fly` first, then come back here.

---

## “Unstick me” (fast reset sequence)

When things feel broken, these are the safest global actions:

1) Press `Esc` (exits pointer lock and often clears confusing focus state)
2) Switch navigation mode to **Orbit** (`O`)
3) Click the canvas once (restore focus)
4) Click **Reset Camera**
5) If multiview is confusing, lock cameras and/or clear snapshots

If the canvas is blank or you see WebGL errors, jump to `a_orientation/02_system_requirements`.

---

## Symptom: “I can’t rotate / pan / zoom anymore”

### Likely causes (ordered)

1) You’re in **Planar** mode (rotation is disabled by design).
2) You’re in **Free-fly** mode and not actively looking (no drag / pointer lock confusion).
3) Pointer lock is active and you expected orbit/planar behavior.
4) The canvas does not have focus (scroll wheel/keys go to the page).
5) A highlight tool is actively capturing interactions.

### How to confirm

- Look at **Compare Views → Navigation → Mode** (Orbit/Planar/Free-fly).
- If the cursor is hidden and moving the mouse rotates the view, pointer lock is active.
- Click once on the canvas and try again (tests focus).

### Fix

1) Press `Esc`.
2) Switch to **Orbit** mode (`O`).
3) Click and drag on the canvas.
4) If you still feel “lost”, click **Reset Camera**.

### Prevention

- Use Planar for 2D, Orbit for 3D; don’t fight the mode.

---

## Symptom: “WASD / shortcuts don’t work”

### Likely causes

1) Focus is in an `input`, `select`, or `textarea` (shortcuts are disabled while typing).
2) The canvas never received focus (page scroll consumes keys).
3) A modal is open (some shortcuts are intentionally ignored).

### How to confirm

- Click the canvas and try `O` / `P` / `G`.
- Press `?` to open the in-app shortcut panel (if it opens, global shortcuts are working).

### Fix

1) Click the canvas.
2) Press `Esc` to dismiss pointer lock or close modals (if relevant).
3) Avoid leaving a dropdown focused (after changing nav mode, the UI attempts to blur it, but not all browsers behave identically).

---

## Symptom: “My cursor disappeared / I’m stuck in capture pointer”

### Likely causes

- Pointer lock is active (free-fly mode, Capture pointer enabled).

### How to confirm

- Mouse movement rotates the camera without dragging.

### Fix

1) Press `Esc` (exits pointer lock).
2) If you want to prevent re-entry, uncheck **Capture pointer** in free-fly controls.

### Prevention

- Only enable Capture pointer when you intentionally want FPS-style looking.

---

## Symptom: “The wrong panel moves in grid compare”

### Likely causes

1) You are interacting with the **focused view** (last clicked), not the view you’re looking at.
2) Cameras are unlocked, so each view has its own camera and switching focus swaps cameras.

### How to confirm

- Click inside the panel you want to control, then drag/scroll again.
- Check whether cameras are locked (Locked Cam toggle).

### Fix

- For simple comparisons, lock cameras.
- For independent cameras, develop the habit: click inside a panel before navigating it.

---

## Symptom: “Navigation badges (Orb/Pan/Fly) don’t show on view badges”

### Likely causes

- Cameras are locked (nav badges only appear when cameras are unlocked).

### How to confirm

- Toggle **Locked Cam** off; badges should gain navigation indicators.

### Fix

- Unlock cameras if you need per-view navigation modes.

---

## Symptom: “Dimension dropdown is missing”

### Likely causes

- The dataset provides only one embedding dimension (common for 3D-only exports).

### How to confirm

- If the view badge never shows `2D`/`1D` and the dropdown is hidden, the dataset likely lacks those embeddings.

### Fix

- Use a dataset with multiple embeddings, or re-export with additional dimensions (data prep issue).

---

## Symptom: “Keep view does nothing”

### Likely causes (ordered)

1) You are in **smoke render mode** (snapshots are points-only).
2) An internal error prevented snapshot creation (rare; check console).

### How to confirm

- Check **Visualization → Render mode**.
- If it’s “Volumetric smoke cloud”, snapshots are blocked.

### Fix

1) Switch render mode to **Points**.
2) Click **Keep view** again.

---

## Symptom: “I can’t switch to smoke mode”

### Likely cause

- You have snapshots (kept views). Smoke mode is intentionally blocked when snapshots exist.

### Fix

- Clear snapshots (Compare Views → Clear), then switch to smoke.

---

## Symptom: “Smoke mode is blank”

### Likely causes (ordered)

1) Cloud density is too low for your current dataset/visibility.
2) Almost all points are filtered out (nothing to render).
3) GPU/driver instability (rare; usually also affects points mode).

### How to confirm

1) Switch back to **Points** mode. If points are visible, data is present.
2) Clear filters/outlier thresholds and try smoke again.

### Fix

- Increase **Cloud density**.
- Reduce filters so something is visible.

---

## Symptom: “Smoke mode is extremely slow”

### Likely causes

- Grid density too high for your GPU.
- Ray quality too high.
- Render resolution too high (close to 1.0× or above).

### How to confirm

- Lower grid density and ray quality; if performance improves immediately, you found the bottleneck.

### Fix (fastest wins)

1) Lower **Grid density**.
2) Lower **Ray quality**.
3) Lower **Render resolution**.

---

## Symptom: “WebGL context lost”

This is almost always GPU memory pressure.

### Fix

- Reload (Cellucid requires a reload to safely reinitialize GPU resources).
- After reload, reduce GPU load (fewer views, points mode, lower quality settings).

Full details and prevention: `a_orientation/02_system_requirements` → “WebGL context lost”.

---

## Symptom: “After switching dimension, everything looks empty”

### Likely causes

- Camera framing is poor for the new embedding.
- Filters removed almost all visible points.

### Fix

1) Click **Reset Camera**.
2) Clear filters/outlier thresholds.
3) Try a different navigation mode (Planar for 2D, Orbit for 3D).

---

## Symptom: “Pointer lock / fullscreen doesn’t work in notebooks”

### Likely cause

- The viewer is embedded in an iframe that does not allow pointer lock/fullscreen.

### Fix

- Use the standalone web app context for those features, or adjust iframe permissions.

---

## Reporting a bug (what to include)

If you’re filing an issue or asking for help, include:

- your browser + version,
- your OS + GPU (from `chrome://gpu` or equivalent),
- whether you were in points or smoke mode,
- whether you had snapshots and whether cameras were locked,
- a screenshot of the relevant sidebar section (navigation/render/multiview),
- any console error text (copy/paste).
