# Camera controls (advanced)

**Audience:** power users, demo builders, and anyone who wants controls to “feel right”  
**Time:** 10–15 minutes  
**What you’ll learn:**
- Every camera/navigation knob in the UI (and when it applies)
- Keyboard shortcuts that affect camera/navigation
- What “Reset Camera” does (and what it doesn’t)
- Persistence: what is saved in sessions vs browser preferences

---

## Where camera/navigation controls live

All navigation/camera controls are in:

- **Compare Views** accordion → **Navigation** section

The controls shown depend on the current navigation **Mode**:

- Orbit controls (keyboard speed, Google Earth-style drag, orbit anchor)
- Planar controls (keyboard speed, zoom-to-cursor, invert axes)
- Free-fly controls (look sensitivity, move speed, pointer lock, invert look)

---

## Orbit controls (what each knob does)

In Orbit mode you rotate around an anchor.

- **Keyboard speed**: scales orbit rotation + zoom speed for keyboard navigation.
- **Google Earth-style drag** (default on): flips horizontal drag direction (A/D and mouse drag feel “grabby”).
- **Show orbit anchor** (default on): shows a compass-like overlay (helps orientation in 3D).

---

## Planar controls

Planar mode is a 2D map-like camera.

- **Keyboard speed**: scales planar pan + zoom speed for keyboard navigation.
- **Zoom to cursor (pinch-style)** (default on): zoom centers on cursor location (not the view center).
- **Invert axes**: flips pan direction (useful on some trackpad/remote setups).

---

## Free-fly controls

Free-fly is first-person movement.

- **Look sensitivity**: how fast mouse motion rotates the view.
- **Move speed**: movement speed through the scene.
- **Invert look axes**: flips the mouse look direction (both X and Y).
- **Capture pointer**: enables pointer lock (FPS-style mouse look without dragging).
  - Exit pointer lock with `Esc`.
  - Pointer lock can be blocked in iframes and some managed browser setups.
- **Projectile shooting** (optional): enables a “click/hold-to-charge” demo interaction.

---

## Keyboard shortcuts that affect camera/navigation

These work unless you’re typing in an input/select.

- `R`: reset camera (camera framing only)
- `O` / `P` / `G`: orbit / planar / free-fly
- `W A S D`: move/rotate/pan (depends on mode)
- `=` or `]`: zoom in (orbit/planar)
- `-` or `[`: zoom out (orbit/planar)
- `Shift`: move faster (all modes)
- `F`: fullscreen
- `H`: toggle sidebar
- `?`: open the in-app “Keyboard Shortcuts” section

:::{note}
If you press keys and nothing happens, click the canvas once (focus issues are the #1 cause). For more, see `c_core_interactions/06_troubleshooting_core_interactions`.
:::

---

## Reset behavior (important nuance)

There are two different “reset” concepts:

1) **Reset Camera button** (in the UI)
   - Resets the camera.
   - Also restores many visualization + navigation UI controls to their initial defaults (captured at app init).

2) **`R` keyboard shortcut**
   - Resets the camera framing only.
   - Does not necessarily reset all visualization sliders/checkboxes.

Practical implication: if you want “start over visually”, use the **button**; if you just want to re-center quickly, use `R`.

---

## Persistence (what gets remembered)

Cellucid uses multiple persistence layers:

- **Session bundles (`.cellucid-session`)**: Save/Load State captures the app state, including camera and many UI control values.
- **Browser preferences (`localStorage`)**: small UI preferences like background/theme can persist across reloads.

Pointer lock state is *not* a persistent preference:

- it is a runtime browser mode,
- it is disabled when you exit free-fly or press `Esc`.

---

## Deep path (exact slider mappings)

This section is for readers who want repeatable “feel” across machines and sessions.

### Look sensitivity (Free-fly)

- UI slider range: 1–30 (default 5)
- Internal mapping: `radians_per_pixel = slider_value * 0.0005`

### Move speed (Free-fly)

- UI slider range: 1–500 (default 100)
- Internal mapping: `units_per_second = slider_value / 100`

### Orbit keyboard speed

- UI slider range: 1–100 (default 40)
- Internal mapping: `multiplier = slider_value / 100`

### Planar keyboard speed

- UI slider range: 1–100 (default 100)
- Internal mapping: maps to a small pan-speed factor in approximately `[0.001, 0.0075]` (higher = faster).

---

## Screenshots (placeholders)

<!-- SCREENSHOT PLACEHOLDER
ID: camera-controls-navigation-panel
Where it appears: Camera controls → Where camera/navigation controls live
Capture:
  - Show Compare Views → Navigation expanded
  - Capture once per mode (Orbit/Planar/Free-fly) OR use one screenshot and call out the mode dropdown
Crop:
  - Include: Navigation controls + a small slice of canvas
Alt text:
  - Navigation controls panel showing the mode dropdown and mode-specific options.
Caption:
  - Navigation controls live under Compare Views and change based on the selected mode.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Navigation controls panel.
:width: 100%

Navigation controls live under Compare Views and change based on the selected mode.
```

---

## Troubleshooting (camera controls)

### Symptom: “WASD moves the page, not the camera”

**Likely causes:**

- the canvas does not have focus (the page is scrolling),
- a dropdown/slider is focused,
- you are typing in an input.

**Fix:**

1) Click the canvas once.
2) Press `Esc` to clear pointer lock/focus if needed.

### Symptom: “Capture pointer doesn’t work”

**Likely causes:**

- you are not in free-fly mode,
- the browser blocked pointer lock (iframe policy or managed settings),
- you didn’t trigger it with a user gesture.

**Fix:**

1) Switch to **Free-fly**.
2) Enable **Capture pointer**.
3) Click the canvas to request pointer lock.

For more symptoms, see `c_core_interactions/06_troubleshooting_core_interactions`.
