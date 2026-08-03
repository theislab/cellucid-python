# Navigation modes (Orbit / Planar / Free-fly)

**Audience:** everyone (new users + power users)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- What each navigation mode is best for
- Exact mouse/trackpad/keyboard controls (as implemented)
- Focus rules in multiview (which view you are actually controlling)
- The “gotchas” that make people think the app is broken

---

## Pick the right mode (fast guidance)

| Mode | Best for | Mental model | When it feels wrong |
|---|---|---|---|
| **Orbit** | 3D exploration | “Rotate around an anchor” | Trying to treat 2D like a map |
| **Planar** | 2D embeddings | “Pan/zoom a flat map” | Expecting 3D rotation |
| **Free-fly** | Immersive movement / demos | “First-person camera” | Using touch-only devices / pointer lock confusion |

Switch modes in **Compare Views → Navigation → Mode** (or use shortcuts `O`/`P`/`G`).

---

## Where the controls live in the UI

Navigation mode controls are in the left sidebar:

- **Compare Views** accordion → **Navigation** section → **Mode** dropdown

The mode options are:

- Orbit
- Planar
- Free-fly (shown as “Free-fly” in the UI, internally stored as `free`)

Keyboard shortcuts (global, unless you’re typing in an input):

- `O` → Orbit
- `P` → Planar
- `G` → Free-fly

You do not have to set the mode at all: the view's dimension picks one for you —
**Orbit** in 3D, **Planar** in 1D and 2D — and keeps doing so every time the
dimension changes. Choosing a mode that differs from that takes the control over
for that view. See {doc}`05_dimension_switching_1d_2d_3d`.

The **Mode** dropdown is followed by a box of options that belong to the chosen
mode only, so the panel looks different in each of the three. A still cannot
show a camera moving — but it can show you exactly which knobs you have, and
that is what people actually get lost looking for.

```{figure} ../../../_static/screenshots/web_app/navigation-controls-orbit.png
:alt: A panel headed NAVIGATION with a small circled i button, a MODE dropdown reading Orbit with the mouse pointer pressing it, and a bordered box below containing a KEYBOARD SPEED slider reading 0.40x, a ticked GOOGLE EARTH-STYLE DRAG checkbox and a ticked SHOW ORBIT ANCHOR checkbox.
:width: 472px

**Orbit** — a keyboard speed slider and two checkboxes.
```

```{figure} ../../../_static/screenshots/web_app/navigation-controls-planar.png
:alt: The same NAVIGATION panel with the MODE dropdown reading Planar and a bordered box below containing a KEYBOARD SPEED slider reading 0.0075x, a ZOOM TO CURSOR (PINCH-STYLE) checkbox and an INVERT AXES checkbox.
:width: 472px

**Planar** — a keyboard speed slider and two different checkboxes. Orbit and
Planar both label their slider `Keyboard speed:`; the value beside it is on a
different scale in each, because one scales rotation and the other pan.
```

```{figure} ../../../_static/screenshots/web_app/navigation-controls-freefly.png
:alt: The same NAVIGATION panel with the MODE dropdown reading Free-fly and a taller bordered box below containing a LOOK SENSITIVITY slider reading 0.05x, a MOVE SPEED slider reading 1.00 u/s, an INVERT LOOK AXES checkbox, a PROJECTILE SHOOTING checkbox and a CAPTURE POINTER checkbox with a small circled i beside it.
:width: 472px

**Free-fly** — two sliders, three checkboxes, and the pointer-lock checkbox
that does not exist in the other two modes.
```

:::{note}
The three boxes are mutually exclusive: only the one belonging to the current
mode is present, so a setting you remember can genuinely vanish when you change
mode. It is not lost — switch back and it returns with its value.
:::

---

## Orbit mode (3D exploration)

Orbit is the default for 3D embeddings: you rotate the camera around a target (“anchor”) and zoom in/out.

### Mouse / trackpad (Orbit)

- **Left-click + drag**: rotate
- **Right-click + drag** *(or Shift + drag)*: pan
- **Scroll wheel / trackpad scroll**: zoom

Orbit has inertia: quick drags can “coast” briefly after mouseup.

### Keyboard (Orbit)

Keyboard navigation works in all modes. In orbit mode:

- `W` / `S`: rotate up/down
- `A` / `D`: rotate left/right (respects “Google Earth-style drag”)
- `=` or `]`: zoom in
- `-` or `[`: zoom out
- `Shift`: move faster

### Touch (Orbit)

- **1 finger drag**: rotate
- **2 finger pinch**: zoom

### Orbit options (UI)

In **Compare Views → Navigation → Orbit**:

- **Keyboard speed**: scales keyboard rotation/zoom speed.
- **Google Earth-style drag** (default on): flips the horizontal drag sign (matches “grab-and-rotate” muscle memory for many users).
- **Show orbit anchor**: shows a compass/anchor overlay in orbit mode.



---

## Planar mode (2D “map” navigation)

Planar mode is designed for 2D embeddings: you pan and zoom like a map. Rotation is disabled.

### Mouse / trackpad (Planar)

- **Click + drag**: pan (always)
- **Scroll wheel / trackpad scroll**: zoom

If **Zoom to cursor (pinch-style)** is enabled (default), the zoom centers around your cursor location rather than the view center.

### Keyboard (Planar)

- `W` / `S`: pan up/down
- `A` / `D`: pan left/right
- `=` or `]`: zoom in
- `-` or `[`: zoom out
- `Shift`: move faster

Planar keyboard speed is controlled by the **Keyboard speed** slider in the Planar controls.

### Touch (Planar)

- **1 finger drag**: pan
- **2 finger pinch**: zoom

### Planar options (UI)

In **Compare Views → Navigation → Planar**:

- **Keyboard speed**: scales keyboard pan/zoom speed.
- **Zoom to cursor (pinch-style)** (default on): zoom around cursor.
- **Invert axes**: flips pan direction (useful for unusual trackpad/remote desktop setups).

---

## Free-fly mode (first-person navigation)

Free-fly is an immersive mode (especially for demos and 3D). You “look around” and move through the scene.

### Mouse (Free-fly)

Two ways to look around:

1) **Click + drag** on the canvas (cursor stays visible)
2) Enable **Capture pointer** (pointer lock): cursor is hidden, and mouse movement rotates the view without dragging

To exit pointer lock, press `Esc`.

:::{note}
Free-fly touch gestures are intentionally disabled in the current implementation (touch events return early in free-fly). If you’re on a touch-only device, use Orbit or Planar instead.
:::

### Keyboard (Free-fly)

- `W` / `S`: move forward/back
- `A` / `D`: strafe left/right
- `Q` / `E`: move down/up
- `Shift`: sprint (move faster)

### Free-fly options (UI)

In **Compare Views → Navigation → Free-fly**:

- **Look sensitivity**: how fast the view rotates per mouse movement
- **Move speed**: movement speed in “scene units per second”
- **Invert look axes**: flips mouse look direction
- **Projectile shooting** (optional demo feature): enables click/hold-to-charge projectiles
- **Capture pointer**: a checkbox, not a button. Ticking it asks the browser for
  pointer lock straight away; if the browser wants a gesture on the canvas
  first, clicking the canvas asks again. Press `Esc` to release.



---

## Multiview focus rules (the #1 “why did it move the wrong view?” issue)

In **Grid compare** layout, you can see multiple views at once.

Key rule:

- The view you last clicked becomes the **focused view** for camera/navigation.

Practical implications:

- If you scroll/drag and “the wrong panel moves”, click inside the panel you intend to control.
- When **cameras are unlocked**, each view can store its own camera state; clicking a different panel can switch you to its saved camera.

---

## Per-view vs global navigation mode (camera lock matters)

### Cameras locked

- All views share the same camera and navigation mode.
- The navigation dropdown changes behavior globally.

### Cameras unlocked

- Each view stores its own camera *and* navigation mode.
- View badges show a small navigation indicator (Orb / Pan / Fly).
- Clicking the nav indicator cycles a view’s nav mode.

This is powerful for comparisons (e.g., one view in planar 2D, another in orbit 3D), but it’s also a common source of confusion—use the badges to verify what mode each view is in.

---

## Edge cases (things that surprise people)

- **Right-click is disabled / triggers context menu**: Cellucid prevents the default context menu on the canvas, but some browser extensions can interfere.
- **Scroll wheel scrolls the page instead of zooming**: click the canvas once so it has focus, then scroll again.
- **Keyboard controls do nothing**: if a dropdown/slider is focused, keys may be captured; click the canvas (or press `Esc`) to return focus.
- **Pointer lock “won’t enable”**: pointer lock requires a user gesture and can be blocked in iframes or by browser policy.

---

## Troubleshooting (navigation)

For a complete symptom catalog, see {doc}`06_troubleshooting_core_interactions`.

### Symptom: “I can’t rotate/pan/zoom anymore”

**Likely causes:**

- pointer lock is active (free-fly) and you expected orbit/planar behavior,
- you are in the wrong navigation mode,
- focus is on an input element (keyboard shortcuts disabled).

**Fix:**

1) Press `Esc` (exits pointer lock and often “un-sticks” focus).
2) Check **Mode** in Navigation controls.
3) Click the canvas once, then try again.
