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

Switch modes in **Navigation → Mode** (or use shortcuts `O`/`P`/`G`). The mode
belongs to the focused view; see
[Each view navigates on its own](#each-view-navigates-on-its-own).

---

## Where the controls live in the UI

Navigation mode controls are in the left sidebar:

- **Navigation** accordion → **Mode** dropdown

The mode options are:

- Orbit
- Planar
- Free-fly (shown as “Free-fly” in the UI, internally stored as `free`)

Keyboard shortcuts (global, unless you’re typing in an input):

- `O` → Orbit
- `P` → Planar
- `G` → Free-fly

You do not have to set the mode at all: the dimension of the embedding a view
shows picks one for it — **Orbit** in 3D, **Planar** in 1D and 2D — every time
the view switches to another dimension. A mode you choose yourself is
remembered for that view and that dimension. See
{doc}`05_dimension_switching_1d_2d_3d`.

The **Mode** dropdown is followed by a box of options that belong to the chosen
mode only, so the panel looks different in each of the three. A still cannot
show a camera moving — but it can show you exactly which knobs you have, and
that is what people actually get lost looking for.

```{figure} ../../../_static/screenshots/web_app/navigation-controls-orbit.png
:alt: The Navigation block with Mode set to Orbit and the pointer pressing the select: a Keyboard speed slider, Google Earth-style drag and Show orbit anchor.
:width: 488px

**Orbit** — a keyboard speed slider and two checkboxes.
```

```{figure} ../../../_static/screenshots/web_app/navigation-controls-planar.png
:alt: The Navigation block with Mode set to Planar: a Keyboard speed slider, Zoom to cursor and Invert axes.
:width: 488px

**Planar** — a keyboard speed slider and two different checkboxes. Orbit and
Planar both label their slider `Keyboard speed`; the value beside it is on a
different scale in each, because one scales rotation and the other pan.
```

```{figure} ../../../_static/screenshots/web_app/navigation-controls-freefly.png
:alt: The Navigation block with Mode set to Free-fly: Look sensitivity and Move speed sliders with their readouts, Invert look axes and Projectile shooting, and Capture pointer with its (i).
:width: 488px

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
- **2 finger drag**: pan (move both fingers together)

With **Select mode** on (its switch is on the selection bar at the bottom of
the plot and in the Highlighting panel), one finger runs the active selection
tool instead of rotating; two fingers still pinch and pan. To rotate between
two selection steps, turn the switch off, rotate, and turn it back on: the
selection in progress stays on the bar, with its **Confirm** and **Cancel**.

### Orbit options (UI)

In **Navigation → Orbit**:

- **Keyboard speed**: scales keyboard rotation/zoom speed.
- **Google Earth-style drag** (default on): flips the horizontal drag sign (matches “grab-and-rotate” muscle memory for many users).
- **Show orbit anchor**: shows a compass/anchor overlay in orbit mode.



---

## Planar mode (2D “map” navigation)

Planar mode is designed for 2D embeddings: you pan and zoom like a map. Rotation is disabled.

### Mouse / trackpad (Planar)

- **Click + drag**: pan (always)
- **Scroll wheel / trackpad scroll**: zoom

If **Zoom to cursor** is enabled (default), the zoom centers around your cursor location rather than the view center.

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
- **2 finger drag**: pan

With **Select mode** on, one finger selects instead of panning; pan with two
fingers.

### Planar options (UI)

In **Navigation → Planar**:

- **Keyboard speed**: scales keyboard pan/zoom speed.
- **Zoom to cursor** (default on): zoom around cursor.
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
Free-fly has no touch navigation: a finger neither looks around nor moves the camera in free-fly. If you’re on a touch-only device, use Orbit or Planar instead. **Select mode** still selects with one finger in free-fly, as `Alt`+click does with a mouse.
:::

### Keyboard (Free-fly)

- `W` / `S`: move forward/back
- `A` / `D`: strafe left/right
- `Q` / `E`: move down/up
- `Shift`: sprint (move faster)

### Free-fly options (UI)

In **Navigation → Free-fly**:

- **Look sensitivity**: how fast the view rotates per mouse movement
- **Move speed**: movement speed in “scene units per second”
- **Invert look axes**: flips mouse look direction
- **Projectile shooting** (optional demo feature): enables click/hold-to-charge projectiles
  that bounce off cells and off the backdrop: the walls of the **Grid**
  background, or the floor of a **Studio** background under a 3D view. With
  neither, or past the edge of the Studio floor, they fall away.
- **Capture pointer**: a checkbox, not a button. Ticking it asks the browser for
  pointer lock straight away; if the browser wants a gesture on the canvas
  first, clicking the canvas asks again. Press `Esc` to release.



---

## Multiview focus rules (the #1 “why did it move the wrong view?” issue)

In **Grid compare** layout, you can see multiple views at once.

Key rule:

- The view you last clicked, tapped, or scrolled over becomes the **focused view** for camera/navigation.

Practical implications:

- Scrolling (wheel zoom) over a panel focuses that panel and zooms it; with **Zoom to cursor** the point under the cursor in that panel stays put.
- If you drag and “the wrong panel moves”, click inside the panel you intend to control.
- Focusing a panel loads its own camera and navigation mode, so the **Mode** dropdown changes to that view's mode.

---

## Each view navigates on its own

Every view keeps its own camera and its own navigation mode:

- A view showing a 3D embedding starts in **Orbit**; a view showing a 1D or 2D
  embedding starts in **Planar**. A 3D view beside a 2D view keeps orbiting
  while the 2D view pans.
- When a view switches to an embedding of another dimension, it takes that
  dimension's mode — unless you chose a mode for that view in that dimension
  before, which it takes instead. A view you set to **Free-fly** in 3D flies
  again whenever it shows a 3D embedding, and pans whenever it shows a 2D one.
- The **Mode** dropdown, its options and the `O`, `P` and `G` keys edit the
  focused view. With several views open, the **Applies to** row above the
  dropdown names it.
- With **Link cameras** on, views that show the same embedding share one camera
  and one mode, so a change reaches all of them, and **Applies to** names them
  ("with view 2"). Views of different embeddings are never linked. See
  {doc}`04_view_layout_live_snapshots_small_multiples`.
- Each view's row in **Compare views** shows its mode (`Orbit`, `Pan` or
  `Fly`); press it to switch that view's mode without focusing it.

A session saves every view's camera and mode.

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
