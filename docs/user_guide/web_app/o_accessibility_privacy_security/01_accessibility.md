# Accessibility

**Audience:** everyone publishing plots + workshop instructors + computational users making reusable figures  
**Time:** 15–40 minutes  
**What you’ll learn:**
- Practical defaults for color/contrast that work for most audiences
- How to check and improve colorblind accessibility *in Cellucid* (without guesswork)
- What keyboard-only navigation can and cannot do today
- How to reduce motion/cognitive load in a dense 3D visualization

---

## What “accessibility” means for a WebGL point cloud

Cellucid has two layers:

1) **The UI layer (HTML)**: sidebar controls, legends, buttons, search boxes, accordions.  
   This layer can be keyboard navigable and screen-reader friendly (when built carefully).
2) **The viewer layer (WebGL canvas)**: the actual millions-of-points plot.  
   A canvas cannot expose individual points to a screen reader — there is no
   element per cell to announce.

What the canvas *does* expose is itself. It is the **first `Tab` stop on the
page**, it sits in a labeled group, and it carries an accessible name and an
accessible description that names its movement keys — so a screen-reader user
reaches the plot first and is told how to drive it, rather than landing on an
unnamed graphic.

So the goal is not “screen reader can read every point”, but:

- readers can **distinguish** what you mean (color/contrast choices),
- you can provide **text alternatives** (legends, labels, captions, methods),
- and the UI is usable with **keyboard and assistive tech** as much as the browser allows.

:::{important}
If you need WCAG-level accessibility for downstream publishing, treat Cellucid as an *exploration tool* plus a *figure generator*.
Do the final accessibility pass on exported figures (color checks, captions, alt text, font size).
:::

---

## Fast path: make an accessible figure in ~10 minutes

If you don’t want to read the full page, do this checklist before exporting a figure:

1) **Prefer perceptually uniform continuous colormaps**  
   Use **Viridis** or **Cividis** when coloring by expression/QC (they remain readable under color vision deficiency and grayscale).
2) **Avoid “rainbow/jet” for quantitative interpretation**  
   Rainbow palettes can create artificial boundaries that look like biology.
3) **Check colorblind simulation before you export**  
   In **Figure export**, use the **Colorblind preview** dropdown (Normal / Deuteranopia / Protanopia / Tritanopia).  
   This catches “red vs green” failure modes quickly.
4) **Choose a background that preserves contrast**  
   If your points look washed out, switch viewer background or theme (dark vs light) and re-check.
5) **Don’t rely on color alone for categories**  
   If categories matter, use a clear legend, meaningful labels, and (when possible) reduce category count by grouping.
6) **Export at readable size**  
   Small PNGs with tiny text are inaccessible. Prefer higher DPI or larger plot size (see {doc}`../k_figure_export/04_quality_knobs_and_best_practices`).

---

## Color accessibility (the biggest win)

### Continuous fields (gene expression, QC, scores)

**Recommended defaults**
- Colormap: **Viridis** (default in Cellucid) or **Cividis**
- Log scale: use only if you can explain it in a caption/methods
- Handle missing values deliberately: Cellucid uses a neutral gray “None” for missing/invalid values

Why this matters:
- Viridis/Cividis are designed so equal numeric steps look like equal visual steps.
- They are significantly more robust under common color vision deficiencies.

Related: {doc}`../d_fields_coloring_legends/03_color_by_behavior` (exact behavior of log scale, “None” gray, and rescaling).

**Edge cases to watch**
- **All zeros + log scale**: if you turn on log scale for a gene that is zero in most cells, many points can become “None” gray (log undefined at 0). That can look like “missing data”.
- **Very narrow numeric range**: if everything looks the same color, use “Rescale to range” so small differences are visible.

---

### Categorical fields (clusters, sample, batch)

A field whose file stores Scanpy's colors for it (`uns['<field>_colors']`)
opens in those colors, whatever they are; see
{doc}`../d_fields_coloring_legends/03_color_by_behavior`. Otherwise Cellucid
assigns categorical colors in three stages, and **none of them is chosen for
color vision deficiency**. This is the most useful thing to know on this page,
so it is worth being precise:

1) **The first 24 categories** come from a vivid palette built from Okabe-Ito,
   IBM Design and Paul Tol colors. It is tuned for distinct hues to normal
   vision, not for color blindness: simulated as a red–green (deuteranopia or
   protanopia) reader sees it (Machado et al., 2009), several of its colors
   look alike — vivid cerulean and iris purple (categories 2 and 10), the
   blues and violets among themselves, bright orange and amber (1 and 16) —
   and already the second and third colors (cerulean and violet) are close.
2) **Beyond 24**, colors are drawn from a general-purpose picker palette.
3) **Beyond that**, colors are generated by walking the color wheel at the
   golden angle. They are guaranteed not to repeat quickly, and nothing more.

So, if your figure must be readable to a colorblind reader:

- Check it in the export preview's **Colorblind preview** (below), which
  simulates the same three deficiencies with the same model.
- Change the colors of the categories that get confused (click a category's
  swatch in the legend), or store a palette chosen for color blindness in the
  file (`uns['<field>_colors']`, which Cellucid uses as it is).
- Keep the categories you rely on few, and label them.

Beyond that, categorical accessibility has limits no palette can fix:

- **Many categories (50, 100, 500…)**: colors become indistinguishable and
  labels become unreadable regardless of how they were chosen.

If you have too many categories:

1) **Collapse to fewer groups** upstream (e.g., coarse cell types instead of clusters).
2) Use categories only as an exploration tool, then export a figure with a smaller subset.
3) Consider using **multiple panels** (snapshots) rather than one legend with 100+ entries.

---

### Colorblind simulation (high value, low effort)

Cellucid includes a **preview-only** colorblind simulation in the Figure Export preview: each option is
the full-severity (dichromat) model of Machado, Oliveira and Fernandes (2009), applied in linear light.

Step-by-step:

1) Open the **Figure export** panel: {doc}`../k_figure_export/index`
2) Click **Preview** (or open the preview section)
3) Find **Colorblind preview**
4) Toggle through:
   - Normal
   - Deuteranopia
   - Protanopia
   - Tritanopia
5) If the plot becomes ambiguous, change your colormap/palette or simplify the categories.

:::{note}
The simulation is meant to catch *obvious* issues quickly. It is not a substitute for a formal accessibility audit of a final figure (especially for print).
:::



---

## Contrast, theme, and background (low vision + print robustness)

Cellucid lets you adjust:
- **Theme** (light/dark UI)
- **Viewer background** (plain white or black, grid, studio, or a painted look
  such as Atelier or Nebula). Every look keeps its large areas outside the
  lightness range of the category palettes, and with the system's *reduce
  motion* setting a look appears without its fade.

Until you pick a background, it follows the theme: switching to dark gives the
dark twin of the current background (Studio (light) becomes Studio (dark),
Atelier becomes Nebula), and back. A background you pick stays when the theme
changes.

Practical guidance:

- If you plan to publish on a white page, preview your figure on a **white background** at least once.
- If you are exploring in dark mode, double-check that your colors remain distinguishable on light backgrounds (and vice versa).

Common failure mode:
- A “pretty” palette on a dark background can become unreadable on white, especially for lighter categories.

---

## Keyboard accessibility (what works today)

Cellucid has global keyboard shortcuts intended to make basic navigation and “getting unstuck” easier.

### Compact `i` help dialogs

Small `i` buttons provide explanatory guidance without hiding the nearby
action label or current-state text.

- Click an `i`, or focus it with {kbd}`Tab` and press {kbd}`Enter` or
  {kbd}`Space`, to open its dialog. Focus moves into the dialog.
- {kbd}`Tab` reaches any links, then closes the dialog and continues to the
  control after the `i`.
- {kbd}`Shift+Tab` or {kbd}`Escape` closes the dialog and returns focus to the
  `i`. Clicking or moving focus elsewhere closes it without moving focus.
- Only one help dialog can be open.

For screen readers, every `i` is a named button with
`aria-controls`, `aria-haspopup="dialog"`, and an accurate `aria-expanded`
state. Its controlled panel is a labeled `role="dialog"`, so the button,
relationship, and open state are announced without removing the visible
control label.

### Shortcuts

{doc}`03_keyboard_shortcuts` lists every key Cellucid answers, from the same
table the in-app **Keyboard shortcuts** panel ({kbd}`?`) and the keys
themselves are made from: the page's shortcuts ({kbd}`?`, {kbd}`/`,
{kbd}`H`, {kbd}`F`, {kbd}`R`, {kbd}`1` {kbd}`2` {kbd}`3`, {kbd}`O`
{kbd}`P` {kbd}`G`, {kbd}`X`, undo and redo), the camera keys, the
selection gestures, the History list's keys and the one way of moving a
panel, a page or a category by keys.

Navigation is **not** pointer-only: {kbd}`W` {kbd}`A` {kbd}`S` {kbd}`D`
orbit, pan or move the camera (one meaning per mode), {kbd}`=` {kbd}`]` and
{kbd}`-` {kbd}`[` zoom in Orbit and Planar, {kbd}`E` {kbd}`Q` move up and
down in Free-fly, which has no keyboard zoom (fly forward with {kbd}`W`).

`Esc` cancels a selection in progress, exactly as the **Cancel** button does,
and with nothing in progress leaves **Select mode**. It is only consumed then,
so the same key still closes a dialog or a menu first.

Who a key belongs to:

- A key is the page's unless what has focus uses it itself: a text field types
  it (and keeps its own undo), a dropdown jumps to the option a typed letter
  starts, and a menu, a list or a dialog keeps the keys pressed in it (the
  welcome card included: `Escape` closes it).
- A slider, a checkbox or a button leaves the letters to the page: after moving
  **Keyboard speed** or ticking a box, {kbd}`W` moves the camera at once. Undo
  and redo work from them too, since they have no undo of their own.

```{figure} ../../../_static/screenshots/web_app/keyboard-shortcuts.png
:alt: The Keyboard shortcuts panel open in the Cellucid sidebar, listing its groups, each with its (i): General (R, F, H, 1 2 3, O P G, ?, /, Ctrl+Z to undo the last step, Shift+Ctrl+Z or Ctrl+Y to redo), Navigation (W A S D and Shift), Zoom (= ] and - [), Free-fly (E Q and Esc), Selection (Alt+click, Alt+drag, Shift+Alt+click, Ctrl+Alt+click, Esc to cancel a selection, Esc to leave Select mode, and X to clear highlights), History list and Moving by keys.
:width: 532px

Press `?` to open this panel in the app. It shows **seven** groups:
**General**, **Navigation**, **Zoom**, **Free-fly**, **Selection**, **History
list** and **Moving by keys**, the same rows as
{doc}`03_keyboard_shortcuts`, both made from one table.
```

### What keyboard cannot do (yet)

The genuine gap is **selection**, not navigation:

- **Lasso and area selection require a pointer.** There is no keyboard equivalent
  for dragging a region over the plot. The pointer does not need a held key,
  though: with **Select mode** on (Highlighting panel or the selection bar on
  the plot), a plain click, tap, or drag runs the active tool, and its
  **Intersect / Add / Subtract** buttons replace `Shift` and `Ctrl`. That serves touch screens and anyone who cannot
  hold a modifier key while pressing.
- **Clicking an individual cell requires a pointer.** Picking is a
  point-and-shoot operation against the drawn scene.
- Dragging a slider quickly is easier with a pointer, though every slider accepts
  `Tab` plus arrow keys and paints a visible focus ring when reached that way.

If you need a keyboard-only workflow:

- navigation, switching embeddings (the **View N shows** picker and keys
  `1`/`2`/`3`), selecting a view (each row in **Compare views** is a button),
  field selection, filtering, page/highlight
  management, docking and the analysis panels are all operable from the
  keyboard, so “set state → export figure” works end to end;
- but use Python/R notebooks for anything that requires selecting cells by
  pointing at them.

### Assistive-technology support that is actually verified

These are asserted by the browser test suite rather than by inspection, so they
are safe to rely on:

- The canvas has an accessible name and a description naming its movement keys,
  and is the first `Tab` stop.
- Every icon-only button resolves to exactly one accessible name; the sidebar
  toggle's name and its expanded state change together.
- While the sidebar is hidden it is `inert`: none of its controls can take
  focus, be reached with `Tab`, or be read by a screen reader until the sidebar
  is shown again.
- No panel header nests a control inside its open/close toggle. A panel's
  undock, dock, copy and close buttons sit beside the header rather than in it,
  are named after their panel (“Undock Visualization panel into a floating
  window”), and come just before the panel's header in `Tab` order. The header
  opens and closes the panel with `Enter` or `Space`; activating one of its
  buttons never does, and after docking or undocking, focus moves to the button
  that reverses the move. While a section floats, the sidebar's **Dock all**
  button (named “Dock all 3 panels back into the sidebar”) docks every one;
  focus then moves to the first of them, and the move is announced.
- Both **Select mode** switches (the Highlighting panel's and the selection
  bar's) are switches (`role="switch"`) named “Select mode” whose
  `aria-checked` state always matches, and each Intersect / Add / Subtract
  choice is a named group with one pressed button, the same one on both. The
  selection bar is a region named “Selection”; its cell count and its note on
  a gesture that changed nothing are announced politely, and every control on
  it is reachable and operable from the keyboard.
- Every control in every analysis panel has an accessible name, and highlight
  page chips are operable with both `Enter` and `Space`, keep focus after
  activation, and announce their state (“120 cells”, “no cells, cannot be
  selected”).
- Every dialog (the welcome, the license, a confirmation, Cite, Share, the
  history's dialogs, the community annotation dialogs) is the browser's own
  modal dialog: it takes focus, the page behind it is out of reach for the
  mouse, the keyboard and a screen reader, and the keys pressed in it never
  reach the page's shortcuts. `Escape` or a click beside the card closes it
  and focus goes back where it was; a dialog opened over another gives focus
  back into the one under it. `Tab` moves through the dialog's controls and,
  past the last, to the browser's own controls and back in.
- Every menu (**Open…**, **My datasets**, the History panel's ⋯, a history
  row's or branch's menu, the history search's **Filter**, a highlight page's
  **Combine**) moves with the arrow keys, `Home` and `End`; `Escape` closes it
  and returns to its button, and `Tab` closes it and moves on from its button
  (`Shift+Tab`, back from it). A menu button states whether its menu is open.
- A color swatch (legend, highlight pages, **Create categorical**, an
  analysis's page chips) draws the focus ring when its color input has
  keyboard focus.
- The app respects the operating system's **reduce motion** setting: interface
  transitions and animations stop, cinematic camera playback and the History's
  guided tour jump to each camera instead of flying, and nothing in the viewer moves
  by itself (see [Motion sensitivity](#motion-sensitivity-and-cognitive-load-3d-can-be-a-lot)
  below). A change of the setting applies at once.
- Text contrast meets WCAG AA (4.5:1) in both themes: every text color token
  is measured against every text surface (`theme-contrast-contract`), and
  every line of sidebar copy, the open colormap menu included, the welcome
  card, the dialogs and the license, is measured as rendered against the
  surface it sits on, in light and in dark.
- Opening a dataset announces its progress: the opening step is a polite
  status region, and notifications are announced the same way.

What is **not** verified, and should not be assumed:

- there is no automated WCAG audit in the test suite, and no formal audit has
  been made, so Cellucid does not claim WCAG conformance;
- rendering progress (a long draw or a smoke build) is not announced.

The site's short accessibility statement, with how to report a barrier, is
[www.cellucid.com/accessibility.html](https://www.cellucid.com/accessibility.html),
linked as **Accessibility** at the foot of the sidebar.

---

## Motion sensitivity and cognitive load (3D can be a lot)

3D navigation, smoke mode, and animated overlays can cause discomfort for some users.

**Start with the operating system setting.** Cellucid honors
`prefers-reduced-motion`: turn on “Reduce motion” (macOS), “Show animations”
off (Windows), or your desktop's equivalent, and nothing moves unless you move
it:

- interface transitions and animations stop;
- cinematic camera playback and the History's guided tour jump straight to each camera
  instead of flying between them;
- background fades complete at once and painted backdrop looks hold still;
- the **vector field overlay** draws a still picture of the flow instead of
  moving particles: streaks traced once along the field, their direction
  (a bright head, a fading tail) and speed (length and color) the same as
  the moving flow shows, the same picture every time
  ({doc}`../i_vector_field_velocity/05_performance_and_quality`);
- highlighted cells keep a steady ring instead of pulsing;
- **Volumetric smoke** stops drifting, whatever its `Animation speed`;
- the camera moves only while you move it: a drag let go fast stops where you
  let go instead of gliding on, a zoom (wheel, keys) arrives at once, and
  free-fly stops when you release its keys.

A still scene then draws nothing until you change something. Switching the
setting while Cellucid is open applies at once, both ways. This costs you
nothing and needs no in-app configuration.

An agent working in Cellucid knows the setting is on: reading the views
(`views_state`) says `reduced_motion`, and a camera path it plays answers that
it went to its end at once (`ended: "reduced-motion"`, 0 seconds of motion),
so it does not tell you a film played.

What you drive yourself still moves while you drive it — orbiting, panning,
zooming, free-fly. For those:

1) Prefer **Planar** navigation over Free-fly.
2) Avoid pointer lock unless you truly need it (exit with `Esc`).
3) Without the system setting, turn the vector field overlay off (or set smoke's `Animation speed` to `Still`) when the motion distracts.
4) Use snapshots (“Keep view”) to build a story as static panels instead of live spinning 3D.

Related:
- {doc}`../c_core_interactions/01_navigation_modes_orbit_planar_free_fly`
- {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke`
- {doc}`../i_vector_field_velocity/index`

---

## Troubleshooting (accessibility)

### Symptom: “I can’t tell categories apart (especially in print)”

**Likely causes (ordered):**
1) Too many categories for a readable legend.
2) Background/theme reduces contrast.
3) Categories include many light colors that wash out on white.

**How to confirm**
- Export a small test PNG and view it at the size it will appear in a paper slide/doc.
- Check the same plot under colorblind simulation (Figure Export preview).

**Fix**
1) Collapse categories (coarser grouping) or export only a subset.
2) Switch viewer/export background to white and re-check.
3) Increase point size and legend font size in the export settings.

**Prevention**
- Design figures with ≤ ~10–20 categorical colors per panel when possible.

---

### Symptom: “Gene expression looks ‘flat’ (everything the same color)”

**Likely causes:**
- The numeric range is very narrow after filtering.
- The color domain is not rescaled to your current slider range.

**Fix**
- Turn on **Rescale to range** (see {doc}`../d_fields_coloring_legends/03_color_by_behavior`).

---

### Symptom: “Keyboard shortcuts don’t work”

**Likely causes (ordered):**
1) The welcome dialog is still up. The keys pressed in a dialog stay in it;
   `Escape` closes it.
2) Keyboard focus is inside a text field or a dropdown (a dropdown's letters
   choose its options). A slider, a checkbox or a button does not hold the keys.
3) Another modal dialog is open (export dialog, voting dialog, etc.).
4) The app is embedded in an iframe that intercepts focus (notebook environments).

**How to confirm**
- Press `Escape` first, then press `?`.
- If that works, a dialog was open: its keys stay in it.

**Fix**
1) Press `Esc` to close the dialog (this also exits pointer lock).
2) Click the canvas or an empty area of the page, then try again.
3) If embedded in a notebook, open the viewer in a new browser tab with
   `cellucid.AnnDataViewer(...)` (see {doc}`../b_data_loading/05_jupyter_tutorial`).

---

### Symptom: “Free-fly makes me dizzy / I get lost”

**Fix**
- Press `o` (Orbit) or `p` (Planar).
- Use Reset camera (see {doc}`../c_core_interactions/02_camera_controls_advanced`).
- Reduce motion: avoid pointer lock; prefer small, deliberate movements.

---

## Next steps

- If you are exporting figures: {doc}`../k_figure_export/index`
- If you are choosing palettes/ranges: {doc}`../d_fields_coloring_legends/03_color_by_behavior`
- If you are teaching workflows: {doc}`../a_orientation/03_quick_tour_60_seconds`
