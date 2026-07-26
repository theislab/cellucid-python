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
   This layer is *not* inherently accessible to screen readers, because it is an image-like canvas.

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
   In Figure Export, use the **Colorblind preview** dropdown (Normal / Deuteranopia / Protanopia / Tritanopia).  
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
- **Very narrow numeric range**: if everything looks the same color, use “Rescale colorbar to slider range” so small differences are visible.

---

### Categorical fields (clusters, sample, batch)

Cellucid’s default categorical palette is designed to be high-contrast and broadly colorblind-safe, but categorical accessibility still has hard limits:

- **Few categories (≤ ~10–20)**: usually fine with a good legend and sensible labels.
- **Many categories (50, 100, 500…)**: no palette is truly accessible; colors repeat and labels become unreadable.

If you have too many categories:

1) **Collapse to fewer groups** upstream (e.g., coarse cell types instead of clusters).
2) Use categories only as an exploration tool, then export a figure with a smaller subset.
3) Consider using **multiple panels** (snapshots) rather than one legend with 100+ entries.

---

### Colorblind simulation (high value, low effort)

Cellucid includes a **preview-only** colorblind simulation in the Figure Export preview.

Step-by-step:

1) Open the Figure Export panel: {doc}`../k_figure_export/index`
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
- **Viewer background** (e.g., white/black/grid modes)

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

### Global shortcuts (works when focus is not inside a text input)

| Key | Action |
|---|---|
| `?` | Open the in-app “Keyboard Shortcuts” section |
| `h` | Toggle sidebar visibility |
| `f` | Toggle fullscreen |
| `1` / `2` / `3` | Switch dimension (when available) |
| `o` | Switch navigation mode to Orbit |
| `p` | Switch navigation mode to Planar |
| `g` | Switch navigation mode to Free-fly |
| `x` | Clear all highlights |

Important focus rule:
- If you are typing in a search box/dropdown (`input`, `textarea`, `select`), shortcuts are intentionally ignored so you can type normally.

```{figure} ../../../_static/screenshots/web_app/keyboard-shortcuts.png
:alt: Keyboard Shortcuts panel listing camera, dimension, navigation, and highlighting shortcuts.
:width: 100%

The Keyboard Shortcuts panel documents camera, dimension, navigation, and highlighting controls.
```

### What keyboard cannot do (yet)

Even with shortcuts, some core interactions still require a pointer device:
- lasso/area selection
- fine-grained orbit/planar navigation on the canvas
- dragging sliders quickly (HTML sliders also accept Tab plus arrow-key input)

If you need a keyboard-only workflow, treat Cellucid as:
- a viewer for “set state → export figure”, and
- use Python/R notebooks for analyses that require precise selection-by-criteria.

---

## Motion sensitivity and cognitive load (3D can be a lot)

3D navigation, smoke mode, and animated overlays can cause discomfort for some users.

If you (or your audience) is motion-sensitive:

1) Prefer **Planar** navigation over Free-fly.
2) Avoid pointer lock unless you truly need it (exit with `Esc`).
3) Disable or minimize animated overlays (e.g., vector field animations) for presentation recordings.
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
- Turn on **Rescale colorbar to slider range** (see {doc}`../d_fields_coloring_legends/03_color_by_behavior`).

---

### Symptom: “Keyboard shortcuts don’t work”

**Likely causes (ordered):**
1) Keyboard focus is inside a text input (search box, dropdown).
2) A modal dialog is open (welcome, export dialog, etc.).
3) The app is embedded in an iframe that intercepts focus (notebook environments).

**How to confirm**
- Click an empty area of the page (or the canvas), then press `?`.

**Fix**
1) Click the canvas, then try again.
2) Close modals with `Esc` (also exits pointer lock).
3) If embedded in a notebook, open the viewer in a new browser tab (see {doc}`../b_data_loading/05_jupyter_tutorial`).

---

### Symptom: “Free-fly makes me dizzy / I get lost”

**Fix**
- Press `o` (Orbit) or `p` (Planar).
- Use Reset Camera (see {doc}`../c_core_interactions/02_camera_controls_advanced`).
- Reduce motion: avoid pointer lock; prefer small, deliberate movements.

---

## Next steps

- If you are exporting figures: {doc}`../k_figure_export/index`
- If you are choosing palettes/ranges: {doc}`../d_fields_coloring_legends/03_color_by_behavior`
- If you are teaching workflows: {doc}`../a_orientation/03_quick_tour_60_seconds`
