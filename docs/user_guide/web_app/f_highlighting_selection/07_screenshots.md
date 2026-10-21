# Verified highlighting captures

Every image below was captured from the running web app against the
**Pancreatic endocrinogenesis (scVelo)** sample (3,696 cells), coloured by
`cell_type`, at a 1440×1000 viewport with a device pixel ratio of 2, in the
light theme. Each one is a named scenario in `docs/_tooling/screenshots/`
(`highlighting-…` in `scenarios.mjs`) and can be regenerated and re-verified
byte for byte with `node capture.mjs check`. The selection gestures are real
gestures: `Alt` held while the mouse presses, drags and releases over the plot.
The pointer glyph marks where the pointer was and whether the button was down,
and every gesture is checked against the step readout the panel wrote for it,
so nothing here is a mock-up.

The 2D embedding is shown in **Planar** navigation because a flat projection
makes the geometry of a gesture legible. Everything on this page behaves the
same way in 1D and 3D; the 3D-specific caveats are in
{doc}`02_selection_tools_document_each_tool`, and the two-panel capture below is
taken in 3D.

## The toolbelt

```{figure} ../../../_static/screenshots/highlighting_selection/choose-highlight-mode.png
:alt: The Highlighting panel before anything is selected: the four Highlight mode buttons — Annotation based (active), KNN drag, Proximity drag and Lasso — ringed, with the pointer on Lasso; below a dashed divider the Select mode switch, off, with the hint "tap, click or drag selects without Alt"; the Annotation based help text; one page tab reading Page 1 (0); and the empty group list reading "No cells highlighted".
:width: 516px

The **Highlighting** accordion before anything is selected. Four modes; below a
divider the **Select mode** switch, which stands in for a held `Alt` on a touch
screen ({ref}`select-mode`); a description of the active mode; one page tab
reading `Page 1 (0)`; and `No cells highlighted`. **Clear** is not drawn at all
while the page is empty.
```

## One lasso, start to finish

```{figure} ../../../_static/screenshots/highlighting_selection/01-draw-lasso.png
:alt: The 2D planar Pancreas embedding mid-lasso: Alt is held, the pointer is pressed at the end of a dotted path drawn most of the way round the upper part of the left lobe, and the cells the path already encloses are drawn in the highlight colour while the step readout has not yet appeared.
:width: 1440px

**Step 1 — draw.** `Alt` is held and the button is down. The enclosed cells are
already repainted in the highlight colour: that is the *preview*, which shows
the set a release would commit under the modifier you are holding. The path
does not have to be closed by hand — releasing closes it.
```

```{figure} ../../../_static/screenshots/highlighting_selection/02-review-lasso-step.png
:alt: The Highlighting panel after the lasso was released, showing the step readout "Step 1: 842 cells" with the help line "Alt to intersect, Shift+Alt to add, Ctrl+Alt to subtract", and the ringed Confirm, undo, redo and Cancel step controls with the pointer on Confirm; the page tab still reads Page 1 (0) and the list still reads "No cells highlighted".
:width: 516px

**Step 2 — read the step.** Releasing over the plot records a step. The readout
names the step number and the candidate count; the step controls appear beneath
it. Nothing is stored yet — `Page 1 (0)` and `No cells highlighted` are still
true.
```

```{figure} ../../../_static/screenshots/highlighting_selection/03-confirm-lasso-group.png
:alt: The Highlighting panel after Confirm: the page tab reads Page 1 842, and the ringed group list shows the summary "842 cells highlighted", a Clear button, and one group row labelled "Lasso (842 cells) (842)" with a checkbox and a remove control.
:width: 516px

**Step 3 — confirm.** **Confirm** turns the candidate set into one group on the
active page. The page tab, the summary line and the group row now all carry the
same count, and **Clear** has appeared.
```

```{figure} ../../../_static/screenshots/highlighting_selection/04-view-highlighted-cells.png
:alt: The whole window after confirming the lasso group: the 842 enclosed cells in the upper part of the left lobe are drawn in the highlight colour, and the Highlighting panel reports 842 cells highlighted with one Lasso group.
:width: 1440px

**Step 4 — the result.** The group survives everything that is not a filter:
change the coloured field, move the camera, switch embedding, and these 842
cells stay members.
```

## The other three tools

```{figure} ../../../_static/screenshots/highlighting_selection/select-category-by-click.png
:alt: Annotation based selection after Alt+click on one cell at the tip of the green Beta arm: all 591 Beta cells are drawn in the highlight colour, and the panel reads "Step 1: 591 cells" with the help line "Alt to replace, Shift+Alt to add, Ctrl+Alt to subtract" above the Confirm, undo, redo and Cancel controls.
:width: 1440px

**Annotation based.** One `Alt+click` on one cell selected all 591 visible cells
of that cell's `cell_type` category, `Beta` — including the ones far from the
pointer. The help line changes to *Alt to replace* on a categorical field,
because intersecting one category with another is almost always empty.
```

```{figure} ../../../_static/screenshots/highlighting_selection/drag-proximity-radius.png
:alt: Proximity drag mid-gesture on the 2D embedding: a dashed circle centred on a dot at the seed cell in the left lobe, the cells inside it drawn in the highlight colour, the pointer pressed on the circle's edge, and the panel still showing the Proximity drag help text.
:width: 1440px

**Proximity drag.** The dashed circle is the radius in *embedding* space, not
screen space; the small dot at its centre is the seeded cell. In 3D the same
gesture is a sphere. The step readout appears only on release, so mid-drag the
panel still shows the tool's ordinary help text.
```

```{figure} ../../../_static/screenshots/highlighting_selection/drag-knn-degree.png
:alt: KNN drag mid-gesture from the same seed cell, dragged upward by the same distance: a small dashed circle labelled 3° around the seed, and the highlighted cells spreading through the left lobe along the neighbour graph rather than filling a disc, with the pointer pressed above the seed.
:width: 1440px

**KNN drag.** The same seed cell and the same drag distance as the proximity
capture above, dragged upward. Compare the shape: the selection follows the
neighbour graph, so it spreads along the manifold instead of filling a disc.
Drag distance sets the degree in whole steps, and the label beside the dashed
circle (`3°`) names the degree reached.
```

## Two pages

```{figure} ../../../_static/screenshots/highlighting_selection/compare-highlight-pages.png
:alt: The Highlighting panel with two page tabs ringed — Page 1 with a blue swatch carrying 842, and the active Page 2 with a red swatch carrying 223, each with a combine (⋈) and a delete (×) control — the group list showing one "Proximity (223 cells) (223)" group, and the pointer on the button that adds another page.
:width: 516px

Pages are alternative workspaces. Each carries its own colour, its own count and
its own group list; only the active one is rendered. The `⋈` control on a tab
combines it with another page by intersection or union.
```

```{figure} ../../../_static/screenshots/highlighting_selection/share-highlight-across-views.png
:alt: Two 3D panels of the same dataset side by side with the sidebar collapsed: the left panel at the default close-up camera, the right panel zoomed out and rotated; the cells confirmed by one lasso drawn in the left panel are drawn in the highlight colour in both.
:width: 1440px

Highlight membership is one list for the whole dataset, so both panels draw the
same members. What each panel draws is still its own business: a cell a panel is
filtering out is not highlighted there. See
{doc}`04_selection_synchronization`.
```

## What the app says when a gesture fails

```{figure} ../../../_static/screenshots/highlighting_selection/read-missed-click-notice.png
:alt: A proximity Alt+click on empty space: the pointer is pressed above the embedding where no cell is drawn, and the ringed line in the Highlighting panel adds "That click did not land on a cell, so there is nothing to select. Zoom in or raise Point size to make cells easier to hit: only the cells the viewer is currently drawing can be clicked." under the Proximity drag help text.
:width: 1440px

A click that hits nothing says so, in one line, under the **Highlight mode**
buttons. That line is the only place these messages appear — no toast, no
history. The full list is in {doc}`06_troubleshooting_highlighting`.
```

```{figure} ../../../_static/screenshots/highlighting_selection/read-no-active-field-notice.png
:alt: The Highlighting panel in Annotation based mode with no coloured field, with the ringed line "Annotation selection needs an active field. Choose a Categorical obs or Continuous obs field under Coloring & Filtering, then Alt+click a cell (or tap it in Select mode)." in place of the usual help text.
:width: 516px

With no coloured field, **Annotation based** has no rule to apply, so it replaces
its usual invitation with an instruction instead of failing silently.
```

```{figure} ../../../_static/screenshots/highlighting_selection/read-annotation-step.png
:alt: The Highlighting panel after an annotation Alt+click on a Beta cell, with the ringed step readout "Step 1: 591 cells" and the categorical help line "Alt to replace, Shift+Alt to add, Ctrl+Alt to subtract".
:width: 516px

The help line under the step readout is tool-specific *and* field-specific. Read
it before you press a modifier: it is the only place the app tells you that
`Alt` means **replace** here and **intersect** everywhere else.
```

## A named page

```{figure} ../../../_static/screenshots/highlighting_selection/highlighting-selected-page.png
:alt: The Highlighting panel with one page renamed "Left lobe" carrying 842 cells, the summary "842 cells highlighted", and one group row labelled "Lasso (842 cells) (842)".
:width: 516px

The confirmed lasso group from above, on a page renamed by double-clicking its
tab. Highlighted cells are stored on a named page, and that page is the explicit
group later analysis and session workflows use. The Python guide reuses this
capture.
```
