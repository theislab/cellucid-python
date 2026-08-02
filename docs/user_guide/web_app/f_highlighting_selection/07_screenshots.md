# Verified highlighting captures

Every image below was captured from the running web app against the
**Pancreatic endocrinogenesis (scVelo)** sample (3,696 cells), at a 1440×1000
viewport with a device pixel ratio of 2, in the light theme. The selection
gestures are real gestures: the pointer glyph marks where the pointer was and
whether the button was down, and each coordinate was checked against the
viewer's own picking before the gesture ran, so nothing here is a mock-up.

The 2D embedding is shown in **Planar** navigation because a flat projection
makes the geometry of a gesture legible. Everything on this page behaves the
same way in 1D and 3D; the 3D-specific caveats are in
{doc}`02_selection_tools_document_each_tool`.

## The toolbelt

```{figure} ../../../_static/screenshots/highlighting_selection/choose-highlight-mode.png
:alt: The Highlighting panel with the four Highlight mode buttons — Annotation based, KNN drag, Proximity drag, Lasso — the pointer on Lasso, and the empty group list reading "No cells highlighted".
:width: 516px

The **Highlighting** accordion before anything is selected. Four modes, a
one-line description of the active mode, one page tab reading `Page 1 (0)`, and
`No cells highlighted`. **Clear** is not drawn at all while the page is empty.
```

## One lasso, start to finish

```{figure} ../../../_static/screenshots/highlighting_selection/01-draw-lasso.png
:alt: The 2D planar Pancreas embedding mid-lasso: Alt is held, the pointer is pressed, and a freehand path encloses the Ductal lobe on the left while the enclosed cells are already drawn in the highlight colour.
:width: 1440px

**Step 1 — draw.** `Alt` is held and the button is down. The enclosed cells are
already repainted in the highlight colour: that is the *preview*, which shows
the set a release would commit under the modifier you are holding.
```

```{figure} ../../../_static/screenshots/highlighting_selection/02-review-lasso-step.png
:alt: The Highlighting panel after the lasso was released, showing the step readout "Step 1: 870 cells" and the Confirm, undo, redo and Cancel step controls, with the pointer on Confirm.
:width: 516px

**Step 2 — read the step.** Releasing over the plot records a step. The readout
names the step number and the candidate count; the step controls appear beneath
it. Nothing is stored yet — `Page 1 (0)` and `No cells highlighted` are still
true.
```

```{figure} ../../../_static/screenshots/highlighting_selection/03-confirm-lasso-group.png
:alt: The Highlighting panel after Confirm, showing the page tab count 870, the summary "870 cells highlighted", a Clear button, and one group row labelled "Lasso (870 cells) (870)" with a checkbox and a remove control.
:width: 516px

**Step 3 — confirm.** **Confirm** turns the candidate set into one group on the
active page. The page tab, the summary line and the group row now all carry the
same count, and **Clear** has appeared.
```

```{figure} ../../../_static/screenshots/highlighting_selection/04-view-highlighted-cells.png
:alt: The whole window after confirming the lasso group, with the enclosed cells drawn in the highlight colour across the embedding and the panel reporting 870 cells highlighted.
:width: 1440px

**Step 4 — the result.** The group survives everything that is not a filter:
change the coloured field, move the camera, switch dimension, and these 870
cells stay members.
```

## The other three tools

```{figure} ../../../_static/screenshots/highlighting_selection/select-category-by-click.png
:alt: Annotation based selection after Alt+click on a single cell: every visible cell of that cell_type category is drawn in the highlight colour across the embedding, and the panel reads "Step 1: 591 cells" with the help line "Alt to replace, Shift+Alt to add, Ctrl+Alt to subtract".
:width: 1440px

**Annotation based.** One `Alt+click` on one cell selected all 591 visible cells
of that cell's `cell_type` category — including the ones on the far side of the
embedding. The help line changes to *Alt to replace* on a categorical field,
because intersecting one category with another is almost always empty.
```

```{figure} ../../../_static/screenshots/highlighting_selection/drag-proximity-radius.png
:alt: Proximity drag mid-gesture, with a dashed circle drawn around the seed cell, the cells inside it in the highlight colour, and the pointer pressed on the circle's edge.
:width: 1440px

**Proximity drag.** The dashed circle is the radius in *embedding* space, not
screen space; the small dot at its centre is the seeded cell. In 3D the same
gesture is a sphere. The step readout appears only on release, so mid-drag the
panel still shows the tool's ordinary help text.
```

```{figure} ../../../_static/screenshots/highlighting_selection/drag-knn-degree.png
:alt: KNN drag mid-gesture, with cells around the seed drawn in the highlight colour spreading along the neighbour graph rather than filling a circle, and the pointer pressed away from the seed.
:width: 1440px

**KNN drag.** Compare the shape with the proximity capture above: the selection
follows the neighbour graph, so it spreads along the manifold instead of filling
a disc. Drag distance sets the degree in whole steps.
```

## Two pages

```{figure} ../../../_static/screenshots/highlighting_selection/compare-highlight-pages.png
:alt: The Highlighting panel with two page tabs, a blue Page 1 carrying 870 and an active red Page 2 carrying 250, the group list showing one Proximity group, and the pointer on the button that adds another page.
:width: 516px

Pages are alternative workspaces. Each carries its own colour, its own count and
its own group list; only the active one is rendered. The `⋈` control on a tab
combines it with another page by intersection or union.
```

```{figure} ../../../_static/screenshots/highlighting_selection/share-highlight-across-views.png
:alt: Two 3D panels of the same dataset at different camera angles with the sidebar collapsed; the cells confirmed by one lasso are drawn in the highlight colour in both panels.
:width: 1440px

Highlight membership is one list for the whole dataset, so both panels draw the
same members. What each panel draws is still its own business: a cell a panel is
filtering out is not highlighted there. See
{doc}`04_selection_synchronization`.
```

## What the app says when a gesture fails

```{figure} ../../../_static/screenshots/highlighting_selection/read-missed-click-notice.png
:alt: A proximity Alt+click on empty space, with the pointer pressed where no cell is drawn and the Highlighting panel showing the notice that the click did not land on a cell and that only the cells the viewer is currently drawing can be clicked.
:width: 1440px

A click that hits nothing says so, in one line, under the **Highlight mode**
buttons. That line is the only place these messages appear — no toast, no
history. The full list is in {doc}`06_troubleshooting_highlighting`.
```

```{figure} ../../../_static/screenshots/highlighting_selection/read-no-active-field-notice.png
:alt: The Highlighting panel in Annotation based mode with no coloured field, showing the line that annotation selection needs an active field chosen under Coloring and Filtering.
:width: 516px

With no coloured field, **Annotation based** has no rule to apply, so it replaces
its usual invitation with an instruction instead of failing silently.
```

```{figure} ../../../_static/screenshots/highlighting_selection/read-annotation-step.png
:alt: The Highlighting panel after an annotation Alt+click, with the step readout and the categorical help line "Alt to replace, Shift+Alt to add, Ctrl+Alt to subtract" ringed.
:width: 516px

The help line under the step readout is tool-specific *and* field-specific. Read
it before you press a modifier: it is the only place the app tells you that
`Alt` means **replace** here and **intersect** everywhere else.
```

## The 246 px capture reused by other sections

```{figure} ../../../_static/screenshots/highlighting_selection/highlighting-selected-page.png
:alt: Highlighting panel showing a named page with highlighted cells.
:width: 246px

Highlighted cells are stored on a named page. That page is the explicit group
used by later analysis and session workflows. This capture is pinned at 246 px
because pages in the Python guide reference it at that width; the sharper
version of the same view is `03-confirm-lasso-group.png` above.
```
