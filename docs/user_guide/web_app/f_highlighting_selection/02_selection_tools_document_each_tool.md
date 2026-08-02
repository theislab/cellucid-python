# Selection tools (annotation, lasso, proximity, KNN)

Selection tools create a **temporary candidate set** (“the current selection”). You then **confirm** that selection to create a new highlight group on the active highlight page.

This page documents each tool, how the modifier keys work, and the tricky parts (especially in 3D and multi-view).

## At A Glance

**Audience**
- Wet lab / beginner: follow the Fast Path to select a cluster and save it as a highlight group.
- Computational users: learn the exact set operations, 2D/3D differences, and performance behavior.
- Developers: learn the concrete inputs/outputs (cell indices) and data requirements (connectivity, active fields).

**Time**
- ~20–40 minutes (longer if you try every tool carefully)

**Prerequisites**
- A dataset loaded with at least one embedding
- For KNN drag: a neighbor graph/edges must be available and loaded

**What you’ll learn**
- Which tool to use for which problem
- How `Alt`, `Shift+Alt`, `Ctrl/Cmd+Alt` change selection behavior
- What “screen-space” vs “embedding-space” selection means
- How multi-step selection works (confirm / undo / redo / cancel)

---

## Quick mental model (read this once)

If you haven’t yet, skim {doc}`01_highlight_mental_model` first.

Key idea:
- You are not “editing a highlight group” live.
- You are building a **candidate set**, step-by-step.
- When you click **Confirm**, Cellucid creates **one new highlight group** whose membership is the candidate set at that moment.

Practical implication:
- If you want one clean group, do all your unions/subtractions/intersections *before* confirming.
- If you confirm multiple times, you will get multiple groups in the list.

---

## Universal rules (apply to all tools)

### Rule 1: selection only sees unfiltered cells
Selection tools only select cells that are visible in the view you are selecting in (filters apply).

If a cell is filtered out:
- you cannot select it with any tool,
- but it can still be part of an existing highlight group from earlier (it’s just not visible right now).

Filters are the *only* thing that narrows what a gesture can pick up. LOD/downsampling
does not: with Level-of-Detail enabled and the camera zoomed out, a lasso, proximity drag,
or KNN growth also takes unfiltered cells the view is not drawing at that moment
({doc}`01_highlight_mental_model`). The `Alt+click` that seeds annotation-based, proximity,
and KNN selection is the exception — it can only land on a cell the view is drawing.

### Rule 2: modifier keys are set operations
Most tools interpret modifiers like this:
- `Alt` = start selection gesture
- `Shift+Alt` = **union** (“add to selection”)
- `Ctrl+Alt` or `Cmd+Alt` = **subtract** (“remove from selection”)

For the **annotation-based** tool on categorical fields, `Alt` behaves like **replace** (not intersect). Details below.

### Rule 3: selections are multi-step (and undoable)
Each gesture typically becomes a “step”.

Step controls (per tool) let you:
- **Confirm**: commit the current candidate set into a new highlight group
- **Undo / Redo**: roll back up to a few recent steps
- **Cancel**: discard the current candidate set

:::{note}
Undo/redo history is intentionally short (to keep memory bounded). If you need a very long “edit history”, confirm intermediate groups and keep them (or save a session).
:::

### Rule 4: a gesture must end on the plot
Release the mouse over the plot. A gesture released anywhere else — over the
floating sidebar, over a notification toast, outside the window — is abandoned
on purpose, because the polygon or radius you finished behind a panel would
select cells the panel was covering. Nothing is committed, no step is recorded,
no undo is spent, and your standing selection is repainted exactly as it was.
The tool explains it under the **Highlight mode** buttons: *“That gesture ended
outside the view, so nothing was selected.”*

### Rule 5: a gesture the app cannot serve says so
Every gesture the app consumes but cannot turn into a step writes one line under
the **Highlight mode** buttons — a click that hit no cell, a cell with no value
on the active field, a gesture that resolved to the selection you already had.
That line is the first thing to read when a tool “does nothing”; the full list
is in {doc}`06_troubleshooting_highlighting`.

### Rule 6: you can switch tools mid-selection
Cellucid keeps one unified candidate set across selection tools.

This enables workflows like:
1) start with **annotation-based** selection to grab “all visible cells in category X”
2) switch to **lasso** and subtract a corner
3) switch to **proximity** and add a small blob
4) confirm

The candidate set can persist across tool switches, but per-tool undo/redo history may reset when switching.

---

## Tool overview (choose quickly)

| Tool | Best for | Space | 2D vs 3D | Requires extra data | Common gotcha |
|---|---|---|---|---|---|
| **Annotation based** | “Select all cells like this one” (category) or “cells in this value range” (continuous) | Field space | Works in both; depends on active field | Active field must be categorical/continuous | Selects only *visible* cells in that category/range |
| **Lasso** | Drawing around a region in the plot | Screen-space polygon (projection) | In 3D it selects by projection | None | In 3D it can include points “behind” your intended region |
| **Proximity drag** | Selecting a local blob by distance | Embedding-space radius | In 2D it behaves like a disk; in 3D like a sphere | None | If you miss-click, center may be wrong (or nothing happens) |
| **KNN drag** | Selecting along the neighbor graph (manifold-aware) | Graph distance (degree) | Independent of camera; still view-filtered | Connectivity/edges must be loaded | “Neighbor graph not available” until edges are loaded |

---

## Tool 1 — Annotation-based selection (Alt+click / Alt+drag)

This is the “select cells similar to this one” tool. It depends on the **active field**.

### Where it gets its meaning
The active field is whatever the UI is currently using for:
- coloring / legends, and/or
- field selection

If there is no active field, annotation-based selection has nothing to operate on.

### Categorical field behavior (Alt = replace)
If the active field is **categorical**:
- `Alt+click` selects **all visible cells in the clicked cell’s category**.
- `Shift+Alt+click` unions that category into the candidate set.
- `Ctrl/Cmd+Alt+click` subtracts that category from the candidate set.

Why `Alt` means “replace” here:
intersecting “category A” with “category B” is usually empty and not what users want.

:::{tip}
If you want to select *all* cells in a category (ignoring current filters), temporarily disable filters first, then do the category selection.
:::

### Continuous field behavior (Alt = intersect)
If the active field is **continuous**:
- `Alt+click` selects a small range around the clicked value (a “default window”).
- `Alt+drag` (vertical drag) expands/contracts a range anchored at the clicked value.
- `Shift+Alt` unions the range results into the candidate set.
- `Ctrl/Cmd+Alt` subtracts the range results from the candidate set.

During `Alt+drag`, Cellucid shows a **range preview label** and preview-highlights what would be selected.

```{figure} ../../../_static/screenshots/highlighting_selection/select-category-by-click.png
:alt: Annotation based selection after Alt+click on a single cell: every visible cell of that cell_type category is drawn in the highlight colour across the embedding, and the panel reads "Step 1: 591 cells" with the help line "Alt to replace, Shift+Alt to add, Ctrl+Alt to subtract".
:width: 1440px

One `Alt+click` on one cell, and every visible cell sharing its `cell_type`
becomes the candidate set — including cells nowhere near the pointer. The count
in the step readout (`591`) is the candidate set; `No cells highlighted` above
it is still true, because nothing has been confirmed.
```

### When to use
- You want “all visible cells in this label” (fast categorical selection).
- You want “cells in this numeric window” without drawing geometry.
- You want to bootstrap a selection before refining it with lasso/proximity/KNN.

### Common pitfalls
- “It selected fewer than expected”: filters are hiding some category members.
- “It selected way too many”: your continuous range is too wide; reduce drag distance or intersect with a lasso.
- “Nothing happens, and the panel says annotation selection needs an active
  field”: no field is coloured at all — choose one under **Coloring &
  Filtering**. That sentence replaces the `Alt+click` invitation exactly while
  the precondition is missing, so the panel and the viewer never disagree about
  whether a gesture can start. More symptoms:
  {doc}`06_troubleshooting_highlighting`.

---

## Tool 2 — Lasso (Alt+draw)

Lasso selects based on a polygon you draw on screen. In 3D, this is a **projection** tool.

### Gesture
1) Choose **Lasso** mode.
2) Hold `Alt` and draw a closed-ish shape.
3) Release mouse → that completes a step.

```{figure} ../../../_static/screenshots/highlighting_selection/01-draw-lasso.png
:alt: The 2D planar Pancreas embedding mid-lasso: Alt is held, the pointer is pressed, and a freehand path encloses the Ductal lobe on the left while the enclosed cells are already drawn in the highlight colour.
:width: 1440px

Mid-draw. The path does not have to be closed by hand — releasing closes it. The
cells already repainted inside the path are the **preview**: the exact set the
release would commit under the modifier currently held.
```

Modifiers during the draw:
- `Alt` = intersect
- `Shift+Alt` = union (add)
- `Ctrl/Cmd+Alt` = subtract (remove)

Releasing records a step and reveals the step controls:

```{figure} ../../../_static/screenshots/highlighting_selection/02-review-lasso-step.png
:alt: The Highlighting panel after the lasso was released, showing the step readout "Step 1: 870 cells" and the Confirm, undo, redo and Cancel step controls, with the pointer on Confirm.
:width: 516px

`Step 1: 870 cells`, then **Confirm**, undo (`↩`), redo (`↪`) and **Cancel**.
Undo history is bounded at five steps; redo is greyed out until you have used
undo at least once.
```

### How it works (2D vs 3D)
- In **2D**, the polygon matches what you see: points whose 2D screen projection lies inside the polygon are selected.
- In **3D**, lasso still selects by **screen projection**. This can include points at different depths that happen to project into your polygon.

This is not a bug; it’s the fundamental ambiguity of selecting 3D points with a 2D gesture.

### The intended 3D workflow (multi-step intersection)
To approximate a “3D volume selection”:
1) Lasso from one angle.
2) Rotate camera substantially.
3) Lasso again with `Alt` (intersect).
4) Repeat until the result is tight.

### Multi-view behavior
In grid multi-view, lasso applies to the specific panel you started the gesture in.
Try to keep the lasso path inside one panel to avoid ambiguity.

### Performance considerations
Lasso preview can be expensive on very large datasets because it may need to evaluate many points repeatedly while you draw.
If lasso feels laggy:
- filter down to fewer visible cells first,
- draw shorter/fewer lasso segments (avoid long slow drags),
- consider proximity or KNN selection.

---

## Tool 3 — Proximity drag (Alt+click, drag radius)

Proximity selection chooses cells within a radius around a 3D point in embedding space.

### Gesture
1) Choose **Proximity drag** mode.
2) Hold `Alt` and click on (or near) a cell to set the center.
3) Drag outward to increase radius (a circle overlay shows the radius).
4) Release mouse → that completes a step.

Modifiers:
- `Alt` = intersect
- `Shift+Alt` = union (add)
- `Ctrl/Cmd+Alt` = subtract (remove)

```{figure} ../../../_static/screenshots/highlighting_selection/drag-proximity-radius.png
:alt: Proximity drag mid-gesture, with a dashed circle drawn around the seed cell, the cells inside it in the highlight colour, and the pointer pressed on the circle's edge.
:width: 1440px

Mid-drag. The dashed circle is the current radius and the dot at its centre is
the seeded cell. Note that the panel still shows the tool's ordinary help text:
proximity publishes its step on release, not while you drag.
```

### What “distance” means here
Distance is computed in embedding coordinates (the same coordinate system used for plotting points).

Practical implications:
- zooming changes how much you drag to achieve a given radius,
- in 3D, it behaves like a sphere; in 2D embeddings (z≈0), it behaves like a disk.

### Miss-click behavior (important)
If you click empty space:
- with no existing candidate set there is no meaningful center, so the gesture
  never starts and the tool says *“That click did not land on a cell, so there
  is nothing to select…”* under the **Highlight mode** buttons.
- with an existing candidate set the tool places a center anyway, where your
  click ray crosses the plane through that selection’s centroid facing the
  camera — that is, *where you clicked, at the selection’s depth*. This is what
  makes “grow/shrink around the current
  selection” work in 3D without seeding on a cell, but the blob is centred on
  your pointer, not on the selection.

### When to use
- You want “all cells near this seed cell” in embedding space.
- You want to add/subtract small blobs cleanly (less projection ambiguity than lasso in 3D).

### Performance considerations
Proximity preview can be expensive for huge datasets because it can scan many points while you drag.
If it lags, reduce visible cells first (filters) or use a smaller radius/fewer drags.

---

## Tool 4 — KNN drag (Alt+click, drag degree)

KNN drag grows a selection outward along the neighbor graph (connectivities).

### Data requirement (non-negotiable)
KNN selection requires a neighbor graph.

Typical symptoms when it’s missing:
- the UI indicates “neighbor graph not available”, or
- KNN mode appears but selecting does nothing.

Practical fix:
- load connectivity/edges for your dataset (often via a “Show edges” toggle), or
- re-export your dataset including connectivities.

### Gesture
1) Choose **KNN drag** mode.
2) Hold `Alt` and click a seed cell.
3) Drag away from the seed to increase degree (degree increases in discrete steps).
4) Release mouse → that completes a step.

Modifiers:
- `Alt` = intersect
- `Shift+Alt` = union (add)
- `Ctrl/Cmd+Alt` = subtract (remove)

```{figure} ../../../_static/screenshots/highlighting_selection/drag-knn-degree.png
:alt: KNN drag mid-gesture, with cells around the seed drawn in the highlight colour spreading along the neighbour graph rather than filling a circle, and the pointer pressed away from the seed.
:width: 1440px

The same seed cell and a comparable drag distance as the proximity capture
above, and a visibly different result: the selection follows the neighbour
graph, so it grows along the manifold rather than filling a disc.
```

### What “degree” means (practical)
Degree `0` = the seed cell only.  
Degree `1` = the seed plus its immediate neighbors.  
Degree `2` = neighbors-of-neighbors, etc.

The step readout names the degree it has reached, e.g. `Step 1: 90 cells (2°
neighbors)`. As degree grows, selection can grow very quickly.

### When to use
- You want “manifold-aware” selection that follows biological neighborhoods rather than Euclidean blobs.
- You want to grow a cluster boundary gradually from a seed.

### Common pitfalls
- Filters can “break” the graph: filtered-out cells are treated as non-selectable, which can change which neighbors are reachable.
- Degree grows fast; use small drags first.

---

## “I selected cells — now what?”

### Option A (recommended): make one clean group
1) Use union/subtract/intersect steps until the candidate set matches what you want.
2) Click **Confirm**.
3) A new highlight group appears in the list (on the active page).

```{figure} ../../../_static/screenshots/highlighting_selection/03-confirm-lasso-group.png
:alt: The Highlighting panel after Confirm, showing the page tab count 870, the summary "870 cells highlighted", a Clear button, and one group row labelled "Lasso (870 cells) (870)" with a checkbox and a remove control.
:width: 516px

After **Confirm**. The group is named for the tool that made it and carries its
count twice — once inside the label the tool wrote, once as the live membership
count beside the checkbox. Groups cannot be renamed or recoloured; see
{doc}`03_highlight_ui`.
```

### Option B: keep intermediate groups
1) Confirm an early coarse selection (Group 1).
2) Continue refining and confirm again (Group 2).
3) Enable/disable groups to control what is currently highlighted (and what analysis uses).

### Option C: use pages for alternative hypotheses
1) Create a new highlight page.
2) Do the same selection workflow again.
3) Switch pages to compare different “versions” of your highlights.

---

## Related pages
- {doc}`01_highlight_mental_model`
- {doc}`03_highlight_ui`
- {doc}`04_selection_synchronization`
- {doc}`05_edge_cases_highlighting`
- {doc}`06_troubleshooting_highlighting`
