# Highlight mental model

Highlights are **annotations on top of visibility**:

- **Filtering** decides which cells are currently visible (per view).
- **Highlighting** marks cells of interest (and can remain meaningful even if those cells are currently filtered out).

If you remember one sentence, make it this:

> Filters control *visibility*. Highlights control *membership*.

## At A Glance

**Audience**
- Wet lab / beginner: learn “selection vs highlight” and how to make a selection persist.
- Computational users: learn the exact set logic (union/intersect/subtract), scope rules, and performance pitfalls.
- Developers: learn what is stored (cell indices), how sessions persist highlights, and why some things are view-scoped.

**Time**
- ~15–25 minutes

**Prerequisites**
- Any dataset loaded (demo is fine)
- Familiarity with basic navigation + fields:
  - {doc}`../c_core_interactions/index`
  - {doc}`../d_fields_coloring_legends/index`
  - {doc}`../e_filtering/index` (recommended)

**What you’ll learn**
- What a highlight page is vs a highlight group
- What “preview highlight” is (and why it disappears)
- Why you cannot select filtered-out cells (by design)
- What persists in sessions vs what is temporary

---

## The mental picture (one diagram)

Think of highlighting as a bookkeeping layer on top of what the view is currently drawing:

```text
Dataset (N cells; fixed order)
│
├─ View A filters → visibility mask A (which cells are interactable in View A)
├─ View B filters → visibility mask B (can differ in snapshots)
│
└─ Highlight pages (global, dataset-scoped)
   ├─ Page 1 (active) → list of highlight groups (each group = a set of cell indices)
   ├─ Page 2
   └─ ...

Renderer (per view):
  draw only visible cells in that view,
  and add a highlight halo for cells that are “highlighted” on the active page.
```

Two key consequences:
- You can have highlighted cells that are *not currently visible*.
- You cannot select cells that a filter has hidden. Level of detail is a
  different thing: when a large dataset is drawn at reduced detail, the cells it
  leaves undrawn are still selectable.

Here is what that looks like once one selection has been confirmed:

```{figure} ../../../_static/screenshots/highlighting_selection/04-view-highlighted-cells.png
:alt: The whole window after confirming the lasso group, with the enclosed cells drawn in the highlight colour across the embedding and the panel reporting 870 cells highlighted.
:width: 1440px

One confirmed group. The panel counts membership (`870 cells highlighted`); the
canvas draws the members it is currently allowed to draw. Those are two
different numbers whenever a filter is active, and the panel says so by
switching to *“N of M highlighted cells visible”*.
```

---

## Definitions (the words people mix up)

### Cell index (what Cellucid stores)
Internally, Cellucid identifies a cell by its **index**: `0..(n_cells-1)` in the loaded dataset’s order.

This makes highlights **dataset-dependent**:
if you reorder/re-export your dataset, “cell #123” may no longer refer to the same biological cell.

:::{important}
Session restore rejects and rolls back the complete operation if the saved
dataset fingerprint does not match the current dataset. It never applies
layout while skipping highlights. This prevents “wrong highlights on the wrong
dataset.”
:::

### Selection candidate set (temporary)
When you use a selection tool, you build up a *candidate set* (“the current selection”).

Key properties:
- It is temporary.
- It can be multi-step (you can refine with multiple gestures).
- It can be confirmed into a highlight group, or canceled.

### Preview highlight (temporary rendering)
While you are mid-selection (especially while dragging), Cellucid shows a **preview highlight**:
it temporarily highlights your current candidate set so you can see “what would be selected if I confirmed now”.

Preview highlight is not saved. It will disappear when you:
- confirm a selection,
- cancel a selection,
- switch highlight pages (which recomputes permanent highlights), or
- reload without restoring a session.

### Highlight group (persistent set of cells)
A highlight group is a named entry in the “Highlighted groups” list.

In the current UI, groups are:
- created by confirming a selection (lasso/proximity/KNN/annotation),
- labeled automatically (e.g. `Lasso (12,345 cells)`),
- enabled/disabled and removable,
- stored as an explicit list of cell indices.

:::{note}
Groups are not “colored labels” yet; Cellucid currently renders highlighted cells using a single highlight style (halo/ring).
Highlight pages have colors in the UI, but highlight groups do not.
:::

### Highlight page (persistent collection of groups)
A highlight page is a container that holds a list of highlight groups.

Why pages exist:
- Wet lab workflow: “my selection for condition A” vs “my selection for condition B”.
- Computational workflow: “hypothesis 1” vs “hypothesis 2”.
- Review workflow: keep multiple alternative selections without destroying prior work.

Only one page is **active** for rendering/highlighting at a time.

---

## The pipeline: visibility first, then highlights

### Step 1: filters decide what is selectable
Selection tools only consider cells that are not filtered out in the view you are selecting in.

That includes:
- standard filters ({doc}`../e_filtering/index`),
- view-specific filters in multi-view snapshots (if snapshots have independent filters).

It does **not** include LOD/downsampling. With **Visualization → Renderer settings →
Level-of-Detail (LOD)** enabled ({doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke`),
a zoomed-out camera makes Cellucid draw only a subset of the points — but a lasso,
proximity drag, or KNN growth still evaluates every unfiltered cell. So a gesture can
select cells the view is not currently drawing, and the group you confirm can be larger
than the dots you drew around. That is deliberate: a selection that changed with the
camera would not be reproducible.

To see how much of a result is actually on screen, read the
“*N* of *M* highlighted cells visible” counter above the group list
({doc}`03_highlight_ui`). That counter *does* account for LOD, so “visible” there can be
a smaller number than “selected”.

One exception: the `Alt+click` that *starts* a gesture picks the cell under the pointer,
and that pick can only land on a cell the view is drawing. Annotation-based, proximity,
and KNN selection all begin that way, so with LOD on and the camera far out a click can
report “no cell” even where cells exist — zoom in to seed, then drag. Lasso has no seed
click and is never affected.

This is why “I can’t lasso those points” means “they are filtered out in this view”,
not “they are too small to be drawn right now”.

### Step 2: highlights decide what is emphasized
Permanent highlights come from:
- all enabled groups on the **active** highlight page.

Preview highlight temporarily adds:
- the current candidate set (while you are mid-selection).

### Step 3: what you see is per view
Each view renders highlights only for cells that are visible in that view.

So in multi-view:
- the *same* highlighted cell may be visible (and highlighted) in View A,
- but invisible (and therefore not highlighted) in View B if View B has stricter filters.

---

## Selection tools are “set operations”

Most selection tools support three operations on your current candidate set:
- **Intersect** (keep only cells in both sets) — “refine / whittle down”
- **Union** (add cells) — “grow”
- **Subtract** (remove cells) — “carve away”

In Cellucid, these operations are controlled by modifier keys held during the gesture:
- `Alt` = start selection gesture
- `Shift+Alt` = union (add)
- `Ctrl+Alt` or `Cmd+Alt` = subtract (remove)

:::{tip}
In 3D, intersection across multiple angles is the most reliable way to isolate a “true 3D region” using 2D tools.
Example: draw a lasso, rotate 90°, draw another lasso with `Alt` to intersect.
:::

For categorical “annotation-based” selection, `Alt` acts like **replace** (not intersect), because intersecting “all cells in category X” with “all cells in category Y” is usually empty and not what users intend.

---

## What persists (and where)

### What persists in a session
When you save a session ({doc}`../l_sessions_sharing/index`), Cellucid can persist:
- highlight pages (names, colors, active page),
- highlight groups (labels/types, enabled flags),
- highlight memberships (cell indices).

Large highlights can make sessions much larger; see {doc}`05_edge_cases_highlighting` for guidance.

### What does not persist automatically
- The temporary candidate set (mid-selection)
- The preview highlight
- Anything if you reload without restoring a session

### Dataset binding (important)
Highlights are dataset-dependent. If the session’s dataset fingerprint doesn’t match your currently loaded dataset, Cellucid restores only dataset-agnostic layout and skips highlight data.

---

## Fast Path (Wet Lab Friendly)

Goal: select a cluster and make it persist.

1) Load any dataset (demo is fine).
2) Open **Highlighting** in the left sidebar.
3) Pick **Lasso**.
4) Hold `Alt` and draw around a cluster.
5) Release mouse → you should see a “Step 1” count and the points preview-highlighted.
6) Click **Confirm** (in the step controls).
7) You now have a highlight group in the list, and the highlight remains even if you change fields.

```{figure} ../../../_static/screenshots/highlighting_selection/choose-highlight-mode.png
:alt: The Highlighting panel with the four Highlight mode buttons — Annotation based, KNN drag, Proximity drag, Lasso — the pointer on Lasso, and the empty group list reading "No cells highlighted".
:width: 516px

Step 3 in the panel. The line under the buttons always describes the mode that
is currently active, so it is the fastest way to confirm you are in the tool you
think you are in.
```

The four steps as they actually look, gesture by gesture, are in
{doc}`07_screenshots`.

Try this once to learn the key idea:
- Apply a filter that hides that cluster → the group still exists, but becomes invisible.
- Remove the filter → the group becomes visible again.

---

## Practical Path (Computational)

### Workflow A: tight 3D selection using intersection
1) Lasso roughly around the region (`Alt`).
2) Rotate to a very different angle.
3) Lasso again with `Alt` (intersect) to reduce false positives from projection.
4) Repeat once more if needed.
5) Confirm.

### Workflow B: graph-aware selection (KNN drag)
1) Enable/show connectivity edges (if available).
2) Use **KNN drag** to grow a selection along the neighbor graph.
3) Union/subtract as needed.
4) Confirm.

### Workflow C: per-condition pages
1) Create/rename pages like `Control`, `Treated`.
2) Do selections and confirm groups inside each page.
3) Use analysis later to compare pages.

### Performance notes (read before selecting 1M cells)
Operations that scale with the number of cells:
- selection preview for lasso/proximity (can scan many cells repeatedly during drags),
- recomputing highlight buffers when you add/remove/enable/disable groups,
- saving sessions (membership arrays are stored in the bundle).

If selection feels laggy, the most effective “first aid” is:
- filter down to fewer visible cells,
- temporarily reduce view count (fewer snapshots),
- avoid long continuous drags that force repeated full scans.

---

## Deep Path (Expert / Developer)

### What is stored
- Highlight membership is stored as arrays of **cell indices**.
- Rendering uses a per-cell `Uint8Array` highlight intensity buffer (binary in practice: `0` or `255`).

### Why the highlight count can say “X of Y visible”
The highlight UI distinguishes:
- **total highlighted** cells (membership), and
- **visible highlighted** cells (after filters and LOD/downsampling in the current view).

If you zoom/pan and see “visible highlighted” change, you may be seeing LOD behavior rather than filters.

### Session storage (high level)
For session bundles, highlight memberships are stored as dataset-dependent chunks and compressed (delta-encoded indices + varint + gzip).
This is efficient for typical selections but can still be large for very big groups.

### Python/Jupyter synchronization (implemented)

Cellucid has a notebook bridge (`cellucid-python`) that supports:
- **Python → UI commands** (e.g. `viewer.highlight_cells(...)`, `viewer.set_color_by(...)`)
- **UI → Python hooks** (e.g. `@viewer.on_selection`, `@viewer.on_click`, `@viewer.on_hover`)
- **No-download session capture**: pull the current `.cellucid-session` into Python as an object and apply it onto an `AnnData`

Recommended patterns:
- If you want to react immediately, use hooks.
- If you want a durable, reproducible artifact (and to mutate an `AnnData`), use a session bundle:

  ```python
  viewer.wait_for_ready()
  bundle = viewer.get_session_bundle(timeout=60)
  adata2 = bundle.apply_to_anndata(
      adata,
      expected_dataset_id="my-study-v1",
      inplace=False,
  )
  ```

Important constraint: selections/highlights in sessions are currently **index-based** (cell identity is the row position), so only apply sessions to an `AnnData` with the same row order.

Developer reference: {doc}`../../python_package/e_jupyter_hooks/16_reference`

---

## Edge cases and troubleshooting
- Edge cases checklist: {doc}`05_edge_cases_highlighting`
- Troubleshooting: {doc}`06_troubleshooting_highlighting`

## Next steps / Related pages
- {doc}`02_selection_tools_document_each_tool`
- {doc}`03_highlight_ui`
- {doc}`04_selection_synchronization`
- {doc}`../e_filtering/index`
- {doc}`../l_sessions_sharing/index`
