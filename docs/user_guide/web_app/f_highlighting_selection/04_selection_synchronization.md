# Selection synchronization (views, pages, filters, Python)

Highlighting feels “simple” until you use multi-view snapshots, change filters, or try to coordinate with Python/Jupyter. This page makes the scope rules explicit so you can predict what will happen.

## At A Glance

**Audience**
- Wet lab / beginner: understand “why did my selection disappear when I clicked another panel?”
- Computational users: understand reproducibility and what is view-scoped vs global.
- Developers: understand which state lives in the viewer vs the app state, and what is session-persisted.

**Time**
- ~15–25 minutes

**Prerequisites**
- Familiarity with views/snapshots:
  - `../c_core_interactions/04_view_layout_live_snapshots_small_multiples`
- Familiarity with the highlight UI:
  - `03_highlight_ui`

**What you’ll learn**
- What syncs across views/snapshots
- What switching highlight pages changes (and what it does not)
- How filters affect selection and highlight visibility
- What to expect from Python/Jupyter integration (and how to diagnose when it’s not active)

---

## The scopes table (bookmark this)

| Thing | Scope | Persists? | Notes |
|---|---|---|---|
| **Filters / visibility mask** | Per view (live + each snapshot may differ) | Session (dataset-dependent) | Filters control which cells are selectable and drawn in that view. |
| **Selection tool mode** (annotation / lasso / proximity / KNN) | Global UI choice | Session (dataset-dependent layout) | Affects how `Alt` gestures are interpreted. |
| **Selection candidate set** (in-progress selection) | Global (shared across tools) | No | Temporary; can be confirmed into a group or canceled. |
| **Preview highlight** (while selecting) | Global, rendered per view | No | Shown only for visible cells in each view. |
| **Highlight pages** | Global (dataset-scoped) | Yes (session, dataset-dependent) | Pages are alternate “workspaces” of highlight groups. |
| **Highlight groups** (inside a page) | Global (dataset-scoped), but page-local | Yes (session, dataset-dependent) | Groups store explicit arrays of cell indices. |
| **Active highlight page** | Global | Yes (session, dataset-dependent) | Switching pages changes which groups contribute to highlighting. |
| **Dataset identity** | Global | N/A | Cell indices are meaningful only for a specific dataset ordering. |

:::{important}
If a session’s dataset fingerprint does not match the currently loaded dataset, dataset-dependent state (including highlights) is intentionally skipped on restore.
:::

---

## Sync between views and snapshots

### What users usually expect
People typically expect:
- “If I highlight cells, they should highlight in every panel.”
- “If I select cells in one view, I should see the corresponding cells in other views.”
- “If I switch focus to a different view panel, I shouldn’t lose my selection.”

### What actually happens (current behavior)

**Permanent highlights (confirmed groups)**
- Are global: once you confirm a selection into a group, those cells are highlighted across all panels.
- Are still subject to visibility: each panel will only show highlights for cells that are visible in that panel (filters apply per view).

**In-progress selection (candidate set + preview highlight)**
- The candidate set is global (shared across selection tools).
- The preview highlight is global, but rendered per view (visible cells only).
- The gesture you perform is tied to the specific panel you start it in (especially for lasso).

**Lasso in multi-view grid**
- Lasso selection is anchored to the panel you started drawing in.
- Keep the lasso path inside one panel to avoid “wrong panel” confusion.

**Proximity/KNN in multi-view**
- The seed click happens in one panel, but the resulting selected indices are global.
- Filters in the target panel can affect which cells are selectable (filtered-out cells are excluded).

### Common confusion patterns (and the mental fix)

**“I selected in panel A but panel B didn’t change.”**
- Most likely: panel B is filtering out those cells (visibility differs).
- Or: you have not confirmed (you are only seeing preview highlight, and the other panel is hiding those cells).

**“My selection disappeared when I clicked another panel.”**
- Often: you changed focus, not state.
- Look for the step controls: if they still show your step count and candidate count, the selection is still in progress.

**“I switched from Grid compare to Edit selected view and things changed.”**
- Layout changes can change what you *notice* (one panel full size vs small multiples), but the underlying highlight page/group state is global.

---

## Sync with filtering and color-by

### Filters and selection
Selection tools only select **visible** cells. This is intentional:
- it prevents selecting “ghost” cells you cannot see,
- it makes selection consistent across multi-view where each panel may have different filters.

Practical workflow tip:
- If your goal is “select all cells in category X across the dataset”, disable filters first.

### Filters and confirmed highlights
Confirmed highlights can include cells that are currently filtered out (because you confirmed them earlier when they were visible).

This is why the highlight count can show:
- “X of Y highlighted cells visible”

### Color-by and selection/highlighting
Changing the active color-by field:
- does not change group membership,
- can change which selection tool behaviors are available for **annotation-based** selection (because that tool depends on the active field’s kind: categorical vs continuous),
- can change how easy it is to visually judge the correctness of your selection.

---

## Sync with highlight pages

### What switching pages does
Switching the active highlight page:
- changes which groups contribute to highlighting,
- changes what appears in the group list,
- changes the highlight count and what is emphasized in the canvas.

### What switching pages does not do
Switching pages does not:
- change filters (visibility masks),
- change the dataset,
- automatically “merge” pages (you must explicitly combine pages via union/intersection).

### Page combine (union / intersection) is a reproducibility tool
Because page combination is a deterministic set operation, it’s a good way to:
- create a derived page that’s stable and shareable (via sessions),
- avoid manual “redo selection” when you want intersections/unions of prior work.

---

## Python/Jupyter synchronization (advanced)

Cellucid has a Jupyter bridge via `cellucid-python`. Conceptually, it supports two directions:

### Python → UI (commands)
Examples from `cellucid-python` include:
- `viewer.highlight_cells([...])`
- `viewer.set_color_by("cell_type")`
- `viewer.set_visibility([...], visible=False)`

These commands identify cells by **cell index**. That means they only make sense when:
- the dataset in the viewer is the one you think it is, and
- indices are aligned with your Python object (AnnData ordering).

### UI → Python (events / hooks)
UI interactions POST events to the Python server (`/_cellucid/events`), so hooks like `@viewer.on_selection` can fire.

Important nuance for selection:
- In notebooks, the **selection event is emitted when you confirm a selection** (i.e. when it becomes a real highlight group),
  not while you are mid-dragging with a preview highlight.
- Hover/click events are high-frequency; hover is debounced.

If you are debugging Python synchronization, start from:
- `cellucid/markdown/HOOKS_DEVELOPMENT.md`

### Common failure modes (Python/Jupyter)
- Wrong `viewerId` (messages routed to the wrong embedded viewer)
- Dataset mismatch (indices don’t refer to the same cells)
- Server unreachable from the viewer (HTTPS/remote notebooks without Jupyter Server Proxy, or remote/HPC without tunneling)
- Embedded UI assets blocked (hosted-asset proxy blocked and no cached copy available)

First-line diagnostic:
- `viewer.debug_connection()` (server probes + ping/pong + frontend debug snapshot + forwarded console errors)

---

## Related pages
- `01_highlight_mental_model`
- `02_selection_tools_document_each_tool`
- `03_highlight_ui`
- `05_edge_cases_highlighting`
- `06_troubleshooting_highlighting`
- `../g_cross_highlighting/index`
