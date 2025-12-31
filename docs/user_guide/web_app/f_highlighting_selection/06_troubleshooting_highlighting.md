# Troubleshooting (highlighting and selection)

Use this page when selection/highlighting does not match your expectation.

Format: **Symptom → likely causes → how to confirm → fix → prevention**.

## Quick triage (start here)

These checks solve most “it’s broken” reports in under a minute:

1) **Filters:** are you filtering out the cells you expect to see?
   - See `../e_filtering/index`
2) **Active page:** are you on the highlight page you think you are on?
   - Switching pages changes what is highlighted.
3) **Group enabled:** is the group checkbox enabled?
   - Disabled groups do not contribute to highlighting.
4) **Visible vs total:** does the highlight count say “X of Y visible”?
   - That usually means highlights exist, but are hidden by filters or LOD.
5) **View panel:** are you selecting in the right panel (live vs snapshot)?
   - See `../c_core_interactions/04_view_layout_live_snapshots_small_multiples`
6) **Modifier keys:** are you holding `Alt` to start a selection gesture?
   - `Shift+Alt` = add; `Ctrl/Cmd+Alt` = subtract.
7) **Tool mode:** are you in the tool you think you are in (Lasso vs Proximity vs KNN vs Annotation based)?

If these do not resolve it, jump to the symptom that matches what you see.

---

## Symptom: “Nothing highlights (I created a group but can’t see it)”

### Likely causes (ordered)
1) The highlighted cells are **filtered out** (visibility = 0).
2) You are on the wrong **highlight page**.
3) The group exists but is **disabled**.
4) The group exists, but **LOD/downsampling** makes it hard to see at the current zoom.
5) You created the group on a different dataset (session/dataset mismatch) and it never restored.

### How to confirm
- Look at the highlight count:
  - if it says “0 of N visible”, the group exists but is not visible right now.
- Temporarily disable filters:
  - if the highlight appears, it was a visibility problem, not a highlight problem.
- Switch pages and watch the group list:
  - if the group list changes, you were on the wrong page.
- Check the group checkbox state:
  - disabled groups appear “inactive”.

### Fix
- Disable/relax filters first.
- Ensure you are on the intended page tab.
- Re-enable the group checkbox.
- Zoom in (and/or increase point size, if available) to make the highlight halo easier to see.

### Prevention
- Name pages clearly (e.g., `Control`, `Treated`, `Hypothesis A`).
- Save a session after creating important highlights (`../k_sessions_sharing/index`).
- When selecting, avoid having “surprise filters” active (especially if using categorical selection).

---

## Symptom: “The highlight count says cells are highlighted, but the plot looks unhighlighted”

### Likely causes (ordered)
1) Highlight style is subtle relative to your color-by/background/point size.
2) Only a small fraction of points are highlighted and you are zoomed out.
3) The highlighted cells are mostly filtered out (so only a few are visible).

### How to confirm
- Switch to a simple color-by (or a plain-looking field) temporarily.
- Zoom in until points are larger on screen.
- Check whether the count says “X of Y visible”.

### Fix
- Zoom in; highlights are easiest to see when points are not subpixel-sized.
- Adjust display settings (point size/background) if your build exposes them.
- Temporarily simplify filters to increase visible highlighted cells.

### Prevention
- When validating a selection, use a plain or high-contrast coloring.
- Capture “before/after” screenshots when documenting a workflow (`07_screenshots`).

---

## Symptom: “Lasso selects the wrong cells”

### Likely causes (ordered)
1) **3D projection ambiguity**: lasso selects by screen projection; depth is ambiguous.
2) You started the lasso in the wrong multi-view panel (grid).
3) Filters changed between steps; visible set differs.
4) You intended union/subtract but used intersect (or vice versa) due to modifier keys.

### How to confirm
- If you are in 3D: rotate the camera and repeat the same lasso. If the selected set changes a lot, you are seeing projection ambiguity.
- If you are in multi-view grid: start the lasso clearly inside the intended panel (not near borders).
- Check the step label (it often indicates whether you are intersecting/adding/subtracting).

### Fix
- In 3D, use the intended workflow:
  1) lasso once,
  2) rotate substantially,
  3) lasso again with `Alt` to **intersect**.
- If you need a true “local region” in 3D, try **Proximity drag** instead of lasso.
- In grid multi-view, switch to “Edit selected view” temporarily to reduce panel ambiguity.

### Prevention
- Treat lasso in 3D as “projected selection”, not volumetric selection.
- Use multi-step intersection for 3D.
- Keep lasso gestures fully inside one panel in grid view.

---

## Symptom: “I can’t start lasso / nothing happens when I draw”

### Likely causes (ordered)
1) You are not holding `Alt` (lasso mode requires `Alt`).
2) You are not in **Lasso** mode (toolbelt mismatch).
3) You are in a navigation interaction state (dragging/rotating) and never entered lasso.
4) Browser/OS intercepts `Alt` in your environment (rare but possible).

### How to confirm
- Click the **Lasso** mode button again.
- Hold `Alt` and look for the cursor change (crosshair) or lasso overlay.
- Try in a different browser window with no extensions.

### Fix
- Use `Alt+drag` (not plain drag).
- If your OS/browser conflicts with `Alt`, try:
  - a different browser (Chrome/Edge/Firefox),
  - disabling extensions that capture keyboard shortcuts.

### Prevention
- Learn the mental split:
  - drag = navigation,
  - `Alt+drag` = selection tool.

---

## Symptom: “Proximity selection selects nothing / selects the wrong region”

### Likely causes (ordered)
1) You clicked empty space, so there was no valid center (especially with no existing candidate set).
2) Your radius is too small (tiny drag) or too large (huge drag).
3) Filters are hiding nearby cells; selection only considers visible cells.
4) You expected screen-space selection; proximity is embedding-space distance.

### How to confirm
- Start by `Alt+click` directly on a clearly visible cell (zoom in).
- Drag slowly and watch the overlay radius.
- Temporarily disable filters and try again.

### Fix
- Zoom in and click on a cell to set a stable center.
- Use small drags first; increase radius gradually.
- If you need geometric control in screen space, use lasso; if you need embedding-space locality, use proximity.

### Prevention
- Use proximity for “local blob around a seed”, not for “draw around this shape”.
- Avoid selecting with proximity when the view is heavily filtered unless that is your intention.

---

## Symptom: “KNN drag says neighbor graph not available (or does nothing)”

### Likely causes (ordered)
1) Your dataset export does not include connectivities/edges.
2) Edges exist but are not loaded yet (often loaded when you enable “Show edges”).
3) You are using a data loading mode that doesn’t provide connectivities.

### How to confirm
- Look for an edges/graph toggle (e.g., “Show edges”) and enable it.
- If there is an edges overlay, confirm you can see edges.
- If you cannot enable edges, assume connectivities are missing.

### Fix
- Re-export your dataset including connectivities (if you control the export pipeline).
- Use lasso/proximity/annotation-based selection if you cannot load a graph.

### Prevention
- For KNN workflows, treat “connectivities present” as a prerequisite and document it explicitly.

---

## Symptom: “Annotation based selection selects too many / too few cells”

### Likely causes (ordered)
1) Filters are hiding part of the category/range (selection only considers visible cells).
2) You are selecting on a different active field than you think.
3) The field has missing values or extreme distribution (continuous).
4) You expected “Alt = intersect” but on categorical fields it behaves like “replace”.

### How to confirm
- Confirm the active field in the UI (legend/field selector).
- Temporarily disable filters and repeat the selection.
- Try a known categorical field with obvious labels (e.g., `cell_type`) to validate semantics.

### Fix
- For categorical:
  - use `Alt` to replace,
  - `Shift+Alt` to add,
  - `Ctrl/Cmd+Alt` to subtract.
- For continuous:
  - use `Alt+drag` to tune the range,
  - intersect with lasso/proximity if you need geometry constraints.

### Prevention
- Don’t do annotation-based selection while actively changing color-by/fields; treat the active field as part of the selection “spec”.

---

## Symptom: “Undo/redo/confirm buttons are disabled or missing”

### Likely causes (ordered)
1) There is no in-progress candidate set (you haven’t completed a selection step yet).
2) You changed tools/pages and cleared the UI state (candidate set may or may not still exist).
3) Candidate set is empty after subtracting everything.

### How to confirm
- Perform one fresh selection step (e.g., a small lasso) and see whether step controls appear.
- If step controls appear but Confirm is disabled, you likely have an empty candidate set.

### Fix
- Build a non-empty candidate set, then Confirm.
- If you got into a weird state, click Cancel and start a fresh selection.

### Prevention
- Confirm or cancel before switching pages.
- Save sessions if you need “checkpointing”.

---

## Symptom: “Selection is delayed / laggy”

### Likely causes (ordered)
1) You are selecting on a very large dataset (CPU-heavy selection preview).
2) You are in multi-view with many panels (more rendering/overhead).
3) You have huge highlight groups already; toggling/recomputing is expensive.
4) Browser/GPU constraints (integrated GPU, low memory, heavy background tabs).

### How to confirm
- Does lag correlate with dragging (lasso/proximity/KNN preview) rather than simple clicks?
- Does lag disappear if you disable filters and reduce to one view panel?
- Does lag correlate with enabling/disabling groups?

### Fix (safe first)
1) Reduce the visible set with filters (counterintuitive but effective): fewer visible cells → fewer selectable candidates.
2) Use fewer snapshot panels (temporarily).
3) Prefer tools that do less repeated preview scanning:
   - use shorter gestures,
   - confirm sooner,
   - avoid long “slow drags”.
4) Reload the page and restore a session if needed.

### Prevention
- For very large datasets, plan for performance:
  - avoid “live preview while dragging for 10 seconds” workflows,
  - keep highlight groups coarse and few,
  - save sessions at milestones.

---

## Symptom: “Highlights disappear when switching views/pages (or after reload)”

### Likely causes (ordered)
1) You switched highlight pages (different groups).
2) You are looking at a snapshot view with different filters (highlights exist but are not visible there).
3) You reloaded without restoring a session (highlights are not auto-persisted in URLs).
4) You restored a session on the wrong dataset (dataset-dependent chunks skipped).

### How to confirm
- Check the active page tab and group list.
- Check whether the highlight count says “X of Y visible”.
- If you restored a session, look for a warning about dataset mismatch.

### Fix
- Switch to the intended highlight page.
- Disable/relax filters in the view you are looking at.
- Restore the correct session on the correct dataset.

### Prevention
- Treat highlights as part of the session artifact: save sessions intentionally.
- Use clear page names and avoid “Page 1 / Page 2” once you have real work.

---

## Symptom: “My highlights didn’t restore from a session”

### Likely causes (ordered)
1) Dataset mismatch: session restore skipped dataset-dependent highlights.
2) Session file is still loading lazy chunks (large highlight memberships can load later).
3) You never saved highlights into the session you loaded (you loaded an older session).

### How to confirm
- Look for a notification about “dataset mismatch”.
- Wait a few seconds: large sessions can apply highlight memberships lazily.
- Check whether pages/groups appear first but counts fill in later (a sign of lazy restore).

### Fix
- Load the session on the dataset it was created from.
- If you must change the dataset, you must recreate highlights (indices won’t match).

### Prevention
- Keep dataset IDs stable and avoid reordering cells between exports if you rely on highlights.
- Save “milestone sessions” before heavy editing.

---

## Symptom: “Python highlight doesn’t show in the UI” (Jupyter integration)

### Likely causes (ordered)
1) You are not actually in Jupyter-embedded mode (no Jupyter bridge active).
2) Viewer ID mismatch (messages go to a different iframe).
3) Dataset/index mismatch (Python indices refer to different cells than the viewer).
4) The particular command/event wiring is not enabled in your app build yet.
5) Browser blocked cross-origin communication (non-loopback origin).

### How to confirm
- Confirm you are using the Jupyter workflow (Cellucid opened with Jupyter-specific URL params).
- In Python, print/log the viewer URL and verify it is the one embedded.
- Try a simple command that is known to be wired (e.g., `set_color_by`) to validate communication.
- If nothing works, check `cellucid/markdown/HOOKS_DEVELOPMENT.md` (developer-facing).

### Fix
- Ensure you are running the server locally/loopback and using the correct viewer URL.
- Ensure the dataset order in Python matches the dataset loaded in the viewer.
- If the feature is not wired in your build, use session-based workflows instead.

### Prevention
- When using Python control, treat “cell index = identity” and keep ordering stable.
- Prefer pre-exported datasets for reproducibility.

---

## When to file a bug report

If you believe you hit a real bug (not one of the expected edge cases), include:
- Browser + OS (and whether you’re on a laptop with integrated GPU)
- Dataset format and rough size (`n_cells`, `n_genes`, whether connectivities exist)
- Whether multi-view snapshots are enabled (and how many)
- Exact selection tool used + exact modifier keys held
- Whether filters were active
- A minimal reproduction recipe (click-by-click)
- Screenshots with callouts:
  - use `07_screenshots` for capture specs

---

## Related pages
- `05_edge_cases_highlighting`
- `07_screenshots`
