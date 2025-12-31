# Quick tour (60 seconds)

**Audience:** first-time users (wet lab + computational)  
**Time:** 1–3 minutes  
**What you’ll do:** load → move → color → keep a snapshot → highlight → export/save

:::{tip}
If you prefer to learn by definition-first, read `a_orientation/04_ui_glossary_terminology` before this tour.
:::

---

## One screenshot “map” (recommended)

This tour references UI areas by name (the app labels). A single annotated screenshot makes the rest of the docs much easier to follow.

<!-- SCREENSHOT PLACEHOLDER
ID: quick-tour-ui-map
Where it appears: Quick tour → One screenshot “map”
Capture:
  - Load any small demo dataset
  - Make sure the left sidebar is open
  - Show: Session accordion, Coloring & Filtering, Compare Views, Navigation controls, Highlighting, Analysis, Figure Export
Crop:
  - Include: full left sidebar + enough of the plot canvas to orient the reader
  - Exclude: browser bookmarks, personal tabs, any private dataset names
Annotations:
  - Add 6–10 numbered callouts matching the tour steps below (1=Session, 2=Coloring, 3=Compare Views, ...)
Alt text:
  - Cellucid UI with the main sidebar sections labeled.
Caption:
  - The main Cellucid sidebar sections used in this quick tour.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for an annotated Cellucid UI map.
:width: 100%

The main Cellucid sidebar sections used in this quick tour.
```

---

## Fast path (click-by-click)

### 1) Load a dataset (demo is fine)

1) Open the **Session** accordion.
2) Under **Sample datasets**, choose any dataset.
3) Wait for the dataset info panel to populate (Cells/Genes/Obs fields).

**What success looks like**

- The canvas shows points (usually on a grid background).
- The dataset info panel shows non‑zero counts.

:::{note}
You can also load your own data from **Local data** (H5AD/Zarr/Prepared) or from a **Remote server**, but those workflows are documented in `b_data_loading/index`.
:::

### 2) Move the camera (don’t skip this)

1) Open **Compare Views** → **Navigation**.
2) Keep the default mode: **Orbit**.
3) Try:
   - click‑drag to rotate,
   - mouse wheel / trackpad scroll to zoom,
   - right‑drag (or Shift‑drag) to pan.

If you are on a trackpad and the interaction feels “wrong”, try **Planar** mode (best for 2D embeddings).

### 3) Color by a field (metadata or gene expression)

1) Open **Coloring & Filtering**.
2) Pick one:
   - **Categorical obs** (clusters, batch, sample), or
   - **Continuous obs** (QC metrics, scores).
3) Look for the legend changing on the left.

**What success looks like**

- The points change color.
- A legend appears/updates (categories or a color scale).

### 4) Keep a snapshot (small multiple)

1) Open **Compare Views** → **Multiview**.
2) Click **Keep view**.
3) You should see a new numbered badge/panel.

This is how you compare hypotheses side‑by‑side without losing your current view.

### 5) Make a highlight (selection)

1) Open **Highlighting**.
2) Choose a mode (e.g., **Lasso**).
3) Draw a selection; a highlight group should appear in the highlight list.

If you want a keyboard-first workflow, the UI supports:

- `Alt` + drag for area selection (and other shortcuts; see the in‑app “Keyboard Shortcuts” accordion).

### 6) Export or save what you did

Pick one:

- **Save State** (Session accordion) to download a `.cellucid-session` bundle you can reopen later.
- **Figure Export** to export a publication‑ready figure of the current view(s).

---

## Troubleshooting (quick tour blockers)

### “The canvas is blank” / “WebGL2 required”

Go to `a_orientation/02_system_requirements`.

### “I can’t rotate/zoom/pan”

Go to `c_core_interactions/06_troubleshooting_core_interactions` (focus, pointer lock, nav mode confusion).

### “Keep view does nothing”

You are likely in smoke mode (multiview is points‑only). See `c_core_interactions/04_view_layout_live_snapshots_small_multiples`.

---

## Next steps

- Learn navigation + multiview properly: `c_core_interactions/index`
- Learn the language of the UI: `a_orientation/04_ui_glossary_terminology`
