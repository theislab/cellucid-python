# Quick tour (60 seconds)

**Audience:** first-time users (wet lab + computational)  
**Time:** 1–3 minutes  
**What you’ll do:** enter → move → color → keep a snapshot → highlight → export/save

:::{tip}
If you prefer to learn by definition-first, read {doc}`04_ui_glossary_terminology` before this tour.
:::

---

## What appears on ordinary startup

```{figure} ../../../_static/screenshots/web_app/welcome-startup.png
:alt: The Cellucid welcome overlay as a first visit opens it, with Choose a dataset and Learn more, over the Suo atlas loaded behind it.
:width: 1440px

Every ordinary bundled or catalog startup opens the welcome overlay.
The Suo dataset remains the sample
catalog default and can finish loading behind it;
dismissing the overlay is not remembered for the next ordinary load. A
startup whose address already names its data (a Python or R viewer, a
`cellucid serve` server, or a `?data=` link) skips the overlay.
```

Choose **Choose a dataset** to open the data controls, or press {kbd}`Escape`
to inspect the loaded default dataset first.

```{figure} ../../../_static/screenshots/web_app/startup-loaded-build.png
:alt: The Suo atlas loaded after the welcome overlay is dismissed, with every sidebar section collapsed so the footer shows, ending in the Version and Build line.
:width: 1440px

Suo loaded after the welcome overlay is dismissed, on the light **Studio**
backdrop, with every sidebar section collapsed so the footer shows: its last
line is the version and the exact web build. Yours will read a later build than this capture — the value changes
with every web release, which is why bug reports should quote your own.
```

---

## One screenshot “map” (recommended)

This tour references UI areas by name (the app labels). One screenshot makes
the rest of the docs much easier to follow.

```{figure} ../../../_static/screenshots/web_app/window-orientation-map.png
:alt: The whole application window with the Pancreas sample loaded, the sidebar showing Coloring & filtering with the embedding coloured by clusters, and the eight-row legend beside the matching colours on the canvas.
:width: 1440px

The whole app on the Pancreas sample. **Left**: one scrollable sidebar of
collapsible accordions — every control in Cellucid is in there. **Right**: the
canvas. The legend swatches and the point colours are the same palette, which
is the single most useful correspondence to internalise: the sidebar and the
picture are always describing the same thing.
```

```{figure} ../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: The whole Cellucid window with the Pancreas sample coloured by cell_type: the sidebar shows Coloring and filtering with the categorical field selected and its eight-row legend, and the 3D embedding fills the canvas in the matching colours.
:width: 1440px

The same sample coloured by `cell_type`, with the other sections collapsed so
**Coloring & filtering** sits at the top: the three field selectors, the
legend's eight rows with their cell counts, and **Filters** reading *Showing
all 3,696 cells*. The published counts of every built-in sample are listed in
{doc}`../b_data_loading/12_sample_dataset_provenance`.
```

---

## Fast path (click-by-click)

### 1) Confirm or replace the default dataset

1) After dismissing the welcome overlay, inspect the **Session** accordion.
2) Keep the loaded Suo sample, or choose another entry under **Sample
   dataset**. **Pancreatic endocrinogenesis (scVelo)** is the standard
   1D/2D/3D velocity example; see
   {doc}`../b_data_loading/10_standard_pancreas_dataset`.
3) Wait for the **Dataset** box at the top of **Session** to name the dataset
   and its cells and genes.
   Every sample opens each dimension on its default embedding; Pancreas also
   carries a 2D PCA embedding beside its 2D UMAP.

**What success looks like**

- The canvas shows points on the default light studio floor.
- The **Dataset** box shows non‑zero cell and gene counts.

:::{note}
You can also open your own data with **Open…** (**File…** for an `.h5ad`,
`.zarr.zip`, `.h5mu`, `.loom`, 10x `.h5` or R `.rds` / `.rda`; **Folder…** for
a Zarr folder, a Cell Ranger `outs/` folder or a MatrixMarket folder), or a
link or a running server through the **Link or server** bar, but those
workflows are documented in
{doc}`../b_data_loading/index`.
:::

### 2) Move the camera

1) Open the **Navigation** accordion.
2) Check the dimension-specific default (each view keeps its own mode):
   - 1D or 2D starts in **Planar**;
   - 3D starts in **Orbit**.
3) In Planar, drag to pan and use the wheel or trackpad to zoom.
4) In Orbit, drag to rotate, use the wheel or trackpad to zoom, and
   right-drag (or Shift-drag) to pan.

The app does not start camera motion by default. **Camera path** starts
collapsed, its top transport is absent until at least two keyframes exist, and
playback begins when you select **Play**. If you explicitly enable
**Autoplay**, a ready path starts immediately and starts again when that
setting is restored from a session.

### 3) Color by a field (metadata or gene expression)

1) Open **Coloring & filtering**.
2) Pick one:
   - **Categorical obs** (clusters, batch, sample),
   - **Continuous obs** (QC metrics, scores), or
   - **Gene expression** (type a gene name into the search box).
3) Look for the legend changing on the left.

```{figure} ../../../_static/screenshots/web_app/field-selectors-three-routes.png
:alt: The Color by heading with its (i), above the three field selectors stacked in order — Categorical obs holding clusters, Continuous obs holding None, and an empty Gene expression search box — each followed by four small icon buttons, with the pointer pressing the first.
:width: 504px

The three ways in, in the order they appear. Only one can be active at a time —
choosing in one clears the other two.
```

**What success looks like**

- The points change color.
- A legend appears/updates (categories or a color scale).

Colouring by a gene looks like this:

```{figure} ../../../_static/screenshots/web_app/window-color-by-gene-ins1.png
:alt: The whole window coloured by expression of the gene Ins1 on a Viridis scale, dark where the gene is absent and bright in one lobe of the embedding, with the continuous legend in the sidebar.
:width: 1440px

`Ins1` on the Pancreas sample. Dark means little or no expression, bright means
high — and the bright region is exactly the beta-cell arm of the embedding,
which is what you would want a marker gene to do.
```

### 4) Keep a snapshot (small multiple)

1) Open **Compare views** → **Views**.
2) Click **Keep view**.
3) You should see a new numbered badge/panel.

This is how you compare hypotheses side‑by‑side without losing your current view.

### 5) Make a highlight (selection)

1) Open **Highlighting**.
2) Choose a mode (e.g., **Lasso**).
3) Draw a selection; a highlight group should appear in the highlight list.

If you want a keyboard-first workflow, the UI supports:

- `Alt` + drag for area selection (and other shortcuts; see the in‑app “Keyboard shortcuts” accordion).

On a phone or tablet, turn on **Select mode** on the selection bar at the
bottom of the plot (or in the Highlighting panel) and draw with one finger
instead; two fingers still zoom and pan. The bar then confirms, cancels or
clears the selection without the menu.

### 6) Export or save what you did

Pick one:

- **Save** under **Session file** (Session accordion) to download a `.cellucid-session` bundle you can reopen later.
- **Figure export** to export a publication‑ready figure of the current view(s).

For a bug report, scroll to the sidebar footer and include its exact **Build**
value. It identifies the web files you exercised and is separate from the
Python package version.

---

## Troubleshooting (quick tour blockers)

### “The canvas is blank” / “WebGL2 required”

Go to {doc}`02_system_requirements`.

### “I can’t rotate/zoom/pan”

Go to {doc}`../c_core_interactions/06_troubleshooting_core_interactions` (focus, pointer lock, nav mode confusion).

### “Keep view does nothing”

You are likely in smoke mode (multiview is points‑only). See {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`.

---

## Next steps

- Learn navigation + multiview properly: {doc}`../c_core_interactions/index`
- Learn the language of the UI: {doc}`04_ui_glossary_terminology`
