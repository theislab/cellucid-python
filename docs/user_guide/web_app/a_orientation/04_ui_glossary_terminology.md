# UI glossary (terminology)

**Audience:** everyone (this is the shared vocabulary)  
**Time:** 10–20 minutes (skim now, return later)  
**What you’ll get:** exact meanings of UI terms so the rest of the docs feel predictable

---

## How to use this page

- If you’re new: read the **bold terms only**.
- If you’re debugging something subtle: read the “Scope / persistence” notes under each section.

:::{tip}
This glossary is written to match the app UI language, even when a “technical” term might be more precise. The goal is: when the docs say “view” or “snapshot”, you know exactly what screen element/state is meant.
:::

---

## UI map (one screenshot you’ll reuse everywhere)

```{figure} ../../../_static/screenshots/web_app/window-orientation-map.png
:alt: The whole application window with the Pancreas sample loaded, the sidebar showing Coloring & filtering with the embedding coloured by clusters, and the eight-row legend beside the matching colours on the canvas.
:width: 1440px

The vocabulary below maps onto this picture. **Sidebar** on the left,
**canvas** on the right, **legend** in the sidebar's *Legend* section, and the
**count line** in the *Filters* section of *Coloring & filtering* telling you
how many of the dataset's cells are currently drawn.
```

---

## Startup and build identity

- **Welcome overlay**: the modal shown on every ordinary bundled or catalog
  startup. Its dismissal is not remembered for the next ordinary load.
  It lists the sample datasets (the default one says it is opening, with its
  progress; choosing one opens it), **Open a file…** and **Paste a link**;
  {kbd}`Escape` dismisses it for the current load. A startup whose address already names its data (a Python
  or R viewer, a `cellucid serve` server, or a `?data=` link) skips it.
- **Suo**: the sample catalog default. It may finish loading behind the welcome
  overlay before you choose another source, and opens in 3D on its default
  embedding.
- **Default background**: **Studio (light)** in a fresh browser profile. A user-selected
  background remains an explicit preference (see [Background](#background)).
- **Version**: the release version printed first in the sidebar footer. The
  web app, the Python package and the R package share it and are released
  together, so one number names the three.
- **Build**: the web build identity printed after the version in the sidebar
  footer. Include this exact value in bug reports so maintainers can identify
  the HTML, JavaScript, CSS, and asset inventory you exercised; it changes on
  every deployment, the version only at a release.

```{figure} ../../../_static/screenshots/web_app/startup-loaded-build.png
:alt: The Suo atlas loaded after the welcome overlay is dismissed, with every sidebar section collapsed so the footer shows, ending in the Version and Build line.
:width: 1440px

Fresh startup after closing the welcome overlay: Suo is loaded on the light
**Studio** backdrop; with every section collapsed, the footer at the bottom of
the sidebar names the version and the build. The capture shows the build it was taken from, not the one you are
running — read yours off your own footer.
```

---

## Compact help (`i`) buttons

- **`i` help button**: opens supplemental guidance for the nearby control. It
  does not replace the control's visible action label or current-state text.
- Click the button, or focus it with {kbd}`Tab` and press {kbd}`Enter` or
  {kbd}`Space`, to open its help dialog. Only one help dialog opens at a time.
- Press {kbd}`Escape`, move focus outside the dialog, or click elsewhere to
  close it.

See {doc}`../o_accessibility_privacy_security/01_accessibility` for the complete
keyboard and screen-reader behavior.

---

## Render mode and maturity tags

- **`Render mode`**: the dropdown that decides how points are drawn. It has
  exactly two options: `Points`, the default, which draws one point per cell;
  and `Volumetric smoke cloud (beta)`, which draws a continuous density cloud,
  coloured by the active field, instead. **“Smoke mode”** in these pages always
  means the second one.
- **Maturity tag**: a small pill beside a control's label naming how finished
  the currently selected option is. It describes the option you chose, not the
  control it sits next to.
- **`Beta`**: the only maturity tag in the app. It sits beside `Render mode`
  and appears only while `Volumetric smoke cloud (beta)` is the selected mode;
  `Points` never wears it, and nothing else in Cellucid is tagged.
- Its hover text reads
  `Volumetric smoke cloud is in beta: it is slower than points, and figure export cannot reproduce it.`

See {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke` for
what alpha means in practice for that mode.

## Background

- **`Background`**: the **Visualization** dropdown that sets what is drawn
  behind the cells. Its twelve options sit in two groups, **Light** (`White`,
  `Grid (light)`, `Studio (light)`, `Atelier (light)`, `Abyss (light)`,
  `Aurora (light)`) and **Dark** (`Black`, `Grid (dark)`, `Studio (dark)`,
  `Nebula (dark)`, `Abyss (dark)`, `Aurora (dark)`); the two lists pair up
  line by line, so `Nebula` is `Atelier`'s dark twin:
  - `White` and `Black`: a plain colour, nothing drawn.
  - `Grid (light)` and `Grid (dark)`: a matplotlib-style grid box
    of ruled walls around the data.
  - `Studio (light)` (the default) and `Studio (dark)`: a calm stage on a neutral grey or
    charcoal. A 3D view in Orbit or Free-fly stands on a **floor**: one plane
    just under the lowest cell, ruled with hairlines that fade out long before
    its edge and with a soft pool of shade (light) or light (dark) under the
    data. It is the same ruled floor from every angle, straight down
    included; it fades only where the camera sees it nearly edge-on, and once
    the camera has gone below it the stage is plain, as under a real floor. A
    1D or 2D view, and any view in Planar navigation, sits on a **dot
    surface** instead, like dotted paper behind the cells. Which of the two a
    view shows follows its embedding's dimension and its navigation mode,
    never the camera angle, so rotating a view never swaps them. Rules and dots
    keep the same on-screen size at every zoom and screen density, and their
    spacing steps by fives as you zoom, crossfading so the pattern never
    crowds or empties.
  - The **looks** are painted, three-dimensional settings. Each turns with
    the camera like the sky around a room, and a nearer layer (a dust shell,
    marine snow, a sea of clouds, the studio floor) slides against it as you
    orbit, pan or zoom, so the stage has depth in 3D and in 2D:
    - `Atelier (light)`: an endless paper studio. The cells cast a soft
      contact shadow on the floor (in 2D, on the paper behind them) that
      follows filters and representation changes; seen from below the floor,
      the paper glows with the shadow showing through.
    - `Nebula (dark)`: deep space, a nebula and a star field whose colours are
      taken from the colours the cells are coloured by.
    - `Abyss (light)` and `Abyss (dark)`: sunlit shallows over a caustic-lit
      sea floor, and the deep ocean with light shafts and glowing marine snow.
    - `Aurora (light)` and `Aurora (dark)`: a dawn sky with clouds and sun
      haze, and a night sky with aurora curtains over a sea of clouds.

    The cells stay the brightest (dark looks) or darkest (light looks) things
    on the stage: every look keeps its large areas outside the lightness range
    of the category palettes, and stars behind the cells are dimmed. A look is
    prepared in the background the first time you choose it (a fraction of a
    second; the stage shows the look's base colour meanwhile, then the look
    fades in) and costs nothing while the scene is still. The sky is prepared
    coarse first and sharpened in steps while the graphics processor keeps up,
    so it appears quickly on every machine and reaches full detail on a
    capable one. With the system's *reduce motion* setting it appears without
    the fade.
- The choice applies to every view of a grid compare, is remembered in the
  browser, and is saved in sessions. Figure export reproduces the grid box,
  the studio backdrop and the looks; see
  {doc}`../k_figure_export/06_edge_cases`.

---

## Data model terms

- **Dataset**: the loaded data bundle (points + metadata + optional gene expression + optional connectivities/vector fields).
- **Cell / point**: one row/observation rendered as a point. (The UI often says “points” because the viewer is generic.)
- **Gene / feature**: one variable/feature; used for gene expression coloring and some analyses.
- **Embedding** (also *representation*): coordinates for points in 1D/2D/3D
  (UMAP/tSNE/PCA/etc). A dataset may carry several embeddings of one dimension
  (the Pancreas sample has a 2D UMAP and a 2D PCA view); each view shows one of
  them, and each dimension has one **default embedding**, which keys `1`, `2`
  and `3` show first.
- **Data mapping**: which stored arrays the dataset is drawn from — the
  embeddings that are turned on and the default of each dimension; under
  **Genes** the modality (of a MuData or a multi-assay object), the expression
  matrix, gene names, values and highly variable genes; the neighbour graph,
  the latent space and the vector fields. Set under **Session › Data mapping**
  ({doc}`../b_data_loading/13_data_mapping`).
- **Best guess**: the default embeddings a file without a UMAP opens on, marked
  *guessed* in the Data mapping checklist.
- **View embedding**: the embedding one view shows, chosen under **Compare
  views → View N shows**; the view's row then reads e.g. `t-SNE · 2D`.
- **Dimension (1D / 2D / 3D)**: the dimension of the embedding a view shows.
- **Connectivity**: a neighbor graph (edges). When enabled, edges are drawn between visible points only.
- **Vector field / velocity overlay**: an optional overlay that visualizes per-cell vectors, bound to the embeddings it was computed on.

---

## Data and loading terms

Every word here appears on a button, a label, or a page in
{doc}`../b_data_loading/index`. Read the short definition; the linked page owns
the detail.

**What the Session panel calls things** (**`Dataset`** comes first; the next
three sit under **`Open data`**)

- **`Dataset`**: the box at the top of the panel that describes the dataset
  that is open, and only that one: its name, its cells and genes, its format
  and host (`h5ad file · this computer`), its description, how much has been
  read for data read in place from a host, and a collapsed **`Details`** list.
  With nothing open it says "No dataset open"
  ({doc}`../b_data_loading/01_loading_options_overview`).
- **`Sample dataset`**: the dropdown listing the catalog Cellucid was
  configured to read. On www.cellucid.com that is the five built-in samples
  ({doc}`../b_data_loading/12_sample_dataset_provenance`); on a copy you host
  yourself it is whatever catalog you pointed it at. While a linked catalog
  that holds several datasets is connected, it lists that catalog's datasets
  instead and its label reads "Dataset from <host>". A dataset opened from a
  file, a single-dataset link or a server is not listed there; while one is
  open the dropdown shows "—".
- **`File or folder`**: the **`Open…`** button, which reads a file or folder
  off your own computer. Its menu has **`File…`** (an `.h5ad`, `.zarr.zip`,
  `.h5mu`, `.loom`, 10x `.h5`, or an R `.rds` / `.rda` / `.RData`) and
  **`Folder…`** (a `.zarr` store folder, a Cell Ranger `outs/` folder or a
  MatrixMarket folder): a browser's file dialog picks files or a folder,
  never both. Dropping a file or folder
  anywhere on the page does the same. Nothing is uploaded.
- **`Link or server`**: one address bar with an **`Open`** button (Enter works
  too). A port (`8765`), `localhost:8765` or `host:port` connects to a running
  `cellucid serve`; a link or id of data that already lives somewhere (a
  CELLxGENE dataset, a GEO series `GSE…`, a Zenodo or Figshare record, a
  Hugging Face, Dropbox, Google Drive or OneDrive link, the `owner/repo/path`
  shorthand of a public GitHub folder, or a web address of a file, a Zarr
  store folder or a dataset catalog folder) opens it. While
  something is connected, one line under the bar reads **`Server`** or
  **`Link`**, the address and **`Disconnect`**
  ({doc}`../b_data_loading/14_open_from_link`).
- **`Session file`**: the **`Save`** and **`Load`** buttons, which download
  and restore a `.cellucid-session` bundle.

**What the data itself is called**

- **Cellucid-ready AnnData**: the AnnData `prepare()` (Python) or
  `cellucid_prepare()` (R) writes — a `.zarr` folder by default, or a
  `.zarr.zip` or `.h5ad` — with a gene-major expression matrix and its name,
  citation and defaults in `uns['cellucid']`. Any AnnData opens; this one
  reads fastest.
  See {doc}`../b_data_loading/07_file_format_expectations`.
- **Dataset catalog**: a folder holding a `datasets.json` file plus the AnnData
  files it lists. Opening the folder's address offers its datasets in the
  dataset picker.
- **`datasets.json`**: the catalog file, listing each dataset's id, its file's
  path inside the folder, and optionally a default session.
- **Dataset id**: the stable identifier that decides whether two loads are “the
  same dataset” for sessions and annotation.
  See {doc}`../b_data_loading/06_dataset_identity_why_it_matters`.
- **AnnData**: the single-cell object Scanpy works with. Its parts that matter
  to Cellucid are `obs` (one row per cell — the metadata you colour by), `var`
  (one row per gene), `obsm` (per-cell coordinate arrays, where embeddings and
  vector fields live), and `obsp` (per-cell-pair matrices, where the neighbor
  graph lives).
- **H5AD**: one AnnData saved as a single HDF5 file (`.h5ad`).
- **Zarr store**: one AnnData saved as a directory of chunked arrays (Zarr v2
  or v3). The browser reads the folder itself (**Open… › Folder…**, or a
  drop) or zipped as one `.zarr.zip` (**Open… › File…**), a link can name a
  store folder on a web host, and
  the Python server and Jupyter serve the directory itself.
- **Range request**: how the viewer reads an `.h5ad` or Zarr store (and the
  other HDF5 formats: `.h5mu`, `.loom`, 10x `.h5`), opened in the browser or
  served by `cellucid serve` and `show_anndata()` — only the bytes a view
  needs are read, so the expression matrix is never held in memory whole. An
  R object (`.rds`, `.rda`) or a MatrixMarket folder is read once from start
  to end instead.
- **Lazy loading**: gene expression is fetched only when you ask for a gene,
  rather than up front. Which paths are truly lazy is tabulated in
  {doc}`../b_data_loading/01_loading_options_overview`.

---

## Fields, coloring, and legends

- **Field**: something you can color by.
  - **Categorical obs**: discrete labels (clusters, batch, sample).
  - **Continuous obs**: numeric values (QC metrics, scores).
  - **Gene expression / var**: expression of a selected gene/feature.
- **Active field**: the currently selected field that drives coloring (per view).
- **Legend**: the UI element that explains the current coloring.
  - Categorical legend shows categories and their colors.
  - Continuous legend shows a numeric color scale.
- **Color-by**: the act of selecting a field to drive coloring.

The two legends look nothing alike, and telling them apart at a glance is the
fastest way to know what kind of field you are looking at:

```{figure} ../../../_static/screenshots/web_app/legend-categorical-clusters.png
:alt: The categorical legend for clusters: the Show all and Hide all buttons above eight rows, each holding a checkbox, a colour swatch, a label, a cell count, a rename pencil and a delete bin.
:width: 488px

**Categorical legend** — checkboxes and swatches.
```

```{figure} ../../../_static/screenshots/web_app/legend-continuous-sscore.png
:alt: The continuous legend for S_score: the Color scale (Viridis) bar with its two end values, Log color scale off and Rescale to range on, each with its (i), and below a rule the Visible range block with Live filtering on, Minimum and Maximum sliders, and Filter and Reset.
:width: 488px

**Continuous legend** — a colour bar and Min/Max sliders. A gene uses this one
too.
```

Scope note:

- In multiview, the “active field” is typically **per view** (each view can show a different coloring).

---

## Views, snapshots, and “small multiples”

These terms are the most important for understanding Cellucid.

- **View**: one panel that renders the dataset with its own state (embedding, camera, navigation mode, coloring, filters); highlights belong to cells and show in every view.
- **Live view**: your default working view (often labeled “All cells” in the view rows).
- **Snapshot view** (also “kept view”): a view created by clicking **Keep view**; used to compare multiple states side-by-side.
- **Multiview**: the overall feature of having multiple views/snapshots.
- **View layout**
  - **Grid compare**: shows all views in a grid.
  - **Edit selected view**: shows only the active view so you can edit it precisely.
- **Active view**: the view you are editing/configuring (the row selected in **Compare views**).
- **Focused view**: in grid mode, the view under your last click; this determines which view the camera controls target.
- **View row**: the row representing a view in **Compare views**: its number, what it shows (`UMAP · 3D`), its field, `linked`, its navigation mode (`Orbit`/`Pan`/`Fly`) and `×`.
- **Link cameras**: views that show the same embedding share one camera and navigation mode; views of different embeddings always keep their own.

Scope note:

- Cellucid stores camera state and navigation mode **per view** and switches to them when you focus a different view.

---

## Navigation and camera terms

- **Orbit**: rotate around an anchor (best for 3D).
- **Planar**: pan/zoom as if looking at a flat map (best for 2D).
- **Free-fly**: first-person navigation (immersive; supports pointer lock).
- **Pointer lock / Capture pointer**: a browser feature that hides the cursor and reports raw mouse movement (needed for FPS-style looking).
- **Orbit anchor**: a visual compass/anchor indicator shown in orbit mode (optional).
- **Reset camera**: returns the camera to a default framing (and can reset related UI controls depending on how you trigger it).
- **Dimension-specific navigation default**: a view showing a 1D or 2D
  embedding starts in Planar; a 3D one in Orbit. A mode chosen for a view is
  remembered for that view in that dimension and takes precedence.
- **Camera path**: a collapsed-at-start panel for user-created keyframes. Its
  transport appears only after two valid keyframes. Playback is stopped by
  default; explicit **Autoplay** starts a ready path immediately and after a
  successful session restore.

---

## Highlighting and selection terms

- **Tool** (under **Highlighting › Selection**): the active interaction tool for selecting points (e.g., lasso).
- **Highlight group**: a named set of highlighted cells (used for comparison/analysis).
- **Highlight page**: a collection of highlight groups (useful when you want multiple “sets” of groups).

Scope note:

- Highlights can be global or per-view depending on the feature; always check the page you’re using (highlighting is documented under {doc}`../f_highlighting_selection/index`).

---

## Analysis and export terms

- **Analysis**: computations driven by the current dataset state (often highlights/groups).
- **Figure export**: export an image/vector figure of the current view(s) suitable for papers.
- **Save / Load** (under **Session file**): save/restore a `.cellucid-session` bundle that captures the application state.

---

## Scope & persistence vocabulary (used throughout the docs)

- **Global**: one value shared across the entire app/dataset (not specific to a view).
- **Per view**: stored separately for each view/snapshot.
- **Session bundle**: a saved `.cellucid-session` file (explicit user action).
- **Local storage**: persistent browser storage used for small preferences (e.g., theme, background, and the antialiasing preference under the key `cellucid_antialias`).
- **Session storage**: storage that is cleared when the tab closes (used for sensitive tokens like GitHub OAuth in community annotation).
- **Device preference**: a setting that describes the machine you are viewing on rather than the picture you made. A session bundle does not record it, and loading a session never changes it. `Antialiasing` is the only one; the theme, by contrast, travels with a session.
- **Automatic (`auto`)**: a stored preference that is a deferral rather than a value. The antialiasing preference has three states — `auto`, `on`, and `off` — and `auto` is the default, and what an absent key means. Cellucid resolves it from the open dataset's cell count and the view's dimension (on below 5,000,000 cells in a 3D view and below 200,000 in a 1D or 2D view, off at or above) and decides again every time a dataset is opened or the view changes dimension. Ticking or unticking the checkbox replaces `auto` with your choice, permanently and in both directions. This is also the second reason a session does not record the setting: a checkbox cannot express `auto`, so restoring a captured tick would turn "let the dataset decide" into a choice you never made.

See {doc}`../l_sessions_sharing/02_what_gets_saved_and_restored` for what a
session bundle carries.

---

## Next steps

- Choose a workflow: {doc}`05_which_workflow_is_for_me_decision_tree`
- Learn navigation + multiview: {doc}`../c_core_interactions/index`
