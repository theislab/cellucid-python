# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.9.2]
<!-- CELLUCID_VERSION -->

Version 0.9.2 is a patch release on 0.9.1.

### Added

- **Web app, data loading**: drop one `.h5ad` file, one `.zarr.zip` archive or
  one prepared export folder onto the viewer to open it. An outline and a
  *Drop to open* card cover the viewer while you drag; the drop is read like
  the matching button and replaces any load in progress. Anything else is
  refused with a message saying what can be dropped, and a file dropped
  anywhere on the page, even while Cellucid is starting, never makes the
  browser leave it.
- **Web app, cell info**: **Show cell info** (Coloring & filtering, on by
  default) shows a card for the cell under the pointer: its value on the
  field the view is coloured by, where that sits on the colour scale, and for
  a category how many of its cells are visible in that view (what selecting it
  would pick); **Show recent fields** adds the cell's labels in other loaded
  fields. It works in every view of a grid; press and hold on a touch screen.
  A lookup never loads a field or redraws the plot (1M cells: 1.4 ms mean).
- **Web app, citing and sharing**: a **Cite** dialog (APA and BibTeX, read
  from `CITATION.cff`, with a placeholder note until the preprint has a DOI),
  and a **Share** dialog for datasets that can be opened from a link: the
  link, an "Open in Cellucid" README badge and an embed snippet.
- **Docs**: new pages on cell info, share links, badges and embeds, how to
  visualize single-cell data interactively, and how Cellucid relates to other
  viewers.
- **Citation**: a placeholder citation for the forthcoming Cellucid bioRxiv
  preprint, the same in every repository. It is written once, as the
  `preferred-citation` of `theislab/cellucid`'s `CITATION.cff`, and copied from
  there to this repository's `CITATION.cff`, the README, and the new
  [Citing Cellucid](https://cellucid.readthedocs.io/en/latest/citing.html) page;
  the test suite fails if a copy drifts. Until the preprint has a DOI, every copy
  says it is a placeholder.
- **Web app**: cells can be highlighted on a touch screen. Every selection tool
  started only from an `Alt` press, which a phone or tablet does not have, so touch
  users could not select at all. A **Select mode** switch in the Highlighting panel
  now stands in for a held `Alt`: while it is on, a one-finger tap or drag runs the
  active tool through the same path an `Alt` press takes, and its Intersect / Add /
  Subtract buttons stand in for `Shift` and `Ctrl`. Two fingers always pinch to zoom
  and, newly, drag to pan; a second finger landing mid-selection abandons that
  gesture rather than committing it. `Esc` with nothing in progress leaves the mode.
  Mouse and keyboard behaviour is unchanged while it is off.

- **Web app**: a **selection bar** at the bottom of the plot carries Select mode's
  controls, so a selection can be made, combined, confirmed, cancelled and cleared
  without opening the menu, which covers the whole screen on a phone. It shows the
  same switch and Intersect / Add / Subtract choice as the panel, the active
  tool's cell count with its Undo, Cancel and Confirm while a selection is in
  progress, the active page's highlighted count with Clear, and a one-line note
  when a tap on empty plot changed nothing. Every control on it is the panel's
  control shown again, so the two always agree. It is always shown on touch
  screens (just the switch while the mode is off), shown to a mouse while the
  mode is on, centred in the part of the plot the sidebar leaves uncovered, and
  notifications stack above it.
- **Web app, backgrounds**: two new viewer backgrounds, **Studio (light)** and
  **Studio (dark)** (neutral grey and charcoal): a 3D view stands on a floor
  of constant-width hairlines that fades out before its edge, the same ruled
  floor from every angle above it (below it, the stage is plain); 1D, 2D and
  Planar views sit in front of dotted paper. Which one a view shows follows
  its dimension and navigation mode, never the camera angle, so rotating a 3D
  view never swaps it to dots. The pattern rescales with zoom without popping, and
  figure export reproduces it in SVG and PNG. **Studio (light) is the new
  default** (the sample datasets open on it too), and the Background menu is
  grouped into Light and Dark.
- **Web app, shader quality**: a fourth quality, **Sculpted (soft 3D
  light)**. Every cell is a soft-lit, glossy ball — a wrapped key light, a
  sky-and-ground ambient, a broad gloss instead of a hard white highlight, and
  darkening toward the edge that keeps touching cells apart — computed in
  linear light so the centre of each ball is exactly its legend colour. With
  antialiasing on, its edges are smooth in 3D views too, without the
  background showing through nearer cells; highlighted cells get the same 3D
  ring as `Full`, held still. It costs about what `Full` does (3–5% more in
  3D views on an Apple M1 Pro, the same in 2D). `Full` stays the default.
- `var_gene_name_columns=` on `prepare()`: further `var` columns offered as
  alternative gene names (for example symbols beside an Ensembl index).
- `embedding_keys=` (`--embedding-1d/2d/3d KEY[:c0,c1]`), `expression_source=`
  (`--expression X|raw|layer:NAME`) and `--gene-id-column` for serving AnnData.
- `cellucid.data_mapping`: the data-mapping model shared with the web app and
  the R package, with a conformance fixture
  (`tests/fixtures/data_mapping_conformance.json`) and a cross-writer byte
  parity fixture (`tests/fixtures/writer_parity_v3.json`).
- A sixth reconciled surface, **Embedding**: the writer proves `embeddings/`
  holds exactly the files `dataset_identity.json` declares.
- `cellucid serve`, `serve_anndata()` and `show_anndata()` with AnnData publish
  the data mapping they serve at `/_cellucid/mapping/candidates.json` (as
  `defaults`, beside every candidate `obsm`, `obsp`, layer and `var` column,
  read from metadata only), and `AnnDataAdapter` gains `get_data_mapping()` and
  `get_mapping_candidates()`. The web app reads it to record the mapping a
  session was saved with.
- **Web app, data mapping**: a **Data mapping** group at the bottom of the
  Session panel lists every embedding the dataset can be shown with as a
  checklist: turn any number on, of any mix of 1D, 2D and 3D, and mark one per
  dimension as its **Default**. For an `.h5ad` or Zarr ZIP opened in the
  browser every numeric `obsm` array is listed, a wider one at any selection
  of 1 to 3 of its columns (`X_pca` columns 2 and 3 as `PCA 3,4`), as many
  selections as you like; a prepared dataset lists the embeddings it was
  exported with, all on. Turning an embedding on loads nothing until a view
  shows it. Choices are staged and applied as one step that loads and checks
  every embedding a view must switch to first; the neighbour graph, gene
  names, expression matrix, latent space and vector fields are listed
  read-only. The data mapping's embeddings are now `{enabled, defaults}`, in
  the web app, `cellucid.data_mapping` and the R package's conformance fixture.
- **Web app, several embeddings at once**: every view shows one of the
  embeddings that are on, of any dimension — two 3D embeddings side by side, a
  3D one beside a 2D one, a t-SNE beside a UMAP — chosen in **View N shows**
  under **Compare views**, grouped by dimension with the default first. Keys
  `1`, `2` and `3` show that dimension's default in the selected view and,
  pressed again, its next embedding ("2D · t-SNE (2 of 3)"). Each view row
  names what its view shows (`UMAP · 3D`), its field, its navigation mode, and
  whether its camera is linked. Views of one embedding share its coordinates
  and GPU buffers; coordinates no view shows are kept in a cache bounded by a
  quarter of the device memory (512 MiB to 2 GiB) and released least recently
  shown first. Up to four views (the live view and three kept views).
- **Web app, per-view cameras**: each view owns its camera and navigation
  mode. **Navigation** is its own accordion and edits the selected view; its
  **Applies to** row names that view and the views linked to it. **Link
  cameras** (on by default) links only views that show the same embedding, and
  is off and disabled while no two views do.
- **Web app, rendering**: points are stored on the GPU in the draw order of
  their spatial tree (by LOD band, then tree leaf, then cell). A reduced level
  of detail draws a prefix of the stored points, frustum culling draws whole
  tree-node ranges (one multi-draw call, at most 256 ranges, with no per-point
  index upload), and full detail draws faster because neighbouring points are
  neighbours in memory: 10M points 18 ms against 68 ms (Chromium), 15 against
  57 (WebKit), 14 against 57 (Firefox); a zoomed-in culled frame 7-8 ms against
  21-31 ms. LOD no longer keeps a second copy of the positions: about 200 MB
  of GPU memory at 10M points instead of 315. The order is adopted a few
  milliseconds a frame when a tree arrives, so it never stalls a frame, and a
  view switched back to a representation it showed keeps its colours without
  an upload. Translucent points in 2D overlap in tree order, the same in every
  view, every level of detail and figure export.
- **Web app**: a notification can carry one action button; the best-guess
  notification uses it (**Show**) to open the Data mapping section.
- **Web app, usage analytics**: `data_load` reports how the embeddings were
  chosen (`mapping_origin`: `default` or `best_guess`) and how many the
  dataset could show (`embedding_candidates`: `1`, `2-3`, `4-9`, `10+`);
  how many it opened with turned on (`enabled_representations`); `feature_used`
  counts `data_mapping` (`embedding`, `embedding_columns`, `default`) and
  `dimension` when the picker switches a view to another dimension. No key or
  embedding name is ever sent.

### Changed

- **Export format version 3.** One export can now carry several embeddings per
  dimension and several sets of gene names. This is a clean break: exports in
  the earlier format are not read, and format-3 exports need a web app that
  reads format 3. Re-export with `prepare()` (or `cellucid_prepare()` from the
  matching R release, which writes the same bytes). `dataset_identity.json`
  has `version: 3`, and `embeddings` is now `{default_dimension, defaults, items}`: every exported
  embedding is an item with a `label`, a `dimension`, the exact payload `file`
  `embeddings/{index}_{dimension}d.bin(.gz)` and its `source` (`key` and the
  0-based `columns` it was read from). Items are ordered by dimension, then id;
  `defaults` names the item each dimension shows first. The `points_*d.bin`
  files at the export root are gone.
- **Vector fields bind to embedding ids.** A field's `files` are keyed by
  embedding id, `vectors/{fieldIndex}_{embeddingIndex}.bin(.gz)`, scaled by that
  embedding's normalization factor. `basis`, `available_dimensions` and the
  per-field `default_dimension` are removed; `label` is the field id.
- **Manifests are `compact_v2`.** A categorical field's element `[5]` is
  `centroidsByEmbedding`, keyed by embedding id. `var_manifest.json` replaces
  `var_gene_id_column` with `gene_name_sets`: the default set (from
  `var_gene_id_column` or the var index) and any alternative sets, stored as
  metadata only, so a reader can switch gene names without reading expression
  again.
- **`prepare()` API:** `X_umap_1d`, `X_umap_2d` and `X_umap_3d` are replaced by
  `embeddings={id: array or Embedding}` (required) and `embedding_defaults`
  (required for every dimension with more than one embedding). The new
  `cellucid.Embedding(values, columns=None, label=None, source_key=None)` reads
  chosen 0-based columns of a wide array, such as `X_pca` columns `[0, 1]`.
  `vector_fields` is now `{field: {embedding id: array}}`.
- **Gene names are made unique instead of refused.** Repeated names become
  `MATR3`, `MATR3-1`, ... by the rule anndata's `var_names_make_unique()` uses,
  and an empty name takes its row's var index name. The same rule runs in the
  Python writer, the R writer and the AnnData server, and what changed is
  printed and recorded in `var_manifest.json`. `gene_identifiers` matches the
  unique names.
- **`cellucid serve` / `serve_anndata()` / `show_anndata()` for AnnData** serve
  identity v3. `serve_connectivity` and `--connectivity` are replaced by
  `knn_key` and `--knn KEY`: the named `obsp` matrix is read and checked when it
  is named, and one that is not symmetric (such as `distances`) is refused with
  that reason instead of being symmetrized.
- When an object has no `X_umap_1d`/`X_umap_2d`/`X_umap_3d` and no 1- to
  3-column `X_umap`, the server no longer stops with rename hints. It makes a
  deterministic best guess from the other `obsm` arrays (t-SNE, force-directed
  layouts, PHATE, spatial coordinates and diffusion maps are preferred over
  latent spaces such as `X_pca`) and reports it; when nothing can be plotted it
  says what `obsm` holds and why each entry was not used.
- Sessions: `datasetFingerprint` must carry `cellOrder` (now naming the
  `embedding` its digest was taken on) and the new `dataMapping`. User-defined
  categorical fields no longer store `centroidsByDim` or `normalizedDims`.
  `apply_cellucid_session_to_anndata()` checks the mapping against the AnnData,
  compares `varCount` on the mapped gene axis, and resolves var-sourced fields
  through the mapped gene names and expression matrix.

- **Web app, data**: the viewer reads export format 3 and only format 3 — from a
  prepared folder, the sample datasets, GitHub, a server, and the notebook —
  and opens `.h5ad` and Zarr files as the same format. Every embedding an
  export carries is loaded only when a view shows it; each dimension opens on
  the embedding `defaults` names. An `.h5ad` or Zarr file's vector field now
  reports its id as its label (`velocity_umap`, not `Velocity (UMAP)`), as the
  format requires, and the reader now refuses any export whose vector-field
  label is not its id, from every source.
- **Web app, data**: an `.h5ad` or Zarr ZIP without a UMAP is no longer
  refused. It opens on the same deterministic best guess `cellucid serve`
  makes (t-SNE, force-directed layouts, PHATE, spatial coordinates, …, and last
  the first two columns of a wide array such as `X_pca`), marked *guessed*,
  with a notification naming what is shown. Only a file with nothing drawable
  is refused, with the message `cellucid serve` prints: "Cellucid found
  nothing to plot…", listing what `obsm` holds and why each entry was not used
  (the error code is `NO_PLOTTABLE_EMBEDDING`). Every embedding is checked
  before it is drawn — NaN, infinity, a value outside 32-bit float range, or
  no coordinate variation is refused naming the first cell — so a one-cell
  file can no longer be opened. A one-dimensional `obsm` array is read as one
  column, and a wide array is read only around the columns a view needs.
- **Web app, vector fields**: the overlay is offered by the embedding a view
  shows, not by its dimension: a field is drawn only in views showing the
  embedding it is bound to, and the info line names those embeddings
  ("Vector fields are bound to … Show one of them in this view to enable.").
- **Web app, sessions**: a session records its data mapping (the embeddings
  that are on and each dimension's default), the embedding its cell order was
  taken on, the embedding every view shows, every view's camera and navigation
  mode, and whether cameras are linked. Loading checks the dataset, its data
  mapping and the cell order of the embedding the session names before
  anything changes on screen, loading that embedding if another is shown — so
  a session saved in 2D now opens while 3D is shown, instead of asking you to
  switch back. A session whose embeddings this dataset can show turns them on
  before its views are restored, and one naming an embedding the dataset does
  not have is refused naming it ("This session shows obsm['X_tsne'] columns
  [0, 1], which this dataset does not have."). A
  session saved by an earlier release is refused with a message saying so.
  User-defined categorical fields no longer store their centroids; they are
  recomputed from the codes.

### Removed

- `vector_fields.missing_embedding_message()`, `resolve_embedding_keys()`,
  `embedding_declaration_hints()` and the other helpers of the one-UMAP-per-
  dimension rule, together with `AnnDataAdapter.get_points_binary()`.

### Fixed

- **Web app, analysis performance**: analysis no longer freezes the page on
  large selections, and every result is the same number it was. Medians,
  quartiles, ranks and sorted copies (Quick, Detailed summaries and tests, CSV
  export, Spearman correlation, Kruskal–Wallis, high-cardinality Wilcoxon in
  differential expression and Marker Genes, gene signatures) come from
  linear-time radix ordering instead of comparison sorts; correlation pairs its
  colour field by merging sorted cell indices instead of a per-cell map; box,
  violin and density plots of large pages compute their statistics on the
  analysis workers; a histogram hands the plot a few values per bin instead of
  every value; and a box plot draws at most 16 copies of an identical outlier
  (all on one spot) instead of one SVG element per cell. At one million cells:
  box plot of a zero-inflated gene 1.1 s → 0.12 s (longest main-thread task
  70 ms), histogram 0.17 s → 0.06 s, Detailed summary 0.10 s → 0.03 s,
  three-page Kruskal–Wallis 0.67 s → 0.09 s, Spearman 0.86 s → 0.14 s,
  high-cardinality marker statistics 0.33 s → 0.06 s per gene. The Plotly
  notification watcher no longer runs a subtree query for every element any
  plot inserts.
- **Web app, data loading**: when the connection to a server is lost, Cellucid
  clears only the dataset that came from it, and not while another dataset you
  chose is still loading; it used to clear over that choice. Recovering from a
  lost connection while its dataset was on screen also failed outright
  ("Remote dataset URL requires an active connection"), because the open
  dataset's address was re-read from the closed connection; it is now the
  address the dataset was loaded from. A server with one dataset that answers
  **Connect** after you applied a session no longer replaces that dataset and
  session; its dataset is listed in the picker instead.
- **Web app, data loading**: a dataset chosen while another is still loading
  now replaces it. At startup the data controls stayed disabled until the
  startup dataset had downloaded, and **Load State** did nothing at all in
  that window; afterwards, a file still being read published over a sample or
  **None** chosen after it. Now the controls open once the dataset catalog has
  been read, and every dataset change (the startup dataset, a sample, **None**,
  a file or folder, a server or GitHub repository with one dataset, the
  benchmark) supersedes the load in flight the moment you act: its downloads
  are cancelled, it never appears late, and the dataset picker shows what is
  loading. A connection that answers after a newer choice keeps its dataset in
  the picker instead of replacing that choice, and a session chosen during a
  load waits for that dataset, says so, and is then applied.
- **Web app, rendering**: switching from **Grid (light)** to **Grid (dark)**
  no longer leaves the scene redrawing on every frame (the grid fade only
  moved upwards, so it never reached the lower target). With reduced motion
  requested by the system, background fades and crossfades are instant.
- **Web app, rendering**: connectivity edges declare their per-vertex
  attribute at location 0. Left to the linker, Firefox on macOS could draw
  them without an enabled attribute 0 and emulate one per vertex.
- **Web app, no dataset**: after choosing **None**, Keep view, figure Export,
  Save and Connect repo are unavailable and say why instead of throwing;
  Differential expression, Gene signature and Marker genes show `Requires gene
  expression data.` (as on a dataset without genes) instead of a page error;
  Analyze performance and Copy situation report ask for data; and Load says no
  dataset is open instead of quoting `0 cells and 0 genes`.
- **Web app, data**: a dataset re-exported at the same address can no longer
  show values from its previous version. Gene and observation payloads were
  reused from the browser cache without asking the server, and payload files
  are numbered, not named: a re-export that added one gene left
  `var/3.values.u8.gz` holding another gene of exactly the same size, so the
  cached file drew the old gene's values under the new gene's name. Every
  dataset file — catalogs, manifests, payloads and saved views — is now
  revalidated with the server before use; an unchanged file costs a short
  "not modified" reply.
- **Web app, interface**: one calmer design. One font (Inter) in sentence case
  with a small type scale; soft borders; compact spacing. Every panel has
  short headings, and every explanation moved from the panel into an (i) tip
  that links to the exact docs page. Navigation is its own panel and shows
  which view it drives when several are open; the benchmark panel is called
  Benchmark; Visualization, Navigation and Compare views start collapsed; the
  Outlier filter sits under Legend; the logo follows the theme. Slider thumbs
  no longer overhang the rows around them, and the legend tick no longer
  touches its swatch.
- **Web app, rendering**: a hover lookup and any click on the page no longer
  redraw the scene for a second.
- **Web app, picking**: hovering and clicking now pick the dot drawn under the
  pointer, at the size it is drawn, instead of the first cell near the view
  ray, which could be a nearer cell beside the pointer. Where dots overlap the
  front-most wins; filtered-out cells and cells the current level of detail
  does not draw are never picked. Cell info, Alt+click selection, proximity
  and KNN all use it, each in the pane under the pointer.
- **Web app, selection**: in a grid of views, Alt+click and Alt+drag in a
  pane other than the focused one now select by that pane's field. They used
  the focused view's field, so a click in a pane colored by `cell_type`
  beside a focused view colored by a continuous field selected a value range
  of the wrong field. The category, the value range, the drag preview and its
  range label now all come from the pressed pane, among the cells it shows;
  view focus still does not move. A pane colored by no field now says so.
- **Web app, selection**: whether an Alt press can start, and the armed
  cursor, now follow the pane under the pointer instead of the focused view,
  so Alt+click works in a colored pane while the focused view has no field,
  and the panel asks for a field only when no pane on screen has one.
- **Web app, cell info**: the card no longer disappears for good. A key whose
  release the page never saw (a Cmd shortcut such as the macOS screenshot
  keys, or a key held while the window lost focus) read as camera motion
  forever, and so did a drag released outside the page; the card waits for
  the camera to settle, so it never came back. Held keys are now only the
  navigation keys, are dropped on blur, and Shift alone no longer counts as
  motion. The card no longer repeats the pane's name in a grid of views.
- **Web app, benchmark**: the presets are now 500K, 1M, 5M, 10M, 20M and 40M.
- **Web app, welcome**: the welcome dialog opens with no choice pre-selected;
  focus rests on the dialog and the first Tab reaches **Choose a dataset**.
  The sidebar footer links GitHub next to Docs and Cite.
- **Web app, onboarding**: closing the welcome dialog, by any path, no longer
  focuses the **Sample dataset** dropdown; **Choose a dataset** opens the
  **Session** panel, and from the keyboard moves focus to its heading.
- **Web app, data loading**: dataset catalogs, identities and manifests are
  revalidated on every load, so for a few minutes after a dataset is updated
  a stale cached copy can no longer clash with the new files (for example a
  sample view failing with a missing control).
- **Web app, pointer**: hover, click and selection read the camera as it is,
  not as the last frame drew it, so a pointer read right after the camera
  moved no longer lands on the previous picture.

- **Web app**: automatic antialiasing now accounts for the view's dimension.
  A 1D or 2D view has no depth test, so every point is shaded and antialiasing
  costs about ten times more there (1M cells: 33.8 ms with it, 3.4 ms without).
  It is now on below 5,000,000 cells in 3D and below 200,000 cells in 1D/2D,
  and the automatic choice is re-decided when the active view switches
  dimension. The status line names both (`Chosen automatically: 3,696 cells in
  3D.`); an explicit choice is never changed.
- **Web app**: less memory per cell, so larger datasets fit. Per-cell visibility
  is one byte (exactly 0 or 1) instead of four, the outlier quantiles are no
  longer copied into the state, categorical codes keep their stored width, and
  a kept view shares the visibility arrays copy-on-write instead of cloning
  them. Loaded fields are released by bytes rather than by count, so browsing
  genes no longer grows memory without bound. A loaded 3D dataset went from
  36 to 22 bytes per cell, and each kept view from 26 to about 6; 150M cells in
  2D and 120M in 3D load and render on a 16 GB machine. Pixels, filters and
  exports are byte-identical. The benchmark's pre-publication GPU estimate now
  follows the dimension (13 B/point in 2D, 17 in 3D) instead of a fixed 28.

- **Web app**: usage analytics rebuilt around one declared schema
  (`assets/js/analytics/schema.js`). It no longer sends dataset names, ids built
  from local file names, exact counts, error text or full page URLs (the page
  view carried `remote=` and `github=` addresses and their tokens); a dataset is
  described by count ranges, dimensions and writer version, with an id only for
  the public samples. Each load is reported once (it was reported twice), timed
  from the click or file pick (local loads excluded the file read), and reported
  when it fails before activation (such failures sent nothing). New events cover
  startup success and failure, feature adoption, and device classes. Every
  payload is validated against GA4's limits, and analytics is never awaited.

- **Web app**: faster, quieter full-detail rendering with identical pixels.
  Hovering the plot or typing into a sidebar field redrew the whole scene for a
  one-second settle window (237 full frames over two seconds of hover at 1M
  cells); only input that changes something now draws. In 2D views — drawn
  without depth testing — the point shaders use a single `discard`, and `Light`
  quality uses a no-depth program that blends a zero-alpha fragment instead of
  discarding: at 1M cells with antialiasing, `Light` 30.7 → 17.0 ms and `Full`
  33.8 → 30.6 ms per frame, byte-identical in 45 compared scenes. Orbit and pan
  inertia, smooth zoom and held-key acceleration now scale with elapsed time, so
  a 120 Hz display no longer moves the camera at double speed, and the render
  loop allocates about half as much per frame.

- **Web app**: with an ad or tracking blocker (uBlock Origin, Brave, AdGuard, iOS
  content blockers — anything shipping EasyPrivacy, which lists
  `/analytics/tracker.js`), Cellucid did not start and showed “Cellucid could not
  load its code…”. The analytics module was a static import, so its refused
  request aborted the whole module graph. Analytics is now loaded with a dynamic
  import and stays off when refused; the dataset-load vocabulary it shared with
  the controls moved to `app/data-load-methods.js`. No other application file
  matches the EasyList, EasyPrivacy, uBlock Origin or AdGuard lists.

- **Web app**: choosing volumetric smoke on a device without float blending
  (Safari on iPhone and iPad, many Android browsers) failed mid-build with
  “Smoke density unavailable: GPU smoke density requires EXT_float_blend for
  exact R32F accumulation.”, and restoring a smoke session there threw. Smoke's
  requirements are now checked up front; such a device gets a plain message
  (“Volumetric smoke isn't available on this device: …”) and stays in Points,
  sessions restore into Points, and other smoke errors no longer show internal
  detail.

- **Web app**: the volumetric smoke cloud is now **beta** rather than alpha. It
  is coloured, bounded and idle-aware after this release's rewrite; the tag beside
  `Render mode:` reads `Beta` in the calmer info tone, and its hover text no
  longer says the mode is still changing.

- **Web app**: analysis plots of large pages draw several times faster, with
  identical results. Plotly's violin density evaluated its kernel once per value
  at every curve point; a local, test-pinned patch to the vendored Plotly
  evaluates it once per run of equal values (expression data is mostly ties) and
  skips terms that are exactly zero, adding the rest in the same order. Box plots
  that draw points hand Plotly exact precomputed statistics; density curves are
  computed once (not twice), reused across cosmetic changes, and split across the
  worker pool when large. At 2,000,000 cells: a gene violin 3.7 → 0.67 s, a
  continuous box 1.15 → 0.26 s, a continuous density 1.39 → 0.56 s. Page values
  keep their own numeric type, the page-data cache is bounded at 256 MiB as well
  as 100 entries, and decodes go to the least-loaded worker.

- **Web app**: on browsers that composite in software (blocklisted GPUs, virtual
  machines, remote desktops), part of the plot vanished beside a floating panel
  whose content scrolled: Chromium culled the canvas under the panel's opaque
  scroller and drew the hole at its mirror image. Each floating window is now its
  own compositing surface.

- **Web app**: multiview and the data path are faster. Kept views drew their
  reduced Level-of-Detail levels through a gather from the full-detail streams;
  they now draw the live view's compact positions with their own compact colors,
  and share the live positions instead of copying them. Measured on an Apple M1
  Pro (Chromium, ANGLE Metal), 10,000,000 cells in 2D with LOD: two views
  15.9 → 6.9 ms a frame, four views 43.6 → 13.5 ms, and each kept view now adds
  38 MB of GPU memory instead of 114 MB. Hiding a category took 225 ms instead
  of 413 ms, a cached recolor 225 ms instead of 321 ms, and loading a prepared
  export 1.30 s instead of 1.43 s (30,000,000 cells: 2.89 s instead of 3.45 s).
  The documentation's claim that a four-view grid costs about as much as one
  view held only for large, overlapping points and is corrected.

- **Web app**: a panel's undock, dock, copy and close buttons no longer sit inside
  its header's open/close toggle. A `<summary>` is itself a button, so each was a
  control nested in a control — announced poorly by screen readers and flagged by
  accessibility checkers. They now stand beside the header, look and sit exactly
  where they did, are named after their panel ("Undock Visualization panel into a
  floating window"), and come just before the header in `Tab` order.

- **Web app**: multiview panes are laid out in the part of the window the sidebar
  leaves uncovered. The grid was cut from the whole canvas, so with the sidebar
  open the left pane was partly hidden under it (at 1440 px with two panes, 430 px
  of the left pane was visible against 720 px of the right one). Every pane now
  has the same size beside the open sidebar, and the layout follows at once when
  the sidebar is closed, opened, resized, or the window or its pixel density
  changes. One layout drives drawing, picking, hover, focus, every selection
  tool, projectiles, the pane title chips, centroid labels, figure export, and
  the selection bar, so what you click is always the pane you see. Wheel zoom
  now focuses the pane under the cursor, free-fly is centred in the uncovered
  area like Orbit and Planar, and the hidden sidebar is `inert`, so `Tab` no
  longer walks into it.

- **Web app**: the sidebar's ✕ works from the first frame. It was wired only once
  the first dataset had loaded, so during a slow first load a phone user, who
  sees the menu full-screen, tapped it and nothing happened.

- **Web app**: smaller layout fixes.
  - On a short screen (a phone on its side), startup notifications covered the
    welcome dialog's buttons; while it is open they now wait beneath it.
  - A long view title in a narrow multiview pane ran into the next pane's
    title; it now ends in an ellipsis inside its own pane, with the full name
    on hover.
  - The benchmark's Frame Time reading wrapped onto three lines. The second
    number now sits on its own line, labelled `CPU` for what it measures: the
    main-thread time spent issuing the frame, not GPU time.
  - Figure export's Background and Font selects now have visible labels, with
    the same trailing colon as every other label.

- **Web app**: in multiview, a finger dragging a pane now focuses that pane first,
  as a mouse press does, instead of moving whichever pane last had focus, and a
  touch drag saves the camera to the focused pane when cameras are unlocked.

- **Web app**: the volumetric smoke cloud is rebuilt and fixed.
  - It drew near-black soot whatever the cells were coloured by. The cloud now
    takes the active field's colours: each region shows the mean display colour
    of the cells there, lit by a key light and sky occlusion pre-computed once
    per build rather than by a shadow ray from every sample.
  - Switching dimension or recolouring (field, colormap, category colour, log
    scale) left the previous volume on screen. Both now rebuild it, and a 2D or
    1D embedding is spread into a soft slab instead of a sheet rays stepped over.
  - `Ray quality` changed almost nothing. It now sets samples per voxel, 0.5 to
    2 (`High`, 1.41, by default), and empty space is skipped at no cost.
  - The cloud shimmered and the GPU never idled: the ray-start jitter moved
    every frame and every frame was redrawn. The jitter is now fixed, and at
    `Animation speed` 0, which reads `Still`, a still view is not redrawn.
  - `Render resolution` reached 2.00x and `Noise detail` 256³, either of which
    could stall one frame long enough to lose the WebGL context. Render
    resolution now runs 0.25x–1.00x (default 0.50x) and a frame marches at most
    2560×1440 rays; noise detail is 32³–128³ (default 128³).
  - `Cloud density`, `Fine detail` and `Light absorption` printed readouts that
    changed with grid and noise resolution, so a recorded value did not
    reproduce a look. Every smoke slider now has one fixed scale.

  The documentation's control tables, recovery order and pipeline description
  are corrected (grid size was never the cause of context loss: the grid is
  capped at 128³, an 8 MiB volume), and the sample datasets' saved views now
  open with the new smoke defaults.

- **Web app, user-defined fields**: a categorical made with
  `Create Categorical`, a category merge or delete, or a field duplicate now
  shows its category labels, keeps them when a view switches dimension (that
  switch threw `Raw 2D positions are unavailable for centroid computation`),
  and gets them back from a saved session exactly where they were; restoring a
  session whose active field is one of these no longer fails with
  `Active categorical field is absent from its rendered selector`.
- **Web app, sessions**: a dataset without genes can be saved to a session;
  capture used to fail with `Dataset fingerprint requires the current var
  field inventory`. Such a dataset now has an empty gene inventory rather than
  none.
- **Web app, camera**: at a low frame rate (a large dataset with level of
  detail off) inertia, zoom and the grid fade now run in real time instead of
  slowing down below 20 fps; only a stall over a quarter second is cut short.
- **Web app, rendering**: switching to another window or tab no longer redraws
  the whole scene for a second; focus changes redraw only when they end a held
  key or a selection in progress.
- **Web app, touch**: on an iPhone or iPad, pressing and holding the plot no
  longer selects the whole page as text or opens the callout menu; the press
  reaches the viewer (Cell info's hold, Select mode, the camera). Safari reads
  only `-webkit-user-select`, so every non-selectable control now declares it
  too.
- **Web app, highlights**: a highlighted selection no longer keeps the whole
  scene redrawing every frame in `Light` and `Ultra-light` shader quality,
  where the ring does not animate; the view goes idle as it does without a
  selection. `Full` still pulses its ring.
- **Web app, fog**: points, highlight rings, connectivity edges, the orbit
  anchor and the velocity overlay share one fog model, evaluated once per
  point instead of per pixel, so every layer of a scene fades at the same rate.
  Fog now thins a point toward whatever is behind it instead of tinting it
  with a fixed fog colour, so fogged points fade correctly over a grid or a
  vignette; over a plain background the picture is unchanged. With fog off a
  dense 10M-cell view draws about 8% faster (Apple M1 Pro).
- **Web app, 2D views**: point edges in 1D and 2D views are smooth (a
  one-pixel antialiased rim) in `Full` and `Light`, instead of a jagged
  staircase that the antialiasing setting could not smooth; `Full` draws 2D
  views as cheaply as `Light`, and figure export draws them the same way.
- **Web app, rendering**: connectivity edges take about a third less GPU time
  (all 6.3M edges of the Suo sample: about 190 to 115–128 ms a frame) and the
  volumetric smoke cloud about 15% less (17.2 to 14.5 ms), on an Apple M1 Pro.
  Every shader now pins a vertex attribute to location 0, which Firefox
  otherwise emulates per vertex. Projectiles no longer allocate new vertex
  arrays every frame.
- **Web app, projectiles**: projectiles bounce off the **Studio** floor under a
  3D view, as they bounce off the grid's walls, and land on it instead of
  falling through. A projectile that had come to rest against a cell no longer
  falls through the grid's floor either.
- **Web app, rendering**: after a change the viewer draws two more frames
  and goes idle, instead of at least twelve frames and a full second; an
  adaptive LOD still stepping toward its level keeps the scene drawing on its
  own account, so nothing is cut short. A fast display no longer draws 60–120
  identical frames after every click, and where a frame takes seconds (a
  software renderer) every change costs seconds less.
- **Web app, smoke**: generating the cloud's noise volumes no longer queues
  GPU work faster than the GPU finishes it. Each batch of slices waits for the
  previous one and shrinks while the GPU lags, so on a slow or software
  renderer the page stays responsive while smoke mode prepares (before, a
  click there could wait half a minute behind the queued noise work).
- **Web app, selection tools**: while a lasso, proximity or KNN selection is
  drawn, a pointer move that leaves the previewed cells unchanged no longer
  re-uploads the highlight or redraws the scene; only a change does.
- **Web app, sessions**: a session saved while a view was coloured by a
  continuous field could not be loaded ("Session UI control … does not
  exist"): it recorded the legend's range sliders, whose ids are made afresh
  on every page. Controls with page-made ids are no longer recorded; the
  colour range and filter they show are saved with the field, as before.
- **Web app, shader quality**: `Full` shading is a little brighter on the
  side of each cell away from the light (ambient 0.5, was 0.4); the lit side
  still shows the colour itself. Highlight rings and centroid markers use the
  same lighting.
- **Web app, Studio backgrounds**: the floor's rules and the dots are lighter
  and finer (major rules no thicker than minor ones), so the stage reads as a
  calm, clean surface behind the data.
- **Web app, rendering**: the plot is drawn at exactly the screen pixels it
  covers. Its drawing buffer used to be rounded down from a fractional layout
  width (a plot beside the sidebar is often half a pixel wider), and the
  browser then rescaled every frame by a fraction of a pixel, which softened
  one-pixel lines and dots.
- **Python sessions**: reading a session now requires each user-defined
  field's `normalizedDims` to name exactly the dimensions in its
  `centroidsByDim`, as the web app writes them.

- A refused upload or event now reaches the client as the status the server sent.
  Every early refusal of a POST -- an expired or unknown session-bundle request, bad
  viewer credentials, a wrong `Content-Type`, an oversized body, an unknown path, a
  foreign `Host` -- was answered before the request body was read, and the socket
  was then closed with that body unread. A socket closed with unread input is reset
  rather than shut down, so a client still sending or not yet reading lost the
  response: Windows reported `WinError 10053`, macOS and Linux a reset or a broken
  pipe. The server now drains a refused request's declared body, within the
  largest accepted upload and for at most five seconds, before it closes. Both
  servers share one `do_POST` that does this.

## [0.9.1] - 2026-10-17

Version 0.9.1 is the PyPI submission release.

### Added

- Scanpy's own `X_umap` is served as written: when no `X_umap_1d`/`X_umap_2d`/`X_umap_3d`
  key exists, a bare `X_umap` of 1, 2, or 3 columns is read at the dimension its column
  count states, without renaming or writing to the caller's object. A dimensional key
  still decides its own dimension and still takes precedence; any other width is refused
  with its shape named. The same rule governs `add_transition_drift_to_obsm` and the
  browser's h5ad and Zarr readers, so one file opens the same way everywhere.
- A long build reports its phases. The direct-AnnData server runs five reported steps, the
  adapter reports obs classification, embedding resolution and which key resolved it,
  vector-field discovery and the graph read, and the manifest and centroid pass is a step
  of its own instead of silence after a success line. `AnnDataAdapter` takes `quiet`.
- `GET /_cellucid/protocol` answers with the wire capabilities the installed package
  accepts, so a web build can ask instead of assuming.
- A notebook learns when a command it sent was refused, rather than reading `None` as
  success.
- Release-contract validation for package, citation, documentation, and downstream recipe
  metadata; installed wheel and sdist gates on Linux, macOS, and Windows before
  publication; reproducible sdist normalization.
- Antialiasing is a live setting, applied on the next frame instead of at context
  creation, and defaults from the cell count until the checkbox is touched.
- A dataset opens at a point size derived from its cell count; the slider now
  reaches 0.050 and every stored session keeps the size it recorded.
- The volumetric smoke render mode is marked `alpha` in the viewer.

### Changed

- **Web app**: the synthetic benchmark generator's custom point-count ceiling is
  raised from 50,000,000 to 250,000,000. The previous value was a declared product
  limit chosen well below the hardware's; the new one is the hardware's, sitting
  just inside the `MAX_TEXTURE_SIZE²` bound on the renderer's interaction-state
  texture (268,435,456 on a context reporting 16,384). A request that the machine
  cannot allocate now fails as an allocation error rather than being refused in
  advance, because the limit that actually binds is the machine's memory — 16 bytes
  per point, so 4 GB of source arrays at the new ceiling — and the control cannot
  know what that is. The point-count input, both validators and the worker boundary
  all read the one constant, so they cannot drift apart.
- The neighbor graph and the vector-field overlay are read only when asked for:
  `serve_connectivity=True` / `--connectivity` and `serve_vector_fields=True` /
  `--vector-fields` on the direct-AnnData path. Each is otherwise never touched,
  and a dataset serving neither publishes exactly what an object holding neither
  publishes. Both are read and validated in full before the socket binds, because
  the manifests the viewer fetches first declare what they contain — and on a
  large object that is the longest part of startup, for a capability most sessions
  never turn on. A 50-neighbor graph over 18 million cells is roughly 900 million
  stored neighbors. A graph asked for and absent is reported; vector fields asked
  for and undeclared serve none, which is ordinary; `vector_field_default` is
  refused unless the fields it selects among are served. The startup report names
  each state separately. `prepare()` is unchanged in Python and R — both already
  took `connectivities=` and `vector_fields=` as data — and now announces all
  three optional capabilities alike when they are not requested.
- Connectivity validation is measurably cheaper: degrees are counted with `np.bincount`
  rather than the unbuffered `np.add.at` (about twenty times faster on that step), one
  coordinate expansion feeds every check, and a graph already in compressed form is
  compared against its own transpose instead of being rebuilt first. Measured on a
  9M-stored-neighbor graph: 1.2x faster end to end and 36% lower peak memory. A graph
  above 50,000,000 stored neighbors logs what it is about to cost.
- The tested Python range is 3.11 through 3.14, and direct AnnData operation targets
  AnnData 0.12, Zarr 3, and numcodecs 0.16.
- PyPI publication uses GitHub OIDC trusted publishing and accepts only a pushed tag
  matching package version metadata exactly.
- Export format: every scientific payload is named by its integer index or a fixed
  neutral name, so no filename carries an observation key, a gene name, or a field id.
- The notebook viewers, the session applier, the data server, and the export writer are
  packages rather than single files. Every public name is imported from where it was, and
  the servers no longer import numpy at startup.
- Adaptive level of detail is bounded by a point budget rather than a fixed ratio.
  A cloud already inside the budget is never reduced, and one above it stops at the
  budget instead of falling to a 44x reduction on any pull-back: at 18.1M cells that
  is 2.4M points rather than 412,000.
- Served gene columns are cached under a byte budget rather than a count of a
  hundred columns, which on an 18.1M-cell object was 7.3 GB of resident memory.

### Fixed

- `cellucid serve --help` no longer fails with `UnicodeEncodeError` on a Windows console
  using a single-byte code page such as CP1252. argparse wrote help, usage, version and
  error text straight to the stream, bypassing the escaping every other console line
  already gets; the CLI's parsers now route that text through the same writer, so any
  character the console cannot show is printed as a backslash escape instead of
  aborting the command. The `serve` help text itself is now plain ASCII.
- A bind address is no longer printed as a URL. `--host 0.0.0.0` used to print
  `http://0.0.0.0:<port>` and try to open it; the banner now prints the loopback origin,
  which is what a browser on the serving machine and the near end of an SSH tunnel both
  use, plus the machine's own name for clients elsewhere. An IPv6 literal bind is
  bracketed as a URL requires.
- An import failure names the condition it is: a distribution absent from the environment
  (`pip install`), one present but too old to hold the name asked of it
  (`pip install --upgrade`), or a standard-library module the running interpreter cannot
  provide -- reported with that interpreter's version, its executable, and the supported
  range, because no package supplies it.
- `cellucid serve` answers an operator mistake with one actionable line instead of a
  traceback, and keeps the traceback for genuine defects.
- Both writers refuse the same numeric input everywhere it can enter. Float32
  underflow (a nonzero value below `2^-149` becoming `0.0`) is now refused for vector
  fields here, and for quantized fields, embeddings, and `latent_space` in the R
  writer, where normalization hid the caller's magnitude from the write-time check.
- A label the browser cannot draw, including an unpaired surrogate, is refused by name.
- `prepare()` reconciles the export root, so it cannot publish coordinates the viewer is
  never told to fetch, and every manifest is checked as written bytes rather than only as
  an in-memory payload.
- Every published payload is little-endian whatever the host is.
- `prepare()` and the web cache both refuse a directory nobody set aside for them,
  including the working directory, the home directory, and paths reached through a link.
- A constant continuous field -- a gene at one level in every cell, or a column a subset
  flattened -- has an encoding instead of failing the export.
- `debug_connection()` reports the dataset probes for a published sample state.
- Publishing a sample's starting view restates its identity. A view is captured
  against a locally served copy of the generation, and the catalog serves the
  same bytes under its own id, so a capture published as-is was refused on every
  open and the sample showed none of its view. The publish step now renames the
  capture -- only after re-deriving the generation's counts and cell-order digest
  from what it ships -- and repoints the catalog digest itself, so the artifact
  and its advertisement cannot part company.
- A refusal that names two identical cell and gene counts as proof of a
  different dataset now names the two dataset identities instead.
- `prepare()` no longer rejects an export because of a gene it was never asked to export.
- A payload that is not entirely finite is refused with counts -- how many NaN,
  infinite, or outside float32, and the first positions -- instead of a bare `500`.
  The server answers `422` with that diagnosis as JSON; the exporter and the R writer
  raise the same report, and the viewer shows it.
- An official sample's published view applies over plain HTTP. Its integrity check used
  `crypto.subtle`, which is secure-context only, so a LAN bind opened the dataset
  uncoloured; the digest is computed the same way with or without WebCrypto.
- `Auto` level of detail follows the camera at every dataset size. Its floor was an
  absolute point budget, unsatisfiable below two million points, so on most datasets
  `Auto` answered full detail at every distance and did nothing; the floor is now the
  smaller of that budget and an 8x reduction cap.
- Level of detail no longer draws square patches. Sorting on a bit-reversed Morton
  code made each level a set of whole lattice cells; sorting on the plain code and
  bit-reversing the rank takes the same fraction from every neighbourhood. At 44x
  reduction: density error 138% to 5%, empty populated regions 42% to none.
- The velocity overlay's `Sync with LOD` read the level scale backwards, running at
  full particle count on the coarsest level and disposing itself at full detail.

### Security

- Jupyter events, session uploads, web asset inventories, and server artifacts use
  authenticated schemas, bounded payloads, and path confinement.
- Both servers validate the `Host` header ahead of routing, file reads, and event
  delivery. A non-loopback bind is an explicit request for network exposure and has no
  authentication of any kind: every user who can reach the port reads the dataset.

### Documentation

- The Python, R, and web guides are organized around the current executable contracts,
  with real browser screenshots for primary workflows.
- Serving from an HPC cluster is documented end to end, including why a compute node
  commonly refuses SSH, the login-node port forward that works instead, why that requires
  `--host 0.0.0.0`, and why the browser address stays `127.0.0.1`.
- Corrected embedding guidance: an all-identical embedding has no finite normalization
  range and is rejected.
- Added a complete standard scVelo recipe for vector fields.

## [0.0.9] - 2026-01-01

This release graduates Cellucid Python out of alpha (still pre-1.0).
The GitHub release was tagged `v0.9.0`, while its Python distribution metadata
declared `0.0.9`. Release 0.9.1 restores one exact version across both systems.

### Added

- AnnData-first viewing and serving with `.h5ad` and `.zarr` support (lazy/backed loading where possible).
- Unified CLI serving with auto-detection for AnnData files, Zarr stores, and
  pre-exported dataset directories.
- Jupyter notebook integration (`CellucidViewer`, `AnnDataViewer`, `show()`, `show_anndata()`) with event hooks and session export.
- Session bundle support (`.cellucid-session`) including `CellucidSessionBundle` and `apply_cellucid_session_to_anndata()` for round-tripping highlights and user-defined fields back into AnnData.
- Multi-dimensional embedding exports (1D/2D/3D) and vector-field overlays (RNA velocity / drift) via `prepare()` and `vector_fields` helpers.
- Hosted web UI proxy mode with on-disk caching helpers (`get_web_cache_dir()`, `clear_web_cache()`).

### Changed

- Export format now includes explicit dataset identity metadata (`dataset_identity.json`) for reproducible sharing and session compatibility checks.
- Reduced export size and improved load performance with optimized manifests, improved connectivity edge export, and optional quantization + gzip compression knobs.

### Security

- Session bundles are treated as untrusted input with bounds checks and dataset mismatch policies when applying to AnnData.

### Documentation

- Major Read the Docs expansion and restructuring (Python package + web app guides), plus new publishing and contributing documentation.

## [0.0.1a0] - 2025

### Added

- Initial alpha release
- AnnData visualization with UMAP embeddings (1D, 2D, 3D)
- Gene expression overlays with sparse matrix support
- Cell metadata coloring (categorical and continuous)
- Interactive filtering and cell selection
- KNN connectivity visualization
- Multiple deployment modes:
  - Local demo (browser-only)
  - Browser file picker (h5ad and exported formats)
  - Server CLI (`cellucid serve`)
  - Python API (`serve()`, `serve_anndata()`)
  - Jupyter integration (`show()`, `show_anndata()`)
- Export functionality for web deployment

[0.9.2]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.2
[0.9.1]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.1
[0.0.9]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.0
[0.0.1a0]: https://github.com/theislab/cellucid-python/releases/tag/v0.0.1a0
