# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.9.2]
<!-- CELLUCID_VERSION -->

Version 0.9.2 is a patch release on 0.9.1.

### Added

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

### Fixed

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
