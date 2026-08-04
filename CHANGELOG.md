# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.9.1]
<!-- CELLUCID_VERSION -->

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

[0.9.1]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.1
[0.0.9]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.0
[0.0.1a0]: https://github.com/theislab/cellucid-python/releases/tag/v0.0.1a0
