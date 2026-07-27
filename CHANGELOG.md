# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.9.1]

The immutable Git tag and PyPI files, when present, are the publication record
for this version; this changelog does not claim a prepublication date.

### Added

- Exact release-contract validation for package, citation, documentation, and
  downstream recipe metadata.
- Installed wheel and source-distribution gates on Linux, macOS, and Windows
  before PyPI publication.
- Reproducible source-distribution normalization with an exact downstream
  recipe SHA-256 gate.
- Reproducible wheel ZIP timestamps through one canonical
  `SOURCE_DATE_EPOCH`.
- An explicit, strictly validated `prepare(created_at=...)` provenance timestamp
  for byte-identical complete export-directory builds.
- Comprehensive current-contract tests for data identity, embeddings, field
  encoding, weighted connectivity, vector fields, sessions, server lifecycle,
  Jupyter messaging, and hosted viewer assets.

### Changed

- The tested Python range is now 3.11 through 3.14.
- Direct AnnData operation now uses one bounded current runtime:
  AnnData 0.12, Zarr 3, and numcodecs 0.16.
- PyPI publication now uses GitHub OIDC trusted publishing and accepts only a
  pushed tag that exactly matches package version metadata and whose commit is
  contained in `origin/main`.
- Server, Jupyter, cache, session, and prepared-export contracts now reject
  invalid or ambiguous inputs before publishing partial state.

### Fixed

- Python 3.14 installation on macOS no longer resolves to the old
  `numcodecs<0.16` build path that failed on Apple silicon.
- Dataset preparation and direct AnnData serving now enforce finite,
  float32-representable scientific values and exact cell/gene alignment.
- Prepared generation, replacement, and server shutdown are atomic and retain
  actionable errors at their owning boundary.
- Cross-platform file names, dataset identifiers, request ranges, and cache
  paths are validated before filesystem or network mutation.
- Prepared exports and direct AnnData responses now use canonical gzip headers,
  eliminating clock and output-path bytes from compressed payloads.

### Security

- Jupyter events, session uploads, web asset inventories, and server artifacts
  use exact authenticated schemas, bounded payloads, and path confinement.

### Documentation

- Reworked the Python, R, and web guides around the current executable
  contracts and added real browser screenshots for primary workflows.
- Corrected embedding guidance: an all-identical embedding is rejected because
  it has no finite normalization range.
- Added a complete standard scVelo Pancreas sample guide covering exact
  catalog selection, 1D/2D/3D navigation, velocity, scientific ownership,
  network payload, provenance, and reproducibility.

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

[Unreleased]: https://github.com/theislab/cellucid-python/compare/v0.9.1...HEAD
[0.9.1]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.1
[0.0.9]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.0
[0.0.1a0]: https://github.com/theislab/cellucid-python/releases/tag/v0.0.1a0
