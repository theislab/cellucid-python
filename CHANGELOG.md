# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.9.1]

Version 0.9.1 is the PyPI submission release.

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
- Session-bundle readers now accept a `cellOrder` record in the manifest's
  `datasetFingerprint`, beside the existing `sourceType`, `datasetId`,
  `cellCount`, and `varCount`. The web app records it so that a dataset
  republished under the same identity with the same cell and gene counts but its
  rows in a different order can no longer silently re-point every saved
  selection at different cells. The field is **accepted, not required**, and its
  digest is never re-derived here: it is taken over exported coordinates, while
  `apply_cellucid_session_to_anndata()` targets an `AnnData` whose row order the
  caller controls, so enforcing that identity belongs to the viewer. Accepting
  both fingerprint shapes is also what keeps existing session files loadable. A
  `cellOrder` that is present is shape-checked strictly: exactly `dimension` and
  `digest`, `dimension` exactly 1, 2, or 3, and `digest` exactly 16 lowercase
  hexadecimal characters.
- An explicit `allowed_hosts=` declaration on `CellucidServer`, `AnnDataServer`,
  `serve`, `serve_anndata`, `CellucidViewer`, `AnnDataViewer`, `show`, and
  `show_anndata`, plus a repeatable `cellucid serve --allowed-host HOST`, naming
  the exact extra `Host` authorities a server answers to. This is what a reverse
  proxy such as `jupyter-server-proxy` requires, because the proxy forwards the
  browser's `Host` verbatim and a proxied request carries no signal that
  separates it from a rebound one. Entries are bare host names — no port, no
  scheme, no wildcard — validated when the server is constructed and matched on
  any port; declaring names also applies the `Host` check to a non-loopback
  bind.

### Changed

- The tested Python range is now 3.11 through 3.14.
- Direct AnnData operation now uses one bounded current runtime:
  AnnData 0.12, Zarr 3, and numcodecs 0.16.
- PyPI publication now uses GitHub OIDC trusted publishing and accepts only a
  pushed tag that exactly matches package version metadata and whose commit is
  contained in `origin/main`.
- Server, Jupyter, cache, session, and prepared-export contracts now reject
  invalid or ambiguous inputs before publishing partial state.
- Every string Cellucid prints verbatim must now read on screen as the value it
  stores. A string category label, `dataset_name`, `dataset_description`,
  `source_name`, `source_url`, and `source_citation` are rejected when they
  carry a control character, one of the zero-width characters `U+200B`,
  `U+2060`, or `U+FEFF`, or leading or trailing whitespace of any kind
  including `U+00A0` NO-BREAK SPACE. An empty category label is rejected, as
  are two labels in one field that a whitespace-collapsing renderer draws
  identically. Nothing is trimmed: trimming would rewrite an annotation the
  caller never asked to change, and would merge `"Liver "` into a separate
  `"Liver"` category and move cells between them. The message names every
  offending label in the field at once and gives the one-line repair. The
  `cellucid` R package enforces the identical rule.
- `dataset_description`, `source_name`, `source_url`, and `source_citation`
  previously accepted padded text in Python while the R writer rejected it;
  the two writers now agree.
- `serve_anndata()` and `show_anndata()` now apply the same category-label rule
  as `prepare()`; the adapter's duplicate label validator was removed in favour
  of the exporter's.
- `generate_datasets_manifest()` now validates a published
  `dataset_identity.json` `description` the same way it validates its `name`.

### Fixed

- `prepare()` no longer rejects an export because of a gene it was never asked
  to export. The portable-filename-component rule and the case-insensitive
  collision rule are about the file a gene is written to, so they now apply to
  the genes `gene_identifiers` selects rather than to every row of `var`. A
  `var` carrying an identifier such as `HLA-DRB1/2`, a Windows device name, or
  a case-only twin of another identifier now exports cleanly whenever
  `gene_identifiers` leaves that gene out, which is what `obs_keys` has always
  done for observation keys and what `cellucid_prepare()` in the R package has
  always done for genes. The rule that is *not* about paths is unchanged and
  still spans the whole `var`: every gene identifier must be a non-empty string
  and must be distinct, because `gene_identifiers` addresses `var` rows by
  identifier and a repeated one resolves to no single row. Identifiers that are
  exported are validated exactly as before.
- Python 3.14 installation on macOS no longer resolves to the old
  `numcodecs<0.16` build path that failed on Apple silicon.
- Dataset preparation and direct AnnData serving now enforce finite,
  float32-representable scientific values and exact cell/gene alignment.
- Prepared generation, replacement, and server shutdown are atomic and retain
  actionable errors at their owning boundary.
- Prepared-directory writers share one persistent exact-target lock with the R
  exporter, reject concurrent publication, and recover after process death
  without leaking descriptors.
- Cross-platform file names, dataset identifiers, request ranges, and cache
  paths are validated before filesystem or network mutation.
- Prepared exports and direct AnnData responses now use canonical gzip headers,
  eliminating clock and output-path bytes from compressed payloads.

### Security

- Jupyter events, session uploads, web asset inventories, and server artifacts
  use exact authenticated schemas, bounded payloads, and path confinement.
- The data servers validate the `Host` header ahead of routing, file reads, and
  event delivery. A loopback bind stops remote clients but not DNS rebinding:
  once an attacker page's name re-resolves to `127.0.0.1` the browser treats the
  server as same-origin, sends no `Origin`, and runs no CORS check. Only one
  well-formed `Host` naming `localhost` or a loopback IP literal, on the port
  actually bound, is served; every other authority receives `421 Misdirected
  Request` with a body that never echoes the supplied value.

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

[0.9.1]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.1
[0.0.9]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.0
[0.0.1a0]: https://github.com/theislab/cellucid-python/releases/tag/v0.0.1a0
