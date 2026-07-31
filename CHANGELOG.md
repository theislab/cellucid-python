# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.9.1]

Version 0.9.1 is the PyPI submission release.

### Fixed — export format: a constant continuous field has an encoding

A gene expressed at one identical level in every exported cell — very often
zero, once an atlas is subset to one lineage — and a continuous `obs` column a
subset flattened are both ordinary scientific data. `prepare()` used to refuse
the whole export when it met one, which cost real features: subsetting Suo to
`LVL0 == "Haematopoeitic_lineage"` leaves four genes detected in none of the
published cells, one of them the nameable `NFIB-AS1`.

compact_v1 now has a named case for them:

- `minValue == maxValue` and every code `0`. The entry keeps its exact shape and
  length — `[index, key, minValue, maxValue]` — so nothing about the manifest
  contract moves.
- The writer takes an explicit branch and never derives a scale, so nothing
  divides by `maxValue - minValue`.
- The reader takes the matching branch and returns `minValue` directly, so the
  constant comes back bit-exact rather than as an approximation of itself.
- Native-double variation finer than float32 resolution collapses to one float32
  value, which is one constant, and is published the same way instead of being
  rejected.

`cellucid_prepare()` in the R package implements the identical case. The only
continuous payload still without an encoding is a categorical field's generated
outlier quantiles when *every* quantile is missing, because no category holds
`centroid_min_points` cells — there is no value to publish at all, and that
failure still names every affected field in one run before any file is written.

### Changed — export format: payload files are named by index

Every scientific payload an export writes is now named by its integer index, or
by a fixed neutral name. No filename anywhere in the tree carries an observation
key, a gene name, or a vector field id:

```text
var/0.values.u8.gz     obs/0.codes.u8.gz     obs/0.outliers.u8.gz
obs/1.values.u8.gz     vectors/0_2d.bin.gz
```

- `_varSchema.pathPattern` and every `_obsSchemas` path pattern now substitute
  `{index}` instead of `{key}`, and each field entry declares its own index as
  element `[0]`:
  - `var` — `[index, name]`, or `[index, name, minValue, maxValue]` quantized;
  - obs continuous — `[index, key]`, or `[index, key, minValue, maxValue]`;
  - obs categorical — `[index, key, categories, dtype, sentinel, centroids]`,
    plus `outlierMinValue`/`outlierMaxValue` when quantized.
  A `var` entry and an obs-continuous entry now share a shape deliberately, so a
  reader must take a field's kind from the array it came from, never from the
  entry's length.
- Within one axis directory the indices are exactly `0 … N-1`, each used once,
  and `obs` shares that one space across `_continuousFields` and
  `_categoricalFields` because both write into `obs/`. `prepare()` asserts this
  against the manifest it has just built, and re-derives every declared payload
  path to compare it with the directory it actually wrote — two fields holding
  one index would write one payload over another and the app would then draw one
  field's values under another field's name, with nothing raised.
- The Jupyter/direct AnnData server serves the same index routes
  (`/var/0.values.f32`, `/obs/0.codes.u8`, `/vectors/0_2d.bin`) and emits
  byte-identical manifest shapes to the batch exporter. Only the exact unpadded
  decimal form resolves; `/var/00.values.f32` is a 404.
- Because an identifier is no longer a path, the filename-portability,
  case-collision, and Windows-reserved-name rules are **removed** from
  observation keys, gene names, and vector field ids. `HLA-DRB1/2`, `CON`,
  `% mito`, `细胞`, and the two distinct fields `Field` and `field` now export
  unchanged. What remains is what an identity is for: non-empty, distinct within
  its axis, and drawable exactly as stored — the same display-text rule every
  string category label already obeyed. `dataset_id` is unaffected: it names a
  real directory and a served URL, so it keeps the full portable-component rule.
- The direct AnnData adapter no longer carries its own copy of that identifier
  validator. It had sanitized unsafe characters into `_` where the batch
  exporter rejected them, so the two producers disagreed about which inputs were
  admissible; both now call the one shared rule.

An index is a position, not a stable name: adding a gene or reordering
`obs_keys` renumbers the payload files. Resolve a payload through the manifest
every time rather than caching or hard-coding a path across exports. Exports are
regenerated from source, so re-export with `force=True`; there is no dual-read
path.

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
  to export. The two gene-identifier rules that survive the move to
  index-named payloads have two different scopes, and each is now checked at
  its own. Being drawable is a property of a name the viewer shows, so it
  covers the genes `gene_identifiers` selects; a `var` row left out reaches no
  manifest and is not checked, exactly as an `obs` column left out of
  `obs_keys` is not. Distinctness spans the whole `var`, because
  `gene_identifiers` addresses `var` rows by identifier and a repeated one
  names no single row. `cellucid_prepare()` in the R package documents the
  identical pair of scopes.
- Python 3.14 installation on macOS no longer resolves to the old
  `numcodecs<0.16` build path that failed on Apple silicon.
- Dataset preparation and direct AnnData serving now enforce finite,
  float32-representable scientific values and exact cell/gene alignment.
- Prepared generation, replacement, and server shutdown are atomic and retain
  actionable errors at their owning boundary.
- Prepared-directory writers share one persistent exact-target lock with the R
  exporter, reject concurrent publication, and recover after process death
  without leaking descriptors.
- Dataset identifiers, request ranges, and cache paths are validated before
  filesystem or network mutation.
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
