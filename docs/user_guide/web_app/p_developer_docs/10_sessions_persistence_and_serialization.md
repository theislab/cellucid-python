# Sessions: persistence and serialization

This page documents the one current `.cellucid-session` writer/reader contract,
transaction boundary, and extension discipline.

**Audience:** contributors changing persistent state, data-source owners, and
maintainers serving official starting states

**Time:** 35–60 minutes
**Prerequisites:** {doc}`05_app_architecture_overview` and
{doc}`06_state_datastate_and_events`

## Design invariants

- A session stores declarative application state, not the prepared dataset.
- The generic reader accepts exactly one current schema and restores the
  complete document all-or-nothing.
- Manifest, chunk profiles, framing, byte lengths, decoded payloads, and the
  complete dataset fingerprint are exact contracts.
- Row indices are only restored against a dataset whose cell ordering is proven
  to be the one they were saved on; an unprovable ordering rejects the restore.
- Eager and lazy are internal ordering classes inside one awaited public
  operation.
- All contributor state commits or registered state rolls back.
- Explicit empty state is serialized; omission never means “keep the
  destination value.”
- Cancellation and replacement remain observable to the caller but do not
  produce false success or red product diagnostics.

## Code map

- Factory and fixed contributor order:
  `cellucid/assets/js/app/session/index.js`
- Orchestrator, profiles, limits, cancellation, and transaction:
  `cellucid/assets/js/app/session/session-serializer.js`
- Transaction/context and exact fingerprint:
  `cellucid/assets/js/app/session/session-context.js`
- Container framing:
  `cellucid/assets/js/app/session/bundle/format.js`
- gzip/DEFLATE codec:
  `cellucid/assets/js/app/session/codecs/gzip.js`
- Feature owners:
  `cellucid/assets/js/app/session/contributors/`
- State capture helpers:
  `cellucid/assets/js/app/state-serializer/`
- UI command boundary:
  `cellucid/assets/js/app/ui/modules/session-controls.js`

## Container and manifest

The container is:

1. `CELLUCID_SESSION\n`;
2. unsigned 32-bit little-endian manifest length;
3. exactly that many UTF-8 JSON bytes; and
4. for every manifest chunk, an unsigned 32-bit stored length plus exact stored
   payload bytes.

The manifest has exactly `createdAt`, `datasetFingerprint`, and `chunks`.
`datasetFingerprint` has exactly `sourceType`, `datasetId`, `cellCount`,
`varCount`, and `cellOrder`. `cellOrder` has exactly `dimension` (the integer
1, 2, or 3) and `digest` (16 lowercase hexadecimal characters). Each chunk has
exactly:

- `id`
- `contributorId`
- `priority`
- `kind`
- `codec`
- `label`
- `datasetDependent`
- `storedBytes`
- `uncompressedBytes`

There is deliberately no `version`, `dependsOn`, compatibility marker, or
extension bag.

## Closed current inventory

Generic Save State always emits the current singleton profiles:

| ID | Owner | Priority | Kind/codec |
|---|---|---|---|
| `core/field-overlays` | `field-overlays` | eager | JSON/gzip |
| `core/state` | `core-state` | eager | JSON/gzip |
| `ui/dockable-layout` | `dockable-layout` | eager | JSON/gzip |
| `analysis/windows` | `analysis-windows` | eager | JSON/gzip |
| `highlights/meta` | `highlights-meta` | eager | JSON/gzip |
| `analysis/cache-inventory` | `analysis-artifacts` | eager | JSON/gzip |
| `cinematic/camera` | `cinematic-camera` | eager | JSON/gzip |

`analysis/cache-inventory` and `cinematic/camera` are mandatory even when
empty. Their emptiness means replace the destination cache/path with empty
state.

Dynamic chunks are also closed:

- `user-defined/codes/<fieldId>`: one binary/gzip chunk per categorical
  user-defined field; eager exactly for target live/snapshot active fields,
  lazy otherwise;
- `highlights/cells/<groupId>`: one lazy binary/gzip membership chunk per
  highlight group; and
- `analysis/artifacts/bulk-gene/<cacheKey>`: lazy binary/gzip chunks exactly
  matching the ordered cache inventory.

Static profile aliases, missing unused codes, extra artifacts, priority lies,
label changes, or reordered contributor groups are invalid.

## Capture ordering

Contributors are registered in `session/index.js`. Capture preserves their
relative order within separate eager and lazy buckets. This lets critical code
columns appear before `core/state`, while inactive columns can remain in the
lazy bucket.

Before writing:

1. every contributor chunk descriptor is exact-validated;
2. JSON or binary payload is encoded;
3. uncompressed and stored limits are enforced;
4. gzip is applied where declared;
5. exact byte counts are recorded; and
6. the resulting manifest is revalidated against the current generic profile.

## Restore pipeline

The public file and URL entry points use one owner-controlled restore:

1. a newer restore supersedes and awaits the older owner;
2. framing streams under stored-byte and aggregate limits;
3. the exact manifest/profile/order and dataset fingerprint validate;
4. small reversible owner snapshots are registered;
5. eager chunks decode/apply in order;
6. lazy chunks decode/apply in the same awaited operation, with macrotask
   yielding where needed;
7. every transaction participant prepares;
8. every participant commits in order;
9. the final UI refresh commits; and
10. the serializer publishes **Session fully restored**.

Any error invokes rollback in reverse participant order and awaits asynchronous
rollback. If rollback itself fails, that failure is preserved as the original
error's cause rather than replacing the primary diagnosis.

## Reversible large-state ownership

Avoid cloning cell-scale arrays:

- field overlays retain and reattach exact prior categorical typed-array
  references;
- the analysis data layer replaces Map/LRU container ownership and can restore
  the exact prior containers;
- in-flight cache writers must be generation-isolated so they cannot cross the
  replacement boundary; and
- cinematic restore snapshots keyframes, navigation/settings, and actual
  stopped/playing/paused runtime state.

The cache inventory verifies exact count and order during transaction prepare.
Cinematic autoplay may begin only after successful commit. Rollback restores
the previous playback state rather than merely copying the autoplay checkbox.

## Bounded gzip

Declared `storedBytes` and `uncompressedBytes` are necessary but not sufficient:
native `DecompressionStream` may materialize one very large output chunk,
especially in WebKit. The codec therefore preflights the complete
single-member DEFLATE structure before constructing the native decompressor.

The preflight:

- accepts stored, fixed-Huffman, and dynamic-Huffman blocks;
- validates trees, distances, back-references, trailer, ISIZE, and exact output
  length;
- rejects concatenated members and trailing data;
- checks cancellation throughout; and
- yields by bounded macrotask intervals so timers and user input can abort large
  valid streams.

Only an exact-length valid member reaches native decompression.

## Dataset identity

`getDatasetFingerprint()` in `session-context.js` builds the record that is
written into the manifest at save and re-derived from live state at restore:

| Field | Derived from | Notes |
|---|---|---|
| `sourceType` | `dataSourceManager.getCurrentSourceType()` | `null` only when no data-source manager owns the loaded state |
| `datasetId` | `dataSourceManager.getCurrentDatasetId()` | same nullability |
| `cellCount` | `state.pointCount` | |
| `varCount` | `state.varData.fields.length` | |
| `cellOrder.dimension` | `state.getViewDimensionLevel('live')` | exactly 1, 2, or 3 |
| `cellOrder.digest` | `digestCellOrder(state.positionsArray)` | 16 lowercase hex characters |

`getDatasetFingerprint()` throws unless `state.positionsArray` is a
`Float32Array` of exactly `cellCount * 3` values, so a fingerprint is never
produced from a partially published dataset.

### Why `cellOrder` exists

Selections, categorical codes, and cached analysis values persist as raw
dataset row indices. The four scalars all survive a row permutation: re-export
a dataset at the same id from re-sorted input and `sourceType`, `datasetId`,
`cellCount`, and `varCount` are unchanged, so every restored index would denote
a different cell with nothing to detect it.

Observation and variable names were the alternative and do not work here — they
are invariant under a row permutation, which is the change being caught.
Positions are the only per-cell payload the viewer always holds in memory, and
a re-ordered export permutes them.

`digestCellOrder()` reads the coordinate buffer as 32-bit words and folds two
32-bit accumulators over it, each updated with a multiply and a rotation per
word, then emits them as one 64-bit hex value. Per-word rotation is what makes
the result depend on the order of the words rather than only on their multiset.
On 842k cells (9.64 MiB of Float32 coordinates) the digest costs 4.3 ms and the
complete fingerprint 5.6 ms. The result is memoized in a `WeakMap` keyed by the
coordinate array itself; that array is replaced, never rewritten, when the
dataset or the displayed embedding changes, so the per-contributor re-derivations
during a capture cost 0.06 ms in total.

`cellOrder.dimension` is part of the record because the 1D, 2D, and 3D
embeddings are separate exported files, normalized independently. A digest is
therefore only comparable within one of them, and the dimension is what lets a
mismatch be attributed to the view on screen instead of to the data.

### Validation and comparison

`assertDatasetFingerprint(value, context)` validates one record from untrusted
input: exact five-key set, nullable strings for the source identity, safe
integers for the counts, and `cellOrder` with exactly `dimension` in 1–3 and a
`/^[0-9a-f]{16}$/` digest. It runs one check before the schema check: a record
whose key set is exactly the four pre-`cellOrder` keys is a session file written
before the record existed, and is refused with
`SESSION_WITHOUT_CELL_IDENTITY_MESSAGE` rather than reported as a schema
violation. `session-serializer.js` calls it from `validateManifest()`, so this
refusal happens before any contributor is consulted.

`datasetFingerprintMatches(a, b)` validates both records and then compares all
six scalars. `describeDatasetFingerprintMismatch(saved, current)` returns `null`
when they match and otherwise the sentence the session controls surface
verbatim. `_restoreFromBundleSource()` throws that string as a `RangeError`.

Four refusals are distinguished, and attributing one to the wrong cause would be
its own integrity failure — telling a user their data was re-ordered when they
merely switched the view would make them distrust a sound dataset:

| Cause | Detected by | Message |
|---|---|---|
| A session file written before `cellOrder` existed | `assertDatasetFingerprint` during manifest validation | says the file predates cell-content recording, so its selections can never be confirmed, and asks for the selections to be re-created and saved again |
| A different dataset | any of the four scalars differs | names the saved and current cell and gene counts and asks for the dataset the session was saved on |
| A different embedding on screen | `cellOrder.dimension` differs | names both dimensions, explains that the check reads the coordinates on screen, and asks for a switch back to the saved dimension |
| A re-ordered dataset | `cellOrder.digest` differs, everything else matches | says the dataset has the same name and the same cell and gene counts but stores its cells in a different order, so every saved selection would mark the wrong cells |

None of these messages carry identifiers or internal vocabulary.

A file predating the record is refused rather than accepted. Accepting it would
mean keeping exactly the unverifiable state the record exists to eliminate,
permanently, for every file already written.

`datasetDependent: false` remains useful profile metadata for the
dockable-layout owner, but it does not authorize layout-only salvage from a
mismatched bundle. The entire operation rejects and rolls back. This applies to
`restorePublishedDefaultState()` too: the published-default profile shares
`validateManifest()` and the same fingerprint comparison, so an advertised
`default.cellucid-session` must carry `cellOrder` and must have been saved on
the dimension the dataset publishes as its default.

The fingerprint is still not a content hash: the digest covers coordinates, not
expression or metadata. Dataset publishers must assign a new identity when
scientific content changes. What `cellOrder` adds is that a re-ordered
republication under an unchanged identity is now detected instead of trusted.

## Official sample profile

Catalog-advertised `default.cellucid-session` uses a separate exact profile:

1. `core/field-overlays`
2. `core/state`
3. `ui/dockable-layout`
4. `analysis/windows`
5. `highlights/meta`

The catalog advertises the state manifest and SHA-256 as a pair. Transport is
strictly smaller and bounded. The loader verifies the manifest and digest and
applies this profile only through `restorePublishedDefaultState()` after the
scientific dataset is published and while the dataset selection generation is
still current.

Do not route this five-chunk artifact through generic Load State. A browser test
for the advertised path must author a real five-chunk fixture; reusing generic
Save State output is invalid by construction.

## Cancellation outcomes

Two exact coded aborts exist:

- user-canceled restore; and
- restore superseded by a newer owner.

Direct serializer APIs reject with the coded `AbortError`, allowing lifecycle
owners to settle correctly. Progress owners dismiss only those exact outcomes.
The Session panel maps them to a no-op and does not log or post failure/success.
An unrelated plain `AbortError` remains a real failure.

## Adding or changing persistence

Persistence changes are schema changes to the only current format:

1. identify the sole feature owner and its rollback boundary;
2. choose eager only when another target chunk needs the state before commit;
3. define exact static or dynamic ID, metadata, payload keys, and size bounds;
4. encode empty replacement state explicitly;
5. add manifest-wide completeness/order checks where a dynamic family depends
   on another chunk;
6. retain cell-scale buffers by ownership/reference where safe; never clone
   them merely for rollback;
7. prove success, late failure, commit failure, async rollback, cancel,
   supersession, empty replacement, metadata mutation, and large-array
   behavior;
8. update this page, the user-facing session chapter, the maintainer README,
   and exact documentation tests in the same change.

Do not introduce a reader for an older shape. Update every current producer,
fixture, consumer, and published artifact together.

## Intentional exclusions

Generic sessions do not contain:

- prepared points, observation tables, or full expression matrices;
- dataset picker/URL/GitHub/Jupyter connection state;
- Community Annotation authentication, votes, or repository workflow;
- Figure Export controls;
- Benchmarking controls;
- notifications, hover, focus, pending selections, or in-flight requests; or
- DOM, worker, network, and WebGL runtime objects.

Next: {doc}`11_analysis_architecture`.
