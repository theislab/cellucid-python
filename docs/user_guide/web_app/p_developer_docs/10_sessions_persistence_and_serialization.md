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
`datasetFingerprint` has exactly `sourceType`, `datasetId`, `cellCount`, and
`varCount`. Each chunk has exactly:

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

All four fingerprint values must match. `datasetDependent: false` remains
useful profile metadata for the dockable-layout owner, but it does not authorize
layout-only salvage from a mismatched bundle. The entire operation rejects and
rolls back.

The fingerprint is not a content hash. Dataset publishers must assign a new
identity when cell/variable order or scientific content changes.

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
