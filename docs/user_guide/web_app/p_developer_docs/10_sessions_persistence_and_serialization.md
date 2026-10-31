# Sessions: persistence and serialization

This page documents the one current `.cellucid-session` writer/reader contract,
transaction boundary, and extension discipline.

**Audience:** contributors changing persistent state, data-source owners, and
maintainers serving official starting states

**Time:** 35–60 minutes
**Prerequisites:** {doc}`05_app_architecture_overview` and
{doc}`06_state_datastate_and_events`

## Design invariants

- A session stores declarative application state, not the dataset.
- The generic reader accepts exactly one session format — the version named
  in the file's signature, currently 3 — and restores the complete document
  all-or-nothing. A file of an earlier format is refused by its signature,
  before anything else is read; there is no migration.
- A session stores choices, not what the dataset already holds: labels, filter
  summaries, field and category positions, and per-cell gene values are
  recomputed or re-read on restore.
- Manifest, chunk profiles, framing, byte lengths, decoded payloads, and the
  complete dataset fingerprint are exact contracts.
- Every per-cell chunk (highlight memberships, user-defined codes and values)
  names cells by **original position**, ascending: cell `i` is the open cell
  with the `i`-th smallest original position (every cell open: position `i`;
  a subsample: `kept[i]`), whatever order the file stores the cells in (a file
  in priority order) or the app opened them in (a subsample's priority order).
  `assets/js/app/session/session-cell-order.js` translates between that order
  and the open order on capture and restore, so Python and R apply a session to
  an object in its original order without knowing any layout.
- Those positions are only restored against a dataset whose cells are proven
  to be the ones they were saved on; an unprovable ordering rejects the restore.
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
- Container framing and format version:
  `cellucid/assets/js/app/session/bundle/format.js`
- Codecs: `cellucid/assets/js/app/session/codecs/` — `gzip.js` (gzip/DEFLATE),
  `varint.js`, `cell-set.js` (highlight memberships), `typed-json.js`
  (analysis settings and results)
- Feature owners:
  `cellucid/assets/js/app/session/contributors/` (among them
  `user-defined-values.js`, the values of a column a notebook pushed, and
  `gene-sets.js`, the gene sets a notebook pushed)
- The cells a session was saved on (the manifest's `cells`) and the order of
  its per-cell chunks: `cellucid/assets/js/app/session/cell-subset-record.js`,
  `cellucid/assets/js/app/session/session-cell-order.js`
- State capture helpers:
  `cellucid/assets/js/app/state-serializer/`
- UI command boundary:
  `cellucid/assets/js/app/ui/modules/session-controls.js`

## Container and manifest

The container is:

1. the 19-byte signature `CELLUCID_SESSION/3\n` (`SESSION_FORMAT_VERSION` in
   `bundle/format.js`; a file with another format version's signature is
   refused naming it, `sessionSignatureRefusal()`: a format-2 file, saved by an
   earlier Cellucid, carries a cell-order digest no later version can check;
   any other file is refused with `NOT_A_SESSION_FILE_MESSAGE`);
2. unsigned 32-bit little-endian manifest length;
3. exactly that many UTF-8 JSON bytes; and
4. for every manifest chunk, an unsigned 32-bit stored length plus exact stored
   payload bytes.

The manifest has exactly `createdAt`, `datasetFingerprint`, `source`,
`requires`, and `chunks`, plus `cells` for a session saved on a subsample
(`cell-subset-record.js`: `{spec, kept, digest}` — the subsample's
specification, the kept original positions as the LEB128 cell-set bytes in
base64, and their digest, {doc}`19_subsampling_specification`; its `k` must
equal the fingerprint's `cellCount`). A manifest without `cells` was saved on
every cell. Loading such a session opens exactly those cells (the dataset is
read again first when it is open with others) before the restore runs. A
session saved while genes were read for a sample of the open cells
(`data/anndata/subsample-plan.js`, `genes`) also has `genes`
(`cell-subset-record.js` `captureGeneSampleRecord`: `{cells, order,
digest}` — how many of the open cells, the first in their order; every
cell's order, `{method: 'blocks', seed, blockSize}`, when every cell was
open, else null; and the cell-set digest of their original positions).
Loading it reads genes for those cells: in the request that opens the
dataset's cells (`genes`, `order`), or, on the cells already open in the
same order, by reading genes for more cells (`readGenesForMoreCells`); with
the cells in another order, the open cells are read again in the session's.
cellucid-python (`CellucidSessionBundle.genes`) and cellucid-r check it as
the web app does. `source` (`session/dataset-source-record.js`) is
where the dataset was opened from — a sample, a data link's locator, a server
address, or a local file's or folder's name (`{kind: 'file', name, dataset}`),
with the dataset id — or `null` (a Python or R viewer); **Load** reads only the manifest
first and, when another dataset (or none) is open and the source can be opened
(not a local file), opens it through the same connection controls a person
would use before the restore runs. `requires` (`collectRequiredFields()` in
`session-context.js`) is `{obs, var}`: the sorted dataset keys of every
observation field and gene a view shows or filters on or a highlight group is
bound to (a renamed field by its original key; user-defined fields excluded).
`datasetFingerprint` has exactly `datasetId`, `name`, `cellCount`, `varCount`,
`cellOrder`, and `dataMapping`; `name` is the dataset's display name or
`null`, used in messages and never compared. `cellOrder` has exactly `embedding`
(the embedding reference of the layout the live view shows), `dimension` (that
layout's dimension, the integer 1, 2, or 3) and `digest` (16 lowercase
hexadecimal characters, the cell-order digest below). `dataMapping` is the dataset's canonical data mapping
(`data/data-mapping.js`), always of `kind` `"anndata"`, with exactly
`embeddings`, `knn`, `geneNames`, `expression`, `log1p`, `highlyVariable`,
`latent` and `vectorFields`, plus `modality` for a dataset of several
modalities (a MuData file; absent otherwise): its embeddings are
`{enabled, defaults}`, the embeddings that are turned on and each dimension's
default, and its vector bindings are `{key, embedding}`. Each chunk has
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

The format version lives in the signature only; there is no `dependsOn`,
compatibility marker, or extension bag.

## Closed current inventory

A generic **Save** always emits the current singleton profiles, in
contributor order:

| ID | Owner | Priority | Kind/codec |
|---|---|---|---|
| `core/field-overlays` | `field-overlays` | eager | JSON/gzip |
| `cinematic/camera` | `cinematic-camera` | eager | JSON/gzip |
| `core/state` | `core-state` | eager | JSON/gzip |
| `ui/dockable-layout` | `dockable-layout` | eager | JSON/gzip |
| `ui/figure-export` | `figure-export` | eager | JSON/gzip |
| `analysis/windows` | `analysis-windows` | eager | binary (typed JSON)/gzip |
| `analysis/docked` | `analysis-docked` | eager | binary (typed JSON)/gzip |
| `highlights/meta` | `highlights-meta` | eager | JSON/gzip |

`cinematic/camera` is mandatory even when empty; its emptiness means replace
the destination path with an empty one. `ui/figure-export` is the record
`exportSessionSettings()` of the figure-export panel returns (versioned within
itself, validated by `figure-export/utils/session-settings.js`).
`analysis/docked` is `comparisonModule.exportSessionState()` — the open mode
and every initialized mode's settings and result — and `analysis/windows` is
`{windows}`, each window's mode, geometry, settings and result.

Dynamic chunks are also closed:

- `user-defined/codes/<fieldId>`: one binary/gzip chunk per categorical
  user-defined field; eager exactly for target live/snapshot active fields,
  lazy otherwise; and
- `highlights/cells/<groupId>`: one lazy binary/gzip membership chunk per
  highlight group, in group order.

Static profile aliases, missing unused codes, priority lies, label changes, or
reordered contributor groups are invalid.

### Payload encodings

- **Highlight cells** (`codecs/cell-set.js`): byte 0 selects sorted-delta
  uvarints (0), a bitset over all cells (1), runs of uvarint gap and length
  (2), a bitset over the cells of a reference category (3) or the whole
  reference category (4); the writer keeps the smallest. 3 and 4 apply to a
  category group on a dataset observation field; they carry the category's
  size and a 32-bit order-independent digest (`digestCellSet()`), and the
  reader recomputes both from the loaded field and refuses a difference.
- **User-defined codes** (`contributors/user-defined-codes.js`): an encoding
  byte, the uvarint cell count, then raw codes at the declared width (0), runs
  of uvarint value and length (1), or codes packed at 1, 2 or 4 bits (2); the
  writer keeps the smallest.
- **Analysis** (`codecs/typed-json.js`): a u32 JSON length, JSON in which the
  k-th typed numeric array is `{"$typed": k}`, then each array as a dtype
  byte, a uvarint length and its little-endian elements. Objects are rebuilt
  with own properties, so a `__proto__` key in a crafted file is data, not a
  prototype.
- **Core state**: a kept view records its id, embedding, camera, active field,
  filters and centroid toggles; the live view its camera and centroid toggles
  beside the top-level embedding, filters and active field. A view's label and
  filter summary are recomputed when it is rebuilt. Highlight groups record
  `fieldKey`/`categoryName`, not positions; so do a categorical filter's hidden
  categories (`visibility: [{categoryName, visible: false}]`), where the
  category named `null` is the cells without a category (the legend's Missing
  row). A name the field no longer has is refused, except on a text or boolean
  column (`discovered`: its categories are the labels its open cells hold):
  there a hidden or recoloured category no open cell holds (a published
  starting view on a subsample) is held, said once ("2 categories of
  "cell_type" in this session aren't among the open cells…"), and applied when
  cells holding it open and the column grows it
  (`state/managers/field/pending-categories.js`). Category renames are by
  label already and apply when the label appears. `#hp-antialias` is never captured
  (a device capability); the theme is captured and restored through its owner
  like any other control.

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
3. the exact manifest/profile/order and dataset fingerprint validate, the
   session's data mapping is applied, and every key of `requires` is checked
   (`describeMissingRequiredFields()` names every missing one);
4. small reversible owner snapshots are registered;
5. eager chunks decode/apply in order;
6. lazy chunks decode/apply in the same awaited operation, yielding a
   macrotask (`utils/yield-macrotask.js`) between chunks — never an animation
   frame, which does not run in a background tab or an off-screen notebook
   iframe;
7. every transaction participant prepares;
8. every participant commits in order;
9. the final UI refresh commits; and
10. the serializer publishes **Session fully restored** by completing the
    restore's own progress notification — the one report of the outcome, as a
    failure ends that same entry as **Loading failed**; a session file you
    chose announces its outcome however quickly it restored. A sample's own
    view (a published default) restores without a word: it is part of
    opening the sample, not a session you asked for. A published sample view
    is the exception on failure: its caller knows the dataset itself loaded and
    reports **Sample view unavailable**, so the progress entry is dismissed
    rather than failed a second time.

Any error invokes rollback in reverse participant order and awaits asynchronous
rollback. If rollback itself fails, that failure is preserved as the original
error's cause rather than replacing the primary diagnosis.

## Reversible large-state ownership

Avoid cloning cell-scale arrays:

- field overlays retain and reattach exact prior categorical typed-array
  references; and
- cinematic restore snapshots keyframes, navigation/settings, and actual
  stopped/playing/paused runtime state.

The core state rebuilds the view graph in one pass
(`state-serializer/multiview.js`): every field, gene and embedding any view
names is requested together and the embeddings are pinned until the views are
built (`preloadSessionViews()` in `contributors/core-state.js`); each kept view
is replayed once on the live view — filters, active field, embedding — and
published; the live view's own state is applied last. Nothing is reset to a
base state between views, because each replay applies a complete state.

Analysis panels, analysis windows and figure-export settings are applied in
the transaction's presentation stage, after every state participant commits.
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
| `datasetId` | `dataSourceManager.getCurrentIdentityId()` | the dataset's identity id, whatever route opened it (a local folder's `local-user:` key is not used); `null` only when no data-source manager owns the loaded state |
| `name` | `dataSourceManager.getCurrentMetadata().name` | trimmed display name or `null`; informational, not compared |
| `cellCount` | `state.pointCount` | |
| `varCount` | the dataset's own var fields in `state.varData.fields` | user-defined fields (for example a duplicated gene) are not counted, so a session that added one still matches the freshly loaded dataset |
| `cellOrder.embedding` | `state.getViewLayoutKey('live')`, named under the data mapping | `{key, columns}`: an `obsm` key with its 0-based columns |
| `cellOrder.dimension` | `state.getViewDimensionLevel('live')` | exactly 1, 2, or 3 |
| `cellOrder.digest` | `state.dimensionManager.getCellOrderDigest(<the live layout>)`: the reader worker's sums of that embedding's stored values (`data/cell-order-digest.js`) | 16 lowercase hex characters |
| `dataMapping` | `state.getDataMapping()` | a canonical data mapping |

`getDatasetFingerprint()` is asynchronous (the digest is the reader worker's)
and throws unless `state.positionsArray` is a `Float32Array` of exactly
`cellCount * cellOrder.dimension` values, so a fingerprint is never produced
from a partially published dataset.

### Cell-order digest

The digest is a function of the cells and their stored coordinates only, so
the same cells of the same file digest alike however they were opened: whole,
grown step by step, raised or lowered by **Apply**, cut in place to a
session's share, in priority order or in the order the file gives them.

- Each cell contributes a 64-bit hash, two 32-bit lanes: MurmurHash3 (x86,
  32-bit) of four little-endian 32-bit words, under seed `0x9e3779b1` (low
  lane) and `0x85ebca6b` (high lane). The words are the cell's **original
  position** ({doc}`19_subsampling_specification`) and the float32 bit
  patterns of its values on axes 0, 1 and 2 of the embedding's `obsm` columns,
  each stored value rounded to float32 (`+0.0`, the zero word, for an axis a
  1D or 2D embedding does not have; every NaN as `0x7fc00000`).
- A set of cells sums its cells' lanes, each modulo 2^32: the order of the
  cells does not matter, and a set's sums are the sums of its parts.
- The digest of `n` cells is `fmix32(high ^ fmix32(n ^ 0x297a2d39))` followed
  by `fmix32(low ^ fmix32(n ^ 0x2c1b3c6d))`, 16 lowercase hexadecimal digits
  (`fmix32` is MurmurHash3's finaliser).

It is never taken over display coordinates: those are centred and scaled by
the bounds of the cells an embedding was first read for, a frame an opening
grown or cut in place keeps, so they depend on how the cells were opened.

The reader worker sums an embedding as it reads it (`EmbeddingPayloads`,
`embedding-payload.js`), in blocks of 4,096 open rows; a whole read is summed
after its payload is handed to the page, 262,144 rows a task, so the page's
next requests are answered between slices (anything that needs the sums
first -- a digest, a growth, an opening cut in place -- takes the rest at once): a growth step adds only
the rows it reads, an opening cut in place keeps the sums of the whole blocks it
keeps and reads the fewer than 4,096 rows of the block it cuts through when the
digest is asked for, and the same cells in another order keep the total. A
digest of an embedding not read in the opening reads its stored values (no
display coordinates). Saving a session therefore folds nothing on the page. A
layout held only in memory (a generated dataset, every cell open) is digested
over the coordinates it holds.

The fingerprint is taken with the open cells in use (`useOpenCells`,
`data/open-cells-work.js`): a growth step or **Apply** waits until the session
is written, so the fingerprint, the `cells` record and every per-cell chunk
describe the same cells.

### Why `cellOrder` exists

Selections and categorical codes persist as raw dataset row indices. The three scalars all survive a row permutation: write
a dataset again at the same id from re-sorted input and `datasetId`, `cellCount`, and
`varCount` are unchanged, so every restored index would denote
a different cell with nothing to detect it.

Observation and variable names were the alternative and do not work here — they
are invariant under a row permutation, which is the change being caught.
Coordinates are what a re-ordered export permutes: with each cell named by its
original position, another cell at a position changes that position's hash.

`cellOrder.embedding` is part of the record because a digest is comparable only
within one embedding's columns. A restore verifies the digest of the very
embedding the session names, over the cells the session names (opened first),
whichever embedding is on screen and without loading it for display, so a
session saved on the 2D t-SNE opens while a 3D UMAP is shown. A capture
re-checks every fingerprint field but the digest between contributors
(`datasetFingerprintStillCurrent()`); the digest stands while the open cells,
the embedding and the dataset do.

### Validation and comparison

`assertDatasetFingerprint(value, context)` validates one record from untrusted
input: exact six-key set, nullable strings for the id and name, safe integers
for the counts, a `dataMapping` that is valid and already canonical, and
`cellOrder` with exactly an embedding reference (`{key, columns}`),
`dimension` in 1–3 and a `/^[0-9a-f]{16}$/` digest. `session-serializer.js`
calls it from `validateManifest()`, so a malformed record is refused before any
contributor is consulted. Files of the earlier format never reach it: their
signature is refused first.

`datasetFingerprintMatches(a, b)` validates both records and then compares
every field but `name`. `verifySessionFingerprint(ctx, saved)` returns `null` when the
saved record can be restored here and otherwise the sentence the session
controls surface verbatim; `_restoreFromBundleSource()` throws that string as a
`RangeError`.

Every refusal is distinguished, and attributing one to the wrong cause would be
its own integrity failure — telling a user their data was re-ordered when the
session only named another embedding would make them distrust a sound dataset:

| Cause | Detected by | Message |
|---|---|---|
| A session of another format version | the signature, in `bundle/reader.js` | `otherSessionFormatMessage()`: "This session was saved by an earlier version of Cellucid (session format 2): the check that it marks the same cells cannot be verified, so it was not opened. Re-create it in this version of Cellucid." (a newer format: "…saved by a newer version of Cellucid (session format 4), which this version cannot read. Open it in that version of Cellucid.") |
| A file without a session signature | the signature, in `bundle/reader.js` | `NOT_A_SESSION_FILE_MESSAGE`: "This file is not a Cellucid session file of format 3 (it does not start with CELLUCID_SESSION/3)." |
| A different dataset | the counts differ, or the counts agree but the dataset id differs | names the saved and current cell and gene counts, or both identities (`"<name>" (id "<id>")`), and asks for the dataset the session was saved on |
| An embedding this dataset cannot show | an enabled embedding of the saved mapping is neither enabled now nor one the dataset can show | "This session shows obsm['X_tsne'] columns [0, 1], which this dataset does not have." (or the reason it cannot be plotted) |
| A different neighbour graph, modality, gene names, expression matrix, log1p, highly variable genes, latent space or vector fields | that group of `dataMapping` differs | names the group and both values |
| An embedding that is gone | `cellOrder.embedding` names no layout of this dataset | names the embedding |
| A re-ordered dataset | `cellOrder.digest` differs, everything else matches | says the dataset has the same name and the same cell and gene counts but its cell coordinates differ, so every saved selection could mark the wrong cells; either the cells are listed in a different order, or the same cells were exported again from a re-computed embedding |
| Fields or genes the session uses are missing | `describeMissingRequiredFields()` after the session's mapping is applied | names every missing observation field and gene at once |
| A category highlight whose category changed | `decodeCellSet()` while the membership chunk restores | "The category this selection was saved relative to no longer holds the same cells." |

A session whose enabled embeddings or defaults differ from the open dataset's,
but which the dataset can show, is not refused: `restoreSessionDataMapping()`
turns its embeddings on beside the current ones, with its defaults, before any
chunk restores a view, and applies the session's mapping exactly at commit,
once the restored views no longer show the extras; its rollback restores the
previous mapping.

None of these messages carry identifiers or internal vocabulary.

`datasetDependent: false` remains useful profile metadata for the
dockable-layout owner, but it does not authorize layout-only salvage from a
mismatched bundle. The entire operation rejects and rolls back. This applies to
`restorePublishedDefaultState()` too: the published-default profile shares
`validateManifest()`, the same fingerprint comparison and the same `requires`
check, so an advertised
`default.cellucid-session` must carry `cellOrder` and `dataMapping`, and its
`cellOrder` must name an embedding the dataset has.

The fingerprint is still not a content hash: the digest covers one embedding's
coordinates, not expression or metadata. Dataset publishers must assign a new identity when
scientific content changes. What `cellOrder` adds is that a re-ordered
republication under an unchanged identity is now detected instead of trusted.

## Official sample profile

Catalog-advertised `default.cellucid-session` uses a separate exact profile
(`PUBLISHED_DEFAULT_STATIC_CHUNK_PROFILE` in `session/published-default.js`):

1. `core/field-overlays`
2. `core/state`
3. `ui/dockable-layout`
4. `ui/figure-export`
5. `analysis/windows`
6. `analysis/docked`
7. `highlights/meta`

Only `cinematic/camera` is dropped when publishing, and only when it is empty;
the trim (`trimToPublishedDefault()`) keeps `requires` and rebinds the
fingerprint, `name` included, to the published identity. The catalog's
`session` entry advertises the file's path, its SHA-256 and its exact byte
length (`bytes`) together. A published default has no size limit:
the read stops past the declared length, each chunk inflates only to its own
declared `uncompressedBytes`, and the loader verifies the manifest, the length
and the digest and
applies this profile only through `restorePublishedDefaultState()` after the
scientific dataset is published and while the dataset selection generation is
still current.

Do not route this seven-chunk artifact through the generic **Load** path. A browser test
for the advertised path must author a real seven-chunk fixture; reusing generic
**Save** output is invalid by construction.

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

- dataset points, observation tables, or full expression matrices;
- connection state beyond the manifest's `source` record (a sample id, a
  link's locator, a server's address, or a local file's name; nothing for
  a Python or R viewer);
- per-cell gene values an analysis read (they are re-read from the dataset);
- the antialiasing preference (a device capability);
- Community Annotation authentication, votes, or repository workflow;
- Benchmarking controls;
- notifications, hover, focus, pending selections, or in-flight requests; or
- DOM, worker, network, and WebGL runtime objects.

Next: {doc}`11_analysis_architecture`.
