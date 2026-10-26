# Reference: current session bundles

**Audience:** power users and developers  
**Time:** 20–40 minutes
**What you’ll learn:**
- the exact `.cellucid-session` framing and manifest;
- which chunks a current user-authored bundle contains; and
- how ordinary sessions differ from integrity-pinned official sample states.

## File and framing

**Save** (under **Session file**) downloads one binary file with the
`.cellucid-session` extension. Its framing is:

1. the 19-byte ASCII signature `CELLUCID_SESSION/2\n` — the file kind and
   its session format version, `2`;
2. `manifestByteLength` as an unsigned 32-bit little-endian integer;
3. exactly that many UTF-8 JSON manifest bytes; and
4. for each manifest entry, a 32-bit little-endian stored length followed by
   exactly that chunk's stored bytes.

A file that does not start with `CELLUCID_SESSION/2\n` is refused before
anything else is read, with: *This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2).*

There may be no duplicate chunks, unlisted payload, truncated payload, or
trailing bytes. Cellucid treats every bundle as untrusted input and validates
the complete framing, inventory, byte lengths, gzip stream, and dataset
identity before reporting success.

## Exact manifest

The root has exactly `createdAt`, `datasetFingerprint`, `source`, `requires`,
and `chunks`:

```json
{
  "createdAt": "2026-07-27T12:00:00.000Z",
  "datasetFingerprint": {
    "datasetId": "suo",
    "name": "Suo et al. 2022",
    "cellCount": 1000,
    "varCount": 2000,
    "cellOrder": {
      "embedding": { "key": "X_umap", "columns": [0, 1] },
      "dimension": 2,
      "digest": "9f2c41a0b7d3e185"
    },
    "dataMapping": {
      "kind": "anndata",
      "embeddings": {
        "enabled": [
          { "key": "X_pca", "columns": [0, 1] },
          { "key": "X_umap", "columns": [0, 1] },
          { "key": "X_umap_3d", "columns": [0, 1, 2] }
        ],
        "defaults": {
          "1d": null,
          "2d": { "key": "X_umap", "columns": [0, 1] },
          "3d": { "key": "X_umap_3d", "columns": [0, 1, 2] }
        }
      },
      "knn": "connectivities",
      "geneNames": "index",
      "expression": "X",
      "log1p": {},
      "highlyVariable": null,
      "latent": null,
      "vectorFields": []
    }
  },
  "source": { "kind": "sample", "dataset": "suo" },
  "requires": { "obs": ["cell_type", "organ"], "var": ["CD3E"] },
  "chunks": []
}
```

- `createdAt` is a canonical ISO timestamp.
- `source` says where the dataset was opened from, so **Load** can open
  it first: `{"kind": "sample", "dataset"}` for a sample of the page's
  catalog, `{"kind": "link", "locator", "dataset"}` for a data link (the
  `locator` is the `?data=` value, see {doc}`../b_data_loading/14_open_from_link`),
  `{"kind": "remote", "url", "dataset"}` for a Cellucid server,
  `{"kind": "file", "name", "dataset"}` for a file or folder opened from this
  computer (`name` is its file or folder name), or `null` for a
  notebook. A local file cannot be opened from the session: **Load** applies
  the session to whatever is open and, when that is another dataset, the
  refusal names the dataset it needs. A published sample view's `source` is
  always `null`. The source is only a pointer: the fingerprint, not the
  source, decides whether the session applies.
- `datasetFingerprint` has exactly `datasetId`, `name`, `cellCount`,
  `varCount`, `cellOrder`, and `dataMapping`. `name` is the dataset's display
  name (or `null`); it is informational — messages name a dataset as
  `"<name>" (id "<id>")` — and is never compared.
- `requires` lists, sorted and without repeats, every dataset observation field
  (`obs`) and gene (`var`) a view shows or filters on or a highlight group is
  bound to, by the dataset's own key (a field renamed in the session by its
  original key). Fields the session creates itself are not listed. After the
  session's data mapping is applied and before any chunk is, Cellucid refuses a
  dataset missing any of them, naming every missing key at once: *This session
  uses the observation fields "…" and the genes "…", which the dataset that is
  open now does not have.*
- `cellOrder` has exactly `embedding`, `dimension`, and `digest`. `embedding` is
  `{"key", "columns"}`: an AnnData `obsm` key with its 0-based columns (a list
  of integers, never null) — the embedding the live view
  showed when the session was saved. `dimension` is its dimension, exactly `1`,
  `2`, or `3`. `digest` is 16 lowercase hexadecimal characters computed over
  every cell coordinate of that embedding in dataset row order.
- `dataMapping` records the arrays the dataset was shown with. Its `kind` is
  always `"anndata"` (any other `kind` is refused: `kind must be "anndata"`),
  and it has exactly `kind`, `embeddings` (`enabled`: the embedding
  references that are turned on, at least one, sorted by key then columns
  without repeats; and `defaults`: the `1d`/`2d`/`3d` default, each an enabled
  reference or null, at least one named), `knn` (an `obsp` key or null),
  `modality` (only for a dataset with several modalities, such as a MuData
  file: the modality's name, whose genes the next four keys then describe),
  `geneNames` (`"index"` or `"column:<name>"`), `expression` (`"X"`, `"raw"`,
  `"layer:<name>"`, or null), `log1p` (for each expression matrix whose
  values are `log1p`, keyed `"X"`, `"raw"` or `"layer:<name>"` in code-point
  order, `{"base": null}` for e or `{"base": <number>}`; a matrix not listed
  is as stored), `highlyVariable` (the true/false var column marking the
  highly variable genes of the expression matrix's axis, or null), `latent`
  (an `obsm` key or null), and
  `vectorFields`: the vector bindings in use, each `{"key", "embedding"}` —
  the vector's `obsm` key and the enabled embedding reference it is drawn
  over. Bindings are sorted by embedding key, columns, then key. The Python package reads the
  same shape; see
  {doc}`../../python_package/e_jupyter_hooks/11_session_bundles_apply_session_to_anndata`.
- `chunks` is the ordered inventory.

Every chunk entry has exactly these nine keys:

```json
{
  "id": "core/state",
  "contributorId": "core-state",
  "priority": "eager",
  "kind": "json",
  "codec": "gzip",
  "label": "Core state",
  "datasetDependent": true,
  "storedBytes": 1234,
  "uncompressedBytes": 5678
}
```

`priority` is either `eager` or `lazy`, `kind` is `json` or `binary`, and
`codec` is `none` or `gzip`. The current chunk object has exactly the nine keys
shown above; an extra key such as `dependsOn` is invalid.

Eager and lazy are deterministic scheduling classes inside one awaited restore:
all eager chunks precede all lazy chunks, and contributor groups retain their
registered order. **Session fully restored** is emitted only after both classes
have validated, applied, and committed.

## Exact dataset identity

Before anything is applied, Cellucid checks, in order: the cell and gene
counts, the dataset id, the data mapping, and the cell-order digest; then,
under the session's own mapping, the `requires` field list. A dataset id,
cell count, variable count, cell-order embedding, or cell-order digest that
differs rejects the whole restore. A data mapping may differ: a session whose
embeddings or other choices (graph, modality, gene names, matrix, …) differ
is restored onto them when the dataset can show or use them (they are applied
before the first chunk, and rolled back with the rest), and refused, naming
what is missing, when it cannot; a group the dataset cannot change must equal
the one it is shown with. The digest is checked on the embedding
`cellOrder.embedding` names, which is loaded for the check whichever embedding
is on screen. Cellucid does not apply a layout-only subset and does not skip
dataset-dependent chunks. Any departure from the fingerprint's six-key shape
is a schema violation.

This strict equality prevents a session's cell-indexed filters, categorical
codes, or highlights from being attached to a different dataset, or to a
re-ordered republication of the same one. The scalars are all preserved when a
dataset is written again at the same id from re-sorted input, so `cellOrder`
is what makes that case detectable. It is still your responsibility not to
reuse the same identity for scientifically changed content that leaves the
coordinates alone.

## Current user-authored inventory

Every bundle an ordinary **Save** writes contains every singleton chunk, in
this order (eager user-defined code chunks, when there are any, sit between the
first two), including an explicit empty camera path:

| Chunk id | Priority | Kind | Dataset dependent | Contains |
|---|---:|---|---:|---|
| `core/field-overlays` | eager | JSON | yes | field rename/delete registries and user-defined field metadata |
| `cinematic/camera` | eager | JSON | yes | exact keyframes, playback settings, and explicit empty path |
| `core/state` | eager | JSON | yes | every sidebar control, the live view's embedding, filters and active field, and the multiview graph: layout, link state, the live view's camera and centroid toggles, and each kept view's embedding, camera, active field, filters, and centroid toggles |
| `ui/dockable-layout` | eager | JSON | no | non-analysis floating-panel state and geometry |
| `ui/figure-export` | eager | JSON | no | every Figure Export setting (size, frame, background, fonts, labels, legend, preview, format, DPI, SVG strategy, all-views) |
| `analysis/windows` | eager | binary | yes | floating analysis windows: mode, geometry, settings, and result |
| `analysis/docked` | eager | binary | yes | the analysis panel in the sidebar: which analysis is open, and each analysis's settings and result |
| `highlights/meta` | eager | JSON | yes | highlight pages and group metadata |

Dynamic binary chunks are closed by those singleton inventories:

| Chunk family | Priority | Required relationship |
|---|---:|---|
| `user-defined/codes/<fieldId>` | eager or lazy | exactly one for every categorical user-defined field; eager exactly when a target live/snapshot active field needs it |
| `highlights/cells/<groupId>` | lazy | exactly one membership payload for every advertised highlight group, in group order |

The built-in singleton profiles and dynamic metadata are exact. Missing,
duplicate, reordered, aliased, unknown, or dishonestly labelled chunks reject
the operation. Empty replacement state is represented by a real singleton
chunk; omission never means “keep whatever is already open.”

What a session records is the choice, not what the dataset already holds: a
kept view's label and filter summary, a highlight group's field and category
positions, and the analysis gene caches are recomputed from the dataset. The
antialiasing setting is a property of the machine drawing the view and is not
saved; the app theme travels with the session.

## Payload encodings

All numbers are little-endian; a *uvarint* is an unsigned LEB128 integer.
Chunks with codec `gzip` are one gzip member around the bytes below.

**Highlight cells** (`highlights/cells/<groupId>`): byte 0 names the encoding,
and **Save** writes whichever is smallest for that selection.

| Byte 0 | Encoding | Body |
|---:|---|---|
| 0 | sorted delta | uvarint count, then that many uvarint gaps (the first is the index itself) |
| 1 | bitset | one bit per dataset cell, `ceil(cellCount / 8)` bytes, cell `i` in bit `i & 7` of byte `i >> 3` |
| 2 | runs | uvarint run count, then per run a uvarint gap from the end of the previous run and a uvarint length |
| 3 | category subset | uvarint category size, u32 category digest, then one bit per cell of the category, in ascending order |
| 4 | whole category | uvarint category size, u32 category digest |

Encodings 3 and 4 are used for a category highlight on a dataset observation
field: the selection is stored relative to the cells of that category, so a
whole category costs a few bytes. On restore the category's cells are
recomputed, and a category whose size or digest differs is refused (*The
category this selection was saved relative to no longer holds the same
cells.*) rather than selecting other cells.

**User-defined codes** (`user-defined/codes/<fieldId>`): a byte naming the
encoding, the uvarint cell count, then raw codes at the field's code width
(`0`), runs of uvarint value and length (`1`), or codes packed at 1, 2 or 4
bits each (`2`, preceded by one byte giving the width). **Save** writes the
smallest.

**Analysis** (`analysis/windows`, `analysis/docked`): typed JSON — a u32 JSON
length, the JSON with every numeric column replaced by `{"$typed": k}`, then
each column as a dtype byte, a uvarint length, and its raw elements. A result
column of 8,192 genes costs its 4- or 8-byte values instead of JSON text.

## Restore order and speed

Restore applies the eager chunks, then the lazy ones, in file order. The core
state rebuilds the views in one pass: every field, gene and embedding any view
shows or filters on is requested in the same step, so the reader can serve the
genes together; each kept view is rebuilt once, and the live view's own state is
applied last. Between lazy chunks the restore yields to the browser with a
macrotask rather than an animation frame, so a session loaded in a background
tab or an off-screen notebook frame finishes instead of waiting for the page
to be shown.

## Atomic restore and cancellation

Restore is one transaction:

1. validate framing, manifest, profiles, dataset identity, and declared bounds;
2. capture the small reversible owner state;
3. apply eager and lazy chunks in order;
4. prepare and commit all contributors; and
5. refresh the final UI before publishing success.

Any validation, decode, contributor, commit, or final-refresh failure rolls
back registered state, including categorical code references, camera-path
playback state, panels, figure-export settings, and the docked and floating
analyses.
Starting a newer restore or pressing **Cancel** aborts the older operation and
dismisses its progress without a false failure or success.

## Official sample state is a separate profile

An official sample may advertise an integrity-pinned default session (the
`session` entry of its catalog, such as `pancreas.cellucid-session`). Cellucid fetches and applies it automatically only
through the catalog capability described in
{doc}`04_official_sample_states`. Its exact profile contains these seven gzip
chunks, in order:

1. `core/field-overlays`
2. `core/state`
3. `ui/dockable-layout`
4. `ui/figure-export`
5. `analysis/windows`
6. `analysis/docked`
7. `highlights/meta`

It deliberately contains no camera path and no cell memberships. Do not download
that internal default and choose **Load**: the ordinary loader requires
the complete user-authored profile and rejects it. To create a portable manual
session for an official sample, choose the sample, wait for its verified
starting view, then use **Save**.

Only the chunk profile differs. The published default is read through the same
manifest validation, fingerprint comparison and `requires` check as any other
bundle, so it too must carry `cellOrder`, and its cell-order embedding must be one the
dataset publishes.

## Implementation pointers

- Orchestration and exact profiles:
  `cellucid/assets/js/app/session/session-serializer.js`
- Bundle framing and format version:
  `cellucid/assets/js/app/session/bundle/format.js`
- Encodings: `cellucid/assets/js/app/session/codecs/` (`cell-set.js`,
  `typed-json.js`, `varint.js`) and
  `cellucid/assets/js/app/session/contributors/user-defined-codes.js`
- Feature contributors:
  `cellucid/assets/js/app/session/contributors/`
- Maintainer inventory:
  `cellucid/assets/js/app/state-serializer/README.md`
