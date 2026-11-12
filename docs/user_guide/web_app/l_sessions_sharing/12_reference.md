# Reference: session files (format 4)

**Audience:** power users and developers  
**Time:** 20–40 minutes
**What you'll learn:**
- the exact `.cellucid-session` framing and manifest;
- which chunks a session file holds, and how its values and arrays are stored;
  and
- how a catalog's starting view differs from a session you save.

## File and framing

**Share…** (under **State › Session**) downloads one binary file with the
`.cellucid-session` extension. It holds the history of the open dataset as
the person chose it (`scope`, below): the trackers and branches ticked, each
a tree of steps, and every value those steps refer to; or the current state
only; or the steps only. Its
framing is:

1. the 19-byte ASCII signature `CELLUCID_SESSION/4\n`: the file kind and
   its session format version, `4`;
2. `manifestByteLength` as an unsigned 32-bit little-endian integer;
3. exactly that many UTF-8 JSON manifest bytes; and
4. for each manifest chunk entry, in order, a 32-bit little-endian stored
   length followed by exactly that chunk's stored bytes.

A frame length is a u32, so the manifest and each stored chunk hold at most
4 GiB − 1 bytes; that is the format's only size limit.

This release accepts exactly one session format, the one in the signature. A
session of another format version is refused before anything else is read,
naming it. One saved by an earlier Cellucid (format 3 held one state in
per-feature chunks, not the history of the steps; formats 1 and 2 before it) is
refused with *This session was saved by an earlier version of Cellucid (session
format 3), which this version no longer reads: a session file now holds the
history of every step (session format 4). It was not opened. Open it in the
version of Cellucid that saved it, or re-create the session in this version.*;
one saved by a newer Cellucid with *This session was saved by a newer version
of Cellucid (session format 5), which this version cannot read. Open it in that
version of Cellucid.* (each with its own format number). Any other file that
does not start with `CELLUCID_SESSION/4\n` is refused with: *This file is not a
Cellucid session file of format 4 (it does not start with
CELLUCID_SESSION/4).*

There may be no duplicate chunks, unlisted payload, truncated payload, or
trailing bytes. Cellucid treats every file as untrusted input and checks the
framing, manifest, chunk inventory, byte lengths, gzip streams, digests,
hashes and dataset identity before reporting success; a file that breaks the
format is refused as a damaged file, the notification's **Details** starting
*This session file is damaged:* (framing faults: *Invalid session file*).

## Exact manifest

The root has exactly `createdAt`, `app`, `actions`, `scope`,
`datasetFingerprint`, `source`, `requires`, `chrome`, `activeTracker`,
`trackers` and `chunks`, plus `cells` and `genes` when the session was saved on a subsample
or with genes read for some of the open cells. `app` is what wrote it (the
Cellucid version and the SHA-256 of its action catalog), `actions` one stamp
per action type and parameters version its steps use, and `scope` what the
save or share chose ({doc}`../../python_package/g_api_reference_coverage/api/sessions`
says what each choice keeps). The complete layout
is the web app's `assets/js/app/session/README.md` ("The file"), in the
`cellucid` repository; this page shows its shape:

```json
{
  "createdAt": "2026-11-03T12:00:00.000Z",
  "app": { "version": "0.9.7", "catalog": "<sha256 of the action catalog>" },
  "actions": [
    { "type": "category.merge", "version": 1, "params": "<sha256 of its parameter schema>" }
  ],
  "scope": {
    "trackers": ["t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10"],
    "branches": [{ "tracker": "t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10", "branches": ["b0"] }],
    "content": "history", "results": "include", "notes": true, "chats": false,
    "journal": false, "times": "exact", "source": "full"
  },
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
  "chrome": { "pressed": { "highlight-mode": "highlight-mode-lasso" } },
  "activeTracker": "t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10",
  "trackers": [
    {
      "id": "t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10",
      "name": "Annotation",
      "dataset": {
        "id": "suo", "name": "Suo et al. 2022",
        "cellOrder": { "embedding": { "key": "X_umap", "columns": [0, 1] }, "dimension": 2, "digest": "9f2c41a0b7d3e185" }
      },
      "createdAt": 1793725004466,
      "by": { "origin": "user", "agent": null, "model": null, "host": null, "person": null, "turn": null },
      "start": { "kind": "opened", "from": null, "report": null, "mapping": null },
      "removed": null,
      "rootId": "s0",
      "headId": "s7",
      "activeBranch": "b0",
      "mainBranch": "b0",
      "paused": null,
      "branches": [
        { "id": "b0", "name": "Main", "base": null, "tip": "s7", "createdAt": 1793725004466,
          "by": { "origin": "user", "agent": null, "model": null, "host": null, "person": null, "turn": null },
          "lane": 0, "locked": false, "archived": false, "notes": [] }
      ],
      "groups": [],
      "steps": ["steps/0/0"],
      "turns": [],
      "notes": [],
      "journal": null,
      "chats": [],
      "nextStep": 8, "nextBranch": 1, "nextGroup": 0, "nextNote": 0, "nextJournal": 0
    }
  ],
  "chunks": []
}
```

- `createdAt` is a canonical ISO timestamp.
- `chrome` is page chrome kept outside the history, which the web app puts back
  after a load: `{"pressed": {<group id>: <button id>}}`, the pressed button of
  each pressed-button group (the selection tool: `{"highlight-mode":
  "highlight-mode-lasso"}`); `{"pressed": {}}` from a file written without a
  page. Every form of file keeps it (history, state, steps). A group id is
  non-empty and names its button by a non-empty id; anything else is refused.
- `source` says where the dataset was opened from, so a load can open
  it first: `{"kind": "sample", "dataset"}` for a sample of the page's
  catalog, `{"kind": "link", "locator", "dataset"}` for a data link (the
  `locator` is the `?data=` value in its public form, never a key, signature
  or token the link carried; see {doc}`../b_data_loading/14_open_from_link`),
  `{"kind": "remote", "url", "dataset"}` for a Cellucid server,
  `{"kind": "file", "name", "dataset"}` for a file or folder opened from this
  computer (`name` is its file or folder name), or `null` for a notebook. A
  local file is opened from the session only where this browser kept its
  handle (Chrome and Edge, {ref}`reopen-local-file`, keyed by the record's
  `dataset`): a load reopens it at once while the browser still lets it be
  read, else offers **Reopen <name>**, whose click asks the browser and then
  loads the session. Anywhere else a load goes onto whatever is open and,
  when that is another dataset, the refusal names the dataset it needs. The
  handle is never part of the session. A starting view's `source` is always
  `null`. The source is only a pointer: the fingerprint, not the source,
  decides whether the session applies.
- `datasetFingerprint` has exactly `datasetId`, `name`, `cellCount`,
  `varCount`, `cellOrder`, and `dataMapping`. `name` is the dataset's display
  name (or `null`); it is informational (messages name a dataset by it) and is
  never compared.
- `requires` lists, sorted and without repeats, every dataset observation field
  (`obs`) and gene (`var`) the shown step uses: what its views are colored by,
  filter on or have changed the colors of, and the fields its highlight
  groups are bound to, by the
  dataset's own key (a field renamed in the session by its original key).
  Fields the session creates itself are not listed. Before anything is put
  back, Cellucid refuses a dataset missing any of them, naming every missing
  key at once: *This session uses the observation fields "…" and the genes
  "…", which the dataset that is open now does not have.* Genes are compared
  only when the session names genes the way the open dataset does (the same
  gene-name column and modality).
- `cells` (present only for a subsample) is the subsample's specification
  (`spec`), its kept original positions (`kept`, cell-set bytes in base64) and
  their digest; a load opens exactly those cells first. A manifest without
  `cells` was saved on every cell.
- `genes` (present only when genes were read for some of the open cells) is
  how many of the open cells, in which order, and their digest; a load reads
  genes for at least those cells again.
- `cellOrder` has exactly `embedding`, `dimension`, and `digest`. `embedding` is
  `{"key", "columns"}`: an AnnData `obsm` key with its 0-based columns (a list
  of integers, never null): the embedding the live view showed when the
  session was saved. `dimension` is its dimension, exactly `1`, `2`, or `3`.
  `digest` is 16 lowercase hexadecimal characters computed over the stored
  values of that embedding's columns for the cells the session was saved on,
  each cell named by its original position (MurmurHash3 lanes summed over the
  cells; specified in
  {doc}`../p_developer_docs/10_sessions_persistence_and_serialization`), so how
  the cells were opened (whole, grown, a share, priority order) does not
  change it.
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
  (an `obsm` key or null), and `vectorFields`: the vector bindings in use,
  each `{"key", "embedding"}`, the vector's `obsm` key and the enabled
  embedding reference it is drawn over. Bindings are sorted by embedding key,
  columns, then key. The Python and R packages read the same shape; see
  {doc}`../../python_package/g_api_reference_coverage/api/sessions`.
- `activeTracker` is the id of the tracker shown; it must be one of
  `trackers`.
- `trackers` lists every tracker: each is the history model's `Tracker`
  (`assets/js/app/history/model.js`, mirrored to Python and R as
  `history_model.schema.json`) with exactly `id`, `name`, `dataset` (the
  dataset it belongs to), `createdAt`, `by` (who started it), `start` (how:
  opened, shown, a replay or an import), `removed` (null unless it is in
  Recently deleted), `rootId` (its first entry, **Opened**), `headId` (the
  step shown), `activeBranch`, `mainBranch`, `paused` (null while it records,
  else since when, who acted and, never in a file, a draft), `branches` (each
  `{id, name, base, tip, createdAt, by, lane, locked, archived, notes}`),
  `groups`, `steps`, `turns`, `notes`, `journal`, `chats` and the counters
  `nextStep`, `nextBranch`, `nextGroup`, `nextNote` and `nextJournal`. Its
  `steps` lists its steps chunks in the order the file holds them, `journal`
  names its journal chunk (or null) and `chats` its conversation chunks (only
  when the share included them).
- `chunks` is the ordered inventory.

## Chunks

Every chunk entry has exactly `id`, `kind`, `priority`, `codec`,
`storedBytes` and `uncompressedBytes`, plus its kind's own keys:

| Chunk id | Kind (own keys) | Priority | Contains |
|---|---|---|---|
| `steps/<n>/<k>` | `steps` (`tracker`) | eager | a run of the n-th tracker's steps, in id order; the tracker lists its steps chunks, read in order as one list |
| `parts/<priority>/<k>` | `parts` | eager for what the values of the step shown name, else lazy | sub-values stored once (below) |
| `array/<digest>` | `array` (`digest`, `encoding`) | eager for what the values of the step shown name, else lazy | one typed array, stored once however many values hold it |
| `values/<priority>/<k>` | `values` | eager: the values of the step shown; lazy: every other value a step or a journal entry refers to | values, as skeletons |
| `journal/<n>` | `journal` (`tracker`) | lazy | the n-th tracker's journal (Recently deleted), when included |
| `chat/<n>/<m>` | `chat` (`tracker`, `conversation`) | lazy | a conversation kept with the n-th tracker, when included |

`priority` is `eager` or `lazy`. Every eager chunk comes before every lazy
one, and within each priority the kinds come in the order of the table
(`steps`, `parts`, `array`, `values`, `journal`, `chat`). `steps`, `parts` and
`values` are JSON cut where their text passes 16 MiB: one list or map in
several chunks, read back as one. `codec` is `gzip` (where gzip saves at least
a sixteenth of the chunk's bytes) or `none` (stored as it is: `storedBytes`
equals `uncompressedBytes`); each chunk is decompressed to exactly its
`uncompressedBytes`. A file holds at least one eager `values` chunk, except a
file of steps only, which holds none; each tracker names `steps` chunks of its
own. Unknown, missing, duplicate, unreferenced or reordered chunks, a value,
part or array nothing refers to and a step referring to a value the file lacks
are refused.

### Steps

Each step is the model's `Step`, with exactly `id`, `parent`, `children`,
`branch`, `kind` (`root`, `action`, `outside` for a change made outside an
action, or `squashed`), `type` (the action), `params` (its parameters, by
name), `v` (the version of its parameters), `code` (its line of pseudocode),
`label`, `by` (`{origin, agent, model, host, person, turn}`: who acted, origin
`user`, `agent`, `host`, `session`, `replay`, `startup` or `outside`), `at` and
`until` (when, in milliseconds since 1970), `merged`, `ms`, `touches` (the
parts it changed), `before` and `after` (each part it changed, mapped to the
content hash of its value), `reads` (the hash of each input it read), `status`
(`applied`, or `unapplicable` for a step that no longer applies), `reason`,
`outcome`, `coalesceKey`, `replayOf`, `squashed` (for a squashed step: why,
how many steps it replaced, who acted and what they carried), `notes`,
`bookmark` and `group`. `params` and `outcome` are skeletons, like values.
Each tracker's steps must form one tree with one root, parents and children
agreeing, each branch one line from its fork step to its tip, and the step
shown on the active branch's path; the full rules are the README's and
`model.js`'s.

The parts a step can change are the sixteen state domains, in the order a load
puts them back: `cells`, `mapping`, `expression`, `field-overlays`, `gene-sets`, `views`,
`active-field`, `filters`, `colors`, `visibility`, `highlights`, `camera`,
`camera-path`, `settings`, `analysis`, `layout`
({doc}`02_what_gets_saved_and_restored` says what each holds). The value of a
part at a step is the `after` value of the nearest step on the way from the
root that changed it, or the root's. A part this version does not have is
left out of every step when the file is opened, and named.

### Values

A `values` chunk is a JSON object mapping each value's content hash (64
lowercase hexadecimal characters) to its **skeleton**: the value as JSON with
every typed array replaced by `{"$array": "<digest>"}`, every number JSON
cannot hold by `{"$number": "NaN" | "Infinity" | "-Infinity" | "-0"}`, and
every object or list whose skeleton text is at least 512 characters by
`{"$part": "<name>"}`, its skeleton stored once in a `parts` chunk under the
SHA-256 of its text (a part names the parts inside it). Values that share a
sub-value (the gene names of every analysis result, a highlight group) store
it once. Every value is checked against its hash once its parts and arrays are
joined back in: a part's name is never trusted. A part that names itself
(through others), and a value whose text would pass 2^29 − 24 characters, are
refused before anything is expanded. There is
no encoder per feature: a part's value is stored as it is. The hash's
canonical form is specified in
{doc}`../p_developer_docs/10_sessions_persistence_and_serialization`.

### Arrays

An array's digest is `<type>:<length>:<sha256 of its bytes>`, for example
`Uint32Array:5:<64 hex characters>`. Its `encoding` is `raw` (its elements,
little-endian) or `cell-set`, used for a strictly increasing `Uint32Array`
(the cells of a set) when that is smaller: a uvarint bound (the last cell plus
one), then the cell-set bytes, whose byte 0 names the encoding:

| Byte 0 | Encoding | Body |
|---:|---|---|
| 0 | sorted delta | uvarint count, then that many uvarint gaps (the first is the index itself) |
| 1 | bitset | one bit per cell below the bound, `ceil(bound / 8)` bytes, cell `i` in bit `i & 7` of byte `i >> 3` |
| 2 | runs | uvarint run count, then per run a uvarint gap from the end of the previous run and a uvarint length |

The writer writes whichever is smallest. Every array is checked against its
digest after it is decoded. A *uvarint* is an unsigned LEB128 integer.

## Exact dataset identity

Before anything is put back, Cellucid checks, in order: the cell and gene
counts, the dataset id, the data mapping, and the cell-order digest; then,
under the session's own mapping, the `requires` field list. A dataset id,
cell count, gene count, cell-order embedding, or cell-order digest that
differs refuses the whole load. A data mapping may differ: a session whose
embeddings or other choices (graph, modality, gene names, matrix, …) differ
is put back onto them when the dataset can show or use them (they are put
back right after the cells, before the views, and rolled back with the rest),
and refused, naming what is missing, when it cannot; a group the dataset
cannot change must equal the one it is shown with. The digest is checked on
the embedding `cellOrder.embedding` names, over the cells the session was
saved on (opened first), whichever embedding is on screen; it is not loaded
for display. Cellucid does not put back the layout alone and does not skip the
parts that name cells. Any departure from the fingerprint's six-key shape is a
format violation.

This strict equality keeps a session's per-cell filters, fields and highlights
from being attached to a different dataset, or to a re-ordered republication
of the same one. The scalars are all preserved when a dataset is written again
at the same id from re-sorted input, so `cellOrder` is what makes that case
detectable. It is still your responsibility not to reuse the same identity for
scientifically changed content that leaves the coordinates alone.

## Load order and outcome

A load is one awaited operation, all or nothing:

1. the session's dataset is opened if needed, with the cells it was saved on;
2. the signature, manifest, fingerprint and `requires` are checked;
3. the eager chunks are read and the step shown is put back, part by part in
   domain order, changing only the parts whose value differs; a part that
   cannot be put back is left as it was shown and the others go on;
4. the lazy chunks are read while that step is shown; a new step waits; and
5. the file's trackers are added beside the history's, and the load ends in
   one notice reading **Session fully restored.**: titled **Loaded "…" beside
   yours**, with **Undo**, when trackers of yours are beside them, else the
   progress notification finishing as **Loaded**.

A part that could not be put back ends the load in **Session restored but
for N parts** instead, naming each part (a part saved by another version of
Cellucid is said as that), with each failure's words under **Details**. Any
other failure, and **Cancel**, puts back what was shown and keeps the
trackers there were. Starting a newer load cancels the older one, which is dismissed
without a false failure or success. A camera path saved with **Autoplay** on
starts once the load is done. The history kept in this browser
({doc}`../m_state_history/03_saved_states_and_autosave`) is read the same way
when you open a saved state.

## A catalog's starting view

A catalog may name a SHA-256-pinned starting view for a dataset (the `session`
entry of its catalog, such as `pancreas.cellucid-session`). Cellucid fetches
it and puts it back by itself only through the catalog, as described in
{doc}`04_official_sample_states`. It is a format-4 file with this profile,
checked on the manifest before any chunk is read:

- exactly one tracker, whose root is its step shown and its only step (no
  step, no branch), with no conversation;
- every chunk eager;
- `source` is `null`, and there is no `cells` and no `genes`.

Once its one step is read, its root must hold neither `cells` nor `layout`. It
is put back without a step, before the dataset's first tracker starts, so the
history starts from it; it applies to any subsample of its dataset (on a
subsample the cell-order digest is not compared, as the catalog pins the
file). Otherwise it goes through the same manifest checks, fingerprint
comparison and `requires` check as any session, so it too must carry
`cellOrder` and `dataMapping`.

A shared session file becomes a starting view only once it is published:
`trimToPublishedDefault` (in `assets/js/app/session/published-default.js` of
the web app, used by `cellucid-datasets/scripts/publish_session_preset.mjs`)
makes the head of the active tracker the one root, drops every step, and
binds the fingerprint to the identity the catalog serves the dataset under. It
refuses a capture whose steps made content (highlights, gene sets, a camera
path, cells hidden one by one, field and category edits).

## Implementation pointers

In the web app repository:

- The format (manifest, chunks, checks): `assets/js/app/session/history-file.js`
- Framing and format version: `assets/js/app/session/bundle/format.js`
- Save and load: `assets/js/app/session/session-serializer.js`
- Dataset fingerprint and `requires`: `assets/js/app/session/session-context.js`
- Starting views: `assets/js/app/session/published-default.js`,
  `assets/js/app/session/catalog-default-session.js`
- Values as skeletons and arrays: `assets/js/app/history/value-parts.js`,
  `assets/js/app/history/value-pool.js`
- Cell-set encodings: `assets/js/app/session/codecs/cell-set.js`
- The state domains: `assets/js/app/state/domains/README.md`
- Maintainer notes: `assets/js/app/session/README.md`
