# Sessions: persistence and serialization

This page documents the one current `.cellucid-session` writer/reader
contract: session format 4, in which the file **is the history** of the open
dataset (its trackers, their steps, and every state value the steps refer to),
its load transaction, the starting views made from it, and the extension
discipline.

**Audience:** contributors changing persistent state, data-source owners, and
maintainers serving official starting views

**Time:** 35–60 minutes
**Prerequisites:** {doc}`29_action_layer_domains_and_history` (actions, state
domains, the history and its value pool) and
{doc}`06_state_datastate_and_events`

## Design invariants

- A session stores the dataset's history, not the dataset: every tracker
  with its step tree, every state domain value its steps refer to (by content
  hash), and, only when the save includes them, the conversations kept with a
  tracker. Loading it restores every tracker and shows the active one's head.
- The generic reader accepts exactly one session format — the version named
  in the file's signature, currently 4 — and restores the complete document
  all-or-nothing. A file of another format version is refused by its
  signature, before anything else is read; there is no migration.
- The file holds state domain values as they are: a value is a JSON skeleton
  plus the typed arrays it names. No feature has an encoder of its own, so a
  new piece of state needs no change here ({doc}`29_action_layer_domains_and_history`,
  "How to add a piece of state").
- Values name fields, genes, categories and views by name and cells by
  **original position** (the cell's row in the dataset as it was written),
  whatever order the cells are open in, so Python and R read a session's
  per-cell state without knowing any layout.
- Those positions are only applied to a dataset whose cells are proven to be
  the ones they were saved on (the dataset fingerprint); an unprovable
  ordering rejects the load.
- Sessions are untrusted input: the manifest, the chunk list, every size,
  every array digest, every value hash and every step tree are exact
  contracts, checked before success.
- A load is one awaited public operation; eager and lazy chunks are an
  internal ordering inside it.
- Cancellation and replacement remain observable to the caller but do not
  produce false success or red product diagnostics.

## Code map

All paths are under `cellucid/assets/js/app/`.

| Module | Owns |
|---|---|
| `session/index.js` | `createSessionSerializer({state, dataSourceManager, getHistory, dispatcher, onLoaded})` |
| `session/session-serializer.js` | `SessionSerializer`: writing (`createSessionBundle`, `downloadSession`), Load (`loadSessionFromFile`, `restoreFromBlob`), a saved state of this browser (`loadSavedHistory`), a session file as a replay plan (`pickSessionFile`, `openSessionFile`), a starting view (`restorePublishedDefaultState`), the manifest alone (`readSessionDataset`), one owned operation at a time |
| `session/share-scope.js` | What a Save or a Share chose (`ShareChoice`, `SAVE_DEFAULTS`, `SHARE_DEFAULTS`, `STATE_DEFAULTS`), made into the trackers a file holds (`shareHistory`) |
| `history/persistence-actions.js` | The catalogue actions that write, list or delete: `history.share` (a session file to download), `history.save` (a saved state of this browser), `history.saved.list`, `history.saved.delete` (one saved state of this browser, by its id; `deletes-saved`), and `history.methods`, `history.pseudocode` |
| `history/anndata-copy-save.js`, `history/exports/anndata-copy.js` | `history.anndata.save` (History ⋯ › **Save AnnData with this state…**; `writes-file`): the state shown, written as a "Current state only" share to name it by its SHA-256, applied to a new Cellucid-ready AnnData `.zarr.zip` of the open dataset by the packages' apply plan ({ref}`history-anndata-save-action`) |
| `session/history-file.js` | The format: manifest and chunks, written (`writeHistoryFile`) and read in two stages (`openHistoryFile`: `readEager`, `readRest`), checked whole (`assertHistoryManifest`), the manifest alone (`readHistoryFileManifest`) |
| `session/bundle/` | The framing: `format.js` (signature, `SESSION_FORMAT_VERSION`, refusal sentences), `reader.js` (streamed, never the whole file in memory), `writer.js` |
| `session/codecs/` | `gzip.js` (bounded native gzip), `cell-set.js` (a set of cells in its smallest encoding), `varint.js` |
| `history/value-parts.js`, `history/value-pool.js` | A value as skeleton plus arrays (`splitValue`, `joinValue`); the content hash (`hashValue`, `typedArrayDigest`) |
| `session/session-context.js` | The dataset fingerprint (`getDatasetFingerprint`, `verifySessionFingerprint`), `requires` (`collectRequiredFields`, `describeMissingRequiredFields`) |
| `session/dataset-source-record.js`, `session/cell-subset-record.js` | Where the dataset opens from (`source`); the cells and genes it was saved on (`cells`, `genes`) |
| `session/published-default.js`, `session/catalog-default-session.js` | Starting views: their profile, the publish step, and the catalog's pinned file |
| `session/view-carry.js` | What a reading of the same dataset with other cells keeps (every domain but the cells, through the domains' owners) |
| `ui/modules/session-controls.js` | Load (State › Session › **Load…** › **Open a session file…**) and opening a saved state, its dataset first |
| `ui/modules/state/share-dialog.js`, `saved-states-dialog.js` | State › Session › **Share…** (what goes in, What's inside, `history.share`) and the **Load…** dialog (saved states; save as history, state or steps, `history.save`; open or replay a session file) |

## Container

A `.cellucid-session` file is:

1. the 19-byte signature `CELLUCID_SESSION/4\n` (`SESSION_FORMAT_VERSION` in
   `bundle/format.js`: the format version is part of the signature);
2. an unsigned 32-bit little-endian manifest length;
3. exactly that many bytes of UTF-8 JSON (the manifest); and
4. for every chunk of the manifest, in manifest order, an unsigned 32-bit
   little-endian stored length and exactly that many stored bytes, then the
   end of the file.

A frame length is a u32, so the manifest and each stored chunk hold at most
4 GiB − 1 bytes (`MAX_FRAME_BYTES`); that is the format's only size limit.
The reader streams the file (`bundle/reader.js`); it never holds it whole.

`sessionSignatureRefusal()` reads the first bytes: a session of another format
version (format 3 and earlier, or a newer one) is refused with
`otherSessionFormatMessage(version)` — "This session was saved by an earlier
version of Cellucid (session format 3), which this version no longer reads: a
session file now holds the history of every step (session format 4). It was
not opened. Open it in the version of Cellucid that saved it, or re-create the
session in this version." (or, for a newer one, "…saved by a newer version of
Cellucid (session format 5), which this version cannot read…") — and any other
file with `NOT_A_SESSION_FILE_MESSAGE`: "This file is not a Cellucid session
file of format 4 (it does not start with CELLUCID_SESSION/4)." A library file
of My datasets (`CELLUCID_LIBRARY/1\n`, the same framing,
`app/library/library-file.js`) is refused with
`LIBRARY_FILE_NOT_A_SESSION_MESSAGE`, which says where it is imported
(Session › Open data › My datasets › All… › Import…). `writeBundle` and
`readBundle` take the file kind's `signature` (`SESSION_SIGNATURE` by
default), so both kinds share one framing.

## Manifest

The manifest has exactly `createdAt`, `app`, `actions`, `scope`,
`datasetFingerprint`, `source`, `requires`, `chrome`, `activeTracker`,
`trackers` and `chunks`, in that order, plus `cells` for a session saved on a subsample and
`genes` for one saved while genes were read for a sample of the open cells.
`assertHistoryManifest` checks it whole before any chunk is read. The layout
is written out in `app/session/README.md` ("The file"); every record is
`app/history/model.js`'s (`historyModelDocument()`, mirrored to Python and R as
`history_model.schema.json`).

| Key | Holds |
|---|---|
| `createdAt` | The save's time, a canonical ISO string (`toISOString()` gives it back unchanged) |
| `app` | `{version, catalog}`: the Cellucid version that wrote it and the SHA-256 of its action catalogue document (`app-stamp.js`) |
| `actions` | `[{type, version, params}]`: one stamp per action type and parameters version its steps use, sorted by type then version, `params` the SHA-256 of that version's parameter schema as canonical JSON |
| `scope` | What the save or share chose (`Scope`: `trackers`, `branches`, `content`, `results`, `notes`, `chats`, `journal`, `times`, `source`; `share-scope.js`): the file holds exactly that |
| `datasetFingerprint` | The dataset it was saved on, exactly six fields ([Dataset identity](#dataset-identity)) |
| `source` | Where that dataset opens from (`dataset-source-record.js`): `{kind: 'sample', dataset}`, `{kind: 'link', locator, dataset}` (the link in its public form: what it carried to grant access is never recorded), `{kind: 'remote', url, dataset}`, `{kind: 'file', name, dataset}` (a file's or folder's name, never a path or a handle), or `null` (a Python or R viewer) |
| `chrome` | Page chrome kept outside the history, put back after a load (`assertPageChrome`): `{pressed: {<group id>: <button id>}}`, the pressed button of each pressed-button group (the selection tool, `{'highlight-mode': 'highlight-mode-lasso'}`); `{pressed: {}}` from a file written without a page. Kept in every form (history, state, steps). An empty group id, or a button not named by a non-empty id, is refused |
| `requires` | `{obs, var}`: the sorted dataset keys of every observation field and gene the state shown at the save uses: a view is coloured by it, filters on it or holds colours for it, or a highlight group is bound to it (a renamed field by its original key; user-defined fields excluded) |
| `cells` | `{spec, kept, digest}` (`cell-subset-record.js`): the subsample's specification, the kept original positions as LEB128 cell-set bytes in base64, and their cell-set digest ({doc}`19_subsampling_specification`); its `k` equals the fingerprint's `cellCount`. Absent: every cell |
| `genes` | `{cells, order, digest}`: how many of the open cells genes were read for (the first in their order), the cell order they were open in (`{method: 'blocks', seed, blockSize}` when every cell was open, else null), and the cell-set digest of their original positions |
| `activeTracker` | The id of the tracker shown |
| `trackers` | Each the history model's `Tracker` (`app/history/model.js`) whose `steps` lists its steps chunks, `journal` names its journal chunk (or null) and `chats` its conversation chunks |
| `chunks` | Every chunk, in file order |

Every chunk record has exactly `id`, `kind`, `priority` (`eager` or `lazy`),
`codec` (`gzip`, or `none` for a chunk stored as it is), `storedBytes`,
`uncompressedBytes`, and its kind's own keys: `tracker` (steps, journal);
`digest`, `encoding` (array); `tracker`, `conversation` (chat).

### Chunks

| Chunk id | Kind | Priority | Contains |
|---|---|---|---|
| `steps/<n>/<k>` | `steps` | eager | A run of the n-th tracker's steps as JSON, in id order (the tracker lists its steps chunks, read in order as one list): `before` and `after` name domain values by hash; `params` and `outcome` are skeletons |
| `parts/<priority>/<k>` | `parts` | eager for what the head's values name, else lazy | `{<name>: <skeleton>}`: sub-values stored once |
| `array/<digest>` | `array` | eager for what the head's values name, else lazy | One typed array, stored once whatever number of values and steps hold it |
| `values/<priority>/<k>` | `values` | eager: the values the active head shows; lazy: every other value a step or a journal entry refers to | `{<hash>: <skeleton>}` |
| `journal/<n>` | `journal` | lazy | The n-th tracker's journal (Recently deleted and changed); written only when the save includes it |
| `chat/<n>/<m>` | `chat` | lazy | One conversation of the n-th tracker, JSON ([Conversations](#conversations)); written only when the save includes conversations |

Every eager chunk comes before any lazy one, and within each priority the kinds
come in the order of the table (`KIND_ORDER` in `history-file.js`). `steps`,
`parts` and `values` are JSON cut where their text passes 16 MiB
(`CHUNK_TEXT_CHARS`): no single string nears V8's longest (2^29 − 24
characters; one JSON string of every lazy value passed it near 70k steps), and
a reader inflates one piece at a time. Each chunk is stored gzip'd when that
saves at least a sixteenth of its bytes, else as it is (`codec: 'none'`,
`storedBytes` equal to `uncompressedBytes`): skeletons shrink about 5-8×,
float arrays and cell-set bytes barely. A step whose parameters name an array
or a part only the lazy chunks hold is read with them still a skeleton, and
given joined by `readRest` (`params`). Rejected before success: an unknown kind
or key, a missing or duplicate chunk, a chunk out of that order, no eager
`values` chunk (or any `values` chunk in a file of steps only), a steps,
journal or chat chunk no tracker names (or another tracker's), an array stored
twice, a chunk that inflates to another size than declared, a value, part or
array nothing refers to, a step that refers to a value the file lacks, and
bytes after the last chunk.

Every tracker is checked whole by the model (`validateSession`, `model.js`):
one root holding no action, parents and children that agree, each branch one
line from its fork step to its tip, groups that are runs of consecutive steps
holding whole runs of a turn, unique branch, group, tracker and bookmark
names, the head on the active branch's path, every id below its counter, a
paused tracker naming who paused it, a squashed step's carried notes and
bookmarks, each journal entry's own rules (a squash's range, a group edit's
memberships), the action stamps exactly those the steps use, and what the
scope promises.

### Conversations

A `chat` chunk holds one conversation the chat bar kept with the tracker
(`history.setChat`; the tracker's `chats` list `[conversation id,
transcript]`, `model.js`). The transcript is the chat's conversation as
`app/agent/chat/conversation.js` `transcriptOf` writes it: `{schema, id,
createdAt, updatedAt, summary, messages}` (`schema` is `TRANSCRIPT_SCHEMA`,
1; a transcript of another schema is refused when continued), each message
one of

- `{role: 'user', text, turn, at, context?, attached?}` (one per message the
  person sent; `turn` names every step the agent took in answer; `context`
  what Cellucid told the model with it; `attached.selection` the selection
  reference: a page's name, its groups and counts, never cells);
- `{role: 'assistant', turn, provider, agent, model, parts, at, epoch, stop,
  detail?, usage?, otherCalls?}` (one per model call, `epoch` null, without
  the providers' reasoning parts, which are opaque and bound to the
  conversation);
- `{role: 'tool', turn, results: [{id, name, title, content, isError,
  step}]}` (`step` the history step a change made, `{id, label}`, or null;
  each picture a tool returned replaced by a line saying it was there);
- `{role: 'note', turn, kind, text, otherCalls?}` (`stopped`, `error`,
  `refusal`, `limit`, `summary`: what the chat told the person, never sent to
  a model).

No key, provider credential or picture is ever in it. A conversation loaded
from a file continues as a conversation of its own (`conversationFrom`); one
kept mid-turn has its calls with no result answered with what became of them
(`closeInterrupted`). With `chats: false` the chunks are left out, and so is
what the turns say of them (`conversation`, `message`). Python and R read the
chunk as JSON and do not interpret it.

### Arrays

An array's digest is `<TypeName>:<length>:<sha256 of its bytes>`
(`typedArrayDigest`). An array chunk is `raw` (its elements, little-endian,
exactly `length × BYTES_PER_ELEMENT` bytes) or `cell-set`: a non-empty,
strictly increasing `Uint32Array` (a highlight group's cells, cells hidden by
hand) stored as a uvarint bound (its last element + 1) followed by the
cell-set codec's bytes, when that is smaller than its raw bytes. The codec
keeps the smallest of sorted-delta uvarints (encoding 0), a bitset (1) and
runs of gap and length (2) (`codecs/cell-set.js`); any other encoding is
refused. Every decoded array is checked against its digest.

### Values: skeletons and content hashes

A value is stored as a **skeleton** (`value-parts.js`): the value as JSON with
every typed array as `{"$array": "<digest>"}`, every number JSON cannot
hold as `{"$number": "NaN" | "Infinity" | "-Infinity" | "-0"}`, and every
object or list whose skeleton text is at least 512 characters
(`PART_MIN_CHARS`) as `{"$part": "<name>"}`, its skeleton stored once in a
`parts` chunk under the SHA-256 of its text (bottom-up: a part names the parts
inside it). A part is storage only: a value is checked against its own hash
once joined, a part's name is never trusted, a part is joined once and shared
by every value that names it (`partJoiner`), and a part that names itself
(through others) or a value whose text would pass 2^29 − 24 characters
(`MAX_VALUE_TEXT`) is refused before anything is expanded. Its key in a
values chunk is its **content hash** (`value-pool.js` `hashValue`): the
SHA-256, 64 lowercase hex digits, of its canonical text, which is JSON with
object keys sorted, `undefined` properties left out, and typed arrays and
those numbers in the same marker forms. So NaN, ±Infinity and -0 never hash
like the text "NaN" or like 0, and two steps show the same state exactly when
their hashes agree. A value's own object may not use the key `$array` or
`$number` or `$part`, so no skeleton is ambiguous; a value is plain data (an object of
another prototype is refused when it is captured).

The reader rebuilds each value from its skeleton (object keys defined as own
properties, so a `__proto__` key in a crafted file is data, not a prototype),
recomputes its hash and refuses a mismatch. cellucid-python
(`src/cellucid/session.py`) and cellucid-r (`R/session-file.R`, `R/json.R`,
`src/number_text.c`) compute the same canonical text and hash, numbers written
as JavaScript's `Number::toString` writes them.

## Saving

`createSessionBundle({choice})` runs one save at a time, after any load under
way, holding the capture settlement (the gene list read; no field interaction
under way). `choice` is a `share-scope.js` choice: `SAVE_DEFAULTS` (all of the
history, Recently deleted included, conversations left out), `SHARE_DEFAULTS`
(the tracker shown with its main and shown branches and the notes; no
conversation or Recently deleted), `STATE_DEFAULTS` (Current state only: no
notes, conversations, authors or times) or one made from them:

1. the actions in flight finish (`dispatcher.idle()`);
2. what changed while recording was paused becomes one squashed step
   (`history.settlePaused`), and `history.reconcile()` keeps a change made
   outside an action as its own step, so the file holds what is shown;
3. the dataset part (`captureSessionHeader`: fingerprint, `source`, `cells`,
   `genes`, `requires`) is taken with the open cells in use (`useOpenCells`,
   `data/open-cells-work.js`): a growth step or **Apply** waits until it is
   taken, and a dataset replaced, remapped or resized meanwhile fails the save
   (`SESSION_CAPTURE_DATASET_CHANGED`, "Save the session again.");
4. `shareHistory` (`share-scope.js`) makes what the choice keeps of
   `history.snapshot()` and its values, and `writeHistoryFile` writes it with
   the `app` and `actions` stamps (`app-stamp.js`), checking what it wrote
   with the reader's own manifest and model checks.

Writing a file is the catalogue action `history.share` (State › Session ›
**Share…**, `ui/modules/state/share-dialog.js`: History, Current state only
or Steps only; branches, notes, conversations, Recently deleted, results,
times, dataset; What's inside lists what the file will hold; an agent's
`history.share` is asked about as `writes-file`). Its parameters left out,
it writes `SHARE_DEFAULTS` (`shareScopeOf`, the catalogue's one rule).
`history.save` keeps the same choice as a saved state of this browser (the
Load dialog). The Python and R hosts' `requestSessionBundle` pass
`SAVE_DEFAULTS`, so conversations are left out (`chats: false`).

### What a choice keeps (`scope`)

The manifest's `scope` records the choice, and the file holds exactly that
(`validateSession`): the trackers and branches named (a branch brings the
steps it starts from); notes and bookmarks only with `notes`; conversations
only with `chats`; the journal only with `journal`.

- `content: 'history'`: the trackers and branches with every value their
  steps name.
- `content: 'state'` (Current state only): one tracker of one branch holding
  the root (the dataset as opened) and one squashed step "Shared state"
  (`cause: 'share'`, no count, no author) shown, with no turn, group, note,
  journal or conversation. Python, R and this reader apply it as any file.
- `content: 'steps'` (Steps only): no `values` chunk; each step keeps its
  record with `null` in place of every parameter that names cells by
  position (`x-cells`) or holds a value per cell (`x-values`), so it replays
  as unmatched ("Its cells are not in this file."); the file opens as a
  replay plan (`history.replay`).
- `results: 'recompute'`: each analysis result is left out (`result: null`
  beside its `settings`), so its value is stored under its own hash and the
  analysis is run again from its settings (always so with `steps`).
- `times: 'dates'` rounds every time down to its UTC day, `'none'` writes
  every time as null (the file's own time 1970-01-01): the records' times,
  the times values hold (a user field's `createdAt`, an analysis result's
  `computedAt`, `share-scope.js` `VALUE_TIMES`) and a conversation's own.
- `source: 'name'` writes `source: null`. Chunks are compressed from paced slices
(`formats/io/inflate.js` `createSlicedInputStream`), so a large value is never
one long task.

## Loading

**Load** (State › Session › **Load…** › **Open a session file…**,
`session-controls.js`) reads only the manifest first (`readSessionDataset`), waits for a dataset change under way, then opens
the session's own dataset when another (or none) is open: from its `source`
through the same connection controls a person would use; a file of this
computer through the handle this browser kept of it (at once when the browser
allows reading it, else after the person's click on **Reopen <name>**); with
no handle, the restore names the dataset it needs. The dataset is then read
with the session's cells (`cells`: exactly the kept cells; none: every cell)
and genes for the cells the session read them for, and the file is loaded
(`restoreFromBlob`). A saved state of this browser takes the same path
(`loadSavedHistory`; State › Session › **Load…**, the Continue offer, My datasets'
Continue).

**Beside, not instead.** A file chosen with Load is **added beside** the
person's trackers (`history.importTrackers`): a tracker whose id is taken
here is issued another (its conversations, Recently deleted and late
parameters follow it), a name taken is suffixed "(shared)", and each
tracker added is journaled (`load`, "Loaded … from <file>"): the notice
after the load offers **Undo**, which removes them to Recently deleted. A
saved state opened from this browser replaces the history (`history.load`).
A file of steps only opens as a replay. The page chrome a file carries is
planned against this page before anything is replaced (`planPageChrome`),
applied once its head is shown, and what this page lacks of it is said. A
file whose dataset fingerprint does not match can still be **replayed** on
the dataset shown (Load… › Replay a session file here…, `openSessionFile`):
it is read whole and checked as a file, its values are never put back, and
its steps are applied again by name in a new tracker
({doc}`29_action_layer_domains_and_history`, "Replay on another dataset").

A file a Python or R host pushes (`load_session` on the viewer
channel: `viewer.load_session()`, `cellucid_load_session()`) is loaded by
`restoreFromBlob` too, while the viewer shows the host's dataset, and
answered with the tracker then shown and its step shown (`session_loaded`).

The load itself is one owned operation (`_runOwnedRestore`):

1. a newer load cancels the one under way and waits for it to settle;
2. the restore settlement is held, and a progress notification starts;
3. the manifest is checked whole (`openHistoryFile`);
4. the fingerprint is verified against the open dataset
   (`verifySessionFingerprint`); each mismatch has its own sentence
   ([Refusals](#refusals));
5. every key of `requires` is checked (`describeMissingRequiredFields` names
   every missing one at once; genes only when the saved mapping names genes the
   way the open one does, else the domain that holds a gene names it);
6. the eager chunks are read: the steps, the parts and arrays the head's
   values name, the head's values; every tree is checked and the head's values
   must all be there;
7. `history.load` refuses a domain this version does not have (naming it),
   keeps the value shown for a domain the file's root does not hold, and shows
   the active head: only the domains whose value differs are put back, in
   domain order (`cells` first, then `mapping`, then `views` before every
   per-view domain);
8. the lazy chunks are read while the head is shown and new actions wait:
   the other parts and arrays, the other values, the journal, the
   conversations; then the file must end, every value a step refers to must be
   held, and no value, part or array may be left unreferenced;
9. the trackers are added beside the history's (a saved state: they replace
   them), the progress notification ends with
   **Session fully restored**, and `onLoaded` runs: the camera path starts
   playing when it was saved with Autoplay on and has two keyframes
   (`cinematic-camera` `startAutoplay`).

Any failure after step 7 began, and **Cancel**, puts back what was shown before
and keeps the trackers there were; nothing is applied piecemeal. A failure is
said once, on the load's own progress notification (marked reported, so the
session controls do not say it again). The autosave keeps a loaded session as
a new saved state, whole.

A session whose data mapping differs from the open dataset's, but which the
dataset can show and use (other embeddings turned on or as defaults, another
expression matrix the dataset holds), is not refused: the `mapping` domain is
put back right after the cells and before every view, so every view the
session names can be shown; a failure later puts the previous mapping back
with everything else.

### Bounded gzip

A chunk's declared `uncompressedBytes` is both its exact output length and the
allocation ceiling. Native `DecompressionStream` may produce one very large
output chunk at once (WebKit inflates everything written to it), so
`codecs/gzip.js`:

- checks the envelope's constant-size parts before constructing the native
  decompressor: at least a header and a trailer, the gzip magic, method 8
  (DEFLATE), no reserved flag bit, and a trailer ISIZE equal to the declared
  size modulo 2^32;
- feeds the input in slices (`createSlicedInputStream`), so the stream is
  cancelled within one slice once its output passes the declared size, and
  the thread is handed back between slices;
- requires exactly the declared length.

What needs the DEFLATE stream decoded (block structure, the CRCs, the trailer
ISIZE against the real output, further members and trailing bytes) is the
native decoder's job: every supported engine rejects each of them, which
`tests/browser/gzip-native-envelope.spec.mjs` pins per engine. A refused
inflate is "Gzip input is not exactly one well-formed member: …".

### Cancellation outcomes

Two coded aborts exist: a load cancelled by the person (**Cancel** on its
progress, `cancelled('user', …)`) and a load superseded by a newer one
(`cancelled('superseded', …)`). Direct serializer APIs reject with the coded
`AbortError` (its `kind` `user` or `superseded`), so lifecycle owners settle
correctly; the progress notification is dismissed for exactly those, and the
session controls map them to a no-op: nothing is logged or said.
`loadSessionFromFile` resolves `false` for them. An unrelated plain
`AbortError` remains a real failure.

## Starting views

An official catalog dataset may advertise one starting view: its catalog
entry's `session` names the file's path, its SHA-256 and its exact byte length
(`catalog-default-session.js`). It is a format-4 file of **one tracker of one
root**, checked on the manifest before any chunk is read
(`assertPublishedDefaultProfile`):

- one tracker whose root is its head and its active tip (no step, no branch),
  no conversation, every chunk eager;
- no `source` (it ships beside its dataset), no `cells` and no `genes` (it is
  captured on every cell and applies to any subsample of its dataset);
- its root holds neither `cells` (the cells open are the person's or the
  link's choice) nor `layout` (the person's own panels, sidebar and theme)
  (`STARTING_VIEW_OMITTED_DOMAINS`, `assertPublishedDefaultRoot`).

It is put back after the dataset is published, without a step, before the
dataset's first tracker starts (`restorePublishedDefaultState` →
`history.showValues`), so the history starts from it. It runs the same
fingerprint verification and the same `requires` check, so it must carry
`cellOrder` and `dataMapping` and be captured in the dataset's default data
mapping; on a subsample of its dataset the cell count compared is the
dataset's and the cell-order digest is not compared (the catalog pins the file
by SHA-256). A starting view succeeds without a word; when it fails, its
progress notification is dismissed and its caller says **Sample view
unavailable** (the dataset itself is open). The same dataset read again on
other cells keeps its history and puts no starting view back.

Publishing (`trimToPublishedDefault(source, publishedIdentity)`) turns a Save
into a starting view: the active tracker's head becomes the one root (named
"Starting view"), every step and conversation goes, `cells` and `layout` are
dropped, and the fingerprint is bound to the identity the catalog serves the
dataset under (`rebindPublishedDefaultFingerprint`, which refuses, naming the
field, a capture taken against other data). It refuses a capture whose steps
made content — highlights, gene sets, a camera path, cells hidden one by one,
field and category edits (`STARTING_VIEW_CONTENT_DOMAINS` must hold at the head
the value they held at the root) — rather than drop it silently.
`cellucid-datasets` imports these functions (`scripts/publish_session_preset.mjs`)
rather than restating the rule. Do not route a starting view through the
generic **Load**: it is applied only by `restorePublishedDefaultState()`.

## The library file (`.cellucid-library`)

My datasets' **All… › Export…** writes the datasets and saved states the
person chose into one file to take to another browser, and **Import…** adds
them (`app/library/library-file.js`: `exportLibrary`, `importLibrary`,
`readLibraryFile`). It is the session bundle's framing under its own
signature `CELLUCID_LIBRARY/1\n`. Its manifest is `{createdAt, app:
{version}, datasets: [<My datasets record>], histories: [{chunk, name, kind,
dataset: {id, name}, createdAt, updatedAt}], chunks: [{id, kind: 'session',
codec: 'none', storedBytes}]}`; each history chunk is one complete
`.cellucid-session` (format 4) of a saved state, read whole, with every
tracker and branch, results, notes, Recently deleted, conversations and exact
times (`LIBRARY_SCOPE`), stored without a second gzip. The export takes, by
default, every dataset and its saved states this version opens. A saved
state's paused draft cannot be in a session file: the export leaves it out,
keeps the tracker paused, and says so.

A library file is untrusted input. The import checks all of it first (the
manifest, every record by `readRecord`, every history by `readHistoryFile`),
then adds only datasets the list does not hold (`library.add` never replaces
one) and keeps each history as a named saved state "<name> (imported)"
beside the others: nothing kept is replaced, and two imports keep two
copies. A session file given to the importer is refused saying it is one
(`SESSION_FILE_NOT_A_LIBRARY_MESSAGE`), and a library file given to Load
(`LIBRARY_FILE_NOT_A_SESSION_MESSAGE`). Tests:
`tests/library-file.test.mjs`, `tests/browser/my-datasets-library-file.spec.mjs`.

## Dataset identity

`getDatasetFingerprint()` in `session-context.js` builds the record written
into the manifest at save and re-derived from the open dataset at load:

| Field | Derived from | Notes |
|---|---|---|
| `datasetId` | `dataSourceManager.getCurrentIdentityId()` | the dataset's identity id, whatever route opened it (a local folder's `local-user:` key is not used); `null` only when no data-source manager owns the loaded state |
| `name` | `dataSourceManager.getCurrentMetadata().name` | trimmed display name or `null`; informational, never compared |
| `cellCount` | `state.pointCount` | |
| `varCount` | the dataset's own var fields in `state.varData.fields` | user-defined fields are not counted; at load, the count of the var axis the saved mapping names (`describeMappingVarCount`: `raw.X` has its own) |
| `cellOrder.embedding` | `state.getViewLayoutKey('live')`, named under the data mapping | `{key, columns}`: an `obsm` key with its 0-based columns |
| `cellOrder.dimension` | `state.getViewDimensionLevel('live')` | exactly 1, 2, or 3 |
| `cellOrder.digest` | `state.dimensionManager.getCellOrderDigest(<the live layout>)`: the reader worker's sums of that embedding's stored values (`data/cell-order-digest.js`) | 16 lowercase hex characters |
| `dataMapping` | `state.getDataMapping()` | a canonical data mapping (`data/data-mapping.js`): the embeddings turned on and each dimension's default, vector bindings, neighbour graph, gene names, expression matrix, `log1p`, highly variable genes, latent space, and `modality` for a dataset of several modalities |

`getDatasetFingerprint()` is asynchronous (the digest is the reader worker's)
and throws unless `state.positionsArray` is a `Float32Array` of exactly
`cellCount * cellOrder.dimension` values, so a fingerprint is never produced
from a partially published dataset. Once the digest is answered, the capture
re-checks every other field (`datasetFingerprintStillCurrent()`).
`assertDatasetFingerprint` validates a record from untrusted input (exact
six-key set, nullable strings for the id and name, safe integers for the
counts, a valid and already canonical `dataMapping`, a `cellOrder` with exactly
`embedding`, `dimension` in 1–3 and a `/^[0-9a-f]{16}$/` digest);
`datasetFingerprintMatches(a, b)` compares every field but `name`.

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
`data/embedding-payload.js`), in blocks of 4,096 open rows
(`CELL_ORDER_BLOCK_ROWS`); a whole read is summed after its payload is handed
to the page, 262,144 rows a task, so the page's next requests are answered
between slices (anything that needs the sums first — a digest, a growth, an
opening cut in place — takes the rest at once): a growth step adds only the
rows it reads, an opening cut in place keeps the sums of the whole blocks it
keeps and reads the fewer than 4,096 rows of the block it cuts through when
the digest is asked for, and the same cells in another order keep the total.
A digest of an embedding not read in the opening reads its stored values. A
layout held only in memory (a generated dataset, every cell open) is digested
over the coordinates it holds. Saving a session folds nothing on the page.

### Why `cellOrder` exists

A session's values name cells by original position. The three scalars all
survive a row permutation: write a dataset again at the same id from re-sorted
input and `datasetId`, `cellCount` and `varCount` are unchanged, so every
restored position would denote a different cell with nothing to detect it.
Observation and variable names do not help either: they are invariant under a
row permutation, which is the change being caught. Coordinates are what a
re-ordered export permutes: with each cell named by its original position,
another cell at a position changes that position's hash.

`cellOrder.embedding` travels with the digest because a digest is comparable
only within one embedding's columns. A load verifies the digest of the very
embedding the session names, over the cells the session names (opened first),
whichever embedding is on screen and without loading it for display, so a
session saved on a 2D t-SNE opens while a 3D UMAP is shown.

The fingerprint is still not a content hash: the digest covers one embedding's
coordinates, not expression or metadata. Dataset publishers must assign a new
identity when scientific content changes; what `cellOrder` adds is that a
re-ordered republication under an unchanged identity is detected instead of
trusted.

### Refusals

`verifySessionFingerprint(ctx, saved)` returns `null` when the saved record can
be loaded here, otherwise the sentence the session controls show verbatim.
Every check runs before any value is put back, and each cause has its own
sentence: telling someone their data was re-ordered when the session only named
another embedding would make them distrust a sound dataset.

| Cause | Detected by | What the person is told |
|---|---|---|
| A session of another format version (format 3 and earlier, or a newer one) | the signature, `bundle/reader.js` | `otherSessionFormatMessage(version)` ([Container](#container)) |
| A library file of My datasets (`.cellucid-library`) | the signature, `bundle/reader.js` | `LIBRARY_FILE_NOT_A_SESSION_MESSAGE`: what it is, and where it is imported |
| A file without a session signature | the signature, `bundle/reader.js` | `NOT_A_SESSION_FILE_MESSAGE` |
| A damaged file (a manifest, chunk, array, value or tree that breaks the format) | `history-file.js` | "This session file is damaged: …", naming what |
| No dataset open | `verifySessionFingerprint()` | to open the dataset the session was saved on (named), then load the session again |
| A different dataset | the counts differ, or the counts agree but the id differs | both cell and gene counts, or both datasets by name, and to open the one the session was saved on |
| An enabled embedding this dataset cannot show; a data-mapping group (neighbour graph, modality, gene names, expression, `log1p`, highly variable genes, latent space, vector fields) whose saved choice this dataset cannot use, or which this dataset cannot change and is shown with another value | `describeDataMappingMismatch()` | the embedding or the group is named, with the reason it cannot be used, or both values |
| `cellOrder.embedding` names no layout of this dataset | `verifySessionFingerprint()` | the embedding and its dimension are named |
| `cellOrder.digest` differs, everything else matches | `verifySessionFingerprint()` | the same name and counts but other cell coordinates: the cells are listed in another order, or were exported again from a re-computed embedding, so every saved selection could mark the wrong cells |
| A field or gene in `requires` is missing | `describeMissingRequiredFields()` | every missing key at once |
| A state domain this version does not have | `history.load` | the domain is named, with "Open it in the version of Cellucid it was saved with." |

None of these messages carries identifiers or internal vocabulary. No
dataset-agnostic part of a refused file (its layout, say) is salvaged: the
whole load rejects and rolls back.

## Adding or changing persistence

A new persisted feature is a state domain's value, or part of one
({doc}`29_action_layer_domains_and_history`): the session file stores it as it
is, with no encoder, chunk or profile of its own. Name its fields, genes,
categories and views by name and its cells by original position, keep typed
arrays immutable, restore it exactly and idempotently, and add capture,
restore, undo and file round-trip coverage
(`tests/state-domains-contract.test.mjs`,
`tests/browser/state-domains-roundtrip.spec.mjs`,
`tests/session-history-file.test.mjs`).

A change to the file itself (the container, the manifest, a chunk kind, the
skeleton or the hash) changes the only current format. Do not introduce a
reader for an older shape. Update every current producer, consumer, fixture
and reference together:

- `app/session/README.md` and this page (`tests/session-docs-current-contract.test.mjs`
  pins the README's wording);
- the shared golden history: `scripts/generate-history-golden.mjs` makes it
  with this app's model and writer (every kind of record: three trackers,
  branches, turns, groups, notes, a bookmark, squashed steps (a paused period
  and a squash in mid-branch with a fork after it), shared parts, the journal,
  a stale result) and writes `tests/fixtures/history-golden.cellucid-session` and
  `history-golden.expected.json` (what reading it and applying it must give,
  the apply plans row for row), byte for byte the same as cellucid-python
  `tests/fixtures/history_golden/` and cellucid-r
  `tests/testthat/fixtures/history_golden/` (also
  `inst/extdata/history_golden.cellucid-session`); its digests are pinned in
  `tests/history-golden-fixture.test.mjs`, cellucid-python
  `tests/test_sessions.py` and cellucid-r `tests/testthat/test-sessions.R`;
- the packages that read a session and apply one of its steps
  (cellucid-python `src/cellucid/session.py`, `history.py` and
  `anndata_session/`, cellucid-r `R/session-file.R`, `R/history.R` and
  `R/session-apply.R`), which read every record and the `cells`,
  `field-overlays`, `highlights` and `gene-sets` values;
- the starting views `cellucid-datasets` publishes, and the user-facing session
  chapter ({doc}`../l_sessions_sharing/index`).

## What a session does not keep

- the dataset itself (a session names where it opens from and assumes it is
  loaded): points, observation tables, expression matrices;
- connection state beyond `source`, and nothing that granted access to a link;
- per-cell gene values an analysis read (they are read from the dataset again);
- the antialiasing preference (a device capability);
- community annotation state (votes, comments, moderation, GitHub sync and
  sign-in);
- benchmarking controls, dataset selection and connection controls;
- notifications, hover, focus, a pointer lock, a selection not confirmed, or
  in-flight requests; DOM, worker, network and WebGL runtime objects;
- the history settings (`history-settings.js`: per browser, never part of a
  session).

## Tests

`tests/session-history-file.test.mjs` (the format: round trip, eager and lazy,
arrays once, conversations only when included, every refusal),
`tests/session-serializer-contract.test.mjs` (save and load, settlements,
cancel and supersession, `requires`, starting views),
`tests/session-starting-view-contract.test.mjs`,
`tests/session-cell-order-identity-contract.test.mjs`,
`tests/session-cell-order-digest-contract.test.mjs`,
`tests/session-gzip-codec.test.mjs`, `tests/history-golden-fixture.test.mjs`,
`tests/session-share-scope.test.mjs` (every choice, the state-only form, the
steps-only form), `tests/session-selection-tool-contract.test.mjs` (page
chrome), `tests/library-file.test.mjs` (the library file's framing),
`tests/session-docs-current-contract.test.mjs`,
`tests/browser/history-share.spec.mjs` (Share and Load on a page),
`tests/browser/history-session-file.spec.mjs` (a saved history loads in
another page at its active step; a format-3 file is refused naming its
version, and nothing changes) and `tests/browser/gzip-native-envelope.spec.mjs`.

Next: {doc}`11_analysis_architecture`.
