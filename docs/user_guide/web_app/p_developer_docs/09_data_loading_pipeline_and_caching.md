# Data loading pipeline and caching

This page documents the **developer-facing data loading pipeline** for the Cellucid web app:

- how datasets are discovered and selected
- how the app derives a dataset **base URL**
- which files are fetched (and in what order)
- where caching happens (HTTP cache, in-memory caches, byte-bounded field residency, IndexedDB caches)
- how remote server / links / file picker / notebook and R viewer modes differ

It is the reference: what the code does now. Why it is built so, and the
measurements behind it, are in the decision records:
{doc}`23_decisions_partial_reading` (what is read of a large file, and its
full description of shares, gene samples and kept indexes),
{doc}`22_decisions_link_model_and_request_planning` (the link, requests to a
web host, bytes kept across readers and visits),
{doc}`21_decisions_decoders_and_codecs` (decoding) and
{doc}`27_decisions_memory_model` (what may be held). The module map of the
read path is {doc}`20_read_path_architecture`.

If you are looking for a user-facing guide, start with:
- {doc}`../b_data_loading/index`

## At a glance

**Audience**
- Computational users: read “Sources” + “Debugging checklist”.
- Developers: read fully before changing data formats or loader behavior.

**Time**
- 30–60 minutes

**Prerequisites**
- {doc}`05_app_architecture_overview`

---

## Mental model: “base URL” + loaders

The frontend separates data loading into two layers:

1) **Source selection**
   - Choose *where* the dataset comes from (the sample catalog, a local file or folder, a remote server, a link, a Python/R viewer).
   - Output: a canonical `{ sourceType, datasetId, baseUrl, metadata }`.

2) **File loading**
   - Given `baseUrl` (`anndata://<dataset id>/`), load the dataset's
     identity, field inventories, embeddings and graph.
   - Output: typed arrays + field loaders attached to `DataState`.

Every dataset is read by the one AnnData reader through one data source
(`anndata-file-source.js`, `AnnDataFileSource`), whatever file it arrived as,
so every base URL is an `anndata://` URL; the format is a fact of the open
dataset (`source.format`, `describeStorage()`), not of its address.
`anndata-provider.js` resolves a route under it to the AnnData
adapter of the dataset that owns it (`base-anndata-adapter.js`), which answers
it from the file's `uns['cellucid']`, `obs`, `var`, `obsm` and `obsp`. Each
route is named for the AnnData element it reads:

| Route | What it answers |
| --- | --- |
| the base itself | the dataset identity (`loadDatasetIdentity`) |
| `obs` | the obs field list (`loadObsFields`) |
| `var`, `var/{g}` | the gene list the file opens with, and gene `g` (`loadGeneList`); read after the first view |
| `obsp` | the neighbor graph the file opens with (`loadConnectivitySummary`, `loadEdges`) |
| `candidates` | the data-mapping candidate catalog (`loadMappingCandidates`) |
| `obsm/{i}/{cols}` | `obsm` entry `i` read at some columns (an embedding) |
| `vectors/{v}/{vcols}/{e}/{ecols}` | a vector field drawn on an embedding |
| `obsp/{p}` | another neighbor graph |
| `var/{axis}/{names}`, `…/{g}` | another gene list: an expression matrix named by a name source |
| `outliers/{l}/{f}` | the outlier quantiles of categorical field `f` in latent space `l` |

Routes name entries by their position in the candidate catalog
(`data/data-mapping-candidates.js`), never by key. A sample catalog or a
linked catalog is a `datasets.json` that lists AnnData files; everything else
about a dataset comes from the file itself.

- UI code should never build a data path itself; it should always flow
  through the dataset's **base URL**.
- Loader code should never assume where the bytes are (“GitHub raw URL…”):
  that is the source's and the provider's job.

```
  local-demo   local-user   remote      link      jupyter
      │            │          │           │           │
      └────────────┴──────┬───┴───────────┴───────────┘
                          ▼
              DataSourceManager                       ← source selection
              data/data-source-manager.js
                          │
                          │ { sourceType, datasetId, baseUrl, metadata }
                          ▼
              data-loaders.js                          ← file loading
                          │
   ┌──────────────────────┼──────────────────────┐
   ▼                      ▼                      ▼
identity               obs / var            embeddings / points
(uns['cellucid'],      inventories          graph (obsp, optional)
 obsm, obsp)
                          │
                          ▼
                     DataState                        ← typed arrays + field residency
```

Everything below `DataSourceManager` sees only a base URL and relative paths.
That is the whole point of the split: a picked file, a linked file and a
served store are all read by the same AnnData reader, and no loader has to
know which.

Code pointers:
- Source selection: `cellucid/assets/js/data/data-source-manager.js`
- Loaders: `cellucid/assets/js/data/data-loaders.js`

---

## Startup path (what happens at boot)

At startup, `cellucid/assets/js/app/main.js` (its start-up sequence; the dataset's
publication and change are `app/dataset/publication.js` and `reload.js`):

1) Registers core sources:
   - `local-user` (file picker)
   - `remote` (cellucid-python server)
   - `jupyter` (only on a viewer page a host program opened:
     `?viewer=true&viewerId=…&viewerToken=…`, `getViewerHostConfig()`)
2) Initializes `DataSourceManager`, which also registers:
   - `local-demo` (the sample catalog: the **Sample dataset** list)
   - `link` (data opened from a link through a provider of
     `data/sources/registry.js`, or any web address)
3) Applies URL overrides (if present):
   - `?remote=...` connects remote server
   - `?data=...` opens a link (a provider's locator such as
     `gh:owner/repo/path`, `cellxgene:<id>`, `geo:GSE…`, `zenodo:<record>`,
     `gdrive:<id>`, or a web address; see
     {doc}`04_configuration_env_vars_and_feature_flags`)
   - exactly one canonical `?sampleCatalog=...` may replace the configured
     sample catalog for that launch
   - `?dataset=...&source=...` selects a specific dataset from a registered source
4) Resolves the startup selection (source, dataset id, load method) without
   loading it, publishes the empty runtime, and builds the interface, so every
   data control is live before any dataset is fetched.

Then the startup selection loads as the first ordinary dataset change (see
{ref}`dataset-change-transactions`), and its loader phase begins:
- open the AnnData file in the reader worker (the structure the first
  view reads, `ANNDATA_OPENING_PATHS`, and `uns['cellucid']`), and read the
  dataset's identity
- load the obs inventory → attach obs field loader
- optionally load the var inventory → attach var field loader
- load the positions of the embedding the dataset opens on (the default of
  its opening dimension, `dimensionManager.getDefaultLayoutKey()`), read
  beside the obs field lists rather than after them; every other enabled
  embedding, including the other dimensions `uns['cellucid']['defaults']`
  names, is read (and its values checked) when a view first shows it (what
  that saves: {doc}`23_decisions_partial_reading`, Measurements)
- optionally describe the neighbor graph (edges load when shown)
- note the latent space the file stores outlier quantiles for
  (`uns['cellucid']['outlierQuantiles']`); no latent space and no quantile is
  read at open

Opening the startup dataset runs as the `dataset` startup step
(`app/startup-phase.js`). Every request it makes bounds itself (see
{ref}`web-host-requests`), so the step's 60 s deadline measures silence, not
duration: it is restarted whenever the transfer log (`data/transfer-log.js`) counts a request or
hears bytes of one still under way (`arriving`), and it fires only when
nothing at all has arrived for 60 s. Its error names the reads still open
(the read activity log) and the request being asked again, with its URL:
`Startup stalled while opening the dataset (its metadata, field lists and the
embedding its first view draws): nothing arrived from the host of dataset
"suo" for 60 s, while reading obsm['X_umap_3d'] columns 0, 1, 2 (last
failure: https://…/obsm/X_umap_3d/3.0: no answer for 35.6 s. Asked 5 times so
far.)`.

Every opening, at startup or later, is measured by the same rule without
being stopped by it (`data/opening-progress.js`). An opening whose step,
layers, reads (requests, bytes, chunks decoded) and transfers have not moved
for 60 s (`STARTUP_STEP_DEADLINE_MS`) says what it waits on in its
notification, or the panel over the empty viewer, and in the console:
`[cellucid:opening] Opening a.h5ad stalled while reading the file's structure:
nothing has changed for 1 min 0 s, while reading the structure of a.h5ad (0 B
of 512 KB arrived).` A wait that is not a read is named by the step its layer
set: the dataset change names staging the first view, stopping the previous
dataset's saved-state restore, waiting for the GPU, allocating its memory and
drawing the first view (`reloadActiveDatasetInPlace`); the file and
connection layers name updating the dataset list. A wait on the person
(`utils/person-wait.js`) is not silence, and while a startup step's deadline
runs (`onStartupStepDeadline` in `app/startup-phase.js`) that deadline is the
one measure: it ends startup and names the stall, so the opening's measure
stands down and the stall is reported once. The opening goes on, and Cancel
stays.

Closing the application (`window._cellucidDispose`, `pagehide`) stops the
dataset change in flight first (`retire()` of the coordinator below) and the
openings the dataset controls began (a file's reader, a connection's
candidate), so its teardown waits only for work told to stop, never for a
read to answer on its own.

A stream reader is read to its end, or canceled, and dropped: no browser code
calls `releaseLock()` (`tests/browser-stream-reader-release-contract.test.mjs`),
because a release can hang WebKit's garbage collector for the rest of the tab
(the decision: {doc}`21_decisions_decoders_and_codecs`).

(dataset-change-transactions)=
## Dataset changes: the newest choice wins

Every change of the displayed dataset is one transaction of the dataset-change
coordinator (`createLatestDatasetReloadCoordinator` in
`app/dataset-reload-outcome.js`): the startup selection, a sample, **None**, a
local file or folder, a connection that opens its only dataset, and
the synthetic benchmark. The rules:

- **A transaction begins when the user acts**, before any of that choice's own
  work — before a file is read, before a sample's first request. Beginning one
  aborts the previous transaction's signal, so the load in flight stops
  fetching (its staged `DimensionManager` is retired, which aborts its
  requests; a local candidate is disconnected, which aborts a file mid-read) and can never publish over the newer choice.
- **Ownership is checked against the newest transaction of all**, not only
  the newest of one kind. The dataset controls own their presentation (the
  select shows the pending choice, or an *Opening …* / *Loading …* stand-in
  for one without a catalog option) only while their intent's transaction
  `isLatest()`; an intent overtaken by the benchmark puts the select back
  itself.
- **A dataset is named as the Dataset box names it**, never by its id (a
  link's is its encoded address, `url:https%3A%2F%2F…`; a provider's names
  its record and file): the opening's label (*Loading atlas…*), the read
  again with other cells (*Reading atlas again…*), a failure (*"atlas" could
  not be opened*), the opening progress's title and the step a stalled
  opening names, and a session's refusal of another dataset. The dataset
  controls learn each name from every listing of every source (those the
  picker does not offer too), from the listing startup chose from
  (`selectStartupDataset`'s `listed`), and from every dataset shown
  (`datasetDisplayName`); the activation carries it (`name`). One no listing
  named yet is *the dataset*, and the opening progress takes its name once
  its source lists it.
- **A connection may or may not change the dataset**, which is known only once
  the server answers. Every user action that decides the dataset advances
  `actionSequence`: each transaction begun, and each session file applied
  (`recordAction()`, which begins no transaction and aborts nothing). The
  address bar's **Open** reads `actionSequence`; a one-dataset activation then calls
  `beginIfNoActionSince(sequence)`, which supersedes loads that were running at
  the click but never a choice made or a session applied after it. A
  connection that lost that race is kept as a catalog. The startup selection
  begins its transaction with `beginStartup()`: it supersedes like any other
  but is no action of the person's, so a link typed before it began still
  replaces it when its server answers.
- **The shown dataset growing does not supersede a newer choice.** A staged
  selection is refused at commit only when another change of *which* dataset
  is shown was committed meanwhile (a selection, its rollback, **None**:
  `DataSourceManager`'s committed-change count). Growth in place
  (`adoptResizedMetadata`, every step of a progressive opening) moves the
  selection revision, which caches keyed to the open cells read, but never
  that count: a link opened while a large sample still grows opens.
- **A session file** chosen during a load waits on `whenSettled()` (the newest
  transaction has called `finish()`) and is then restored onto whatever
  dataset that produced.
- **A lost remote connection** clears only its own dataset, and only once
  `settled` is true: a newer choice still loading replaces the lost dataset
  itself, and recovery clears it only if that choice failed and the lost
  dataset is still active. The manager's `getCurrentBaseUrl()` is the base URL
  the selection committed, never re-derived from a transport that may have
  closed.
- **A drop onto the viewer** (`ui/modules/dataset-drop.js`) is classified
  inside the `drop` event — one file whose name a format of
  `data/input-formats.js` names, or one directory entry from
  `webkitGetAsEntry()`, listed with each file's `webkitRelativePath` set —
  and handed to the same local reader and intent as the file inputs. A
  folder is listed by `data/local-folder.js` (`readFolderEntry`; a folder's
  handle by `readFolderHandle`): each `entry.file()` and `readEntries` batch
  is a round trip to the browser's file system, so they run concurrently,
  as many in flight as pays (`ConcurrencySearch`: from one, doubled while
  the calls answered per millisecond rise by at least half, then the last
  level that paid, the rule `formats/io/link.js` uses for a multiplexed
  link), every directory at once and one directory's batches in turn, the
  files returned sorted by path; the intent's signal and **Cancel** stop the
  listing, and an entry that cannot be read fails it naming the entry. Where
  the browser has `getAsFileSystemHandle()` (Chromium), the item's handle is
  asked for in the same event and kept once the dataset opened
  (`app/library/local-file-handles.js`). The
  classic script `ui/core/file-drop-guard.js` cancels every file drag at the
  document from the first paint, so a drop never navigates, and offers it as
  `none` until the drop target sets `data-dataset-drop="ready"` on the root.
- **A published dataset's own continuations** (UI synchronization, the
  published sample view) belong to its publication, not to the request: a
  newer request that fails must not interrupt the dataset still on screen.
- **Startup** loads its selection through `ui.loadStartupDataset`. Its network
  awaits are deadline-bound startup steps; if it is superseded, startup
  publishes `ready` without it; only a failure nothing superseded becomes the
  terminal startup failure.
- **Closing the application** retires the coordinator (`retire()`): the
  transaction in flight is aborted as a newer choice would abort it, and no
  transaction begins afterwards.

A superseded load is silent: it dismisses its own progress, announces nothing,
and leaves the outcome to the choice that replaced it.

---

## Sources (where datasets can come from)

All sources implement a similar interface:
- `isAvailable()`
- `listDatasets()`
- `getMetadata(datasetId)`
- `getBaseUrl(datasetId)`

### `local-demo` (the sample catalog)

Code:
- `cellucid/assets/js/data/local-demo-source.js`
- `cellucid/assets/js/data/sources/anndata-catalog.js` (the `datasets.json`
  reader shared with catalog links)

Behavior:
- Reads `datasets.json` (schema 2) in the sample catalog folder
  (`DATA_CONFIG.SAMPLE_CATALOG_URL`) and lists its datasets. Listing reads each
  file's metadata and `uns` only (`describe` in the reader worker), never its
  data; a dataset's file is opened when it is chosen.
- Listing waits for any dataset of the catalog being opened, and lists that
  dataset from its opening (`_listedFromOpen`): at startup the default
  sample's store is read once, and the other files' metadata is read only
  after the dataset on screen has what its first view needs, so the listing
  never competes with it.
- `DATA_CONFIG.SAMPLE_CATALOG_URL` resolves in this order:
  1) exactly one `?sampleCatalog=...` query field
  2) the one `<meta name="cellucid-sample-catalog" content="...">` in
     `index.html` (`https://theislab.github.io/cellucid-datasets/anndata/` in
     production)
- Both must be one canonical absolute HTTP(S) folder URL with a trailing `/`
  and no query or fragment; anything else is a startup configuration error.
- A dataset whose catalog entry names a `session` gets that starting view
  (a format-4 file of one tracker of one root, pinned by its byte length and
  SHA-256) put back after it is published, without a step, before its history
  starts (`app/session/catalog-default-session.js`, `published-default.js`;
  {doc}`10_sessions_persistence_and_serialization`, and
  {doc}`../l_sessions_sharing/04_official_sample_states`).

Common failure mode:
- `datasets.json` missing or invalid → the sample list is “not available”,
  and the message names what is wrong (for a catalog without `schema`:
  `Invalid datasets.json: the catalog must contain exactly: datasets,
  default, schema (missing: schema).`).

Cross-origin transport contract:
- `datasets.json` and every starting view are fetched by
  `fetchCatalogFile` and revalidated with the host (`CATALOG_FILE_CACHE`);
  each AnnData file is read in place with HTTP range requests.
- The catalog host must permit the viewer origin with CORS.
- Credentials are omitted; an HTTP, CORS, URL, or payload failure terminates
  that dataset load.

### `link` (data opened from a link)

Code:
- `cellucid/assets/js/data/link-source.js` (the source)
- `cellucid/assets/js/data/sources/registry.js` (providers and `?data=` locators)
- `cellucid/assets/js/data/sources/providers/*.js` (one module per provider:
  `github`, `cellxgene`, `geo`, `zenodo`, `huggingface`, `figshare`,
  `dropbox`, `gdrive`, `onedrive`, `http`)
- `cellucid/assets/js/data/sources/auth/` and `data/sources/access.js` (the
  Google Drive and OneDrive sign-ins and how the reader worker is handed a
  read address; see {doc}`04_configuration_env_vars_and_feature_flags`)
- `cellucid/assets/js/data/sources/anndata-catalog.js` (a `datasets.json`
  catalog at a folder URL)
- `cellucid/assets/js/data/sources/remote-file.js` (one file on a web host:
  its format, exact size and identity)
- `cellucid/assets/js/formats/` and `cellucid/assets/js/data/anndata/` (reading
  AnnData files in place; see [Reading AnnData directly](#reading-anndata-directly))

Behavior:
- The user pastes a link or id, or the page has `?data=<locator>`.
- The registry picks the provider that recognizes it. A provider only says
  where the bytes are: a dataset catalog folder (a folder URL holding
  `datasets.json`; GitHub resolves to
  `https://raw.githubusercontent.com/owner/repo/<branch>/path/`), one file of
  a format of `data/input-formats.js` with its exact size, the files that
  make one dataset together (a GEO sample's MatrixMarket triplet), or a Zarr
  store folder (read through its consolidated metadata).
- A catalog then behaves like `local-demo` (the same `anndata-catalog.js`
  reader; each listed file opened when chosen, its starting view applied).
  A file or Zarr folder is opened by the AnnData
  reader worker in place, with parallel HTTP range requests, and served
  through `anndata://` like a local file (see
  [Reading AnnData directly](#reading-anndata-directly)). How it is read
  follows the link the reader measures while opening: over a link as fast as
  a disk (the same machine or network) a cell-major (CSR) matrix is indexed
  by gene at once, in the background, as a local file is; over a slower one
  it is indexed when the first gene is shown, so opening downloads only what
  the first view draws, and a gene's values are fetched as its entries need
  them.
- Requests carry only a single `Range: bytes=a-b` header, which needs no CORS
  preflight, and bypass the browser's HTTP cache (`cache: 'no-store'`, not a
  header, so no preflight either), because Chromium lets one request of a URL
  through its cache at a time. File sizes come from `HEAD` or the provider,
  because many hosts do not expose `Content-Range`.
- An `.h5ad` opens one round trip after its link is followed: its size's
  `HEAD` goes out with a read of its first 512 KiB, and the reader starts from
  those bytes, which a Cellucid-ready `.h5ad` lays out to hold everything
  opening reads (the specification's "Where things sit in the file"). How the
  probe and the head are made and what each saves:
  {doc}`22_decisions_link_model_and_request_planning`.

Common failure modes:
- `datasets.json` not found (wrong path)
- the host does not allow cross-origin reads (CORS), or ignores range requests
- a CELLxGENE dataset newer than the latest Census release
- GitHub raw blocked by corporate environments

### `local-user` (browser file picker)

Code:
- `cellucid/assets/js/data/local-user-source.js`
- `cellucid/assets/js/app/ui/modules/open-data-menu.js` (the **Open…** menu)
- `cellucid/assets/js/app/ui/modules/dataset-drop.js` (drag and drop)
- `cellucid/assets/js/data/local-folder.js` (a dropped or chosen folder's files, listed concurrently)
- `cellucid/assets/js/app/library/local-file-handles.js` (handles kept to reopen, Chromium)

Behavior:
- Lets the user pick, or drop on the viewer:
  - one file with **Open… › File…** (one file input,
    `#user-data-file-input`, accepting the single-file endings of the format
    registry, `data/input-formats.js` — `.h5ad,.zip,.loom,.h5,.h5mu,.rds.gz,`
    `.rds,.rdata.gz,.rda.gz,.rdata,.rda,.h5seurat` — set from
    `filePickerAccept()`; `localFileLayout()` takes it by its name, and
    anything else is refused naming the file)
  - one folder with **Open… › Folder…** (a directory input,
    `#user-data-folder-input`; every file below it, each with its
    `webkitRelativePath`): a Zarr store, a Cell Ranger `outs/` folder or a
    MatrixMarket folder, told apart by the registry's folder signatures
  - where the page is given handles (`localHandlesAvailable`: Chromium's
    `showOpenFilePicker` and `showDirectoryPicker`, a page not in a frame,
    IndexedDB), **File…** and **Folder…** call those pickers instead of the
    inputs (the same endings, `filePickerAccept()`); the handle, or a dropped
    item's (`getAsFileSystemHandle`), is kept once the dataset opened, in
    IndexedDB `cellucid-local-files` under the dataset id (`keepLocalFile`:
    the handle, name, kind, a file's size and `lastModified` or a folder's
    file count and the SHA-256 of every file's path, size and
    `lastModified`, a line each in path order, when opened). The tab notes the local dataset it shows in
    `sessionStorage` (`rememberTabDataset` on every dataset change); after a
    reload `dataset-connections.js` offers *Reopen <name>* for it, and
    `session-controls.js` reopens a `{kind: "file"}` session source from its
    handle (at once when `queryPermission` is `granted`, else on *Reopen
    <name>*, then applies the session); My datasets opens it again with one
    click (`ui/modules/my-datasets/`). The reopen (`reopenLocalFile`) runs
    in the click: `readKeptLocalFile` asks `requestPermission({mode:
    'read'})` unless granted, reads the file or lists the folder, compares
    what it held, and fails explained (`details.reopen`: `refused`, `gone`,
    `changed`) saying to choose it with **Open…**; a gone or changed handle
    is forgotten
- Each is opened as a candidate by the AnnData reader worker and adopted only
  when it opened completely; it is served through `anndata://` URLs. What it
  is is decided by the reader from its content (a file's first bytes), so a
  misnamed file opens as what it is. Its dataset id is
  `uns['cellucid']['id']`, else `<format>_<name>_<structure>`
  (`h5ad_<file name>_<structure>`, `zarr_<store name>_<structure>`:
  `defaultDatasetIdOf`, the structure being eight hexadecimal digits of the
  file's cell count, gene count and gene names, so files of one name holding
  other data are other datasets; the id names the dataset, never its
  version, which kept stores file beside it); its name is
  `uns['cellucid']['name']` (or `uns['title']`), else the file name without
  its format's ending.

Common failure modes:
- denied or policy-blocked local-file permissions
- a folder no format's signature files mark (a Zarr store's `.zgroup`,
  `.zarray` or `zarr.json`, a Cell Ranger `filtered_feature_bc_matrix.h5`, a
  `matrix.mtx[.gz]`): refused, naming what can be opened; in Safari, whose
  folder dialog leaves out hidden files, a Zarr v2 store is refused with the
  advice to drop the folder or open its `.zarr.zip`
  (`describeHiddenMetadataProblem`)
- a malformed `uns['cellucid']` or non-finite scientific payloads
- invalid/multi-root Zarr ZIPs

#### Browser memory (no fixed ceilings)

No reader refuses data for its size. Every array is declared exactly by the
file's own array metadata (shape × dtype, one chunk at a time), checked
against the transfer and the compression or ZIP envelope, and an array the
browser cannot allocate is refused as **too large** naming the bytes
(`allocateTypedArray` in the readers; an inflate whose output cannot be
allocated is the same refusal). The connectivity graph has no
working-set ceiling either (see
{doc}`07_rendering_pipeline_webgl_and_performance_notes`).

Two real limits remain, both the browser's own: a JSON document (a catalog, a
session manifest) longer than the engine's longest string cannot be parsed
(the message names the file), and one typed array cannot exceed the engine's buffer maximum, which is why the
gene-column index is held in segments (`data/sparse-column-index.js`). Memory
budgets (`data/memory-budget.js`) only decide what is kept for reuse; they
never refuse. If you add a reader, follow the same pattern: bound by the size
the data declares, and name the bytes when an allocation fails.

### `remote` (cellucid-python server)

Code:
- `cellucid/assets/js/data/remote-source.js`

High-level behavior:
- Connects to a server URL (e.g. `http://127.0.0.1:8765`).
- Validates health and server info:
  - `GET /_cellucid/health`
  - `GET /_cellucid/info`
- Enumerates datasets (`GET /_cellucid/datasets`).
- Reads the dataset's file or store from the server with HTTP range requests.
- Opens no WebSocket and polls nothing: health and info are asked when the
  page connects. The viewer channel is the `jupyter` source's alone.

Why remote is important:
- A Cellucid server (`cellucid serve` on any format of
  `data/input-formats.js`, `serve_anndata()`, `show_anndata()`, R's
  `cellucid_serve()` / `cellucid_view()`) computes nothing. Its catalog entry
  (`_cellucid/datasets`) carries `anndata: {file: {url, kind, layout, size?},
  defaults}`: the store's route under `_cellucid/anndata/` (named by the
  dataset id), and the launch defaults the user set. `served-anndata.js` opens
  that route with the same AnnData reader worker a linked file uses (see
  [Reading AnnData directly](#reading-anndata-directly)), under the server's
  dataset id and name, and passes the defaults to the adapter, where they
  override the file's own `uns['cellucid']['defaults']`.
- The server answers the store's bytes with exact `Range` / `206` /
  `Content-Range` / `HEAD` / `ETag` responses over HTTP/1.1 keep-alive (a
  thread per connection, `TCP_NODELAY`, the body sent by the kernel with
  `sendfile`), and a Zarr store folder key by key (404 for an absent key,
  consolidated metadata generated in memory when the store has none). Store
  bytes carry `Cache-Control: no-cache` with a strong `ETag` (size,
  modification time, status-change time, inode, read from the open file — on
  Windows NTFS's change time and file index, read from the file's handle; R
  reads the same from the open file in C), which `/_cellucid/datasets` also states beside the
  file's size, read from the file at each request: the version of the file's
  identity (`served-anndata.js`), so its byte cache, a kept gene index and a
  converted store last across visits and a rewritten file is a new version,
  and the byte source refuses an answer with another `ETag` than the one
  listed. The browser keeps what it read and revalidates
  it on a reload (`304 Not Modified`, nothing sent again), and a range asked
  with `If-Range` naming an older version is answered with the whole current
  file, so kept bytes of a changed file are never mixed with new ones.
- The browser measures the link to the server like any other (see
  [Reading AnnData directly](#reading-anndata-directly)): `cellucid serve` on
  the same machine reads like a local file, the same server behind a remote
  proxy like a file on the internet, and a notebook widget at the speed of
  the notebook's own connection.

Common failure modes:
- mixed-content blocking (https app trying to fetch non-local http server)
- missing CORS headers
- server schema/version mismatch

### `jupyter` (a viewer a host program opened)

Code:
- `cellucid/assets/js/data/viewer-host-source.js` (the source)
- `cellucid/assets/js/data/viewer-channel.js` and `data/viewer-protocol.js`
  (the channel and its rules); `data/viewer-embed.js` (the notebook widget's
  transport)
- `cellucid/assets/js/app/viewer-host/results.js` (what the person made,
  published to the host as `result` events) and
  `app/registries/derived-results.js` (the gene sets and scores the viewer
  holds besides its fields: pushed sets, the shown analyses' sets and scores)
- Protocol architecture: {doc}`../../python_package/e_jupyter_hooks/05_architecture_viewer_channel`

Behavior:
- A browser-tab viewer (Python `AnnDataViewer`, `serve_anndata()`, R
  `cellucid_view()` / `cellucid_serve()` outside a notebook) is a page its
  host serves, opened with `?viewer=true&viewerId=…&viewerToken=…`. Its
  `/_cellucid/info` names the program (`runtime`), which the Dataset box
  shows. Commands from the host and events from the viewer travel on one
  WebSocket to the same server (`/_cellucid/ws`), authenticated with the
  viewer id and token, with no polling. The host's pings keep an idle socket
  open; a dropped socket is reopened once, and one health request then tells
  a stopped host from a proxy that does not pass WebSockets.
- A notebook widget (Python `show_anndata()` / `cellucid.Widget`, R
  `cellucid_view()` in IRkernel) frames the app from its web source with
  `?viewer=embedded`. The widget hands the page its viewer id and token in a
  message posted to the app's origin only, over a private `MessagePort`
  (accepted only from the frame's embedder), and every request and the
  viewer channel then run over that port and the kernel's comm: no server,
  port or proxy. Above the transport -- the documents, the store reads, the
  protocol -- everything is the same as in a browser tab.
- Either way the dataset is read like a server's (`served-anndata.js`), and
  every cell crossing the channel is an original position (a row of the
  dataset's `obs`), whatever share of the cells is open.

---

(reading-anndata-directly)=
## Reading AnnData directly

Every dataset, whatever format it arrived in (`data/input-formats.js`) and
whether it was opened from the file picker, a drop, a link or a server, is
read by Cellucid's own readers, in the AnnData reader worker
(`data/anndata-reader-worker.js`): nothing is copied into memory first, and
only what a view asks for is read. The code is layered like a package; each
layer depends only on the ones below it.

**What a file is** is decided in one place, the input format registry
(`data/input-formats.js`): each format's id (`h5ad`, `zarr`; the id a server
states as `anndata.file.kind` and its health `format`, and the prefix of the
id of a dataset that states none), the name endings of each layout it arrives
in (`file`: an `.h5ad`; `zip`: a `.zarr.zip`; `folder`: a `.zarr` store), the
container signatures (HDF5's at byte 0, 512, 1024 or 2048; ZIP's) and the
marker files of a store folder. Every entry point asks it: the File… picker's
`accept` list, a drop, the address bar's routing, the link providers, the
dataset catalog, the served-store validator, and the reader worker. A name
decides whether an entry point takes a file; once its bytes are at hand,
its content decides what it is. `tests/fixtures/input-formats.json` is the
registry's document byte for byte, and cellucid-python (`cellucid serve`) and
cellucid-r read the same rules from their copies.

**Opening** has one entry, `openDatasetStore({source, expect}, options)` in
`data/anndata/open.js`. The source is one file `{file}` (a `File` or a remote
file), a picked folder `{files, rootName}`, or a store folder on a web host
`{url, name}`; `expect` is what a server or a catalog says it is
(`{format, layout}`). It opens the container (`formats/hdf5`,
`formats/zarr`), then `presentStore()` asks the layout views of the formats
that container can hold (`data/anndata/views/index.js`) whether they
recognize the root: a view presents its format as an AnnData store root
without copying it (`views/view-nodes.js`, the nodes `older-layouts.js` is
made of), and a root no view recognizes is AnnData as stored. A format that is
one stream is turned into an AnnData Zarr store by its transcoder
(`data/anndata/transcode/index.js`) and opened with the Zarr reader. The
views are Loom, a 10x `.h5` and a Cell Ranger `outs/` folder, and MuData
(`.h5mu`, MuData Zarr); the transcoders read R objects (`.rds`, `.rda`) and
MatrixMarket folders into the origin-private file system (memory where the
browser has none), holding a block of entries at a time, never the matrix
(`transcode/entries.js`): a Seurat v5 `scale.data` layer is never stored,
split layers are joined in one pass, and X is laid out by gene
(`csc_matrix`, `byFeature`: buckets of consecutive genes of at most 4 M
entries, dealt in one pass and laid out one at a time), so a gene is one
contiguous read as in a Cellucid-ready AnnData; the other layers keep the
form they arrived in. Every per-cell array is written with its cells in
priority order under seed 0 (`transcode/priority.js`): `obs` with
`cellucid_original_position`, `obsm`, both axes of `obsp`, every layer's
rows (`majorInOrder`, bounded like `byFeature`), each gene's rows ascending
with `varm['cellucid_prefix_counts']`, each modality's rows by their cells'
priority, and `uns['cellucid']['priority']`, so a subsample of a converted
file opens as a prefix (the subsampling specification's priority-ordered
files); an uncompressed R file, whose arrays are read in place and never
copied, keeps the object's order. The store is kept across visits under the source's identity
(`transcode/space.js`, `session.js`): the same file opened again opens the
store without reading the file. The opening says up front that such a file
cannot be read in parts, how much is read and where it is kept
(`data/whole-read-notice.js`). A file whose host serves only whole files is
downloaded into the same storage as it arrives, never into memory, and kept
(`data/whole-download.js`). A large gzip stream, and a MatrixMarket file's
entry lines, are decoded off the reader thread (`formats/rds/gunzip-worker.js`,
`formats/text/matrix-market-entries.js`; {doc}`21_decisions_decoders_and_codecs`).

| Layer | Code | What it does |
|---|---|---|
| Bytes | `formats/io/` | Byte sources (a local `File`, HTTP range requests, memory) and **one byte reader per file**: reads asked for in the same task are planned together, neighbors closer than the merge gap become one request, requests are split at the request size so large reads run in parallel (`concurrency` at once), requests waiting for a slot take their neighbors with them when round trips cost more than bytes, each byte in flight is fetched once, metadata is cached in aligned blocks (16 KiB from a web host), and a read's `AbortSignal` cancels it. The merge gap, request size and concurrency follow **the link as measured** (`link.js`). A host that content-codes a file even in a range (`Content-Encoding`, hidden from the page) is told by the probe (`probeHttpResource`), and the file is downloaded whole. `inflate.js` decompresses DEFLATE with the platform's `DecompressionStream` or fflate, `decoder-choice.js` measures which decoder of a codec family is fastest on this machine, and `decode-pool.js` decodes chunks on decode workers (`decode-worker.js`), so a background pass over a matrix leaves the reader thread free; `runOnDecodeWorker` gives the same workers the reader thread's work on decoded values. `activity.js` counts what each operation fetched and decoded, and which shared loads it reads for (`keys`: a job joining one adopts the read, `read-activity-log.js`). How requests are planned, and why: {doc}`22_decisions_link_model_and_request_planning`. Which decoder decodes a chunk: {doc}`21_decisions_decoders_and_codecs`. |
| Formats | `formats/hdf5/`, `formats/zarr/`, `formats/zip/` (and, for the transcoders, `formats/rds/`: R's serialization, gzip native, bzip2 and xz in JavaScript, zstd in WebAssembly; `formats/text/`: MatrixMarket) | Cellucid's own HDF5 reader (every superblock, object header, group, attribute and chunk-index version; deflate, shuffle, Fletcher-32, LZF and hdf5plugin's Blosc, LZ4, Zstandard and bitshuffle), Zarr reader (v2 and v3 with sharding; gzip, zlib, Blosc, LZ4 and Zstandard; a picked folder, a ZIP, or a folder on a web host via consolidated metadata) and ZIP reader. Zstandard, LZ4 and Blosc are decoded by their own libraries compiled to WebAssembly (`assets/wasm/codecs.wasm`, built by `scripts/build-wasm-codecs.mjs`). Each opens a file as the **array store** of `formats/store.js`: a store's tree is read as it is navigated into, a group at a time (`group.load(paths)`: the groups on the way and the node at each path to a depth, each object read once however many readings ask, a tree level's reads in flight together), and navigation is synchronous over what is loaded (a group navigated before it is loaded throws `StoreNotLoadedError`). Opening an HDF5 file reads its root and what its format opens with (`ANNDATA_OPENING_PATHS` in `data/anndata/elements.js`: every element's own node -- the group to depth 2: each dataframe column, each `obsm`, `varm`, `obsp`, `varp`, `layers` and `uns` entry, `raw` likewise -- and every expression matrix and `uns['cellucid']` whole; a layout view's `openingPaths`), named to `openHdf5` so they load with the root and before any opening step runs; what lies inside an element (a categorical's codes and categories, a nullable column's values and mask, a graph's arrays, nested `uns` entries) loads when a reading needs it. A Zarr store's metadata is described in one place (consolidated) or read in one batch, so it is loaded whole at open. Array values are read on demand, fetching only the chunks a selection touches, in parallel, each decoded once (a per-file decoded-chunk cache serves neighboring reads); of an unfiltered HDF5 chunk or a Zarr chunk stored as its values' bytes, only the run of bytes a selection needs is fetched (`partialChunks`). A column of variable-length strings reads the rows asked for alone (`readRows`), and a categorical column with more categories than open cells is coded by the categories its open cells hold ({ref}`partial-reading-strings`). |
| AnnData | `data/anndata/` | **One** AnnData reader for both formats: the element rules of the current layout (`elements.js`, `dataframe.js`, `sparse.js`; an element that breaks one is an `AnnDataFormatError`, which opening says to the person as the file's problem, with what to do), files from anndata 0.6/0.7 and SeuratDisk presented in the current layout (`older-layouts.js`), slabbed reads with read-ahead and exact conversions (`values.js`), gene columns and cell rows per layout (`matrix-access.js`), and opening plus error explanations (`open.js`, `views/`, `transcode/`). The reader class is `reader.js` (its state and lifecycle) with its method groups in `reader/`. |
| Worker | `data/anndata-reader-worker.js`, `data/anndata-reader-client.js` | The reader's methods over `postMessage`, results transferred; live `activity` and `transfer` events. |

**Genes and cells by layout.** Only what one gene (or cell) needs is read:

| Matrix | A gene (column) | A cell (row) |
|---|---|---|
| CSC | its own entries: indices and values fetched together | the cells' entries of every gene, read as a cell subset (each gene's entries once, keeping the cells'; a priority-ordered file opened as a prefix reads each gene's prefix alone), never an index of every stored row |
| CSR, over a link as fast as a disk (a local file, the same machine or network) | indexed by gene once in the background, at open, indices and values into a column store (in memory, or on disk when it does not fit); a gene asked for meanwhile is answered from the same pass, and every later one is a copy out of the store | its own entries |
| CSR, over a slower link | indexed when the first gene is shown, from the **indices** alone (kept in memory when they fit the `columnIndex` budget; otherwise a second pass fills a column store, values and all, on disk when the site's storage quota has room); a gene's values are then fetched as the ranges its entries occupy; once the ranges already read and those asked for together would cost as much as every value (a request costs a round trip of the measured link), every value is read once, after which every gene is answered from memory (never more than twice the cheaper choice in hindsight) | its own entries |
| Dense, chunked in strips of few genes | the chunks of its strip (kept for its neighbors) | the chunks of its row strip |
| Dense, in one piece or in wide chunks | indexed by gene once, streaming | its row |

A CSC matrix's column pointers are read while the dataset opens, beside the
obs and var indexes, so the first gene shown is one round trip, not two.

**Reading a large file in parts.** What a gene, a share of the cells or a
cell subset reads is described once, in the "What it does now" of
{doc}`23_decisions_partial_reading`:

| Topic | Where |
|---|---|
| Gene indexes per matrix, the first gene's cost said before it is read (`AnnDataReader.geneReadPlan`), a stopped or failed index pass | {ref}`partial-reading-gene-indexes` |
| The gene sample (`geneSampleOf`) and the page's top-up (`DataState.topUpGeneSample`) | {ref}`partial-reading-gene-sample` |
| A subsample's genes filled in segments, partial answers, the two forms of a gene | {ref}`partial-reading-genes-arriving` |
| Some genes in some cells (`getMatrixCellSubset`), priced requests, gene delivery to the analysis workers | {ref}`partial-reading-cell-subsets` |
| A subsample's neighbor graph | {ref}`partial-reading-edges` |
| The open file read again with other cells, grown or cut in place, progressive opening | {ref}`partial-reading-reading-again` |

### Graphs and embeddings arrive finished

For a neighbor graph the reader worker checks the `obsp` matrix as its lines
stream in (bounds, order, a diagonal left out, exact symmetry of topology and
weight) and builds its forward adjacency: each edge once, at its smaller
cell, with its partner and its `float32` relative strength, the forward
offsets, and the graph's largest weight and degree (`ConnectivityEdges`,
`data/connectivity-edges.js`). That is 8 bytes an edge and 4 a cell, the only
copy: its three arrays are transferred to the page, never copied, and the
build is part of the working set (`'graph build'`, `memory-budget.js`) while
it runs, so every thread's caches make room for it (float64 weights are held
whole until the strongest is known, 8 bytes an edge more; what a graph at
scale holds: {doc}`27_decisions_memory_model`). Which edges a pane shows is
derived from it on the page; nothing is stored in draw order. Each upper
entry of a line is an edge kept with its partner and exact weight; each lower
entry must match the next due edge of its partner's line, so one pass in
storage order checks the whole graph. The page checks only the payload's
shape against the graph summary. The adapter keeps only the neighbor graph
asked for last: it is the payload the viewer shows, held by reference, so
keeping it costs nothing, and a graph no longer shown is released and read
again if it is shown again. A neighbor graph is not read at open: a direct
dataset records that it has one and counts its edges when they are first
read.

The k-NN tool reads the rows it expands through from the file, by cell, in
the reader worker (`getConnectivityRows`, through the chunks that hold them;
a frontier whose rows lie in every chunk streams the graph instead, as a
subsample's edges do), so no thread holds every cell's neighbor list. While a
view shows the graph's edges the page holds them anyway, and a frontier is
then answered from them (`forwardEdgeRows`: two passes over every edge, the
same neighbors) whenever that is faster than the file by the rates both ways
measured on this page (`knnRowReader` in `app/ui/modules/connectivity-controls.js`). An index built ahead
of need (the cells of a CSC matrix, asked for before) does not take turns
with an extraction someone waits for. Edges at scale:
{doc}`24_decisions_graph_edges_at_scale`.

For an embedding the worker reads the selected `obsm` columns at their stored
precision, checks every value, normalizes them once with Cellucid's rule (or
fits a file's own coordinates to the view) and sends `float32` positions with
their display transform and per-axis bounds; nothing on the page scans them
again. Vector fields arrive scaled by their embedding's transform, with their
largest magnitude, read concurrently with the positions.

### Values, strings and older releases

- **Values.** Gene values, obs columns and nullable integers are float32, as
  the views draw them: every stored type converts as numpy's
  `astype(float32)` does (64-bit integers to the nearest float32, a float64
  beyond the float32 range to an infinity), so one outlier never makes a
  column unreadable. Boolean matrices read as 0 and 1. Embeddings and vector
  fields are coordinates: a value float32 cannot hold is refused there.
  Integer category labels beyond 2^53 are kept as their exact decimal text.
- **String obs columns.** A string or categorical column whose categories are
  not listed in the file's metadata (a direct AnnData string column) is listed
  with its categories unknown; they are found when the column is first read,
  so opening never reads a column's values.
- **Older AnnData releases** are read through views (`older-layouts.js`), so
  the element rules exist once.
- **Memory.** Every cache that decides what stays resident takes a disjoint
  share of one allowance from `data/memory-budget.js` (see
  {ref}`memory-budgets`); a refused allocation in any thread halves the
  allowance everywhere ({doc}`27_decisions_memory_model`).
- **Requests from the page.** The gene reads of one task are one request per
  matrix and form, each gene read once however many callers asked, and every
  read takes an `AbortSignal` that stops the worker's reads
  ({ref}`link-model-shared-reads`).
  **A shared read is canceled only when every reader waiting on it withdrew**
  (`utils/abort.js`; {ref}`link-model-shared-reads`).

### A Zarr folder on a web host that hides `_` and `.` paths

Jekyll (GitHub Pages without a `.nojekyll` file) serves no path segment that
starts with `_` or `.`. A Zarr v3 store still opens there (its consolidated
metadata is in the root `zarr.json`), but every chunk under anndata's
`obs/_index` and `var/_index` answers 404, which a Zarr reader must take as
an absent chunk: cell and gene names would read as empty, with no error. So
when a store is opened from consolidated metadata on a web host, the metadata
document of the first `_`/`.` node on each such path (`obs/_index/zarr.json`,
or `.zarray`/`.zgroup` in Zarr v2) is fetched once, by the first read of an
array under it and beside that read's own chunk requests, which wait for its
answer too: a host that serves everything costs one parallel request per such
node read and no round trip in sequence, and opening or listing a store
fetches nothing for it. A host that answers it as absent is refused, naming
the fix: an empty `.nojekyll` file at the root of the site, or the store
opened through `raw.githubusercontent.com` or a `gh:` link. (A Zarr v2 store
on such a host already fails at its `.zmetadata`.)

A picked or dropped `.zarr` folder is read file by file
(`formats/zarr/stores/directory.js`), each chunk a
`File.slice().arrayBuffer()`, and a `.zarr.zip` through the byte reader as
one file (its last 256 KiB read at open); how each is read:
{ref}`link-model-zarr-folders`.

### Seeing what is read

Every reader operation (opening, a field, a batch of genes, an embedding, a
gene index, the graph) runs under a `ReadActivity`. It knows how much is
left: the byte reader counts each request when it is planned and when it
arrives, and an operation that streams a whole array says its stored size up
front, so a running operation reports `x of y` bytes, its rate and the time
left (from steps where it counts them, as an index pass does). The worker
reports each live (at most every 250 ms) and when it ends;
`data/read-activity-log.js` keeps them for the open dataset and writes a
running line every 2 s for anything longer than a second, and one line when
it ends:

```text
[cellucid:read] Sparse X, indexed by gene — 98 MB of 215 MB (45%), 112 MB/s, ~2 s left, 41 requests, 5210 chunks decoded to 512 MB, 1.1 s (reading)
[cellucid:read] 2 genes of X — 137 KB read, 9 requests, 1 from cache, 16 chunks decoded to 192 KB, 6 ms (done)
```

`app/read-progress-notifications.js` shows the reads a person waits for that
nothing else shows as one progress notification, in plain words ("Reading
field cell_type and the neighbor graph · 4.1 MB of 9.3 MB · ~3 s left"), on
screen once the reads have run a second (`NOTICE_AFTER_MS` in
`app/notification-center.js`: the one rule for every card of work in
progress) and gone when they end. A failed read only leaves that card: what
asked for it (the opening, the field or gene, the job) reports the failure,
once, and the read's own words are in its `[cellucid:read]` warning line. Not
shown there: an opening's reads (the opening progress), the cells
added to a growing opening (the **Dataset** box), a job's reads (the job's
card) and background work (an index built ahead of need). The reader's own
words for each read stay in these `[cellucid:read]` console lines.

A remote file's running totals also feed the transfer log
(`data/transfer-log.js`, `[cellucid:transfer]` lines), which the **Dataset**
box (`app/ui/modules/dataset-card.js`) shows as one line, "5.5 MB read in 43
requests", beside where the dataset is stored (`data/anndata-storage.js`:
format, container, host or this computer, size).

### Verification and benchmarks

`tests/hdf5-reader-conformance.test.mjs` reads every object of
an HDF5 corpus (`scripts/generate-hdf5-conformance-corpus.py`) and every `.h5ad`
fixture with this reader and with h5wasm (the HDF5 C library, a pinned
development dependency) and requires the same answer;
`tests/zarr-format-conformance.test.mjs` does the same against values
zarr-python wrote (`scripts/generate-zarr-conformance-corpus.py`); the AnnData
corpus (`tests/fixtures/anndata/`, made by `scripts/anndata-corpus/`) holds
573 small files written by every anndata release line from 0.6 to 0.12 (and
SeuratDisk-style, h5sparse and h5py `libver='latest'` files), each with the
values anndata read back; `tests/anndata-corpus-conformance.test.mjs` reads
every one through the AnnData reader, `.h5ad` locally and over HTTP, Zarr as
an archive and as a picked folder, and compares every obs column, gene and
cell, embedding, graph and title (two deliberate differences, where
Cellucid reads more than anndata did, are listed in the test). The other
formats have corpora of their own: `tests/formats-corpus-conformance.test.mjs`
(Loom, 10x `.h5` and `outs/`, MuData; `scripts/formats-corpus/`),
`tests/rds-corpus-conformance.test.mjs` (Seurat and SingleCellExperiment
objects written by R, `scripts/rds-corpus/`, and the MatrixMarket folders of
`scripts/mtx-corpus/`). The R object corpus is written again with the
current CRAN and Bioconductor packages once a month
(`.github/workflows/r-corpus-drift.yml`, `scripts/rds-corpus/drift.sh`).

Benchmarks: `npm run bench:hdf5` (the HDF5 reader against h5wasm, on local
files), `npm run bench:anndata` (`scripts/benchmark-anndata.mjs`: the AnnData
reader against an earlier revision's, on generated `.h5ad` and Zarr v2/v3
files, local and over a simulated link; `--protocol http/1.1|h2` simulates a
host that answers six requests at once or all of them), `npm run bench:formats`
(`scripts/benchmark-formats.mjs`: every other format, Loom, a Cell Ranger
`outs/` folder, a MatrixMarket folder, Seurat and SingleCellExperiment
`.rds`/`.rda`, `.h5mu` and MuData Zarr, against the same data as a
Cellucid-ready `.h5ad`: first view, first gene, peak memory, and a stream
format's pass against plain inflate of its own bytes) and
`scripts/benchmark-real-links.mjs` (the whole app on real hosts, against an
earlier revision; `--corpus` takes the real-link corpus below). A corpus of
real public links of every source and format
(`scripts/efficiency-floor/public-links.json`: CELLxGENE from 2.5k cells to
the 2.3M-cell lung atlas, Zenodo, GEO `.h5ad`, `.h5ad.gz`, 10x `.h5`,
MatrixMarket, `.rds.gz` and text matrices, the Single Cell Expression
Atlas, Figshare, Hugging Face, GitHub, Dropbox, Google Drive, public S3, the
sample catalog) is opened by `scripts/efficiency-floor/public-links.mjs`
through the app's own path in Node over the real network (a browser's six
connections to an HTTP/1.1 host and its name): each step's time, bytes,
requests and share of the file (open and first view, a coloring field,
the gene list, genes), the host's behavior (HTTP version, ranges, CORS and
the headers it exposes, content coding, validators, rate limits), and every
request a second visit sends, which should be the version check alone.

---

(web-host-requests)=
## Requests to a web host: retries, stalls and failures

Every request to a web host, whatever it fetches, goes through one retry and
stall implementation (`fetchOnce` and `withRetries`, `formats/io/byte-source.js`),
and none is left waiting for ever. In short (every rule, its reason and
its measurements: {ref}`link-model-web-host-requests`):

- **Transient failures are asked again** (HTTP 408, 429 and 5xx but 501 and
  505, a 403 that names a rate limit, a failed or dropped connection), with
  waits grown from the link's retransmission timeout; a request is reported
  as failed once asked again at least three times over 100 seconds
  (`RETRY_AT_LEAST`, `RETRY_FOR_MS`) on a read its host already delivered
  to, after three repeats otherwise.
- **A wait the host names** (`Retry-After`, a spent request budget with its
  reset) is waited for exactly and is no failure; a served store's `503`
  "being prepared" shows as **Waiting for the host**.
- **A host's request budget paces and prices every request to it**
  (`formats/io/hosts.js`), a host that says it is busy is asked less at once,
  and a host that blocks this network is given a pace, kept across visits.
- **A stalled request is aborted and asked again at once** (the link's
  `ResponseClock`), except one to a notebook's kernel; a request queued
  behind another reader's is waiting, not stalled.
- **A file that changed is refused**, and so is an answer of other bytes
  (`changed`, `range-ignored`).
- **A failure that will stay is reported at once** (`worthAskingAgain`); a
  host that stops answering is left alone, not asked by every read, and a
  read that was delivered to before waits a refusal or an offline computer
  out.
- **A failure says whose problem it is** (`hostFaultSentence`, a host-side
  status said as the host's), and names the request: `https://…/suo.zarr/obsm/X_umap_3d/3.0:
  the host could not answer (HTTP 503). Asked 8 times.`
- **A sign-in never follows a redirect to another host.**
- **A whole download** resumes from the byte it stopped at where the host
  serves ranges.
- Each repeat is reported to the transfer log while it happens (`retrying`;
  the Dataset box's one "Retrying a read" notification), and bytes of a
  request still under way as they arrive (`arriving`).
- Cancellation: the caller's signal aborts the request in flight and any wait
  to ask again at once, and nothing is asked after it.

`tests/http-request-resilience-contract.test.mjs` holds the cases against a
fake `fetch`, and `tests/host-budget-contract.test.mjs` every budget header
form, a 429 and a 403 rate limit, a served store's 503 "being prepared", the
per-host stall watch, a redirect's headers and the probe's two versions.

## Loader layer (what is read)

Most of the “what do we read?” logic lives in:
- `cellucid/assets/js/data/data-loaders.js` (it accepts only `anndata://`
  URLs)
- `cellucid/assets/js/data/base-anndata-adapter.js` (what each path means in
  the file)

What a dataset generation reads (`loadDatasetGeneration`):

- the **identity**: dataset id, name, description, citation and source (from
  `uns['cellucid']`), counts (the gene count from the var axis's length, no
  gene name read), embeddings (`obsm`), vector fields, and whether there is a
  graph (`obsp`). Its id must be the one that was selected.
- the **obs inventory**: the `obs` columns, categorical or continuous, and the
  outlier quantile route of each categorical one.
- the **graph description** (when the file has a graph in `obsp`): edges are
  read when first shown.

No gene name is read before the first view. The dataset is published with an
empty var inventory, its gene list pending (`DataState.beginGeneListing`),
and the data mapping it opens with: its `geneNames` is the var index, or the
column a launch named (`geneIdColumn`); nothing guesses another. Right after
the first view, `DataState.listPublishedGenes` reads the gene list (`var`):
the names the mapping says, and nothing else (no other var column). The list
is checked against the identity's gene count and installed, and
`genes:listed` is emitted. Until it is installed, the gene search, the
analyses that need genes and their gene pickers show that names are being
read; sessions, h5ad exports, mapping changes and view carries wait for it
(`whenGenesListed`). A var column named in Data mapping is read when the
mapping is applied.

Then, on demand: an embedding's positions when a view first shows it, an obs
column when it is colored, a gene's values when it is chosen, and a
categorical field's outlier quantiles on its outlier slider's first move below
`100%` (never when the field is picked; `state/managers/outlier-quantiles.js`):
the field's stored array when the file stores them for the latent space in use
and they still fit its categories, otherwise computed in the reader worker
from the latent space (`latentOutlierQuantiles`), which is read once and kept
within the memory budget for every other field. Either is one cancellable job
with progress (`cancelOutlierQuantiles`).

---

## Caching layers (what is cached where)

Cellucid uses multiple caches; knowing which one you’re hitting makes debugging much faster.

| Layer | What it keeps | Lives, and what bounds it | Code | More |
|---|---|---|---|---|
| 1. HTTP cache (browser) | Catalog files (`datasets.json`, starting views), revalidated before every use; nothing of a dataset's bytes from a web host | The browser's | `data/data-source.js` (`fetchCatalogFile`, `CATALOG_FILE_CACHE`), `formats/zarr/stores/http.js` | below |
| 2. In-memory caches in the data layer | Dataset metadata per source (`DataSourceManager`); a catalog's `datasets.json` and each file's description; the reader's loaded tree of an open file and its decoded chunks (`formats/store.js`, `formats/io/decoded-budget.js`); the open store, reused by a reading of other cells in the same worker (`sourceKey`); the bytes of remote files every reader worker fetched, and a Zarr folder's objects (the page's byte cache) | The page (the byte cache), the reader worker (the rest); the `decodedBlocks` and `fetchedRanges` memory shares | `data/data-source-manager.js`, `formats/store.js`, `data/anndata-reader-worker.js`, `data/byte-range-cache.js`, `formats/io/range-cache.js` | {ref}`partial-reading-reading-again`, {ref}`link-model-kept-bytes` |
| 3. `DataState` field residency | Loaded obs fields and gene vectors, inactive ones least recently used released | The page; a share of the memory allowance, all inactive ones when the page is hidden | `state/managers/field/residency.js`, `data/memory-budget.js` | below |
| 4. Community annotation caches | What each Pull read, by Git blob id; the last Pull's commit and tree ids (IndexedDB); the working annotation session (localStorage) | Across visits | `app/community-annotations/file-cache.js` | below |
| 5. Marker cache | Marker-discovery results, with their statistic arrays | Memory, for the dataset load they were computed on | `app/analysis/genes-panel/marker-cache.js` | below |
| 6. Origin-private file system | Converted stores, whole downloads, the downloaded parts of remote files, kept gene indexes | Across visits, under the file's version; checked as read | `data/anndata/transcode/space.js`, `formats/io/kept-parts.js`, `data/anndata/kept-column-index.js` | below |
| 7. Host link memory | What each host delivered (round trips, rate), and a blocking host's pace | `localStorage`, across visits, 32 hosts | `data/host-link-memory.js` | below |

### 1) HTTP cache (browser)

Every file a dataset catalog names — `datasets.json` and a starting view —
is fetched by `fetchCatalogFile` with `{ cache: 'no-cache' }`
(`CATALOG_FILE_CACHE` in `data/data-source.js`); a caller cannot choose
another mode. The browser keeps its copy but asks the host before every use,
and an unchanged file costs a conditional request answered with 304.
Reusing a cached catalog without asking would be wrong, not merely stale: a
starting view is pinned by the catalog's digest, and hosts such as GitHub
Pages cache for minutes, so a new catalog paired with an old session would be
refused.

AnnData bytes are read with range requests. A Cellucid server answers them
with `Cache-Control: no-cache` and a strong `ETag`, and a range asked with
`If-Range` naming an older version is answered with the whole current file,
so bytes of a changed file are never mixed with old ones.

Every request into a Zarr store on a web host bypasses the HTTP cache
(`cache: 'no-store'`, `formats/zarr/stores/http.js`), as every range read of a
file does: WebKit does not finish reading chunk bodies it serves from its
HTTP cache to a worker (the decision and how it was found:
{doc}`22_decisions_link_model_and_request_planning`). The suite asserts the
cache mode (`http-request-resilience-contract`). Bypassing the cache also
means a browser never pairs a store's metadata with chunks it cached from an
earlier version of the store; a reload reads the store again.

A query string cannot version a store instead: GitHub Pages' CDN ignores it
(a never-seen `?v=` on a store file is answered `x-cache: HIT`), and it would
break presigned bucket URLs, whose signature covers the query. For up to ten
minutes after a store is rebuilt, the CDN may still answer some of its
objects from before the rebuild; a chunk that does not fit the new metadata
is refused with a decode or size error naming its URL, never waited on.

Debug tip:
- DevTools → Network → “Disable cache” while you debug.

### 2) In-memory caches in the data layer

- `DataSourceManager` caches dataset metadata per source.
- A catalog source keeps the `datasets.json` it read and each file's
  description; the reader keeps what it loaded of an open file's tree and its
  decoded chunks (`formats/store.js`).
- Opening another cell subset of the open file (the subsample control,
  `app/ui/modules/subsample-control.js`; a session's recorded cells, a
  progressive opening's steps) happens in the same reader worker
  (`AnnDataFileSource.reopenCells` → `AnnDataReaderClient.reopenCells`, a new
  adapter on the same client) and reads through the store already open
  (`anndata-reader-worker.js`, `sourceKey`) when the source names the same
  version of the same bytes; the open cells are grown or cut in place
  (`app/open-cells-resize.js`) when one opening's cells are a prefix of the
  other's. What is kept, when, and what it saves:
  {ref}`partial-reading-reading-again`.
- The page keeps the bytes of remote files its reader workers fetched, for
  every reader worker, and the objects of a Zarr folder on a web host under a
  prefix naming the store's version, so a dataset opened again in a new worker
  fetches nothing the old one fetched ({ref}`link-model-kept-bytes`: modules,
  bound and order of dropping).

### 3) `DataState` field residency (bounded in bytes)

`DataState` keeps loaded field arrays — obs fields and gene vectors — and
releases inactive ones least recently used once they exceed their budget, a
share of the memory allowance (`state/managers/field/residency.js`,
`data/memory-budget.js`), and all of them when the page is hidden. The active field, fields with filter
state, and fields still loading are always kept. Categorical codes
stay in the type the file stored (uint8 stays uint8). An inactive gene whose
column has a sparse form (nearly every gene: `data/column-forms.js`) keeps
that form alone, whatever the budget — its entries instead of four bytes a
cell — and is made again from it, bit for bit and
without a read, when it is shown (`compactedColumn`); so browsing genes holds
one gene's column at a time. Why it is bounded in bytes:
{doc}`27_decisions_memory_model`.

### 4) Community annotation caches (IndexedDB + the session in localStorage)

Community annotation stores:
- what each file a Pull read holds for the open dataset, keyed by its Git
  blob id, or why the Pull left the file out (IndexedDB); the blob ids a scope
  knows are the store's keys, so "which files changed?" needs no separate index
- the last completed Pull's commit and tree ids (IndexedDB), so a Pull that
  finds the same `annotations` tree stops after one or two requests
- the working annotation session itself (localStorage), which holds votes,
  comments, suggestions and profile — dataset-derived content, not preferences
  — and only the blob ids of the files it publishes onto (its own user file,
  the config and the merges file)

Code:
- `cellucid/assets/js/app/community-annotations/file-cache.js`

### 5) Marker cache (in memory)

Marker-discovery results, including computed statistic arrays, are cached in
memory for the dataset load they were computed on (the cache identity names
the dataset id, its selection revision and its gene axis, and is never
persistent): an AnnData file records no write time that would prove a stored
result still matches its values, so nothing is written to IndexedDB.

Code:
- `cellucid/assets/js/app/analysis/genes-panel/marker-cache.js`

### 6) Origin-private file system (kept across visits)

The reader worker keeps what is costly to read again in the browser's
origin-private file system: an R object or MatrixMarket folder converted
once, a file a host serves only whole, the parts of remote files its reads
downloaded (`formats/io/kept-parts.js`, written as journals as the bytes
arrive: {ref}`link-model-kept-bytes`), and the gene index of a
matrix stored cell by cell, built over the open cells and found again by a
later opening of the same cells of the same file (what each run is checked
against before use, and what that guarantees:
{doc}`19_subsampling_specification`, *Kept across visits*). Directories,
manifests, lifetimes and room-making are listed under *Origin-private file
system (OPFS)* in {doc}`04_configuration_env_vars_and_feature_flags`.

Everything kept is checked as it is read: a converted store, a kept
download and downloaded parts carry a CRC-32 of every 16 KiB page (damaged
parts are fetched again), a kept gene index a CRC-32 of
each gene's piece (`formats/io/page-checks.js`); damaged bytes fail the read
as damaged, the copy is removed and the next opening reads the file again.
Safari's workers have the origin-private file system as Chromium's and
Firefox's do; only a private window (any engine) has none, and then what is
converted stays in memory for the tab.

**Session › Stored** (`app/ui/modules/kept-stores-line.js`) sums up what is
kept; **Manage…** opens the Stored dialog
(`app/ui/modules/kept-stores-dialog.js`):

| Row | What it lists | Code |
|---|---|---|
| A dataset | every item kept of it (`space.js` `listKeptStores`: a converted store, a whole download, a file's downloaded parts, with the bytes they spare the next opening, `saves`, or a gene index per file, matrix and cell order), named by the manifest's `about`, its file's name and host from `keptAbout`, never a link, grouped by dataset (`file`); the rows chosen, a dataset's, or all, that no tab holds, are removed (`removeKeptStores`) | `data/anndata/transcode/space.js` |
| *Saved states* | the histories this browser keeps (`app/history/autosave.js` `listSavedStates`, IndexedDB `cellucid-history`): a row per saved state, its dataset, whether it is the history kept as it was made or a save under a name, its steps, size and when it was last written; ids `saved:<id>`, chosen like the kept items, removed with `removeSavedStates`, heard through `onSavedStatesChange` ({doc}`29_action_layer_domains_and_history`) | `app/history/autosave.js` |
| My datasets | one row (`app/library/my-datasets.js` `listMyDatasets`, IndexedDB `cellucid-library`, with the handles of files to reopen, `app/library/local-file-handles.js` `listLocalFiles`, IndexedDB `cellucid-local-files`): how many datasets and files to reopen, **Show…**, which opens My datasets' own dialog (each dataset is listed there and only there), and **Remove**, which forgets the list and the handles (`forgetMyDatasets`) | `app/library/my-datasets.js` |
| Connection speeds | the hosts whose connection speed is remembered (`data/host-link-memory.js` `rememberedHosts`, in `localStorage`, layer 7), with its own **Remove** (`forgetHostLinks`) | `data/host-link-memory.js` |
| Decoder choices | the decoders measured fastest in this browser (`formats/io/decoder-memory.js` `countDecoderDecisions`, `forgetDecoderDecisions`) | `formats/io/decoder-memory.js` |

**Remove all** forgets them all. The Stored row names the datasets something
is kept of and the bytes, and counts the saved states and connection speeds
as "more" (`summarizeKeptStores`: *2 datasets · 540 MB · 3 more*, its title
*2 datasets · 540 MB · 2 saved states · 1 connection speed*; with no dataset
kept, the one kind kept or how many items; it does not count My datasets,
which has its own control). Listing and removal run one at a time in a page:
each asks whether a tab holds a directory by requesting its lock if free,
which, while it lasts, is a holder to any other such request. The dialog's
foot (*Where it's kept*: *Holds* files, saved states, My datasets, connection
speeds) comes from `keptStoragePlace`: whether the window gives the site an
origin private file system at all (asked for, never inferred from the engine:
`navigator.storage.getDirectory()` answers or refuses),
`navigator.storage.estimate()` (usage and what the quota leaves, labeled as
the browser's figures) and `navigator.storage.persisted()` (*Kept when space
is low*); `keepKeptStoresWhenSpaceIsLow` calls `navigator.storage.persist()`
only on the person's click (**Keep it when space is low**). Without storage
the Stored row says *Kept for this tab only* and the dialog lists no dataset.

#### When a dataset changes

What is kept of a file, of every kind, is filed under its **version**, one
rule for all of them (`space.js` `persistsAcrossVisits`):

- a name of the content where the source has one (a CELLxGENE Census
  build, a checksum, a commit or Git LFS object id, a Drive revision, a
  OneDrive cTag);
- else the host's validators: the file's exact size with the
  `Last-Modified` date and the strong `ETag` the host sends
  (`writeVersion`: `written:<size> modified:<date> etag:<tag>`; a weak
  `W/` tag adds nothing);
- a file served by `cellucid serve` or `cellucid_serve()`: its strong
  `ETag` (size, times, inode);
- a file on this computer: its size and modification time;
- none of these (`size:<n>`): nothing is kept for a later opening (what an
  opening converts or indexes of it is used by that opening and its
  readings again with other cells only), and nothing kept is removed: an
  opening that cannot tell its version (GitHub past its API's limit) only
  uses nothing (`space.js` `forgetKept`).

Every opening asks the source for its version before anything kept or
cached of it is used: the probe's `HEAD` past the HTTP cache (`no-store`),
the provider's API, or a server's dataset list, read at that opening, never
a listing remembered from an earlier one. An opening from the source's own
description makes that check itself, in the same opening: the probe (a
`HEAD` and the head read), or the provider's answer that names the version
(a record, a listing, an item). A reader opened again from an earlier
opening's description (the dataset no longer shown and then read again, a
page back from the back/forward cache: `resumed`,
`anndata-reader-client.js`) asks it explicitly first, whatever is kept
(`kept-parts.js` `confirmFileVersion`, one small request a file): a provider
that alone can tell (Google Drive's public reads: the item's checksum;
OneDrive: the item's cTag) is asked; a file served only whole has its
validators asked with a `HEAD`; any other file has its first byte read from
the host, past every kept copy and the page's memory, under the validators
it was opened with. A version that no longer holds forgets every kept copy
of the file and fails the opening as changed ("open it again"); a reading
again with other cells within one opening (Apply, a growth step) keeps the
version its opening checked.

Before reading, `open.js` (`forgetOtherVersionsOf`) forgets what is kept of
the file's other versions, every kind (`space.js` `forgetOtherVersions`:
every directory of the same address, `addressOf`, whose `fileKey` is not
this version's loses its manifest and is removed once no tab holds it), and
the page's memory of their bytes (`range-cache.js`
`forgetOtherVersionsInMemory`). A file that changes on its host while it is
open is refused by its byte source (`changed`: "open it again to read the
new version"); the kept parts then forget every version of the file, every
kind (`forgetKeptFile`), and the page's memory drops it
(`forgetFileInMemory`).

A dataset is one address: GitHub and Hugging Face links name their ref in
it (`gh:owner/repo@<ref>/path`, `hf:repo@<revision>/path`), so two refs are
two datasets and never replace each other's kept copies; a link's address
is taken in its public form (`data/sources/link-credentials.js`
`publicLink`: the address and only the query parameters on an allow-list of
those that name a file or its version and grant nothing, such as
`versionId`, `generation`, `alt` or `download`; never a user name,
password, fragment or any other parameter). What a link *is*
(`linkIdentity`: the public form plus a short digest of the left-out
parameters that are not known credentials, so `download?uuid=a` and
`download?uuid=b` stay two datasets) is the one rule for a dataset id, a
kept store's key, the page's memory of bytes and a served store's key; known
credential and signing parameters are in neither, so a link signed again is
the same dataset.

A reading of the open dataset with other cells (Apply, a growth step, a
session's cells) keeps what is shown only when it is of the same version
of the same bytes: the reader worker took over the store of the reading
before (`anndata-reader-worker.js` `sameVersion`; a remote file by its
version, a folder at a content-named version by it, any other folder while
its metadata answers unchanged). A file with no
version is the same bytes within the opening that read it, so its growth
steps and Apply take over its store and read only the added cells; a
reader worker started again (a page restored from the back/forward cache)
is another opening, and what is shown is read anew
(`AnnDataFileSource.reopenCells`). Old values are never shown beside new
ones.

Every kept directory is stamped with the code fingerprint of its kind
(converted stores, gene indexes, files downloaded whole, the downloaded
parts of files), written by the build into
`transcode/kept-code-fingerprints.js` (`scripts/kept-code-fingerprints.mjs`
names each kind's entry modules; comments and layout do not count). A kind
covers every module its entries reach through imports and, inside a class
composed from several modules (the AnnData reader, `utils/apply-mixin.js`),
the module of every method a covered one calls through `this`: the gene
index's covers `reader/kept-indexes.js`, `reader/gene-segments.js` and the
row pointers they write (`reader/row-pointers.js`). A directory
stamped by other code is removed, never read, and what it held is made
again on the next opening. A file downloaded whole and a file's parts are
the file's own bytes, so their fingerprints cover only `kept-format.js`,
`page-checks.js`, `sources/identity.js` and `sources/link-credentials.js`
(the parts' also `kept-journal.js`): a release that changes readers,
decoders or transport keeps every download.

| Source | Its version | Asked at each opening by |
|---|---|---|
| A web link, S3 and its mirrors, Dropbox, GEO | size, `Last-Modified`, strong `ETag` | the probe (`HEAD`, `no-store`) |
| CELLxGENE | the Census build (its files are a folder of that build) | the release list (`no-cache`) |
| Zenodo, Figshare | the version record and the file's checksum | the record (`no-cache`) |
| Hugging Face | the file's LFS hash or blob id, read at the commit that last wrote it | `paths-info` for the path, beside the revision's commit (`no-cache`); never the repository's listing |
| GitHub | the commit that last wrote the file, or its LFS object id; for a Zarr folder, the commit that last wrote anything below it, read at that commit and filed under the link's address (past the API's 60 requests an hour: `size:` or no version, nothing kept across visits); a dataset of a catalog on GitHub or Hugging Face is opened by its provider as a link to it is (`sources/anndata-catalog.js` `providerSourceOf`: the same identity and version whichever way it is reached; its listing asks no provider) | the commits API |
| Google Drive, OneDrive | the md5 or revision (signed in, the revision the item names is what is read; public reads, whose key may not read a revision, are kept only once the item confirms it again); the cTag | the Drive metadata; the Graph item |
| `cellucid serve`, `cellucid_serve()`, a notebook | the file's strong `ETag` (its key names the server and the dataset); no parts are kept (the file is on this machine or in the kernel) | the server's dataset list (`no-store`) |
| A folder picked from this computer (a Zarr store, MatrixMarket, Cell Ranger `outs/`) | every file's path, size and modification time | the browser's file list |
| A few files of one dataset by link (MatrixMarket, `outs/`) | each file's own version | one round of probes, all at once |
| A Zarr store folder read at a commit (Hugging Face) | its git tree id, and each kept object's own validator | `paths-info` |
| A Zarr store folder on a plain host (GitHub Pages, S3, GitHub by branch) | none across openings: kept for its opening only | — |
| A `.zarr.zip` | as a file by link | the probe |

What this cannot notice: a file replaced by other content of exactly the
same size within the same second (`Last-Modified`'s resolution), on a host
that sends no `ETag` and whose source names no content. A GitHub link
opened past the API's 60 requests an hour has no version (`size:`) and
opens without anything kept.

Code:
- `cellucid/assets/js/data/anndata/transcode/space.js`
- `cellucid/assets/js/formats/io/kept-parts.js`
- `cellucid/assets/js/data/anndata/kept-column-index.js`
- `cellucid/assets/js/formats/io/page-checks.js`
- `cellucid/assets/js/app/ui/modules/kept-stores-line.js`

### 7) Host link memory (localStorage)

An opening sizes its gene sample and a progressive first view at the speed
of its link, and its first requests are small, so what it measures is a
lower bound. So what each host delivered is kept across visits
(`data/host-link-memory.js`, `localStorage` `cellucid:host-link:v1`): by
origin (`scheme://host[:port]`, never a path, a dataset, a query or a token),
the round trip, the round trip under load, the rate its reads sustained, the
bytes they rest on and when they were measured; at most `HOST_LINK_RECORDS`
(32) hosts, the least recently measured dropped. The reader prices with
`pricingLink`: the opening's own measurement once it has measured its
bandwidth or rests on as many bytes as the record, the record before that. A
host on this machine (`onThisMachine`: a loopback origin) is priced as this
machine at rest (`thisMachinePricing`). `rememberedHosts` lists the hosts,
and `forgetHostLinks` removes them (the Stored dialog's *Connection speeds
remembered* row, and its **Remove all**). How a record is blended and used,
and why: {ref}`link-model-remembered-links`.

Code:
- `cellucid/assets/js/data/host-link-memory.js`

:::{important}
Layers 4 and 6 persist data derived from the user’s dataset, and layer 7
the hosts data was read from (layer 5 lives only in memory), so they belong
in the privacy model, not just the performance story. Anything you add here
that outlives a reload must be documented in
{doc}`../o_accessibility_privacy_security/02_privacy_model`.
:::

---

## Debugging checklist (data loading)

When something “doesn’t load”, do the following in order:

1) **Identify the source type**
   - local-demo vs local-user vs remote vs link vs jupyter

2) **Find the base URL**
   - In console:
     - `window._cellucidState` exists after bootstrap
     - The dataset base URL (`anndata://<id>/`: the address the data loaders use for the dataset's parts, not a fetchable one) is held by the dataset publication (`app/dataset/publication.js`)

3) **Watch the Network tab**
   - Look for:
     - `404` (wrong path)
     - CORS errors (no `Access-Control-Allow-Origin`)
     - HTML returned for a JSON/binary request (“Unexpected token `<`”)

4) **Confirm the file opens**
   - The `[cellucid:read]` console lines name each read operation and its
     outcome; an open that fails says what the file lacks (for example
     nothing in `obsm` to plot) or which `uns['cellucid']` key is wrong.

5) **Confirm point buffer sizes**
   - After load, check:
     - `window._cellucidState.pointCount`
     - `window._cellucidState.positionsArray.length === pointCount * window._cellucidState.positionsDimension`
       (positions hold the embedding's own components: 1, 2 or 3 per cell)

6) **Check decompression support**
   - Cellucid requires native `DecompressionStream` (gzip and deflate) to read
     compressed chunks and session files.
   - A browser without that API receives a terminal capability error before the data request.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Dataset loads forever”

Likely causes (ordered):
1) The host stopped answering: requests that stall or fail are asked again
   ({ref}`web-host-requests`), each repeat a `[cellucid:transfer]` warning
   naming its URL; the startup error names it once nothing has arrived for
   60 s.
2) A request is blocked (CORS/mixed content): it fails as a failed
   connection, is asked again, then reported naming its URL.
3) A large read on a slow or throttled host (the `[cellucid:read]` lines show
   its progress).

How to confirm:
- Network tab: which request is pending?
- Console: enable debug and look for the last log line.

Fix:
- For remote datasets: validate CORS + protocol (https vs http).
- For GitHub datasets: verify raw URLs are reachable from the target environment.
- For a served AnnData store: watch the `_cellucid/anndata/...` range
  requests in the Network tab; the `[cellucid:read]` console lines name the
  operation that is still reading.

### Symptom: a host answers with a web page

A host that answers a path it has nothing at with its own page, as `200
text/html` (a single-page site: Cloudflare Pages without a `404.html`, a
Netlify or Vercel rewrite to `index.html`), or a link that answers with a
sign-in or preview page, is told apart by the answer's type and its first
bytes (`formats/io/byte-source.js` `isWebPageAnswer`, `looksLikeWebPage`):

- a Zarr object answered with a page is absent, as a 404 says
  (`formats/zarr/stores/http.js`): a chunk zarr-python left out (it holds
  only the fill value) reads as its fill value, and a folder with no store at
  it says so;
- a catalog file (`datasets.json`) answered with a page is not found
  (`fetchJson`), never a network failure to try again;
- a data file whose link answers with a page, or a page saved under a data
  file's name, is said to be a web page, with what to do (`httpByteSource`,
  `data/anndata/open.js`), not "not a valid HDF5 file".

Fix: fetch the URL in a new tab; a data file or Zarr object must come back as
its bytes, and a missing one as a `404`.

### Symptom: a load reported “invalid format” but the file is fine

Read the message again — this class of failure is split four ways, and
“invalid format” means only the last one:

| Message contains | Cause |
|---|---|
| “was canceled” | An abort — a newer selection superseded this load. Re-thrown untouched so it is never relabeled. |
| “is longer than the longest text this browser can hold” | A JSON document (a catalog, Zarr metadata) exceeds the engine's longest string |
| “response body transfer failed” | The connection dropped mid-body |
| “must contain valid JSON” | Genuinely malformed content |

If you add a new bounded or cancellable read, preserve this: check for an abort
first and re-throw it unchanged, distinguish the range failure from the parse
failure, and give each its own error code; one message for all of them would
leave these failures undiagnosable.

### Symptom: “Gzip decompression not supported”

Likely cause:
- The browser does not provide native `DecompressionStream('gzip')`.

Fix:
- Update to the current stable release of Chrome, Edge, Firefox, or Safari.

---

Next: {doc}`10_sessions_persistence_and_serialization` (how loaded state is saved/restored).
