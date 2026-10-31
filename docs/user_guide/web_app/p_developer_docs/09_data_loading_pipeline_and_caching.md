# Data loading pipeline and caching

This page documents the **developer-facing data loading pipeline** for the Cellucid web app:

- how datasets are discovered and selected
- how the app derives a dataset **base URL**
- which files are fetched (and in what order)
- where caching happens (HTTP cache, in-memory caches, byte-bounded field residency, IndexedDB caches)
- how remote server / links / file picker / notebook and R viewer modes differ

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
| `obsp` | the neighbour graph the file opens with (`loadConnectivitySummary`, `loadEdges`) |
| `candidates` | the data-mapping candidate catalog (`loadMappingCandidates`) |
| `obsm/{i}/{cols}` | `obsm` entry `i` read at some columns (an embedding) |
| `vectors/{v}/{vcols}/{e}/{ecols}` | a vector field drawn on an embedding |
| `obsp/{p}` | another neighbour graph |
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

At startup, `cellucid/assets/js/app/main.js`:

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
  names, is read (and its values checked) when a view first shows it. Before
  its first render, Suo reads 5.5 MB in 27 requests (its metadata and its
  3D embedding), not the 10.9 MB in 139 requests that reading all three
  embeddings costs; over HTTP/2 at 80 ms and 2 MB/s its first render takes
  3.8 s in WebKit (6.9 s reading all three)
- optionally describe the neighbour graph (edges load when shown)
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

A stream reader is read to its end, or cancelled, and dropped: no browser code
calls `releaseLock()` (`tests/browser-stream-reader-release-contract.test.mjs`).
WebKit (314030@main) runs the release's JavaScript while it holds a lock its
garbage collector takes when it visits the reader; a collection that starts
inside the release never finishes, and every later collection in the tab,
the page's included, waits behind it. The inflation of a compressed chunk on
the reader thread (`formats/io/inflate.js`) released its reader, so in WebKit
a file could stay "Opening …" for good and the page then hung on close.

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
- **A connection may or may not change the dataset**, which is known only once
  the server answers. Every user action that decides the dataset advances
  `actionSequence`: each transaction begun, and each session file applied
  (`recordAction()`, which begins no transaction and aborts nothing). The
  address bar's **Open** reads `actionSequence`; a one-dataset activation then calls
  `beginIfNoActionSince(sequence)`, which supersedes loads that were running at
  the click but never a choice made or a session applied after it. A
  connection that lost that race is kept as a catalog.
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
  (`app/local-file-handles.js`). The
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
- A dataset whose catalog entry names a `session` gets that default session
  applied after it is published (`app/session/catalog-default-session.js`;
  see {doc}`../l_sessions_sharing/04_official_sample_states`).

Common failure mode:
- `datasets.json` missing or invalid → the sample list is “not available”,
  and the message names what is wrong (for a catalog without `schema`:
  `Invalid datasets.json: the catalog must contain exactly: datasets,
  default, schema (missing: schema).`).

Cross-origin transport contract:
- `datasets.json` and every default session are fetched by
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
- The registry picks the provider that recognises it. A provider only says
  where the bytes are: a dataset catalog folder (a folder URL holding
  `datasets.json`; GitHub resolves to
  `https://raw.githubusercontent.com/owner/repo/<branch>/path/`), one file of
  a format of `data/input-formats.js` with its exact size, the files that
  make one dataset together (a GEO sample's MatrixMarket triplet), or a Zarr
  store folder (read through its consolidated metadata).
- A catalog then behaves like `local-demo` (the same `anndata-catalog.js`
  reader; each listed file opened when chosen, its default session applied).
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
  header, so no preflight either): Chromium lets one request of a URL
  through its cache at a time, so ranges of one file asked for together went
  out one round trip after another (five 100-byte ranges of one URL: five
  round trips in Chromium, one in Firefox and WebKit; one in all three with
  `no-store`). File sizes come from `HEAD` or the provider, because many hosts
  do not expose `Content-Range`.
- An `.h5ad` opens one round trip after its link is followed: the `HEAD` that
  asks its size goes out with a read of its first 512 KiB (`HEAD_BYTES`,
  `probeHttpResource`; GitHub's provider sends it beside its own `HEAD`), and
  the reader starts from those bytes (`createByteReader`'s `head`, which also
  measures the link). A Cellucid-ready `.h5ad` keeps everything opening reads
  there (the specification's "Where things sit in the file"): Suo over 80 ms
  at 50 MB/s, the probe and opening 327 → 98 ms (7 requests in three dependent
  waves → the `HEAD` and the head in one), and in the app in Chromium the
  first embedding request 0.81 → 0.22 s after the link is followed, the
  dataset ready 1.9 → 1.3 s. A global heap collection read through bytes the
  reader holds is read only as far as they go (`cachedEnd`): reading on for a
  round trip's worth would cost the round trip it saves.

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
- `cellucid/assets/js/app/local-file-handles.js` (handles kept to reopen, Chromium)

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
    <name>*, then applies the session). The reopen (`reopenLocalFile`) runs
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
- `cellucid/assets/js/app/viewer-host-results.js` (what the person made,
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
recognise the root: a view presents its format as an AnnData store root
without copying it (`views/view-nodes.js`, the nodes `older-layouts.js` is
made of), and a root no view recognises is AnnData as stored. A format that is
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
(`data/whole-download.js`). A large gzip stream (1 MiB or more) is
decompressed on a worker of its own (`formats/rds/gunzip-worker.js`), since
`DecompressionStream` inflates on the thread that calls it; a MatrixMarket
file's entry lines go to the decode workers in blocks of about 1 MiB and come
back as typed arrays (`formats/text/matrix-market-entries.js`), taken in file
order and collected a block at a time.

| Layer | Code | What it does |
|---|---|---|
| Bytes | `formats/io/` | Byte sources (a local `File`, HTTP range requests, memory) and **one byte reader per file**: reads asked for in the same task are planned together, neighbours closer than the merge gap become one request, requests are split at the request size so large reads run in parallel (`concurrency` at once: six per HTTP/1.1 host, the browser's limit), requests waiting for a slot take their neighbours with them when round trips cost more than bytes, each byte in flight is fetched once (a read of bytes a request already carries waits for that request, and no plan fetches a gap across bytes in flight; a request that fetched a gap no read asked for is kept a while within the cache's budget, so the rest of a heap collection read in part is answered from it), metadata is cached in aligned blocks (16 KiB from a web host: an anndata-written `.h5ad` keeps each element's headers within a few KiB beside its data, so larger blocks fetched mostly data; a 1M-cell `.h5ad` opened at 80 ms/2 MB/s in 1.4 s instead of 4.2 s), and a read's `AbortSignal` cancels it. A host that content-codes a file even in a range (`Content-Encoding`, hidden from the page; NCBI does so for GEO's `.rds`) is told by the probe (`probeHttpResource`: the ranged read fails while the whole file answers), and the file is downloaded whole. The merge gap and request size follow **the link as measured** (`link.js`, see below). `inflate.js` decompresses with the browser's native `DecompressionStream`, to an exact or bounded size, feeding it the input in 32 KiB slices so the size bound holds while it decodes (WebKit's decompressor inflates whatever is written to it at once: a chunk written whole expanded a 60 MB bomb past its 1 KiB declared size), or, for a chunk of known decoded size, with fflate 0.8.2's inflater (`assets/external/fflate-0.8.2/`, vendored unchanged and pinned by digest), into an output one byte larger than the declared size and checked as the stream checks it (Adler-32 or CRC-32 and size); fflate cannot seek, so `deflate-seek.js` keeps the streams read in parts. Which of the two inflates is measured, not named (`decoder-choice.js`; Zstandard, LZ4 and Blosc choose between WebAssembly and JavaScript the same way): per codec family and size class of the decoded chunk. Until a size class is decided, every chunk is decoded with the engine's prior (`ENGINE_PRIORS`: the decoder measured fastest on that engine for a read; Firefox: fflate and the JavaScript Zstandard, LZ4 and Blosc decoders; Chromium: the JavaScript Zstandard and Blosc decoders, whose WebAssembly is fast only once V8 has optimised it in each worker, which a gene's read does not reach (10M cells, first gene zstd 175–244 ms in JavaScript against 386–678 cold WebAssembly); WebKit: the references), never with an untested decoder, so the first gene of a first visit costs what the best known decoder costs. The decoders are measured only once the pool has had no other work for a second (`EXPLORE_AFTER_IDLE_MS`: a read's chunks come in waves, and a measuring message between two would delay the next): it keeps copies of a few chunks of each undecided size class (one per worker; chunks that only travel in batches are not kept and keep the prior) and sends them, one message a worker at a time, in phases (eight per worker, every one of a phase with one decoder, so each is measured under its own copies' contention on the other workers), first one warm-up phase per decoder (the engines compile it tier by tier), the prior first, then two phases each in turn, every phase counted, the warm-up too, so a decoder fast only once optimised is charged its cold cost, which a person's short reads pay; a read's chunk sent meanwhile voids the phase under way, which runs again later. A message's time is its worker's busy time for it, measured in the pool's thread (from its sending, or the worker's previous answer if still busy, to its answer: waiting for work is not counted, travel is); the decoder with the most decoded bytes per busy millisecond is kept. Decisions are kept in this browser (`decoder-memory.js`, IndexedDB `cellucid-decoders`, under a version naming the WebAssembly module's digest, the decoder lists, the user agent and the core count; any change measures again), so only a first visit measures; the Stored dialog lists and forgets them (*Decoder choices measured*). A decoder other than the reference is checked once against it on one chunk of the measuring messages (never on a read's path: a read's chunk is decoded once) and dropped if it differs or throws. Node (four worker threads, real chunks, each decoder fixed in the pool against the probe free): the probe kept the pool's fastest in all six cases measured (h5py-sized DEFLATE native 621 vs fflate 572 MB/s, zlib 538 KB 1,429 vs 454, Zstandard WebAssembly 2,726 vs JavaScript 1,209, LZ4 8,045 vs 2,518, Blosc-LZ4 4,498 vs 3,002, Blosc-zlib native 1,769 vs fflate 340). A stream's fixed cost contends across threads in some engines (Firefox's eight workers inflated h5py's 34 KB chunks at 160 MB/s in total, one at 300), which a choice made on one core or by engine name gets wrong. `decode-pool.js` decodes chunks on decode workers (`decode-worker.js`, nested in the reader worker: one started when a file opens, the others as chunks need them, up to as many as the cores allow, less one for the reader thread and one for the page, and as many as the `decodeWorkers` memory share holds), small chunks that wait together going to a worker in one message, a background pass's chunks behind every other read's, so a background pass over a matrix leaves the reader thread free; where no worker can be made, chunks are decoded in the reader thread. The same workers take the reader thread's work on decoded values (`runOnDecodeWorker`: a function of a module that depends on nothing else, its arguments transferred, with a context each worker is sent once), as a cell subset's runs of genes. `activity.js` counts what each operation fetched and decoded, and which shared loads it reads for (`keys`: a job joining one adopts the read, `read-activity-log.js`). |
| Formats | `formats/hdf5/`, `formats/zarr/`, `formats/zip/` (and, for the transcoders, `formats/rds/`: R's serialization, gzip native, bzip2 and xz in JavaScript, zstd in WebAssembly; `formats/text/`: MatrixMarket) | Cellucid's own HDF5 reader (every superblock, object header, group, attribute and chunk-index version; deflate, shuffle, Fletcher-32, LZF and hdf5plugin's Blosc, LZ4, Zstandard and bitshuffle), Zarr reader (v2 and v3 with sharding; gzip, zlib, Blosc, LZ4 and Zstandard; a picked folder, a ZIP, or a folder on a web host via consolidated metadata) and ZIP reader. Zstandard, LZ4 and Blosc (Zarr's codecs and hdf5plugin's filters alike) are decoded by their own libraries compiled to WebAssembly: zstd 1.5.6, lz4 1.10.0 and c-blosc 1.21.6 (with the blosclz and bitshuffle it carries) in one module, `assets/wasm/codecs.wasm`, built by `scripts/build-wasm-codecs.mjs` from the release tarballs pinned by SHA-256 with zig 0.13.0 (the build is byte-for-byte reproducible; `tests/wasm-codecs-contract.test.mjs` pins the module's digest and rebuilds it where the build cache holds the inputs), licences in `assets/wasm/LICENSES/`. `formats/io/wasm-codecs.js` loads it once per thread at the address `assets/workers/scripts.js` gives (`?v=` its bytes' version, kept by the service worker), decodes into exactly the size each format declares, and lets go of an instance whose memory grew past a decode worker's budget. A Blosc frame of zlib streams (what `prepare()` and `cellucid_prepare()` write for a shuffled Zarr v3 array) is read by c-blosc too, but its streams are inflated by the page's DEFLATE decoder (below), each into its place, and its blocks byte-unshuffled in JavaScript (`bloscZlib`): no zlib is built into the module. The HDF5 bitshuffle filter (32008: bit-transposed blocks, stored as they are or each LZ4- or zstd-compressed, the elements after the last block as they are) is undone in one call (`bitshuffleChunk`), and an R file's zstd stream piece by piece (`zstdStream`, in an instance of its own). The JavaScript decoders written for this reader (`formats/zarr/codecs/javascript/`) stay as each family's other candidate: the decode pool measures both (`decoder-choice.js`, below) and keeps the faster on the machine, some engines' JIT running them faster than WebAssembly. Measured on one core against them (Node/V8, real float32 chunks re-encoded by numcodecs, byte-identical): Zstandard 345 → 717 MB/s, LZ4 664 → 2,695, Blosc-LZ4 841 → 1,224, Blosc bit shuffle 236 → 1,108, blosclz 710 → 968, Blosc-zlib 614–648 → 638–685. Each opens a file as the **array store** of `formats/store.js`: a store's tree is read as it is navigated into, a group at a time (`group.load(paths)`: the groups on the way and the node at each path to a depth, each object read once however many readings ask, a tree level's reads in flight together), and navigation is synchronous over what is loaded (a group navigated before it is loaded throws `StoreNotLoadedError`). Opening an HDF5 file reads its root and what its format opens with (`ANNDATA_OPENING_PATHS` in `data/anndata/elements.js`: every element's own node -- the group to depth 2: each dataframe column, each `obsm`, `varm`, `obsp`, `varp`, `layers` and `uns` entry, `raw` likewise -- and every expression matrix and `uns['cellucid']` whole; a layout view's `openingPaths`), named to `openHdf5` so they load with the root and before any opening step runs; what lies inside an element (a categorical's codes and categories, a nullable column's values and mask, a graph's arrays, nested `uns` entries) loads when a reading needs it (about half of a CELLxGENE file's objects at open: HLCA's 259). A Zarr store's metadata is described in one place (consolidated) or read in one batch, so it is loaded whole at open. Array values are read on demand, fetching only the chunks a selection touches, in parallel, each decoded once (a per-file decoded-chunk cache serves neighbouring reads); of an unfiltered HDF5 chunk (anndata's default sparse `.h5ad`), or a Zarr chunk stored as its values' bytes (no compressor or filter, C order: what `cellucid serve` and `cellucid_serve` answer), only the run of bytes a selection needs is fetched (`partialChunks`; a few columns of some rows of such a chunk read those rows' bytes), so 1,000 cells of a 1M-cell CSR `.h5ad` over 30 ms/5 MB/s read 4.5 MB in 12 s, not 121 MB in 25 s. A one-dimensional column of variable-length strings (an obs index, a text column) reads the rows asked for alone (`readRows`): of an HDF5 column the references of the chunks they lie in, then of each global heap collection those rows point into the prefix up to the furthest of their objects (the chunk's other references, and its neighbours', name the objects before it), a collection let go once the rows pointing into it are decoded; of a Zarr column the chunks they lie in, held as stored (bytes and an offset an item) with only the rows asked for decoded. A share of such a column decodes and holds its own strings, not every string of the chunks it touches. A categorical obs column whose dictionary has more categories than there are open cells, and more than one request of the link reads (its stored elements past the link's `maxRequestBytes`; a smaller dictionary is read whole, every category listed, the unused ones too), is coded by the categories its open cells hold (`AnnDataReader._categoryCoding`, `_readPresentCategorical`): its codes are read, then the labels and `uns['<key>_colors']` entries at the stored positions they name alone (`readCategoriesAt`, the dictionary gathered by those positions through `gatheredArray`) -- unless the rows read are enough to name a category in every chunk of a chunked dictionary (more than c · ln c rows for c chunks), when the dictionary and its colours are read whole beside the codes, in one plan, and the present ones picked (no round trip waits for the codes; no byte or request more than reading the dictionary); the field's categories are those in stored order, its codes positions among them, and the worker hands the page `present: {count, stored, colors}` (the adapter's `_adoptPresentCategories`, `presentCategoryPalette`: each category's colour is the file's or Cellucid's at its stored position). The coding is decided on the field's first read and kept: a clear of the reader's cache keeps it, an opening grown or cut in place carries it (`_takeCarried`), so added rows are coded with every category found plus the ones they add, in stored order, and the page grows the field by label (`category-growth.js`); the page marks such a field `discovered` (a saved category it lacks is one no open cell holds, `pending-categories.js`) with `storedCategoryCount` for the legend's line. |
| AnnData | `data/anndata/` | **One** AnnData reader for both formats: the element rules of the current layout (`elements.js`, `dataframe.js`, `sparse.js`), files from anndata 0.6/0.7 and SeuratDisk presented in the current layout (`older-layouts.js`), slabbed reads with read-ahead and exact conversions (`values.js`), gene columns and cell rows per layout (`matrix-access.js`), and opening plus error explanations (`open.js`, `views/`, `transcode/`). |
| Worker | `data/anndata-reader-worker.js`, `data/anndata-reader-client.js` | The reader's methods over `postMessage`, results transferred; live `activity` and `transfer` events. |

**Genes and cells by layout.** Only what one gene (or cell) needs is read:

| Matrix | A gene (column) | A cell (row) |
|---|---|---|
| CSC | its own entries: indices and values fetched together | the cells' entries of every gene, read as a cell subset (each gene's entries once, keeping the cells'; a priority-ordered file opened as a prefix reads each gene's prefix alone), never an index of every stored row |
| CSR, over a link as fast as a disk (a local file, the same machine or network) | indexed by gene once in the background, at open, indices and values into a column store (in memory, or on disk when it does not fit); a gene asked for meanwhile is answered from the same pass, and every later one is a copy out of the store | its own entries |
| CSR, over a slower link | indexed when the first gene is shown, from the **indices** alone (kept in memory when they fit the `columnIndex` budget; otherwise a second pass fills a column store, values and all, on disk when the site's storage quota has room); a gene's values are then fetched as the ranges its entries occupy; once the ranges already read and those asked for together would cost as much as every value (a request costs a round trip of the measured link), every value is read once, after which every gene is answered from memory (never more than twice the cheaper choice in hindsight) | its own entries |
| Dense, chunked in strips of few genes | the chunks of its strip (kept for its neighbours) | the chunks of its row strip |
| Dense, in one piece or in wide chunks | indexed by gene once, streaming | its row |

A CSC matrix's column pointers are read while the dataset opens, beside the
obs and var indexes, so the first gene shown is one round trip, not two.

**Gene indexes per matrix.** Switching the expression matrix (X, a layer,
`raw.X`, in the data mapping) keeps each built gene index while every kept
index fits the `columnIndex` memory budget, the most recently read first:
switching back to a matrix already indexed reads nothing. An index still
being built when another matrix is chosen stops (the genes waiting on it are
still answered). A sparse index keeps 8 bytes an entry (a uint32 cell and a
float32 value: 414 MB for Suo's 51.7 million entries); building it keeps the
first pass's indices for the second in the smallest type that holds a gene
position (2 bytes for up to 65,536 genes, 103 MB for Suo instead of 207 MB),
and the values a waiting gene made the first pass read, so no slab is
decoded twice.

**What the first gene costs, said before it is read.** How a matrix will be
indexed is decided from its metadata before anything is read
(`sparseIndexPlan`, the store factory's `plan`; `geneReadCost` in
`data/anndata/matrix-access.js`, which a subsampled opening weighs too), and
`AnnDataReader.geneReadPlan(element)` reports it: whether the matrix is
stored cell by cell, the stored bytes the first gene reads, the index it
builds and where that is kept (memory, disk, or nowhere, with the reason),
and the measured rate. The gene's loading notification states it in one
line, with the reasons as fact rows in its details (`geneReadPlanFacts`,
`data/gene-read-plan.js`; the cause from the layout through one helper,
`geneLayoutReason`/`geneLayoutFacts`, which the auto-fit notice and the
index pass's read activity use too), and has Cancel: a build asked for by
genes stops when the last of them withdraws — one round trip of the
measured link later, so a gene picked instead of the withdrawn one joins
the pass under way with its reads in flight — and its queued reads are
never sent. A stopped pass keeps the slabs it counted (`_pass1`: counts,
kept indices, row pointers, values read whole for genes), and the next gene
goes on from the first slab not counted; pointer windows end at fixed
multiples, so the slabs are the ones a pass from row 0 makes (pass 2 finds
them by number). A failed pass starts over. A matrix no store can hold is
scanned once per batch of waiting genes, each scan its own read activity
naming why (`onUnindexedScan`). The kept rows' row pointers a subsample's
first gene needs are read as one batch through the store's decoded chunks,
as a read activity of their own (`_readyRows`), and only the rows the gene
is read for: a gene sample's first cells (`readRowPointers` given `rows`,
then `rowsMatrix`: the sample's rows as a CSR matrix of their own;
`_openRowsMatrix`), never every open cell's (an 11.4M-cell CELLxGENE
atlas's 2,011-cell sample read its whole `indptr` before, 31 MB in 4,379
requests: the writer lays the pointer chunks between the matrix's, so none
merge); a stored matrix opened whole reads its first `rows + 1` pointers.
The rows' pointers are read in one pass over the rows ascending
(`CellRows.ascending`, made once), pointer ranges up to three pointers
apart read and held as one (a range of its own holds as many numbers): a
50 % uniform share's millions of pointer ranges are a few, and no run of
the rows in open order is placed in them until something reads the rows in
open order (`RowPointers.runIndexFor`). A share's pass in store order
(rows spread through the store, `_inStoreOrder`) answers in the open rows'
own file order (`CellRows.byRow`, made once, its rows the rows ascending),
and the gene plan of such a pass makes no matrix of the rows in open order
(`_segmentProgress`); the kept values and the share of index chunks they
lie in are one walk forward (`keptValueStats`), and pricing a read lists no
piece (`planBoxes` `pieces: false`). A 4M-cell CSR `.h5ad`'s 50 % uniform
share's first gene at 20 ms / 50 MB/s took 5.06 s against every cell's
4.83 s before these, 4.69 s after. Spread rows to the last open row whose
own pass is more work than every cell's, for the same bytes (they lie in
every chunk), read through **every stored row's index restricted to them**
(`_restrictsWhole`): every cell's pass over the stored matrix as stored (the
same reads and decoded chunks), kept in this browser as every cell's index;
no kept row's pointer is read (`_readyGeneRows`) and no entry gathered. A gene is placed from it on the open rows (`restrictedAccess`:
each value that is not zero by the open cells' lookup; once built, from the
gene's entries alone). The work is weighed at rates measured as passes run
on the source (`PassRates`: the index's own work per entry, and the view's
-- its kept rows' pointers, rows put in file order, entries gathered -- per
entry a pass over it indexes), carried to each opening of the same store;
before any is measured, a share keeping at least half of the entries reads
the whole. On a link faster than decoding the 50 % share's own pass took
about 1.4x every cell's first gene; through the whole it is every cell's.
Pointers read once are kept together and copied from by every later read,
piece by piece (a range only part of which is held reads only what it
lacks; a grown sample reads only its added rows'), and a read of fewer cells
than are open by their rows (a cell export, cells past the sample) reads
only theirs (`_namedRows`). Over a link, values stored in
chunks read whole (compressed or filtered) are read beside every slab's
indices while a gene waits, and kept, when they fit beside the kept indices
(`sparseIndexPlan` `buysValues`): the first gene was priced with them
(`firstGeneReadCost`), and no later gene reads a value again; a raw array's
values are read per gene by range. A raise reads the pointer
chunks its added rows share with the open ones from those decoded chunks,
as every gathered read of a subsample does (`readKeptRowPointers`,
`inducedGraphLines`, a modality's cell map: `chunks`). A subsample's neighbour
graph (the induced subgraph of the open cells) streams into its edges a slab
at a time, as a whole graph does; it is never held whole. It is read in the
stored graph's row order, never in the order the cells were opened (a
uniform share's open order is scattered through the file: read that way,
every slab spanned the whole graph and decoded its chunks again — a
2.5M-cell uniform share of a 5M-cell `.h5ad` took 19.9 s for its edges
against 2.0 s for every cell's), by the cheaper of two reads
(`AnnDataReader._graphReadPlan`, priced as a matrix's rows are,
`firstGeneReadCost`): the whole graph streamed once from the first open row
to the last, keeping the open rows, whenever that reads no more bytes than
the open rows alone (they lie in every chunk: a uniform share of a graph in
chunks of thousands of cells), else the open rows' pointers and entries
only (`readKeptRowPointers`, `gatheredArray`, a block share). Each line's
column is its row's rank among the open rows in stored order, so the lines
ascend and the edges are built and checked exactly as a whole graph's, then
placed on the open cells (`ConnectivityLines` `relabel`,
`relabelledEdges`): strengths and partners placed in one walk when the four
edge arrays held at once take no more than every stored line's edges take
to build (`wholeBuildBytes`), else strengths placed and let go before the
partners, about 1.5 times the edges' bytes, counted in `'graph build'`. Each
edge is counted at its smaller cell as its line arrives (the builder's
`cellOffsets`, while the walk waits for the next slab), each placement walk
is a typed loop of its own (no closure an edge, no count after the last
line), the open rows' stored order is the open cells' own (`CellRows`
`byRow`), and a stream's row pointers are asked for before the lines'
order and ranks are made.
Over a link slower than a disk, slabs are read ahead toward filling every
request the link runs at once (`_inFlightEntries`), within the decoded
chunks' share of memory (`_graphReadAhead`); from a disk, one ahead. A MuData modality's graph is read the same way in its
own rows (`placedGraphLines`, `readPlacedGraphRows`). Measured (Node, local
5M-cell CSR `.h5ad`, k = 15): the 2.5M uniform share's edges 19.9 → 2.2 s
(every cell 1.9–2.5 s), 0.5M 3.7 → 1.2 s; edges bit-identical.

**The gene sample prices the first gene.** The map opens every cell the
device holds (`autoFitCount`: the working set, and the gene indexes every
open cell holds); a cell-major matrix's genes are read for a **gene sample**
(`geneSampleOf` in `data/anndata/subsample-plan.js`): the first `cells` open
cells in their order, every cell of a file in its own order opened in its
`blocks` subsample's order (`cellOrder`) so the sample is a nested prefix,
and a share spread through such a file (`uniform`, `stratified`, `spatial`)
opened in that same order with its cells filtered (`shareInCellOrder`,
`priorityOrder` of the cell order's `blocks` over the share's cells): its
sample is its cells in every cell's first blocks, a subset of the same
reads. In its own priority order each of its first cells lay in a chunk of
its own (1 % of a 10M-cell CSR `.h5ad` over an emulated S3 link: 71 cells
in 10.7 s against every cell's 1,572,326 in 5.5 s). The share's cells are
unchanged; the order's rung is recorded (`cellOrder`), a request reading
the same cells again names it (`order`, `requestOfOpenCells`), the kept
gene index files it with the share's specification (`openOrderSpecOf`,
`orderKey`), and `blockRung` looks for a kept index under the share's order
at each rung (`keptOrder`). A `blocks` share keeps its own order; a share
of a priority-ordered file opens in file order (its cells by file row: a
uniform sample of them at every prefix, in the file's first rows); a share
whose every open cell is read anyway keeps its priority order, as every
cell then keeps file order. A share's sample (and the page's top-up after
it) never makes more bytes or requests than every cell's sample in the
same time (`wholeSampleCost`, `geneTimeFit`'s `ceiling`,
`geneSampleCost`'s `geneTime.spentRequests`): on Tabula Sapiens a 1 %
share in the cell order, bounded by bytes alone, read 3.4 MB in 434
requests for 83 cells where every cell's sample made 351.
Over a measured link `firstGeneReadCost` prices the index pass for the
sample's rows from the layout alone — the stored chunks those rows lie in
(each row taken to hold the mean stored values) for the indices (and
values, unless read per gene) and the pointers, the requests their runs
make, and the chunk-index reads that find them (`chunkIndexReads` of the
store contract: an HDF5 version 1 chunk B-tree's nodes not read yet, each
counted once, a level a round trip; a sharded Zarr array's shard indexes),
in dependent steps (`steps`: the pointers' index levels, their chunks, then
the entries' index levels and chunks) — and `linkReadMs` turns that into
time (a round trip per wave of requests the link runs at once in each step,
or a round trip a step and the bytes at its rate). The round trip is the
link's under load (`link.js` `loadedLatencyMs`, `loadedRoundTripMs`: the
median of the latest requests' times less their bytes' transfer at their
share of the bandwidth, of requests made on a connection already set up:
one that opened its connection waited for its handshakes too, and an
opening's first requests priced an export's round trip at 883 ms against
365 ms as it ran), not its shortest: over S3 a request took 0.35–0.5 s
with 42 in flight against 230 ms idle, and an export priced at the shortest
said half its time. A request costs a plan no less than the host takes one
under load (`loadedRequestMs`: with every slot busy, the time between
completions less their bytes' transfer at the fastest delivery seen; part
of `requestCostMs`, so of `requestMs`, every merge gap and `linkReadMs`):
a path may complete only 70–88 small requests a second whatever the
connections (CELLxGENE's S3 mirrors from Europe with 42 in flight: 14 ms a
request where a round trip over the connections said 6 ms), so a read's
floor is the longer of its requests at that rate and its bytes at the
bandwidth, and plans merge through the gaps that pay at that price. The
other side of that floor: on a host whose throughput is per connection
(S3: the full HLCA's 703.5 MB of edges took 172 s as some 2,300 requests
of 64–256 KB and 252–296 s as 320–600 requests of several MB, about 7 of
42 connections busy), one long request moves only what its connection
does, so a plan's long ranges are split over every connection
(`byte-reader.js` `requestPieceBytes`: the plan's bytes over
`shapingConcurrency`, never below one connection's round trip of bytes,
`LinkMeter.connectionRoundTripBytes`, nor below `mergeGapBytes` or a block;
at most `maxRequestBytes`). One connection is measured only by a
transfer-bound request that was alone on the link from start to end, and
plans are split only when that delivery is below 0.8 of the most the link
has delivered (`CONNECTION_BOUND_SHARE`), shown by four transfer-bound
request ends, never one (`LINK_EVIDENCE`: one delivery timed over a stretch
little longer than the clock's tick can be any rate, and on a local web host
it split every later plan of the page: a 1 % share's `.h5ad` export of a 10M
file went from 473 to 3,135 requests, 6.2 to 18.2 s): a request sharing the link
delivers its share of a busy link, and taking that for a connection's rate
split each next plan finer (a 10M-cell export over an emulated 30 MB/s
HTTP/2 link went from 255 to 7,556 requests of 16–40 KB, 18 % slower, the
same bytes); on a host whose one connection carries the link (most web
hosts and CDNs) nothing is split. The round trip (the shortest of the latest requests
that were mostly a round trip) shapes reads (merge gaps, request
sizes); every price of a read's time (the gene sample, Read for all cells,
an export, a progressive first view, an analysis run's read-ahead) takes
the loaded one. The rows are priced as
runs (`PricedRows`): a `blocks` order a block at a time
(`subsampleBlockRuns`), a subsample's rows ordered once; `geneTimeFit` takes
the most cells within half of `GENE_READ_TARGET_MS` (`FIRST_GENE_SHARE`;
20 s, a response-time target, not a platform limit; the fit's `budgetMs`
is the whole 20 s), `geneIndexFit` the most whose gene index fits the
`columnIndex` share of what the open cells leave, and the fewer is the
sample, recorded as `genes` `{cells, of, fit}` (`fit` one of the two, with
`limit` `geneTime` or `columnIndex`). The reader indexes a matrix in
segments up to `genes.cells` only (`prefix-column-index.js`, `readRows`): a
gene's column is NaN on the other open cells, a sparse column has no entry
there, and `getMatrixCellSubset` reads a cell past it by its row (its own
values, never the index's absence of it as a zero).
`AnnDataReader.growGeneSample` reads it for more cells in segments doubling
from the sample (`growthSegmentBounds`: the first as many cells as the
sample, about the first gene's time; `geneSampleCost` prices them first,
and when the gene index of every open cell is kept nowhere here names the
most cells it is kept for, `index.fits`, which the legend then offers),
keeping what is built when stopped; a request names the sample to reopen
with (`genes`: a count, or `'all'`) and every cell's order (`order`).
After the first gene of a sample the gene time sized, the page reads genes
for more cells in what is left of `budgetMs` (`DataState.topUpGeneSample`,
`app/state/managers/field/gene-sample.js`): each step through
`readGenesForMoreCells`, sized from the reads' measured times (`cellsWithin`:
a fixed time and a time a cell from the last two reads of different sizes;
after the first alone, only a step at least as large, at its rate; the
first gene's cells answered by a kept gene index, `fit.kept`, not counted),
capped by `geneSampleCost`'s `index.fits`, never by a price, and for a
share by the rule its first gene was fitted by: its sample transfers no more
over the gene time than every cell's would (`geneSampleCost`'s `geneTime`
option: `wholeSampleCost`, every cell's sample at the gene time over the
link measured after the first gene, less what the sample transferred and
requested so far; `ceiling` its bytes and requests, `target` the cells
lowered to them), since a scattered share's steps each
fetch nearly every stored chunk again; on a later visit the cells whose
gene index was kept (`fit.kept`) count as read at the rate the first gene
measured on the others, so the gene time buys no more cells in all than it
would have read; a gene asked for never waits behind a step
(`yieldGeneSampleTopUp`, called as every gene load starts: the step under
way stops at once, the reader trimming its unbuilt segments on the abort,
and the steps go on after the gene within the time left); a step adding
under a quarter of the cells read is not made; a person's Read for all cells
takes over (the top-up stops, keeping what it read); the sample keeps its
`fit` (still the gene time's: the reader's own, `growGeneSample`'s `topUp`;
a person's count has none). Each opening of the same bytes carries the
cells the earlier ones opened, with their cell order and sample as they
were when replaced (`AnnDataReader._geneSamplesSeen`,
`noteGeneSampleSeen`): the same cells opened again with a request naming
no sample — a share lowered and raised back — take that order (a gene
index kept under it is found), a person's count or every cell read as they
were, and a gene-time sample fitted again over the link as measured now but
floored at the cells it read and never every cell when it read fewer
(`geneSampleOf`'s `seen`/`seenFitted`); a share's `wholeSampleCost` is every
cell's fitted sample as it was read. Fitted afresh alone, a 50M-cell store's
sample of 4.49M became every cell (markers 15.4 → 62.3 s, a reader worker of
8.7 GB), a 2M-cell `.h5ad` over an emulated 40 MB/s link 1,690,543 → 352,680
cells in another block order; reused alone, E-ANND-5's opening sample of 68
cells (priced on the link's first small reads, its top-up cut short by the
share change) came back as 68 where the measured link reads 13,290. A progressive opening's steps carry the
fitted sample with its fit (`cell-subset-policy.js` `geneSampleRequestOf`:
`genes: {cells, fit}`; a count alone is a person's or a session's choice,
`fit` null), so an atlas grown to every cell is topped up too; before, its
last step's request named the count only and the top-up never ran after a
progressive opening. A price made at the opening rests on the
link's first, small reads: on Census S3 the opening's resolve measured
0.96–1.1 MB/s where the first gene then measured 11 MB/s.
`geneSampleOf` answers before any gene or row pointer is read (the open
cells' view of a CSR matrix is described by its shape), so the page knows
the sample at open and every analysis asks it of the source first.

**A subsample's genes fill in as its cells are read.** A subsample opens its
cells in an order whose every run of first cells is a sample of them (its
priority order, the file's cell order, or file order: `subsample-plan.js`). Its cell-major matrix is indexed in segments of that
order (`data/anndata/prefix-column-index.js`): the first about a
thirty-second of the cells (at least one slab of the index's work), each next
one as large as all before it, one after another, sharing one memory budget
(each segment's store placed with what the built ones left,
`withMemoryBudget`). A segment's `indices` and `data` are the subsample's
gathered entries cut to its rows, read as one gather of the stored arrays
(`gatheredArray` composes the two sets of runs, `composedRuns`), so a gene of
k cells holds memory of their entries and the decoded-chunk budget, never of
the file: read through the subsample's estimated chunks instead, each slab
was tens of thousands of reads at once, and a gene of 3,401 cells of HLCA
(4 billion stored entries) held gigabytes in the reader worker for 20 MB of
entries (`tests/subsample-gene-read-memory.test.mjs`). Their builds report into one read activity
(`wholeGeneIndexActivity`): it expects the stored bytes of every unbuilt
segment as it starts and counts each segment's progress as its share of
their rows, so its bytes and time left are the whole gene's. `AnnDataReader.getMatrixColumn(…, {onPartial})` answers
a gene after each segment with its values on the cells indexed so far and NaN
on the rest; the reader worker posts those as `{event: 'partial', column}`
before the reply (`getMatrixColumn`'s `progressive` argument), as the cells
read so far in their sparse form when that is smaller (below), so a partial
answer moves a gene's entries, not a whole column (40 MB at 10M cells, 200
MB at 50M); the client hands them to the callers that asked (`onPartial`),
and the page draws them as colours only (`DataState.showPartialVarField`:
the gene is not loaded until it is whole, so filters, analyses and sessions
never see partial values; Cancel puts the colours of the field shown back).
`geneReadPlan` says what remains to read and when the first colours come
(`progressive`).

**Two forms of a gene.** A view colours by a gene's dense column (float32,
one value per cell), made in the reader worker and transferred, so the page
never expands one. Beside it the worker makes, in one pass, its sparse form
when that is smaller (`data/column-forms.js`: `{length, indices, values,
readRows}`, cells from `readRows` on not read, a gene sample's NaN tail),
and the client notes it beside the column (`noteColumnForm`): the colour
kernels (`state/managers/continuous-colors.js`) and the range scan read the
form, a word fill and the gene's entries, instead of two passes over every
cell on the main thread (at 10M cells, a gene shown 66–70 ms → 6–10 ms;
a partial answer 146 → 6 ms), and a gene grown or cut in place keeps its
form (`resizedForm`). A column is never written in place, so its form stays
exact while it lives. A computation reads the same gene as its non-zero entries,
`{length, indices: Uint32Array, values: Float32Array}`
(`getMatrixColumnsSparse`; `DataLayer.readGeneColumn` →
`DataState.readVarSparseColumn` → the var loader's `readSparseColumn` →
`getGeneSparseValues`): the dense column's values bit for bit (a duplicate
coordinate summed in stored order, a zero left out), about 45 times fewer
bytes for a median Suo gene (0.2 % non-zero), and nothing installed in field
state. The analysis kernels take either form (`compute/gene-statistics.js`).
Every typed array a worker reply carries owns its whole buffer (a view is
copied out first), so replies are transferred, never cloned.

**Genes an analysis computes never reach the page.** Marker discovery,
differential expression, the gene-moments pass, the dot plot and UCell
compute in the analysis workers, and the reader delivers their genes to those
workers directly (`data/gene-delivery.js`): a run opens a gene route
(`app/analysis/data/gene-route.js`), one `MessageChannel` per analysis worker,
the reader holding one end of each and each worker the other; the page asks
for a delivery (`DataLayer.deliverGeneColumns` /
`deliverGeneCellSubset` → `DataState.deliverVarSparseColumns` → the var
loader's `deliverSparseColumns` → `AnnDataReaderClient.deliverMatrixColumnsSparse`)
and hears back only the genes' bytes, then hands each worker a task naming
its delivery. Relayed through the page, a run's genes had cost its main
thread more than the messages: WebKit counts every array a worker sends the
page as the page's own allocation, so a run over every gene drove the page's
garbage collector, whole-heap collections included, while it ran. A gene only
the page holds (a user-defined copy) is read here and sent with its task.

**Some genes in some cells.** Differential expression between two named
groups reads one gene-major block of the groups' cells
(`getMatrixCellSubset`): each gene's entries, with each entry's position in
the cells asked for (ascending) and its value. It is read by whichever route
decodes or copies fewer bytes (`readCellSubset` in `matrix-access.js`): the
cells' rows (CSR or dense: every chunk the rows touch, since a chunk is
decoded whole; nearby rows are read as one range, entries filtered into the
block without a per-row expansion) or the genes' columns (a CSC matrix's
ranges, or a gene index's entries in memory plus what is left of building
it). Few cells of a large matrix go by cell, many cells, or any once the gene
index is built, by gene. Suo, all 8,192 genes, Node, local `.h5ad`, index
built: the two largest cell types (103,736 cells) 8.8 s → 0.69 s (by gene),
two types of 16,376 cells 0.92 → 0.30 s, one type of 61 cells 16 → 15 ms (by
cell); the run's peak memory 3.6 → 2.0 GB.

**A request is priced.** Every read of the cells' rows is priced at no less
than a request of the store's link (`requestPriceOf` in `matrix-access.js`,
from `formats/io/link.js`: its merge gap, and a round trip's worth of bytes,
or the merge gap at every request the link runs at once until it is
measured), in three places: rows nearer than the merge gap are read as one
range even where an array is read in part (an uncompressed or contiguous
`.h5ad`, whose chunks do not join rows), a read in flight counts as at least
a request's worth of entries (`pipelined`), so no more requests are in flight
than the in-flight entries hold at that price, and the route's cost
(`readCost`) adds a request's worth for every range. Unpriced, a two-group
subset of a 10-million-cell uncompressed `.h5ad` (two categories, a million
cells scattered through the file, about 10 entries each) put a read per cell
in flight at once: a million requests, each with its promises, plan and
pieces on the reader's heap, which grew to its 4 GB limit within ten seconds
and took the tab down; at a million cells the same read held 1 GB of heap for
5.5 s. Priced, the same 10M-cell differential expression finishes in 2.8 s
with the reader's heap at 4 MB (Chromium, a synthetic uncompressed CSR
`.h5ad` of 200 genes and 10 values a cell; at 1M cells 5.8 → 0.8 s). CSC gene ranges and dense rows are
priced the same way. Nor does an
analysis ask for every gene at once: two groups' blocks are read as many
genes at a time as the run's memory holds ({doc}`11_analysis_architecture`).

A CSC matrix read by gene reads consecutive genes as one range of the
indices and one of the values (every gene of Suo, 51.8 million entries: 26
reads of each array instead of 8,192), several ranges in flight while the
entries already read are kept, and keeps each gene's entries in one pass
without a branch on which cells the subset holds; only a gene stored with a
duplicate coordinate or unsorted rows goes through `sparseColumnOf` first.
Suo, the two largest cell types (118,905 cells), every gene, local `.h5ad`,
Chromium: 0.82 s → 0.36 s for the read alone (decoding the 1,580 chunks it
touches takes 121 ms on the decode pool; the rest is the reader thread
keeping 51.8 million entries and copying them out of their chunks), 1.16 →
0.62 s in the benchmark's sequence of steps. Each range of genes is then kept
on a decode worker (`keepSubsetRun` in `subset-entries.js`, sent with
`runOnDecodeWorker`: the range's arrays are transferred, the cells' positions
are sent to each worker once as a context), while the reader thread reads and
decodes the next: the rebuilt Suo (5,103 genes), the same cells, local file
in Chromium, 0.31 s → 0.21 s, against 0.15 s for reading the whole of `X`
(2.1× → 1.4×), byte for byte the block the reader thread keeps.

**Decoding keeps up.** An `.h5ad` read fetches its chunks in groups of
neighbours, as many bytes in flight (stored and decoded) as twice what keeps
the link full, twice the largest group a read has fetched, and twice what
keeps every decompression slot busy with the largest chunk seen
(`formats/hdf5/context.js`), so the next group is fetched while one decodes:
every value of Suo's expression matrix 210 → 175 ms (121 ms to decode the
chunks all at once), its neighbour graph 110 → 100 ms (54 ms).

**Graphs and embeddings arrive finished.** For a neighbour graph the reader
worker checks the `obsp` matrix as its lines stream in (bounds, order, a
diagonal left out, exact symmetry of topology and weight) and builds its
forward adjacency: each edge once, at its smaller cell, with its partner and
its `float32` relative strength, the forward offsets, and the graph's largest
weight and degree (`ConnectivityEdges`, `data/connectivity-edges.js`). That is
8 bytes an edge and 4 a cell, the only copy: its three arrays are transferred
to the page, never copied, and the build is part of the working set
(`'graph build'`, `memory-budget.js`) while it runs, so every thread's caches
make room for it. A k = 15 graph of 50M cells (350M edges) holds 3.0 GB and
peaks at 3.2 GB while it is built (float32 weights; float64 weights are held
whole until the strongest is known, 8 bytes an edge more). Which edges a
pane shows is derived from it on the page; nothing is stored in draw order.
The k-NN tool reads the rows it expands through from the file, by cell, in
the reader worker (`getConnectivityRows`, through the chunks that hold them;
a frontier whose rows lie in every chunk streams the graph instead, as a
subsample's edges do), so no thread holds every cell's neighbour list. While
a view shows the graph's edges the page holds them anyway, and a frontier is
then answered from them (`forwardEdgeRows`: two passes over every edge, the
same neighbours) whenever that is faster than the file by the rates both
ways measured on this page (`knnRowReader` in `app/main.js`; 20M edges, a
65,536-cell frontier: 0.18 s, against 1.7 s from a 10M-cell file). Each upper entry of a line is
an edge kept with its partner and exact weight; each lower entry must match
the next due edge of its partner's line, so one pass in storage order checks
the whole graph.
An index built ahead of need (the cells of a CSC matrix, asked for before)
does not take turns with an extraction someone waits for: Suo's graph after
1,000 cell rows, Chromium, 0.80 → 0.47 s. The page checks only the
payload's shape against the graph summary. The adapter keeps only the
neighbour graph asked for last: it is the payload the viewer shows, held by
reference, so keeping it costs nothing, and a graph no longer shown is
released and read again if it is shown again. For an embedding the worker reads the selected `obsm`
columns at their stored precision, checks every value, normalizes them once
with Cellucid's rule (or fits a file's own coordinates to the view) and sends
`float32` positions with their display transform and per-axis bounds; nothing
on the page scans them again. Vector fields arrive scaled by their
embedding's transform, with their largest magnitude, read concurrently with
the positions.

**Memory budgets.** Every cache that decides what stays resident (the gene
indexes with their row pointers, positions of layouts no view shows,
inactive fields, the analysis, the adapter, the latent space held for
outlier quantiles, decoded chunks, fetched bytes, decode workers, undo
history) takes a disjoint share of what the drawn cells leave of one
allowance from `data/memory-budget.js` (see {ref}`memory-budgets`): half the
estimated device, divided with the other Cellucid pages of the origin, the
shares summing to exactly one. Each registers what it keeps, in whichever
thread it lives, and is trimmed when the allowance or the working set
changes. A refused allocation in any thread halves the allowance everywhere.

**Requests from the page.** The reader client sends the gene reads asked for
in one task as one request per matrix and form (dense or sparse), each gene
read once however many callers asked; a second caller of the same gene gets
its own copy. Every read takes an `AbortSignal`: aborting rejects the caller
at once and sends the worker `{cancel: id}`, whose reads then stop (a request
the network has not sent is never sent; a gene index other reads still need
goes on). A batched gene request is cancelled once every caller in it has
cancelled. The same rule holds below it: **a shared read is cancelled only
when every reader waiting on it withdrew** (`utils/abort.js`). Reads
merged into one request (neighbouring ranges of a file, a Zarr store's
chunks), a read two openings of the same file both wait on, and a read the
page's byte cache hands to another reader each stop for one reader at once
-- that reader's wait is rejected -- and the request itself is aborted only
when no reader is left waiting. The worker's open reply and its transfer events say how many
reads the source answers at once (unbounded for a local file, the link's
background concurrency for a web host), which the page's download scheduler
follows for long runs (`setSourceRequestConcurrency`).

**String obs columns.** A string or categorical column whose categories are
not listed in the file's metadata (a direct AnnData string column) is listed
with its categories unknown; they are found when the column is first read,
so opening never reads a column's values.

**The link.** Nothing asks where a file is, except that a picked Zarr folder
is this machine's disk (see *A Zarr folder from disk* below). Every byte
reader (and the request layer of a Zarr folder on a web host) measures its
link as it reads
(`formats/io/link.js`) over its latest requests (`RECENT_REQUESTS`, 16),
since a connection changes while it is read: the round trip is the shortest
of the latest requests whose time was mostly a round trip, the bandwidth the
fastest of the latest deliveries, over one request or over the stretch of
requests in flight it ended (the windowed minimum and maximum TCP BBR models
a path by). A link that gets faster is followed at its first faster request,
one that gets slower once a window of requests saw nothing faster, and one
odd request (a stall) moves neither. A delivery bounded by its round trip
(a request's time under three round trips: small reads, a subsample's kept
rows) is only a lower bound of the bandwidth: it may raise it, never lower
it (BBR's application-limited samples).
Each request is counted half a tick of the clock longer than measured
(`clockTickMs`: `performance.now()` ticks in whole milliseconds in WebKit's
and Firefox's workers, 0.1 ms in Chromium's), so the shortest of many lies in
the middle of where the true one does and a request to a server on the same
machine, ending on the tick it began on, does not make the round trip 0.
Until it is measured, a local file is assumed to be a disk and a web host
nothing. From these:

- each plan of reads asked for together fetches the gaps between its ranges
  smallest first, as many as make it quickest as the link prices it: its
  requests' round trips (`requests / concurrency` × latency) plus its
  bytes over the bandwidth. Over a host whose protocol is not known yet (six
  requests at a time over HTTP/1.1, a hundred over HTTP/2) the plan taken is
  the one nearest the best under either: merging for six fetched 2.6× the
  bytes over HTTP/2 (99 genes of Suo: 9.3 MB for 3.6), merging for a hundred
  costs round trips over HTTP/1.1. 99 genes of Suo's `.h5ad` at 80 ms /
  50 MB/s, the protocol not yet known: over HTTP/2 316 → 186 ms (11.3 →
  4.6 MB; one round trip and the bytes take 151 ms), over HTTP/1.1 688 →
  780 ms (the quickest plan of single-range requests, simulated, takes
  621 ms: six connections, not bytes, set that floor). The request size is one round trip's
  worth of the whole link (bandwidth × latency, never less than the
  source's prior): the bandwidth measured is shared by every request in
  flight, so a request that long already keeps the link busy for its own
  round trip, and a longer one only serialises its bytes on one connection
  (24 connections at 50 MB/s and 150 ms asked 180 MB requests before; 7.5 MB
  now);
- requests waiting for a slot are merged with their waiting neighbours,
  whichever reads asked for them, while their round trips would take longer
  than their bytes and a gap is no wider than the waiting requests' own plan
  would fetch, and a request that fetched gaps is kept (within the metadata
  cache budget) so the next wave of a scattered read finds its bytes;
- the sparse cost model prices a request at one measured round trip;
- **concurrency** starts at the host's HTTP/2 stream limit: 100, the
  minimum RFC 9113 recommends and the value GitHub announces. A browser keeps
  only six connections to an HTTP/1.1 host, and the page often cannot see
  which protocol a host speaks (GitHub sends no `Timing-Allow-Origin`, so
  `nextHopProtocol` is empty), so the link learns it: Resource Timing settles
  it where the host shows it; otherwise, over HTTP/1.1 the seventh request of
  a burst completes a whole round trip after the first, so seven completing
  within half a round trip prove multiplexing, and six completing together
  followed by half a round trip of silence, twice, mean HTTP/1.1 (the reader
  drops to six). A request belongs to the burst it began in even when its
  round trip is timed from later (a read that first asked the page's byte
  cache, `range-cache.js`): timing those from the host's answer had put
  every completion outside its burst, so no host was ever taken for
  HTTP/1.1 this way. On a host proven to multiplex the concurrency keeps doubling
  while the completion rate rises by at least half what doubling could give.
  Until a host is proven to multiplex, reads are split as for six requests
  (`shapingConcurrency`) and their gaps priced under both six and the host's
  streams (above), requests beyond six start only while those
  in flight carry less than two round trips' worth of bytes, and background
  work (a gene index built ahead of need) never holds more than five, so on an
  HTTP/1.1 host a read a view waits for still gets a connection. A local file
  needs none of this: it is read ahead from the start, whatever one small read
  measures. Zarr folder on a GitHub-like host (simulated, 80 ms / 50 MB/s):
  open 307 → 203 ms, 57 obs fields 475 → 199 ms, 99 genes 332 → 196 ms;
  `.h5ad` 1,000 cell rows 5.98 → 1.51 s;
- a chunk B-tree of an HDF5 file (`.h5ad`, `.h5mu`, `.loom`, a 10x `.h5`) is
  read from its root as reads need it: the levels below the root are read
  with the table only while a level is one wave of the link (its nodes
  transfer in less than a round trip and go out at once), every other node
  -- internal or leaf -- when a read first needs the chunks below it, a level
  a round trip. What a dataset stores is estimated from the nodes read so far
  (every key holds the size of the first chunk below it), so pricing a read
  costs no request. Before, every leaf's parent was read with the table: on
  an 11.4M-cell CELLxGENE atlas (`X/indices` in 2.26M chunks) that was 1,099
  requests and 13.5 s over S3 before the first view; now its opening reads
  58 index and structure requests (5.2 s), 1,281 → 240 requests for the
  whole opening. Once a tree's leaves are all known, the rest of them follow
  in the background when they transfer in less than a round trip, so a later
  gene is one round trip, not two. The reads
  asked for in one task (the genes of one call, the arrays of one sparse
  gene) wait for the leaves of all of them together (`chunk-index.js`,
  `indexBurst`), so their chunks go out in one task and the byte reader
  plans them as one set instead of wave by wave as each read's own leaves
  land (100 Loom genes over 150 ms / 5 MB/s: 11.1 → 7.8 s for every cell,
  2.1 → 1.6 s at 1 %). A sharded Zarr v3 array does the same with its shard
  indexes (`formats/zarr/array.js`): a read asks for every shard index it
  needs at once and for the inner chunks of all its shards once they are
  in. A read that
  goes on where the last one of the same array ended (rows in file order: a
  cell-major matrix read slab by slab, a subsample's kept rows wave by wave)
  has the leaves after it read ahead in the background, as many as it
  needed, asked for as its own leaves land, in the same task as its chunks:
  a leaf sits among the chunks it indexes, so the next one comes in their
  request when the gap is worth it, and each next read waits one round
  trip, for its chunks, not two; a read elsewhere reads nothing ahead. Background work waits for a
  leaf already being read ahead instead of asking for it again. What an
  array stores (`storedBytes`, which cost models and progress read before a
  read) is known without reading a leaf: a leaf's parent holds the key of
  the leaf's first chunk, its stored size included, so a chunk whose leaf is
  not read yet counts that size (exact for unfiltered chunks; within a few
  per cent for gzip ones, whose chunks of one array compress alike), and a
  chunk whose leaf is read counts its own. The 2.3M-cell lung atlas on
  CELLxGENE (`X` CSR, 470,000 gzip chunks of `indices`, 235,000 of `data`)
  read all 12,382 leaves of both, 2 KB each and 850 KB apart, one request
  each, before its first gene could start: 11,706 requests in the first 10
  minutes in Chromium; now none. A sharded Zarr v3 array on a web host is
  answered the same way: the shard indexes in hand count exactly, one more
  is read, and the other shards' inner chunks count the mean stored size of
  the indexes in hand (anndata's automatic shards hold about 11.5M values,
  so a matrix of 4 billion values has some 350 shards: one request, not 350);
  a folder or ZIP, which knows its shards' sizes, reads every index;
- **background work is coalesced too**: the reads of a gene index built
  over a subsample's kept rows (background work the first gene waits for)
  are merged across their gaps by the same price as any read, and while
  requests wait for a slot the gaps are priced at twice the measured
  bandwidth (`LinkMeter.pricingBytesPerMs`): a link measured only while its
  connections were full of small requests measured its round trips, not its
  bandwidth, and priced at that no gap would ever pay. On an atlas-shaped
  file over 150 ms / 50 MB/s (24 connections), 20 % blocks: first gene
  10.6 → 6.9 s, 1,617 → 1,013 requests (more bytes, fewer round trips);
- **mirrors**: a public Amazon S3 object (an unsigned URL of a regional S3
  endpoint: the CELLxGENE Census bucket every CELLxGENE link reads, or any
  S3 link) is read through every address S3 answers it at
  (`formats/io/s3.js`: virtual-hosted and path-style, each also on the
  dual-stack endpoint and on the region's dashed endpoint `s3-<region>`, and
  virtual-hosted on the global endpoint; one `ETag`, CORS open on each,
  measured on the Census bucket in us-west-2). S3 speaks only HTTP/1.1, and
  a browser keeps six connections to one host, so a read of many small
  requests waited on six round trips at a time; spread over seven hosts it
  runs on 42 (a Region without a dashed endpoint does not resolve it: that
  mirror is dropped at its second failure). The `.h5ad` export of 21,285
  scattered cells of CELLxGENE's HLCA core (Node, the app's reader, six
  connections a host): 4 hosts 493 s, 7 hosts 383 s; read as one stream of
  parts (`readMatrixRows`) 271 s, every connection busy throughout
  (`httpByteSource` `mirrors`, and the same for a
  Zarr folder on S3, `httpStore` `mirrors`). Each request goes to the
  address with the fewest requests asked of it; a mirror that fails, or
  answers another status, size or `ETag` than the file's own would, is
  dropped for the rest of the read; one that stalls or cannot be reached is
  dropped at its second failure in a row; one that says it is busy (a 503
  every address of a bucket gives alike) is kept (`formats/io/mirrors.js`,
  shared by both). Either way the request is asked of the file's own URL at
  once, so a mirror adds connections and never a failure. The lung
  atlas's first gene over a link of about 2 MB/s and 0.3–1 s to S3 (Node,
  six connections a host as a browser keeps): in its first 90 s, 34 → 135 MB
  with mirrors alone; with both changes, in its first 120 s, 2,500 requests
  of index nodes (12 MB, none of the matrix) → 940 requests and 200 MB of
  the matrix, the link's full rate;
- one `abort` listener per `AbortSignal` (`utils/abort.js`), however
  many reads wait on it: thousands of reads under one signal (the rows of a
  large selection) register in constant time each, and an abort still
  rejects every one at once;
- **read-ahead**: a link that delivers at least 64 MB/s (`READ_AHEAD_BYTES_PER_MS`,
  about the rate at which the reader indexes stored bytes, measured) is not
  what a read waits for, so the file is read as a local one: the gene index
  starts when the file opens, and every value is read with it. Below that,
  and until the link is measured, nothing is read before a view needs it.

Measured with `scripts/benchmark-anndata.mjs` (`--protocol http/1.1|h2`
simulates a host that answers six requests at once or all of them; Node, so
decoding shares the reader thread) on Suo (562k cells × 8,192 genes), and with
the same calls against `cellucid serve` on the same machine. Before → now:

| Store and link | Step | Before | Now |
|---|---|---|---|
| CSR `.h5ad`, `cellucid serve` on localhost, a gene chosen 3 s after opening | first gene / 100 genes | 1,322 / 574 ms | 3 / 22 ms |
| CSR `.h5ad`, 20 ms / 100 MB/s | 99 genes / 1,000 cells | 682 / 4,448 ms | 13 / 1,902 ms |
| CSR `.h5ad`, 80 ms / 50 MB/s | 1,000 cells | 16.7 s | 7.1 s |
| Cellucid-ready Zarr folder, 80 ms / 50 MB/s, HTTP/2 host | 99 genes / 1,000 cells | 1.13 / 21.7 s | 0.29 / 3.7 s |
| Cellucid-ready Zarr folder, 150 ms / 5 MB/s, HTTP/2 host | 99 genes / 1,000 cells | 2.3 / 43.6 s | 1.2 / 25.5 s |
| Cellucid-ready `.h5ad`, 80 ms / 50 MB/s, HTTP/2 host | 99 genes | 1.08 s | 0.33 s |

Over an HTTP/1.1 host a Zarr folder stays where six requests at a time put
it (1,000 cells: 22 s at 80 ms): each chunk is its own object, and nothing
merges two objects into one request.

**A Zarr folder on a web host that hides `_` and `.` paths.** Jekyll (GitHub
Pages without a `.nojekyll` file) serves no path segment that starts with `_`
or `.`. A Zarr v3 store still opens there (its consolidated metadata is in the
root `zarr.json`), but every chunk under anndata's `obs/_index` and
`var/_index` answers 404, which a Zarr reader must take as an absent chunk:
cell and gene names would read as empty, with no error. So when a store is
opened from consolidated metadata on a web host, the metadata document of the
first `_`/`.` node on each such path (`obs/_index/zarr.json`, or
`.zarray`/`.zgroup` in Zarr v2) is fetched once, by the first read of an array
under it and beside that read's own chunk requests, which wait for its answer
too: a host that serves everything costs one parallel request per such node
read and no round trip in sequence, and opening or listing a store fetches
nothing for it. A host that answers it
as absent is refused, naming the fix: an empty `.nojekyll` file at the root of
the site, or the store opened through `raw.githubusercontent.com` or a `gh:`
link. (A Zarr v2 store on such a host already fails at its `.zmetadata`.)

**A Zarr folder from disk.** A picked or dropped `.zarr` folder is read file
by file (`formats/zarr/stores/directory.js`), each chunk a
`File.slice().arrayBuffer()`:

- its link is this machine's disk, a local file's policy, and is **never
  measured**: a folder read is one file each, whose time is the file system's
  cost of opening a file plus however long the reader thread is busy before it
  takes the bytes, not the disk's speed. Measured from its first small files
  (a 21 KB `.zmetadata` in 1 ms), a disk delivering gigabytes a second fell
  below the 64 MB/s read-ahead line in most opens, and a CSC matrix's cell rows
  were indexed without their values (1,000 rows of Suo, Node: 2.3–2.7 s
  instead of 1.4–1.6 s);
- 16 file reads are in flight at once, a view's ahead of background work's
  (an index built ahead of need never takes the last slot): enough that chunk
  files arrive faster than the decode pool consumes them (Chromium, Apple M1
  Pro SSD: 1,580 chunk files of 63 KB arrive at 640 MB/s eight at a time,
  800 MB/s at 16, while the decode pool consumes about 670 MB/s of stored
  bytes); more only buffers files ahead of the decoders;
- the inner chunks of a shard (Zarr v3 `sharding_indexed`) a read needs are
  read in one slice per run of neighbours, up to 1 MiB: merging saves a read's
  fixed cost, a longer slice only delays decoding its first chunk (1,000 cell
  rows of a sharded Suo: 744 ms chunk by chunk, 519 ms merged to 1 MiB; the
  graph's 1 MB inner chunks 554 ms, against 687 ms merged to 32 MiB);
- copying a chunk into a read walks the fewest axes the layout allows (axes
  both sides lay out as one are merged, axes of extent 1 dropped): an
  embedding chunked by rows is one block move per chunk instead of one per row
  (24-column `scanvi`, 562k cells: 69 → 59 ms in Chromium; 127 → 109 ms in
  Node);
- the queue of waiting reads (and of waiting chunk decodes) is kept per
  activity, so starting the next one never scans every waiting read: 25,000
  queued background reads cost 2.4 s of queue work before, 14 ms now.

A `.zarr.zip` is read through the byte reader as one file. Opening reads its
last 256 KiB (the footer and central directory), and every entry lying in
those bytes is served from them: an archive no longer than that is read once,
in one request, locally or over HTTP (the 39 KB demo archive
`synthetic-cell-types-2d.zarr.zip`, opened and every array read: 30 requests
and 56 KB before, now 1 request and 39 KB). Every entry is checked against its
CRC-32, now computed sixteen bytes at a time from 32-bit words (2.2 against
1.4 GB/s on the reader thread).

Measured in a Chromium worker (Playwright, headless; the same calls as
`scripts/benchmark-anndata.mjs`, best of three, Apple M1 Pro, files in the
page cache) on Suo written by `cellucid.prepare()` (Zarr v2 zlib 6, v3 gzip 6,
v3 sharded, and v2 with 4 KB chunks, 25,260 files) and by anndata (Zarr v2
Blosc-LZ4, v3 Zstandard, CSR float64), whole runs in milliseconds, before →
now:

| Store | Total | 1,000 cells | Graph | Embedding (24 columns) |
|---|---|---|---|---|
| `.h5ad` (for comparison) | 1,424 → 1,195 | 490 → 485 | 657 → 425 | 59 → 63 |
| Cellucid-ready `.zarr.zip` | 1,576 → 1,165 | 596 → 526 | 611 → 343 | 138 → 78 |
| Cellucid-ready folder, Zarr v2 | 1,236 → 1,082 | 551 → 581 | 440 → 255 | 69 → 59 |
| Cellucid-ready folder, Zarr v3 | 1,197 → 1,027 | 530 → 514 | 410 → 261 | 71 → 52 |
| Cellucid-ready folder, Zarr v3 sharded | 1,308 → 996 | 590 → 488 | 420 → 286 | 71 → 51 |
| Cellucid-ready folder, 4 KB chunks | 4,366 → 3,955 | 3,528 → 3,349 | 482 → 272 | 76 → 57 |
| anndata `.zarr.zip` | 2,220 → 1,780 | 246 → 174 | 551 → 282 | 174 → 141 |
| anndata folder, Zarr v2 | 1,894 → 1,670 | 138 → 138 | 474 → 257 | 85 → 105 |
| anndata folder, Zarr v3 | 2,460 → 2,149 | 229 → 232 | 555 → 357 | 110 → 123 |

(The graph and the `.h5ad` row also move with changes to the reader outside
`formats/zarr` made in the same release.) A Cellucid-ready folder is now
faster in total than the same data as `.h5ad` or `.zarr.zip` in a browser,
and within 1.2× of the faster of the two at every step that takes more than
a few milliseconds, except the first obs field (39 ms against 19 ms from the
`.h5ad`; about 34 ms from the `.zarr.zip`, mostly starting the decode
workers). A folder of 25,260 files of 4 KB stays bound by the file count:
1,000 cells take 3.3 s against 0.5 s, the cost of opening every file plus
decoding chunks too small to send to a decode worker on the reader thread.
An anndata store's matrix is CSR, so its gene index is built in the
background from the moment the folder opens, as for a local `.h5ad` (it was
not before, whenever the folder's link had been mis-measured): the first obs
field read meanwhile takes 39–55 ms instead of 19 ms, while the index
competes for the reader thread.

Node is not a browser here: a `File` from `fs.openAsBlob` opens and reads a
file on the reading thread, about 0.35 ms per file, so in Node a folder's
1,580 chunk files of one cell-row read cost about half a second more than the
same bytes in one `.zarr.zip` (1,000 cells: 1.34 s against 0.72 s); in
Chromium the browser reads files off the worker's thread.

Reads of many genes or cells at once (`getMatrixColumns`, `readMatrixRows`, or
many single reads issued together, as analysis does) go out as one set of
coalesced, parallel requests. `readMatrixRows` (a cell export) reads every
cell's row in one call: the rows' pointers once, then parts of about a slab
of entries in the cells' order, each part asked for before the one before it
is handed out (`partial` messages from the reader worker), so the link stays
busy across parts and the caller holds one part at a time; a read that
arrived frees its request's share of the in-flight budget at once
(`matrix-access.js` `pipelined`), so a slow read at the head never holds the
reads behind it below the link's concurrency. An index built ahead of need is background
work: its requests wait while any read a view waits for is queued and never
take the last request slot, so colouring by a field or drawing an embedding
is not held up by it; once a gene waits for the index, the index is no
longer background.

**Values.** Gene values, obs columns and nullable integers are float32, as
the views draw them: every stored type converts as numpy's
`astype(float32)` does (64-bit integers to the nearest float32, a float64
beyond the float32 range to an infinity), so one outlier never makes a
column unreadable. Boolean matrices read as 0 and 1. Embeddings and vector
fields are coordinates: a value float32 cannot hold is refused there.
Integer category labels beyond 2^53 are kept as their exact decimal text.

**Older AnnData releases** are read through views (`older-layouts.js`), so the
element rules exist once. **A neighbour graph** is not read at open: a direct
dataset records that it has one and counts its edges when they are first
read.

**Seeing what is read.** Every reader operation (opening, a field, a batch of
genes, an embedding, a gene index, the graph) runs under a `ReadActivity`.
It knows how much is left: the byte reader counts each request when it is
planned and when it arrives, and an operation that streams a whole array
says its stored size up front, so a running operation reports `x of y`
bytes, its rate and the time left (from steps where it counts them, as an
index pass does). The worker reports each live (at most every 250 ms) and
when it ends; `data/read-activity-log.js` keeps them for the open dataset
and writes a running line every 2 s for anything longer than a second, and
one line when it ends:

```text
[cellucid:read] Sparse X, indexed by gene — 98 MB of 215 MB (45%), 112 MB/s, ~2 s left, 41 requests, 5210 chunks decoded to 512 MB, 1.1 s (reading)
[cellucid:read] 2 genes of X — 137 KB read, 9 requests, 1 from cache, 16 chunks decoded to 192 KB, 6 ms (done)
```

`app/read-progress-notifications.js` shows the reads a person waits for that
nothing else shows as one progress notification, in plain words ("Reading
field cell_type and the neighbour graph · 4.1 MB of 9.3 MB · ~3 s left"), on
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

**Verification.** `tests/hdf5-reader-conformance.test.mjs` reads every object of
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
`scripts/mtx-corpus/`). Benchmarks: `npm run bench:hdf5` (the
HDF5 reader against h5wasm, on local files), `npm run bench:anndata`
(`scripts/benchmark-anndata.mjs`: the AnnData reader against an earlier
revision's, on generated `.h5ad` and Zarr v2/v3 files, local and over a
simulated link), `npm run bench:formats` (`scripts/benchmark-formats.mjs`:
every other format, Loom, a Cell Ranger `outs/` folder, a MatrixMarket
folder, Seurat and SingleCellExperiment `.rds`/`.rda`, `.h5mu` and MuData Zarr,
against the same data as a Cellucid-ready `.h5ad`: first view, first gene,
peak memory, and a stream format's pass against plain inflate of its own
bytes) and `scripts/benchmark-real-links.mjs` (the whole app on real hosts,
against an earlier revision; `--corpus` takes the real-link corpus below). A
corpus of real public links of every source and format
(`scripts/efficiency-floor/public-links.json`: CELLxGENE from 2.5k cells to
the 2.3M-cell lung atlas, Zenodo, GEO `.h5ad`, `.h5ad.gz`, 10x `.h5`,
MatrixMarket, `.rds.gz` and text matrices, the Single Cell Expression
Atlas, Figshare, Hugging Face, GitHub, Dropbox, Google Drive, public S3, the
sample catalog) is opened by `scripts/efficiency-floor/public-links.mjs`
through the app's own path in Node over the real network (a browser's six
connections to an HTTP/1.1 host and its name): each step's time, bytes,
requests and share of the file (open and first view, a colouring field,
the gene list, genes), the host's behaviour (HTTP version, ranges, CORS and
the headers it exposes, content coding, validators, rate limits), and every
request a second visit sends, which should be the version check alone. The R object corpus is written again with the
current CRAN and Bioconductor packages once a month
(`.github/workflows/r-corpus-drift.yml`, `scripts/rds-corpus/drift.sh`).

---

(web-host-requests)=
## Requests to a web host: retries, stalls and failures

Every request to a web host (a byte range of an `.h5ad` or `.zarr.zip`, a
Zarr store's objects, the size probe, catalog files and provider APIs) goes
through `fetchOnce` and `withRetries` in `formats/io/byte-source.js`
(`fetchWhenNotBusy` in `data/data-source.js` for catalog files and provider
APIs). None is left waiting for ever:

- **Transient failures are asked again**, with waits that are the link's
  own measure, as TCP retransmits: the first is the link's retransmission
  timeout (`ResponseClock.timeoutMs`, below: at least one second), each next
  one twice the last (RFC 6298 5.5), each stretched by up to a half by a
  draw from the request's own URL and range (`retryWaitMs`, `spreadOf`:
  requests a busy host turned away together come back apart, and the same
  request always waits the same), never shorter than a `Retry-After` the
  page may read. A request is reported as failed once it was asked again at
  least three times and its waits add up to 100 seconds (`RETRY_AT_LEAST`,
  `RETRY_FOR_MS`: RFC 1122 4.2.3.5's R1 and R2, the least a TCP connection
  is given before it is dropped). That patience is for a read its host
  already delivered to (`withRetries` `delivered`: a file whose ranges, a
  folder whose objects arrived, a download resumed part way): a busy or
  failing answer to a request whose host has delivered nothing to it yet
  (a catalog file, a provider's API, a file's or a folder's opening) is
  reported once asked again three times (R1 alone, a few seconds), so a
  server answering with errors is said to the person, with **Try again**,
  instead of holding a startup step past its bound; a request that gets no
  answer at all keeps the whole patience either way. Transient: HTTP 408,
  429 and 5xx (but 501 and 505), a
  403 that says it is a rate limit (Google's APIs name `rateLimitExceeded`
  in its JSON; GitHub's API sets `X-RateLimit-Remaining: 0`,
  `isLimitAnswer`), and a connection that fails or drops. GitHub Pages
  answers a burst of requests with an occasional 503 that carries no CORS
  headers, which every browser reports to the page as a failed connection,
  so a failed connection is transient too.
- **A wait the host names is waited for exactly, and is no failure.** A
  `Retry-After` of a second or more that the page may read (the host lists
  it in `Access-Control-Expose-Headers`), or a request budget with nothing
  left and its reset stated (below), is the host saying when: the request is
  asked again then, however often the host says so, and it counts toward
  neither the repeats nor the 100 seconds. A `503` with `Retry-After` is the
  server preparing what was asked for (`cellucid serve` laying a matrix out
  gene-major answers so, with how far it is in a `text/plain` line, `Retry-
  After` one second until it has a pace): it holds that request only, never
  the server's other keys, and the read's activity carries the server's line
  (`ReadActivity.waiting`), which the gene's notification shows ("Reading
  gene CD3E: X is being laid out gene-major for the viewer (37% done, about
  12 s left).") and the opening progress shows as **Waiting for the host**
  (the transfer log's `waiting`), never as a stall or a retry.
- **A host's request budget paces every request to it**
  (`formats/io/hosts.js`, one book of hosts per realm). Where a host states
  its budget on its answers -- `X-RateLimit-Limit/Remaining/Reset` (Zenodo:
  133 a minute, every range read counted, the reset in epoch seconds;
  GitHub's API: 60 an hour without a sign-in), the IETF `RateLimit-*`
  fields, or their structured form (`RateLimit: "api";r=499;t=227` with
  `RateLimit-Policy`, Hugging Face's) -- requests are sent while what remains
  exceeds the requests on their way, and wait for the next answer, or with
  none on its way for the reset, when it does not. A `429` with a readable
  `Retry-After` holds every request to the host until then, and the host is
  asked fewer at once after it (halved from what was on its way, one more per
  round it answers). A block -- a host that answered this site (most of the
  minute's requests) and then sends nothing a page may read to two rounds in
  a row (Zenodo's and
  Google's "unusual traffic" pages, which Zenodo sent after about 550 range
  reads in a few minutes, within its stated 133 a minute) -- sets the host's
  **pace** (`learnPace`; a `429` is an answer, handled as above, and sets
  none): requests to it are spaced so a minute holds at most
  half the requests sent in the minute before the block, a further block
  halving again, whatever its headers allow; `requestCostMs` prices each
  request at no less, so plans read through larger gaps in fewer requests.
  The pace is told to the other realms, kept per origin by the page
  (`host-link-memory.js` `rememberHostPace`, `cellucid:host-pace:v1`) and
  seeded into the page's and the reader worker's books when the host is
  opened again (`seedPace`), so a later visit never sends the burst that met
  the block. It recovers slowly (owner, 2026-10-30): a host that blocked
  once and has refused nothing for a day (`PACE_RECOVERY_AFTER_MS`) starts
  each later visit a step looser (`visitPace`: a quarter more requests a
  minute, once a visit); a block met at a loosened pace snaps the kept pace
  back to the one the first block taught and marks the host as one that
  blocks again, which then stays paced, as does one that blocked twice. What one realm learns of a host (its budget, a wait it
  named) is told to the others of this site -- the page, its reader workers,
  other tabs -- on a `BroadcastChannel`, since the host counts them all. A
  request whose answer is optional (GitHub's commit of a file, without which
  the file is read at its branch) is never held for a spent budget. Measured
  on Zenodo, 150 range reads of one file at once, from a window with 133
  left: all 150 read in 62 s, the 35 refused (sent before the first answer
  told the budget) waited for the reset, none failed.
- **A host's budget prices every request to it** (`hosts.js`
  `requestCostMs`, `link.js` `requestMs`): its interval, the window over the
  limit (Zenodo: 60 s / 133, 0.45 s; the window learnt from the longest
  reset stated with a whole window), in the share of the window already
  spent -- nearly nothing while most of it remains, so early reads stay
  exact, and the whole interval once nothing remains for what is asked or a
  request has waited for the budget -- however fast the host answers. The
  byte reader's plans (`planMergeGap`), the requests a waiting one takes
  along, the link's merge gap and the cost models (`linkReadMs`: auto-fit,
  the progressive plan) price a request at no less, so as the budget runs
  down gaps worth that long at the link's rate are read through; and the
  metadata blocks around a cached read (an HDF5 header, heap or chunk-index
  node) are read with it, half of what the link moves in that interval each
  side, so the file's structure
  comes in a few large requests. A host that states no budget (or no
  window) is priced as before.
- **A host that says it is busy is asked less at once.** Each 429 (or a
  403 that names a rate limit, `isLimitAnswer`) halves the requests sent to
  that host at once, from those on their way, for every reader and store of
  this site that asks it (`formats/io/hosts.js`, `slowDown`), and they grow
  back by one for every round of requests it answers (additive increase,
  multiplicative decrease). A 5xx is retried, not counted as a limit.
- **A stalled request is aborted and asked again at once**, except one to
  a notebook's kernel (`formats/io/host-fetch.js`): the kernel answers its
  requests one after another, and none while a cell runs, so asking again
  would only queue a copy; its silence is reported ("still waiting") and
  the request waited for. The link's
  `ResponseClock` (`formats/io/link.js`) applies TCP's retransmission
  timeout (RFC 6298) to the time requests wait for their answers: the timeout
  is the smoothed wait plus four times its variation, never under one second,
  and three seconds before the first answer, which waits for the connection
  to be made too (`FIRST_ANSWER_TIMEOUT_MS`, RFC 6298 5.7 and RFC 1122
  4.2.3.1: a Dropbox share link's first byte takes about 2 s, S3's up to
  2.4 s, measured), doubled by every stall and computed afresh by the next
  answer. A request is stalled when neither it nor any
  other request of the link, nor any request to its host anywhere in this
  site, has received anything (an answer, or bytes of a body) for a whole
  timeout: a browser keeps six connections to an HTTP/1.1 host for every
  reader and tab of the site, so a request queued behind another reader's
  (a catalog listing several files of one host, GEO's three MatrixMarket
  files, an old and a new reader worker of one file) is waiting, not stalled,
  while the host delivers (`hosts.js` counts each host's deliveries, and
  tells the other realms of the site at most every half second, half the
  shortest timeout). A range of known length is also
  given the time its bytes take at the measured bandwidth, since a host that
  buffers a range (a proxy, a CDN fetching the object from its origin)
  answers only once it holds all of it. Bodies are read as streams, so a
  body that stops arriving is caught as well as an answer that never comes.
  The file's head, read beside the probe of its size, seeds the timeout with
  how long the host took to answer it (and a file whose size its provider
  gave without a probe, from Zenodo or Hugging Face, from how long the host
  took to answer the provider's own request, `hosts.js` `rememberAnswer`): a host slower than a second to answer
  (S3's colder storage classes take 1–3 s for a first byte, measured on the
  CELLxGENE bucket) no longer has every first read taken for stalled and
  asked twice.
- **A file that changed between its description and its first read is
  refused**: the `ETag` or, from a host that sends none, the `Last-Modified`
  the probe saw (`describeRemoteFile` carries both; with the size they are
  the file's version, `writeVersion`, unless its provider names its content)
  must be what every answer says, a whole download's included
  (`downloadRemoteFile`). The probe's `HEAD` passes the HTTP cache by
  (`no-store`), so the version is the host's at this opening; it goes out
  beside the head read, and a file whose two answers name two versions
  (their `ETag`s, weak marks aside, or `Last-Modified` dates, where both
  state one) was replaced between them: it is probed again, and refused as
  `changed` when the next two disagree too. Nothing keeps a
  description past the opening that made it: a link opened again asks its
  provider again (`link-source.js`), and a dataset catalog describes a file
  at each opening that reads it.
- **An answer of other bytes is refused, never read.** Where the host lets
  the page see `Content-Range` (`cellucid serve`, `cellucid_serve()`, most
  object stores), a 206 must name exactly the range asked for (another range:
  `range-ignored`) of a file of the size it was opened at (another size:
  `changed`); and an answer whose `ETag` differs from the file's, or, from a
  host that sends no `ETag`, whose `Last-Modified` differs, means the file
  changed while it was open (`changed`: "open it again to read the new
  version"), whether the answer came from the file's URL, a redirect or a
  mirror. A Zarr folder's ranges (a shard's index and inner chunks) are held
  to the same `Content-Range` check.
- **A failure that will stay is reported at once** (`worthAskingAgain`).
  The signs that a host will answer later keep the whole schedule: a
  `Retry-After`, a status that says the host or its gateway is busy (408,
  429, 502, 503, 504: S3's `503 SlowDown` has the same body and no
  `Retry-After` every time), a connection dropped in the middle of a body, a
  stall. Any other error status without `Retry-After` that repeats exactly
  (the same status and `Content-Length`: Zenodo answers a missing record
  with the same 500 every time) is reported after its second answer, and a request that got no
  answer at all (no connection, or a host that does not let this site read
  it) after it was asked once more; the momentary 503 GitHub Pages sends
  without CORS headers is answered on that second request.
- **A host that stops answering is left alone, not asked by every read**
  (`hosts.js`). A request that gets no answer a page may read is a lost
  request. A file's range reads count their losses by address (its URL;
  Drive's file and its details are two addresses, as Google holds back only
  the file), requests on their way together as one round: a second lost
  round in a row (the retries of the first) leaves the address unasked for
  RFC 6298's retransmission timeout, one second, doubled for every further
  lost round up to its 60 s maximum. Meanwhile its reads fail at once,
  unsent (`HostNotAnswering`, `sent: false`), and then one request asks
  again; the first readable answer ends it. A read of a file whose ranges
  (or a folder whose objects) were delivered before waits such a refusal
  out instead of failing (`withRetries` `delivered`,
  `refusesForNow`): asked again after a back-off from the link's
  retransmission timeout doubled to 60 s, uncounted, its activity told
  "`<host>` stopped answering this site, which it answered before … Asking
  again in N s; what was read is kept.", until it is answered or no longer
  wanted; a long read (a graph's edges, genes for all cells) then goes on
  where it was (CELLxGENE's S3 refused this machine part way through
  HLCA's graph twice on 2026-10-30, and the whole read failed). A file that
  never delivered a range fails as before (CORS, a host that cannot serve
  ranges). A provider's
  own requests (a record, a file's details) are asked as the person asked,
  each twice. Zenodo answered this machine with an HTML "unusual traffic"
  page (HTTP 403) without CORS headers from 2026-10-28 on, after many
  requests from one network: 40 reads of one file then cost 8 requests to
  the refusing host instead of 74. The failure says what the page knows:
  from a host that answered this site before (`answeredBefore`), it is not
  CORS, so "`<host>` stopped answering this site, which it answered
  before: it may be refusing requests from this network for now (hosts do
  so for a while after many requests from one network), or failing, or this
  computer lost its connection. Try again later."; a provider's first
  request without an answer says "`<provider>` (`<host>`) did not answer
  this site. It lets other websites read it, so this computer's connection
  may be down, `<provider>` may be failing, or it may be refusing requests
  from this network for now … Try again later."
- **A failure that outlasts its repeats names the request**, with the
  host's own words when it sent a `text/plain` line:
  `https://…/suo.zarr/obsm/X_umap_3d/3.0: the host could not answer (HTTP
  503). Asked 8 times.` A failure that is not transient (404 on a file, a
  short or changed answer, a range ignored) is reported at once.
- **A sign-in never follows a redirect to another host.** A redirected URL
  is reused until it fails or its stated lifetime ends; reused on another
  host than the one asked, it is asked with `Range` alone, so a token given
  for Google Drive or Microsoft Graph never reaches a CDN it redirected to
  (a browser drops `Authorization` on a cross-origin redirect it follows
  itself).
- **A whole download** (a host that serves no part of a file, or a
  compressed copy read from its start) is asked as any request is, its body
  watched for silence as it is written to storage (`fetchOnce`'s `stream`):
  a busy answer is asked again, a host that stops sending fails it as
  stalled instead of holding the opening, and a body cut short goes on from
  the byte it stopped at (`Range: bytes=<received>-`, under the `ETag` and
  date it started with) where the host serves ranges, so a drop on an
  unsteady link never starts a large compressed copy over; only from a host
  that ignores ranges or content-codes the file is it downloaded again from
  its start, for as long as a request would be asked again. Each repeat is
  said in the Transfers log.
- Each repeat is reported to the transfer log while it happens (`retrying`,
  naming the URL; the console's `[cellucid:transfer]` warning; the Dataset
  box's one "Retrying a read" notification, closed when bytes arrive again),
  and bytes of a
  request still under way are reported as they arrive (`arriving`, at most
  every 250 ms), so a read of one large request shows progress before it
  completes. The bytes also count into the read's activity as they arrive
  (the `[cellucid:read]` progress, its rate and time left), and an attempt
  that fails has its bytes taken back, so what is fetched never passes what
  was planned before it has arrived.
- Cancellation is unchanged: the caller's signal aborts the request in
  flight and any wait to ask again at once, and nothing is asked after it. A
  range read's planned bytes are counted once, however often its request was
  made.

Measured in WebKit, Chromium and Firefox against Suo served with injected
faults (the app's own requests: a 503 without CORS headers, a 429 with
`Retry-After: 1`, a request never answered, a body that stops after its first
bytes, a dropped connection): each recovers and the dataset opens; a request
that is never answered on any attempt ends the startup step with the error
above. `tests/http-request-resilience-contract.test.mjs` holds the same cases
against a fake `fetch`, and `tests/host-budget-contract.test.mjs` every
budget header form, a 429 and a 403 rate limit, a served store's 503 "being
prepared", the per-host stall watch, a redirect's headers and the probe's two
versions.

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
column when it is coloured, a gene's values when it is chosen, and a
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

### 1) HTTP cache (browser)

Every file a dataset catalog names — `datasets.json` and a default session —
is fetched by `fetchCatalogFile` with `{ cache: 'no-cache' }`
(`CATALOG_FILE_CACHE` in `data/data-source.js`); a caller cannot choose
another mode. The browser keeps its copy but asks the host before every use,
and an unchanged file costs a conditional request answered with 304.

Reusing a cached catalog without asking would be wrong, not merely stale: a
default session is pinned by the catalog's digest, and hosts such as GitHub
Pages cache for minutes, so a new catalog paired with an old session would be
refused.

AnnData bytes are read with range requests. A Cellucid server answers them
with `Cache-Control: no-cache` and a strong `ETag`, and a range asked with
`If-Range` naming an older version is answered with the whole current file,
so bytes of a changed file are never mixed with old ones.

Every request into a Zarr store on a web host bypasses the HTTP cache
(`cache: 'no-store'`, `formats/zarr/stores/http.js`), as every range read of a
file does. WebKit (Safari, and Playwright's WebKit) never finished reading the
bodies of chunks it served from its HTTP cache to the AnnData reader worker:
Suo opened from the sample catalog loaded once, then every reload, with the
store unchanged, waited for those bodies until the startup step's deadline
(measured 2026-10-23 with a persistent WebKit profile and the store served
with GitHub Pages' `Cache-Control: max-age=600`: the default cache mode and
`no-cache` revalidation both hung on every reload; `no-store` opened it three
reloads in a row). A small store (the browser fixtures) does not trigger it,
so the suite asserts the cache mode (`http-request-resilience-contract`)
rather than reproducing the hang. Bypassing the cache also means a browser
never pairs a store's metadata with chunks it cached from an earlier version
of the store. A reload reads the store again (Suo's first view: 5.5 MB).

A query string cannot version a store instead: GitHub Pages' CDN ignores it
(a never-seen `?v=` on a store file is answered `x-cache: HIT`), and it would
break presigned bucket URLs, whose signature covers the query. For up to ten
minutes after a store is rebuilt, the CDN may still answer some of its
objects from before the rebuild; a chunk that does not fit the new metadata
is refused with a decode or size error naming its URL, never waited on.

Debug tip:
- DevTools → Network → “Disable cache” while you debug.

### 2) In-memory caches in the data layer

Examples:
- `DataSourceManager` caches dataset metadata per source.
- A catalog source keeps the `datasets.json` it read and each file's
  description; the reader keeps what it loaded of an open file's tree and its
  decoded chunks (`formats/store.js`).
- Opening another cell subset of the open file (the subsample control,
  `app/ui/modules/subsample-control.js`; a
  session's recorded cells) happens in the same reader worker
  (`AnnDataFileSource.reopenCells` → `AnnDataReaderClient.reopenCells`, a
  new adapter on the same client) and reads through the store already open
  (`anndata-reader-worker.js`, `sourceKey`): its structure, chunk indexes,
  cached metadata blocks, decoded chunks and measured link stay, and only
  the cells change, so a new subset costs no request before its own cells'
  bytes (a CSR gzip `.h5ad` over 80 ms: the 31 requests and 0.8 s of its
  structure; the lung atlas: its 130 requests and 10–30 s). The bytes are
  the same when the source names the same version: a local file by its
  name, size and modification time, a remote file by its identity's key and
  version (`ETag`, provider revision) or its URL, size and `ETag`; a picked
  folder (a Zarr store, a Cell Ranger `outs/` or MatrixMarket folder) or a
  set of files by every file's path and that key, so a file added, removed
  or rewritten is another source; a Zarr store folder on a web host read at
  a version naming its content (a GitHub folder at a commit: the provider's
  key and commit) by that key and version, whatever the host lets the page
  see of its headers (GitHub's raw host exposes none to a page, so its
  validators read as absent: every reading of a GitHub folder, a growth step
  included, was opened anew and the page stopped growing at its first view);
  any other store folder on a web host (a served store, a folder link) by its URL while
  the validators of its consolidated metadata (`ETag` and `Last-Modified`,
  each the host names) are the ones it was opened with, asked again with one
  `HEAD` (`httpStore.unchanged`). That requires the metadata to move with
  any rewrite of the store: every Zarr writer writes the consolidated
  metadata again when it writes the store, so its write time
  (`Last-Modified`), or a tag made from it (`cellucid serve`'s `ETag`: size,
  times, inode), changes. A host naming only a tag made from the bytes
  shows a rewrite only when the metadata's bytes change: one keeping every
  shape and attribute and changing only values is not seen. Another version, a folder whose host names no validator, or a
  set of files one of which is known only by its size opens anew. The new reader also takes
  over (`AnnDataReader.open`'s `carry`) what the open one read that does not
  depend on the cells — CSR row pointers, a CSC matrix's column pointers —
  and the gene index of the open cells both share as a prefix: a subsample's
  cells open in priority order (`subsample-plan.js`, `priorityOrder`), so a
  larger `uniform` or `blocks` subsample under the same seed begins with the
  smaller one, and its index (built in segments of that order,
  `anndata/prefix-column-index.js`) is extended for the new cells only; a
  smaller one cuts it and reads nothing. A failed or cancelled re-reading
  reads the cells open before again (carried, so it fetches nothing) and
  leaves the dataset on screen as it was.
- **The open cells grown or cut in place** (`app/open-cells-resize.js`).
  Every reading of the open dataset names its own cells (a cell subset
  request passed to `reopenCells`; the module-level request of
  `cell-subset-policy.js` is a link's, for the next opening only): Apply
  during a progressive opening's step, and the earlier cells read back
  after a cancelled one, never take or overwrite another's. Before
  reading, the work under way on the open cells ends on them
  (`data/open-cells-work.js`: field, gene, embedding and vector field
  loads, the gene list read after the first view, a data-mapping change,
  cancellable jobs such as an analysis run or outlier quantiles; a run's
  result arrives and is marked stale; and the operations that hold the open
  cells in use until they finish, `useOpenCells`: saving a session, the
  notebook's session request and the CSV and `.h5ad` cell exports); the
  change begins in the same turn as it sees no work under way
  (`beginCellsChangeOnceSettled`); from then until the commit,
  work about to start waits and starts on the cells after it (a field
  picked meanwhile loads on them). When the new opening's first cells are all of the open ones, or its cells
  are all among their first (`sharedOpenCells`, `anndata-file-source.js`:
  `reopenCells` answers `resized: {from, to}`), nothing shown is read
  again: the new adapter does not read the first view's embedding
  (`initialize({firstView: false})`), the reader keeps every embedding in the
  frame the cells were first shown in (`EmbeddingPayloads.adoptTransforms`)
  and the cell-order sums of the rows both openings share
  (`adoptCellOrderSums`, the session's cell check:
  {doc}`10_sessions_persistence_and_serialization`),
  and the page reads only rows `[from, to)` of what it holds — every
  resident layout (`loadCoordinateRows`), every loaded obs field and gene
  (the loaders' `readRows`; `BaseAnnDataAdapter`, "Rows of the open cells").
  The reader answers `rows` reading only what those rows need where the
  layout allows: an obs column or embedding's row range, a cell-major
  matrix's index segments the rows lie in, a priority-ordered file's
  gene-major `X` between the prefix-count levels around the rows
  (`_placedColumnRows`, `_prefixEntryRows`). Everything the commit needs
  (the identity and priority of the new count, the live view's layout, the
  length of every added row) is checked before it changes anything; an
  embedding whose added rows hold a value that cannot be drawn is let go,
  as an opening of those cells refuses it (when it is the one shown, Apply
  shows the cells as a new reading and a progressive opening stops at the
  cells open). Then, in one synchronous step, the dimension manager's
  layouts (`resizeCells`) and every per-cell structure of `DataState`
  (`state/managers/open-cells.js`, which tabulates each) grow or are cut,
  and the viewer takes the new live points with the camera kept
  (`resizeLiveCells`). A failure anywhere before it reads the earlier cells
  back with their own request. Kept views are retired and made again from
  their record (the renderer sizes them by the count); a selection under
  way stays when cells are added and is retired when cells are cut; the
  neighbour graph, the vector field overlay and outlier quantiles in use
  are read again on the new cells; analyses mark their results stale. What
  is derived from the file takes the added cells as opening them would: a
  user-defined category made from a column (a copy, edits, a consensus)
  labels them as it labels its other cells of the same source category,
  and a category or range highlight group takes those in its category or
  range that are visible; labels and groups made by hand stay as they were.
  Each per-cell structure is bit for bit what opening at the new count
  directly holds after the same actions
  (`tests/open-cells-resize-state.test.mjs`). A dataset read again with
  other cells gives its newly open cells the same (`cell-subset-carry.js`,
  `takeInNewCells`): a category or range group keeps its type and rule on
  both paths. The same step grows a
  progressive opening (`anndata/progressive-open.js`: a first view sized
  by one rule for every source, a picked file as a web host -- the
  measured link (a picked file's is its reads'), the stored bytes of the
  embedding chunks its cells lie in (`firstViewCost`), their decoding at
  the rate chunks are measured to decode in the reader (`decodingRate`;
  before any chunk has been decoded there, the first open cell's chunk is
  read to measure it when a read costs no round trip worth pricing, a
  read the first view makes anyway) and the page's measured cost of
  showing a cell, `data/opening-cost.js`, taken only when it arrives at
  least twice as soon as every cell; then doubling steps,
  `cellSubset.growth`). Every size the opening resolves -- the cells asked
  for, the first view, each step, and the cells asked for again at the
  last -- is selected through one subsample cache the reader holds for the
  opening and carries to each step (`subsample.js` `createSubsampleCache`):
  the cells (a `blocks` subsample's blocks, a stratum's cells) are ranked in
  one pass for the largest size, each smaller size is a prefix of that rank
  order and its priority order the prefix itself, and a priority-ordered
  file's uniform steps are prefixes of its computed row order, so a step
  costs its own cells, not the file's (10M cells in Node: a 50 % uniform
  growth's resolutions 1.3 s → 0.2 s; a priority-ordered file opened whole,
  its steps 570 → 49 ms). A subsample's cells lie across the whole file, so
  every step's added cells lie in chunks earlier steps decoded: a
  subsample's view reads every array decoded a chunk at a time through
  whole decoded chunks the store keeps (`anndata/decoded-chunks.js`,
  carried from each opening of the store to the next, let go with the
  store), within the one decoded-chunk budget the store's own cache shares
  (`formats/io/decoded-budget.js`: together they never hold more than the
  `decodedBlocks` share; chunks stay until something needs their room, the
  allowance or the working set changes, or the page is hidden:
  `trimDecodedChunks`), and a step, or a column read again, copies its
  rows out of them instead of decoding them again (an
  HLCA-shaped `.h5ad`, 1 % blocks of 19 cells: the last step decoded 33
  embedding chunks and 2 of `cluster`, not 84 and 46). A gene index's pass
  in store order over a CSR matrix whose decoded chunks pass that budget
  reads through a window of its own instead (`AnnDataReader._scanChunks`:
  the two chunks of `indices` and of `data` a slab shares with the next):
  each pass reads each chunk once, in order, and kept by recency the chunks
  are gone before a later pass comes back to them, so keeping them only
  held memory (a 50 % uniform share of 5M cells: 132 MB above every
  cell's). Describing the open cells' CSR matrix makes none of its arrays
  (`_requireMatrixElement` reads its shape; the stored matrix was checked
  at opening): its arrays in open order are its kept rows' entry runs in
  open order, 2.5M runs and 140 MB for that share, which a pass in store
  order never uses. Row runs hold 32-bit positions wherever every one fits
  (`runsOfPositions`). Measured (Node, the reader worker's handler, 5M-cell
  CSR `.h5ad` over an emulated 20 ms link, first gene with `genes: 'all'`):
  50 % uniform 3.98 s, peak ArrayBuffers 922 MB, RSS 1,582 MB → 3.86 s,
  711 MB, 1,383 MB; every cell 3.49 s, 599 MB, 1,119 MB. Each read of an
  `obsm` array names the open cells it read and what for in the read
  activity (`open cells 7,435–20,000 (embedding, cells added)`). Only a
  `uniform` or `blocks` subsample or a priority-ordered file opened whole
  is priced at all: the chunks of any other are never looked at. The cost
  of showing a cell is what the prediction needs: each opening's time from
  its cells in the page (their GPU stores made next) to the first frame
  that draws them (`recordFirstFrame`, from `reloadActiveDatasetInPlace`;
  a frame that drew none of them, or cells replaced first, records
  nothing). An opening also costs what showing any dataset costs (about
  25 ms in Chromium, against 9 ns a cell: 10K cells took 27 ms, 5M 70 ms),
  so the measure is the time a cell adds, the least-squares slope of the
  records' time against their cells (through the origin while they hold
  one size, or when the line would give an opening a negative cost), each
  earlier record counting half as much as the next, so no single slow
  moment keeps sizing later openings; kept in the browser under
  `cellucid:opening-cost:v2` as the records' decayed sums. An in-place resize is not recorded: the GPU
  keeps the cells kept and takes only those added, while every open
  cell's colours and filters are made again, so its time is neither its
  added cells' cost on an opening nor its open cells'. A subsample of a file
  in its own order whose first-view embedding is chunked much coarser than
  its blocks grows from few chunks instead (`progressive-open.js`,
  `chunkMajorOrder`, `chunkedPlan`): the same kept cells, taken chunk by
  chunk in a seeded random chunk order (each chunk's cells in priority
  order), each step a `chunks` request (`subsample-plan.js`, `chunkStep`)
  ending at a chunk boundary so each chunk is read once; then one
  `reordered` reading puts the same cells in priority order, every per-cell
  array moved with its cell (`DimensionManager.reorderCells`,
  `DataState.reorderOpenCells`), the frame kept. It is planned only over a
  measured link (any source's) and when its first view is predicted at least twice as soon
  as the priority order's (or as every cell's, when that is not
  progressive), each ordering priced by the embedding chunks its first view
  touches -- and, as every growth of a share, only when its cells are all
  on the map no later than one reading of them would put them
  (`grownLater`: the same bytes in more readings, a loaded round trip each,
  every cell shown again at the end), so over a priced link a share opens
  whole (50 % uniform of 10M cells on a local web host: grown from chunks,
  all 5M at 4.5 s; whole, 2.8 s). Every cell of a file in its own order
  never grows, and opens in its default subsample's cell order only when
  its gene sample brings the first gene sooner by more than the order
  costs the opening end to end -- made and timed at opening, its runs
  found a block at a time (`blockOrderRuns`), and each of the
  `ORDER_OPENING_PASSES` passes over the open cells it adds (its copy in
  the reply, the page's cell identity, the first view's embedding and
  codes gathered and scattered) priced at that time (`orderOpeningMs`;
  checked on a browser trace: made in 11 ms, 60-70 ms more) --
  (`subsample-plan.js` `orderPays`, the
  sample's priced first gene `geneTimeFit` `sampleMs` against every cell's
  `neededMs`): a sample of 9.84M of 10M cells over a local web host saved
  tens of milliseconds while the order cost the opening hundreds and made
  it grow in two steps (all 10M at 3.4 s against 1.5 s in file order). The plan is made after the store's structure is read, and
  opening reads no gene name (the var index is read when the names are first
  asked for, after the first view), so nothing else is priced as the first
  view's: a
  blocks share, whose first cells in priority order lie in few chunks
  already, keeps priority order; a uniform share of a coarsely chunked file
  grows from few chunks. A share opened in an order its smaller sizes are
  not prefixes of -- one spread through a file in its own order, opened in
  the file's cell order (`cellOrder`), or a share of a priority-ordered file
  opened in file order -- grows instead through the first cells of that
  order (`prefix` requests, `subsample-plan.js` `prefixStep`, carrying the
  `order` and the fitted sample): each step's first cells and gene sample
  are the share's, so each grows the one before in place (`sharedOpenCells`)
  and the last is the share itself, also in place; the chunk-major plan is
  weighed against it as against a priority order. A step's cells are
  recorded with the subsample's specification of their count, which does
  not give them again (`reproduced: false`; a prefix step's `part`), so a
  session saved during the growth opens exactly them (in the same order:
  the request names the cell order). A later visit's first view is priced
  by each chunk's own bytes to fetch (`firstViewCost` per chunk once
  anything is kept, `storedFor` for decoding a kept chunk), and it never
  opens fewer cells than those whose chunks are all kept
  (`progressivePlan`); a folder on a web host read at a content version
  (a GitHub folder at a commit) is taken over by every reading again
  without asking its host (`sourceKey`), so its growth steps stay in place:
  GitHub's raw host shows a page no validator, and every step read other
  cells than those shown (Suo: growth stopped at 32,768 of 561,947 cells,
  and its revisit at 750).
  An analysis's first pass on the first cells of priority order reads those
  cells' part of each gene the same way (`getMatrixCellSubset`).
- The page keeps the bytes of remote files its reader workers fetched
  (`data/byte-range-cache.js`, bounded by its `fetchedRanges` memory share,
  least recently used ranges dropped first), and every reader worker reaches
  it through a `MessagePort` posted when the worker starts
  (`formats/io/range-cache.js`), closed when the reader closes or the page is
  left (`AnnDataReaderClient`); the cache itself lives as long as the page. A remote file's reads ask the page for the
  bytes it holds first and fetch only the rest, which the page then keeps
  (except bytes only a cell export reads, `ReadActivity.readOnce`: read once
  for a file written out, they would only push out what other reads use, at
  a copy and a message each -- 812 MB for 21,285 cells of HLCA), and finds
  by binary search, under a key naming the file's bytes exactly (its identity's key and
  version, or its URL, size and `ETag`; a file known by no validator is not
  kept). A reader asks the page once a file which ranges it holds of it,
  adds the ranges it hands the page itself, and asks for a read's bytes only
  when the read touches them: a read of bytes the page does not hold goes to
  the host at once (on HLCA, 99 of the first 245 reads had waited a round
  trip each for a page busy opening the file). A read that touches them
  waits for the page's answer no longer than one round trip
  of its link (before the link is measured, no longer than the host is
  given to answer), since the page answers from its main thread, which a
  long task may hold; then it fetches, and the late answer is dropped. A
  cancelled read stops waiting at once, and the request's round trip is
  timed from when it went to the host. So a dataset opened again in a new worker (a worker the page ended
  and started again, the same file in another reader) fetches nothing the
  old one fetched: on an atlas-shaped CSR `.h5ad` over 150 ms / 50 MB/s,
  re-opening at 10 % after 5 % blocks took 1.1 s and 12 requests before and
  none now, and its first gene 6.0 s / 81 MB before and 2.7 s / 34 MB now.
  A read answered in part from the page measures no link. Pieces the page
  holds split a read into the parts it does not, each a request of its own;
  a held piece smaller than the gap the link makes worth fetching to save a
  request (`planMergeGap`, the byte reader's own pricing, once the link is
  measured) is fetched again with its neighbours instead: 100 genes of a
  dense `.h5ad` read after a few others, 150 ms / 5 MB/s over HTTP/1.1,
  took 208 requests and 6.1 s while the earlier genes' chunks split the
  read, and take 22 requests and 2.9 s (its floor 3.0 s), fetching 4 % of
  the bytes again.
- A Zarr folder on a web host (every object a file of its own: the sample
  datasets, GitHub, Hugging Face, S3 and GCS folders, a Python or R
  server's store) keeps its objects in the same page memory, for the tab,
  under a prefix naming the store's version: its address (the link's
  identity, never its credentials) and the version its provider names its
  content by (a GitHub commit, a Hugging Face git tree: shared even where
  the host shows no validator, as GitHub's raw host does not), else the
  validator its host gave its consolidated metadata (`range-cache.js`
  `storeObjectsKey`); a folder named by neither keeps nothing. Its metadata is read from the host at each
  opening (from this browser's kept copy only under a version that names
  its content, below), so a rewritten store is read under a new prefix,
  and its old objects only age out. A reader asks the page once a store which objects it holds
  (not once an object), then asks for an object only when the page holds
  it, within the same one round trip; whole objects, ranges of shards and
  shard indexes are kept apart (`storeObjects`). So an embedding or field
  read again, the first chunk the counts probe read, and every opening of
  the store in another reader (a share raised or lowered) fetch nothing the
  tab fetched before.
- Beneath the page's memory, what a reader fetches of a remote file is kept
  in this browser across visits: its **downloaded parts**
  (`formats/io/kept-parts.js`, `withKeptParts`, in the reader worker). Every
  range fetched is kept as it arrives, whatever the visit's length: each
  reader keeps a **journal** (`transcode/space.js` `openKeptJournal`), one
  directory whose manifest is written first, so it is kept from its first
  byte; each piece is appended to its `d.bin` through a synchronous access
  handle, then a record of where it is of the file (`j.bin`: framed by its
  length and CRC-32, with the piece's page checksums). A write is the
  operating system's once it returns, so a tab closed the moment its reads
  ended, a worker ended from outside or a browser that crashes loses at most
  the piece being written; a machine that loses power may lose what the disk
  had not yet stored, whose records then do not read whole or do not match
  their checksums and are not used (`parseKeptJournal` stops at the first).
  The page asks its reader workers to have the disk hold what they wrote
  when it is hidden or left (`flushKeptParts` on `visibilitychange` and
  `pagehide`; a reader closed, or a page left, asks its worker to leave,
  `{leave: true}`, which stops its reads, flushes and ends itself). A
  journal is closed, never appended to again, once its reader's fetches go
  quiet for the link's retransmission timeout, when the page is hidden or
  left and when the reader closes (the next piece begins a new one): a
  journal still being written (its writer holds its writing lock: another
  tab, or the previous generation of the dataset in another worker during
  a read) is left to the next opening, because a file another worker still
  appends to cannot be read in place. A file whose reads name no version (Google
  Drive's public reads) records a confirmation once its provider confirms
  the version it opened as still holds (asked once its reads go quiet for
  the link's retransmission timeout and when the page is hidden,
  `sources/access.js` `version`), and a later opening uses only the pieces
  a confirmation follows. Every piece has a CRC-32 of every 16 KiB page
  (`page-checks.js`); every journal is filed under `{parts: identity, size}`: the file's
  identity as this opening established it (layer 6, *When a dataset
  changes*) and its exact size. A read is answered from the reader's own
  memory (cached metadata blocks, the surplus of merged requests: within
  its budget, plus the requests a read still running made, held until that
  read ends, `ReadActivity.finished`, and counted under `fetchedRanges`, so
  a heap collection fetched inside a reference chunk's request is never
  fetched again before its references are decoded; what all of them hold is
  bounded by the memory model's `fetchedRanges` budget, the least recently
  used going first, `configureByteReaders`, since a long read such as a cell
  export asks for each region once and held every request it made until it
  ended: 210 of 216 fetched MB for a 21,285-cell export of a 10M-cell file)
  first,
  then from the kept parts (a read of this machine's disk, waited for only
  once, before the reader's first plan: `kept.ready`; the kept ranges are
  held in order, `SortedRanges`, so finding them and every lookup cost a
  sort and a search whatever their number: a first visit whose analyses
  read 200,000 ranges opens again in about a second), then from the page's
  memory, and only what none holds is fetched. The byte reader's plans treat
  kept bytes as present (`byte-reader.js`): no request asks for them, two
  ranges are never merged across kept bytes, a partly kept read fetches only
  its missing pieces, and a metadata read whose block holds kept bytes is
  read exactly (a whole block would fetch them again).
  `ByteReader.bytesToFetch(start, end)` says what a range still costs
  (neither kept nor cached; `keptReady` resolves once the kept parts are
  known), and `storedBytes(selection, {toFetch: true})` (`formats/store.js`)
  what an array read would fetch: the cost models price with them (the
  first gene's notice, auto-fit by gene time, a progressive opening's first
  view). A file a transcoder reads whole (an R object, a MatrixMarket file)
  keeps no parts: its converted store is what is kept of it. Nothing is kept
  as parts of a source on this machine or a notebook kernel
  (`kept-parts.js` `sameMachine`: a loopback address — `localhost`,
  `127.x`, `::1`, `cellucid serve`, a tunnel — or the widget's host,
  `EMBEDDED_HOST_ORIGIN`), nor of a store a Python or R server serves
  (`served-anndata.js` `keptParts: false`; a store served from an object in
  memory is named anew every run): a copy spares no download. A Zarr store
  folder on a web host keeps its objects (`formats/zarr/stores/http.js`
  `keeping`: whole objects and a shard's ranges, each with the validator the
  host answered it with; under such a version nothing of it changes, so its
  consolidated metadata is kept too, and the objects its host answered as
  absent, which a later visit does not ask for again: its only request is
  the provider's version) only under a version that names its content: a
  Hugging Face folder, named by its git tree (`hf-tree:<tree id>`, the
  folder's `oid` the Hub's `paths-info` answers: a hash of every file below
  it by its own content name, so a commit elsewhere in the repository leaves
  it), threaded from the provider through `link-source.js` and
  `loadRemoteFolder`; a GitHub folder, read at the commit that last wrote
  anything below it (`commit:<sha>`, one request to the commits API) and
  filed under its link's address (`key`: the ref the link names, so another
  commit's kept objects are stale). While a folder is open, an object
  answered with another validator than before (in this opening, or kept) is
  refused as changed. A folder on a
  plain host (GitHub Pages, S3, GitHub by branch when its API refuses)
  keeps nothing across openings (identity `{kind: 'url'}`, its opening's
  only): its consolidated
  metadata's validators say when the metadata was written, not that every
  object is as it was (`rclone`/`gsutil rsync` and in-place writes rewrite
  objects under the same metadata). When anything
  is kept of a file at its address (`space.js` `keptAtAddress`, from the
  directory names alone), its probe asks only for the version (a `HEAD`),
  not its first 512 KiB, which then come from what is kept when the version
  is the same (`data/sources/remote-file.js`). Google Drive and OneDrive
  read their file's head beside their own API request (it proves the file
  can be read without a sign-in, and resolves the download link): one byte
  of it when something is kept (`headBytesFor`). Figshare reads its head
  after its item names the file, one byte of it when something is kept; a
  link naming the file by its id reads the head beside the item, whole (512
  KiB transfer in less than the round trip a read after the item would
  cost).

### 3) `DataState` field residency (bounded in bytes)

`DataState` keeps loaded field arrays — obs fields and gene vectors — and
releases inactive ones least recently used once they exceed their budget, a
share of the memory allowance (`state/managers/field/residency.js`,
`data/memory-budget.js`), and all of them when the page is hidden. The active field, fields with filter
state, and fields still loading are always kept. This is what stops “scroll
through genes for 10 minutes → tab uses 10 GB”; the earlier entry-count bound
did not, because each field object kept its arrays alive. Categorical codes
stay in the type the file stored (uint8 stays uint8). An inactive gene whose
column has a sparse form (nearly every gene: `data/column-forms.js`) keeps
that form alone, whatever the budget — its entries instead of four bytes a
cell (40 MB at 10M cells) — and is made again from it, bit for bit and
without a read, when it is shown (`compactedColumn`); so browsing genes holds
one gene's column at a time.

### 4) Community annotation caches (IndexedDB + localStorage index)

Community annotation stores:
- raw downloaded JSON file content (IndexedDB)
- a small path→sha index (localStorage) for “which files changed?”
- the working annotation session itself (localStorage), which holds votes,
  comments, suggestions and profile — dataset-derived content, not preferences

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
downloaded (layer 2, `formats/io/kept-parts.js`), and the gene index of a
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
converted stays in memory for the tab. **Session › Open data › Stored**
(`app/ui/modules/kept-stores-line.js`) sums up what is kept; **Manage…**
opens the Stored dialog (`app/ui/modules/kept-stores-dialog.js`), which lists
every item (`space.js` `listKeptStores`: a converted store, a whole download,
a file's downloaded parts, with the bytes they spare the next opening
(`saves`), or a gene index per file, matrix and cell order, named by the manifest's
`about`, its file's name and host from `keptAbout`, never a link) grouped by
dataset (`file`), and removes the rows chosen, a dataset's, or all, that no
tab holds (`removeKeptStores`). Under the datasets, *Unsaved highlights*
lists the highlight drafts (`app/highlight-draft.js` `listHighlightDrafts`,
IndexedDB `cellucid-highlight-drafts`, a row per dataset by the name the
draft records, its size and when it was saved; ids `draft:<dataset id>`,
removed with `removeHighlightDrafts`, heard through
`onHighlightDraftsChange`), then *Files to reopen* the kept handles
(`app/local-file-handles.js` `listLocalFiles`, IndexedDB
`cellucid-local-files`; ids `local:<dataset id>`, forgotten with
`forgetLocalFiles`, heard through `onLocalFilesChange`). Both stores are
`utils/indexed-db-store.js` `indexedDbKeyValueStore`. Then one row lists the hosts
whose connection speed is remembered (`data/host-link-memory.js`
`rememberedHosts`, in `localStorage`, layer 7) with its own **Remove**
(`forgetHostLinks`); **Remove all** forgets them too, and
the Stored row counts them (*2 datasets · 540 MB · 1 highlights draft ·
1 file to reopen · 2 connection speeds*; with no dataset kept, *Nothing
kept · 2 connection speeds*). *Where it's kept* names what the site's
storage holds (*Holds*: files, highlight drafts, files to reopen,
connection speeds). Listing and removal run one at a time in a
page: each asks whether a tab holds a directory by requesting its lock if
free, which, while it lasts, is a holder to any other such request. The
dialog's foot (*Where it's kept*) comes from `keptStoragePlace`: whether the
window gives the site an origin private file system at all (asked for, never
inferred from the engine: `navigator.storage.getDirectory()` answers or
refuses), `navigator.storage.estimate()` (usage and what the quota leaves,
labelled as the browser's figures) and `navigator.storage.persisted()`
(*Kept when space is low*); `keepKeptStoresWhenSpaceIsLow` calls
`navigator.storage.persist()` only on the person's click (**Keep it when
space is low**). Without storage the Stored row says *Kept for this tab only* and the dialog lists
no dataset.

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

Every opening asks the source for its version before anything kept of it
is used. An opening from the source's own description makes that check
itself, in the same opening: the probe (a `HEAD` and the head read), or the
provider's answer that names the version (a record, a listing, an item). A
reader opened again from an earlier opening's description (the dataset no
longer shown and then read again, a page back from the back/forward cache:
`resumed`, `anndata-reader-client.js`) asks it explicitly first, whatever is
kept (`kept-parts.js` `confirmFileVersion`, one small request a file): a
provider that alone can tell (Google Drive's public reads: the item's
checksum; OneDrive: the item's cTag) is asked; a file served only whole has
its validators asked with a `HEAD`; any other file has its first byte read
from the host, past every kept copy and the page's memory, under the
validators it was opened with. A version that no longer holds forgets every
kept copy of the file and fails the opening as changed ("open it again"); a
reading again with other cells within one opening (Apply, a growth step)
keeps the version its opening checked.

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

The raw bytes kept (downloads and downloaded parts) have their own revision,
`KEPT_BYTES_REVISION`, apart from converted stores and gene indexes: they
are the file's own bytes, so it moves only when how they are laid out and
checked changes, held by golden bytes of fixed inputs
(`tests/kept-parts-contract.test.mjs`); a release that changes readers or
decoders keeps every download.

Every opening asks the source for the version before anything kept or
cached is used: the probe's `HEAD` past the HTTP cache (`no-store`), the
provider's API, or a server's dataset list, read at that opening, never a
listing remembered from an earlier one. Before reading, `open.js`
(`forgetOtherVersionsOf`) forgets what is kept of the file's other
versions, every kind (`space.js` `forgetOtherVersions`: every directory of
the same address, `addressOf`, whose `fileKey` is not this version's loses
its manifest and is removed once no tab holds it), and the page's memory of
their bytes (`range-cache.js` `forgetOtherVersionsInMemory`). A file that
changes on its host while it is open is refused by its byte source
(`changed`: "open it again to read the new version"); the kept parts then
forget every version of the file, every kind (`forgetKeptFile`), and the
page's memory drops it (`forgetFileInMemory`).

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

An opening sizes its gene sample (by gene time,
`subsample-plan.js` `geneTimeFit`) and a progressive first view at the
speed of its link. Its first requests are small, a round trip each, so what
it measures is a lower bound (HLCA's opening measured 244 ms and about
640 KB/s; its gene reads then ran at 1.5–3.4 MB/s). So what each host
delivered is kept across visits (`data/host-link-memory.js`, `localStorage`
`cellucid:host-link:v1`): by origin (`scheme://host[:port]`, never a path,
a dataset, a query or a token), the round trip, the round trip under load
(`LinkMeter.loadedLatencyMs`: a revisit's first reads are priced at it,
scaled to the opening's round trip and never below its own, until the
opening rests on as many bytes), the rate, the bytes they rest on
(`LinkMeter.deliveredBytes`) and when they were measured. The rate
is what the reads sustained: bytes delivered over the time any request was
in flight (`LinkMeter.busyMs`), round trips and all, never the link's fastest
moment (`bytesPerMs`), which overstates a long read: HLCA's fastest delivery
was 6.2 MB/s while its gene reads sustained 1.6 MB/s, and sizing from the
fastest gave a later opening's first gene three times its target. A visit's
measurement is blended into the record, each weighed by its bytes, the
record at half its weight (`rememberHostLink`); at most
`HOST_LINK_RECORDS` (32) hosts are kept, the least recently measured
dropped. The page records it from the reader's link measurements and hands
the host's record to each opening (`hostLink`); the reader prices with
`pricingLink`: the opening's own measurement once it has measured its
bandwidth with a delivery bounded by transfer (`bandwidthMeasured`: a slow
day is priced as slow within a few requests, never at a remembered fast
rate), or once it rests on as many bytes
as the record, before that the record's rate scaled by the ratio of the two
round trips (a longer path delivers less) and never below the opening's
own rate. The opening's round trip is one measured on a connection already
made (`LinkMeter.warmLatencyMs`: a request begun after another completed),
never its first requests', whose time holds the connection's set-up; until
it has one, the record's rate stands as it was measured. `rememberedHosts` lists the hosts, and `forgetHostLinks` removes
them (the Stored dialog's *Connection speeds remembered* row, and its
**Remove all**).

A host on this machine (`onThisMachine`: a loopback origin -- `localhost`,
`127.0.0.0/8`, `[::1]`: `cellucid serve`, Python's and R's served stores, a
web server run here) is marked in the seed (`hostLinkSeed` `thisMachine`),
and `pricingLink` prices it as this machine at rest, whatever its link
measured (`thisMachinePricing`): at no less than a local file's prior (a
disk, `LOCAL_POLICY.bytesPerMs`, 1 GB/s), at its shortest round trip
(`latencyMs`, never its round trip under load) and at no time a request
beyond that (`requestCostMs` 0). What its reads measure past that is the
machine's own load, never a network: its opening's metadata reads are each
mostly their fixed cost (a 10M-cell `.h5ad` on a loopback host in Chromium
measured 45 MB/s from 0.9 MB of them, a Zarr v3 store 9 MB/s from 51 KB),
and priced at that its first gene's read of every cell (5-7 s in fact) was
priced at 15-37 s, so the first gene was read for a sample of 22-68 % of
the cells; with the machine busy (another browser starting beside it), one
of those reads was slowed past three round trips and counted as a delivery
bounded by transfer (`bandwidthMeasured`) at 0.1-13 MB/s, or the requests
waited and measured a round trip under load of 863 ms (1.5-7 ms at rest),
and the first gene was read for 19,779 to 1,947,096 of 10M cells (5 of 6
such openings). A local file's own link is priced the same way: its rate
by the meter (`LinkMeter` `_pricedBytesPerMs`, carried by
`snapshot().bytesPerMs`), the rest by `pricingLink`, which its snapshot
asks for (`snapshot().onThisMachine`: its prior is a disk). A tunnel to a
remote server (`ssh -L` to a loopback port) is priced so too. A notebook kernel's host
(`EMBEDDED_HOST_ORIGIN`) is not marked: a kernel may run across the internet
(Colab, a JupyterHub).

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
     - The dataset base URL is logged by `main.js` when debug is enabled

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

### Symptom: “Unexpected token `<`” during loading

Likely cause:
- The server returned an HTML document for a JSON request (`datasets.json`,
  or a Zarr metadata key).

Fix:
- Fetch the URL in a new tab; confirm it’s real JSON and returns `404` when missing.

### Symptom: a load reported “invalid format” but the file is fine

Read the message again — this class of failure is now split four ways, and
“invalid format” means only the last one:

| Message contains | Cause |
|---|---|
| “was cancelled” | An abort — a newer selection superseded this load. Re-thrown untouched so it is never relabelled. |
| “is longer than the longest text this browser can hold” | A JSON document (a catalog, Zarr metadata) exceeds the engine's longest string |
| “response body transfer failed” | The connection dropped mid-body |
| “must contain valid JSON” | Genuinely malformed content |

If you add a new bounded or cancellable read, preserve this: check for an abort
first and re-throw it unchanged, distinguish the range failure from the parse
failure, and give each its own error code. Collapsing them back into one message
is what made these failures undiagnosable before.

### Symptom: “Gzip decompression not supported”

Likely cause:
- The browser does not provide native `DecompressionStream('gzip')`.

Fix:
- Update to the current stable release of Chrome, Edge, Firefox, or Safari.

---

Next: {doc}`10_sessions_persistence_and_serialization` (how loaded state is saved/restored).
