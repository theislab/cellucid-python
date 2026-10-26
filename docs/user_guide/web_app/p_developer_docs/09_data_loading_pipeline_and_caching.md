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
| `var`, `var/{g}` | the gene list the file opens with, and gene `g` (`loadGeneList`) |
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
- open the AnnData file in the reader worker (its node tree and
  `uns['cellucid']`), and read the dataset's identity
- load the obs inventory → attach obs field loader
- optionally load the var inventory → attach var field loader
- load the positions of the embedding the dataset opens on (the default of
  its opening dimension, `dimensionManager.getDefaultLayoutKey()`), read
  beside the obs field lists rather than after them; every other enabled
  embedding, including the other dimensions `uns['cellucid']['defaults']`
  names, is read (and its values checked) when a view first shows it. Before
  its first render, Suo now reads 5.5 MB in 27 requests (its metadata and its
  3D embedding). It used to read all three embeddings, 10.9 MB in 139
  requests, one after another. Over HTTP/2 at 80 ms and 2 MB/s its first
  render went from 6.9 to 3.8 s in WebKit; the 0.9.4 app took 4.3 s with the
  prepared export
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
  and handed to the same local reader and intent as the file inputs. The
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
- Each is opened as a candidate by the AnnData reader worker and adopted only
  when it opened completely; it is served through `anndata://` URLs. What it
  is is decided by the reader from its content (a file's first bytes), so a
  misnamed file opens as what it is. Its dataset id is
  `uns['cellucid']['id']`, else `<format>_<name>` (`h5ad_<file name>`,
  `zarr_<store name>`: `defaultDatasetIdOf`); its name is
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
  modification time, inode): the browser keeps what it read and revalidates
  it on a reload (`304 Not Modified`, nothing sent again), and a range asked
  with `If-Range` naming an older version is answered with the whole current
  file, so kept bytes of a changed file are never mixed with new ones.
- The browser measures the link to the server like any other (see
  [Reading AnnData directly](#reading-anndata-directly)): `cellucid serve` or a
  notebook on the same machine reads like a local file, the same server
  behind a remote JupyterHub proxy like a file on the internet.

Common failure modes:
- mixed-content blocking (https app trying to fetch non-local http server)
- missing CORS headers
- server schema/version mismatch

### `jupyter` (a viewer a host program opened)

Code:
- `cellucid/assets/js/data/viewer-host-source.js` (the source)
- `cellucid/assets/js/data/viewer-channel.js` and `data/viewer-protocol.js`
  (the channel and its rules)
- Protocol architecture: {doc}`../../python_package/e_jupyter_hooks/05_architecture_viewer_channel`

Behavior:
- A host program (Python `show_anndata()` in a notebook or a script, R
  `cellucid_view()`) serves one dataset and opens the page it serves with
  `?viewer=true&viewerId=…&viewerToken=…`; the page is the same whether it
  sits in a notebook output frame or a browser tab. Its `/_cellucid/info`
  names the program (`runtime`), which the Dataset box shows.
- The dataset is read like a server's (`served-anndata.js`). Commands from
  the host and events from the viewer travel on one WebSocket to the same
  server (`/_cellucid/ws`), authenticated with the viewer id and token: no
  `postMessage` to a parent window and no polling. The host's pings keep an
  idle socket open; a dropped socket is reopened once, and one health
  request then tells a stopped host from a proxy that does not pass
  WebSockets.

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
split layers are joined in one pass, and a matrix that arrives by cell or by
feature is kept as it arrived. A large gzip stream (1 MiB or more) is
decompressed on a worker of its own (`formats/rds/gunzip-worker.js`), since
`DecompressionStream` inflates on the thread that calls it; a MatrixMarket
file's entry lines go to the decode workers in blocks of about 1 MiB and come
back as typed arrays (`formats/text/matrix-market-entries.js`), taken in file
order and collected a block at a time.

| Layer | Code | What it does |
|---|---|---|
| Bytes | `formats/io/` | Byte sources (a local `File`, HTTP range requests, memory) and **one byte reader per file**: reads asked for in the same task are planned together, neighbours closer than the merge gap become one request, requests are split at the request size so large reads run in parallel (`concurrency` at once: six per HTTP/1.1 host, the browser's limit), requests waiting for a slot take their neighbours with them when round trips cost more than bytes, metadata is cached in aligned blocks (16 KiB from a web host: an anndata-written `.h5ad` keeps each element's headers within a few KiB beside its data, so larger blocks fetched mostly data; a 1M-cell `.h5ad` opened at 80 ms/2 MB/s in 1.4 s instead of 4.2 s), and a read's `AbortSignal` cancels it. A host that content-codes a file even in a range (`Content-Encoding`, hidden from the page; NCBI does so for GEO's `.rds`) is told by the probe (`probeHttpResource`: the ranged read fails while the whole file answers), and the file is downloaded whole. The merge gap and request size follow **the link as measured** (`link.js`, see below). `inflate.js` decompresses with the browser's native `DecompressionStream`, to an exact or bounded size, feeding it the input in 32 KiB slices so the size bound holds while it decodes (WebKit's decompressor inflates whatever is written to it at once: a chunk written whole expanded a 60 MB bomb past its 1 KiB declared size). `decode-pool.js` decodes chunks on decode workers (`decode-worker.js`, nested in the reader worker, started when a file opens: as many as the cores allow, less one for the reader thread and one for the page, and as many as the `decodeWorkers` memory share holds), small chunks that wait together going to a worker in one message, a background pass's chunks behind every other read's, so a background pass over a matrix leaves the reader thread free; where no worker can be made, chunks are decoded in the reader thread. The same workers take the reader thread's work on decoded values (`runOnDecodeWorker`: a function of a module that depends on nothing else, its arguments transferred, with a context each worker is sent once), as a cell subset's runs of genes. `activity.js` counts what each operation fetched and decoded, and which shared loads it reads for (`keys`: a job joining one adopts the read, `read-activity-log.js`). |
| Formats | `formats/hdf5/`, `formats/zarr/`, `formats/zip/` (and, for the transcoders, `formats/rds/`: R's serialization, gzip native, bzip2 and xz in JavaScript; `formats/text/`: MatrixMarket) | Cellucid's own HDF5 reader (every superblock, object header, group, attribute and chunk-index version; deflate, shuffle, Fletcher-32 and LZF), Zarr reader (v2 and v3 with sharding; gzip, zlib, Blosc, LZ4 and Zstandard; a picked folder, a ZIP, or a folder on a web host via consolidated metadata) and ZIP reader. Each opens a file as the **array store** of `formats/store.js`: the whole node tree (names, shapes, element types, attributes) is read at open, with a tree level's reads in flight together, so navigation is synchronous; array values are read on demand, fetching only the chunks a selection touches, in parallel, each decoded once (a per-file decoded-chunk cache serves neighbouring reads); of an unfiltered HDF5 chunk (anndata's default sparse `.h5ad`) only the run of bytes a selection needs is fetched (`partialChunks`), so 1,000 cells of a 1M-cell CSR `.h5ad` over 30 ms/5 MB/s read 4.5 MB in 12 s, not 121 MB in 25 s. |
| AnnData | `data/anndata/` | **One** AnnData reader for both formats: the element rules of the current layout (`elements.js`, `dataframe.js`, `sparse.js`), files from anndata 0.6/0.7 and SeuratDisk presented in the current layout (`older-layouts.js`), slabbed reads with read-ahead and exact conversions (`values.js`), gene columns and cell rows per layout (`matrix-access.js`), and opening plus error explanations (`open.js`, `views/`, `transcode/`). |
| Worker | `data/anndata-reader-worker.js`, `data/anndata-reader-client.js` | The reader's methods over `postMessage`, results transferred; live `activity` and `transfer` events. |

**Genes and cells by layout.** Only what one gene (or cell) needs is read:

| Matrix | A gene (column) | A cell (row) |
|---|---|---|
| CSC | its own entries: indices and values fetched together | indexed by cell once (the index of the transposed matrix) |
| CSR, over a link as fast as a disk (a local file, the same machine or network) | indexed by gene once in the background, at open, indices and values into a column store (in memory, or on disk when it does not fit); a gene asked for meanwhile is answered from the same pass, and every later one is a copy out of the store | its own entries |
| CSR, over a slower link | indexed when the first gene is shown, from the **indices** alone; a gene's values are then fetched as the ranges its entries occupy; once the ranges already read and those asked for together would cost as much as every value (a request costs a round trip of the measured link), every value is read once, after which every gene is answered from memory (never more than twice the cheaper choice in hindsight) | its own entries |
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

**Two forms of a gene.** A view colours by a gene's dense column (float32,
one value per cell), made in the reader worker and transferred, so the page
never expands one. A computation reads the same gene as its non-zero entries,
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
worker checks the `obsp` matrix in storage order (bounds, order, zero
diagonal, exact symmetry) and builds its edges as 12-byte records (both
cells and the `float32` relative strength) already in the renderer's fixed
random draw order, every cell's neighbours as one ascending list (the k-NN
tool's adjacency and the renderer's visible-cell walk), the forward offsets,
and the graph's largest weight and degree: about 20 bytes an edge, nothing
ordered or copied on the page. Checking and extracting are two passes in
storage order: the first checks and counts every entry and notes where each
line's upper half starts; the second walks each half once, placing each lower
entry where its partner's edge is and matching each upper entry against it.
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
indexes, positions of layouts no view shows, inactive fields, the analysis,
the adapter, the latent space held for outlier quantiles, decoded chunk
blocks, decode workers, undo history) takes a disjoint share of one allowance
from `data/memory-budget.js` (see {ref}`memory-budgets`): half the estimated
device, the shares summing to exactly one. A refused allocation halves the
allowance of the thread that saw it.

**Requests from the page.** The reader client sends the gene reads asked for
in one task as one request per matrix and form (dense or sparse), each gene
read once however many callers asked; a second caller of the same gene gets
its own copy. Every read takes an `AbortSignal`: aborting rejects the caller
at once and sends the worker `{cancel: id}`, whose reads then stop (a request
the network has not sent is never sent; a gene index other reads still need
goes on). A batched gene request is cancelled once every caller in it has
cancelled. The worker's open reply and its transfer events say how many
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
(`formats/io/link.js`): the round trip is the shortest request, the
bandwidth the fastest delivery seen, over one request or over a busy period.
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
  621 ms: six connections, not bytes, set that floor). The request size carries a round trip's
  worth for every request in flight (bandwidth × latency × concurrency,
  never less than the source's prior);
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
  drops to six). On a host proven to multiplex the concurrency keeps doubling
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
- a chunk B-tree of an `.h5ad` is read leaf by leaf as reads need it, and the
  rest of its leaves follow in the background when they transfer in less
  than a round trip, so a later gene is one round trip, not two;
- one `abort` listener per `AbortSignal` (`formats/io/abort.js`), however
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

Reads of many genes or cells at once (`getMatrixColumns`, `getMatrixRows`, or
many single reads issued together, as analysis does) go out as one set of
coalesced, parallel requests. An index built ahead of need is background
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

`app/read-progress-notifications.js` shows any read longer than 0.7 s
(except opening, which the startup view shows) as a progress notification
with the same figures, closed with what it read and how long it took.

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
against an earlier revision). The R object corpus is written again with the
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

- **Transient failures are asked again**, with growing waits
  (`RETRY_DELAYS_MS`: 250 ms to 8 s, six repeats), or as long as a
  `Retry-After` asks (seconds or an HTTP date; a cross-origin page sees it
  only where the host exposes it): HTTP 408, 429 and 5xx (but 501 and 505),
  and a connection that fails or drops. GitHub Pages answers a burst of
  requests with an occasional 503 that carries no CORS headers, which every
  browser reports to the page as a failed connection, so a failed connection
  is transient too.
- **A stalled request is aborted and asked again at once.** The link's
  `ResponseClock` (`formats/io/link.js`) applies TCP's retransmission
  timeout (RFC 6298) to the time requests wait for their answers: the timeout
  is the smoothed wait plus four times its variation, never under one second
  and one second before the first answer, doubled by every stall and computed
  afresh by the next answer. A request is stalled when neither it nor any
  other request of the link has received anything (an answer, or bytes of a
  body) for a whole timeout: while the link delivers bytes, a request queued
  behind them (in a browser's HTTP/1.1 queue, or behind an HTTP/2 host's
  other streams) is waiting, not stalled. A range of known length is also
  given the time its bytes take at the measured bandwidth, since a host that
  buffers a range (a proxy, a CDN fetching the object from its origin)
  answers only once it holds all of it. Bodies are read as streams, so a
  body that stops arriving is caught as well as an answer that never comes.
- **A failure that will stay is reported at once** (`worthAskingAgain`).
  The signs that a host will answer later keep the whole schedule: a
  `Retry-After`, a connection dropped in the middle of a body, a stall. An
  error status without `Retry-After` that repeats exactly (the same status
  and `Content-Length`: Zenodo answers a missing record with the same 500
  every time) is reported after its second answer, and a request that got no
  answer at all (no connection, or a host that does not let this site read
  it) after it was asked once more; the momentary 503 GitHub Pages sends
  without CORS headers is answered on that second request.
- **A failure that outlasts its repeats names the request**:
  `https://…/suo.zarr/obsm/X_umap_3d/3.0: the host could not answer (HTTP
  503). Asked 7 times.` A failure that is not transient (404 on a file, a
  short or changed answer, a range ignored) is reported at once.
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
against a fake `fetch`.

## Loader layer (what is read)

Most of the “what do we read?” logic lives in:
- `cellucid/assets/js/data/data-loaders.js` (it accepts only `anndata://`
  URLs)
- `cellucid/assets/js/data/base-anndata-adapter.js` (what each path means in
  the file)

What a dataset generation reads (`loadDatasetGeneration`):

- the **identity**: dataset id, name, description, citation and source (from
  `uns['cellucid']`), counts, embeddings (`obsm`), vector fields, and whether
  there is a graph (`obsp`). Its id must be the one that was selected.
- the **obs inventory**: the `obs` columns, categorical or continuous, and the
  outlier quantile route of each categorical one.
- the **var inventory** (when the file has genes): the gene names of the chosen
  gene-name column (`var`).
- the **graph description** (when the file has a graph in `obsp`): edges are
  read when first shown.

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
  description; the reader keeps an open file's node tree and decoded chunks
  (`formats/store.js`).

### 3) `DataState` field residency (bounded in bytes)

`DataState` keeps loaded field arrays — obs fields and gene vectors — and
releases inactive ones least recently used once they exceed their budget, a
share of the device memory the browser reports (`state/managers/field/residency.js`,
`data/memory-budget.js`). The active field, fields with filter
state, and fields still loading are always kept. This is what stops “scroll
through genes for 10 minutes → tab uses 10 GB”; the earlier entry-count bound
did not, because each field object kept its arrays alive. Categorical codes
stay in the type the file stored (uint8 stays uint8).

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

:::{important}
Layer 4 persists data derived from the user’s dataset (layer 5 lives only in
memory), so it belongs in the privacy model, not just the performance story. Anything you add here that
outlives a reload must be documented in
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
