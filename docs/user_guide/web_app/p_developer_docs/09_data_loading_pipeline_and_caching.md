# Data loading pipeline and caching

This page documents the **developer-facing data loading pipeline** for the Cellucid web app:

- how datasets are discovered and selected
- how the app derives a dataset **base URL**
- which files are fetched (and in what order)
- where caching happens (HTTP cache, in-memory caches, LRU field caches, IndexedDB caches)
- how remote server / GitHub / file picker / Jupyter modes differ

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
   - Choose *where* the dataset comes from (demo exports, local folder, remote server, GitHub repo, Jupyter).
   - Output: a canonical `{ sourceType, datasetId, baseUrl, metadata }`.

2) **File loading**
   - Given `baseUrl`, load manifests + binary assets (obs, embeddings, connectivity, etc.).
   - Output: typed arrays + field loaders attached to `DataState`.

The base URL is the critical bridge:
- UI code should never hardcode demo-exports paths; it should always flow through an exports **base URL**.
- Loader code should never assume “GitHub raw URL…”.

Code pointers:
- Source selection: `cellucid/assets/js/data/data-source-manager.js`
- Loaders: `cellucid/assets/js/data/data-loaders.js`

---

## Startup path (what happens at boot)

At startup, `cellucid/assets/js/app/main.js`:

1) Registers core sources:
   - `local-user` (file picker)
   - `remote` (cellucid-python server)
   - `jupyter` (only if in iframe context)
2) Initializes `DataSourceManager`, which also registers:
   - `local-demo` (static exports)
   - `github-repo` (GitHub-hosted exports)
3) Applies URL overrides (if present):
   - `?remote=...` connects remote server
   - `?github=...` connects GitHub exports
   - exactly one canonical `?exportsBaseUrl=...` may replace the configured
     demo exports base URL for that launch
   - `?dataset=...&source=...` selects a specific dataset from a registered source
4) Chooses the active dataset and sets `EXPORT_BASE_URL` to the active `baseUrl`.

Then the loader phase begins:
- load obs manifest → attach obs field loader
- optionally load var manifest → attach var field loader
- load dataset identity + embeddings metadata
- load points positions (for the active dimension)
- optionally load connectivity manifest + edges

---

## Sources (where datasets can come from)

All sources implement a similar interface:
- `isAvailable()`
- `listDatasets()`
- `getMetadata(datasetId)`
- `getBaseUrl(datasetId)`

### `local-demo` (exports base URL)

Code:
- `cellucid/assets/js/data/local-demo-source.js`

Behavior:
- Reads `datasets.json` under the exports base URL (`DATA_CONFIG.EXPORTS_BASE_URL`).
- `DATA_CONFIG.EXPORTS_BASE_URL` resolves in this order:
  1) exactly one `?exportsBaseUrl=...` query field
  2) `<meta name="cellucid-exports-base-url" content="...">` in `index.html`
- Uses `dataset_identity.json` inside each dataset folder as the required metadata anchor.

Common failure mode:
- `datasets.json` missing or invalid → local demo source is “not available”.

Cross-origin transport contract:
- Every sample request is one canonical absolute URL below the configured
  exports root.
- The exports host must permit the viewer origin with CORS.
- Credentials are omitted and redirects are rejected; an HTTP, CORS, URL, or
  payload failure terminates that dataset load.

### `github-repo` (public GitHub-hosted exports)

Code:
- `cellucid/assets/js/data/github-data-source.js`

Behavior:
- User enters a path like `owner/repo/path`.
- The app resolves it to a raw URL like:
  - `https://raw.githubusercontent.com/owner/repo/<branch>/path/`
- It then behaves like `local-demo` (loads `datasets.json` + `dataset_identity.json`).

Common failure modes:
- `datasets.json` not found (wrong path)
- GitHub raw blocked by corporate environments
- unexpected redirects / rate limits (rare for public raw URLs)

### `local-user` (browser file picker)

Code:
- `cellucid/assets/js/data/local-user-source.js`

Behavior:
- Lets the user pick:
  - one prepared export directory with **Prepared**
  - one `.h5ad` file with **H5AD**
  - one `.zarr.zip` or `.zip` archive with **Zarr ZIP**
- Internally, the source resolves reads to “local-user URLs” (not real HTTP URLs) and the loaders know how to fetch them.

Common failure modes:
- denied or policy-blocked local-file permissions
- H5AD larger than the 512 MiB browser limit
- malformed current AnnData identity or non-finite scientific payloads
- invalid/multi-root Zarr ZIPs or archive/chunk memory limits

### `remote` (cellucid-python server)

Code:
- `cellucid/assets/js/data/remote-source.js`

High-level behavior:
- Connects to a server URL (e.g. `http://127.0.0.1:8765`).
- Validates health and server info:
  - `GET /_cellucid/health`
  - `GET /_cellucid/info`
- Enumerates datasets (endpoint is server-defined; the client caches dataset paths).
- Loads dataset files via HTTP from the server base URL.
- Opens a WebSocket for live updates when the server advertises a WebSocket URL.

Why remote is important:
- H5AD remains read-only-backed on the Python side while the browser requests
  fields and genes on demand.
- Zarr is materialized eagerly by `anndata.read_zarr`; browser requests remain
  on demand, but the Python process must have memory for the store.

Common failure modes:
- mixed-content blocking (https app trying to fetch non-local http server)
- missing CORS headers
- server schema/version mismatch

### `jupyter` (embedded viewer)

Code:
- `cellucid/assets/js/data/jupyter-source.js`
- Protocol architecture: {doc}`../../python_package/e_jupyter_hooks/05_architecture_message_routing_http_vs_postmessage`

Behavior:
- In Jupyter, Cellucid often runs inside an iframe and exchanges messages/events with a local server.
- The same loader architecture is used, but the source and event pathways differ.

---

## Loader layer (what files are fetched)

Most of the “what do we fetch?” logic lives in:
- `cellucid/assets/js/data/data-loaders.js`

Common artifacts:

- `dataset_identity.json`
  - dataset id, name, stats, and “what exists” metadata (embeddings, fields, vector fields, etc.)

- `obs_manifest.json`
  - obs field inventory (categorical/continuous), including any precomputed metadata needed for fast UI.

- `var_manifest.json` (optional)
  - gene-expression availability and gene index metadata.

- embeddings / points buffers
  - used to build positions arrays for the viewer.

- connectivity artifacts (optional)
  - `connectivity_manifest.json`
  - edge arrays or edge textures depending on format

:::{important}
Cellucid expects gzipped binary assets (e.g. `.bin.gz`) to be served as *raw gzipped bytes*.
The loader detects `.gz` by filename and decompresses in the browser using `DecompressionStream`.
If a server “helpfully” applies `Content-Encoding: gzip` to these files, you can get double-decompression failures.
See {doc}`03_build_run_and_deployment`.
:::

---

## Caching layers (what is cached where)

Cellucid uses multiple caches; knowing which one you’re hitting makes debugging much faster.

### 1) HTTP cache (browser)

Some loaders use fetch options like:
- `{ cache: 'force-cache' }`

This can speed up repeated loads, but can also make debugging confusing if you expect a changed file to reload immediately.

Debug tip:
- DevTools → Network → “Disable cache” while you debug.

### 2) In-memory caches in the data layer

Examples:
- `DataSourceManager` caches dataset metadata per source.
- Some sources cache manifests (`datasets.json`, metadata) after first load.

### 3) `DataState` field caches (bounded LRU)

`DataState` maintains LRU caches for loaded field arrays:
- obs fields (many small arrays)
- var fields (gene vectors; fewer but larger)

These caches are designed to prevent “scroll through genes for 10 minutes → tab uses 10GB”.

### 4) Community annotation caches (IndexedDB + localStorage index)

Community annotation stores:
- raw downloaded JSON file content (IndexedDB)
- a small path→sha index (localStorage) for “which files changed?”

Code:
- `cellucid/assets/js/app/community-annotations/file-cache.js`

---

## Debugging checklist (data loading)

When something “doesn’t load”, do the following in order:

1) **Identify the source type**
   - local-demo vs local-user vs remote vs github-repo vs jupyter

2) **Find the base URL**
   - In console:
     - `window._cellucidState` exists after bootstrap
     - The dataset base URL is logged by `main.js` when debug is enabled

3) **Watch the Network tab**
   - Look for:
     - `404` (wrong path)
     - CORS errors (no `Access-Control-Allow-Origin`)
     - HTML returned for a JSON/binary request (“Unexpected token `<`”)

4) **Confirm manifest parsing**
   - `obs_manifest.json` is required for most UI behavior.
   - If it is missing or invalid, loading terminates with that manifest error.

5) **Confirm point buffer sizes**
   - After load, check:
     - `window._cellucidState.pointCount`
     - `window._cellucidState.positionsArray.length === pointCount * 3`

6) **Check gzip support**
   - Cellucid requires native `DecompressionStream('gzip')` for `.gz` files.
   - A browser without that API receives a terminal capability error before the data request.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Dataset loads forever”

Likely causes (ordered):
1) A required manifest request is hanging or blocked (CORS/mixed content).
2) A binary file fetch is large and the server is slow (or throttled).

How to confirm:
- Network tab: which request is pending?
- Console: enable debug and look for the last log line.

Fix:
- For remote datasets: validate CORS + protocol (https vs http).
- For GitHub datasets: verify raw URLs are reachable from the target environment.
- For large H5AD files: choose remote server mode for read-only-backed access.
- For Zarr directories: choose server mode only when the Python process can
  materialize the store.

### Symptom: “Unexpected token `<`” during loading

Likely cause:
- The server returned an HTML document for a JSON request.

Fix:
- Fetch the URL in a new tab; confirm it’s real JSON and returns `404` when missing.

### Symptom: “Gzip decompression not supported”

Likely cause:
- The browser does not provide native `DecompressionStream('gzip')`.

Fix:
- Update to the current stable release of Chrome, Edge, Firefox, or Safari.
- Re-export without compression only when gzip transport is not required.

---

Next: {doc}`10_sessions_persistence_and_serialization` (how loaded state is saved/restored).
