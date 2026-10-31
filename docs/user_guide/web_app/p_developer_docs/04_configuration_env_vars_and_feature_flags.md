# Configuration, environment variables, and feature flags

Cellucid is primarily a **static browser app**, so most configuration is done via:

- **URL parameters** (deep-linking to datasets / connections / annotation repos)
- **browser storage** (preferences like theme; safe caches; debug toggles)
- **explicit global overrides** (developer-only `window.__…` flags for local testing)

There is no general “dotenv environment variable” configuration for the frontend, because there is no build step required to run it.

:::{note}
Privacy/security expectations for storage are documented for users in:
{doc}`../o_accessibility_privacy_security/02_privacy_model`.
This page focuses on **developer-facing** “what keys exist and where they are read/written”.
:::

---

## URL parameters (deep links)

### Data loading params (dataset/source selection)

Handled by:
- `cellucid/assets/js/app/main.js`
- `cellucid/assets/js/app/url-state.js` (keeps the URL in sync with UI connections)

| Param | Example | Meaning |
|---|---|---|
| `dataset` | `?dataset=suo` | Selects a local-demo dataset by id; on a page served by the Python server (which auto-connects to its own origin) a bare `dataset` names a server dataset instead |
| `source` | `?source=local-demo&dataset=suo` | Selects the source type for `dataset` (default: `local-demo`). Picking a sample writes this explicit form, so the URL reloads as the sample on any origin |
| `remote` | `?remote=https://data.example.org&dataset=study-a` | Connects to a remote data server; `dataset` is required when its manifest contains multiple datasets |
| `anndata` | `?remote=http://127.0.0.1:8765&anndata=true` | Says the source is a Python server serving one dataset (`cellucid serve`, `serve_anndata()`; any format of `data/input-formats.js`); startup requires a server source and refuses a static origin. The browser reads the file or store the server's `_cellucid/datasets` names |
| `data` | `?data=gh:owner/repo/path&dataset=study-a` | Opens a link through a provider of `data/sources/registry.js`: `gh:owner/repo/path` (GitHub), `cellxgene:<dataset id>`, `geo:<GSE or GSM>[/<file>]`, `zenodo:<record>[/<file>]`, `hf:<user>/<repo>[@<revision>][/<path>]` (Hugging Face), `figshare:<item>[/<file>]`, `gdrive:<id>` (Google Drive, sign-in), `onedrive:<drive>/<item>` (OneDrive, sign-in), a Dropbox link, or a web address of a file, a store folder or a dataset catalog folder. `dataset` is required when the link holds several datasets ({doc}`../b_data_loading/14_open_from_link`) |
| `sampleCatalog` | `?sampleCatalog=http://localhost:9000/catalog/` | Replaces the configured sample catalog (`<meta name="cellucid-sample-catalog">`) for this page: one absolute folder URL ending in `/` that holds a `datasets.json` (`data/data-source.js`) |
| `cells`, `sample`, `seed`, `block`, `by`, `alloc`, `min`, `embedding`, `grid` | `?data=…&cells=10%&sample=blocks&seed=7` | The subsample the dataset opens with (`app/subsample-url.js`): `cells` a share (`10%`) or an exact count (`250000`), `cells=100%` every cell; `sample` `uniform`, `blocks`, `stratified` (with `by=<obs field>`, a categorical or text field, optional `alloc=balanced`, `min=<n>`) or `spatial` (optional `embedding=<obsm key>`, `grid=<n>`), by default the file's own (`defaultSubsampleMethod`: `uniform` for a file in a verified priority order, `blocks` for any other); `seed` a uint32 (default the file's priority seed, else 0); `block` cells per block (default from the file's layout, `defaultBlockSize`); `cells` at least 2 when a count. Overrides fitting to the device. Parsed strictly: a malformed value or a parameter of another method is refused by name and the dataset opens as if none were given. Written back as the exact count, method, seed and settings of a subsample chosen in the Session panel, never one fitted to the device; never the cells themselves ({doc}`../b_data_loading/21_subsampled_loading`) |
| `viewer=embedded` | `?viewer=embedded` (nothing else) | The page a notebook widget frames (Python `show_anndata()`, R `cellucid_view()` in Jupyter): it receives the viewer's id and token and a `MessagePort` from its parent frame in a message, never in its address, and reads the data and runs the viewer channel over that port (`data/viewer-embed.js`) |
| `viewer`, `viewerId`, `viewerToken` | `?viewer=true&viewerId=…&viewerToken=…` | The page of a browser-tab viewer a host program opened (Python `AnnDataViewer`, R `cellucid_view()` outside a notebook), on that program's own server; the viewer connects its channel (`/_cellucid/ws`) with the id and token. Written by the host, never by the app (`data/viewer-host-source.js`, `data/viewer-protocol.js`) |

Notes:
- A remote or linked dataset list containing exactly one dataset selects that sole
  dataset when `dataset` is omitted. A list containing multiple datasets
  requires `dataset=<exact-dataset-id>`; Cellucid does not choose an arbitrary
  catalog entry.
- **Mixed content:** when the viewer origin is HTTPS (e.g. `https://www.cellucid.com`), browsers will block `remote=http://127.0.0.1:<port>` (HTTPS page fetching HTTP). For a local `cellucid serve`, open the exact server-backed Viewer URL it prints (`http://127.0.0.1:<port>/?anndata=true`), served by `cellucid-python` itself.
- For **local-user** (browser file picker), the URL is intentionally kept clean (no local paths in the URL).
- `url-state.js` uses `history.replaceState()` so UI changes don’t spam browser history.
- Publishing a new data source replaces every data-source param (`dataset`, `source`, `remote`, `data`, `anndata`) and every subsample param; `anndata=true` describes the Python AnnData server the page was opened from, so it is not carried over to another source.

### Community annotation param

Handled by:
- `cellucid/assets/js/app/url-state.js`
- `cellucid/assets/js/app/community-annotations/REPO_SETUP.md` (behavior + sharing)

| Param | Example | Meaning |
|---|---|---|
| `annotations` | `?annotations=owner/repo` | Pre-selects the annotation repo for the current dataset |
| `annotations` | `?annotations=owner/repo@branch` | Same, but pins to a branch |

:::{important}
Tokens are **never** stored in the URL.
Only the repo reference may be stored (and can also be persisted in `localStorage` for convenience).
:::

### Debug / diagnostics params

| Param | Example | Meaning |
|---|---|---|
| `debug` | `?debug=1` | Enables analysis-module debug logging (see `cellucid/assets/js/app/analysis/shared/debug-utils.js`) |

:::{note}
There is no analytics debug parameter. Analytics is decided entirely by the
page's hostname in `cellucid/assets/js/app/ui/core/ga-init.js`: it initialises on
`cellucid.com`, `www.cellucid.com` and `theislab.github.io`, and nowhere else.
The resulting `window.cellucidAnalyticsEnabled` flag is defined non-writable and
non-configurable, so it cannot be flipped from the console or by a query
parameter — by design. Check it to find out whether the current origin reports.
:::

---

## Browser storage keys (preferences + safe caches)

### `localStorage` keys

| Key | Purpose | Where used |
|---|---|---|
| `CELLUCID_DEBUG` | Enables app-wide debug logger | `cellucid/assets/js/utils/debug.js` |
| `cellucid_theme` | Theme preference (`light`/`dark`) | `cellucid/assets/js/utils/theme-manager.js`, `cellucid/assets/js/app/ui/core/theme-init.js` |
| `cellucid_viewer_background` | Viewer background mode (`white`, `grid`, `studio`, `atelier`, `abyss`, `aurora`, `black`, `grid-dark`, `studio-dark`, `nebula`, `abyss-dark`, `aurora-dark`; the enum is `cellucid/assets/js/rendering/viewer-backgrounds.js`) | `cellucid/assets/js/app/ui/core/theme-init.js`, viewer/render controls |
| `cellucid_antialias` | Antialiasing preference — `auto` (also what an absent key means), `on`, or `off`. `auto` resolves from the dataset's cell count on every dataset publication | `cellucid/assets/js/app/ui/core/antialias-preference.js`, `cellucid/assets/js/app/ui/modules/render-controls.js` |
| `cellucid_last_quote_index` | Welcome modal “quote rotation” state | `cellucid/assets/js/app/ui/onboarding/welcome-modal.js` |
| `debug` | Analysis-module debug toggle | `cellucid/assets/js/app/analysis/shared/debug-utils.js` |
| `cellucid:cellxgene-census-release` | The CELLxGENE Census release this browser read last, so a `cellxgene:` link probes its file in that release while the release list loads | `cellucid/assets/js/data/sources/providers/cellxgene.js` |
| `cellucid:opening-cost:v2` | What showing a cell costs this page, measured (cell counts and milliseconds of earlier openings), so a progressive opening sizes its first view; kept in memory and here as a convenience | `cellucid/assets/js/data/opening-cost.js` |
| `cellucid:analysis-cost:v2` | What each analysis kind costs per cell on this device (`fixedMs + perCellMs × cells`), so a run is sized before it starts | `cellucid/assets/js/app/analysis/shared/progressive-cells.js` |
| `cellucid:host-link:v1` | What each host delivered on earlier visits, by origin (`scheme://host[:port]`: never a path, a dataset, a query or a token): the round trip, the round trip under load (what a read's time is priced at), the rate its reads sustained (bytes over the time a request was in flight), the bytes they rest on, when last measured. At most 32 hosts, the least recently measured dropped. An opening sizes its first gene from it until its own reads measure the link's bandwidth or rest on as many bytes ({doc}`09_data_loading_pipeline_and_caching`, layer 7); cleared from the Stored dialog (`forgetHostLinks`) | `cellucid/assets/js/data/host-link-memory.js` |
| `cellucid:host-pace:v1` | The pace each host that blocked this network was found to bear, by origin (never a path, a dataset, a query or a token): the least milliseconds between two requests (half the requests of the minute that met the block: a host that answered most of them and then sent nothing a page may read to two rounds, as Zenodo's "unusual traffic" page), the pace that block taught (`learntMs`), when it was and how many blocks the host has met. At most 32 hosts. A host that blocked once and has not refused for a day starts each later visit a step looser (a quarter more requests a minute, `PACE_RECOVERY_STEP`); a refusal at a loosened pace snaps it back to the learnt pace, and a host that blocks again (or twice in a visit) stays paced. An opening's requests to the host start at it, from the page and its reader workers ({doc}`09_data_loading_pipeline_and_caching`); cleared from the Stored dialog with the hosts' speeds (`forgetHostLinks`) | `cellucid/assets/js/data/host-link-memory.js`, `cellucid/assets/js/formats/io/hosts.js` |

Community annotation (localStorage-backed):

| Key / prefix | Purpose | Notes |
|---|---|---|
| `cellucid:community-annotations:repo-map` | “datasetId + user → chosen annotation repo” mapping, with its branch mode | Convenience only; no secrets |
| `cellucid:community-annotations:auto-pull:v1` | Auto-pull preference | Convenience only |
| `cellucid:community-annotations:lock:*` | Cross-tab “single editor per scope” lease locks | Safety feature to prevent silent data loss |
| `cellucid:community-annotations:cache:*:files:shas` | SHA index for cached GitHub files | Content is in IndexedDB; this is only the small index |
| `cellucid:community-annotations:cache:*:session` | **The working annotation session**: profile, votes, comments, suggestions, deleted suggestions, moderation merges, last sync time | **Dataset-derived content, not a preference.** Scope is datasetId + repo@branch + GitHub user id |

:::{warning}
The `:session` entry is written to **`localStorage`**, so it survives closing the
tab. The helper that builds its key is named `toSessionStorageKey`, which refers
to the *annotation session*, not to `sessionStorage`. Do not infer a lifetime
from that name; earlier documentation did and got it wrong.

Because it holds annotation content rather than preferences, it belongs in the
privacy model: {doc}`../o_accessibility_privacy_security/02_privacy_model`.
:::

### IndexedDB databases

| Database | Store | Contents |
|---|---|---|
| `cellucid_community_annotation_file_cache` | `files` | Raw annotation JSON downloaded from GitHub (`app/community-annotations/file-cache.js`) |
| `cellucid-decoders` | `decisions` | One record: which decoder decodes each codec family and chunk size class fastest in this browser (`deflate:15` → `fflate`), under a version naming the WebAssembly module's digest, the decoder lists, the user agent and the core count (`formats/io/decoder-memory.js`); read by the reader thread's decode pool so only a first visit measures; listed and forgotten by the Stored dialog |
| `cellucid-highlight-drafts` | (one key–value store) | Per dataset id, the highlight pages and groups made on it, each group's cells as file rows, and the dataset's name (`app/highlight-draft.js`); listed and removed by the Stored dialog |
| `cellucid-local-files` | `handles` | Per dataset id, the browser's handle of a file or folder opened from this computer, its name, kind, size and modification time (a folder's file count and digest), and when it was opened; Chromium-based browsers only (`app/local-file-handles.js`); listed and removed by the Stored dialog |

Marker-discovery results are cached in memory for the
dataset load they were computed on and never written to IndexedDB
(`cellucid/assets/js/app/analysis/genes-panel/marker-cache.js`). When IndexedDB
is unavailable, Community Annotation Pull fails with a message and leaves the
annotation view as it was.

Code: `cellucid/assets/js/app/community-annotations/file-cache.js`.

### Origin-private file system (OPFS)

The AnnData reader worker writes **dataset-derived** bytes to the
origin-private file system, never anything else:

| Directory | Contents | Lifetime |
|---|---|---|
| `cellucid-transcoded/<digest>-<address>-<random>` | The AnnData store an R object (`.rds`, `.rda`) or a MatrixMarket folder is converted into as it is read whole once, or a file a host serves only whole as downloaded (`data/anndata/transcode/space.js`); `<digest>` starts the SHA-256 of the source's identity (`cacheKey`: a link's signature and token query parameters left out), `<address>` the first 12 hex digits of the SHA-256 of its address (the identity without its version, `addressOf`; `local` for a file on this computer), `manifest.json` (written last) holds the whole digest (never the identity: no address, id or credential at rest), `fileKey` (the digest of the file's identity) and `address` (the digest of its address, null for a local file), `about` (the file's name and host, `keptAbout`, for the Stored dialog: never a link, path or query), the store revision, whether it is kept for later visits (else the opening it belongs to), every file's size and the CRC-32 of each file's page-checksum table (`<file>.sums`: a CRC-32 of every 16 KiB of the file, taken as it was written, `formats/io/page-checks.js`), `used` when it was last opened | Kept across visits under the identity's version: a content name where the source has one (a Census build, a checksum, a commit or LFS oid, a Drive revision, a OneDrive cTag), else the host's validators (`writeVersion`: `written:<size> modified:<Last-Modified> etag:<strong ETag>`, each the host sends), a served file's strong ETag, a local file's size and modification time; with none of them (`size:<n>`) nothing is kept across openings, and such an opening removes nothing kept. A dataset's address names a GitHub or Hugging Face ref (`gh:owner/repo@<ref>/path`, `hf:repo@<revision>/path`), and a link is taken without its credentials (`data/sources/link-credentials.js`). Opening any version of a remote file first forgets what is kept of its other versions, every kind (`forgetOtherVersions`: the manifest removed, the directory once no tab holds it); a file that changes while open forgets every version (`forgetKeptFile`). Checked against its manifest and its checksum tables before it is reopened, and every read of it checked page by page (a damaged page fails the read as damaged, the manifest is removed, so the directory goes once nobody holds it and the next open reads the file again), removed when work needs the room under the real quota (after every kept gene index run and downloaded part nobody holds, least recently used first), or from the Stored dialog (Session › Open data › Stored › Manage…: `removeKeptStores`). One of a file with no version is used only by the opening that wrote it and that opening's readings again with other cells (the reader worker names a tab per opening), never by a later opening, in this tab or another; it is kept while a tab holds it (a Web Lock) and removed when none does (at startup, before a new conversion, after a cancelled one). Memory instead where the browser has no OPFS, or a page cannot write it (Safari before 26, for downloads) |
| `cellucid-transcoded/<digest>-<address>-<random>` (`kind: column-index`) | One run of open rows of a cell-major matrix's gene index, kept for a later opening of the same cells (`data/anndata/kept-column-index.js`): `<digest>` starts the SHA-256 of its group (the file's identity, the matrix, its shape, the cells' method, seed and settings, never their count); `d.bin` holds the index as built (entry positions; rows and values, also for a structure built over a link once every value of its rows is held, bought with the first gene or read whole since, when the run is kept again with them and the positions-only run it replaces forgotten; or a dense column store) with a CRC-32 of each gene's piece of each, the run's file rows and row pointers; the manifest its first open row, row and entry counts, the index revision (`KEPT_INDEX_REVISION`), the digest of the file's identity (`source`: the Session panel counts a file once), `fileKey`, `address`, `about` (the file's name and host, the matrix and the cells' order, for the Stored dialog) | Kept across visits under the file's version (the rule above), else for its opening. Checked before use: manifest, file rows, row pointers against the file's (a mismatch removes it), column starts and checksum tables. Each gene's piece is checked against its CRC-32 as it is read; a damaged run is removed and built and kept again. Makes room only from other runs and downloaded parts; removed before converted stores when work needs room; the Stored dialog removes it (all of one file's runs for one cell order are one item) |
| `cellucid-transcoded/<digest>-<address>-<random>` (`kind: parts`, `journal: true`) | The parts of a remote file (or the objects of a Zarr store folder on a web host) one reader fetched, kept as they arrived so a later opening of the same version downloads none of them (`formats/io/kept-parts.js`, a journal: `space.js` `openKeptJournal`): `<digest>` starts the SHA-256 of `{parts: identity, size}`; the manifest is written **first** (kept from the first byte on) with the raw-bytes revision (`KEPT_BYTES_REVISION`, also a whole download's), `size` (the file's exact size; null for a folder), `confirms` (a source whose reads name no version), `source`, `fileKey`, `address`, `about` (the file's name and host); `d.bin` holds the bytes as they arrived, appended through a synchronous access handle; `j.bin` the records, each appended after its bytes and framed by its length and CRC-32: a range (object key, `''` for a file; start in the object; offset and length in `d.bin`; the CRC-32 of every 16 KiB page of the piece), what an object's host said of it (its size and validator, the ETag or Last-Modified), an object of a folder its host answered as absent (404, 410: not asked for again under that version; a folder's consolidated metadata is kept like any object), or a confirmation that the provider still vouches for every range before it. A tab closed at once or a crash loses at most the piece being written; a record cut short or not matching its CRC-32 ends the journal there. The writer holds an exclusive lock `<directory>#writing` while it writes and closes the journal for good once its reads go quiet for the link's retransmission timeout, when the page is hidden or left, and when its reader closes (the next piece begins a new journal); another reader reads a journal only once that lock is free | Kept across visits under the file's version (the rule above); a file with no version keeps no parts, nor does a source on this machine or a notebook kernel (a loopback address, `cellucid serve`, the widget's host: `sameMachine`) or a store a Python or R server serves (`keptParts: false`), nor a Zarr store folder whose version names no content (a folder on a plain host). Checked before use: manifest, revision, every record whole and every range inside the file and `d.bin`, only the confirmed ranges of a source that needs confirming; every read checks the pages it touches (a damaged page forgets the journal and the bytes are fetched again). An object of a folder (read at a commit) fetched again with another validator than the one kept is a changed store: the read is refused and every kept copy forgotten. Makes room only from other parts and gene index runs, least recently used first; removed before converted stores when work needs room; the page holds each journal its tab used for the tab's life; the Stored dialog lists one item per file (*Downloaded parts*, with the bytes it saves) |
| `cellucid-column-index/<directory>` | A gene-column index of a row-major (CSR) matrix, spilled to disk when it outgrows its memory budget (`data/column-store-placement.js`); when the quota leaves less room than it needs, the kept copies no tab holds are removed for it first (gene index runs, then converted stores and downloads) | Held by a Web Lock while its reader lives; removed at startup and before a new index is written when no tab holds it |

### `sessionStorage` keys

Session storage is used for **short-lived secrets and per-tab identity**.

| Key | Purpose | Notes |
|---|---|---|
| `cellucid:github-app-auth:session` | GitHub access token, its expiry **and** the signed-in user payload, in one record | Community annotations; a GitHub App user token that expires after 8 hours, forgotten then; Sign out revokes it through the worker (`/auth/revoke`) |
| `cellucid:google-auth:session` | Google access token (`drive.file`), its expiry and granted scope | Google Drive files; about an hour, renewed with one click; Sign out revokes it at `oauth2.googleapis.com/revoke` |
| `cellucid:microsoft-auth:session` | Microsoft access token, its expiry and granted scopes (never the refresh token) | OneDrive files; about an hour. The refresh token lives in the page's memory only and renews the access without a click for up to 24 hours within the page; after a reload one click signs in again. Microsoft offers no revocation |
| `cellucid:community-annotations:tab-id:v1` | Unique tab id for scope locks | Used by scope-lock lease mechanism |

`sessionStorage` holds the app's three sign-ins, each a short-lived access
token; closing the tab destroys them. Browsers write `sessionStorage` to the
profile on disk for session restore, which is why nothing long-lived (no
refresh token) is kept there.

#### The sign-in token model

Every sign-in goes through one module,
`cellucid/assets/js/data/sources/auth/token-store.js`: GitHub
(`app/community-annotations/github-auth.js`), Google (`data/sources/auth/google.js`)
and Microsoft (`data/sources/auth/microsoft.js`). The rules, each held by a
contract test (`tests/source-sign-in-contract.test.mjs`):

- An access token lives in the tab's `sessionStorage` and in memory: it
  survives a reload of the tab until it expires, and is gone when the tab
  closes. Only short-lived access tokens are stored: Microsoft's refresh token
  stays in memory, and GitHub's never leaves the worker. Sign out revokes the
  token where the provider allows it (GitHub through the worker, Google at its
  revocation endpoint) and says so in a notification when it could not
  (Microsoft, an unreachable provider, a worker deployed before
  `/auth/revoke`). Nothing writes a token to `localStorage`, IndexedDB,
  cookies, a URL, a session file, figure provenance, the console, analytics
  or an error message.
- No client secret exists. Google uses the token model of Google Identity
  Services (a popup, no refresh token); Microsoft uses the authorization code
  flow with PKCE for a single-page application, redeemed in the page. The
  click opens the popup empty, on the page's own origin, and the page then
  sends it to Microsoft: Safari hands a popup opened at
  `login.microsoftonline.com` to a Mac's single sign-on extension (Microsoft's,
  from Company Portal), which signs in in a hidden window and closes it. The
  popup returns to `assets/auth/callback.html`, which posts the code over a
  `BroadcastChannel` and closes once the page answers that it has it; a popup
  found closed before then offers the step again, and the page keeps
  listening until the step ends (a popup cut off by a
  Cross-Origin-Opener-Policy reads as closed while still open). The public client IDs, API key, project
  number and the origins registered with Google (`registeredOrigins`) are in
  `assets/js/data/sources/provider-config.js`; empty values mean "not
  configured on this site yet". Where Google or Microsoft cannot answer the
  page (another origin, a plain `http://` page of another computer, a viewer
  embedded in another site, `data/sources/auth/page-context.js`), Cellucid
  says so before any request or popup.
- A file read behind a sign-in names only `{provider, ref}` in its resource
  (`data/sources/access.js`). The AnnData reader worker asks the page for the
  URL and headers to read from (`{event: 'source-access'}`), holds them in
  memory, and shares one renewal across all its parallel reads
  (`httpByteSource`'s `access` in `formats/io/byte-source.js`). Sessions and
  `?data=` links record the provider and file ids (`gdrive:<id>`,
  `onedrive:<drive>/<item>`), never a token or a download address.
- A step the person must take (sign in, sign in again, grant more access)
  waits in `data/sources/auth/sign-in.js` and is shown under the address bar
  and as a notification (`app/source-sign-in-view.js`); its button runs the
  sign-in inside the click, because a popup opens only from one.

---

## Global overrides (developer-only)

These are intended for local development and troubleshooting only.

| Global | Meaning | Where used |
|---|---|---|
| `window.CELLUCID_DEBUG = true` | Alternative way to enable the main debug logger | `cellucid/assets/js/utils/debug.js` |
| `window.__CELLUCID_DEBUG__ = true` | Enables analysis-module debug logging | `cellucid/assets/js/app/analysis/shared/debug-utils.js` |
| `window.__CELLUCID_GITHUB_WORKER_ORIGIN__` | Override GitHub auth worker origin (local dev only) | `cellucid/assets/js/app/community-annotations/github-auth.js` |
| `window._cellucidViewer`, `window._cellucidState` | The live viewer and state, for console debugging | `cellucid/assets/js/app/main.js` |
| `window._cellucidBenchmarkHarness` | The configuration-matrix benchmark harness, published when the Benchmark panel is first opened | `cellucid/assets/js/app/main.js` |
| `window.cellucidAnalyticsEnabled` | Read-only: whether this origin reports analytics | `cellucid/assets/js/app/ui/core/ga-init.js` |

:::{important}
For security, `window.__CELLUCID_GITHUB_WORKER_ORIGIN__` is rejected on non-local hosts unless it matches the compiled-in default origin.
This is designed to prevent accidental token exfiltration via an untrusted proxy.
:::

---

## Troubleshooting (when “configuration” causes weird behavior)

### Symptom: “The app behaves strangely after I changed something weeks ago”

Likely causes:
- Stale `localStorage` preferences or caches.

How to confirm:
- DevTools → Application → Local Storage / Session Storage.

Fix (safe reset):
- Clear `localStorage` keys that start with `cellucid` or `CELLUCID_` (theme/debug).
- For community annotation debugging, you may also clear `cellucid:community-annotations:*` keys.

### Symptom: “Community annotation sign-in works locally but fails on production”

Likely causes:
- A `window.__CELLUCID_GITHUB_WORKER_ORIGIN__` override is being used on a non-local host (blocked by design).

Fix:
- Remove the override and use the compiled-in worker origin, or update the compiled-in default for your deployment (see `cellucid/assets/js/app/community-annotations/REPO_SETUP.md`).

---

Next: {doc}`05_app_architecture_overview`.
