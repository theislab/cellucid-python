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
| `cellucid_community_annotation_file_cache` | `files` | Raw annotation JSON downloaded from GitHub |

It is the only database. Marker-discovery results are cached in memory for the
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
| `cellucid-transcoded/<digest>-<random>` | The AnnData store an R object (`.rds`, `.rda`) or a MatrixMarket folder is converted into as it is read whole once, or a file a host serves only whole as downloaded (`data/anndata/transcode/space.js`); `<digest>` starts the SHA-256 of the source's identity (`cacheKey`: a link's signature and token query parameters left out), `manifest.json` (written last) holds the whole digest (never the identity: no address, id or credential at rest), the store revision, whether it is kept for later visits (else the tab it belongs to), every file's size and the CRC-32 of each file's page-checksum table (`<file>.sums`: a CRC-32 of every 16 KiB of the file, taken as it was written, `formats/io/page-checks.js`), `used` when it was last opened | Kept across visits when the identity names a version (a local file's modification time, a link's readable ETag or provider revision, else its readable `Last-Modified`, the `ETag` a viewer host lists a served file with): checked against its manifest and its checksum tables before it is reopened, and every read of it checked page by page (a damaged page fails the read as damaged, the manifest is removed, so the directory goes once nobody holds it and the next open reads the file again), removed when work needs the room under the real quota (after every kept gene index run nobody holds, least recently used first), or by **Clear** (Session › Open data › Stored). One known only by its size (a link with no version) is opened again by its own tab only, kept while that tab holds it (a Web Lock) and removed when none does (at startup, before a new conversion, after a cancelled one). Memory instead where the browser has no OPFS, or a page cannot write it (Safari before 26, for downloads) |
| `cellucid-transcoded/<digest>-<random>` (`kind: column-index`) | One run of open rows of a cell-major matrix's gene index, kept for a later opening of the same cells (`data/anndata/kept-column-index.js`): `<digest>` starts the SHA-256 of its group (the file's identity, the matrix, its shape, the cells' method, seed and settings, never their count); `d.bin` holds the index as built (entry positions, or rows and values, or a dense column store) with a CRC-32 of each gene's piece of each, the run's file rows and row pointers; the manifest its first open row, row and entry counts, the index revision (`KEPT_INDEX_REVISION`), the digest of the file's identity (`source`: the Session panel counts a file once) and, for a file whose identity does not name its content, its fingerprint (one seeded row per doubling band of its rows, with the CRC-32 of its entries) | Kept across visits when the file's identity names a version, else for its tab. Checked before use: manifest, file rows, row pointers against the file's (a mismatch removes it), column starts and checksum tables; the fingerprint re-read from the file unless the identity names the content (a provider checksum or hash, a commit, a Census release). Each gene's piece is checked against its CRC-32 as it is read; a damaged run is removed and built and kept again. Makes room only from other runs; removed before converted stores when work needs room; **Clear** removes it |
| `cellucid-column-index/<directory>` | A gene-column index of a row-major (CSR) matrix, spilled to disk when it outgrows its memory budget (`data/column-store-placement.js`); when the quota leaves less room than it needs, the kept copies no tab holds are removed for it first (gene index runs, then converted stores and downloads) | Held by a Web Lock while its reader lives; removed at startup and before a new index is written when no tab holds it |

### `sessionStorage` keys

Session storage is used for **short-lived secrets and per-tab identity**.

| Key | Purpose | Notes |
|---|---|---|
| `cellucid:github-app-auth:session` | GitHub access token **and** the signed-in user payload, in one record | Community annotations |
| `cellucid:google-auth:session` | Google access token (`drive.file`), its expiry and granted scope | Google Drive files; about an hour, renewed with one click |
| `cellucid:microsoft-auth:session` | Microsoft access token, refresh token, their expiries and granted scopes | OneDrive files; renewed without a click for 24 hours |
| `cellucid:community-annotations:tab-id:v1` | Unique tab id for scope locks | Used by scope-lock lease mechanism |

`sessionStorage` holds the app's three sign-ins, and closing the tab destroys
them.

#### The sign-in token model

Every sign-in goes through one module,
`cellucid/assets/js/data/sources/auth/token-store.js`: GitHub
(`app/community-annotations/github-auth.js`), Google (`data/sources/auth/google.js`)
and Microsoft (`data/sources/auth/microsoft.js`). The rules, each held by a
contract test (`tests/source-sign-in-contract.test.mjs`):

- A token lives in the tab's `sessionStorage` and in memory: it survives a
  reload of the tab and is gone when the tab closes. Nothing writes one to
  `localStorage`, IndexedDB, cookies, a URL, a session file, figure
  provenance, the console, analytics or an error message.
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
