# Build, run, and deployment

This page explains what “build” means for Cellucid, how to run it reliably, and what to watch for when deploying the web app (or hosting datasets).

## Key idea: Cellucid is a static ES-module web app

Cellucid is shipped as **static files** (HTML/CSS/JS + data assets).
There is no required bundler step to run it locally:

- `cellucid/index.html` starts the app from its boot script
  (`<script id="cellucid-page-boot">`, CELLUCID_PAGE_BOOT), which imports the
  entry once the document is parsed.
- On this machine (`localhost`, `127.0.0.1`, `[::1]`: a development server,
  the browser tests' server, the Python and R packages' servers) the entry is
  `assets/js/app/main.js`, and the browser loads the ES modules directly from
  the server. Every script URL names its bytes (`main.js?v=<first 12 hex
  digits of its SHA-256>`); an import map in `index.html` does the same for
  every import. A static server ignores the query, so a module you edit loads
  before you have run the generator.
- Anywhere else (a deployed site, a notebook served through a proxy) the entry
  is the **page bundle**: the same module graph built by esbuild (pinned in
  `package.json`) into `assets/page/`, minified with function and class names
  kept, the parts loaded on demand (`import()`: the analysis panel, the
  benchmark, the figure renderers) in chunks of their own, every file named by
  its content (`<name>.<content hash>.js`). The boot script asks for the entry
  and every file it imports at startup at once, from the head
  (`modulepreload`). A module that reads its own address (`import.meta`) is not
  bundled but imported at its versioned URL (`PAGE_EXTERNALS` in
  `scripts/build-page-bundle.mjs`): the worker scripts, the worker-start
  registry and `app/asset-url.js`.
- Workers follow the same rule: their sources on this machine, a bundle per
  worker in `assets/workers/` anywhere else.
- Stylesheets name what they import the same way: every `@import` and
  `url()` in `assets/css/` carries the version of its file (`?v=`), and each
  page preloads every stylesheet its stylesheet imports (one round trip, not
  one per `@import` level).
- The generator writes all of it (checklist item 5); the bundles are
  committed with the sources, so a host needs no build of its own.

That makes development simpler, but it also means **your server configuration matters**.

---

## Run locally (static server)

Recommended minimal server:

```bash
cd cellucid
python -m http.server 8000
```

Then open `http://localhost:8000`.

Why this is required:
- Browsers restrict ES module imports and `fetch()` behavior under `file://`.
- A real HTTP server makes your behavior match real deployments.

---

## Deploy the web app (static hosting)

Any static host works (GitHub Pages, Netlify, S3+CloudFront, Cloudflare Pages, nginx, Apache), as long as it satisfies the constraints below.

### Deployment checklist (web app files)

1) **Correct MIME types for JS modules**
   - `.js` must be served as JavaScript (e.g. `text/javascript`).
   - If your host serves `.js` as `text/plain`, module imports will fail.

2) **Don’t break module paths**
   - The app uses relative module imports under `assets/js/…`.
   - If you deploy under a subpath, verify that relative URLs still resolve.

3) **CSP expectations**
   - `cellucid/index.html` defines CSP in a `<meta http-equiv="Content-Security-Policy">`.
   - The policy is `default-src 'self'`. Every JavaScript dependency is vendored
     under `assets/external/` and served from the app's own origin; the only
     external script origins allowed are Google's tag manager and, loaded only
     when a Google Drive file needs them, Google's sign-in and picker
     (`accounts.google.com`, `apis.google.com`; the picker's frame comes from
     `docs.google.com`). Adding a CDN
     dependency requires editing this policy, which is the intended friction.
   - `script-src` grants `'wasm-unsafe-eval'` and nothing else that evaluates
     code: the readers compile the WebAssembly decoders
     (`assets/wasm/codecs.wasm`: zstd 1.5.6, lz4 1.10.0 and c-blosc 1.21.6,
     built from their pinned release tarballs by
     `scripts/build-wasm-codecs.mjs`), and JavaScript evaluation
     (`'unsafe-eval'`) stays refused. Only `index.html` runs the readers, so
     only it grants it; `tests/data-loading-contract.test.mjs` and
     `tests/browser/content-security-policy.spec.mjs` fail if `'unsafe-eval'`
     is granted or `'wasm-unsafe-eval'` is missing.
   - If you move CSP into an HTTP header (recommended for production), keep it consistent with required endpoints (datasets, link providers, Google/Microsoft sign-in, GitHub worker, analytics).
   - If you edit the inline JSON-LD `<script type="application/ld+json">`, the CSP hash in `script-src` must be updated or the browser will block it.

4) **Know what your hostname turns on**
   - Analytics initialises **only** on `cellucid.com`, `www.cellucid.com` and
     `theislab.github.io`. A deployment on any other hostname — including a
     staging copy of the production site — reports nothing, by design.
   - If you deploy an internal instance for sensitive data, that is the property
     you are relying on. Document it for your reviewers:
     {doc}`../o_accessibility_privacy_security/02_privacy_model`.

5) **A deploy never pairs files of two builds from a cache**
   - GitHub Pages answers every file with `Cache-Control: max-age=600`, and
     the headers cannot be changed; a static server without cache headers
     (`python -m http.server`) lets the browser guess a freshness of its own.
     A page could then run a module its browser had cached from the build
     before a deploy beside one fetched from the build after it, and stop on
     "does not provide an export named …".
   - So every script is requested at a URL that names its bytes:
     `<path>?v=<first 12 hex digits of its SHA-256>`. `index.html` carries an
     import map listing every `assets/**/*.js` at its versioned URL (a
     module's imports, relative or not, are resolved to a URL and then looked
     up there), and each `<script src>` of the shipped pages carries the same
     `?v=`. A file whose bytes did not change keeps its URL, so it stays
     cached across deploys; a changed one has a URL no cache has seen.
     `scripts/generate-web-asset-inventory.mjs` writes the map, the versions
     and the `script-src` hashes of each page's inline scripts (the map and the
     JSON-LD block) with the inventory; `npm test` refuses a tree where any is
     stale (`tests/web-asset-inventory-contract.test.mjs`), and
     `tests/browser/deploy-skew.spec.mjs` simulates deploys in all three
     engines.
   - **Workers.** No browser gives a worker an import map (measured
     2026-10-25 in Chromium, Firefox and WebKit), so a worker started from its
     source would ask for each of its imports by plain path, and could meet a
     cached file of the previous build. So the generator also builds every
     worker (the AnnData reader, decode, gunzip, spatial-tree, point-raster,
     figure-export, analysis and benchmark workers, and the two modules the
     decode workers import by URL) with esbuild (pinned in `package.json`)
     into one file, `assets/workers/<name>.<first 12 hex digits of its
     SHA-256>.js`, with a source map to the sources beside it, and writes
     `assets/workers/scripts.js`, which names each. A page served under any
     name but `localhost`, `127.0.0.1` or `[::1]` starts every worker from its
     bundle: one request at an address no other build uses, so a worker never
     mixes builds. A server on this machine (a development server, the
     browser tests' server, the Python and R packages' servers) serves the
     files as they are, and there workers start from their sources: nothing
     is built while developing, and `assets/js/workers/worker-script.js` is
     the only place the two differ (the URL a worker starts from). The
     bundles are committed beside the sources (GitHub Pages builds nothing),
     listed in the inventory with their source maps (so the Python and R web
     caches download and verify them like every other file) and allowed by
     the policy's `worker-src 'self'`; `npm test` refuses a bundle that is not
     exactly what the sources build (`tests/worker-bundles-contract.test.mjs`).
     A page that outlived its build asks for a bundle the host no longer has:
     the worker fails to start and the page says **Cellucid has been updated**
     with **Reload**, as for a part loaded on demand.
   - **The page bundle.** The same reasoning applies to the page's own
     graph on a deployed site: `assets/page/<name>.<esbuild content hash>.js`
     (built by `scripts/build-page-bundle.mjs`, written by the generator
     after the worker bundles, listed in the inventory with its source maps;
     `npm test` refuses a bundle that is not exactly what the sources build:
     `tests/page-bundle-contract.test.mjs`). The import map lists only the
     files whose name does not already name their bytes.
   - **The code cache.** GitHub Pages' `max-age=600` makes a visit ten
     minutes after the last revalidate every file. A deployed page registers
     `service-worker.js` (a root file, since a service worker controls only
     pages below its folder) once it has loaded, with its build id
     (`?build=`), and names the files it loaded. The worker keeps every
     same-origin file under `assets/` whose address names its bytes (a `?v=`
     file only when its bytes hash to its version, since the CDN may answer a
     new version's address with old bytes in the minutes after a deploy; a
     bundle file by its name) and answers later requests for them from the
     cache. Nothing else goes through it (where the browser supports it,
     `addRoutes` sends everything else straight to the network; Firefox 153
     has no `addRoutes`, so there those requests reach the worker, which
     leaves them to the network). Caches are per
     build; the current build's and the newest other's are kept. A deploy
     cannot make it serve stale code: a kept file's address names its bytes,
     so a new build's changed file is a new address; the new page registers
     the worker under its own `?build=`, the new worker takes over at once
     (`skipWaiting`) and, as it activates, removes every build's cache but its
     own and the newest other's (a tab still open on the previous build finds
     its files there). It is never
     registered on this machine. Specs run without service workers
     (`playwright.config.mjs`); `tests/browser/code-cache.spec.mjs` and
     `tests/code-cache-contract.test.mjs` cover it. Playwright's Firefox runs
     only the spec's first visit: under its automation protocol (Juggler) a
     navigation the worker leaves to the network commits without a
     controller, so its next visit never starts under the worker. Measured
     2026-10-31 with a 20-line worker as well as the app: the same Firefox 153
     build, run with no automation (a page that drove its own second visit)
     or through Playwright's WebDriver BiDi channel (`channel: 'moz-firefox'`,
     same executable), starts the next visit under the worker and took none
     of the 22 startup files and the stylesheet from the network (Chromium
     and WebKit none either; Firefox under Juggler all 23). A Juggler page is controlled only when
     the worker answers the navigation itself, which the code cache does not
     do for a test harness's sake: every page load would wait on the worker.
     Firefox reports `workerStart` 0 for a file its worker answered, so the
     spec's witness for "answered by the worker" holds in Chromium and WebKit
     only.
   - What a content address cannot reach, and what the app does about it:
     - **A page that outlived its build.** GitHub Pages' CDN ignores the query
       string (a never seen `?v=` on an existing file is answered
       `x-cache: HIT`, measured 2026-10-23), and the origin keeps only the
       current files. A tab open across a deploy, or an `index.html` still in
       the browser's cache for up to ten minutes, asks for its own build's
       URLs and receives the next build's bytes. At startup, a module that
       fails to download, parse, link or run shows the **Cellucid could not
       load** card (`ui/core/code-load-watch.js`), which asks the server for
       its build id and says **Cellucid has been updated** when that is not
       the page's; **Reload** loads one build. A part loaded on demand later
       (`app/app-code.js` `loadAppCode`: the benchmark panel, the figure
       export renderers, Plotly) that fails says the same in a notification
       with **Reload**.
   - Whether GitHub Pages purges its CDN on a deploy is not measured. If it
     does not, a page loaded within ten minutes of a deploy can receive a
     file's previous bytes at its new URL; a mismatch that fails to link is
     shown as above.
   - Dataset stores are not affected in the browser: every request into a
     Zarr store bypasses the HTTP cache (see
     {doc}`09_data_loading_pipeline_and_caching`, "HTTP cache").

6) **Catch-all route rewrites**
   - Cellucid is a single page, but it is not using a framework router that needs wildcard rewrites.
   - Do not serve `index.html` for missing data files; that turns a missing-file
     response into a misleading JSON parse error.

---

## Host datasets safely

Many Cellucid workflows involve hosting datasets separately from the app (e.g. a dataset catalog on GitHub Pages, an institutional HTTP server, or a CDN).
This section highlights server behaviors that commonly break loading.

### 1) CORS: allow the app origin to fetch your dataset files

If you load a dataset hosted at `https://example.org/my_datasets/` from `https://www.cellucid.com`, your dataset host must send CORS headers like:
- `Access-Control-Allow-Origin: https://www.cellucid.com`

For a development or internal server, you may choose to allow `*` (public datasets only).

This requirement has no alternate browser-side transport: if the dataset host
does not authorize the app origin, host the files on an endpoint whose CORS
policy you control.

### 2) Serve dataset bytes as they are, and honour `Range`

Cellucid reads an HDF5 file (`.h5ad`, `.h5mu`, `.loom`, 10x `.h5`), a Zarr ZIP
or a Zarr store's chunks in place with HTTP range requests
(`Range: bytes=a-b`), and decompresses chunks itself (using
`DecompressionStream`). R objects and MatrixMarket files are read once from
the start (`access: 'stream'` in `data/input-formats.js`).

That means:
- The host must answer a range request with `206 Partial Content` and exactly
  those bytes; a host that ignores `Range` makes every read download the whole
  file.
- The host should **not** apply `Content-Encoding: gzip` to these files: byte
  ranges then address the transformed stream, not the file.
- The file's size must be known (`HEAD` with `Content-Length`, or the
  provider's own listing).

### 3) Avoid returning HTML for missing JSON/binary files

Some static hosts return a friendly HTML page for 404s.
From the app’s perspective, that looks like:
- “JSON parse failed” (`datasets.json`, Zarr metadata)
- “Unexpected token `<`”

A Zarr store folder on a web host is read key by key, and an absent key must be
a real `404`: Zarr reads a missing chunk as fill values.

Always verify that a missing asset is a real `404` (not `200` with HTML).

### 4) `.cellucid-session` files and proxies

Session bundles are binary `.cellucid-session` files.
Some proxies compress or transform unknown binary types.

If you host sessions:
- Prefer serving them as raw bytes (no transformations).
- If your host compresses them, be aware `Content-Length` may not reflect decompressed size (Cellucid treats it as a hint, but extreme misconfiguration can still break loads).

See {doc}`10_sessions_persistence_and_serialization`.

---

## Troubleshooting (deployment failures)

### Symptom: “Module script failed to load”

Likely causes (ordered):
1) Wrong MIME type for `.js`.
2) You deployed under a subpath and module import URLs are wrong.
3) A caching layer is serving a module of another build (checklist item 5).

How to confirm:
- DevTools → Network → click the failing `.js` request → check Response Headers and preview.

Fix:
- Adjust static host configuration to serve `.js` as JavaScript.
- A card or notification saying **Cellucid has been updated**: choose
  **Reload** (checklist item 5). If a module still does not fit after a
  reload in your own checkout, run `node scripts/generate-web-asset-inventory.mjs`:
  an `index.html` whose import map names old versions keeps asking for them,
  and the browser may answer them from its cache.

### Symptom: “Dataset loads on my machine but not for others”

Likely causes:
- CORS headers missing on your dataset host.
- You accidentally relied on local filesystem access or a local server.

Fix:
- Host the dataset on a real HTTP server with explicit CORS.
- Test from a clean browser profile/incognito.

---

Next: {doc}`04_configuration_env_vars_and_feature_flags`.
