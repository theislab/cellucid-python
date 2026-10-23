# Build, run, and deployment

This page explains what “build” means for Cellucid, how to run it reliably, and what to watch for when deploying the web app (or hosting datasets).

## Key idea: Cellucid is a static ES-module web app

Cellucid is shipped as **static files** (HTML/CSS/JS + data assets).
There is no required bundler step to run it locally:

- `cellucid/index.html` loads `assets/js/app/main.js` via `<script type="module">`.
- The browser loads ES modules directly from the server.

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
     external script origin allowed is Google's tag manager. Adding a CDN
     dependency requires editing this policy, which is the intended friction.
   - `script-src` grants no code evaluation: neither `'unsafe-eval'` nor
     `'wasm-unsafe-eval'`. The app compiles no WebAssembly (`.h5ad` and Zarr
     files are read by its own JavaScript readers), and
     `tests/data-loading-contract.test.mjs` fails if either is granted.
   - If you move CSP into an HTTP header (recommended for production), keep it consistent with required endpoints (datasets, GitHub worker, analytics).
   - If you edit the inline JSON-LD `<script type="application/ld+json">`, the CSP hash in `script-src` must be updated or the browser will block it.

4) **Know what your hostname turns on**
   - Analytics initialises **only** on `cellucid.com`, `www.cellucid.com` and
     `theislab.github.io`. A deployment on any other hostname — including a
     staging copy of the production site — reports nothing, by design.
   - If you deploy an internal instance for sensitive data, that is the property
     you are relying on. Document it for your reviewers:
     {doc}`../o_accessibility_privacy_security/02_privacy_model`.

5) **Catch-all route rewrites**
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

### 2) Serve AnnData bytes as they are, and honour `Range`

Cellucid reads an `.h5ad`, a Zarr ZIP or a Zarr store's chunks in place with
HTTP range requests (`Range: bytes=a-b`), and decompresses chunks itself
(using `DecompressionStream`).

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
3) A caching layer is serving stale JS modules.

How to confirm:
- DevTools → Network → click the failing `.js` request → check Response Headers and preview.

Fix:
- Adjust static host configuration to serve `.js` as JavaScript.
- Bust caches (new filenames, or cache-control headers) for JS.

### Symptom: “Dataset loads on my machine but not for others”

Likely causes:
- CORS headers missing on your dataset host.
- You accidentally relied on local filesystem access or a local server.

Fix:
- Host the dataset on a real HTTP server with explicit CORS.
- Test from a clean browser profile/incognito.

---

Next: {doc}`04_configuration_env_vars_and_feature_flags`.
