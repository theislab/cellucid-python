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
   - `cellucid/index.html` contains a dev-phase CSP in a `<meta http-equiv="Content-Security-Policy">`.
   - If you move CSP into an HTTP header (recommended for production), keep it consistent with required endpoints (datasets, GitHub worker, analytics).
   - If you edit the inline JSON-LD `<script type="application/ld+json">`, the CSP hash in `script-src` must be updated or the browser will block it.

4) **Service worker / SPA fallbacks**
   - Cellucid is a single page, but it is not using a framework router that needs wildcard rewrites.
   - Avoid “SPA fallback” rewrites that accidentally serve `index.html` for data files (it makes fetches look like JSON parse errors).

---

## Host datasets (exports) safely

Many Cellucid workflows involve hosting datasets separately from the app (e.g. GitHub-hosted exports, an institutional HTTP server, or a CDN).
This section highlights server behaviors that commonly break loading.

### 1) CORS: allow the app origin to fetch your dataset files

If you load a dataset hosted at `https://example.org/my_export/` from `https://www.cellucid.com`, your dataset host must send CORS headers like:
- `Access-Control-Allow-Origin: https://www.cellucid.com`

For a development or internal server, you may choose to allow `*` (public datasets only).

Exception (GitHub Pages):
- GitHub Pages does not allow custom CORS headers.
- Cellucid can still load GitHub Pages exports via the `cellucid-datasets/bridge.html` iframe bridge (no CORS headers required).

### 2) Don’t “helpfully” auto-decompress `.gz` exports

Cellucid’s data loader treats URLs ending in `.gz` as **gzip-compressed files** and transparently decompresses them in the browser (using `DecompressionStream`).

That means:
- The `.gz` response body must be the *actual gzipped bytes*.
- Your server should **not** apply `Content-Encoding: gzip` to these already-compressed files.

If a server auto-decompresses `.bin.gz` (or adds `Content-Encoding: gzip`), the browser will deliver already-decompressed bytes and Cellucid will attempt to decompress again → load failure.

### 3) Avoid returning HTML for missing JSON/binary files

Some static hosts return a friendly HTML page for 404s.
From the app’s perspective, that looks like:
- “JSON parse failed”
- “Invalid binary length”
- “Unexpected token `<`”

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
