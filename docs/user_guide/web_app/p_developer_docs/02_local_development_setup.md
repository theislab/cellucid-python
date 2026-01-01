# Local development setup

This page gives a **repeatable, low-friction** way to run the Cellucid web app locally for debugging and development.

Because the app is shipped as **native ES modules** (no bundler required), the setup is usually just “run a static server and reload”.

## At a glance

**Audience**
- Wet lab / non-technical: you can use this to reproduce bugs and capture screenshots, without writing code.
- Computational users: use this to run on the same machine as `cellucid-python` server mode.
- Developers: use this for rapid iteration; then read {doc}`13_debugging_playbook`.

**Time**
- 2–5 minutes for “hello world” local run
- 10–20 minutes if you also want to test remote server + Jupyter workflows

**Prerequisites**
- A modern desktop browser with WebGL2 (Chrome/Edge/Firefox recommended)
- Python 3 (for a simple local static server), or any equivalent local server

---

## Step 0: Know the two “local” modes

There are two distinct local workflows:

1) **Local web app**: you run the frontend from this workspace (`cellucid/`), and you load demo or local files.
2) **Local data server**: you run `cellucid-python` to serve a dataset, and the frontend connects to it (either from your local web app or from `https://www.cellucid.com`).

They are often used together, but debugging is faster when you can clearly say which one you’re using.

---

## Step 1: Run the Cellucid web app locally

From the workspace root:

```bash
cd cellucid
python -m http.server 8000
```

Then open:
- `http://localhost:8000`

:::{important}
Do **not** open `cellucid/index.html` via `file://`.
Browsers will block ES module imports and cross-origin fetches in ways that look like “random failures”.
Always use an HTTP server.
:::

---

## Step 2: Turn on debug logging (recommended)

In DevTools → Console:

```js
localStorage.setItem('CELLUCID_DEBUG', 'true');
location.reload();
```

To disable:

```js
localStorage.removeItem('CELLUCID_DEBUG');
location.reload();
```

What this does:
- Enables the lightweight logger in `cellucid/assets/js/utils/debug.js`.
- You’ll see prefixed logs like `[Cellucid] ...` across data/state/UI flows.

---

## Step 3: Load a dataset (choose one)

### Option A: Use built-in demo datasets (fastest)

If the app can fetch `datasets.json` from the configured exports base URL, the dataset dropdown will populate from the **local-demo** source.

Configure the exports base URL via:
- `<meta name="cellucid-exports-base-url" ...>` in `cellucid/index.html`, or
- `?exportsBaseUrl=...` as a runtime override (useful in development).

This is best for:
- UI work
- rendering/LOD work
- session/export work

### Option B: Use the browser file picker (local-user)

Use the “Browse local data…” UI to select:
- a prepared `exports/` folder, or
- a `.h5ad` file, or
- a `.zarr` directory

This is best for:
- “real user” workflows with no Python server running
- testing permission/UX issues around file pickers

### Option C: Connect to a local Python server (recommended for large h5ad/zarr)

Run (in another terminal) from your `cellucid-python` environment:

```bash
cellucid serve /path/to/dataset.h5ad
# or
cellucid serve /path/to/dataset.zarr
```

Then in the Cellucid UI, use the “Remote server” connection and point it at:
- `http://127.0.0.1:<port>`

This is best for:
- large datasets (lazy loading)
- reproducible debugging (server logs + client logs)

:::{note}
If you are using the hosted app (`https://www.cellucid.com`) and a local server, you are doing an **https → http** connection.
Most modern browsers allow `https://…` pages to talk to `http://localhost` / `http://127.0.0.1`, but if your environment blocks it, run the web app locally (`http://localhost:8000`) and connect to the server from there.
:::

### Option D: Jupyter embedded viewer (hooks + events)

If you are developing the Python↔frontend hooks system:
- Read the Python docs section: {doc}`../../python_package/e_jupyter_hooks/index`
- And the frontend development guide: `cellucid/markdown/HOOKS_DEVELOPMENT.md`

---

## Step 4: Make local iteration fast (recommended defaults)

### Disable cache while debugging

In Chrome/Edge DevTools:
- Network tab → check “Disable cache”
- Keep DevTools open while reloading

This avoids “I fixed it but it still reproduces” confusion caused by cached ES modules.

### Keep reproducible datasets handy

For bugs that involve state/persistence, keep a small “debug dataset” around:
- ~5k–50k cells (fast reloads)
- at least one categorical field with many categories
- at least one continuous field with outliers and NaNs
- optionally one vector field overlay

This helps you reproduce edge cases quickly without waiting for large loads.

---

## Troubleshooting

### Symptom: “Blank page / canvas is black”

Likely causes (ordered):
1) WebGL2 not available (browser/GPU issue).
2) A JS exception during startup prevented viewer initialization.
3) The app cannot fetch required files due to server/CORS/missing paths.

How to confirm:
- DevTools → Console: look for the first error after “Starting…”.
- DevTools → Network: reload and look for `404` or blocked requests.

Fix:
- Try a different browser (Chrome/Edge first).
- Ensure you are serving from the `cellucid/` directory.
- Enable debug logs and reload (Step 2).

### Symptom: “Module script failed to load”

Likely causes:
- You opened via `file://`, or the server is sending incorrect MIME types for `.js`.

Fix:
- Use `python -m http.server` from `cellucid/` (Step 1).

### Symptom: “Remote server connect fails / CORS error”

Likely causes:
- The server is not running.
- The server does not send CORS headers for the Cellucid origin you’re using.
- Corporate proxy/ad blocker interferes with localhost requests.

How to confirm:
- Open the remote server URL in a new tab and see if it responds.
- Check DevTools → Network for failed `fetch` requests and CORS messages.

Fix:
- Use the locally hosted app (`http://localhost:8000`) and connect from there.
- Or configure the server origin/CORS (Python side).

---

Next: {doc}`03_build_run_and_deployment`.
