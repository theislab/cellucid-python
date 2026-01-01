# Privacy, security, and offline vs online

**Audience:** everyone (especially anyone handling sensitive data)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- When your data stays local vs when it can be exposed
- What the Python server is doing with `cellucid.com` (hosted UI assets)
- How to make Cellucid predictable in offline / firewalled environments

---

## Mental model (one sentence)

Cellucid is a web app UI, but with `cellucid-python` your **dataset is served by your Python process** — the main privacy/security decisions are about *where you host the server* and *whether the viewer UI assets are cached for offline use*.

---

## “Does my dataset get uploaded to cellucid.com?”

### Notebook/server modes (most `cellucid-python` usage)

By default:
- your dataset is served by a local Python server (usually `127.0.0.1:<port>`),
- the viewer UI is served from the same origin via a hosted-asset proxy + local cache,
- the browser loads data from your local server.

Your dataset is **not uploaded** anywhere automatically.

### Web app file picker mode (no Python server)

If you open the Cellucid web app and use the browser file picker to select an export folder:
- the browser reads files locally,
- nothing is uploaded unless you explicitly host/share the export folder yourself.

Start here: {doc}`../../web_app/b_data_loading/index`

---

## Security model of the Python server (what can be exposed)

### Default: localhost-only

`CellucidServer` and notebook viewers bind to `127.0.0.1` by default.

That means:
- only the same machine can reach the server,
- people on your network cannot access it unless you deliberately bind to a public interface or tunnel it.

### If you bind to a non-localhost host (be careful)

If you run a server on:
- `0.0.0.0` (all interfaces), or
- a public IP / hostname,

then **anyone who can reach that host/port can fetch your dataset files**.

Practical implications:
- this is not an authenticated server by default,
- treat it like “public to whoever can reach it”.

Safe sharing patterns:
- SSH port forwarding (keeps the server bound to localhost on the remote machine)
- VPN + firewall rules
- hosting behind an authenticated reverse proxy (advanced)

---

## CORS and why it exists

The Python server sets CORS headers so the viewer can load data in different deployment modes.

Two common scenarios:
- Viewer UI served from `https://www.cellucid.com` but data served from your server (cross-origin).
- Viewer UI served from your local server (same-origin; still safe to keep consistent headers).

If you see browser errors mentioning CORS, it usually means:
- the browser could not reach the server at all (network/proxy),
- or you are mixing HTTPS and HTTP in a way the browser blocks (mixed content).

Start debugging at: {doc}`08_debugging_mental_model_where_to_look`.

---

## Offline vs online: hosted-asset proxy and the UI cache

### Why the viewer UI is “downloaded once”

In notebook/server modes, Cellucid serves the web app UI from your Python server so:
- the viewer and dataset share the same origin (avoids mixed-content),
- your notebook embed is more reliable across environments.

If the UI is not packaged locally, the server runs in **hosted-asset proxy mode**:
- it fetches `index.html` and `/assets/...` from `https://www.cellucid.com`,
- caches them on disk,
- and then serves them locally.

### What this means for offline environments

- First run in a new environment may require internet access to populate the cache.
- After the cache is populated, you can often run fully offline (the server serves cached UI assets).

### Where the cache lives

The cache directory is:
- `CELLUCID_WEB_PROXY_CACHE_DIR` if set, else
- a temp directory (platform-dependent).

In Python:

```python
from cellucid import get_web_cache_dir
get_web_cache_dir()
```

### Prefetch the UI intentionally (recommended for workshops / offline HPC)

Run once while online:

```python
from cellucid import show_anndata

viewer = show_anndata("data.h5ad", auto_open=False)
viewer.ensure_web_ui_cached(force=False, show_progress=True)
viewer.stop()
```

Or clear and re-download (useful when debugging stale caches):

```python
from cellucid import clear_web_cache
clear_web_cache()
```

```{important}
If you rely on offline use, set `CELLUCID_WEB_PROXY_CACHE_DIR` to a persistent, writable location.
Otherwise a system temp directory may be cleared between sessions.
```

---

## Privacy implications of session bundles

A `.cellucid-session` file is not “just UI state” — it can encode:
- highlight membership lists (which cells you selected),
- user-defined labels and categories,
- potentially other feature-specific state as the format evolves.

Treat session bundles as data artifacts:
- don’t share them publicly if the underlying dataset is sensitive,
- version them alongside the matching dataset export folder.

---

## Screenshot placeholder: offline UI-unavailable error page

If you want to help beginners, capture the UI error page shown when the UI cannot be fetched and no cache exists.

<!-- SCREENSHOT PLACEHOLDER
ID: offline-web-ui-unavailable
Suggested filename: web_app/offline-web-ui-unavailable.png
Where it appears: Python Package → Concepts → Privacy/offline → Offline vs online
Capture:
  - UI location: browser tab or notebook iframe
  - State prerequisites: run Cellucid in a fresh environment with no cached UI while offline
  - Action to reach state: start a viewer; it fails to load the UI and shows an error page
Crop:
  - Include: the main error message + suggested fixes (run once online, set cache dir)
  - Exclude: private dataset names and local paths if shown
Alt text:
  - Error page indicating the Cellucid viewer UI could not be loaded and suggesting caching steps.
Caption:
  - The Python server proxies and caches the web UI assets; when offline with no cache, you’ll see this page—run once online or set a persistent cache directory.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the offline viewer UI unavailable page.
:width: 100%

When offline and the viewer UI is not cached, the Python server cannot serve the web app UI until you run once online or configure a persistent cache.
```

---

## Troubleshooting (privacy/offline/security)

### Symptom: “I’m offline and the viewer UI won’t load”

Likely causes:
- no cached copy of the viewer UI exists yet,
- cache directory is not writable or not persistent.

Fix:
- run once while online to populate the cache, or
- set `CELLUCID_WEB_PROXY_CACHE_DIR` to a persistent location and prefetch the UI.

### Symptom: “Someone else on my network can see my dataset”

Likely causes:
- you bound the server to `0.0.0.0` or a public interface,
- your machine firewall allows inbound traffic to the port.

Fix:
- bind to `127.0.0.1`,
- use SSH tunneling for remote access,
- add firewall rules.

---

## Next steps

- Performance/scaling: {doc}`07_performance_mental_model_and_scaling`
- Debugging checklist: {doc}`08_debugging_mental_model_where_to_look`
