# Configuration, environment variables, and logging

This page lists the **configuration knobs** available in `cellucid-python`: CLI flags, function arguments, environment variables, and logging behavior.

---

## Configuration layers (highest priority wins)

In general, configuration flows like this:

1) **Function arguments** (e.g. `serve(..., host=..., port=...)`)
2) **CLI flags** (e.g. `cellucid serve ... --host ... --port ...`)
3) **Environment variables** (only for a small set of notebook/web-cache behaviors)
4) **Defaults** in code (`DEFAULT_HOST`, `DEFAULT_PORT`)

---

## Defaults

Defined in `cellucid-python/src/cellucid/_server_base.py`:

- `DEFAULT_HOST = "127.0.0.1"` (local only by default)
- `DEFAULT_PORT = 8765`

For remote access, explicitly bind:
- `--host 0.0.0.0` (and use SSH tunneling or a firewall)

---

## Environment variables

### `CELLUCID_WEB_PROXY_CACHE_DIR`

Controls where the **hosted-asset proxy** stores cached web UI files (downloaded from `https://www.cellucid.com`).

Why it matters:
- offline/airgapped environments,
- persistent caching across notebook restarts,
- avoiding re-downloading large UI bundles repeatedly.

Example:

```bash
export CELLUCID_WEB_PROXY_CACHE_DIR="$HOME/.cache/cellucid-web"
```

Related APIs:
- `cellucid.get_web_cache_dir()`
- `cellucid.clear_web_cache()`

### `CELLUCID_CLIENT_SERVER_URL`

Overrides the URL that the **browser** should use to reach the Python server.

This is primarily for remote notebook environments where:
- the kernel runs on a remote VM,
- and `http://127.0.0.1:<port>` is not reachable from your browser.

If you have an HTTPS reverse proxy for the server (or a stable tunnel URL), set:

```bash
export CELLUCID_CLIENT_SERVER_URL="https://<your-proxy-host>/proxy/<port>"
```

The notebook embedder will use this as the base for `viewer.viewer_url`.

See details in: {doc}`10_jupyter_embedding_architecture`

---

## Logging

### CLI logging

The CLI uses Python logging with two user-facing knobs:

- `--verbose`: `DEBUG` logs (more detail, stack traces)
- `--quiet`: suppresses most informational prints

### Library logging

When you import `cellucid` as a library, logging is “normal Python logging”:
- by default you may see little/no output,
- configure logging in your application/notebook if you want debug output.

Example:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

---

## CORS / origin policy (important for security)

Servers add CORS headers that allow:
- `localhost` / `127.0.0.1`
- the canonical hosted viewer origin (`https://www.cellucid.com`)
- `https://cellucid.com` and `https://*.cellucid.com`

This is implemented in:
- `cellucid-python/src/cellucid/_server_base.py` (`CORSMixin._get_allowed_origin`)

If you run a server on a public host, do not assume CORS is a full security boundary.
Treat the server as trusted/local unless you have a hard threat model and additional protections (see {doc}`17_security_privacy_and_networking`).

---

## Troubleshooting

### Symptom: “Notebook says mixed content is blocked”

Cause:
- the notebook is served over HTTPS,
- but the viewer tries to load `http://127.0.0.1:<port>`.

Fix options (ordered):
1) Use `jupyter-server-proxy` (recommended).
2) Use a tunnel/reverse proxy and set `CELLUCID_CLIENT_SERVER_URL`.

See: {doc}`10_jupyter_embedding_architecture`

### Symptom: “Viewer UI can’t be loaded (offline)”

Cause:
- hosted UI couldn’t be fetched,
- and no cached UI exists.

Fix:
- run once while online to populate the cache,
- ensure `CELLUCID_WEB_PROXY_CACHE_DIR` is writable and persistent,
- or clear a corrupted cache via `cellucid.clear_web_cache()`.
