# Security: CORS, origins, and mixed content

This page explains the browser security constraints that most often break notebook embedding:

- **Origins**: where the iframe thinks it is running (`scheme://host:port`)
- **CORS**: which origins are allowed to call `/_cellucid/events`
- **Mixed content**: HTTPS notebooks blocking HTTP loopback servers

If your viewer shows a blank iframe, “Failed to fetch”, or “notebook proxy required”, this is the page to read.

## Why this matters for hooks

Hooks depend on the viewer being able to send:

```text
POST /_cellucid/events
```

If the browser blocks that request (mixed content, CORS, or proxy policies), your Python callbacks will never fire.

## CORS policy (Python server)

The Cellucid server includes a small origin allowlist.

As implemented in `cellucid-python/src/cellucid/_server_base.py` (`CORSMixin._get_allowed_origin`), the server may allow:
- `http(s)://localhost:<port>` and `http(s)://127.0.0.1:<port>`
- the canonical hosted viewer origin `https://www.cellucid.com`
- `https://cellucid.com` and `https://*.cellucid.com`

Everything else is rejected (no `Access-Control-Allow-Origin`).

```{note}
In the recommended notebook mode, the viewer UI is served from the same origin as the dataset server (hosted-asset proxy),
so most requests are same-origin and do not rely on CORS at all.
```

## Mixed content (the #1 failure mode)

### The problem

If your notebook frontend is served over HTTPS (very common on JupyterHub / hosted environments),
the browser may block an iframe that tries to load or fetch from:

```text
http://127.0.0.1:<port>
```

This is “active mixed content” and browsers treat it as a security risk.

Symptoms:
- the iframe stays blank,
- the app UI loads but data fetches fail,
- browser devtools shows “Mixed Content” or “Blocked” errors,
- hooks don’t fire (`/_cellucid/events` never reaches the server).

### Cellucid’s solution: serve UI from the server origin + use notebook proxies when needed

Cellucid notebook embeds:
1. Serve the viewer UI from the same origin as the dataset server (hosted-asset proxy).
2. When a direct HTTP loopback URL is unlikely to work, embed via **Jupyter Server Proxy**:

```text
https://<notebook-origin>/<base>/proxy/<port>/?jupyter=true&viewerId=...&viewerToken=...
```

The embed script probes `/_cellucid/health` before selecting the proxy URL.

### If you see “notebook proxy required”

That message means:
- the notebook is HTTPS or remote,
- direct loopback is blocked/unreachable,
- and the proxy URL probe failed.

Fix options:
- Install/enable `jupyter-server-proxy` on the notebook server (recommended).
- Use SSH port forwarding if the kernel is remote and your browser is local.
- Advanced: set `CELLUCID_CLIENT_SERVER_URL` to a browser-reachable HTTPS URL for the server.

## Environment variables that affect connectivity

### `CELLUCID_CLIENT_SERVER_URL`

Overrides the URL the browser should use to reach the server.

Use cases:
- remote kernels + custom reverse proxies,
- environments where loopback/proxy auto-detection fails.

### `CELLUCID_WEB_PROXY_CACHE_DIR`

Controls where hosted viewer UI assets are cached.

Use cases:
- offline/air-gapped environments (cache once, reuse),
- persistent cache on shared machines,
- debugging stale cache behavior.

## How to debug in practice (recommended checklist)

1. Run in Python:
   ```python
   viewer.debug_connection()
   ```
2. In a browser tab, open:
   - `<viewer.server_url>/_cellucid/health`
3. In browser devtools → Network:
   - check requests to `/_cellucid/events` (should be 200 with JSON body)
   - check mixed-content/CORS errors

## Next steps

- Full troubleshooting guide: {doc}`14_troubleshooting_hooks`
- Security model overview: {doc}`12_security_model`

