# Security: CORS, origins, and mixed content

This page explains the browser security constraints that most often break notebook embedding:

- **Origins**: where the iframe thinks it is running (`scheme://host:port`)
- **CORS**: which origins are allowed to call `/_cellucid/events`
- **Mixed content**: HTTPS notebooks blocking HTTP loopback servers

If your viewer shows a blank iframe, “Failed to fetch”, or a mixed-content
error, this is the page to read.

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
In notebook mode, the exact verified viewer generation is served from the same
origin as the dataset server,
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

### Cellucid’s solution: serve the UI from the server origin

Cellucid notebook embeds:
1. establish and serve the verified viewer generation from the same origin as
   the dataset server, and
2. use the exact `client_server_url=` supplied by the caller, or the bound
   server URL when that argument is omitted.

### If direct loopback is blocked

This means:
- the notebook is HTTPS or remote,
- and the browser cannot load the kernel machine’s HTTP loopback URL.

Fix options:
- expose the Cellucid port through Jupyter Server Proxy or another HTTPS
  reverse proxy, then pass that exact HTTP(S) base as `client_server_url=`, or
- use SSH port forwarding and pass the local forwarded base.

Cellucid does not probe, infer, or select a proxy URL.

## Explicit connectivity arguments

### `client_server_url`

Sets the exact URL the browser uses to reach the server. It is for remote
kernels, explicit reverse proxies, and notebook frontends that cannot reach
the bound server URL.

### `web_cache_dir`

Sets the directory where the verified viewer generation is published.

Use cases:
- selecting a controlled writable cache location,
- verifying an existing generation with `force=False`,
- clearing and re-establishing the current source generation.

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
