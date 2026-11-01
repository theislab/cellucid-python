# Security: CORS, origins, and mixed content

## In a notebook: none of this applies

A notebook widget (`show_anndata()`, R's `cellucid_view()` in Jupyter) frames
the app from `https://www.cellucid.com` (or your `web_source_url`) and carries
every read and message over the notebook's own connection. The page never
contacts the kernel's machine, so there is no loopback address, no mixed
content, no CORS and no `Origin` check for the kernel to make. The one thing
the widget checks is that the page's messages come from its own frame and
from the app's origin ({doc}`12_security_model`).

## Browser-tab viewers

The rest of this page is about a viewer in a browser tab of its own
({class}`~cellucid.AnnDataViewer`, `serve_anndata`, `cellucid serve`), which
reaches a server on this machine:

- **Origins**: where the page thinks it is running (`scheme://host:port`)
- **CORS** and the **`Origin` check**: which pages may read the server's
  routes, and which may open the viewer channel (`/_cellucid/ws`)
- **Mixed content**: an HTTPS page fetching from an HTTP loopback server

If such a viewer shows “Failed to fetch” or a mixed-content error, this is
the page to read.

### Why this matters for hooks

Hooks and commands depend on the viewer being able to open its channel:

```text
GET /_cellucid/ws   (Upgrade: websocket)
```

If the browser or a proxy blocks it (mixed content, a proxy that does not pass
WebSockets, an `Origin` the server does not know), the viewer says it cannot
open its channel and your Python callbacks never fire.

## CORS policy (Python server)

The Cellucid server includes a small origin allowlist.

As implemented in `cellucid-python/src/cellucid/_server_base.py` (`CORSMixin._get_allowed_origin`), the server may allow:
- `http(s)://localhost:<port>` and `http(s)://127.0.0.1:<port>`
- `https://www.cellucid.com` and `https://cellucid.com` (the hosted app's origins, exactly)

Everything else is rejected (no `Access-Control-Allow-Origin`).

## The viewer channel's `Origin` check

CORS does not apply to WebSockets, so the channel checks the page's `Origin`
itself (`cellucid._viewer_channel.origin_is_allowed`). It accepts a page on:
- the same host the request's `Host` header names (the server's own page),
- a loopback name (`localhost`, `127.0.0.1`, `::1`), or
- a host declared in `allowed_hosts=` (`--allowed-host`), such as the
  JupyterHub host a proxy serves the page from.

Anything else is refused with 403 before the upgrade. The viewer's token is
checked after that, on the first frame.

```{note}
A browser-tab viewer is served from the same origin as its dataset server, so
most requests are same-origin and do not rely on CORS at all.
```

## Mixed content (the #1 failure mode)

### The problem

If a page served over HTTPS (the hosted app, a hub's proxy route) tries to
fetch from:

```text
http://127.0.0.1:<port>
```

This is “active mixed content” and browsers treat it as a security risk.

Symptoms:
- the app UI loads but data fetches fail,
- browser devtools shows “Mixed Content” or “Blocked” errors,
- hooks don’t fire (the viewer channel never opens).

### Cellucid’s solution: serve the UI from the server origin

A browser-tab viewer:
1. establishes and serves the verified viewer generation from the same origin
   as the dataset server, and
2. uses the exact `client_server_url=` supplied by the caller, or the bound
   server URL when that argument is omitted.

### If the browser cannot reach the server

Fix options:
- in a notebook, use `show_anndata()` instead (no address involved);
- use SSH port forwarding and open the forwarded local address; or
- expose the Cellucid port through an HTTPS reverse proxy, pass that exact
  HTTP(S) base as `client_server_url=`, and its host as `allowed_hosts=`.

Cellucid does not probe, infer, or select a proxy URL.

## Explicit connectivity arguments

### `client_server_url`

Sets the exact URL the browser uses to reach the server. It is for remote
kernels, explicit reverse proxies, and notebook frontends that cannot reach
the bound server URL. The viewer channel uses the same address and path,
`ws://` for `http://` and `wss://` for `https://`.

### `allowed_hosts`

Extra host names the server answers to (`Host` header) and accepts pages
from (`Origin` header). Name the proxy's host here when `client_server_url`
points at a reverse proxy (a hub's proxy route).

### `web_cache_dir`

Sets the directory where the verified viewer generation is published.

Use cases:
- selecting a controlled writable cache location,
- verifying an existing generation with `force=False`,
- clearing and re-establishing the current source generation.

## How to debug in practice (recommended checklist)

1. Run in Python (a browser-tab viewer):
   ```python
   viewer.debug_connection()
   ```
2. In a browser tab, open:
   - `<viewer.server_url>/_cellucid/health`
3. In browser devtools → Network:
   - filter `WS`: the `/_cellucid/ws` connection should be `101 Switching
     Protocols`, with `hello` and `welcome` as its first frames
   - check mixed-content/CORS errors, and a 403 (page origin not accepted) or
     421 (host not accepted) on the upgrade

## Next steps

- Full troubleshooting guide: {doc}`14_troubleshooting_hooks`
- Security model overview: {doc}`12_security_model`
