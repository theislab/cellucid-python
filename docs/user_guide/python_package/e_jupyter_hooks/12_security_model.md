# Security model

This page explains the security assumptions behind notebook embedding + hooks, and how to run Cellucid safely in local and shared environments.

If you primarily need origin and HTTPS details, see {doc}`13_security_cors_origins_and_mixed_content`.

## Threat model (what Cellucid assumes)

Default assumption for notebook usage:
- You are running a notebook you trust.
- The Python kernel and the browser are under your control.

Within that model, the viewer opens one authenticated channel to Python, on
which Python sends commands and the viewer sends events back.

## The notebook widget (`show_anndata`)

A notebook widget opens **no port and no server**. Three hops carry
everything, each authenticated by what carries it:

- **page ⇄ widget**: the widget's front end accepts the page's `ready` only
  from its own frame and from the app's origin (`web_source_url`), and
  answers with a private `MessagePort` and the viewer's id and token, posted
  to that exact origin. The token never appears in an address, a widget trait
  (which a notebook may save) or a log.
- **widget ⇄ kernel**: the notebook's own connection (a Jupyter comm), as for
  every widget; whoever can use your notebook server can use your kernel.
- **in the kernel**: the same viewer channel rules as below (`hello` with the
  token, closed field sets, a newer channel replacing an older one).

What leaves the kernel: the dataset's bytes the viewer reads, only to the
page in the widget's output. The app itself comes from `web_source_url`
(`https://www.cellucid.com` by default) and runs in your browser; your data is
not uploaded anywhere. A copy of the app you host can be named with
`web_source_url=`.

## A browser-tab viewer (`AnnDataViewer`, `serve_anndata`, `cellucid serve`)

Here the viewer reaches Python over HTTP and a WebSocket to a server on this
machine, bound to `127.0.0.1` unless you choose otherwise.

### The data server

The Cellucid server exposes the dataset (the bytes of one file or folder under
`/_cellucid/anndata/<dataset_id>…`) plus:
- `/_cellucid/health`
- `/_cellucid/info`
- `/_cellucid/datasets`
- `/_cellucid/protocol` (where the viewer channel is and the commands and
  events this `cellucid` version accepts — the same answer for every viewer and
  every dataset, so it exposes nothing about your data)
- `/_cellucid/ws` (the viewer channel, a WebSocket; nothing passes on it before
  a viewer's token is checked)

It answers `GET` and `HEAD` (and the CORS preflight `OPTIONS`); a `POST` is
refused (405).

`AnnDataViewer`'s server is always bound to `127.0.0.1` (it has no `host`
argument), so it is reachable only from the same machine; a browser elsewhere
reaches it through a tunnel or proxy you set up.

### The viewer UI assets

Each browser-tab server startup downloads the declared web generation from
its configured source, verifies it completely, and publishes it locally
before the server binds.

Security implications:
- Your dataset is **not uploaded** to cellucid.com by this mechanism.
- Every browser-tab server startup needs outbound access to its configured
  `web_source_url`; an existing generation directory is not an offline
  substitute.

## Authentication and integrity

### The viewer channel (`/_cellucid/ws`)

- The upgrade passes the server's `Host` check (421 otherwise) and an `Origin`
  check: the page must be this server's own host, a loopback name, or a host
  in `allowed_hosts` (403 otherwise). WebSockets are not covered by CORS, so
  this check is what keeps a page on another site from opening a channel.
- The first frame must be `hello` with the viewer's `viewerId` and
  `viewerToken`; the server compares the token in constant time and closes an
  unknown viewer or a wrong token with 4401 before anything else passes.
- After `welcome` the socket is that one viewer's: commands and events carry
  no credentials, and a newer channel of the same viewer replaces the older
  one (4409).
- Every command and event is checked against its exact closed field set; a
  frame that breaks the rules closes the channel (1008) before any hook runs.

## High-risk configurations (avoid these unless you know what you’re doing)

### Binding the server to `0.0.0.0` (LAN exposure)

If you run `cellucid serve --host 0.0.0.0` (or equivalent), anyone on the network who can reach your machine can potentially access:
- your dataset endpoints and the dataset list;
- the viewer channel is still closed to anyone without a viewer's token, and to
  pages on other sites by the `Origin` check.

If you need remote access, prefer:
- SSH tunneling, or
- a properly secured reverse proxy with authentication.

### Shared JupyterHub

On shared infrastructure, ensure notebook server access is authenticated. A
notebook widget needs no port or proxy route; if you serve a browser-tab
viewer through a hub's proxy route, declare the hub's host name with
`allowed_hosts=` and keep the route behind the hub's authentication.

## Data privacy notes

- Hooks events can contain:
  - selected cell indices,
  - hover/click indices,
  - the exact ready and diagnostic fields documented in
    {doc}`07_frontend_to_python_events`.
- Session bundles can contain user-defined labels and highlight membership.

Treat session bundles as sensitive data:
- store them like you would store curated annotations,
- avoid committing private identifiers into public repos.

## Recommended safe defaults

- In a notebook, use the widget (`show_anndata`): nothing to bind or expose.
- Keep a browser-tab server bound to `127.0.0.1` unless you need remote
  access; for remote access, use SSH port forwarding.
- Do not run untrusted notebooks: their code runs in your kernel, where it can
  read any viewer's token.
- Clear a browser-tab server's web UI cache if you suspect stale assets or want to remove cached copies:
  - `viewer.clear_web_cache()` or `cellucid.clear_web_cache()`

## Next steps

- CORS/origins/mixed-content: {doc}`13_security_cors_origins_and_mixed_content`
- Troubleshooting: {doc}`14_troubleshooting_hooks`
