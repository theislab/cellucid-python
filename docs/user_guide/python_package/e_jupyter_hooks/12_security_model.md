# Security model

This page explains the security assumptions behind notebook embedding + hooks, and how to run Cellucid safely in local and shared environments.

If you primarily need HTTPS/mixed-content details, see {doc}`13_security_cors_origins_and_mixed_content`.

## Threat model (what Cellucid assumes)

Default assumption for notebook usage:
- You are running a notebook you trust.
- The Python kernel and the browser are under your control.
- The Cellucid server is bound to `127.0.0.1` (localhost) and is not exposed to the network.

Within that model, it’s reasonable to allow the viewer iframe to open one
authenticated WebSocket to the server (`/_cellucid/ws`), on which Python sends
commands and the viewer sends events back.

## What is “exposed” and to whom?

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

In notebook mode, the server is always bound to `127.0.0.1` (`AnnDataViewer`
has no `host` argument), so it is reachable only from the same machine; a
browser elsewhere reaches it through a tunnel or proxy you set up.

### The viewer UI assets

In notebook mode, each viewer startup downloads the declared web generation
from its configured source, verifies it completely, and publishes it locally
before the server binds.

Security implications:
- Your dataset is **not uploaded** to cellucid.com by this mechanism.
- Every viewer-serving startup needs outbound access to its configured
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
- The notebook output holds no script, only the iframe, so nothing runs in
  the notebook page itself.

## High-risk configurations (avoid these unless you know what you’re doing)

### Binding the server to `0.0.0.0` (LAN exposure)

If you run `cellucid serve --host 0.0.0.0` (or equivalent), anyone on the network who can reach your machine can potentially access:
- your dataset endpoints and the dataset list;
- the viewer channel is still closed to anyone without a viewer's token, and to
  pages on other sites by the `Origin` check.

If you need remote access, prefer:
- SSH tunneling, or
- a properly secured reverse proxy with authentication.

### Running on shared JupyterHub without proxy hardening

On shared infrastructure, ensure:
- notebook server access is authenticated,
- ports are not exposed directly,
- jupyter-server-proxy is configured appropriately (if used).

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

- Keep the server bound to `127.0.0.1` unless you need remote access.
- If you need remote access, use SSH port forwarding.
- Do not run untrusted notebooks: their code runs in your kernel, where it can
  read any viewer's token.
- Clear web UI cache if you suspect stale assets or want to remove cached copies:
  - `viewer.clear_web_cache()` or `cellucid.clear_web_cache()`

## Next steps

- CORS/origins/mixed-content: {doc}`13_security_cors_origins_and_mixed_content`
- Troubleshooting: {doc}`14_troubleshooting_hooks`
