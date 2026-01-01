# Security model

This page explains the security assumptions behind notebook embedding + hooks, and how to run Cellucid safely in local and shared environments.

If you primarily need HTTPS/mixed-content details, see {doc}`13_security_cors_origins_and_mixed_content`.

## Threat model (what Cellucid assumes)

Default assumption for notebook usage:
- You are running a notebook you trust.
- The Python kernel and the browser are under your control.
- The Cellucid server is bound to `127.0.0.1` (localhost) and is not exposed to the network.

Within that model, it’s reasonable to allow:
- the iframe to POST events back to the server (`/_cellucid/events`)
- Python to send commands into the iframe (`postMessage`)

## What is “exposed” and to whom?

### The data server

The Cellucid server exposes dataset endpoints (exported or AnnData-backed) plus:
- `/_cellucid/health`
- `/_cellucid/info`
- `/_cellucid/datasets`
- `/_cellucid/events` (POST)
- `/_cellucid/session_bundle` (POST)

In notebook mode, the server is started with `host="127.0.0.1"` by default, so it should only be reachable from the same machine.

### The viewer UI assets

In notebook mode, the server may download the viewer UI from `https://www.cellucid.com` and cache it on disk (hosted-asset proxy).

Security implications:
- Your dataset is **not uploaded** to cellucid.com by this mechanism.
- Your environment does need outbound HTTPS access to cellucid.com (at least once, unless cached).

## Authentication and integrity

### Python → iframe commands (`postMessage`)

- Commands are authenticated using a per-viewer token (`viewerToken`).
- The iframe ignores commands that don’t include the expected token.

This avoids relying on strict `postMessage` origin allowlisting, which is brittle in notebook frontends that proxy/transform origins.

### Iframe → Python events (`/_cellucid/events`)

- Events are HTTP POSTs to the local server.
- The server routes by `viewerId`.
- There is no per-event cryptographic auth token today; safety relies on:
  - localhost binding by default, and
  - CORS origin checks (see {doc}`13_security_cors_origins_and_mixed_content`).

## High-risk configurations (avoid these unless you know what you’re doing)

### Binding the server to `0.0.0.0` (LAN exposure)

If you run `cellucid serve --host 0.0.0.0` (or equivalent), anyone on the network who can reach your machine can potentially access:
- your dataset endpoints,
- hooks endpoints,
- and session bundle upload endpoints.

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
  - potentially other metadata if custom events are added.
- Session bundles can contain user-defined labels and highlight membership.

Treat session bundles as sensitive data:
- store them like you would store curated annotations,
- avoid committing private identifiers into public repos.

## Recommended safe defaults

- Keep the server bound to `127.0.0.1` unless you need remote access.
- If you need remote access, use SSH port forwarding.
- Do not run untrusted notebooks that embed arbitrary iframes and can message your browser context.
- Clear web UI cache if you suspect stale assets or want to remove cached copies:
  - `viewer.clear_web_cache()` or `cellucid.clear_web_cache()`

## Next steps

- CORS/origins/mixed-content: {doc}`13_security_cors_origins_and_mixed_content`
- Troubleshooting: {doc}`14_troubleshooting_hooks`

