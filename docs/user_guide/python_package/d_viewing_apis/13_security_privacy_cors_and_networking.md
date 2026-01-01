# Security, privacy, CORS, and networking

This page documents the security/privacy model of Cellucid viewing workflows, especially:
- what is (and is not) sent over the network,
- why `--host 0.0.0.0` is risky,
- and how CORS behaves.

## At a glance

**Audience**
- Anyone working with sensitive data (patient, unpublished, internal)
- Anyone serving on shared networks

## What network traffic happens in the recommended workflow?

### Your dataset

In the recommended workflows (`cellucid serve`, `show`, `show_anndata`):
- your **dataset bytes are served by your local Python server**,
- and your **browser downloads the dataset from that local server**.

Cellucid does not “upload your dataset to cellucid.com” as part of this workflow.

### The viewer UI (hosted-asset proxy)

The Python server downloads the **viewer UI assets** from:

```text
https://www.cellucid.com
```

and caches them locally. This is done so the UI can be served from the same origin as your dataset server.

Implications:
- the first run may require internet access,
- offline use is supported *after* assets are cached,
- the cached assets are UI files (JS/CSS), not your data.

Configure cache location with `CELLUCID_WEB_PROXY_CACHE_DIR`.

## The biggest risk: exposing the dataset server

### Default behavior (safe)

By default, the server binds to:

```text
127.0.0.1
```

Meaning: only your machine can access it.

### Risky behavior: `--host 0.0.0.0`

Binding to `0.0.0.0` makes the server reachable from other machines on the network (subject to firewall rules).

This can expose:
- embeddings and cell metadata (obs fields),
- gene expression values,
- and any other exported artifacts the server can serve.

```{warning}
CORS does not “secure” a server. CORS is a browser policy; any machine that can reach your server can still download data using non-browser clients (curl, Python, etc.).
Treat the server as a dataset-serving HTTP endpoint.
```

Recommended alternative for remote access: {doc}`12_remote_servers_ssh_tunneling_and_cloud` (SSH tunnel).

## CORS behavior (what it is and why it exists)

Cellucid’s viewer UI may run on a different origin than your dataset server in some workflows.
To make browser requests work, the server sends CORS headers.

High-level behavior:
- The server only allows a limited set of origins (loopback + Cellucid’s hosted origin).
- Requests from random websites should not be able to read your dataset from a browser tab.

## Notebook-specific security notes

Notebook mode uses two channels:

1) Python → frontend commands via `postMessage` (includes a per-viewer `viewerToken`)
2) Frontend → Python events via `POST /_cellucid/events` (identified by `viewerId`)

Practical guidance:
- keep the server bound to localhost unless you have a strong reason,
- don’t embed the viewer inside untrusted web pages,
- and treat session bundles as potentially sensitive (they can encode selections, highlights, and analysis state).

## Troubleshooting security/networking issues

- Remote access (SSH tunnels): {doc}`12_remote_servers_ssh_tunneling_and_cloud`
- Viewer UI cache/offline issues: {doc}`15_troubleshooting_viewing`
