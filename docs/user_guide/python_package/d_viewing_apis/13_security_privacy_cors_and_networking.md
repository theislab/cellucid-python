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

### The verified viewer UI generation

The Python server downloads the **viewer UI assets** from:

```text
https://www.cellucid.com
```

and publishes them as one verified local generation. This lets the UI and
dataset API use the same origin.

Implications:
- every viewer-serving startup requires access to the configured source,
- the complete declared generation is verified before the server starts,
- source failure is reported instead of substituting a previous generation,
- the local files are UI assets, not your dataset.

Configure the local generation directory with `--web-cache-dir PATH` or the
`web_cache_dir=` Python argument.

## The biggest risk: exposing the dataset server

### Default behavior (safe)

By default, the server binds to:

```text
127.0.0.1
```

Meaning: only your machine can access it.

A loopback bind is still not the whole boundary. Any page you have open can
re-point its own DNS name at `127.0.0.1` (DNS rebinding), after which the
browser treats this server as same-origin and applies no CORS check. Cellucid
closes that by validating `Host` on every request before routing: only one
well-formed authority naming `localhost` or a loopback IP literal on the bound
port is served, and everything else gets `421 Misdirected Request`.

Behind a reverse proxy (`jupyter-server-proxy`, a Colab proxy port), the
forwarded `Host` is the proxy's, so it must be declared:
`--allowed-host hub.example.org` or `allowed_hosts=["hub.example.org"]`. One
bare host name per entry, no port, no scheme, no wildcard. Details:
{doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`.

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
- Viewer source/generation issues: {doc}`15_troubleshooting_viewing`
