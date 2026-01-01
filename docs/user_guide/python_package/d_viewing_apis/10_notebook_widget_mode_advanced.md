# Notebook / widget mode (advanced)

This page explains the notebook embedding details for `show(...)` and `show_anndata(...)`, including:
- why HTTPS notebooks can break “simple” localhost iframes,
- how the embed chooses a working URL,
- and how to debug the browser ↔ kernel connection.

If you want the minimal notebook quickstart, go to {doc}`06_jupyter_show_and_show_anndata_quickstart`.

## At a glance

**Audience**
- JupyterHub / cloud notebook users
- anyone seeing “proxy required” or mixed-content issues

## The core problem: “where is localhost?”

In notebook mode, Python runs a local server and the browser loads an iframe pointing at that server.

This is easy when:
- the notebook server is on your laptop, and
- the notebook is served over HTTP, and
- the browser can reach `http://127.0.0.1:<port>`.

It gets tricky when:
- the notebook is served over **HTTPS** (JupyterHub, cloud),
- the kernel is remote (HPC),
- or the notebook frontend rewrites/proxies URLs.

## How Cellucid chooses the iframe URL

Cellucid computes a “direct” URL (what the server thinks it is):

```text
http://127.0.0.1:<port>/?jupyter=true&viewerId=...&viewerToken=...
```

Then, inside the notebook frontend, the embed script may replace that with a safer proxy URL:

- If the notebook environment supports **jupyter-server-proxy**, it tries:

```text
<notebook-origin>/<baseUrl>/proxy/<port>/?jupyter=true&viewerId=...&viewerToken=...
```

- If you are in Colab, it uses Colab’s HTTPS port proxy.
- If you set `CELLUCID_CLIENT_SERVER_URL`, that override is used as the base.

If none of these work in an HTTPS/remote context, the iframe shows a message explaining that a proxy is required.

## Recommended solution: install/enable `jupyter-server-proxy`

For JupyterHub and many remote notebooks, `jupyter-server-proxy` is the easiest and safest fix.

Typical install (local environments):

```bash
pip install jupyter-server-proxy
```

```{note}
In managed JupyterHub deployments you might not control installs; ask your admin to enable it.
```

## Advanced solution: `CELLUCID_CLIENT_SERVER_URL`

Set this environment variable in the **kernel environment** to override what URL the browser should use to reach the server.

Examples:

- If you have a reverse proxy that exposes the Cellucid server over HTTPS:

```bash
export CELLUCID_CLIENT_SERVER_URL=https://your.domain.example/cellucid
```

- If you are using an SSH tunnel and an HTTP notebook (no mixed content), you might set:

```bash
export CELLUCID_CLIENT_SERVER_URL=http://127.0.0.1:8765
```

```{warning}
If your notebook is served over HTTPS, embedding an `http://...` iframe can be blocked by the browser (mixed content).
Prefer an HTTPS proxy URL (jupyter-server-proxy, Colab proxy, or your own HTTPS reverse proxy).
```

## Debugging notebook connectivity (high-signal)

In a notebook, run:

```python
viewer.debug_connection()
```

This report includes:
- server health probes (`/_cellucid/health`, `/_cellucid/info`)
- whether `jupyter-server-proxy` is installed
- cache status for the hosted-asset proxy UI
- recent frontend console warnings/errors (forwarded to Python when available)

Also inspect:

```python
print(viewer.server_url)
print(viewer.viewer_url)
```

## When all else fails: open the viewer manually

Even if embedding is blocked, the server can still work in a normal browser tab:

1) Start the viewer (still creates the server):

```python
viewer = show_anndata("data.h5ad")
```

2) Copy and open:

```python
print(viewer.viewer_url)
```

## Troubleshooting

- Embed shows “proxy required” → install/enable jupyter-server-proxy or set `CELLUCID_CLIENT_SERVER_URL`.
- Embed shows “viewer UI unavailable” → hosted-asset proxy cache/network issue: {doc}`15_troubleshooting_viewing`.
