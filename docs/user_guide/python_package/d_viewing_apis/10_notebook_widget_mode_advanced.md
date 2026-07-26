# Notebook / widget mode (advanced)

This page explains the notebook embedding details for `show(...)` and `show_anndata(...)`, including:
- why HTTPS notebooks can break “simple” localhost iframes,
- how to supply the one browser-reachable URL,
- and how to debug the browser ↔ kernel connection.

If you want the minimal notebook quickstart, go to {doc}`06_jupyter_show_and_show_anndata_quickstart`.

## At a glance

**Audience**
- JupyterHub / cloud notebook users
- anyone seeing a blank iframe, mixed-content error, or unreachable loopback URL

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

Without `client_server_url=`, Cellucid uses the exact URL of its bound server:

```text
http://127.0.0.1:<port>/?jupyter=true&viewerId=...&viewerToken=...
```

With `client_server_url=`, Cellucid uses that exact absolute HTTP(S) base.
The argument must have no credentials, query, fragment, trailing slash, or
surrounding whitespace.

Cellucid does not probe for Jupyter Server Proxy, call a Colab proxy API, or
rewrite the URL in the notebook frontend. This keeps routing explicit and makes
connectivity failures reproducible.

## Jupyter Server Proxy

For JupyterHub and many remote notebooks, Jupyter Server Proxy can expose the
Cellucid port over the notebook’s HTTPS origin.

```{note}
`jupyter-server-proxy` is already a core Cellucid dependency. A Jupyter
administrator must still enable and configure its server route; installing the
package does not make Cellucid discover or select a proxy URL.
```

Determine the proxy base that reaches the selected port and pass it when
constructing the viewer:

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(
    adata,
    port=8765,
    dataset_name="Example",
    dataset_id="example",
    client_server_url="https://notebooks.example/user/alice/proxy/8765",
)
```

For an SSH tunnel viewed from an HTTP notebook:

```python
viewer = AnnDataViewer(
    adata,
    port=8765,
    dataset_name="Example",
    dataset_id="example",
    client_server_url="http://127.0.0.1:8765",
)
```

```{warning}
If your notebook is served over HTTPS, embedding an `http://...` iframe can be blocked by the browser (mixed content).
Use an HTTPS proxy base for an HTTPS notebook.
```

## Debugging notebook connectivity (high-signal)

In a notebook, run:

```python
viewer.debug_connection()
```

This report includes:
- server health probes (`/_cellucid/health`, `/_cellucid/info`)
- the configured and browser-facing server URLs
- the verified local viewer-generation status
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
viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

2) Copy and open:

```python
print(viewer.viewer_url)
```

## Troubleshooting

- Mixed-content or unreachable iframe → expose the port through an appropriate
  proxy/tunnel and pass its exact base as `client_server_url=`.
- Web-generation startup error → diagnose the source, inventory, declared
  object, or generation directory: {doc}`15_troubleshooting_viewing`.
