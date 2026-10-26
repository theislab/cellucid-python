# Notebook / widget mode (advanced)

This page explains the notebook embedding details for `show_anndata(...)`, including:
- why HTTPS notebooks can break “simple” localhost iframes,
- how to supply the one browser-reachable URL,
- what the constructor is doing during a long, silent cell,
- how to choose the neighbor graph, which a notebook viewer reads only when a
  view draws it,
- and how to debug the browser ↔ kernel connection.

If you want the minimal notebook quickstart, go to {doc}`06_jupyter_show_anndata_quickstart`.

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
http://127.0.0.1:<port>/?viewer=true&viewerId=...&viewerToken=...&anndata=true
```

`<port>` is 8765 unless another viewer or program holds it, in which case a
viewer created without `port=` takes the next free port and prints a note; a
`port=` you pass is used exactly or refused.

With `client_server_url=`, Cellucid uses that exact absolute HTTP(S) base.
The argument must have no credentials, query, fragment, trailing slash, or
surrounding whitespace.

Cellucid does not probe for Jupyter Server Proxy, call a Colab proxy API, or
rewrite the URL in the notebook frontend. This keeps routing explicit and makes
connectivity failures reproducible.

## A notebook viewer starts quietly

`AnnDataViewer` and `show_anndata(...)` run their server with output suppressed,
so the numbered startup report a terminal prints — `[1/4] Detecting format`
through `[4/4] Starting server` — never appears in a notebook. None of those
steps reads an array: the server reads the store's metadata, checks the launch
choices you passed, and binds; the viewer then reads what each view needs. Two
things can make the constructor take longer:

| Phase | When it runs |
|---|---|
| Writing a temporary `.h5ad` (`.h5mu`) | an in-memory `AnnData` (`MuData`), once; or an AnnData store with an array the browser cannot decode |
| Viewer generation and bind | every start; slowest on the first run on a machine, or after a cleared cache |

To see the report, run the same object through a terminal server once —
`serve_anndata(..., quiet=False)` prints the four steps.

## The neighbor graph

A notebook viewer reads the neighbor graph only when a view draws it. Without
`knn_key` it uses `obsp['connectivities']` when the store has it; `knn_key`
names another `obsp` matrix. The parameter is keyword-only and accepted by both
{class}`~cellucid.AnnDataViewer` and {func}`~cellucid.show_anndata`:

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(
    adata,
    port=8765,
    dataset_name="Example",
    dataset_id="example",
    knn_key="connectivities",
)
```

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    knn_key="connectivities",
)
```

```{important}
Naming a matrix the object's `obsp` does not hold raises during construction
rather than producing a graph-less viewer. Compute it first with
`sc.pp.neighbors(adata)`, or leave the parameter unset.
```

The rules are in {doc}`08_anndata_mode_show_anndata_and_serve_anndata`.

## Jupyter Server Proxy

For JupyterHub and many remote notebooks, Jupyter Server Proxy can expose the
Cellucid port over the notebook’s HTTPS origin.

```{note}
`jupyter-server-proxy` is already a core Cellucid dependency. A Jupyter
administrator must still enable and configure its server route; installing the
package does not make Cellucid discover or select a proxy URL.
```

Determine the proxy base that reaches the selected port and pass it — together
with the proxy’s host name — when constructing the viewer:

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(
    adata,
    port=8765,
    dataset_name="Example",
    dataset_id="example",
    client_server_url="https://notebooks.example/user/alice/proxy/8765",
    allowed_hosts=["notebooks.example"],
)
```

`allowed_hosts=` is required here, not optional. Jupyter Server Proxy forwards
the browser’s `Host` header verbatim, so the loopback-bound Cellucid server
receives `notebooks.example` rather than `127.0.0.1` and refuses it with
`421 Misdirected Request` — the same refusal that stops DNS rebinding. A proxied
request is indistinguishable from a rebound one, so the proxy can only be
declared, never detected. Each entry is one bare host name with no port, scheme
or wildcard; see
{doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`.

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
- `anndata_store_probes`, keyed by every exact listed dataset id
- the configured and browser-facing server URLs
- the verified local viewer-generation status
- `channel_connected`: whether the viewer's WebSocket channel is open
- `frontend_roundtrip` and `frontend_debug_snapshot`: a ping over the channel
  and the page's own view of its URL and origin (each says why it failed when
  the channel is not connected)
- recent accepted-event counts grouped by exact event type, and the commands
  the viewer refused (`recent_events["command_errors"]`)

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
- The cell has produced no output for minutes → you are inside the constructor.
  Nothing is printed there; see “A notebook viewer starts quietly” above.
- The edge overlay is missing although `adata.obsp['connectivities']` exists →
  the graph is read only when a view draws edges; turn the edges on in the
  viewer, or choose the graph under **Session › Data mapping**.
- Construction raised “knn_key names obsp['connectivities'],
  which this store does not have.” → run `sc.pp.neighbors(adata)` first, or
  construct without the parameter.
- The viewer loads but says it “cannot open its channel”, or hooks and commands
  never arrive → a proxy between the page and the server does not pass
  WebSockets (`/_cellucid/ws`); `viewer.debug_connection()` shows
  `channel_connected: False`.
- Construction raised “Cellucid found nothing to plot” → no `obsm` array can be
  plotted. The message lists every `obsm` entry and why it was not used; name a
  usable one with `embedding_keys={"2d": "<key>"}`.
  {doc}`08_anndata_mode_show_anndata_and_serve_anndata` has the complete rules.
