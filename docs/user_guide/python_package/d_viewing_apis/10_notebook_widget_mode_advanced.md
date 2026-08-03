# Notebook / widget mode (advanced)

This page explains the notebook embedding details for `show(...)` and `show_anndata(...)`, including:
- why HTTPS notebooks can break “simple” localhost iframes,
- how to supply the one browser-reachable URL,
- what the constructor is doing during a long, silent cell,
- how to ask for the neighbor graph, which a notebook viewer does not read by
  default,
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

## A notebook viewer builds quietly

`AnnDataViewer` and `show_anndata(...)` run their server with output suppressed,
so the numbered startup report a terminal prints — `[1/5] Detecting format`
through `[5/5] Starting server` — never appears in a notebook. On a small object
that is invisible. On a very large one it means the cell runs for minutes with
no output at all, and the constructor is where every one of those minutes is
spent.

The phases are the same ones the terminal names, in this order:

| Phase | What it does | What makes it slow |
|---|---|---|
| Obs columns | classifies every served `obs` column, so an unservable dtype fails now rather than inside a later request | many columns; wide categorical fields |
| Embeddings: resolving obsm keys | resolves `X_umap_1d` / `X_umap_2d` / `X_umap_3d`, or a bare `X_umap` at the dimension its column count states | validating and normalizing each embedding over all cells |
| Vector fields | scans `obsm` for `<field>_umap_<dim>d` declarations | one pass per declared field |
| Connectivity | reads and validates `obsp['connectivities']` — **only if you asked for it** | the largest phase by far, when it runs at all |
| Manifests and centroids | one centroid per categorical field per available dimension, then a cross-check of every served route | number of categories × number of dimensions |
| Viewer generation and bind | establishes the configured web generation and binds the port | first run on a machine, or a cleared cache |

To see the same work reported, run the object through a terminal server once —
`serve_anndata(..., quiet=False)` prints the five-step report, and what it says
about your object is what the notebook constructor is doing silently.

## Connectivity is opt-in here too

A notebook viewer does not read `adata.obsp['connectivities']` unless you ask.
The parameter is keyword-only, `False` by default, and accepted by both
{class}`~cellucid.AnnDataViewer` and {func}`~cellucid.show_anndata`:

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(
    adata,
    port=8765,
    dataset_name="Example",
    dataset_id="example",
    serve_connectivity=True,
)
```

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    serve_connectivity=True,
)
```

With it off, the matrix is never read, `dataset_identity.json` reports
`stats.has_connectivity` as `false` and `stats.n_edges` as `null`, and
`connectivity_manifest.json`, `connectivity/edges.src.bin`,
`connectivity/edges.dst.bin`, and `connectivity/edges.weights.f64.bin` all
answer `404`. The overlay is simply not offered.

With it on, the whole graph is read and validated inside the constructor, before
the server binds — the manifest the viewer fetches first declares the edge count
and the neighbor maximum, so there is nothing to defer. That is the phase most
likely to turn a quiet cell into a multi-minute one: a 50-neighbor graph over
millions of cells is hundreds of millions of stored neighbors.

```{important}
Asking for connectivity on an object whose `obsp` holds no `connectivities`
matrix raises during construction rather than producing a graph-less viewer.
Compute it first with `sc.pp.neighbors(adata)`, or leave the parameter alone.
```

The reasoning, the arithmetic, and the terminal transcript are in
{doc}`08_anndata_mode_show_anndata_and_serve_anndata` and
{doc}`../b_concepts_mental_models/07_performance_mental_model_and_scaling`.

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
- `dataset_identity_probes`, keyed by every exact server-declared dataset id
- the configured and browser-facing server URLs
- the verified local viewer-generation status
- recent accepted-event counts grouped by exact event type

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
  Nothing is printed there; see “A notebook viewer builds quietly” above, and
  drop `serve_connectivity=True` if this session will not draw edges.
- The edge overlay is missing although `adata.obsp['connectivities']` exists →
  the graph is opt-in. Rebuild the viewer with `serve_connectivity=True`.
- Construction raised “Connectivity was asked for, and adata.obsp has no
  'connectivities' matrix to serve” → run `sc.pp.neighbors(adata)` first, or
  construct without the parameter.
- Construction raised “No UMAP embedding Cellucid can read was found in
  adata.obsm” → the object declares none of `X_umap_1d`, `X_umap_2d`, or
  `X_umap_3d`, and either has no bare `X_umap` or has one whose width is not 1,
  2, or 3 columns. The message lists the keys present and names that shape;
  {doc}`08_anndata_mode_show_anndata_and_serve_anndata` has the complete rules.
