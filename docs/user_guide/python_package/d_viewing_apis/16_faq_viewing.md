# FAQ (viewing)

## Do I need to run `prepare()`?

No. Every viewing entry point (`show_anndata`, `serve_anndata`,
`cellucid serve`, the web app's file buttons) opens any AnnData. `prepare()`
writes a Cellucid-ready AnnData: still a standard AnnData, but gene-major,
chunked for range reads, and naming itself — the best choice to share, host,
or open a large dataset often.

See:
- {doc}`../c_data_preparation_api/index`
- {doc}`08_anndata_mode_show_anndata_and_serve_anndata`

## What’s the difference between `serve_anndata` and `show_anndata`?

- `serve_anndata(...)` / `cellucid serve ...`: starts a server and you use the viewer in a normal browser tab.
- `show_anndata(...)`: starts a server and embeds the viewer in a notebook output cell (iframe).

## Does Cellucid upload my data to cellucid.com?

In the recommended workflows (server/notebook modes), your data is served by your local Python server and loaded by your browser from that server.

The Python server may download the **viewer UI assets** from `https://www.cellucid.com` (and cache them locally), but that is UI code, not your dataset.

Security details: {doc}`13_security_privacy_cors_and_networking`.

## Can I use Cellucid offline?

Not for a startup that serves the viewer UI. `cellucid serve`,
`serve_anndata(...)` and `show_anndata(...)` establish the configured exact web generation at each
startup and report source failure instead of substituting a cached generation.

`cellucid serve ... --no-web-ui` is the explicit data-endpoint-only mode and
does not fetch viewer assets.

Cache details: {doc}`09_server_mode_advanced`.

## What happens if port 8765 is already in use?

`serve_anndata(...)` and `cellucid serve` bind the requested port
exactly and report an error if it is unavailable. Pass `port=0` or `--port 0`
to request one operating-system-assigned port.

Always use the printed Viewer URL.

## Why is a remote notebook iframe blank or unreachable?

Usually because:
- the notebook is served over HTTPS, and/or
- the kernel/server is remote (JupyterHub/cloud), so `127.0.0.1` is not reachable from your browser.

Fix: expose the Cellucid port through an HTTPS proxy (for example Jupyter
Server Proxy), pass that exact browser base as `client_server_url=`, and declare
the proxy's host name as `allowed_hosts=[...]`. Without the second argument the
proxied requests are refused with `421 Misdirected Request`, because the proxy
forwards the browser's `Host` and Cellucid answers only authorities that name
it.

Details: {doc}`10_notebook_widget_mode_advanced`.

## Should I use `.h5ad` or Zarr for large datasets?

Both work the same way in Python server mode: the viewer reads the store in
place, with range requests, and the server holds no matrix in memory. A Zarr
store can be a folder (served key by key) or zipped as one `.zip`. What matters
more than the container is the matrix layout: a column-major (CSC) expression
matrix reads a gene directly, while a row-major (CSR) one is indexed by gene
once, the first time a gene is read.

`prepare()` writes the expression matrix column-major, whichever container you
choose. To host a dataset on the web, prefer one file (`.h5ad` or `.zarr.zip`)
to a `.zarr` folder: the browser merges a file's chunk reads into few requests
({doc}`../c_data_preparation_api/10_performance_and_storage`).

Performance guide: {doc}`14_performance_scaling_and_lazy_loading`.

## Is there an R package?

Yes. `cellucid-r` writes the same Cellucid-ready AnnData (`.zarr` or
`.zarr.zip`) from Seurat, SingleCellExperiment, matrices, and data frames. Open
those stores in the web app, serve them with `cellucid serve`, or host them as
static files. Server processes, notebook embedding, and programmatic
viewer control are provided by the Python package.

See {doc}`../../r_package/index`.

## Where do I learn the web app UI (filters, analysis, figure export)?

Those are documented in the web app user guide:

- {doc}`../../web_app/index`
