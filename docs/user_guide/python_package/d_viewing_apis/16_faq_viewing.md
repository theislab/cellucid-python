# FAQ (viewing)

## Do I need to run `prepare()`?

No. Every viewing entry point (`show_anndata`, `serve_anndata`,
`cellucid serve`, the web app's file buttons) opens any AnnData, and every
other format Cellucid reads ({doc}`02_loading_routes`). `prepare()`
writes a Cellucid-ready AnnData: still a standard AnnData, but gene-major,
chunked for range reads, and naming itself — the best choice to share, host,
or open a large dataset often.

See:
- {doc}`../c_data_preparation_api/index`
- {doc}`08_anndata_mode_show_anndata_and_serve_anndata`

## What’s the difference between `serve_anndata` and `show_anndata`?

- `serve_anndata(...)` / `cellucid serve ...`: starts a server and you use the viewer in a normal browser tab.
- `show_anndata(...)`: shows the viewer as a widget in a notebook cell's output; no server, no port. The widget is also the object for hooks, commands and `viewer.pull(adata)`.
- `AnnDataViewer(...)`: the same viewer object for a browser tab of its own (a script, an SSH tunnel); it starts a server.

## Does Cellucid upload my data to cellucid.com?

No. In a notebook, the app (loaded by your browser from
`https://www.cellucid.com`) reads your data from your kernel through the
notebook's own connection. In server mode, your data is served by your local
Python server and loaded by your browser from that server, and the server
downloads the **viewer UI assets** from `https://www.cellucid.com` (and caches
them locally). Either way that is UI code, not your dataset.

Security details: {doc}`13_security_privacy_cors_and_networking`.

## Can I use Cellucid offline?

Not for a startup that serves the viewer UI. `cellucid serve`,
`serve_anndata(...)` and `AnnDataViewer(...)` establish the configured exact
web generation at each startup and report source failure instead of
substituting a cached generation. `show_anndata(...)` needs the browser to
load the app from `https://www.cellucid.com` (or a copy you host, named with
`web_source_url=`).

`cellucid serve ... --no-web-ui` is the explicit data-endpoint-only mode and
does not fetch viewer assets.

Cache details: {doc}`09_server_mode_advanced`.

## What happens if port 8765 is already in use?

`serve_anndata(...)` and `cellucid serve` bind the requested port
exactly and report an error if it is unavailable. Pass `port=0` or `--port 0`
to request one operating-system-assigned port. `AnnDataViewer(...)` without
`port=` takes the next free port instead and prints a note: Google Drive files
open only on 8765 (community annotation's GitHub sign-in works only on
www.cellucid.com). A notebook
widget (`show_anndata(...)`) uses no port.

Always use the printed Viewer URL.

## Does it work in a remote notebook (JupyterHub, VS Code Remote, Colab)?

Yes, with the same call: `show_anndata()` reaches the kernel through the
notebook's own connection, so no port, tunnel or proxy is involved. JupyterLab
is tested by the project; Notebook 7, VS Code, Colab, marimo and Databricks
are expected to work. Details: {doc}`10_notebook_widget_mode_advanced`.

A browser-tab viewer (`AnnDataViewer`, `serve_anndata`) on a remote machine is
not reachable at the kernel's `127.0.0.1`: use an SSH tunnel, or a reverse
proxy whose base you pass as `client_server_url=` and whose host you declare as
`allowed_hosts=[...]` (without it the proxied requests are refused with
`421 Misdirected Request`). See {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

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
`.zarr.zip`, and MuData Zarr from several assays) from Seurat,
SingleCellExperiment, matrices, and data frames. Open those stores in the web
app, serve them with `cellucid serve`, or host them as static files. R also has
its own viewer and server (`cellucid_view()`, `cellucid_serve()`) with the same
commands, events and session capture on the same viewer channel.

See {doc}`../../r_package/index`.

## Where do I learn the web app UI (filters, analysis, figure export)?

Those are documented in the web app user guide:

- {doc}`../../web_app/index`
