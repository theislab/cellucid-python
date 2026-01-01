# FAQ (viewing)

## Do I need to run `prepare()`?

No, but it’s often the best choice.

- Use **AnnData mode** (`show_anndata`, `serve_anndata`) when you want convenience in analysis workflows.
- Use **export mode** (`prepare` → `show/serve`) when you want speed, reproducibility, and easy sharing.

See:
- {doc}`07_exported_directory_mode_show_and_serve`
- {doc}`08_anndata_mode_show_anndata_and_serve_anndata`

## What’s the difference between `serve` and `show`?

- `serve(...)` / `cellucid serve ...`: starts a server and you use the viewer in a normal browser tab.
- `show(...)`: starts a server and embeds the viewer in a notebook output cell (iframe).

## Does Cellucid upload my data to cellucid.com?

In the recommended workflows (server/notebook modes), your data is served by your local Python server and loaded by your browser from that server.

The Python server may download the **viewer UI assets** from `https://www.cellucid.com` (and cache them locally), but that is UI code, not your dataset.

Security details: {doc}`13_security_privacy_cors_and_networking`.

## Can I use Cellucid offline?

Yes, after the viewer UI assets have been cached at least once.

If you are offline and have no cached UI, you’ll see a “viewer UI unavailable” page with next steps.

Cache details: {doc}`09_server_mode_advanced`.

## Why did my server pick a different port than 8765?

If `8765` is already in use, Cellucid auto-selects the next free port and prints it.

Always use the printed Viewer URL.

## Why does the notebook embed say “proxy required”?

Usually because:
- the notebook is served over HTTPS, and/or
- the kernel/server is remote (JupyterHub/cloud), so `127.0.0.1` is not reachable from your browser.

Fix: enable `jupyter-server-proxy` (recommended) or set `CELLUCID_CLIENT_SERVER_URL`.

Details: {doc}`10_notebook_widget_mode_advanced`.

## Should I use `.h5ad` or `.zarr` for large datasets?

Both can work well in Python server mode:
- `.h5ad` with backed mode is convenient and widely used.
- `.zarr` is naturally chunked and often performs well for lazy access patterns.

If you’re repeatedly viewing or sharing, exports are usually best.

Performance guide: {doc}`14_performance_scaling_and_lazy_loading`.

## Is there an R package?

`cellucid-R` is planned but not ready yet. Today, the supported “viewing” entry points are:
- CLI (`cellucid serve`)
- Python (`serve`, `serve_anndata`, `show`, `show_anndata`)

## Where do I learn the web app UI (filters, analysis, figure export)?

Those are documented in the web app user guide:

- {doc}`../../web_app/index`
