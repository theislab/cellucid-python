# Viewing methods overview

This page is the “map” of how you can open Cellucid against your data using `cellucid-python`.

If you just want to pick the right command, go to {doc}`03_choose_your_workflow_decision_tree`.

## The mental model (1 minute)

Cellucid always needs two things:

1) the **viewer UI** (HTML/JS) that runs in your browser, and
2) a **dataset server** that can answer “give me points / fields / genes / vectors” requests.

`cellucid-python` starts a **local HTTP server** that serves *both* of those from the same origin (via the hosted-asset proxy), so you can open:

```text
http://127.0.0.1:8765/
```

and the UI can load your data without cross-origin issues.

## At a glance

**Audience**
- Wet lab / non-technical: use the CLI quickstart and stop at “What success looks like”.
- Computational: pay attention to data formats (`exported` vs `.h5ad` vs `.zarr`) and backed mode.
- Power users: pay attention to SSH tunnels, HTTPS notebooks, and caching.

**Prerequisites**
- `pip install cellucid`
- One data input: export directory, `.h5ad`, `.zarr`, or `AnnData`

## The four core entry points

| You want to… | Best entry point | You pass… | Typical place you run it |
|---|---|---|---|
| Open the viewer from the terminal | `cellucid serve …` | a path | terminal |
| Start a server from a Python script | `serve(…)` / `serve_anndata(…)` | export dir or AnnData-ish | Python (notebook optional) |
| Embed the viewer in a notebook | `show(…)` | export dir | Jupyter / VSCode notebooks |
| Embed the viewer for AnnData | `show_anndata(…)` | `.h5ad` / `.zarr` / `AnnData` | Jupyter / VSCode notebooks |

### Server mode vs notebook mode

- **Server mode** (`cellucid serve`, `serve`, `serve_anndata`) is a normal “keep this terminal running” server.
- **Notebook mode** (`show`, `show_anndata`) starts the same kind of server, but also renders an **iframe** in the notebook and gives you a `viewer` object for hooks/commands.

## Exported directory vs AnnData direct

| Data input | How you get it | Recommended APIs | Why you’d choose it |
|---|---|---|---|
| Exported directory | `cellucid.prepare(...)` | `show(dir)` / `serve(dir)` / `cellucid serve dir` | fastest + most reproducible + easiest to share |
| `.h5ad` path | AnnData HDF5 file | `show_anndata("…h5ad")` / `serve_anndata("…h5ad")` / `cellucid serve …h5ad` | convenient + scalable via backed mode |
| `.zarr` path | AnnData Zarr store | `show_anndata("…zarr")` / `serve_anndata("…zarr")` / `cellucid serve …zarr` | convenient + naturally chunked/lazy |
| In-memory `AnnData` | you already loaded it | `show_anndata(adata)` / `serve_anndata(adata)` | great for small/medium; risky for huge matrices |

## Where the “14 loading options” fit

The docs use a canonical list of 14 ways to load data (web app only + server + notebook).

- Full matrix: {doc}`02_the_14_loading_options_breakdown`
- Web-app perspective: {doc}`../../web_app/b_data_loading/01_loading_options_overview`

## Common “gotchas” (read once)

- **Large dataset?** Prefer `exported` + `serve/show` or `.h5ad/.zarr` + `serve_anndata/show_anndata` (backed mode).
- **Remote machine?** Don’t bind to `0.0.0.0` unless you mean to. Use an SSH tunnel first: {doc}`12_remote_servers_ssh_tunneling_and_cloud`.
- **HTTPS notebook?** Direct `http://127.0.0.1:<port>` iframes can be blocked. Use `jupyter-server-proxy` or set `CELLUCID_CLIENT_SERVER_URL`: {doc}`10_notebook_widget_mode_advanced`.
- **Offline?** The viewer UI may need to be cached once while online. See {doc}`09_server_mode_advanced` and {doc}`15_troubleshooting_viewing`.

## Next steps

- Pick a workflow: {doc}`03_choose_your_workflow_decision_tree`
- Terminal flow: {doc}`04_cli_cellucid_serve_quickstart`
- Notebook flow: {doc}`06_jupyter_show_and_show_anndata_quickstart`
