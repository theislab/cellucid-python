# Viewer lifecycle, cleanup, ports, and multiple viewers

This page is about running Cellucid viewers/servers safely over long sessions:
- cleaning up background servers,
- avoiding port conflicts,
- and understanding what happens when you create multiple viewers.

## At a glance

**Audience**
- Notebook users (re-run cells create new viewers)
- Browser-tab viewer users (most common source of “why is my port in use?”)
- Anyone viewing an in-memory `AnnData` or `MuData` (each viewer holds a
  store of it, with the memory-mapped files a cell-major `X` is laid out in, or
  a temporary `.h5mu`)

## The three lifecycles you might be in

### 1) CLI server lifecycle

- Start: `cellucid serve ...`
- Stop: **Ctrl+C**

If you close the terminal or the process dies, the viewer stops working.

### 2) Python server lifecycle (non-notebook)

- Start (blocking): `serve_anndata(...)`
- Start (non-blocking): `AnnDataServer(...).start_background()`
- Stop: `server.stop()`

### 3) Notebook viewer lifecycle

`show_anndata(...)` returns a notebook widget (`cellucid.Widget`) that shows
the viewer in the cell's output and answers its reads through the notebook's
own connection. It starts no server and opens no port. Stop it explicitly
when you’re done:

```python
viewer.stop()
```

### 4) Browser-tab viewer lifecycle

`AnnDataViewer(...)` starts a server in the background and gives you the
tab's address (`viewer.viewer_url`). Stop it with `viewer.stop()`.

## Why `viewer.stop()` matters

In notebooks, it’s very common to:
- re-run cells,
- create multiple viewer objects,
- and forget that each one holds its dataset (and, for a browser-tab viewer,
  its own server).

Calling `viewer.stop()`:
- freezes the view on its last frame and stops answering it,
- for a browser-tab viewer, stops the server and closes its port,
- and releases what it served: the store of an in-memory `AnnData` (or of a
  file not in priority order) with the files its cell-major `X` was laid out
  in, and the temporary `.h5mu` an in-memory `MuData` was written to (a file is
  read in place and left untouched).

```{note}
Cellucid registers an interpreter-exit cleanup handler that logs shutdown
failures. Call `viewer.stop()` explicitly so your workflow can observe and
handle any failure before interpreter teardown.
```

## Ports: defaults, auto-selection, and “port already in use”

### Defaults

- A notebook widget (`show_anndata`) uses no port.
- CLI/server functions and classes default to the exact port `8765`.
- `AnnDataViewer` defaults to `port=None`: port 8765 (the one local port
  Google sign-in accepts), else the next free port.

### What happens if the port is already in use?

- CLI/server mode binds the requested port exactly and raises if it is unavailable.
- `AnnDataViewer` with `port=None` takes the next free port after 8765 and
  print a note that Google sign-in works only on 8765; a port given is used or
  refused.

**Rule:** use the CLI/server banner for blocking server mode and
`viewer.viewer_url` for `AnnDataViewer` objects; never assume port `8765`.

### Choosing a fixed port (SSH tunnels, Google Drive)

`show_anndata(...)` opens no port. If you need a browser-tab viewer on a
fixed port, instantiate `AnnDataViewer`. Port `8765` is also the one local
port on which Google Drive files open in a browser-tab viewer:
Google answers only the addresses Cellucid registered, port included
({doc}`../../web_app/b_data_loading/16_private_drive_onedrive`).

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(
    "data.h5ad",
    port=8765,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer.display()
```

## Multiple viewers: what’s safe and what’s risky

### Safe patterns

- One viewer per notebook session; stop it before starting another.
- Multiple viewers for side-by-side comparisons (small/medium datasets); each
  notebook widget is independent, each browser-tab viewer has its own port.

### Risky patterns (common footguns)

- Creating many viewers in a loop without calling `stop()` → many open datasets
  (and, for browser-tab viewers, many servers and sockets), possible
  file-handle exhaustion.
- Using `show_anndata(adata, dataset_name=..., dataset_id=...)` on a huge
  in-memory `AnnData` repeatedly → each call lays its cell-major `X` out again
  into files of its own until its viewer is stopped; stop the viewers you no
  longer use.

## Cleanup strategies in notebooks

### Manual “always stop before re-run”

Keep the viewer in a single variable and stop it before creating a new one:

```python
viewer.stop()

viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

### If you lost the reference

If you lost the `viewer` object reference, your safest reset is often:
- `from cellucid.jupyter import cleanup_all; cleanup_all()` (stops every viewer
  of this kernel), or
- restart the notebook kernel.

## Troubleshooting

- Port conflicts, servers stuck, and “viewer won’t load” issues: {doc}`15_troubleshooting_viewing`
- Remote usage and fixed-port workflows: {doc}`12_remote_servers_ssh_tunneling_and_cloud`
