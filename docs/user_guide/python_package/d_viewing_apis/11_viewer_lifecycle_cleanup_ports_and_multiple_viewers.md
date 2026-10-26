# Viewer lifecycle, cleanup, ports, and multiple viewers

This page is about running Cellucid viewers/servers safely over long sessions:
- cleaning up background servers,
- avoiding port conflicts,
- and understanding what happens when you create multiple viewers.

## At a glance

**Audience**
- Notebook users (most common source of “why is my port in use?”)
- Anyone viewing an in-memory `AnnData` or `MuData` (each viewer writes a
  temporary `.h5ad` or `.h5mu`)

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

`show_anndata(...)` returns a viewer object that:
- starts a server in the background, and
- embeds an iframe that points at that server.

Stop it explicitly when you’re done:

```python
viewer.stop()
```

## Why `viewer.stop()` matters

In notebooks, it’s very common to:
- re-run cells,
- create multiple viewer objects,
- and forget that each viewer started its own server.

Calling `viewer.stop()`:
- stops the server,
- closes its port,
- and removes the temporary `.h5ad` (`.h5mu`) an in-memory `AnnData`
  (`MuData`) was written to (a path is served in place and left untouched).

```{note}
Cellucid registers an interpreter-exit cleanup handler that logs shutdown
failures. Call `viewer.stop()` explicitly so your workflow can observe and
handle any failure before interpreter teardown.
```

## Ports: defaults, auto-selection, and “port already in use”

### Defaults

- CLI/server functions and classes default to the exact port `8765`.
- Notebook viewer classes default to `port=None`: port 8765 (the one local
  port Google sign-in accepts), else the next free port.

### What happens if the port is already in use?

- CLI/server mode binds the requested port exactly and raises if it is unavailable.
- Notebook viewers with `port=None` take the next free port after 8765 and
  print a note that Google sign-in works only on 8765; a port given is used or
  refused.

**Rule:** use the CLI/server banner for blocking server mode and
`viewer.viewer_url` for notebook viewer objects; never assume port `8765`.

### Choosing a fixed port (SSH tunnels, Google Drive)

The convenience function `show_anndata(...)` does not expose `port` directly.
If you need a fixed port, instantiate the viewer class. Port `8765` is also
the one local port on which Google Drive files open in a notebook viewer:
Google answers only the addresses Cellucid registered, port included
({doc}`../../web_app/b_data_loading/16_private_drive_onedrive`).

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(
    "data.h5ad",
    port=8765,
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer.display()
```

## Multiple viewers: what’s safe and what’s risky

### Safe patterns

- One viewer per notebook session; stop it before starting another.
- Multiple viewers on different ports for side-by-side comparisons (small/medium datasets).

### Risky patterns (common footguns)

- Creating many viewers in a loop without calling `stop()` → lots of servers, lots of open sockets, possible file-handle exhaustion.
- Using `show_anndata(adata, dataset_name=..., dataset_id=...)` on a huge
  in-memory `AnnData` repeatedly → each call writes another temporary `.h5ad`
  until its viewer is stopped; pass the file's path instead when it is on disk.

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
- restart the notebook kernel, or
- find the port and kill the process (advanced).

## Troubleshooting

- Port conflicts, servers stuck, and “viewer won’t load” issues: {doc}`15_troubleshooting_viewing`
- Remote usage and fixed-port workflows: {doc}`12_remote_servers_ssh_tunneling_and_cloud`
