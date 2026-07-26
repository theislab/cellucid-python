# Viewer lifecycle, cleanup, ports, and multiple viewers

This page is about running Cellucid viewers/servers safely over long sessions:
- cleaning up background servers,
- avoiding port conflicts,
- and understanding what happens when you create multiple viewers.

## At a glance

**Audience**
- Notebook users (most common source of “why is my port in use?”)
- Anyone serving `.h5ad` in backed mode (file handles matter)

## The three lifecycles you might be in

### 1) CLI server lifecycle

- Start: `cellucid serve ...`
- Stop: **Ctrl+C**

If you close the terminal or the process dies, the viewer stops working.

### 2) Python server lifecycle (non-notebook)

- Start (blocking): `serve(...)` / `serve_anndata(...)`
- Start (non-blocking): `CellucidServer(...).start_background()`, `AnnDataServer(...).start_background()`
- Stop: `server.stop()`

### 3) Notebook viewer lifecycle

`show(...)` / `show_anndata(...)` returns a viewer object that:
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
- closes AnnData file handles (important for `.h5ad` backed mode),
- and frees memory used by internal caches.

```{note}
Cellucid registers an interpreter-exit cleanup handler that logs shutdown
failures. Call `viewer.stop()` explicitly so your workflow can observe and
handle any failure before interpreter teardown.
```

## Ports: defaults, auto-selection, and “port already in use”

### Defaults

- CLI/server functions and classes default to the exact port `8765`.
- Notebook viewer classes default to `port=None`, which requests one
  operating-system-assigned port.

### What happens if the port is already in use?

- CLI/server mode binds the requested port exactly and raises if it is unavailable.
- Notebook viewers use an operating-system-selected free port when `port=None`.

**Rule:** use the CLI/server banner for blocking server mode and
`viewer.viewer_url` for notebook viewer objects; never assume port `8765`.

### Choosing a fixed port (useful for SSH tunnels)

The convenience functions `show(...)` and `show_anndata(...)` do not expose `port` directly.
If you need a fixed port, instantiate the viewer class:

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

```python
from cellucid import CellucidViewer

viewer = CellucidViewer("/path/to/export_dir", port=8765, height=600)
viewer.display()
```

## Multiple viewers: what’s safe and what’s risky

### Safe patterns

- One viewer per notebook session; stop it before starting another.
- Multiple viewers on different ports for side-by-side comparisons (small/medium datasets).

### Risky patterns (common footguns)

- Creating many viewers in a loop without calling `stop()` → lots of servers, lots of open sockets, possible file-handle exhaustion.
- Using `show_anndata(adata, dataset_name=..., dataset_id=...)` on a huge
  in-memory `AnnData` repeatedly → repeated conversions/caches can spike memory.

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
