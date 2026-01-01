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
Cellucid registers a best-effort cleanup handler at interpreter exit, but you should not rely on it during an interactive notebook session.
```

## Ports: defaults, auto-selection, and “port already in use”

### Default port

The default starting port is `8765`.

### What happens if the port is already in use?

- In CLI/server mode, Cellucid checks if the port is free. If not, it auto-selects the next available port and prints a message.
- In notebook mode, the viewer picks a free port automatically (starting at `8765`).

**Rule:** always use the printed URL (or `viewer.viewer_url`), not “I assume it’s 8765”.

### Choosing a fixed port (useful for SSH tunnels)

The convenience functions `show(...)` and `show_anndata(...)` do not expose `port` directly.
If you need a fixed port, instantiate the viewer class:

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer("data.h5ad", port=8765, height=600)
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
- Using `show_anndata(adata)` on a huge in-memory `AnnData` repeatedly → repeated conversions/caches can spike memory.

## Cleanup strategies in notebooks

### Manual “always stop before re-run”

Keep the viewer in a single variable and stop it before creating a new one:

```python
try:
    viewer.stop()
except Exception:
    pass

viewer = show_anndata("data.h5ad")
```

### If you lost the reference

If you lost the `viewer` object reference, your safest reset is often:
- restart the notebook kernel, or
- find the port and kill the process (advanced).

## Troubleshooting

- Port conflicts, servers stuck, and “viewer won’t load” issues: {doc}`15_troubleshooting_viewing`
- Remote usage and fixed-port workflows: {doc}`12_remote_servers_ssh_tunneling_and_cloud`
