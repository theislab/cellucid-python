# Compatibility matrix

**Audience:** notebook users, IT-restricted environments, power users  
**Time:** 5–10 minutes  
**Goal:** know what works where (and what to do when it doesn’t)

Cellucid has two ways to show a dataset from Python, and they have different
needs:

- **In a notebook cell** — {func}`~cellucid.show_anndata` returns a
  {class}`~cellucid.jupyter.Widget` (an [anywidget](https://anywidget.dev)).
  The widget frames the Cellucid app, which the browser loads from
  `https://www.cellucid.com`, and connects it to the kernel through the
  notebook's own connection. No port is opened and no proxy is involved.
- **In a browser tab of its own** — {class}`~cellucid.AnnDataViewer`,
  {func}`~cellucid.serve_anndata` and `cellucid serve` start a small server on
  this machine (port 8765 by default) that serves the app and the data. The
  browser must reach that port: on the same machine, through an SSH tunnel, or
  through a reverse proxy.

## Notebook environments

The canonical list, with what “tested” means, is
{doc}`../e_jupyter_hooks/03_supported_environments_matrix`. In short:

| Environment | `show_anndata()` | Notes |
|---|---|---|
| JupyterLab 4 (local or remote, HTTP or HTTPS) | tested | The project's browser suite runs it end to end. |
| Jupyter Notebook 7 | expected, not tested | Same front end as JupyterLab. |
| VS Code notebooks (local, Remote-SSH, containers) | expected, not tested | VS Code renders anywidgets; nothing to forward, since no port is used. |
| Google Colab | expected, not tested | Colab renders anywidgets; `!pip install cellucid` in the runtime. |
| marimo | expected, not tested | Make the widget the cell's output (end the cell with `viewer`). |
| Databricks | expected, not tested | |
| JupyterHub, Binder, cloud notebooks | expected, not tested | No proxy or port setup. |
| A plain Python script | no notebook | Use {class}`~cellucid.AnnDataViewer`, `serve_anndata()` or `cellucid serve`. |

What a notebook needs:

1. **Widget support.** JupyterLab and Notebook 7 need `anywidget` in the
   environment that runs the *notebook server* (installing `cellucid`
   installs it in the kernel's environment, the same one in most setups).
   VS Code, Colab, marimo and Databricks bring their own.
2. **The browser can load the app** from `web_source_url`
   (`https://www.cellucid.com` by default). The browser loads it, not the
   kernel: a compute node without internet is fine, a laptop without internet
   is not. An institution can host a copy and pass `web_source_url=`.

The widget's messages reach Python **between cells**, as every widget's do:
hooks fire and `viewer.state` updates while the kernel is idle, and a call that
waits for the viewer in the same cell (`wait_for_ready()`, `wait_for_event()`,
a session request's `result()`) raises at once with what to do instead;
commands such as `push()` and `set_data_mapping()` never wait. Details: {doc}`../e_jupyter_hooks/03_supported_environments_matrix`.

## Browser-tab viewers (scripts, SSH, HPC)

| Situation | Use | What the browser needs |
|---|---|---|
| Same machine | `cellucid serve data.h5ad` or `AnnDataViewer(...)` | `http://127.0.0.1:8765` |
| Remote machine you can SSH to | the same, with a fixed `port=` | `ssh -N -L 8765:127.0.0.1:8765 you@host` |
| HPC compute node behind a login node | `serve_anndata(..., host="0.0.0.0")` / `cellucid serve --host 0.0.0.0` | `ssh -N -L 8765:compute-node:8765 you@login-node` |
| Behind a reverse proxy (a hub's proxy route) | `AnnDataViewer(..., client_server_url=..., allowed_hosts=[...])`, `serve_anndata(..., allowed_hosts=[...])` | the proxy must pass WebSockets |

Recipes: {doc}`../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud` and
{doc}`../d_viewing_apis/17_hpc_slurm_and_compute_node_serving`.

### Network requirements

| Feature | Network access needed | Notes |
|---|---|---|
| `show_anndata(...)` (notebook widget) | The **browser** loads the app from `web_source_url` | The kernel downloads nothing. |
| `cellucid serve data.h5ad`, `AnnDataViewer`, `serve_anndata` | Yes, at startup (the **kernel/server** machine) | The server downloads, byte-verifies and publishes the app's exact build before it binds. |
| `cellucid serve ... --no-web-ui` | No | Data endpoints only; open them from the hosted app. |
| Hosted web app (`www.cellucid.com`) | Yes | The browser loads the website directly. |
| Hooks/events after a viewer has started | No further download | They travel on the viewer channel. |

## Edge cases (things that look like “bugs”)

- **Ad blockers / strict corporate CSP** can stop the browser from loading
  `https://www.cellucid.com` inside a notebook output; the output says the app
  did not load. Host a copy and pass `web_source_url=`.
- **Port collisions** (browser-tab viewers only): `AnnDataViewer` without
  `port=` takes 8765 or the next free port and prints a note (Google sign-in
  works only on 8765); an explicit busy port fails; `port=0` asks the
  operating system for one.
- **Showing one widget twice**: the newest output takes the viewer over; the
  older one says it was opened elsewhere.
- **Kernel restart**: the output keeps its last view and says the notebook
  stopped it; re-run the cell.
- **Remote access**: binding `--host 0.0.0.0` exposes a browser-tab server to
  the network; do this only on a trusted network.

## Troubleshooting (compatibility)

### Symptom: the notebook output is empty or says it cannot show widgets

**Likely causes**
- `anywidget` is not installed where the notebook *server* runs.
- The notebook front end does not render widgets.

**Fix**
- `pip install anywidget` in the notebook server's environment and reload the
  page; or use JupyterLab, Notebook 7, VS Code, Colab, marimo or Databricks.

### Symptom: the output says the app did not load

**Likely cause:** the browser cannot reach `https://www.cellucid.com` (network
policy, offline laptop).

**Fix:** allow the address, or host a copy of the app and pass
`web_source_url=`.

### Symptom: hooks don’t fire (`on_selection`, `on_click`, …)

**Likely causes**
- The viewer has not loaded yet: it loads once the cell that showed it ends.
- An old viewer is still registered (re-run cells, stale kernel state).

**How to confirm** (in a new cell):

```python
viewer.debug_connection()   # channel_connected, connected_outputs, recent_events
```

**Fix**
- Restart the kernel, or call `viewer.stop()` on old viewers.

## Next steps

- End-to-end quick starts: {doc}`04_quick_start_3_levels`
- Bidirectional hooks deep dive: {doc}`../e_jupyter_hooks/index`
