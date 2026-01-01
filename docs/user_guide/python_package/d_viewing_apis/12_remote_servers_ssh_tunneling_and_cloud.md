# Remote servers, SSH tunneling, and cloud

This page is for when your data is **not on your laptop**:
- HPC clusters,
- lab servers,
- cloud VMs,
- Docker containers.

The safest default is almost always: **keep Cellucid bound to `127.0.0.1` on the remote machine and use an SSH tunnel.**

## At a glance

**Audience**
- Anyone with remote data (HPC/cloud)
- Anyone tempted to use `--host 0.0.0.0`

## Recommended workflow: SSH tunnel (no public exposure)

### Step 1 — Start the server on the remote machine

On the remote machine (SSH session):

```bash
cellucid serve /path/to/data.h5ad --no-browser
```

Keep it bound to `127.0.0.1` (default).

### Step 2 — Forward the port to your laptop

On your laptop:

```bash
ssh -L 8765:localhost:8765 user@remote-host
```

Leave this SSH session open.

### Step 3 — Open the viewer locally

In your local browser:

```text
http://127.0.0.1:8765/
```

```{note}
This works because your browser connects to *your laptop* at `localhost:8765`, and SSH forwards that traffic to the remote machine’s `localhost:8765`.
```

## Common variations

### Remote port is not 8765

If the server prints a different port (or you set one):

```bash
# Remote
cellucid serve /path/to/data.h5ad --port 9000 --no-browser
```

```bash
# Local
ssh -L 9000:localhost:9000 user@remote-host
```

Then open:

```text
http://127.0.0.1:9000/
```

### HPC: jump hosts and compute nodes

If you must hop through a login node, you may need `ProxyJump`:

```bash
ssh -J user@login-node user@compute-node -L 8765:localhost:8765
```

Exact commands vary by cluster policy; the key idea is always “forward a port from the compute node to your laptop”.

## Remote notebooks (JupyterHub / cloud notebooks)

Notebook embedding can fail for two common reasons:

1) **HTTPS mixed content**: your notebook is on `https://…` but the viewer server is `http://…`.
2) **Remote loopback**: the kernel’s `127.0.0.1` is not your browser’s `127.0.0.1`.

Best fixes:

- enable `jupyter-server-proxy` so the iframe URL becomes HTTPS and browser-reachable (recommended)
- or open the viewer in a separate tab using an SSH tunnel (works everywhere)

Deep dive: {doc}`10_notebook_widget_mode_advanced`.

## Docker / containers (advanced)

If you run Cellucid inside a container and want to access it from the host machine:

1) bind the server to `0.0.0.0` inside the container, and
2) publish the port with Docker.

Example:

```bash
cellucid serve /data/export_dir --host 0.0.0.0 --port 8765 --no-browser
```

```{warning}
Binding to `0.0.0.0` exposes the dataset server to the network that can reach it. Prefer SSH tunnels whenever possible.
```

## Troubleshooting

- If the page won’t load at `http://127.0.0.1:<port>/`, start with `/_cellucid/health` and `/_cellucid/info` (see {doc}`15_troubleshooting_viewing`).
- If you’re in a notebook and see “proxy required”, see {doc}`10_notebook_widget_mode_advanced`.
