# Remote servers, SSH tunneling, and cloud

This page is for when your data is **not on your laptop**:
- HPC clusters,
- lab servers,
- cloud VMs,
- Docker containers.

Two facts decide every recipe on this page:

1. **A tunnel forwards a port to one machine.** Which machine that is decides
   whether the server may stay on `127.0.0.1`.
2. **A socket bound to `127.0.0.1` accepts connections from its own machine
   only.** Nothing else can reach it, including another machine in the same
   cluster.

When your laptop can open an SSH connection **directly to the machine holding
the data**, keep Cellucid on `127.0.0.1` and tunnel to it. When it cannot --
the usual case on an HPC cluster, where only the login node accepts your SSH
connection -- the tunnel terminates somewhere else, and the server has to accept
a connection from that other machine. That is what `--host 0.0.0.0` is for, and
it is covered in full in
{doc}`17_hpc_slurm_and_compute_node_serving`.

## At a glance

**Audience**
- Anyone with remote data (HPC/cloud)
- Anyone tempted to use `--host 0.0.0.0`
- Anyone whose `ssh` to a compute node is refused

**Decide in one question:** can your laptop `ssh` straight to the machine that
will run `cellucid serve`?

| Answer | Bind | Tunnel |
| --- | --- | --- |
| Yes -- lab server, cloud VM, container host | `127.0.0.1` (default) | `ssh -N -L 8765:127.0.0.1:8765 you@remote-host` |
| No -- only a login node accepts my SSH | `0.0.0.0` on the compute node | `ssh -N -L 8765:compute-node:8765 you@login-node` |

## Recommended workflow: SSH tunnel (no public exposure)

### Step 1 — Start the server on the remote machine

On the remote machine (SSH session):

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset --no-browser
```

Keep it bound to `127.0.0.1` (default). `--no-browser` because the remote
machine has no display to open one on.

### Step 2 — Forward the port to your laptop

On your laptop:

```bash
ssh -N -L 8765:127.0.0.1:8765 you@remote-host
```

Leave this SSH session open. `-N` asks for no remote shell, since the only
thing this connection carries is the forwarded port.

Read the forward as three parts: **listen on 8765 on my laptop**, and for each
connection, ask `remote-host` to connect to **`127.0.0.1` port 8765**. That
second address is resolved *by `remote-host`*, which is why the loopback bind is
reachable here and not from anywhere else.

### Step 3 — Open the viewer locally

In your local browser:

```text
http://127.0.0.1:8765/?anndata=true
```

```{note}
This works because your browser connects to *your laptop* at `localhost:8765`, and SSH forwards that traffic to the remote machine’s `localhost:8765`.
Use the exact Viewer URL printed by the remote command. The H5AD command above
prints `?anndata=true`; a prepared export prints `?source=remote`.
```

## Common variations

### Remote port is not 8765

If the server prints a different port (or you set one):

```bash
# Remote
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset --port 9000 --no-browser
```

```bash
# Local
ssh -N -L 9000:127.0.0.1:9000 you@remote-host
```

Then open:

```text
http://127.0.0.1:9000/?anndata=true
```

Pin the port explicitly whenever you tunnel. `--port 0` asks the operating
system for any free port, which is useful on a shared machine but leaves you
forwarding a number you cannot know until after the server prints it.

### HPC: jump hosts and compute nodes

On a cluster, the machine that runs your job is usually **not** the machine that
accepts your SSH connection. Two shapes exist, and they need different commands.

**The cluster lets you SSH to the compute node.** Then `ProxyJump` relays your
own SSH connection through the login node, you authenticate to the compute node
as usual, and the loopback bind is enough:

```bash
ssh -J you@login-node -N -L 8765:127.0.0.1:8765 you@compute-node
```

**The cluster refuses SSH to the compute node.** This is common: many sites
accept SSH on login nodes only, or admit a compute-node connection only while
you hold a running job there, so the attempt ends in
`Permission denied (publickey,keyboard-interactive)` or a timeout. There is then
no SSH session to the compute node at all, and no `ProxyJump` command can make
one. Forward to the compute node *from the login node* instead:

```bash
# On the compute node, inside your job
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --host 0.0.0.0 --no-browser

# On your laptop
ssh -N -L 8765:compute-node:8765 you@login-node
```

The forward target `compute-node:8765` is resolved **by the login node**, which
then opens an ordinary TCP connection to it. That connection is not SSH and
carries no credentials, so it needs no account on the compute node at all -- but
it arrives from another machine, which a socket bound to `127.0.0.1` will not
accept. Hence `--host 0.0.0.0`.

```{important}
The bind is `0.0.0.0`; the address you open is still `http://127.0.0.1:8765/…`.
`--host` decides which interfaces the server accepts connections on. The browser
address is decided by where the tunnel ends, and it ends on your own machine.
`0.0.0.0` means "every interface" and is not an address anything can open.
```

This is a property of the cluster, not of Cellucid: the Helmholtz Munich / ICB
cluster and many others accept SSH on login nodes only, so the tunnel has nowhere
else to terminate. The same shape appears with Docker containers and VMs, where
the loopback interface belongs to the container or the guest rather than to the
machine opening the connection.

```{warning}
`--host 0.0.0.0` accepts connections from every machine that can route to that
port, and Cellucid has no password, token, or any other authentication. On a
shared cluster, every user who can reach the port can read the embeddings, the
obs columns, and the expression values. Serve only while you are looking, and
prefer a loopback bind whenever your client can reach the machine directly.
```

The full treatment -- finding the node your job is on, Slurm job scripts,
surviving a dropped connection, port collisions on shared nodes, and the whole
failure table -- is in {doc}`17_hpc_slurm_and_compute_node_serving`.

### What the banner prints when you bind a wildcard

`0.0.0.0` is a statement about which interfaces accept connections. It is not an
address anything can connect *to*, so Cellucid never prints it inside a URL. A
wildcard bind reports the loopback origin, which is what a tunnel and a browser
on the serving machine both use, followed by the machine's own name for clients
elsewhere on the network:

```text
════════════════════════════════════════════════════════════
  CELLUCID SERVER RUNNING
════════════════════════════════════════════════════════════

  Local URL:    http://127.0.0.1:8765
  Viewer URL:   http://127.0.0.1:8765/?anndata=true

  Bound to every network interface. From another machine:
  Machine URL:  http://compute-node:8765
  Viewer URL:   http://compute-node:8765/?anndata=true

  Press Ctrl+C to stop

════════════════════════════════════════════════════════════
```

Through a tunnel, use the **loopback** Viewer URL: the tunnel's near end is on
your own laptop, so `127.0.0.1` is the correct host there whatever the server
bound to.

## Remote notebooks (JupyterHub / cloud notebooks)

Notebook embedding can fail for two common reasons:

1) **HTTPS mixed content**: your notebook is on `https://…` but the viewer server is `http://…`.
2) **Remote loopback**: the kernel’s `127.0.0.1` is not your browser’s `127.0.0.1`.

Best fixes:

- configure a browser-reachable HTTPS route for one fixed Cellucid port, pass
  its exact base as `client_server_url=`, and declare the route's host name as
  `allowed_hosts=[...]` (a proxy forwards the browser's `Host`, and an
  undeclared authority is refused with `421 Misdirected Request`)
- or open the viewer in a separate tab using an SSH tunnel (works everywhere)

A kernel running on an HPC compute node has the same two-hop problem as the CLI,
with one extra constraint: the embedded viewer always binds its data server to
`127.0.0.1`, and neither `AnnDataViewer` nor `show_anndata` takes a `host`
argument. So the embedded iframe cannot be reached through a login-node forward.
Two things do work there:

- forward the **notebook's own** port instead, so the browser, the kernel, and
  the viewer's loopback all sit on the same machine as far as the connection is
  concerned; or
- start a separate server that can accept the login node's connection, and open
  it in its own browser tab:

  ```python
  from cellucid import serve_anndata

  serve_anndata(
      adata,
      dataset_name="My dataset",
      dataset_id="my-dataset",
      host="0.0.0.0",
      open_browser=False,
  )
  ```

See {doc}`17_hpc_slurm_and_compute_node_serving`.

Deep dive: {doc}`10_notebook_widget_mode_advanced`.

## Docker / containers (advanced)

If you run Cellucid inside a container and want to access it from the host machine:

1) bind the server to `0.0.0.0` inside the container, and
2) publish the port with Docker.

Example:

```bash
cellucid serve /data/export_dir --host 0.0.0.0 --port 8765 --no-browser
```

A container's loopback interface belongs to the container, so a published port
cannot reach a server bound to `127.0.0.1` inside it -- the same reason a
compute node needs the wildcard bind.

```{warning}
Binding to `0.0.0.0` exposes the dataset server to the network that can reach it. Prefer SSH tunnels whenever possible.
```

## Reading a failure

| What you see | What it means | What to do |
| --- | --- | --- |
| `ssh: connect to host … port 22: Operation timed out` | No route, or a firewall dropping the packets | Tunnel through a machine that can reach it |
| `ssh: connect to host … port 22: Connection refused` | The packets arrive; nothing is listening on that port | Check the port; the host is reachable |
| `Permission denied (publickey,keyboard-interactive)` | The host answered and refused your credentials | On a cluster, forward through the login node instead |
| `channel 2: open failed: connect failed: Connection refused` | The tunnel works; the far end has nothing listening | The server is not running, or bound to loopback while the forward names the node |
| Browser cannot open `http://0.0.0.0:8765` | The wildcard is a bind address, not a destination | Open the loopback Viewer URL from the banner |
| Viewer loads, dataset does not | The tunnel reaches a different server | Check the port; something else may hold it |

## Troubleshooting

- If the printed Viewer URL will not load, test the same origin's
  `/_cellucid/health` and `/_cellucid/info` endpoints (see
  {doc}`15_troubleshooting_viewing`).
- If a notebook iframe is blank, mixed-content blocked, or unable to reach
  loopback, see {doc}`10_notebook_widget_mode_advanced`.
- For cluster-specific failures -- a refused compute node, a job that moved, a
  port another user holds -- see
  {doc}`../i_troubleshooting_index/07_hpc_and_remote_access_issues`.
