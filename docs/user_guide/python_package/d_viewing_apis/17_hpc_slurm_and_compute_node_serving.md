# HPC clusters: Slurm, compute nodes, and `--host 0.0.0.0`

This page is the complete recipe for viewing a dataset that lives on an HPC
cluster: the data is on a shared filesystem, the work runs on a compute node
allocated by a scheduler, and the only machine that accepts your SSH connection
is a login node.

It is written for the shape almost every academic cluster has. Names differ --
`login`, `submit`, `head`, `frontend` for the entry machine; `cpusrv…`,
`node…`, `gpu…` for the workers -- but the topology is the same.

## At a glance

**Audience**
- Anyone whose data sits on cluster storage
- Anyone whose `ssh compute-node` is refused
- Anyone who ran `cellucid serve` on a node and could not reach it

**Prerequisites**
- Cellucid installed in an environment on the cluster
- The ability to start an interactive job (`srun`, `salloc`) or submit a batch one
- SSH access to the cluster's login node

**Fast path**

```bash
# 1 — on the login node: get an interactive job and note where it lands
srun --pty --cpus-per-task=4 --mem=32G --time=04:00:00 bash
hostname                       # -> the compute node your job holds

# 2 — on the compute node: serve on every interface, inside tmux
tmux new -s cellucid
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --host 0.0.0.0 --port 8765 --no-browser

# 3 — on your laptop: forward through the login node to that compute node
ssh -N -L 8765:compute-node:8765 you@login-node

# 4 — in your browser
#     http://127.0.0.1:8765/?anndata=true
```

:::{important}
**`--host 0.0.0.0`, and the address is still `127.0.0.1`.**

These two look contradictory and are not. They answer different questions:

- `--host 0.0.0.0` decides **where the server listens**. It has to accept a
  connection that arrives from the login node, which is a different machine, and
  a socket bound to `127.0.0.1` accepts its own machine only.
- `http://127.0.0.1:8765/…` is **what you type in the browser**. The tunnel's
  near end is a port on *your own laptop*, so from the browser's point of view
  the server is local, whatever the far end bound to.

`0.0.0.0` is never an address you open. It means "every interface" and resolves
nowhere. The banner prints the loopback URL for exactly this reason.
:::

## Why this is needed here, and where else it applies

The two-hop forward is not a Cellucid design choice. It is forced by how the
cluster is configured, and the configuration that forces it is common:

- **Helmholtz Munich / ICB.** Compute nodes carry private addresses and do not
  accept a direct SSH connection from outside the cluster, so a tunnel can only
  terminate on a login node. This is the setup the recipe on this page was
  written against.
- **Many other academic and national clusters** apply the same policy, sometimes
  admitting SSH to a compute node only while you hold a running job there.
- **The same shape appears off-cluster**: a Docker container's loopback belongs
  to the container, and a VM's loopback belongs to the VM. Whenever the process
  that opens the second connection is on a different machine from the server, the
  server needs a non-loopback bind.

If your laptop *can* SSH straight to the machine holding the data -- a lab
server, a cloud VM -- none of this applies: keep the loopback bind and use the
direct tunnel in {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

## Why the ordinary recipe does not work here

The usual remote recipe -- bind to `127.0.0.1`, tunnel straight to the machine --
assumes your laptop can open an SSH connection to the machine holding the data.
On a cluster it usually cannot:

- Compute nodes normally carry **private addresses** that are not routable from
  outside the cluster network. The connection times out.
- Many sites run **no sshd on compute nodes at all**, or admit a connection only
  while you hold a running job on that node. The connection is refused with
  `Permission denied (publickey,keyboard-interactive)`.
- Your cluster home directory may not even be the same one the login node uses,
  so installing a key does not necessarily help.

`ssh -J login-node compute-node` does not rescue this. `ProxyJump` still opens a
real SSH session **to the compute node**, using the login node only as a relay,
so everything the compute node refuses, it still refuses.

## What actually works

Forward to the compute node **from the login node**:

```bash
ssh -N -L 8765:compute-node:8765 you@login-node
```

Read the three parts of `-L 8765:compute-node:8765`:

1. Listen on port `8765` **on your laptop**.
2. For each connection, ask the machine you authenticated to -- the **login
   node** -- to open a connection to
3. `compute-node` port `8765`, resolved **from the login node's point of view**.

Step 3 is an ordinary TCP connection, not SSH. It carries no credentials, needs
no account on the compute node, and is unaffected by whatever the compute node's
sshd allows. This is the same mechanism that makes `--ip 0.0.0.0` the standard
advice for Jupyter on a cluster.

Because that connection arrives **from another machine**, a server bound to
`127.0.0.1` will refuse it: loopback accepts its own machine only. So the server
must accept connections on every interface:

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --host 0.0.0.0 --no-browser
```

:::{warning}
`--host 0.0.0.0` has **no authentication of any kind**. Cellucid has no token in
its URL and no password: every user who can route to that port can read the
embeddings, every obs column, and every gene's expression. On a shared cluster
that is every other user on it. Serve only while you are looking, stop the
server when you are done, and prefer a loopback bind whenever your client can
reach the machine directly. See
{doc}`13_security_privacy_cors_and_networking`.
:::

## Step by step

### 1. Get a compute node, and know its name

Never run a server on a login node: those machines are shared by everyone and
usually forbid sustained work. Ask the scheduler for an allocation.

Interactive, for exploring:

```bash
srun --pty --cpus-per-task=4 --mem=32G --time=04:00:00 bash
```

Or a persistent allocation you can attach to more than once:

```bash
salloc --cpus-per-task=4 --mem=32G --time=08:00:00
```

The node name is what your tunnel has to name. From inside the job:

```bash
hostname
echo "$SLURMD_NODENAME"
```

From the login node, at any time:

```bash
squeue -u "$USER" -o "%.10i %.12P %.10j %.8T %.10M %R"
```

The last column is the node your job holds. **This changes every time you get a
new allocation**, so re-read it rather than reusing yesterday's name.

:::{important}
A tunnel names one node. When your job ends and a new one starts elsewhere, the
old tunnel silently points at a machine that is no longer serving anything.
Check `squeue` again after every new allocation.
:::

### 2. Choose a port nobody else is using

Compute nodes are shared, and 8765 is a popular number. Pick something less
contended, and check before you bind:

```bash
ss -ltn | grep 12765 || echo "12765 is free on this node"
```

If the port is taken, `cellucid serve` says so and names the flag to change it.

### 3. Serve, inside a terminal that survives

An interactive job dies with its SSH session. Run the server inside `tmux` (or
`screen`) so a dropped connection does not kill a server that took minutes to
start:

```bash
tmux new -s cellucid
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --host 0.0.0.0 --port 12765 --no-browser
# detach with Ctrl-b d, re-attach later with: tmux attach -t cellucid
```

Startup reports five numbered steps, and on a large object most of the time is
spent inside step 2:

```text
[1/5] Detecting format...
      Path: /path/to/data.h5ad
      Format: h5ad
      ✓ Format detected

[2/5] Loading AnnData...
      Mode: read-only backed h5ad
      Obs columns: classifying 34
      Embeddings: resolving obsm keys
      Embeddings: 2D from obsm['X_umap']
      ✓ File opened

[3/5] Analyzing dataset...
      Cells: 18,142,044
      Genes: 5,001
      Embeddings: 2D
      Obs fields: 22 categorical, 12 continuous
      Vector fields: no
      Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)
      ✓ Analysis complete

[4/5] Building manifests...
      Centroids: one per categorical field and dimension
      ✓ Manifests built

[5/5] Starting server...
      ✓ Server ready
```

### 4. Forward from your laptop

```bash
ssh -N -L 12765:compute-node:12765 you@login-node
```

Leave it running. To make it one command, put the whole thing in
`~/.ssh/config`:

```text
Host cluster-view
    HostName login-node.example.org
    User you
    LocalForward 12765 compute-node:12765
```

then `ssh -N cluster-view`. Update `compute-node` whenever your job moves.

### 5. Open the viewer

Open the **loopback** Viewer URL the banner printed:

```text
http://127.0.0.1:12765/?anndata=true
```

Yes -- `127.0.0.1`, even though you started the server with `--host 0.0.0.0`.
The bind decided which interfaces the server accepts connections on; the browser
address is decided by where the tunnel ends, and it ends on your own laptop. Do
not substitute the compute node's name, and never type `0.0.0.0`: it names every
interface and resolves nowhere.

The banner's `Machine URL` block is the exception -- those URLs are for a client
already inside the cluster network, reaching the node directly without a tunnel.

## Batch jobs

For a long session, submit a job that records where it landed:

```bash
#!/bin/bash
#SBATCH --job-name=cellucid
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=cellucid-%j.log

echo "serving on: $(hostname)"
cellucid serve /path/to/data.h5ad \
    --dataset-name "My dataset" --dataset-id my-dataset \
    --host 0.0.0.0 --port 12765 --no-browser
```

Submit it, then read the node out of the log:

```bash
sbatch serve.sh
grep "serving on" cellucid-*.log
```

The job holds the allocation for its whole walltime, so the server outlives any
SSH session -- and it stops when the walltime ends, which is the honest way to
bound how long an unauthenticated port stays open.

## Large objects on cluster storage

Two things dominate the cost of serving a very large `.h5ad` from shared
storage, and both are now under your control:

- **The neighbor graph is not read unless you ask for it.** A 50-neighbor graph
  over 18 million cells is roughly 900 million stored neighbors; deduplicating
  it into edges, checking it is symmetric, and ordering it takes minutes and
  several times the graph's own memory. `cellucid serve` does none of that by
  default. Pass `--connectivity` (or `serve_connectivity=True`) only when you
  intend to draw the graph. See
  {doc}`08_anndata_mode_show_anndata_and_serve_anndata`.
- **Every gene query reads from the shared filesystem.** Serving a backed
  `.h5ad` over a network filesystem pays that on every request. For repeated
  sessions, run `prepare()` once into an export directory and serve that
  instead; a prepared export needs no `--dataset-name` or `--dataset-id`,
  because it carries its own identity. See
  {doc}`../c_data_preparation_api/01_prepare_export_overview`.

Startup also needs outbound network access to establish the verified web build,
which some compute nodes do not have. `--no-web-ui` serves the data endpoints
without it; see
{doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`.

## Failure table

| Symptom | Cause | Fix |
| --- | --- | --- |
| `ssh: connect to host compute-node port 22: Operation timed out` | The node's address is not routable from outside the cluster | Forward through the login node instead |
| `Permission denied (publickey,keyboard-interactive)` on the compute node | The node refuses SSH, or admits it only with a running job | Forward through the login node; no SSH session to the node is needed |
| Tunnel connects, browser shows nothing | Server bound to `127.0.0.1` while the forward names the node | Add `--host 0.0.0.0` |
| `channel 2: open failed: connect failed: Connection refused` | Nothing is listening on that port on that node | Check the server is still running and on the port you forwarded |
| Browser cannot open `http://0.0.0.0:12765` | The wildcard is a bind address, not a destination | Open the loopback Viewer URL from the banner |
| It worked yesterday, not today | The job moved to another node | `squeue -u "$USER"` and re-point the tunnel |
| `Port 12765 is already in use` | Another user or an earlier server holds it on a shared node | Choose another port, or `--port 0` and read the printed one |
| Startup aborts fetching the web build | The node has no outbound access | `--no-web-ui`, or serve from a node that has it |
| Startup silent for minutes | A very large object is being read | Watch the five steps; step 2 is the long one |

## Related pages

- {doc}`12_remote_servers_ssh_tunneling_and_cloud` -- the general remote recipe
- {doc}`13_security_privacy_cors_and_networking` -- what a non-loopback bind exposes
- {doc}`08_anndata_mode_show_anndata_and_serve_anndata` -- the direct-AnnData path in full
- {doc}`14_performance_scaling_and_lazy_loading` -- where startup time goes
- {doc}`../i_troubleshooting_index/07_hpc_and_remote_access_issues` -- symptom-first triage
