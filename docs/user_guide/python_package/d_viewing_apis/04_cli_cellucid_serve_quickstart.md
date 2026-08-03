# CLI: `cellucid serve` quickstart

`cellucid serve` is the fastest way to open the Cellucid viewer from a terminal.

It:
- detects your input format (`exported` folder vs `.h5ad` vs `.zarr`),
- starts a small local HTTP server,
- prints a **Viewer URL** you open in your browser,
- and keeps serving data until you stop it.

## At a glance

**Audience**
- Wet lab / non-technical: copy/paste the commands; you don’t need to write Python.
- Computational: focus on backed mode (`.h5ad`) and SSH tunnel workflows.

**Time**
- Local machine: ~5–10 minutes
- Remote/HPC with SSH tunnel: ~15–30 minutes

**Prerequisites**
- `pip install cellucid`
- One of: export directory, `.h5ad`, `.zarr`

## Step-by-step (fast path)

### Step 0 — Install

```bash
pip install cellucid
```

### Step 1 — Start the server (auto-detects format)

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
# or:
# cellucid serve /path/to/data.zarr --dataset-name "My dataset" --dataset-id my-dataset
# cellucid serve /path/to/export_dir
```

You should see a banner like:

```text
CELLUCID SERVER RUNNING
Local URL:  http://127.0.0.1:8765
Viewer URL: http://127.0.0.1:8765/?anndata=true
```

The direct AnnData command above prints `?anndata=true`. A prepared export
prints `?source=remote`. Copy the printed **Viewer URL** exactly; the bare
Local URL is the server origin for diagnostics, not the launch URL.

### Step 2 — Open the viewer in your browser

Copy the **Viewer URL** into your browser address bar.

```{note}
You do not need to open `cellucid.com` manually. Before binding, the Python
server establishes its exact published web generation and then serves it from
the dataset-server origin.
```

### Step 3 — Keep the terminal running

If you close the terminal (or stop the process), the viewer will stop working.

### Step 4 — Stop the server

Press **Ctrl+C** in the terminal.

## What startup prints

Direct AnnData input (`.h5ad` / `.zarr`) runs **five numbered steps**. A prepared
export runs three — `Validating dataset`, `Loading dataset info`,
`Starting server` — because it has no AnnData object to open and its manifests
are already written on disk.

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
      Cells: 1,462,702
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

Step 2 reports the adapter as it is built — obs columns, embeddings, vector
fields, and connectivity when it was asked for — so a long startup names the
phase it is in instead of sitting silent. Step 4 is the manifest and centroid
phase; it is its own numbered step rather than unreported work after a success
line. On a large object nearly all of the wall time is inside step 2.

The `Embeddings:` line names the array each dimension came from, and
`Connectivity:` has three states — served, present but not asked for, or absent.
Both rules are in
{doc}`08_anndata_mode_show_anndata_and_serve_anndata`.

## What the CLI auto-detects (important for debugging)

`cellucid serve <data_path>` decides what you meant using these rules:

- A regular file must use the exact `.h5ad` suffix.
- A Zarr v2 directory must contain both `.zgroup` and `.zattrs` as valid root
  metadata and must not also contain `zarr.json`.
- A directory name ending in `.zarr` is not sufficient by itself, and a valid
  Zarr store does not need that suffix.
- Any other directory that is one complete prepared dataset or contains
  complete prepared dataset subdirectories is treated as **exported**.
- Otherwise → error

If detection fails, the error message will tell you what it expected.

## Most useful flags

Run `cellucid serve --help` anytime.

### Networking / ports

- `--port, -p <int>`: change the port (useful if `8765` is in use)
- `--host, -H <host>`:
  - default: `127.0.0.1`, which no other machine can reach
  - `0.0.0.0`: accept connections on every interface. This is what a compute
    node needs when the tunnel terminates on a different machine, such as a
    cluster login node
  - a non-loopback bind has **no authentication**: every user who can reach the
    port reads the dataset (see {doc}`13_security_privacy_cors_and_networking`)

#### What the banner prints for a wildcard bind

`0.0.0.0`, `::`, and an empty host all name *every* interface. None of them is
an address a browser can open, so none of them is ever printed as a URL. A
server bound to one of them prints the loopback origin as its **Local URL** and
**Viewer URL**, and adds a second block naming the machine itself when the
machine has a hostname that means something outside it:

```text
  Local URL:    http://127.0.0.1:8765
  Viewer URL:   http://127.0.0.1:8765/?anndata=true

  Bound to every network interface. From another machine:
  Machine URL:  http://compute-node.example.org:8765
  Viewer URL:   http://compute-node.example.org:8765/?anndata=true
```

Through an SSH tunnel, open the **loopback** Viewer URL: the tunnel's near end
is on your own laptop, so `127.0.0.1` is right there whatever the server bound
to. The `Machine URL` block is for clients already on the same network. A bind
to an IPv6 literal is bracketed as a URL requires, so `--host ::` prints
`http://[::1]:8765` and `--host ::1` prints the same origin.

### Browser behavior / output

- `--no-browser`: don’t auto-open a browser tab
- `--quiet, -q`: suppress progress banners
- `--verbose, -v`: debug logging (useful when reporting bugs)

### AnnData-specific knobs (only for `.h5ad`/`.zarr`)

- `--latent-key KEY`: choose the latent space in `adata.obsm` used for some derived quantities (e.g. outlier quantiles / centroids)
- `--dataset-name NAME`: exact display name, required for direct AnnData input
- `--dataset-id ID`: exact stable identifier, required for direct AnnData input
- `--obs-key KEY`: serve only this `adata.obs` column; repeat once per column.
  Every column is served when it is omitted. Use it to leave out a column
  Cellucid cannot serve, such as a `datetime64` collection date
- `--vector-field-default FIELD_ID`: exact default field ID; required when
  direct AnnData declares more than one vector field
- `--connectivity`: also serve the `obsp['connectivities']` neighbor graph.
  **Off by default** — see below

Every flag in this group applies to direct AnnData input only. Supplying one
alongside a prepared export directory is an error naming the flag, because an
export already declares its identity, its columns, and its artifacts in
`dataset_identity.json`.

### The graph and the overlay are opt-in

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My study" --dataset-id my-study-v1 \
    --connectivity
```

`--vector-fields` does the same for the `obsm` `<field>_umap_<n>d` overlay
arrays. Both are off by default, because each is read and validated in full
before the server binds and most sessions turn neither overlay on.
AnnData input only: the flag is refused for a prepared export, which already
declares the artifacts it holds.

## Remote server (HPC / cloud): recommended SSH tunnel recipe

This is the safest way to use Cellucid when your data is on a remote machine.

### On the remote machine

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset --no-browser
```

Keep the default host (`127.0.0.1`). Leave this running.

### On your laptop

```bash
ssh -L 8765:localhost:8765 user@remote-host
```

Then open:

```text
http://127.0.0.1:8765/?anndata=true
```

### Serving from an HPC compute node

On a cluster, the tunnel usually cannot terminate on the compute node itself:
you forward through the login node, and the login node opens the second
connection. The server therefore has to accept a connection that does not come
from its own machine, which is what `--host 0.0.0.0` means.

```bash
# on the compute node
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --host 0.0.0.0 --no-browser
```

```bash
# on your laptop
ssh -N -L 8765:<compute-node>:8765 <username>@<login-node>
```

Then open the **loopback** Viewer URL the banner printed, not the `Machine URL`
block:

```text
http://127.0.0.1:8765/?anndata=true
```

```{important}
The bind is `0.0.0.0`; the address you open is still `127.0.0.1`. `--host`
decides which interfaces the server accepts on, and the browser address is
decided by where the tunnel ends — your own laptop. `0.0.0.0` resolves nowhere.
```

Whether you need this at all depends on the cluster, not on Cellucid: the
Helmholtz Munich / ICB cluster and many other sites accept SSH on login nodes
only, so a tunnel has nowhere else to terminate. If your laptop can SSH straight
to the machine holding the data, keep the loopback bind instead.

Full guide (edge cases, firewalls, JupyterHub): {doc}`12_remote_servers_ssh_tunneling_and_cloud`.
Compute nodes, Slurm batch jobs, `tmux`, and the tunnel failure table:
{doc}`17_hpc_slurm_and_compute_node_serving`.

## Debugging endpoints (high signal)

If something looks wrong, open these in a browser (replace `<port>`):

- `http://127.0.0.1:<port>/_cellucid/health` (is the server alive?)
- `http://127.0.0.1:<port>/_cellucid/info` (server config + version)
- `http://127.0.0.1:<port>/_cellucid/datasets` (what datasets does it think exist?)
- `http://127.0.0.1:<port>/_cellucid/protocol` (which wire capabilities does this
  `cellucid` accept?)
- `http://127.0.0.1:<port>/dataset_identity.json` (dataset metadata)
- `http://127.0.0.1:<port>/connectivity_manifest.json` (edge count and neighbor
  maximum — 404 unless the server was started with `--connectivity`)

## Edge cases (common)

- **Port already in use**: the requested port fails. Choose another explicit
  port, or use `--port 0` to request one operating-system-assigned port.
- **`--connectivity` on an object with no graph**: startup fails naming
  `obsp['connectivities']` instead of serving a dataset with no edges.
- **`--connectivity` (or any other AnnData flag) on an export directory**:
  startup fails naming the flag; remove it and run the command again.
- **A bare `X_umap` that is not 1, 2, or 3 columns wide**: startup fails and the
  message names the array's shape.
- **Viewer source blocked**: startup fails instead of substituting a local
  generation; see {doc}`15_troubleshooting_viewing`.
- **Export folder is incomplete**: startup rejects missing required metadata,
  point files, and declared artifacts before binding the server; validate your
  export: {doc}`07_exported_directory_mode_show_and_serve`.

## Troubleshooting

Start here: {doc}`15_troubleshooting_viewing`.
