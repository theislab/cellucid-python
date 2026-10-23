# CLI: `cellucid serve` quickstart

`cellucid serve` is the fastest way to open the Cellucid viewer from a terminal.

It:
- detects your input format (`.h5ad` file vs Zarr store folder vs zipped Zarr store),
- starts a small local HTTP server,
- prints a **Viewer URL** you open in your browser,
- and keeps serving data until you stop it.

## At a glance

**Audience**
- Wet lab / non-technical: copy/paste the commands; you don’t need to write Python.
- Computational: focus on direct AnnData input (`.h5ad`, Zarr) and SSH tunnel workflows.

**Time**
- Local machine: ~5–10 minutes
- Remote/HPC with SSH tunnel: ~15–30 minutes

**Prerequisites**
- `pip install cellucid`
- One of: `.h5ad`, Zarr store (folder or `.zip`)

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
```

The viewer reads the file itself, in place; the server serves its bytes. `--dataset-name` and `--dataset-id` default to the store's own
`uns['cellucid']['name']` and `['id']`, and are required when it has none. A
file written by `cellucid.prepare()` declares both, so
`cellucid serve pbmc_cellucid.h5ad` is the whole command for one.

You should see a banner like:

```text
CELLUCID SERVER RUNNING
Local URL:  http://127.0.0.1:8765
Viewer URL: http://127.0.0.1:8765/?anndata=true
```

Copy the printed **Viewer URL** (ending in `?anndata=true`) exactly; the bare
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

Startup runs **four numbered steps**.

```text
[1/4] Detecting format...
      Path: /path/to/pbmc_demo.h5ad
      Format: h5ad file
      ✓ Format detected

[2/4] Reading metadata...
      Cells: 4,096
      Genes: 512
      obsm: X_umap
      obs fields: 3
      obsp: connectivities
      layers: none
      ✓ Metadata read

[3/4] Checking launch choices...
      ✓ Launch choices checked

[4/4] Starting server...
      ✓ Server ready
```

Step 2 lists what the store holds, from its metadata; no array is read. Step 3
checks each launch choice you set (`--embedding-2d`, `--knn`, `--expression`,
…) against that metadata and lists it; a name the store does not have stops
the server here, before it binds. A choice you leave unset is taken from the
store's `uns['cellucid']['defaults']`, else the viewer's default rule — see
{doc}`08_anndata_mode_show_anndata_and_serve_anndata`.

## What the CLI auto-detects (important for debugging)

`cellucid serve <data_path>` decides what you meant using these rules:

- A regular file must be an `.h5ad`, or a `.zip` holding a Zarr store (a store
  folder written by `adata.write_zarr`, zipped).
- A Zarr store folder must hold `zarr.json` (Zarr v3) or `.zgroup` (Zarr v2) at
  its root.
- A directory name ending in `.zarr` is not sufficient by itself, and a valid
  Zarr store does not need that suffix.
- Anything else → error: `Unable to detect format for: <path>. Expected an
  .h5ad file, a Zarr store folder, or a zipped Zarr store (.zip).`

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

### Dataset and launch choices

- `--latent-key KEY`: the `adata.obsm` latent space (scVI, PCA) the outlier
  filter measures in, never a 2D/3D layout. Without it the store's own
  `latentKey` is used, else none: the filter then waits for one chosen in
  **Data mapping** ({doc}`../../web_app/e_filtering/02_outlier_filtering_per_active_field`)
- `--dataset-name NAME`: exact display name; defaults to the store's
  `uns['cellucid']['name']`, required when it has none
- `--dataset-id ID`: exact stable identifier; defaults to the store's
  `uns['cellucid']['id']`, required when it has none
- `--obs-key KEY`: show only this `adata.obs` column; repeat once per column, in
  order. Every column is shown when it is omitted. Use it to leave out a column
  Cellucid cannot show, such as a `datetime64` collection date
- `--vector-field-default FIELD_ID`: exact vector field id the overlay opens on
- `--embedding-1d/2d/3d KEY[:COLUMNS]`: the `obsm` array plotted in that
  dimension
- `--normalize-embeddings` / `--no-normalize-embeddings`: scale every embedding
  into the viewer's coordinate range (the default), or show the coordinates the
  store holds
- `--centroid-min-points N`: the fewest cells a category needs for the outlier
  filter to measure it (10 when neither this nor the store says otherwise)
- `--knn KEY`: the `obsp` neighbor graph the viewer draws, read only when a view
  draws it
- `--gene-id-column COLUMN`, `--expression X|raw|layer:NAME`: gene names and
  the expression matrix

A flag overrides the store's own `uns['cellucid']['defaults']` for this
server only; the file is not changed.

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
- `http://127.0.0.1:<port>/_cellucid/datasets` (the dataset id and name, the
  store's route under `_cellucid/anndata/`, and the launch choices)
- `http://127.0.0.1:<port>/_cellucid/protocol` (which wire capabilities does this
  `cellucid` accept?)

## Edge cases (common)

- **Port already in use**: the requested port fails. Choose another explicit
  port, or use `--port 0` to request one operating-system-assigned port.
- **`--knn connectivities` on an object with no graph**: startup fails in step 3
  naming `obsp['connectivities']` instead of serving a dataset with no edges.
- **No `uns['cellucid']` identity and no `--dataset-name`/`--dataset-id`**:
  startup fails asking for both.
- **Viewer source blocked**: startup fails instead of substituting a local
  generation; see {doc}`15_troubleshooting_viewing`.

## Troubleshooting

Start here: {doc}`15_troubleshooting_viewing`.
