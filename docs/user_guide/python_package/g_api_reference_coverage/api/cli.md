# CLI (`cellucid`)

The `cellucid` command-line interface is a **single entry point** for common workflows.

Right now, the main user-facing command is:
- the `serve` subcommand, which serves any file or folder Cellucid reads

If you prefer Python APIs, see {doc}`server` and {doc}`jupyter`.

---

## Fast path (beginner-friendly)

```bash
# Serve an AnnData file or store (format detected)
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
cellucid serve /path/to/data.zarr --dataset-name "My dataset" --dataset-id my-dataset
cellucid serve /path/to/data.zarr.zip --dataset-name "My dataset" --dataset-id my-dataset
```

---

## Command reference

### `cellucid --version`

```bash
cellucid --version
```

### `cellucid serve`

```bash
cellucid serve <data_path> [--port 8765] [--host 127.0.0.1] [--no-browser] [--quiet] [--verbose] [dataset and launch options]
```

`<data_path>` can be:
- `.h5ad` file
- Zarr store folder (Zarr v2 with `.zgroup`, or v3 with `zarr.json`), AnnData
  or MuData
- zipped Zarr store (`.zip`, for example `data.zarr.zip`)
- `.h5mu` file
- `.loom` file
- 10x Cell Ranger `.h5` file, or a Cell Ranger `outs/` folder
- MatrixMarket folder (`matrix.mtx[.gz]` with `barcodes` and `features` or
  `genes` lists, also under a shared prefix such as GEO's `GSM4567_`)
- R object: `.rds`, `.rda`, `.RData`, gzipped or not (`.rds.gz`)

An h5Seurat file is not opened, so it is refused before anything is served,
with the route that converts it: in R, read it with
`SeuratDisk::LoadH5Seurat()` and save it with `saveRDS()`, then serve the
`.rds` (or show the object with `cellucid::cellucid_view()`).

A file written by `cellucid.prepare()` declares its own name and id, so
`cellucid serve pbmc_cellucid.h5ad` needs no `--dataset-name`/`--dataset-id`;
an R object, Loom, 10x or MatrixMarket dataset always needs both. Every path
is served as it is on disk; the viewer reads it itself, in place with HTTP
range requests where the format allows it, once from the start for an R object
or a MatrixMarket folder.

#### Options

| Option | Meaning | Notes |
|---|---|---|
| `--port, -p` | Port to bind | Use a different port if 8765 is busy; `--port 0` takes any free port |
| `--host, -H` | Host to bind to | Default `127.0.0.1`, which no other machine can reach. `0.0.0.0` accepts connections on every interface, which is what a compute node needs when a tunnel terminates on a different machine such as a cluster login node. A non-loopback bind has no authentication: every user who can reach the port reads the dataset |
| `--allowed-host HOST` | Extra exact `Host` name this server answers to, and accepts notebook viewer pages from; repeat per name | Needed only behind a reverse proxy; one bare DNS name or IP address, no port, scheme, or wildcard |
| `--no-browser` | Don’t auto-open a browser tab | Helpful on remote servers |
| `--quiet, -q` | Minimal output | Mutually exclusive with `--verbose` |
| `--verbose, -v` | Debug logging | Useful for troubleshooting; also prints the traceback for any failure |
| `--no-web-ui` | Serve the data endpoints without establishing the web application | For machines with no outbound network access |
| `--web-source-url URL` | Exact origin publishing `cellucid-web-assets.json` | |
| `--web-cache-dir PATH` | Directory for the active verified web build | |
| `--latent-key KEY` | The `obsm` latent space (scVI, PCA) the outlier filter measures in | never a 2D/3D layout; without it the file's own `latentKey`, else none (the filter waits for one in Data mapping) |
| `--centroid-min-points N` | The fewest cells a category needs for the outlier filter to measure it | without it the file's own `centroidMinPoints`, else 10; stored quantiles written for another number are computed again in the browser |
| `--dataset-name NAME` | Set the display name | Defaults to the store's `uns['cellucid']['name']`, required without it |
| `--dataset-id ID` | Set the stable dataset identifier | Defaults to the store's `uns['cellucid']['id']`, required without it |
| `--obs-key KEY` | Show this `obs` column; repeat per column, in order | every column is shown when omitted |
| `--vector-field-default FIELD_ID` | The vector field the overlay opens on: `<field>` of an `obsm['<field>_<basis>']` array drawn over `X_<basis>` (`velocity` for `velocity_umap`) | |
| `--modality NAME` | The modality a multimodal dataset (a MuData) opens on | without it the file's own `modality`, else the first gene-expression modality; refused for an AnnData, which has one feature space |
| `--embedding-1d`, `--embedding-2d`, `--embedding-3d KEY[:COLUMNS]` | The `obsm` array (and 0-based columns) of a dimension slot | |
| `--normalize-embeddings`, `--no-normalize-embeddings` | Scale every embedding into the viewer's coordinate range, or show the coordinates the store holds | without either, the file's own `normalizeEmbeddings`, else scaled |
| `--knn KEY` | The `obsp` matrix drawn as the neighbour graph, usually `connectivities` | read only when a view draws it |
| `--gene-id-column COLUMN` | The `var` column (`raw.var` with `--expression raw`) whose values name the genes | the var index when omitted |
| `--expression X\|raw\|layer:NAME` | Where gene expression is read | |

Every AnnData launch option is checked against the store's metadata before the
server binds — a named `obsm`, `obsp`, layer, or column must exist — and handed
to the viewer as a launch default. An option left out is taken from the store's
own `uns['cellucid']['defaults']`, else the viewer applies its default rule
({doc}`../../../web_app/b_data_loading/13_data_mapping`). For a format only the
browser reads (an R object, Loom, a 10x file or folder, a MatrixMarket folder),
the server checks each option's form only (an embedding's columns against its
dimension), and the viewer checks it against the data when it opens it.

A flag overrides the store's own `uns['cellucid']['defaults']` for this server
only; the file is not changed.

---

## Startup output

Startup runs **four numbered steps**: `[1/4] Detecting format`,
`[2/4] Reading metadata`, `[3/4] Checking launch choices`, `[4/4] Starting
server`. The transcript is in {doc}`../../../web_app/b_data_loading/04_server_tutorial`.

Steps 1 to 3 read metadata only, never an array, so they take about the same
time on an atlas as on a small file. Step 2 prints the cell and gene counts and
the `obsm`, `obs`, `obsp`, and `layers` entries of an AnnData or MuData (for a
format only the browser reads it prints `Contents: read by the browser when it
opens the data`); step 3 prints each launch choice it accepted.

### `Embeddings`: which `obsm` key each dimension came from

`--embedding-1d`, `--embedding-2d`, and `--embedding-3d` name the array of a slot,
optionally with 0-based columns (`--embedding-2d X_pca:0,1`). Every other slot
follows the default rule: `X_umap_1d`, `X_umap_2d`, and `X_umap_3d` each fill
their own slot, and when an object declares none of them, the first 1- to
3-column UMAP under a conventional name — `X_umap` (what `sc.tl.umap()` writes),
`X_UMAP` (a SingleCellExperiment's), `X_wnn.umap` or `X_ref.umap` (a Seurat
object's) — fills the slot of its own column count. When neither
fills any slot, the viewer makes a deterministic best guess from the other
`obsm` arrays and reports it as guessed. Nothing is renamed and the object is
not modified. When nothing can be plotted, the viewer reports `Cellucid found
nothing to plot…`, listing every `obsm` entry and why it was not used.

## `--knn KEY`

Names the `obsp` matrix the viewer draws as the neighbour graph. Without it the
viewer draws `obsp['connectivities']` when the store has one. Either way the
graph is read the first time a view draws edges, never at startup.

Naming a matrix the store does not carry is an error before the bind, not
silence — compute it with `sc.pp.neighbors(adata)` first, or drop the flag. A
matrix that is not symmetric, such as `obsp['distances']`, is refused by the
viewer with that reason, never symmetrized.

The Python equivalent is the keyword-only `knn_key="connectivities"` on
{func}`~cellucid.serve_anndata`, {func}`~cellucid.show_anndata`,
{class}`~cellucid.AnnDataViewer`, and {class}`~cellucid.AnnDataServer`; see
{doc}`server`. {func}`~cellucid.prepare` stores one graph,
`obsp['connectivities']`, from its own `connectivities=` argument; the viewer's
default rule finds it without `--knn`.

---

## Practical workflows

### Remote server (HPC) via SSH tunnel

On the remote machine:

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset --no-browser --host 127.0.0.1 --port 8765
```

On your laptop:

```bash
ssh -L 8765:127.0.0.1:8765 user@remote
```

Then open:
- `http://127.0.0.1:8765/?anndata=true`

The Viewer URL ends in `?anndata=true`. Always copy the exact Viewer URL from
the server banner.

### Serving from an HPC compute node

A tunnel to a cluster terminates on the login node, and the login node opens the
second connection to the compute node. The server therefore has to accept a
connection that does not come from its own machine:

```bash
# on the compute node
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --host 0.0.0.0 --no-browser
```

```bash
# on your laptop
ssh -N -L 8765:<compute-node>:8765 <username>@<login-node>
```

Open the **loopback** Viewer URL the banner printed. Full treatment — Slurm
batch jobs, `tmux`, port choice, large objects on shared storage, and the tunnel
failure table: {doc}`../../d_viewing_apis/17_hpc_slurm_and_compute_node_serving`.
Tunnels in general, including firewalls, JupyterHub, and cloud VMs:
{doc}`../../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`.

### What the banner prints for a wildcard bind

`0.0.0.0`, `::`, and an empty host name every interface rather than one address,
and no browser can open any of them, so none is ever printed as a URL. The
banner reports the loopback origin as **Local URL** and **Viewer URL**, and adds
a second block when the machine has a hostname that means something outside it:

A wildcard bind prints the loopback origin plus the machine's own name; see {doc}`../../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`.

An IPv6 literal is bracketed as a URL requires, so a `--host ::` bind prints
`http://[::1]:8765`. Through a tunnel, the loopback URL is the one to open: the
tunnel's near end is on your own laptop. The `Machine URL` block is for clients
already on the same network, and it is absent when the machine's hostname
resolves only to loopback.

---

## Exit codes

- `0`: success (or help text)
- `1`: error (e.g., path not found, format not detected)
- `130`: interrupted by user (Ctrl+C)

## What a failure prints

The last line always begins `Error:` and carries the whole message — what
failed and what to change. Conditions you can correct (a missing flag, a taken
port, an unbindable host, an input that is not a dataset, an optional package
that is not installed) print that line and nothing else.

A Python traceback above the `Error:` line means Cellucid itself failed; that
line then names the issue tracker and asks you to include the traceback. Add
`--verbose` to see the traceback for any failure.

### A missing import names the condition it actually is

An import that fails is not one condition, and the advice differs:

| Condition | What the `Error:` line says |
| --- | --- |
| The distribution is not installed in this environment | `pip install <distribution>`, using the name that installs it rather than the name it imports as |
| The module imported, but a name Cellucid needs is not in it | `pip install --upgrade <distribution>` — installing the same version again would repeat the failure |
| A standard-library module this interpreter cannot provide | The running Python version, the interpreter's own path, and the Python range this Cellucid declares. No `pip` command, because no package installs a standard-library module |

---

## Troubleshooting

### Symptom: “Unable to detect format”
Fix:
- Provide a path that is:
  - a file whose name a format names (`.h5ad`, `.zip`, `.h5mu`, `.loom`, `.h5`,
    `.rds[.gz]`, `.rda[.gz]`, `.RData[.gz]`), or
  - a Zarr store folder (with `.zgroup` or `zarr.json` at its root), or
  - a Cell Ranger `outs/` or MatrixMarket folder.

### Symptom: “This is an h5Seurat file (SeuratDisk), which Cellucid does not open”
Fix:
- In R: `seu <- SeuratDisk::LoadH5Seurat(path)`, then
  `saveRDS(seu, "object.rds")`, and serve the `.rds`.

### Symptom: “Permission denied” / “Address already in use”
Fix:
- Choose a different port: `--port 9000`.

---

### Symptom: “zsh: command not found: cellucid”
Fix:
- Ensure the package is installed in the current environment: `pip install cellucid`
- If you’re inside a virtualenv/conda env, activate it before running `cellucid`.
- The equivalent module invocation is `python -m cellucid.cli --version`.

---

### Symptom: “The server is reachable but the browser didn’t open”
Fix:
- Remove `--no-browser`.
- Or open the printed URL manually.

---

### Symptom: “I need to share this with someone else on my network”
Fix (with caution):
- Bind to `--host 0.0.0.0` and ensure firewall rules are appropriate.
- A non-loopback bind has no authentication: everyone who can reach the port
  reads the dataset.
- Prefer SSH tunneling if you only need access from your own machine.

---

### Symptom: “The browser cannot open `http://0.0.0.0:8765`”
Fix:
- `0.0.0.0` is a bind address, not a destination. Open the Viewer URL the banner
  printed — loopback for a tunnel, the `Machine URL` block from another machine
  on the same network.

---

### Symptom: “The graph is not drawn”
Fix:
- Read the `obsp:` line of step 2. A graph stored under a key other than
  `connectivities` is drawn when you name it with `--knn KEY`.
- No `obsp` entry at all means the object has no graph: compute one with
  `sc.pp.neighbors(adata)` before serving.

---

## See also

- {doc}`server` for the Python server APIs (`serve_anndata`, `AnnDataServer`)
