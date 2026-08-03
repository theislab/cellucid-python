# CLI (`cellucid`)

The `cellucid` command-line interface is a **single entry point** for common workflows.

Right now, the main user-facing command is:
- the `serve` subcommand, which auto-detects direct AnnData input or an exported dataset directory

If you prefer Python APIs, see {doc}`server` and {doc}`jupyter`.

---

## Fast path (beginner-friendly)

```bash
# Serve anything (auto-detected)
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
cellucid serve /path/to/data.zarr --dataset-name "My dataset" --dataset-id my-dataset
cellucid serve /path/to/exported_dataset
```

---

## Command reference

### `cellucid --version`

```bash
cellucid --version
```

### `cellucid serve`

```bash
cellucid serve <data_path> [--port 8765] [--host 127.0.0.1] [--no-browser] [--quiet] [--verbose] [--vector-fields] [--connectivity]
```

`<data_path>` can be:
- `.h5ad` file
- `.zarr` directory
- exported dataset directory (contains `dataset_identity.json`)

#### Options

| Option | Meaning | Notes |
|---|---|---|
| `--port, -p` | Port to bind | Use a different port if 8765 is busy; `--port 0` takes any free port |
| `--host, -H` | Host to bind to | Default `127.0.0.1`, which no other machine can reach. `0.0.0.0` accepts connections on every interface, which is what a compute node needs when a tunnel terminates on a different machine such as a cluster login node. A non-loopback bind has no authentication: every user who can reach the port reads the dataset |
| `--allowed-host HOST` | Extra exact `Host` name this server answers to; repeat per name | Needed only behind a reverse proxy; one bare DNS name or IP address, no port, scheme, or wildcard |
| `--no-browser` | Don’t auto-open a browser tab | Helpful on remote servers |
| `--quiet, -q` | Minimal output | Mutually exclusive with `--verbose` |
| `--verbose, -v` | Debug logging | Useful for troubleshooting; also prints the traceback for any failure |
| `--no-web-ui` | Serve the data endpoints without establishing the web application | For machines with no outbound network access |
| `--web-source-url URL` | Exact origin publishing `cellucid-web-assets.json` | |
| `--web-cache-dir PATH` | Directory for the active verified web build | |
| `--latent-key KEY` | Choose the exact latent key in `obsm` | AnnData-only |
| `--dataset-name NAME` | Set the display name | Required for direct AnnData input |
| `--dataset-id ID` | Set the stable dataset identifier | Required for direct AnnData input |
| `--obs-key KEY` | Serve only this `obs` column; repeat per column | AnnData-only; every column is served when omitted |
| `--vector-field-default FIELD_ID` | Choose the initial vector field | AnnData-only; required when more than one exists |
| `--vector-fields` | Serve the `obsm` `<field>_umap_<n>d` overlay arrays | AnnData-only; **off by default** because every declared field is read and checked before the server binds |
| `--connectivity` | Serve the `obsp['connectivities']` neighbor graph | AnnData-only; **off by default** because reading the graph is the longest part of starting a large object |

The AnnData-only flags are rejected for a prepared export directory. The error
names every flag that does not apply, because an export already declares its
identity, its columns, and its artifacts in `dataset_identity.json`.

---

## Startup output

Direct AnnData input runs **five numbered steps**:

Startup prints five numbered steps; the transcript is in {doc}`../../../web_app/b_data_loading/04_server_tutorial`.

Step 2 reports the adapter as it is built, so a slow load names the phase it is
in. Step 4 is the manifest and centroid phase, reported as its own step. When
`--connectivity` is passed, step 2 gains two more lines — `Connectivity: reading
obsp['connectivities']`, then the edge count and neighbor maximum it found.

A prepared export directory runs **three** steps instead — `Validating dataset`,
`Loading dataset info`, `Starting server` — because it has no AnnData object to
open and its manifests are already written.

### `Embeddings`: which `obsm` key each dimension came from

`X_umap_1d`, `X_umap_2d`, and `X_umap_3d` each name their own dimension and
always decide it. When an object declares none of them, a bare `X_umap` — what
`sc.tl.umap()` writes — is read at the dimension its own column count states:
one column as 1D, two as 2D, three as 3D. Nothing is renamed and the object is
not modified. A bare `X_umap` of any other width is refused with a message
naming its shape, and an explicit dimensional key anywhere in `obsm` means the
bare key never joins the resolved set.

### `Connectivity`: three states

| Printed in step 3 | Meaning |
| --- | --- |
| `yes (N edges)` | The graph was read, validated, and is served |
| `not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)` | The object carries a graph and this run did not ask for it |
| `no` | The object carries no `obsp['connectivities']` |

---

## `--vector-fields`

Serves the `<field>_umap_<n>d` arrays in `obsm` as an animated overlay. Without
it, `obsm` is not scanned for them, `dataset_identity.json` declares no
`vector_fields`, and every `/vectors/...` route answers `404`. `--vector-field-default`
selects among served fields, so it is refused without this flag.

An object that declares no field in the `<field>_umap_<n>d` grammar serves none,
which is an ordinary result rather than an error — unlike `--connectivity`, which
names one exact matrix and therefore reports its absence.

## `--connectivity`

Off by default. Without it, `adata.obsp['connectivities']` is never read,
`stats.has_connectivity` in `dataset_identity.json` is `false`, `stats.n_edges`
is `null`, and `connectivity_manifest.json`, `connectivity/edges.src.bin`,
`connectivity/edges.dst.bin`, and `connectivity/edges.weights.f64.bin` all
return 404.

With it, the whole graph is read, validated, and ordered **before the server
binds**, because the manifest the viewer fetches first declares the edge count
and the neighbor maximum; the weights route then carries one Float64 weight per
edge. The wait is real: a 50-neighbor graph over millions of cells is hundreds
of millions of stored neighbors, and deduplicating it into a symmetric, ordered
edge list costs minutes and several times the graph's own memory. At 50,000,000
stored neighbors or more the command logs a warning naming that cost before it
begins validating. Most sessions never draw the graph, which is why the default
is off.

Asking for a graph an object does not carry is an error naming
`obsp['connectivities']`, not silence — compute it with `sc.pp.neighbors(adata)`
first, or drop the flag.

The Python equivalent is the keyword-only `serve_connectivity=True` on
{func}`~cellucid.serve_anndata`, {func}`~cellucid.show_anndata`,
{class}`~cellucid.AnnDataViewer`, and the adapter; see {doc}`server`.
{func}`~cellucid.prepare` is a separate decision made at export time through its
own `connectivities=` argument, and a server reading an export never revisits
it.

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

This H5AD command prints `?anndata=true`; a prepared export prints
`?source=remote`. Always copy the exact Viewer URL from the server banner.

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
  - a real `.h5ad` file, or
  - a real `.zarr` directory, or
  - an exported folder containing `dataset_identity.json`.

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

### Symptom: “The graph is not drawn / `connectivity_manifest.json` returns 404”
Fix:
- Restart with `--connectivity`. The neighbor graph is off by default.
- If startup then fails naming `obsp['connectivities']`, the object has no
  graph: compute one with `sc.pp.neighbors(adata)` before serving.

---

### Symptom: “Startup takes minutes on a large `.h5ad`”
Fix:
- Read the step numbers. Nearly all of the time is inside `[2/5] Loading
  AnnData`, and with `--connectivity` most of that is the edge list.
- Drop `--connectivity` for sessions that never draw the graph.
- For repeated sessions, run {func}`~cellucid.prepare` once and serve the export
  directory instead; see {doc}`export`.

---

## See also

- {doc}`server` for Python server APIs (`serve`, `serve_anndata`)
