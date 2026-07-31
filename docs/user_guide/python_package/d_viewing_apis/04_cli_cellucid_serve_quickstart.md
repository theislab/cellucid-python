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
  - default: `127.0.0.1` (local only; safest)
  - `0.0.0.0`: expose on your network (use only if you understand the risk; see {doc}`13_security_privacy_cors_and_networking`)

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

Full guide (edge cases, firewalls, JupyterHub): {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

## Debugging endpoints (high signal)

If something looks wrong, open these in a browser (replace `<port>`):

- `http://127.0.0.1:<port>/_cellucid/health` (is the server alive?)
- `http://127.0.0.1:<port>/_cellucid/info` (server config + version)
- `http://127.0.0.1:<port>/_cellucid/datasets` (what datasets does it think exist?)
- `http://127.0.0.1:<port>/dataset_identity.json` (dataset metadata)

## Edge cases (common)

- **Port already in use**: the requested port fails. Choose another explicit
  port, or use `--port 0` to request one operating-system-assigned port.
- **Viewer source blocked**: startup fails instead of substituting a local
  generation; see {doc}`15_troubleshooting_viewing`.
- **Export folder is incomplete**: startup rejects missing required metadata,
  point files, and declared artifacts before binding the server; validate your
  export: {doc}`07_exported_directory_mode_show_and_serve`.

## Troubleshooting

Start here: {doc}`15_troubleshooting_viewing`.
