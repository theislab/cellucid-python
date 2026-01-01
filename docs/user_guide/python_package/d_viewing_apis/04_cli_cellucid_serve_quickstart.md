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
cellucid serve /path/to/data.h5ad
# or:
# cellucid serve /path/to/data.zarr
# cellucid serve /path/to/export_dir
```

You should see a banner like:

```text
CELLUCID SERVER RUNNING
Local URL:  http://127.0.0.1:8765
Viewer URL: http://127.0.0.1:8765/
```

<!-- SCREENSHOT PLACEHOLDER
ID: python-cli-server-banner
Suggested filename: data_loading/python_cli_01_server-banner.png
Where it appears: Python Package Guide → Viewing APIs → CLI quickstart → Step 1
Capture:
  - UI location: terminal window
  - State prerequisites: server started successfully
  - Action to reach state: run `cellucid serve ...`
Crop:
  - Include: the printed Viewer URL and the “Press Ctrl+C to stop” line
  - Exclude: usernames, hostnames, private paths, shell history
Redact:
  - Replace file paths with generic `/path/to/data.h5ad`
Alt text:
  - Terminal output showing the Cellucid server banner and viewer URL.
Caption:
  - The viewer URL printed by the CLI is what you open in your browser; keep the terminal running while you use Cellucid.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Cellucid CLI server banner.
:width: 100%

CLI banner showing the local server URL and the browser viewer URL.
```

### Step 2 — Open the viewer in your browser

Copy the **Viewer URL** into your browser address bar.

```{note}
You do not need to visit `cellucid.com` manually for this workflow. The server serves the viewer UI itself (via the hosted-asset proxy).
```

### Step 3 — Keep the terminal running

If you close the terminal (or stop the process), the viewer will stop working.

### Step 4 — Stop the server

Press **Ctrl+C** in the terminal.

## What the CLI auto-detects (important for debugging)

`cellucid serve <data_path>` decides what you meant using these rules:

- If `<data_path>` ends in `.h5ad` and exists → treat as **h5ad**
- If `<data_path>` ends in `.zarr` and exists → treat as **zarr**
- If `<data_path>` is a directory that looks like an export (or contains exported subfolders) → treat as **exported**
- If `<data_path>` is a directory that contains `.zattrs` or `.zgroup` → treat as **zarr**
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

- `--no-backed`: force loading the entire AnnData into memory
  - this is **not recommended** for large datasets
  - it can defeat lazy loading and blow up RAM usage
- `--latent-key KEY`: choose the latent space in `adata.obsm` used for some derived quantities (e.g. outlier quantiles / centroids)

## Remote server (HPC / cloud): recommended SSH tunnel recipe

This is the safest way to use Cellucid when your data is on a remote machine.

### On the remote machine

```bash
cellucid serve /path/to/data.h5ad --no-browser
```

Keep the default host (`127.0.0.1`). Leave this running.

### On your laptop

```bash
ssh -L 8765:localhost:8765 user@remote-host
```

Then open:

```text
http://127.0.0.1:8765/
```

Full guide (edge cases, firewalls, JupyterHub): {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

## Debugging endpoints (high signal)

If something looks wrong, open these in a browser (replace `<port>`):

- `http://127.0.0.1:<port>/_cellucid/health` (is the server alive?)
- `http://127.0.0.1:<port>/_cellucid/info` (server config + version)
- `http://127.0.0.1:<port>/_cellucid/datasets` (what datasets does it think exist?)
- `http://127.0.0.1:<port>/dataset_identity.json` (dataset metadata)

## Edge cases (common)

- **Port already in use**: the server may pick a different port; always copy the printed Viewer URL.
- **Offline / blocked internet**: the viewer UI is cached from `https://www.cellucid.com`; see {doc}`15_troubleshooting_viewing` if you see a “viewer UI unavailable” page.
- **Export folder is incomplete**: the server can start but the viewer may fail when it requests missing files; validate your export: {doc}`07_exported_directory_mode_show_and_serve`.

## Troubleshooting

Start here: {doc}`15_troubleshooting_viewing`.
