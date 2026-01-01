# CLI (`cellucid`)

The `cellucid` command-line interface is a **single entry point** for common workflows.

Right now, the main user-facing command is:
- `cellucid serve` (auto-detects `.h5ad`, `.zarr`, or an exported dataset directory)

If you prefer Python APIs, see {doc}`server` and {doc}`jupyter`.

---

## Fast path (beginner-friendly)

```bash
# Serve anything (auto-detected)
cellucid serve /path/to/data.h5ad
cellucid serve /path/to/data.zarr
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
cellucid serve <data_path> [--port 8765] [--host 127.0.0.1] [--no-browser] [--quiet] [--verbose]
```

`<data_path>` can be:
- `.h5ad` file
- `.zarr` directory
- exported dataset directory (contains `dataset_identity.json`)

#### Options

| Option | Meaning | Notes |
|---|---|---|
| `--port, -p` | Port to bind | Use a different port if 8765 is busy |
| `--host, -H` | Host/IP to bind | Use `0.0.0.0` only if you want LAN access |
| `--no-browser` | Don’t auto-open a browser tab | Helpful on remote servers |
| `--quiet, -q` | Minimal output | |
| `--verbose, -v` | Debug logging | Useful for troubleshooting |
| `--no-backed` | Load AnnData fully into memory | AnnData-only; disables lazy/backed mode |
| `--latent-key KEY` | Choose latent key in `obsm` | AnnData-only; otherwise auto-detected |

---

## Practical workflows

### Remote server (HPC) via SSH tunnel

On the remote machine:

```bash
cellucid serve /path/to/data.h5ad --no-browser --host 127.0.0.1 --port 8765
```

On your laptop:

```bash
ssh -L 8765:127.0.0.1:8765 user@remote
```

Then open:
- `http://127.0.0.1:8765/`

---

## Exit codes

- `0`: success (or help text)
- `1`: error (e.g., path not found, format not detected)
- `130`: interrupted by user (Ctrl+C)

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
- As a fallback: `python -m cellucid.cli --version`

---

### Symptom: “The server is reachable but the browser didn’t open”
Fix:
- Remove `--no-browser`.
- Or open the printed URL manually.

---

### Symptom: “I need to share this with someone else on my network”
Fix (with caution):
- Bind to `--host 0.0.0.0` and ensure firewall rules are appropriate.
- Prefer SSH tunneling if you only need access from your own machine.

---

## See also

- {doc}`server` for Python server APIs (`serve`, `serve_anndata`)
