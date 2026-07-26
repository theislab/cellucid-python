# CLI architecture and commands

This page documents the `cellucid` command-line interface: what it does today, how it is structured internally, and how to extend it safely.

---

## CLI goals (design constraints)

The CLI should be:

- **Fast to start** (even on machines without heavy scientific stacks installed).
- **Predictable** (consistent flags, consistent output).
- **Low-friction** for users who just want a viewer URL.

That’s why:
- public imports are lazy (`cellucid.__init__` uses `__getattr__`),
- CLI subcommands import heavy deps only when needed.

---

## Commands (current)

### `cellucid serve`

`serve` is the unified command that **auto-detects** what you gave it:

- `.h5ad` → starts `AnnDataServer` with read-only backed access
- `.zarr` directory → starts `AnnDataServer` after an eager Python load
- directory with `dataset_identity.json` → starts `CellucidServer` (static files)

Examples:

```bash
# Export folder
cellucid serve ./my_export

# AnnData file
cellucid serve ./data.h5ad --dataset-name "My dataset" --dataset-id my-dataset

# Zarr store
cellucid serve ./data.zarr --dataset-name "My dataset" --dataset-id my-dataset
```

Common flags:

- `--port/-p`: bind that exact port (`0` requests one operating-system-assigned port)
- `--host/-H`: bind host (use `0.0.0.0` for remote access)
- `--no-browser`: don’t auto-open a browser tab
- `--quiet/-q`: less output
- `--verbose/-v`: debug logging

AnnData-only flags:

- `--latent-key KEY`: choose latent space for outlier quantiles (AnnData mode)
- `--dataset-name NAME`: exact display name (required in AnnData mode)
- `--dataset-id ID`: exact stable identifier (required in AnnData mode)
- `--vector-field-default FIELD_ID`: exact default field ID (required when
  direct AnnData contains more than one field)

---

## How `serve` auto-detection works (internals)

Implementation lives in:
- `cellucid-python/src/cellucid/cli.py`

Detection rules:

1) A regular file is `h5ad` only when it has the exact `.h5ad` suffix.
2) A directory is `zarr` only when it is a complete Zarr v2 root (both valid
   `.zgroup` and `.zattrs`, with no `zarr.json`).
3) A `.zarr` directory name alone does not establish the format.
4) A non-Zarr directory is `exported` when `_list_exported_datasets(...)`
   discovers one or more complete prepared datasets, including an exports root
   with dataset subdirectories.
5) Else → `unknown` (error).

Then:
- `exported` → call `cellucid.server.serve(...)`
- `h5ad`/`zarr` → call `cellucid.anndata_server.serve_anndata(...)`

---

## Troubleshooting

### Symptom: “`cellucid serve` says ‘Unable to detect format’”

Confirm:
- your path exists,
- a file uses the exact `.h5ad` extension,
- a directory is a complete Zarr v2/v3 root as described above,
- or the directory contains one or more complete prepared datasets.

### Symptom: “The port is busy”

The server reports the bind failure for the exact requested port.

Stop the process using that port, choose another exact port, or request an
operating-system-assigned port with `--port 0`:

```bash
cellucid serve ./my_export --port 0
```

### Symptom: “The browser opens but can’t reach the server”

Common causes:
- you bound to `127.0.0.1` but you’re trying to access remotely,
- you’re on a remote machine without SSH tunneling,
- you’re in an HTTPS notebook environment (mixed-content).

Start with:
- {doc}`09_server_mode_architecture_endpoints_and_security`
- {doc}`10_jupyter_embedding_architecture`
