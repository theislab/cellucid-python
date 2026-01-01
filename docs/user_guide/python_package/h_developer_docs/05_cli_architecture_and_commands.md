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

- `.h5ad` → starts `AnnDataServer` (lazy loading in backed mode by default)
- `.zarr` directory → starts `AnnDataServer`
- directory with `dataset_identity.json` → starts `CellucidServer` (static files)

Examples:

```bash
# Export folder
cellucid serve ./my_export

# AnnData file (lazy)
cellucid serve ./data.h5ad

# Zarr store
cellucid serve ./data.zarr
```

Common flags:

- `--port/-p`: choose a port (falls back to the next free port if occupied)
- `--host/-H`: bind host (use `0.0.0.0` for remote access)
- `--no-browser`: don’t auto-open a browser tab
- `--quiet/-q`: less output
- `--verbose/-v`: debug logging

AnnData-only flags:

- `--no-backed`: load entire AnnData into memory (debugging only; can be huge)
- `--latent-key KEY`: choose latent space for outlier quantiles (AnnData mode)

---

## How `serve` auto-detection works (internals)

Implementation lives in:
- `cellucid-python/src/cellucid/cli.py`

Detection rules (simplified):

1) If path ends with `.h5ad` → `h5ad`
2) If path ends with `.zarr` → `zarr`
3) If path is a directory and contains `dataset_identity.json` → `exported`
4) Else → `unknown` (error)

Then:
- `exported` → call `cellucid.server.serve(...)`
- `h5ad`/`zarr` → call `cellucid.anndata_server.serve_anndata(...)`

---

## Adding a new CLI command (recommended workflow)

If you add a new subcommand (e.g., `cellucid prepare`, `cellucid cache`, `cellucid session`), do it in a way that preserves the CLI constraints.

Checklist:

1) Add a new `_add_<cmd>_subparser(...)` in `src/cellucid/cli.py`.
2) Keep the parser creation lightweight (no heavy imports).
3) In the command handler (`_run_<cmd>`), import heavy modules *inside* the function.
4) Add `--quiet/--verbose` behavior consistent with the rest of the CLI.
5) Add docs (this file + user guide page if user-facing).
6) Add tests if behavior is non-trivial (see {doc}`13_testing_and_ci`).

---

## Troubleshooting

### Symptom: “`cellucid serve` says ‘Unable to detect format’”

Confirm:
- your path exists,
- you gave the correct extension (`.h5ad`, `.zarr`),
- or the directory is a valid export folder (has `dataset_identity.json`).

### Symptom: “The port is busy”

The server will attempt to pick the next free port automatically.

If you need a specific port, stop the process using it or choose another:

```bash
cellucid serve ./my_export --port 9000
```

### Symptom: “The browser opens but can’t reach the server”

Common causes:
- you bound to `127.0.0.1` but you’re trying to access remotely,
- you’re on a remote machine without SSH tunneling,
- you’re in an HTTPS notebook environment (mixed-content).

Start with:
- {doc}`09_server_mode_architecture_endpoints_and_security`
- {doc}`10_jupyter_embedding_architecture`
