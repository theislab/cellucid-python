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

`serve` **detects** which AnnData layout you gave it:

- `.h5ad` file → starts `AnnDataServer`, which serves the file as it is
- Zarr store folder, or zipped Zarr store (`.zip`) → starts `AnnDataServer`,
  which serves the store as it is

The viewer reads an AnnData store itself, with range requests.

Examples:

```bash
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

- `--latent-key KEY`: `obsm` latent space (scVI, PCA) the outlier filter
  measures in
- `--dataset-name NAME`: exact display name (default: the store's
  `uns['cellucid']['name']`; required without it)
- `--dataset-id ID`: exact stable identifier (default: the store's
  `uns['cellucid']['id']`; required without it)
- `--obs-key KEY` (repeatable), `--embedding-1d/2d/3d KEY[:COLUMNS]`,
  `--knn KEY`, `--gene-id-column COLUMN`, `--expression X|raw|layer:NAME`,
  `--vector-field-default FIELD_ID`, `--[no-]normalize-embeddings`,
  `--centroid-min-points N`: launch choices, each checked against the
  store's metadata before the server binds and handed to the viewer as a launch
  default; one left out comes from `uns['cellucid']['defaults']`, else the
  viewer's default rule

Direct AnnData startup prints four steps: `[1/4] Detecting format`,
`[2/4] Reading metadata`, `[3/4] Checking launch choices`,
`[4/4] Starting server`.

---

## How `serve` auto-detection works (internals)

Implementation lives in:
- `cellucid-python/src/cellucid/cli.py`

Detection rules:

`_detect_data_format` asks `_anndata_files.classify_anndata_path`, which
applies the input format registry the web app reads too (the packaged
`input_formats.json`, byte for byte the web app's
`tests/fixtures/input-formats.json`; `tests/test_input_formats.py` pins its
digest). The detected format is the registry's id, which the server states
as the store's `kind` and its health `format`:

1) A regular file is taken when its name ends as a format names it (`.h5ad`,
   `.zip`). It is then what its first bytes say: the HDF5 signature (at byte
   0, 512, 1024 or 2048) makes it `h5ad`, the ZIP signature a zipped `zarr`
   store, whatever its name; the name decides only when the bytes carry no
   signature. A `.zip` must hold a Zarr store (a `.zgroup`, `.zmetadata` or
   `zarr.json`).
2) A directory is `zarr` when its root is a Zarr v3 (`zarr.json`) or v2
   (`.zgroup`) group. A `.zarr` directory name alone does not establish the
   format, and a valid store does not need it.
3) Else → `unknown` (error: `Unable to detect format for: <path>. Expected an
   .h5ad file, a Zarr store folder, or a zipped Zarr store (.zip).`).

Then `h5ad`/`zarr` → `cellucid.anndata_server.serve_anndata(...)`.

---

## How failures are reported (internals)

`main()` classifies every exception that escapes a subcommand into exactly one
of two kinds, and the kind decides whether a traceback is printed.

**Operator conditions** are things the person at the terminal can correct. They
are reported as a single `Error:` line naming what failed and what to change,
with no traceback, and exit `1`. The classification is by exception type, listed
in `cli._OPERATOR_ERROR_TYPES`:

| Type | Reaches the CLI as |
|---|---|
| `ValueError` | a missing identity, an undetectable input, a launch choice the store does not have, an invalid `--allowed-host` |
| `TypeError` | an argument of the wrong type |
| `KeyError` | a named `obs` column, gene, or `obsm` key that is not in the data |
| `RuntimeError` | no browser to open, a viewer-source generation that did not return HTTP 200 |
| `OSError` | a taken port, an unbindable host, a privileged port, an unreadable file — including the `socket.gaierror` and `urllib` `URLError` subclasses |
| `ImportError` | an optional dependency such as `zarr` or `anndata` that is not installed |

These are the types the package's own
{doc}`../g_api_reference_coverage/02_error_messages_and_exceptions_document_patterns`
taxonomy raises deliberately, each already carrying an actionable message; the
CLI adds the flag-level fix for the operating-system conditions, which arrive
with an `errno` and no advice.

**Defects** are everything else — `AttributeError`, `IndexError`,
`AssertionError`, `NameError` and the rest. Cellucid did something it did not
intend, so the full traceback is written to `stderr` unconditionally, even under
`--quiet`, and the `Error:` line names the issue tracker and asks for the
traceback. Losing that traceback would make real bugs unreportable.

Three properties hold for both kinds and are covered by
`tests/test_cli_error_reporting_contract.py`:

- `stdout` is flushed before anything is written to `stderr`, so the `Error:`
  line is last even when the command is piped into a file or a pager;
- the message is collapsed to one line, so “read the last `Error:` line” is
  literally true;
- `--verbose` adds the traceback to an operator condition too, for a bug report
  that turns out to have been misclassified.

The traceback is written with `traceback.print_exception`, not through a logger,
so it does not depend on whether `--quiet` skipped `logging.basicConfig` or on
how an embedding process configured logging.

---

## Troubleshooting

### Symptom: “`cellucid serve` says ‘Unable to detect format’”

Confirm:
- your path exists,
- a file uses the exact `.h5ad` extension, or is a `.zip` holding a Zarr store,
- a directory is a Zarr v2/v3 root as described above.

### Symptom: “The port is busy”

The server reports the bind failure for the exact requested port.

Stop the process using that port, choose another exact port, or request an
operating-system-assigned port with `--port 0`:

```bash
cellucid serve ./data.h5ad --dataset-name "My dataset" --dataset-id my-dataset --port 0
```

### Symptom: “The browser opens but can’t reach the server”

Common causes:
- you bound to `127.0.0.1` but you’re trying to access remotely,
- you’re on a remote machine without SSH tunneling,
- you’re in an HTTPS notebook environment (mixed-content).

Start with:
- {doc}`09_server_mode_architecture_endpoints_and_security`
- {doc}`10_jupyter_embedding_architecture`
