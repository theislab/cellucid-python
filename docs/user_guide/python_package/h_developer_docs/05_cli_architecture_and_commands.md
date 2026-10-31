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

`serve` **detects** the format of the file or folder you gave it, by the input
format registry (see the internals below), and starts `AnnDataServer`, which
serves it as it is on disk:

- `.h5ad`, AnnData Zarr (folder or `.zip`), `.h5mu`, MuData Zarr, `.loom`, a
  10x `.h5` file or Cell Ranger `outs/` folder — read by the viewer in place,
  with range requests;
- a MatrixMarket folder or an R object (`.rds`, `.rda`/`.RData`, gzipped or
  not) — read by the viewer once from the start.

The viewer reads every one itself; the server computes nothing.

Examples:

```bash
# AnnData file
cellucid serve ./data.h5ad --dataset-name "My dataset" --dataset-id my-dataset

# Zarr store
cellucid serve ./data.zarr --dataset-name "My dataset" --dataset-id my-dataset
```

Common flags:

- `--port/-p`: bind that exact port (default `8765`; `0` requests one operating-system-assigned port)
- `--host/-H`: bind host (default `127.0.0.1`; use `0.0.0.0` for remote access)
- `--allowed-host HOST` (repeatable): an extra exact `Host` name to answer to, behind a reverse proxy
- `--cache-bytes BYTES`: what a store served in priority order keeps of what it makes (default a 64th of the machine's memory)
- `--no-browser`: don’t auto-open a browser tab
- `--quiet/-q`: less output
- `--verbose/-v`: debug logging (mutually exclusive with `--quiet`)
- `--no-web-ui`: serve the data routes without establishing the web build
- `--web-source-url URL`, `--web-cache-dir PATH`: where the verified web build comes from and is kept

Dataset flags:

- `--latent-key KEY`: `obsm` latent space (scVI, PCA) the outlier filter
  measures in
- `--dataset-name NAME`: exact display name (default: the store's
  `uns['cellucid']['name']`, else the file's name)
- `--dataset-id ID`: exact stable identifier (default: the store's
  `uns['cellucid']['id']`, else the file's name without its extension)
- `--obs-key KEY` (repeatable), `--embedding-1d/2d/3d KEY[:COLUMNS]`,
  `--knn KEY`, `--gene-id-column COLUMN`, `--expression X|raw|layer:NAME`,
  `--vector-field-default FIELD_ID`, `--[no-]normalize-embeddings`,
  `--centroid-min-points N`, `--modality NAME`: launch choices, each checked
  against the store's metadata before the server binds (for an AnnData or
  MuData; a format only the browser reads is checked when the viewer opens it)
  and handed to the viewer as a launch default; one left out comes from
  `uns['cellucid']['defaults']`, else the viewer's default rule

Startup prints four steps: `[1/4] Detecting format`,
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
   `.zarr.zip`/`.zip`, `.h5mu`, `.loom`, `.h5`, `.rds[.gz]`,
   `.rda`/`.RData[.gz]`, ...). It is then what its first bytes say: the HDF5
   signature (at byte 0, 512, 1024 or 2048) makes it an HDF5 format and the
   ZIP signature a zipped Zarr store, whatever its name; an HDF5 root marked
   `encoding-type` `anndata` is `h5ad` and `mudata` is `h5mu`; otherwise the
   name decides (`loom`, `tenx-h5`, `rds`, `rda`). A `.zip` must hold a Zarr
   store (a `.zgroup`, `.zmetadata` or `zarr.json`), `mudata-zarr` when its
   root says MuData.
2) A directory is a Zarr store when its root is a Zarr v3 (`zarr.json`) or v2
   (`.zgroup`) group (`zarr`, or `mudata-zarr` by its root attributes); else
   what its files make it (`folder_format_of_paths`, the web app's
   `folderFormatOfPaths`): a Cell Ranger `outs/` (`tenx-outs`) or a
   MatrixMarket folder (`mtx`, also with GEO's `GSM…_` prefixes).
3) Else → `unknown` (error: `Unable to detect format for: <path>.` followed by
   the registry's reason, such as which file endings or folder contents it
   reads).

Then every format → `cellucid.anndata_server.serve_anndata(...)`.

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
- a file's name ends as a registry format names it (`.h5ad`, `.zarr.zip`/`.zip`
  holding a Zarr store, `.h5mu`, `.loom`, `.h5`, `.rds[.gz]`,
  `.rda`/`.RData[.gz]`),
- a directory is a Zarr v2/v3 root, a Cell Ranger `outs/` or a MatrixMarket
  folder, as described above.

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
