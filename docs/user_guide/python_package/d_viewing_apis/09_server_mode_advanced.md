# Server mode (advanced)

This page is for readers who want to understand:
- what the server is actually doing,
- how the verified local web generation works,
- what endpoints exist (for debugging),
- and how to operate servers safely (local vs remote).

If you just want “copy/paste and go”, start with {doc}`04_cli_cellucid_serve_quickstart`.

## At a glance

**Audience**
- Computational power users (HPC, SSH tunnels)
- Developers (debugging, integrations)

## Two server implementations (exported vs AnnData)

Cellucid runs one of two servers depending on what you are viewing:

### Exported dataset server (`CellucidServer`)

- serves static files from a directory (fast + cacheable)
- viewer URL looks like:

```text
http://127.0.0.1:<port>/?source=remote
```

The exact `source=remote` marker declares prepared-server launch intent before
the first paint. Open `server.viewer_url` rather than constructing a bare URL;
ordinary sample-catalog onboarding is intentionally not presented for this
explicitly user-served source. This URL opens the server's
`/_cellucid/datasets` catalog. A unique catalog selects its sole dataset; a
multi-dataset catalog requires an exact dataset-id selection. Neither
`CellucidServer.viewer_url` nor the viewer guesses the first catalog entry.

### AnnData server (`AnnDataServer`)

- serves “virtual” Cellucid-format files computed from AnnData
- viewer URL includes:

```text
http://127.0.0.1:<port>/?anndata=true
```

## Verified local UI generation

By default, both servers publish one verified local UI generation:

- before binding the public server, Cellucid establishes every file declared by
  the configured source's `cellucid-web-assets.json`;
- this includes `index.html`, root browser metadata, and `/assets/*`;
- the complete staged generation is byte-verified and then published atomically.

This keeps the viewer UI and dataset API on the same origin, which avoids mixed content and cross-origin problems.
Source or inventory failure stops startup; an older cache is not substituted.

### Cache directory

Default: a temporary directory (platform-dependent).

Override with:

```bash
cellucid serve /path/to/data --web-cache-dir /path/to/cache
```

Python callers pass the same choice explicitly:

```python
from cellucid import serve
serve("/path/to/data", web_cache_dir="/path/to/cache")
```

### Generation publication

The inventory declares the build ID, every path, exact MIME type, byte length,
and SHA-256 hash. Cellucid downloads the complete declared generation to a
sibling staging directory and verifies it before atomic publication. If
download, verification, or publication fails, startup propagates the error and
does not serve another generation.

### Clearing the cache (manual)

```python
from cellucid import clear_web_cache, get_web_cache_dir

print(get_web_cache_dir())
clear_web_cache()
```

## Debug endpoints (exported + AnnData)

These endpoints are useful for “is the server alive” and “what does it think it’s serving?”:

- `/_cellucid/health`
- `/_cellucid/info`
- `/_cellucid/datasets`
- `/_cellucid/protocol`

### Example: health probe

```text
http://127.0.0.1:<port>/_cellucid/health
```

Exported servers return a small “ok + version” payload.
AnnData servers also include:
- `format` (`h5ad`, `zarr`, `in-memory`)
- `is_backed`
- `n_cells`, `n_genes`

### Example: what datasets are visible?

```text
http://127.0.0.1:<port>/_cellucid/datasets
```

This is especially useful when serving a directory that contains multiple export subfolders.

### Example: which wire capabilities does this `cellucid` accept?

```text
http://127.0.0.1:<port>/_cellucid/protocol
```

Two lists — `events` and `commands` — describing what the notebook side of the
protocol accepts. The web build reads this before emitting an event that an
older `cellucid` would reject; it is also the quickest way to confirm which
version of the protocol a running server speaks.

## Frontend ↔ Python communication endpoints

Even in “server mode”, the server includes endpoints used by notebook integrations and hooks:

- `POST /_cellucid/events` (frontend → Python events)
- `POST /_cellucid/session_bundle` (frontend uploads session bundle bytes for notebook workflows)

If you are not using notebook hooks, you can ignore these.

## Serving multiple exported datasets (advanced)

You can run:

```bash
cellucid serve /path/to/exports_root
```

and the server will list subdirectories that look like datasets.

Recommended layout:

```text
exports_root/
  dataset_a/
    dataset_identity.json
    obs_manifest.json
    points_3d.bin.gz
    ...
  dataset_b/
    dataset_identity.json
    obs_manifest.json
    points_2d.bin.gz
    ...
```

## Logging and diagnostics

CLI:

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset -v
```

This enables debug logging for:
- format detection,
- request routing,
- and error traces.

## Security model (summary)

- Default host: `127.0.0.1` (local only).
- Binding to `0.0.0.0` exposes the server on your network.
- CORS is restricted to loopback origins and Cellucid’s hosted origin.

Full discussion: {doc}`13_security_privacy_cors_and_networking`.

## Troubleshooting

- Web-generation startup error → source/inventory/cache-directory diagnosis:
  {doc}`15_troubleshooting_viewing`
- Remote access patterns (SSH) → {doc}`12_remote_servers_ssh_tunneling_and_cloud`
