# Server mode (advanced)

This page is for readers who want to understand:
- what the server is actually doing,
- how the hosted-asset proxy + cache works,
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
http://127.0.0.1:<port>/
```

### AnnData server (`AnnDataServer`)

- serves “virtual” Cellucid-format files computed from AnnData
- viewer URL includes:

```text
http://127.0.0.1:<port>/?anndata=true
```

## Hosted-asset proxy (UI) + cache

By default, both servers run in **hosted-asset proxy mode**:

- when the browser requests `/` or `/index.html`,
  - the server fetches `https://www.cellucid.com/index.html` and `/assets/*` (if needed),
  - caches them under a local cache directory,
  - and serves them as if they were local files.

This keeps the viewer UI and dataset API on the same origin, which avoids mixed content and cross-origin problems.

### Cache directory

Default: a temporary directory (platform-dependent).

Override with:

```bash
export CELLUCID_WEB_PROXY_CACHE_DIR=/path/to/persistent/cache
```

### Cache invalidation

The server looks for a build stamp in `index.html`:

```html
<meta name="cellucid-web-build-id" content="...">
```

If the build ID changes, the cache is purged automatically (best-effort) to avoid stale assets.

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
cellucid serve /path/to/data.h5ad -v
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

- “Viewer UI unavailable” page → cache/network problem: {doc}`15_troubleshooting_viewing`
- Remote access patterns (SSH) → {doc}`12_remote_servers_ssh_tunneling_and_cloud`
