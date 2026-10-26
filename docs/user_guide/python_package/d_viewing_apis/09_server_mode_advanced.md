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

## The server (`AnnDataServer`)

Every Python entry point — `cellucid serve`, {func}`~cellucid.serve_anndata`,
{func}`~cellucid.show_anndata` — runs one server class:


- hands one dataset to the viewer, which reads it itself — in place with HTTP
  range requests where the format allows it, once from the start for an R
  object or a MatrixMarket folder; the server computes nothing per request
- serves the bytes of a file under `/_cellucid/anndata/<dataset_id><ending>`,
  keeping the file's own ending (`.h5ad`, `.zarr.zip`, `.h5mu`, `.loom`,
  `.rds.gz`, …); a Zarr store folder key by key under
  `/_cellucid/anndata/<dataset_id>.zarr/<key>`; and a folder of files (a
  MatrixMarket folder, a Cell Ranger `outs/`) file by file under
  `/_cellucid/anndata/<dataset_id>/<file>`, its files listed in
  `/_cellucid/datasets` because a browser cannot list a folder. `Range`
  requests are answered exactly (`206 Partial Content`, `Content-Range`,
  `HEAD`, `ETag`, `Accept-Ranges`). The route is named by the dataset id, so
  the local file name is never published
- an in-memory `AnnData` is written once to a temporary `.h5ad` (a `MuData` to
  a temporary Cellucid-ready `.h5mu`) and served the same way; the file is
  removed when the server stops
- threaded, with HTTP/1.1 keep-alive, because a browser reads a file with many
  range requests; bytes are sent by the kernel (`sendfile`), and a reload of an
  unchanged file revalidates what the browser kept (`ETag`, `304`) instead of
  downloading it again
- viewer URL includes:

```text
http://127.0.0.1:<port>/?anndata=true
```

The exact `anndata=true` marker declares the user-served launch before the
first paint. Open `server.viewer_url` rather than constructing a bare URL;
ordinary sample-catalog onboarding is intentionally not presented for this
explicitly user-served source. The viewer opens the one dataset the server
lists at `/_cellucid/datasets`.

## Startup steps and what each one reports

`[1/4]` Detecting format · `[2/4]` Reading metadata · `[3/4]` Checking launch
choices · `[4/4]` Starting server.

No array is read before the server binds. It reads the store's metadata (cell
and gene counts, the names of the stored elements), checks each launch choice
you set against it, and binds; the viewer then reads what each view needs. The
transcript is in {doc}`../../web_app/b_data_loading/04_server_tutorial`.

Step 3 lists each launch choice as the viewer receives it. A choice left unset
is not listed: the viewer takes it from the store's
`uns['cellucid']['defaults']`, else its default rule. Full rules:
{doc}`08_anndata_mode_show_anndata_and_serve_anndata`.

## Where the server binds, and what the banner prints

`0.0.0.0`, `::`, and an empty host are not addresses. They tell the operating
system to accept connections on every interface, and nothing can connect *to*
them, so Cellucid never renders one inside a URL. A wildcard-bound server
therefore prints the loopback origin as `Local URL` / `Viewer URL` and reports
the addresses another machine needs separately:

A wildcard bind prints the loopback origin plus the machine's own name; see {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

The same three properties are readable from Python on the server:

| Property | Value |
| --- | --- |
| `server.url` | the origin this machine can open — the loopback origin for a wildcard bind, the bound host otherwise |
| `server.viewer_url` | `server.url` plus the launch query (`/?anndata=true`) |
| `server.network_urls` | `(label, origin)` pairs for the **Bound to every network interface** block; empty for a loopback bind, and empty for a wildcard bind on a machine with no non-loopback name |

Reading the block:

- It appears only for a wildcard bind. A loopback bind has no other machine to
  name.
- It appears only when this machine's own hostname resolves to something other
  than loopback. A hostname that resolves to loopback names nothing outside the
  machine, so no origin is reported rather than one that cannot travel.
- An IPv6 literal is bracketed as a URL requires: `http://[::1]:8765`,
  not `http://::1:8765`.

:::{note}
A viewer reached through an SSH tunnel always uses the **loopback** origin. The
tunnel terminates on your own machine, so the URL you open there is
`http://127.0.0.1:<port>/?anndata=true`, whatever host the far-side server
bound to. The wildcard bind exists so
that a second hop — a login node forwarding into a compute node — can connect at
all; it is not the address you type. See
{doc}`12_remote_servers_ssh_tunneling_and_cloud`.
:::

## Verified local UI generation

By default, the server publishes one verified local UI generation:

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
from cellucid import serve_anndata
serve_anndata("/path/to/data.h5ad", web_cache_dir="/path/to/cache")
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

## Debug endpoints

These endpoints are useful for “is the server alive” and “what does it think it’s serving?”:

- `/_cellucid/health`
- `/_cellucid/info`
- `/_cellucid/datasets`
- `/_cellucid/protocol`
- `/_cellucid/anndata/<dataset_id>.h5ad` (or the file's own ending,
  `.zarr/<key>`, or `<dataset_id>/<file>` for a folder of files): the store's
  bytes
- `/_cellucid/ws`: the viewer channel (a WebSocket), below

### Example: health probe

```text
http://127.0.0.1:<port>/_cellucid/health
```

The server returns:

```json
{"status": "ok", "type": "anndata", "version": "<cellucid version>", "format": "h5ad", "n_cells": 4096, "n_genes": 512}
```

`format` is the input format id of the served store (`h5ad`, `zarr`, `h5mu`,
`mudata-zarr`, `loom`, `tenx-h5`, `tenx-outs`, `rds`, `rda`, `mtx`); an
in-memory `AnnData` reports `h5ad` and an in-memory `MuData` `h5mu`, because
they are served as temporary files. `n_cells` and `n_genes` are `null` for a
format only the browser reads (an R object, Loom, 10x, MatrixMarket).
`/_cellucid/info` carries `version`, `type`, `format`, `n_cells`, `n_genes`,
plus `host`, `port`, and `runtime`: the program serving the viewer (`python`,
or `jupyter` for a notebook viewer), which the viewer names.

### Example: what datasets are visible?

```text
http://127.0.0.1:<port>/_cellucid/datasets
```

The server lists its one dataset with an `anndata` member naming the store and
the launch choices you set:

```json
{"datasets": [{"id": "pbmc-demo", "name": "PBMC demo", "anndata": {"file": {"url": "_cellucid/anndata/pbmc-demo.h5ad", "kind": "h5ad", "layout": "file", "size": 2799976}, "defaults": {"knn": "connectivities"}}}]}
```

`url` is relative to the server's base URL; `kind` is the input format id (as
`format` above); `layout` is `file`, `zip`, or `folder`; `size` is given for a
file and absent for a folder, and a folder of files (MatrixMarket, `outs/`)
lists its `files` with each one's `path` and `size`. `defaults` holds only the
choices you set.

### Example: which wire capabilities does this `cellucid` accept?

```text
http://127.0.0.1:<port>/_cellucid/protocol
```

The viewer channel's route (`channel`, `_cellucid/ws`) and two lists —
`events` and `commands` — describing what the notebook side of the protocol
accepts. The web build reads this before emitting an event that an
older `cellucid` would reject; it is also the quickest way to confirm which
version of the protocol a running server speaks.

## The neighbor graph

The server has no connectivity routes: the viewer reads the
`obsp` graph from the store itself, the first time a view draws edges. Without
`knn_key` (`--knn`) it uses `obsp['connectivities']` when the store has it;
with it, the graph named. Naming a graph the store does not carry stops the
server before it binds.

## Frontend ↔ Python communication endpoints

Even in “server mode”, the server answers the viewer channel used by notebook
integrations and hooks:

- `GET /_cellucid/ws` (a WebSocket): commands to the viewer, events,
  selections and session bundles back, opened only by a registered viewer with
  its token

A server started by `cellucid serve` or `serve_anndata(...)` registers no
viewer, so no channel opens; it answers no `POST` (405). If you are not using
notebook hooks, you can ignore this.

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
- Binding to `0.0.0.0` exposes the server on your network. A non-loopback bind
  has **no authentication**: every user who can reach the port reads the
  dataset. `--host`’s own help text says so, and says why you would still do it
  — a compute node has to accept a connection from the login node when the
  tunnel terminates there.
- The wildcard is never shown inside a URL; the banner prints the loopback
  origin plus a separate **Bound to every network interface** block.
- CORS is restricted to loopback origins and Cellucid’s hosted origin.

Full discussion: {doc}`13_security_privacy_cors_and_networking`.

## Troubleshooting

- Web-generation startup error → source/inventory/cache-directory diagnosis:
  {doc}`15_troubleshooting_viewing`
- Remote access patterns (SSH) → {doc}`12_remote_servers_ssh_tunneling_and_cloud`
