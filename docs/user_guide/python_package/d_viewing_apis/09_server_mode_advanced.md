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

## Startup steps and what each one reports

The two servers report different numbers of steps, because they do different
amounts of work before the socket opens.

| Server | Steps |
| --- | --- |
| Exported (`CellucidServer`, `serve`) | `[1/3]` Validating dataset · `[2/3]` Loading dataset info · `[3/3]` Starting server |
| Direct AnnData (`AnnDataServer`, `serve_anndata`, `cellucid serve` on an AnnData path) | `[1/5]` Detecting format · `[2/5]` Loading AnnData · `[3/5]` Analyzing dataset · `[4/5]` Building manifests · `[5/5]` Starting server |

The direct-AnnData path reports **five** steps. Step 2 reports the adapter build
as it happens instead of running silently behind one success line, and the
manifest/centroid phase — previously unannounced work after that success line —
is now its own reported step:

Startup prints five numbered steps; the transcript is in {doc}`../../web_app/b_data_loading/04_server_tutorial`.

`Embeddings: 2D from obsm['X_umap']` names the key that was resolved for each
dimension. `X_umap_1d`, `X_umap_2d`, and `X_umap_3d` each name their own
dimension; an object declaring none of them and carrying a bare `X_umap` has it
read at the dimension its column count states, which is why the report names
both the dimension and the key it came from. Full rule:
{doc}`05_python_serve_and_serve_anndata_quickstart`.

`Connectivity` in the step-3 dataset report has three states, not two:

```text
      Connectivity: yes (74,771 edges)
      Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)
      Connectivity: no
```

The middle line is a graph this object carries that this server was not asked to
serve. Reporting `no` there would read as missing data. Deciding between the
three costs one `in` test against `obsp` — a key lookup, not a read of the
matrix — so naming the state costs nothing that opting out was meant to save.

The prepared-export server is unchanged: an export's manifests already hold the
counts, so three steps still cover it.

## Where the server binds, and what the banner prints

`0.0.0.0`, `::`, and an empty host are not addresses. They tell the operating
system to accept connections on every interface, and nothing can connect *to*
them, so Cellucid never renders one inside a URL. A wildcard-bound server
therefore prints the loopback origin as `Local URL` / `Viewer URL` and reports
the addresses another machine needs separately:

A wildcard bind prints the loopback origin plus the machine's own name; see {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

The same three properties are readable from Python on both server classes:

| Property | Value |
| --- | --- |
| `server.url` | the origin this machine can open — the loopback origin for a wildcard bind, the bound host otherwise |
| `server.viewer_url` | `server.url` plus the launch query (`/?anndata=true` or `/?source=remote`) |
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
`http://127.0.0.1:<port>/?anndata=true` (or `/?source=remote` for a prepared
export), whatever host the far-side server bound to. The wildcard bind exists so
that a second hop — a login node forwarding into a compute node — can connect at
all; it is not the address you type. See
{doc}`12_remote_servers_ssh_tunneling_and_cloud`.
:::

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

## Connectivity routes on the direct-AnnData server (opt-in)

On the direct-AnnData server the neighbor graph is served only when you ask for
it: `serve_connectivity=True` on `AnnDataAdapter`, `AnnDataAdapter.from_file`,
`AnnDataServer`, `serve_anndata`, `AnnDataViewer`, or `show_anndata`, or
`--connectivity` on the command line. The default is off, because building the
edge list is the single longest part of opening a large object — a 50-neighbor
graph over millions of cells is hundreds of millions of stored neighbors — and
most sessions never draw the graph.

| Route | Off (default) | On |
| --- | --- | --- |
| `/connectivity_manifest.json` | `404 No connectivity data` | the weighted edge manifest |
| `/connectivity/edges.src.bin` | `404` | source indices |
| `/connectivity/edges.dst.bin` | `404` | destination indices |
| `/connectivity/edges.weights.f64.bin` | `404` | Float64 weights, one per edge |

With connectivity off, `obsp['connectivities']` is never read, and
`dataset_identity.json` reports `stats.has_connectivity` as `false` with
`stats.n_edges` as `null`. With it on, the whole graph is read and validated
**before** the server binds, because the manifest the viewer fetches first
declares `n_edges` and `max_neighbors`, and the payload routes have to agree
with it byte for byte. Asking for connectivity that `adata.obsp` does not carry
is an error at construction, not a quiet omission.

Two things this does *not* change:

- `prepare()` was already opt-in through its own `connectivities=` argument and
  is untouched.
- The prepared-export server is untouched. An export either holds
  `connectivity_manifest.json` and the three `connectivity/edges.*` artifacts or
  it does not; there is nothing to decide at serve time.

Debug it the same way as any other route — `curl -I` the manifest and read the
status:

```bash
curl -sI http://127.0.0.1:8765/connectivity_manifest.json | head -1
```

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
