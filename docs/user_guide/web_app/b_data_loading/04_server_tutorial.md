# Server Mode (CLI + Python)

Server mode opens `.h5ad` read-only-backed. It also supports `.zarr`, which is
loaded eagerly and therefore must fit server memory.

You run a small local Python server that:
- keeps `.h5ad` matrix access read-only-backed
- materializes `.zarr` input eagerly
- serves only what the viewer needs, on demand
- establishes and serves the exact viewer generation on the same origin

Then you open the exact **Viewer URL printed by the server**:
`/?source=remote` for prepared data or `/?anndata=true` for direct AnnData.

This tutorial covers options **#6–#11** from the “14 loading options” list.

## At A Glance

**Audience**
- Wet lab / beginner: follow the copy/paste commands; you do not need to write code.
- Computational users: focus on lazy loading (`backed` mode), memory, and SSH tunnel workflows.

**Time**
- Local server (same machine): ~5–10 minutes
- Remote server (SSH tunnel): ~15–30 minutes

**Prerequisites**
- `pip install cellucid`
- Your data in one of these forms:
  - pre-exported folder from `prepare()`
  - `.h5ad` file
  - `.zarr` directory

## Security Model (Read Once)

- By default, the server binds to `127.0.0.1` (localhost).
  - This means **only your machine** can access it.
- If you bind to `0.0.0.0`, other machines on your network may be able to access it.
  - Do this only if you understand the risk and you trust your network.

**Best practice for remote machines:**
- Keep the server bound to `127.0.0.1` on the remote machine.
- Use an **SSH tunnel** so you still access it via `http://localhost:...` on your laptop.

This is safer and also avoids browser mixed-content issues.

## Fast Path (CLI)

1) Start the server:

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
```

2) Copy the exact **Viewer URL** printed by the command. With the default port,
   direct AnnData prints:

```text
http://127.0.0.1:8765/?anndata=true
```

3) Keep the terminal running while you use the viewer.
4) Stop the server with **Ctrl+C**.

Real prepared-data startup, against a copy of the standard Pancreas export:

```text
$ cellucid serve ./pancreas_export --no-browser
cellucid v0.9.1

[1/3] Validating dataset...
      Path: /path/to/pancreas_export
      ✓ Dataset valid

[2/3] Loading dataset info...
      Cells: 3,696
      Genes: 3,753
      Connectivity: yes
      ✓ Dataset loaded

[3/3] Starting server...
      Viewer UI generation: establishing exact configured source
Cellucid web UI cache: 100%|██████████████████| 486/486 [00:03<00:00, 155.88file/s]
      ✓ Viewer UI generation established
      ✓ Server ready

════════════════════════════════════════════════════════════
  CELLUCID SERVER RUNNING
════════════════════════════════════════════════════════════

  Local URL:    http://127.0.0.1:8765
  Viewer URL:   http://127.0.0.1:8765/?source=remote

  Press Ctrl+C to stop

════════════════════════════════════════════════════════════
```

Three numbered steps run in order: the directory is validated, its
3,696 cells / 3,753 genes / connectivity are read out of the manifests, and the
486-file web cache is established before the socket opens. The banner is the
part you need — **copy the printed Viewer URL** rather than retyping it.

Two lines vary with your machine rather than with your data. The `Path:` line
prints the **absolute** directory Cellucid resolved, so it will show your own
location, not the `/path/to/…` placeholder above; check it when you are not
certain which copy of an export you just served. The progress bar is sized to
your terminal width, and its rate depends on your network.

:::{note}
The transcript above was produced with the extra flag
`--web-source-url http://127.0.0.1:4183`, which pins the viewer generation to a
checked-out copy of the web app instead of fetching it from
`https://www.cellucid.com`. Everything else — the steps, the counts, the file
count, and the URLs — is what the plain command prints.
:::

### What success looks like in the browser

```{figure} ../../../_static/screenshots/server/open-served-dataset.png
:alt: Cellucid with a dataset loaded from a local Python server; the Session panel names the dataset, reports 4K cells and 33K edges, the Remote server field holds the loopback URL, and Reconnect and Disconnect buttons have replaced Connect.
:width: 1440px

The same server, opened at its printed Viewer URL. Two things tell you the
connection is live rather than a leftover sample: the **Remote server** field
holds the server's address, and a **Disconnect** button has appeared beneath it.
The dataset summary above comes from the server, not from the sample catalog.
```

:::{important}
**The banner prints two URLs, and they are not interchangeable.**

- **Viewer URL** (`http://127.0.0.1:8765/?source=remote`, or `…/?anndata=true`
  for direct AnnData) goes in the **browser address bar**. Opening it is the
  whole workflow; the field below fills itself in.
- `Remote server:` — the text field in the Session panel — takes the bare
  origin the banner prints as **Local URL**: `http://127.0.0.1:8765`. No
  trailing slash, no `?…` query, no credentials, no `#fragment`. The field's
  own placeholder, `http://localhost:8765`, is the exact shape. Pasting the
  Viewer URL there is refused, because of the query string.

You only need the field when you are pointing an already-open Cellucid tab at a
server, for example after a tunnel comes up. Otherwise open the Viewer URL and
skip it.
:::

## Option #6 — Serve a Pre-exported Folder (Best Performance)

Use this if you already ran `prepare()`.

```bash
cellucid serve /path/to/export_dir
```

Notes:
- This is usually the fastest experience.
- The path may be one complete prepared dataset or an exports root containing
  `datasets.json` and multiple complete dataset directories.
- The CLI detects only a complete current export: valid
  `dataset_identity.json`, `obs_manifest.json`, at least one non-empty exact
  points artifact, and every manifest-declared artifact.

## Option #7/#8 — Serve `.h5ad` or `.zarr` Directly (Auto-detected)

```bash
# h5ad
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset

# zarr (directory)
cellucid serve /path/to/data.zarr --dataset-name "My dataset" --dataset-id my-dataset
```

Why this is good:
- `.h5ad` uses read-only backed access instead of loading the full matrix.
- `.zarr` is loaded eagerly, with gene values requested by the browser on demand.

Direct AnnData input runs five numbered steps instead of three, and prints a
different Viewer URL — `/?anndata=true`, not `/?source=remote`:

```text
$ cellucid serve pbmc_demo.h5ad \
      --dataset-name "PBMC demo" --dataset-id pbmc-demo --no-browser --port 8766
cellucid v0.9.1

Importing dependencies (anndata, numpy, scipy)... done

[1/5] Detecting format...
      Path: pbmc_demo.h5ad
      Format: h5ad
      ✓ Format detected

[2/5] Loading AnnData...
      Mode: read-only backed h5ad
      Obs columns: classifying 3
      Embeddings: resolving obsm keys
      Embeddings: 2D from obsm['X_umap']
      ✓ File opened

[3/5] Analyzing dataset...
      Cells: 4,096
      Genes: 512
      Embeddings: 2D
      Obs fields: 2 categorical, 1 continuous
      Vector fields: no
      Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)
      ✓ Analysis complete

[4/5] Building manifests...
      Centroids: one per categorical field and dimension
      ✓ Manifests built

[5/5] Starting server...
      Viewer UI generation: establishing exact configured source
Cellucid web UI cache: 100%|██████████| 493/493 [01:22<00:00,  5.96file/s]
      ✓ Viewer UI generation established
      ✓ Server ready

════════════════════════════════════════════════════════════
  CELLUCID SERVER RUNNING
════════════════════════════════════════════════════════════

  Local URL:    http://127.0.0.1:8766
  Viewer URL:   http://127.0.0.1:8766/?anndata=true

  Press Ctrl+C to stop

════════════════════════════════════════════════════════════
```

Real direct-H5AD startup against a small `.h5ad` holding 4,096 cells, 512 genes,
one `obsm['X_umap']`, and one `obsp['connectivities']`. Nothing about the shape
of the output depends on that size: an atlas prints the same five steps, only
slower.

The numbered steps above are Cellucid's own output. AnnData's library logging
interleaves a few extra `INFO` lines — including the absolute path the file
resolved to — which are omitted here because they repeat what the steps already
report.

Two lines vary with your machine rather than with your data. The `Path:` line
echoes the path you typed, and the progress bar is sized to your terminal width
and paced by your network. Its file count is the size of the published web
generation on the day you run it, so `493/493` above and `486/486` in the
prepared transcript are both real captures of different generations; yours will
be whatever is current.

### What each of the five steps tells you

| Step | What it does | What to read it for |
| --- | --- | --- |
| `[1/5] Detecting format` | Classifies the path as `h5ad`, `zarr`, or a prepared export | The `Format:` line confirms Cellucid agreed with you about what the path is |
| `[2/5] Loading AnnData` | Opens the object and builds the adapter | `Mode: read-only backed h5ad` means the whole matrix is never read into memory; the sub-lines report each phase of the build as it happens |
| `[3/5] Analyzing dataset` | Summarizes what will be served | The sanity check worth reading before you open the browser |
| `[4/5] Building manifests` | Builds and cross-checks every served manifest and every categorical centroid | This phase used to run silently after step 3's success line; it is now reported, because on a large object it is not instant |
| `[5/5] Starting server` | Establishes the exact viewer generation and binds the socket | The banner, and the **Viewer URL** you copy |

Step 2's sub-lines are the ones that matter on a big object, because that is
where the minutes go. `Obs columns` names how many columns are being classified,
and `Embeddings` names the obsm keys as they resolve. The two optional
capabilities print here only when you asked for them: `Vector fields: scanning
obsm` then `Vector fields: N served`, and `Connectivity: reading
obsp['connectivities']` then its edge and neighbor counts. A startup that looks
stalled is stalled inside whichever of those printed last.

Step 3 is the sanity check: if the embeddings, obs-field counts, or connectivity
line disagree with what you expect, fix the file rather than the viewer.

### Which `obsm` keys become embeddings

- `X_umap_1d`, `X_umap_2d`, and `X_umap_3d` each name their own dimension, and
  each one always decides that dimension.
- When an object declares none of those, a bare `X_umap` — the key
  `sc.tl.umap()` writes — is read at the dimension its own column count states:
  1, 2, or 3 columns become the 1D, 2D, or 3D embedding. That is what
  `Embeddings: 2D from obsm['X_umap']` in the transcript above is reporting. No
  key is renamed and your object is not modified.
- A bare `X_umap` of any other width is refused, and the message names its
  shape. Ten columns is a latent space someone named after a plot, and no
  dimension of it is the picture you meant to draw.
- If any dimensional key is present anywhere in `obsm`, the writer declared
  dimensions deliberately, so the bare `X_umap` does not join that set.

```{tip}
Step 2 is the fastest way to answer “does Cellucid see my UMAP?”. An `.h5ad`
that carries no readable embedding fails there, on the terminal, before a
browser is involved — and the failure line lists the `obsm` keys it did find.
```

### Connectivity is opt-in on this path

The neighbor graph is **not** served unless you ask for it. Building the edge
list is the single longest part of startup on a large object — a 50-neighbor
graph over millions of cells is hundreds of millions of stored neighbors — and
most sessions never draw the graph, so `--connectivity` (CLI) and
`serve_connectivity=True` (Python) turn it on:

```bash
cellucid serve /path/to/data.h5ad \
  --dataset-name "My dataset" --dataset-id my-dataset \
  --connectivity
```

With the flag, step 2 reports the read and step 3 reports the count:

```text
[2/5] Loading AnnData...
      Mode: read-only backed h5ad
      Obs columns: classifying 3
      Embeddings: resolving obsm keys
      Embeddings: 2D from obsm['X_umap']
      Connectivity: reading obsp['connectivities']
      Connectivity: 35,167 edges, 26 neighbors at most
      ✓ File opened

[3/5] Analyzing dataset...
      Cells: 4,096
      Genes: 512
      Embeddings: 2D
      Obs fields: 2 categorical, 1 continuous
      Vector fields: no
      Connectivity: yes (35,167 edges)
      ✓ Analysis complete
```

The whole graph is read and validated *before* the socket opens, because the
manifest the viewer fetches first has to declare the edge count and the neighbor
maximum. There is no way to pay that cost later.

The `Connectivity` line in step 3 therefore has three states, not two:

| Line | Meaning |
| --- | --- |
| `Connectivity: yes (35,167 edges)` | The graph was asked for, read, and validated |
| `Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)` | The object carries a graph and this run did not ask for it |
| `Connectivity: no` | The object has no `obsp['connectivities']` at all |

What the flag changes on the wire:

| | `--connectivity` omitted | `--connectivity` passed |
| --- | --- | --- |
| `obsp['connectivities']` | never read | read and validated before the bind |
| `stats.has_connectivity` in `dataset_identity.json` | `false` | `true` |
| `stats.n_edges` | `null` | the validated edge count |
| `connectivity_manifest.json` | 404 | served |
| `connectivity/edges.src.bin` | 404 | served |
| `connectivity/edges.dst.bin` | 404 | served |
| `connectivity/edges.weights.f64.bin` | 404 | served, one Float64 weight per edge |

```{important}
This default applies to the direct-AnnData serve path only. `prepare()` is a
separate decision made at export time through its own `connectivities=`
argument, and a prepared export either holds the connectivity artifacts or does
not — serving one never rebuilds a graph, so there is nothing to opt into.
```

### Vector fields are opt-in on this path too

The velocity/drift overlay takes exactly the same decision as connectivity, for
the same reason: every declared `<field>_umap_<n>d` array in `obsm` is read and
checked before the socket opens. `--vector-fields` (CLI) and
`serve_vector_fields=True` (Python) turn it on:

```bash
cellucid serve /path/to/data.h5ad \
  --dataset-name "My dataset" --dataset-id my-dataset \
  --vector-fields
```

With the flag, step 2 reports the scan (`Vector fields: scanning obsm`, then
`Vector fields: 1 served`) and step 3 names the fields it will serve.

The `Vector fields` line in step 3 has the same three states:

| Line | Meaning |
| --- | --- |
| `Vector fields: yes (velocity_umap)` | The overlay was asked for, read, and validated; the served field IDs are listed |
| `Vector fields: not served (obsm declares vector fields; pass serve_vector_fields=True, or --vector-fields, to draw them)` | The object carries vectors and this run did not ask for them |
| `Vector fields: no` | `obsm` declares no `<field>_umap_<n>d` key at all |

What the flag changes on the wire:

| | `--vector-fields` omitted | `--vector-fields` passed |
| --- | --- | --- |
| `obsm` vector arrays | never read | read and validated before the bind |
| `vector_fields` block in `dataset_identity.json` | absent | present, listing every served field |
| `vectors/<index>_<dim>d.bin` | 404 | served |

`--vector-field-default` chooses among *served* fields, so passing it on its own
raises rather than doing nothing. And as with connectivity, this default applies
to the direct-AnnData path only: a prepared export written by
`prepare(vector_fields=...)` already holds the vectors under `vectors/` and
declares them in `dataset_identity.json`, so there is nothing to opt into.

Quick verification (high signal):
- Open `http://127.0.0.1:8765/dataset_identity.json` and search for
  `vector_fields`. Absent means either “not asked for” or “not present” — step 3
  of the startup output distinguishes the two.

If the block is there but the overlay dropdown is empty, it is a dimension
mismatch: you have 2D vectors and are viewing 3D, or the reverse.

See:
- {doc}`../i_vector_field_velocity/index`
- {doc}`../i_vector_field_velocity/01_what_vector_fields_are_user_facing`

## CLI Options (What They Mean)

Run `cellucid serve --help` anytime. Its complete output:

```text
$ cellucid serve --help
usage: cellucid serve [-h] [--port PORT] [--host HOST] [--allowed-host HOST]
                      [--no-browser] [--quiet | --verbose] [--no-web-ui]
                      [--web-source-url URL] [--web-cache-dir PATH]
                      [--latent-key KEY] [--dataset-name NAME]
                      [--dataset-id ID] [--obs-key KEY]
                      [--vector-field-default FIELD_ID] [--vector-fields]
                      [--connectivity]
                      data_path

Serve data for visualization. Automatically detects format.

positional arguments:
  data_path             Path to h5ad file, zarr directory, or pre-exported
                        dataset

options:
  -h, --help            show this help message and exit
  --port PORT, -p PORT  Port to serve on (default: 8765)
  --host HOST, -H HOST  Host to bind to (default: 127.0.0.1, which no other
                        machine can reach). Use 0.0.0.0 to accept connections
                        on every interface, which is what a compute node needs
                        when a tunnel terminates on a different machine such
                        as a cluster login node. A non-loopback bind has no
                        authentication: every user who can reach the port
                        reads the dataset
  --allowed-host HOST   Extra exact Host name this server answers to; repeat
                        once per name. Needed only behind a reverse proxy such
                        as jupyter-server-proxy. One bare DNS name or IP
                        address, with no port, scheme or wildcard
  --no-browser          Don't open browser automatically
  --quiet, -q           Suppress info messages
  --verbose, -v         Enable verbose/debug logging
  --no-web-ui           Serve scientific data endpoints without establishing
                        the web application
  --web-source-url URL  Exact origin publishing cellucid-web-assets.json
  --web-cache-dir PATH  Directory for the active verified web build
  --latent-key KEY      Explicit obsm latent-space key for AnnData centroid
                        outliers
  --dataset-name NAME   Explicit dataset name (required for h5ad and zarr)
  --dataset-id ID       Explicit stable dataset identifier (required for h5ad
                        and zarr)
  --obs-key KEY         Exact AnnData obs column to serve; repeat once per
                        column. Every column is served when this is omitted.
                        Naming the columns leaves out one Cellucid cannot
                        serve, such as a datetime column
  --vector-field-default FIELD_ID
                        Exact default field id when AnnData contains multiple
                        UMAP vector fields
  --vector-fields       Serve the obsm '<field>_umap_<n>d' vector arrays as an
                        animated overlay. Off by default: every declared field
                        is read and checked before the server binds
  --connectivity        Serve the obsp['connectivities'] neighbor graph. Off
                        by default: reading it costs time and memory
                        proportional to the stored neighbor count, which on a
                        large dataset is the longest part of startup

Auto-detection:
    - .h5ad files → served directly via AnnData
    - .zarr directories → served directly via AnnData
    - Export directories (or exports roots with multiple datasets) → served as pre-exported data

Examples:
    # Serve an h5ad file (auto-detected)
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example

    # Serve a zarr store (auto-detected)
    cellucid serve /path/to/data.zarr --dataset-name Example --dataset-id example

    # Serve pre-exported data (auto-detected)
    cellucid serve /path/to/exported_dataset

    # Serve only named obs columns, leaving out one Cellucid cannot serve
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \
        --obs-key cell_type --obs-key total_counts

    # Serve on a different port
    cellucid serve /path/to/data --port 9000

    # Behind jupyter-server-proxy, name the proxy's host explicitly
    cellucid serve /path/to/data --allowed-host hub.example.org

    # Serve the neighbor graph and the velocity overlay, both off by default
    # because each is read and checked in full before the server binds
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \
        --connectivity --vector-fields

    # For SSH tunnel access from a remote server:
    # On the server: cellucid serve /path/to/data
    # On local machine: ssh -N -L 8765:127.0.0.1:8765 you@server
    # Then open the exact Viewer URL printed by Cellucid.

    # On an HPC compute node, where a tunnel can only terminate on the login node:
    # On the compute node: cellucid serve /path/to/data --host 0.0.0.0 --no-browser
    # On local machine:    ssh -N -L 8765:compute-node:8765 you@login-node
    # The login node opens the second connection, so the server has to accept it
    # on every interface. Then open the exact Viewer URL printed by Cellucid.
```

One thing varies with your machine: `argparse` wraps the option column to your
terminal width, so a narrower window rewraps the help text. The transcript above
was taken at 80 columns. The `Auto-detection` and `Examples` sections are printed
verbatim and never rewrap. The flags below are the ones worth understanding
before you run anything.

Key flags:

- `--port, -p`:
  - Change port if `8765` is in use. `--port 0` asks the operating system for
    an available port; copy the URL Cellucid prints.

- `--host, -H`:
  - Default is `127.0.0.1`, which no other machine can reach.
  - `0.0.0.0` accepts connections on every interface. That is what a compute
    node needs when the tunnel has to terminate on a different machine, such as
    a cluster login node.
  - A non-loopback bind has **no authentication**: every user who can reach the
    port reads the dataset.
  - A server bound to a wildcard never prints `http://0.0.0.0:<port>`, because
    no browser can open it. The compute-node recipe under **Remote Server
    Access** below shows what it prints instead.

- `--allowed-host HOST`:
  - Extra exact `Host` name this server answers to; repeat once per name.
  - Required only behind a reverse proxy such as `jupyter-server-proxy`, which
    forwards the browser’s `Host` unchanged. Without it those requests return
    `421 Misdirected Request` — the refusal that stops DNS rebinding.
  - One bare DNS name or IP address per entry: no port, no scheme, no wildcard.
    It matches whatever port the proxy publishes.

- `--no-browser`:
  - Don’t auto-open a browser tab.

- `--latent-key`:
  - Explicit AnnData `obsm` key used as latent space for categorical centroid
    outlier calculations.

- `--dataset-name` and `--dataset-id`:
  - Required for direct H5AD and Zarr input; invalid for prepared input.

- `--obs-key KEY`:
  - Exact AnnData `obs` column to serve; repeat once per column. Every column is
    served when you omit it.
  - Name the columns when one of them is a type Cellucid cannot serve — a
    datetime column is the usual case — so the rest still load.

- `--vector-fields`:
  - Also serve the `obsm` `<field>_umap_<n>d` arrays as the animated velocity
    overlay. **Off by default**, because every declared field is read and
    checked before the server binds.
  - Direct AnnData input only, exactly like `--connectivity`. Without it the
    overlay is simply absent, and step 3 of startup says so.

- `--vector-field-default`:
  - Required only when direct AnnData exposes multiple vector-field IDs and
    you need to choose the exact default.
  - It selects among *served* fields, so passing it without `--vector-fields`
    is an error rather than a no-op: there is nothing to choose between.

- `--connectivity`:
  - Also serve the `obsp['connectivities']` neighbor graph. **Off by default**,
    because reading it costs time and memory proportional to the stored
    neighbor count, and that is the longest part of startup on a large object.
  - Direct AnnData input only. Like every other AnnData-only flag, it is
    rejected on a prepared export, which already declares in
    `dataset_identity.json` whether it carries connectivity.
  - Passing it when `obsp['connectivities']` is absent is an error, not
    silence.

- `--quiet, -q` and `--verbose, -v`:
  - Mutually exclusive output levels.

- `--no-web-ui`:
  - Serve scientific endpoints without establishing or serving the web
    application. Pair it with `--no-browser`; this is an API-consumer mode, not
    a browser-viewing workflow.

- `--web-source-url` and `--web-cache-dir`:
  - Advanced controls for establishing the exact verified web build. Leave
    them at their defaults unless you operate an audited Cellucid web origin
    and cache.

## Option #9/#10/#11 — Python API Server Mode

You can start the same blocking servers from Python (useful for scripts).

- `serve(export_dir)` serves a pre-exported folder.
- `serve_anndata(data, dataset_name=..., dataset_id=...)` serves `.h5ad`,
  `.zarr`, or an in-memory `AnnData`.

```python
from cellucid import serve

serve(
    "/path/to/export_dir",
    port=8765,
    host="127.0.0.1",
    open_browser=True,
)
```

```python
from cellucid import serve_anndata

serve_anndata(
    "/path/to/data.h5ad",
    port=8765,
    host="127.0.0.1",
    open_browser=True,
    dataset_name="My dataset",
    dataset_id="my-dataset",
)
```

Both convenience calls block until you interrupt them. Use the class APIs when
another part of your program must keep control of the process.

The Python path takes the same two opt-in decisions as the CLI, through the
keyword-only `serve_connectivity` and `serve_vector_fields` parameters. Both
default to `False` on `serve_anndata()`, `show_anndata()`, `AnnDataViewer`,
`AnnDataAdapter`, `AnnDataAdapter.from_file()`, and `AnnDataServer`:

```python
from cellucid import serve_anndata

serve_anndata(
    "/path/to/data.h5ad",
    port=8765,
    host="127.0.0.1",
    open_browser=True,
    dataset_name="My dataset",
    dataset_id="my-dataset",
    serve_connectivity=True,
    serve_vector_fields=True,
)
```

`serve()` has neither parameter: a prepared export either holds the
connectivity and vector artifacts or does not.

### Stopping the server from Python

If you use the class-based API (advanced), you can stop the server programmatically:

```python
from cellucid.anndata_server import AnnDataServer

server = AnnDataServer(
    "data.h5ad",
    open_browser=False,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
server.start_background()

# ... interact in the browser ...

server.stop()
```

For most users, **Ctrl+C** in the terminal is simplest.

## Remote Server Access (SSH Tunnel Workflow)

This is the recommended way to use Cellucid when your data lives on a remote machine (HPC, lab server).

### Step 1 — Start the server on the remote machine

On the remote machine:

```bash
cellucid serve /path/to/data.h5ad \
  --no-browser \
  --dataset-name "My study" \
  --dataset-id my-study-v1
```

Keep it bound to `127.0.0.1` (default).

### Step 2 — Create an SSH tunnel from your laptop

On your laptop:

```bash
ssh -L 8765:localhost:8765 user@remote-host
```

Leave that SSH session open.

### Step 3 — Open Cellucid locally

```text
http://127.0.0.1:8765/?anndata=true
```

Use the exact path printed by the remote command: `/?anndata=true` for direct
H5AD/Zarr or `/?source=remote` for a prepared dataset or catalog.

Why this works well:
- Your browser talks only to `localhost`.
- You avoid exposing the server to the public internet.
- You avoid HTTPS→HTTP mixed-content blocking because the viewer UI and data API are served from the same origin.

### HPC: serving from a compute node

The three steps above assume your laptop can open an SSH connection straight to
the machine running the server. On a cluster it usually cannot: you reach a
login node, and the job — with the data and the memory — runs on a compute node
that only the login node can talk to. The tunnel therefore terminates on the
login node, and the server has to accept a connection arriving from a different
machine.

That is the one case where `--host 0.0.0.0` is the right answer rather than a
shortcut.

1) In the job script, on the compute node:

```bash
cellucid serve /scratch/you/atlas.h5ad \
  --dataset-name "My study" --dataset-id my-study-v1 \
  --host 0.0.0.0 --no-browser
```

2) From your laptop, forwarding through the login node to the compute node the
   scheduler gave you:

```bash
ssh -N -L 8765:compute-node-42:8765 you@login-node
```

3) Open the **Viewer URL** the job printed, with `127.0.0.1` as the host,
   because your end of the tunnel is local.

A wildcard bind changes the banner, and this is what you should expect to find
in the job's output file:

```text
════════════════════════════════════════════════════════════
  CELLUCID SERVER RUNNING
════════════════════════════════════════════════════════════

  Local URL:    http://127.0.0.1:8765
  Viewer URL:   http://127.0.0.1:8765/?anndata=true

  Bound to every network interface. From another machine:
  Machine URL:  http://compute-node-42:8765
  Viewer URL:   http://compute-node-42:8765/?anndata=true

  Press Ctrl+C to stop

════════════════════════════════════════════════════════════
```

Read that banner as two halves. The **Local URL** and **Viewer URL** at the top
are the loopback origin, which is the correct address for anything reaching the
server through a tunnel — including your browser at the laptop end of the
`ssh -L` above. `http://0.0.0.0:<port>` is never printed, because the wildcard
is a statement about which interfaces accept connections and not an address a
browser can open. The **Bound to every network interface** block appears only on
a wildcard bind and names the machine's own hostname, so it tells you the
`compute-node-42` half of the `ssh -L` argument without your having to look it
up. (Both hostnames above stand in for whatever your scheduler allocated, in the
same way the `/path/to/…` placeholders elsewhere on this page do.) An IPv6
literal bind is bracketed, so the URL stays a URL.

:::{warning}
A wildcard bind has no authentication. On a shared cluster, anyone who can
reach that port on that node can read the dataset. Bind to `0.0.0.0` only for
as long as the job runs, and prefer `127.0.0.1` whenever the tunnel can
terminate on the same machine as the server.
:::

For the full remote picture — jump hosts, `ProxyJump`, VS Code and JupyterHub
port forwarding, and cloud instances — see
{doc}`../../python_package/d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`.
For scheduler-specific job scripts and the compute-node hostname mechanics, see
{doc}`../../python_package/d_viewing_apis/17_hpc_slurm_and_compute_node_serving`.

## What a failed start looks like

A failed start ends with **one line beginning `Error:`**, and that line is the
whole message: it names what went wrong and what to change. There is no Python
traceback to wade through, because none of these is a Cellucid bug — they are
things you can fix from the command line. Read the last line.

### You forgot `--dataset-name` / `--dataset-id`

```text
$ cellucid serve lung_atlas.h5ad --no-browser
cellucid v0.9.1
Error: --dataset-name and --dataset-id are required when serving h5ad or zarr data. An .h5ad file and a Zarr store carry no Cellucid identity of their own, so name the dataset on the command line, for example: --dataset-name "My study" --dataset-id my-study-v1
```

Add both flags and run it again. A prepared export needs neither, because its
`dataset_identity.json` already carries them. Supplying only one is reported the
same way, naming just the flag still missing.

### The port is taken

```text
$ cellucid serve ./pancreas_export --no-browser
cellucid v0.9.1

[1/3] Validating dataset...
      ✓ Dataset valid

[2/3] Loading dataset info...
      ✓ Dataset loaded

[3/3] Starting server...
Error: Port 8765 is already in use on 127.0.0.1, so Cellucid could not start its server there. Another program is listening on that port, often an earlier Cellucid that was never stopped. Serve on a different port with --port 9000, or let the operating system pick a free one with --port 0 and use the Viewer URL it prints.
```

The dataset validated and loaded; only the socket failed. Nothing was wrong with
your data — another process holds the port.

### You asked for connectivity and the object has no graph

```text
$ cellucid serve pbmc_demo.h5ad \
      --dataset-name "PBMC demo" --dataset-id pbmc-demo --no-browser --connectivity
cellucid v0.9.1

Importing dependencies (anndata, numpy, scipy)... done

[1/5] Detecting format...
      Path: pbmc_demo.h5ad
      Format: h5ad
      ✓ Format detected

[2/5] Loading AnnData...
      Mode: read-only backed h5ad
      Obs columns: classifying 3
      Embeddings: resolving obsm keys
      Embeddings: 2D from obsm['X_umap']
Error: Connectivity was asked for, and adata.obsp has no 'connectivities' matrix to serve. Compute the neighbor graph with sc.pp.neighbors(adata) before serving, or serve without asking for connectivity.
```

Asking for a graph that is not there is an error rather than silence. Either
compute the neighbors first or drop `--connectivity`; every other capability of
the dataset serves fine without it.

### The bare `X_umap` is not 1, 2, or 3 columns wide

```text
$ cellucid serve wide_umap.h5ad \
      --dataset-name "PBMC demo" --dataset-id pbmc-demo --no-browser
cellucid v0.9.1

Importing dependencies (anndata, numpy, scipy)... done

[1/5] Detecting format...
      Path: wide_umap.h5ad
      Format: h5ad
      ✓ Format detected

[2/5] Loading AnnData...
      Mode: read-only backed h5ad
      Obs columns: classifying 3
      Embeddings: resolving obsm keys
Error: No UMAP embedding Cellucid can read was found in adata.obsm. Cellucid reads the exact keys 'X_umap_1d', 'X_umap_2d', and 'X_umap_3d', and reads a bare 'X_umap' as the dimension its own column count states. Available obsm keys: ['X_umap']. adata.obsm['X_umap'] has shape (4096, 10), and Cellucid draws 1, 2, or 3 dimensions, so that array names a dimension no viewer renders. Assign the columns you mean to draw under the key for their own count.
```

The message names the shape it found, which is the fastest way to see what
happened: a ten-column array under `X_umap` is a latent space named after a
plot. Write the columns you actually want to draw to `X_umap_2d` (or `_1d` /
`_3d`) and serve again.

### A missing package, an old package, and a Python that is too old

`cellucid serve` distinguishes three import conditions instead of advising
`pip install` for all of them:

| Condition | What the message says |
| --- | --- |
| The distribution is not installed in this environment | Names the package and gives `pip install <name>` |
| The module is installed, but too old to hold the name Cellucid asked of it | Says the installed version is not one this Cellucid can use, and gives `pip install --upgrade <name>` |
| The module belongs to the standard library and this interpreter cannot import it | Names the Python version, the interpreter's own path, and the Python range Cellucid requires — and offers no `pip` command, because no package installs a standard-library module |

The third case is the one that used to send people in circles: no amount of
installing fixes a module the running interpreter is supposed to provide, so
the message points at the interpreter instead.

:::{tip}
If you ever do see a Python traceback, that is Cellucid's own failure, not
yours. The final `Error:` line says so and gives the address to report it to;
copy the traceback into the report.
:::

## Edge Cases

- **Port already in use**:
  - Startup gets as far as binding the socket and then fails with
    `Error: Port 8765 is already in use on 127.0.0.1, …`.
  - Choose another port explicitly or use `--port 0`, which asks the operating
    system for a free one; the banner then prints the port it actually got, for
    example `http://127.0.0.1:64088/?source=remote`. Copy that URL.

- **`--host` names an address this machine does not have**:
  - The bind fails with `Error: Host '203.0.113.9' is not an address this
    machine can serve from. …`
  - Use `--host 127.0.0.1` (the default) for this computer only, or
    `--host 0.0.0.0` to accept connections from other machines.

- **You bound to a wildcard and want the address to type**:
  - `--host 0.0.0.0`, `--host ::`, and an empty host all bind every interface,
    and none of them is an address a browser can open.
  - The banner prints the loopback origin as **Local URL** / **Viewer URL** —
    correct for the tunnelled or same-machine case — and adds a
    `Bound to every network interface. From another machine:` block with this
    machine's own hostname URLs when it has one.

- **The graph is in the file but no edges are drawn**:
  - Step 3 says
    `Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)`.
  - Restart with `--connectivity` (CLI) or `serve_connectivity=True` (Python).
    Connectivity is off by default on the direct-AnnData path because building
    the edge list is the longest part of startup on a large object.
  - A prepared export is unaffected: it either ships the connectivity artifacts
    or it does not.

- **Windows firewall prompt**:
  - If you allow public network access accidentally, other machines may reach your server.

- **Large in-memory AnnData or Zarr input**:
  - Can cause high RAM use.

- **Mixed content**:
  - Pointing an HTTPS Cellucid page at an `http://` server is refused by
    Cellucid itself, before the browser gets a chance to block the fetch:
    `An HTTPS Cellucid page requires an explicit HTTPS remote server URL`.
  - Always open the exact local Viewer URL printed by the server
    (`/?source=remote` for prepared data or `/?anndata=true` for direct
    AnnData), which serves the verified UI generation.
  - Prefer an SSH tunnel when the kernel/server is remote.

- **Vector fields only exist in one dimension**:
  - If you have `*_umap_2d` vectors but no `*_umap_3d`, the overlay dropdown will be empty in 3D.
  - Switch to 2D, or provide the missing dimension.

## Troubleshooting (Massive)

### Symptom: “Port already in use”

**Likely causes**
- Another service is using that port.

**How to confirm**
- The final line reads `Error: Port 8765 is already in use on 127.0.0.1, …`.
- Try a different port: `--port 9000`.

**Fix**
- Run:

  ```bash
  cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset --port 9000
  ```

- Then copy the exact printed Viewer URL. For this command it should be:

  ```text
  http://127.0.0.1:9000/?anndata=true
  ```

---

### Symptom: “Cellucid says it can’t connect to the remote server”

**Likely causes**
1) You typed the wrong URL (port mismatch).
2) You pasted the **Viewer URL** into the `Remote server:` field. That field
   takes the bare origin (`http://127.0.0.1:8765`) and refuses a trailing
   slash, a query string, or a fragment.
3) The server is bound to `127.0.0.1` but you are trying to access it from another machine.
4) The Cellucid page is on HTTPS and the server address is `http://`, which
   Cellucid refuses outright.

**How to confirm**
- Open the server URL directly in your browser:

  ```text
  http://127.0.0.1:8765/_cellucid/health
  ```

  (If the server is running locally, you should get a small JSON response.)
- For an exports root, inspect the validated catalog separately at
  `http://127.0.0.1:8765/_cellucid/datasets`.

**Fix**
- Use an SSH tunnel for remote machines.
- Ensure the server is actually running and the URL matches.

---

### Symptom: “It connects, but genes are missing / gene search returns nothing”

**Likely causes**
- Your AnnData has no expression matrix (`adata.X` empty) or no `var`.
- The names you are searching for are in a different `var` column from the one the genes were named with.

**How to confirm**
- In Python: `print(adata.X is None)`, `print(adata.var.head())`.

**Fix**
- Use `show_anndata(..., gene_id_column="...")` in Jupyter.
- Or export via `prepare(var_gene_id_column="...")`.

---

### Symptom: “Vector field overlay toggle is disabled / no fields appear”

**Likely causes (ordered)**
1) The server was started without `--vector-fields` / `serve_vector_fields=True`,
   which is the default on the direct-AnnData path.
2) The dataset contains no vectors (common).
3) Vector fields exist, but are not named using the expected convention (`*_umap_2d`, `*_umap_3d`).
4) Dimension mismatch: vectors exist for 2D but you’re viewing 3D (or vice versa).

**How to confirm**
- Read step 3 of the startup output. `not served (…)` is cause 1 and `no` is
  cause 2 or 3; the wording is quoted in full under **Vector fields are opt-in
  on this path too** above.
- Open `http://127.0.0.1:8765/dataset_identity.json` and check for a `vector_fields` block.
- If using AnnData, list `obsm` keys and look for `velocity_umap_2d`-style entries.

**Fix**
- Cause 1: restart with `--vector-fields`, or `serve_vector_fields=True` from
  Python.
- Rename or regenerate vector fields using the exact keys documented in
  {doc}`../i_vector_field_velocity/01_what_vector_fields_are_user_facing`.
- Switch the viewer to the dimension that has vectors.

For overlay UI behavior and deeper debugging, see:
- {doc}`../i_vector_field_velocity/index`
- {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

---

### Symptom: “The connectivity/edge overlay has nothing to draw”

**Likely causes (ordered)**
1) The server was started without `--connectivity` / `serve_connectivity=True`,
   which is the default on the direct-AnnData path.
2) The object genuinely has no `obsp['connectivities']`.
3) You are serving a prepared export that was written without a
   `connectivities=` argument to `prepare()`.

**How to confirm**
- Read step 3 of the startup output. `not served (…)` is cause 1 and `no` is
  cause 2; the wording is quoted in full under
  **Connectivity is opt-in on this path** above.
- Open `http://127.0.0.1:8765/connectivity_manifest.json`. A 404 means the
  edges are not being served, and so are
  `connectivity/edges.src.bin`, `connectivity/edges.dst.bin`, and
  `connectivity/edges.weights.f64.bin`.

**Fix**
- Cause 1: restart with `--connectivity`, or `serve_connectivity=True` from
  Python. The graph is read and validated before the socket opens, so expect
  startup to take longer.
- Cause 2: compute the graph first, for example with `sc.pp.neighbors(adata)`.
- Cause 3: re-export with `prepare(connectivities=...)`.

---

### Symptom: “It’s extremely slow”

**Likely causes**
- Large dataset + dense X
- Running over a high-latency SSH tunnel
- Connectivity was asked for on a very large graph

**Fix**
- Prefer pre-exported data for best performance.
- If using SSH, run the server close to the data (on the same machine as the file) and tunnel to it.
- Leave `--connectivity` off unless you are going to draw the graph. A
  50-neighbor graph over millions of cells is hundreds of millions of stored
  neighbors, and every one of them is read and validated before the server
  binds. Past 50,000,000 stored neighbors, the validator logs a warning naming
  the cost before it starts, so a long silence there is expected rather than a
  hang.

## Next Steps

- Dataset identity (sessions/sharing): {doc}`06_dataset_identity_why_it_matters`
- Full troubleshooting matrix: {doc}`08_troubleshooting_data_loading`
- Want notebook embedding + programmatic control? → {doc}`05_jupyter_tutorial`
- Want browser-only loading without any server? → {doc}`03_browser_file_picker_tutorial`
- Want to publish a static public catalog instead? →
  {doc}`11_custom_dataset_repository`

Will it be usable at your scale? These two answer that, and are worth reading
before you hand the URL to a lab:

- What the server keeps in memory, what it reads lazily, and what backed mode
  costs → {doc}`../../python_package/d_viewing_apis/14_performance_scaling_and_lazy_loading`
- What to turn on and off in the viewer at tens of millions of cells →
  {doc}`../n_benchmarking_performance/03_large_dataset_best_practices`
