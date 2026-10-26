# Server Mode (CLI + Python)

Server mode hands an `.h5ad` file or a Zarr store — or any other file or
folder Cellucid reads (MuData, Loom, 10x, MatrixMarket, an R object) — to the
viewer, which reads it itself: the browser's own HDF5/Zarr reader fetches only
the bytes a view needs, with HTTP range requests. Nothing is converted first.

You run a small local Python server that:
- serves the store's bytes exactly as they are on disk (an `.h5ad` file, a Zarr
  store folder, a zipped Zarr store, or any other format of
  {doc}`07_file_format_expectations`), and computes nothing per request
- names the dataset and states the launch choices you made on the command line
- establishes and serves the exact viewer generation on the same origin

Then you open the exact **Viewer URL printed by the server** (`/?anndata=true`).

This tutorial covers the routes where a Cellucid server serves the data: the
terminal and the Python API (the full list: {doc}`01_loading_options_overview`).

## At A Glance

**Audience**
- Wet lab / beginner: follow the copy/paste commands; you do not need to write code.
- Computational users: focus on what is read when, launch choices, and SSH tunnel workflows.

**Time**
- Local server (same machine): ~5–10 minutes
- Remote server (SSH tunnel): ~15–30 minutes

**Prerequisites**
- `pip install cellucid`
- Your data in one of these forms:
  - `.h5ad` file
  - Zarr store folder (Zarr v2 with `.zgroup`, or v3 with `zarr.json`)
  - zipped Zarr store (`.zip`, for example `data.zarr.zip`)
  - or any other format Cellucid reads: `.h5mu` or MuData Zarr, `.loom`, a 10x
    `.h5` or `outs/` folder, a MatrixMarket folder, an `.rds` / `.rda` /
    `.RData` (gzipped or not)

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

   An AnnData that already declares its Cellucid identity in
   `uns['cellucid']['name']` and `uns['cellucid']['id']` needs neither flag.

2) Copy the exact **Viewer URL** printed by the command. With the default port
   it is:

```text
http://127.0.0.1:8765/?anndata=true
```

3) Keep the terminal running while you use the viewer.
4) Stop the server with **Ctrl+C**.

The viewer does not watch the server: if the server stops while a viewer is
open, what is already shown stays, and the next read it needs (a gene, a
field, the neighbour graph) is asked again a few times and then fails with a
notification. Start the server again and open its address again to continue.

### What success looks like in the browser

```{figure} ../../../_static/screenshots/server/open-served-dataset.png
:alt: Cellucid opened at the Viewer URL of cellucid serve, coloured by clusters: the Dataset box names the Pancreas sample with its cells and genes, "Zarr v2 folder" and the server's loopback address, and how much was read; the Link or server bar holds that address with one line under it reading Server, the address and Disconnect.
:width: 1440px

The same server, opened at its printed Viewer URL. Two things tell you the
connection is live rather than a leftover sample: the **Link or server** bar
holds the server's address, and a line under it reads **Server**, the address
and **Disconnect**. The **Dataset** box above describes the server's dataset,
with the format and the server as its host (here the Pancreas sample's Zarr
folder: `Zarr v2 folder · 127.0.0.1:<port>`; an `.h5ad` reads `h5ad file`),
not the sample catalog.
```

:::{important}
**The banner prints two URLs, and they are not interchangeable.**

- **Viewer URL** (`http://127.0.0.1:8765/?anndata=true`) goes in the
  **browser address bar**. Opening it is the
  whole workflow; the **Link or server** bar fills itself in.
- **Link or server** — the address bar under **Open data** in the Session
  panel — takes the server's address: the port alone (`8765`),
  `localhost:8765` or `127.0.0.1:8765` (its placeholder reads
  `Link or localhost:8765`), or the bare origin the banner prints as
  **Local URL**, `http://127.0.0.1:8765`. The Viewer URL carries a `?…` query,
  and the bar opens an address with a query as a link to data, not as a
  server, so paste the Local URL instead. See
  [The address bar](14_open_from_link.md#the-address-bar) for every rule.

You only need the bar when you are pointing an already-open Cellucid tab at a
server, for example after a tunnel comes up. Otherwise open the Viewer URL and
skip it.
:::

## Serve `.h5ad` or Zarr (Auto-detected)

```bash
# h5ad
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset

# Zarr store folder
cellucid serve /path/to/data.zarr --dataset-name "My dataset" --dataset-id my-dataset

# zipped Zarr store
cellucid serve /path/to/data.zarr.zip --dataset-name "My dataset" --dataset-id my-dataset
```

What happens:
- The file or folder is served as it is on disk: never copied, never converted,
  never loaded into the server's memory.
- The viewer reads it itself with HTTP range requests, as each view needs it:
  the obs columns and embeddings when the dataset opens, a gene when you color
  by it, the neighbour graph when a view draws it, a vector field when the
  overlay shows it.
- The store is published under a route named by the dataset id
  (`/_cellucid/anndata/<dataset_id>.h5ad`, `.zarr.zip`, or `.zarr/<key>`), so
  your local file name never leaves the machine.

The server runs four numbered steps, then prints its URLs:

```text
$ cellucid serve pbmc_demo.h5ad \
      --dataset-name "PBMC demo" --dataset-id pbmc-demo --no-browser --port 8766
cellucid v0.9.5

[1/4] Detecting format...
      Path: /path/to/pbmc_demo.h5ad
      Format: AnnData (.h5ad)
      ✓ Format detected

[2/4] Reading metadata...
      Cells: 4,096
      Genes: 512
      obsm: X_umap
      obs fields: 3
      obsp: connectivities
      layers: none
      ✓ Metadata read

[3/4] Checking launch choices...
      ✓ Launch choices checked

[4/4] Starting server...
      Viewer UI generation: establishing exact configured source
Cellucid web UI cache: 100%|██████████| 590/590 [02:31<00:00,  3.90file/s]
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
three `obs` columns, one `obsm['X_umap']`, and one `obsp['connectivities']`.
Steps 1 to 3 read metadata only — names, shapes, and the codecs arrays are
stored with — never an array, so an atlas prints the same steps in about the
same time.

Two lines vary with your machine rather than with your data. The `Path:` line
prints the **absolute** path Cellucid resolved, so it shows your own location,
not the `/path/to/…` placeholder above. The progress bar is sized to your
terminal width and paced by your network; its file count is the size of the
published web generation on the day you run it, so yours will be whatever is
current.

### What each of the four steps tells you

| Step | What it does | What to read it for |
| --- | --- | --- |
| `[1/4] Detecting format` | Classifies the path by the input format registry (`AnnData (.h5ad)`, `AnnData Zarr, Zarr v2 folder`, `MuData (.h5mu)`, `Loom (.loom)`, …) and reads the dataset's identity | The `Format:` line confirms Cellucid agreed with you about what the path is |
| `[2/4] Reading metadata` | Reads the cell and gene counts and the names of the `obsm`, `obs`, `obsp`, and `layers` entries | The sanity check worth reading before you open the browser |
| `[3/4] Checking launch choices` | Checks every launch option you passed against that metadata, and prints each one it will hand the viewer | A misspelled key fails here, on the terminal, before a socket opens |
| `[4/4] Starting server` | Establishes the exact viewer generation and binds the socket | The banner, and the **Viewer URL** you copy |

Two more lines appear only when they apply:

- A Zarr store folder without consolidated metadata prints
  `Consolidated metadata: built in memory (.zmetadata)` (Zarr v2) or
  `(zarr.json)` (Zarr v3) in step 3. A browser cannot list a folder, so the
  server builds the document once from the store's own metadata files and
  answers it in place of the missing one. Nothing is written to your store.
- A store holding an array the browser cannot decode — an HDF5 filter other
  than gzip, shuffle, Fletcher-32, or LZF (for example Blosc, Zstandard, szip,
  or scaleoffset), or a Zarr codec such as bz2, lzma, or zfpy — prints
  `Browser cannot decode:` with the arrays and their codecs, then
  `Serving: a decodable copy, written once to a temporary .h5ad`. The store is
  read once with anndata (and its codec plugins), and the browser reads the
  copy. The copy is removed when the server stops.

### Launch choices: what the viewer opens on

Every AnnData option is optional. Each one you pass is checked against the
store's metadata in step 3 — the named `obsm`, `obsp`, layer, or column must
exist — and handed to the viewer as a launch default. An option you leave out
is taken from the store's own `uns['cellucid']['defaults']` when it declares
one, else the viewer applies its default rule. The rules, and how to change any
of them while the dataset is open, are on {doc}`13_data_mapping`.

| CLI flag | Python keyword | Chooses |
| --- | --- | --- |
| `--embedding-1d`, `--embedding-2d`, `--embedding-3d KEY[:COLUMNS]` | `embedding_keys={"2d": "X_tsne"}` or `{"2d": ("X_pca", [0, 1])}` | The `obsm` array (and 0-based columns) of a dimension slot |
| `--knn KEY` | `knn_key` | The `obsp` matrix drawn as the neighbour graph |
| `--expression X\|raw\|layer:NAME` | `expression_source` | Where gene expression is read |
| `--gene-id-column COLUMN` | `gene_id_column` | The `var` column (`raw.var` with `raw`) whose values name the genes |
| `--obs-key KEY` (repeat) | `obs_keys` | The `obs` columns shown, in this order |
| `--modality NAME` | `modality` | The modality a MuData dataset opens on ({doc}`19_multimodal_mudata`) |
| `--latent-key KEY` | `latent_key` | The `obsm` latent space (scVI, PCA) the outlier filter measures in |
| `--vector-field-default FIELD_ID` | `vector_field_default` | The vector field the overlay opens on |
| `--normalize-embeddings`, `--no-normalize-embeddings` | `normalize_embeddings` | Whether embeddings are scaled into the viewer's coordinate range (the default) or shown at the coordinates the store holds |
| `--centroid-min-points N` | `centroid_min_points` | The fewest cells a category needs for the outlier filter to measure it (10 when left out) |

Step 3 prints each choice it accepted, one per line:

```text
$ cellucid serve pbmc_demo.h5ad --dataset-name "PBMC demo" --dataset-id pbmc-demo \
      --no-browser --knn connectivities --embedding-2d X_umap:1,0
...
[3/4] Checking launch choices...
      embedding2d: {'key': 'X_umap', 'columns': [1, 0]}
      knn: connectivities
      ✓ Launch choices checked
```

Which `obsm` keys become embeddings:

- `--embedding-1d`, `--embedding-2d`, and `--embedding-3d` name the `obsm` array
  of a slot, optionally with 0-based columns: `--embedding-2d X_pca:0,1`.
- Every slot you do not name follows the default rule: `X_umap_1d`, `X_umap_2d`,
  and `X_umap_3d` each fill their own slot, and when an object declares none of
  those, a bare `X_umap` — the key `sc.tl.umap()` writes — fills the slot of its
  own column count (1, 2, or 3). No key is renamed and your object is not
  modified.
- Only when the default rule fills no slot does the viewer make a deterministic
  best guess from the other `obsm` arrays — UMAP-like names first, latent
  spaces such as PCA last — and report it as guessed. A ten-column array with no
  narrower candidate is shown in 2D at columns 0 and 1.

No launch choice costs startup time. `--knn connectivities` names the graph
the viewer draws; the viewer reads it the first time a view draws edges, and
`--vector-field-default` names the field the overlay opens on, read when the
overlay shows it.

A store written by `cellucid.prepare()` carries its own launch choices in
`uns['cellucid']['defaults']`; a choice you name on the command line overrides
the store's.

### What the server publishes

The server answers:

| Route | Answer |
| --- | --- |
| `/_cellucid/datasets` | The dataset, where its store is, and the launch choices you set |
| `/_cellucid/anndata/<dataset_id>.h5ad` (or `.zarr.zip`, `.h5mu`, `.loom`, `.h5`, `.rds.gz`, … — the file's own ending) | The file's bytes, with `Range`, `206 Partial Content`, `Content-Range`, `HEAD`, and `ETag` answered exactly |
| `/_cellucid/anndata/<dataset_id>.zarr/<key>` | One key of a Zarr store folder; an absent key is a 404, which a Zarr reader takes as "absent" |
| `/_cellucid/anndata/<dataset_id>/<path>` | One file of a MatrixMarket or Cell Ranger `outs/` folder, which `/_cellucid/datasets` lists with its size |
| `/_cellucid/health` | `{"status": "ok", "type": "anndata", "version", "format", "n_cells", "n_genes"}`; the counts are `null` for a format Python does not read (an R object, Loom, 10x, MatrixMarket) |
| `/_cellucid/info` | The same facts plus the bound `host` and `port` |
| `/_cellucid/protocol` | Where the notebook viewer channel is, and the commands and events it carries |
| `/_cellucid/ws` | The notebook viewer channel, a WebSocket (see {doc}`05_jupyter_tutorial`) |

For the demo above, `/_cellucid/datasets` is:

```json
{"datasets": [{"id": "pbmc-demo", "name": "PBMC demo", "anndata": {"file": {"url": "_cellucid/anndata/pbmc-demo.h5ad", "kind": "h5ad", "layout": "file", "size": 2799976}, "defaults": {"embedding2d": {"key": "X_umap", "columns": [1, 0]}, "knn": "connectivities"}}}]}
```

`url` is relative to the server's base URL, because the server cannot know the
address a browser reaches it by (a tunnel, a proxy). `size` is given for a file
and absent for a folder. `defaults` holds only the choices you set. The server
is threaded and keeps connections alive (HTTP/1.1), because a browser reads one
file with many range requests.

```{note}
The viewer reads a served store with the same reader it uses for an `.h5ad` or
Zarr store you open from your disk or from a link, so what it accepts, what it
reports, and every data-mapping control behave the same way on all three
routes.
```

## CLI Options (What They Mean)

Run `cellucid serve --help` anytime. Its complete output:

```text
$ cellucid serve --help
usage: cellucid serve [-h] [--port PORT] [--host HOST] [--allowed-host HOST]
                      [--no-browser] [--quiet | --verbose] [--no-web-ui]
                      [--web-source-url URL] [--web-cache-dir PATH]
                      [--latent-key KEY] [--centroid-min-points N]
                      [--dataset-name NAME] [--dataset-id ID] [--obs-key KEY]
                      [--modality NAME] [--vector-field-default FIELD_ID]
                      [--embedding-1d KEY[:COLUMNS]]
                      [--embedding-2d KEY[:COLUMNS]]
                      [--embedding-3d KEY[:COLUMNS]]
                      [--normalize-embeddings | --no-normalize-embeddings]
                      [--knn KEY] [--gene-id-column COLUMN]
                      [--expression X|raw|layer:NAME]
                      data_path

Serve a dataset file or folder to the viewer, which reads it itself.

positional arguments:
  data_path             Path to a file or folder Cellucid reads: .h5ad,
                        AnnData or MuData Zarr (folder or .zip), .h5mu, .loom,
                        a 10x .h5 or outs/ folder, a MatrixMarket folder, or
                        an R object (.rds, .rda, .RData, gzipped or not)

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
  --latent-key KEY      obsm latent space (scVI, PCA) the outlier filter
                        measures in
  --centroid-min-points N
                        Fewest cells a category needs for the outlier filter
                        to measure it (default: the store's own, else 10)
  --dataset-name NAME   Dataset name (default: the store's
                        uns['cellucid']['name'], required without it)
  --dataset-id ID       Stable dataset identifier (default: the store's
                        uns['cellucid']['id'], required without it)
  --obs-key KEY         Exact AnnData obs column the viewer shows; repeat once
                        per column, in order. Every column is shown when this
                        is omitted
  --modality NAME       Modality a multimodal dataset (a MuData) opens on; the
                        default one when omitted
  --vector-field-default FIELD_ID
                        Vector field the overlay opens on: <field> of an obsm
                        '<field>_<basis>' array drawn over 'X_<basis>'
  --embedding-1d KEY[:COLUMNS]
                        AnnData obsm array to plot in 1D, optionally with its
                        0-based columns, e.g. X_pca:0. Overrides the default
                        (X_umap_<n>d or X_umap) for this dimension
  --embedding-2d KEY[:COLUMNS]
                        AnnData obsm array to plot in 2D, optionally with its
                        0-based columns, e.g. X_pca:0,1. Overrides the default
                        (X_umap_<n>d or X_umap) for this dimension
  --embedding-3d KEY[:COLUMNS]
                        AnnData obsm array to plot in 3D, optionally with its
                        0-based columns, e.g. X_pca:0,1,2. Overrides the
                        default (X_umap_<n>d or X_umap) for this dimension
  --normalize-embeddings, --no-normalize-embeddings
                        Scale every embedding into the viewer's coordinate
                        range; the no- form shows the coordinates as stored
                        (default: the store's own, else scaled)
  --knn KEY             AnnData obsp matrix the viewer draws as the neighbor
                        graph, usually 'connectivities'; read only when a view
                        draws it
  --gene-id-column COLUMN
                        AnnData var column (raw.var with --expression raw)
                        whose values name the genes; the var index is used
                        when omitted
  --expression X|raw|layer:NAME
                        Where AnnData gene expression is read: X, raw, or
                        layer:NAME (default: X when present, else raw or a
                        layer, reported as guessed)

Formats (detected, by the input format registry the web app reads):
    - .h5ad files; AnnData Zarr store folders (v2 or v3) and zipped stores (.zip)
    - .h5mu files and MuData Zarr stores
    - .loom files; 10x Cell Ranger .h5 files and outs/ folders
    - MatrixMarket folders (matrix.mtx.gz, or GEO's GSM..._matrix.mtx.gz naming)
    - R objects: .rds, .rda, .RData, gzipped or not (.rds.gz from GEO)
    The browser reads every one itself; the server serves the bytes. A store
    written by cellucid.prepare() names itself (uns['cellucid']), so
    --dataset-name and --dataset-id are needed only for other data.

Examples:
    # Serve an h5ad file (auto-detected)
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example

    # Serve a zarr store (auto-detected)
    cellucid serve /path/to/data.zarr --dataset-name Example --dataset-id example

    # Serve a Cellucid-ready AnnData written by cellucid.prepare()
    cellucid serve /path/to/my_dataset.zarr

    # Serve only named obs columns, leaving out one Cellucid cannot serve
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \
        --obs-key cell_type --obs-key total_counts

    # Serve on a different port
    cellucid serve /path/to/data --port 9000

    # Behind jupyter-server-proxy, name the proxy's host explicitly
    cellucid serve /path/to/data --allowed-host hub.example.org

    # Open on a chosen neighbor graph and velocity field
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \
        --knn connectivities --vector-field-default velocity

    # Choose what is plotted: t-SNE in 2D, PCA components 1 and 2 in a second
    # run, and expression from a layer named by the gene-symbol column
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \
        --embedding-2d X_tsne --expression layer:log1p --gene-id-column gene_symbols
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \
        --embedding-2d X_pca:0,1

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
was taken at 80 columns. The `Formats` and `Examples` sections are printed
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
  - The `obsm` latent space (scVI, PCA) the outlier filter measures in, never
    a 2D/3D layout. Without it the store's own `latentKey` is used, else none:
    the outlier filter then waits for one chosen in **Data mapping**.
  - The viewer reads the outlier quantiles the store holds for that latent
    space (`prepare(latent_key=...)` stores them), or computes them in the
    browser on the slider's first move, reading the latent space once; see
    {doc}`../e_filtering/02_outlier_filtering_per_active_field`.

- `--dataset-name` and `--dataset-id`:
  - Name the dataset, whatever its format. Each defaults to the
    store's own `uns['cellucid']['name']` / `uns['cellucid']['id']`, and is
    required when the store declares none.
  - Sessions and community annotations are keyed by the dataset id, so keep it
    stable across runs.

- `--obs-key KEY`:
  - Exact AnnData `obs` column the viewer shows; repeat once per column, in
    order. Every column is shown when you omit it.

- `--vector-field-default FIELD_ID`:
  - The vector field the velocity overlay opens on: the `<field>` of an
    `obsm['<field>_<basis>']` array drawn over `obsm['X_<basis>']` (`velocity`
    for `velocity_umap`, scVelo's name for velocity on `X_umap`). Naming a field
    the store does not declare is an error.
  - The viewer reads a vector field only when the overlay shows it.

- `--knn KEY`:
  - The `obsp` matrix the viewer draws as the neighbour graph, usually
    `connectivities`. The viewer reads it the first time a view draws edges, so
    it costs nothing at startup.
  - Passing a key the store does not have is an error, not silence.

- `--centroid-min-points N`:
  - The fewest cells a category needs for the outlier filter to measure it;
    smaller categories are left unfiltered. Without it the store's own
    `centroidMinPoints` is used, else 10.
  - Stored outlier quantiles are used only when they were written for the same
    number (`prepare(centroid_min_points=...)`); otherwise the browser computes
    them on the slider's first move.

- `--normalize-embeddings` / `--no-normalize-embeddings`:
  - Whether the viewer scales every embedding into its coordinate range (the
    default) or shows the coordinates the store holds, so figure axes,
    centroids and vectors keep the file's own units. Without either, the
    store's own `normalizeEmbeddings` is used, else scaling.

- `--embedding-1d`, `--embedding-2d`, `--embedding-3d`, `--expression`,
  `--gene-id-column`:
  - The launch choices in the table under **Launch choices** above.

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

## Python API Server Mode

You can start the same blocking server from Python (useful for scripts).

`serve_anndata(data, dataset_name=..., dataset_id=...)` serves an `.h5ad`
file, a Zarr store folder, a zipped Zarr store, any other file or folder
Cellucid reads, or an in-memory `AnnData` or `MuData`. An in-memory object is
written once to a temporary `.h5ad` or `.h5mu` (sparse matrices column-major,
uncompressed, so one gene is one contiguous read) that is removed when the
server stops or Python exits.

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

The convenience call blocks until you interrupt it. Use the class APIs when
another part of your program must keep control of the process.

The Python path takes the same launch choices as the CLI, as keyword-only
parameters that all default to `None` on `serve_anndata()`, `show_anndata()`,
`AnnDataViewer`, and `AnnDataServer`:

```python
from cellucid import serve_anndata

serve_anndata(
    "/path/to/data.h5ad",
    port=8765,
    host="127.0.0.1",
    open_browser=True,
    dataset_name="My dataset",
    dataset_id="my-dataset",
    knn_key="connectivities",
    vector_field_default="velocity",
)
```

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

Use the exact URL printed by the remote command.

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

### The store declares no identity and you named none

```text
$ cellucid serve lung_atlas.h5ad --no-browser
cellucid v0.9.5

[1/4] Detecting format...
      Path: /path/to/lung_atlas.h5ad
      Format: AnnData (.h5ad)
Error: dataset_name and dataset_id are required: this AnnData declares no Cellucid identity of its own (uns['cellucid']['id'] and ['name']). Name the dataset, for example dataset_name="My study", dataset_id="my-study-v1" (--dataset-name/--dataset-id on the command line).
```

Add `--dataset-name` and `--dataset-id` and run it again, or write the identity
into the object once (`adata.uns["cellucid"] = {"name": "My study", "id":
"my-study-v1"}`) so every later run picks it up. A file written by
`cellucid.prepare()` needs neither, because it carries both. When the store
declares one of the two, only the other is reported missing.

### The port is taken

```text
$ cellucid serve pbmc_demo.h5ad --dataset-name "PBMC demo" --dataset-id pbmc-demo --no-browser
…
[4/4] Starting server...
Error: Port 8765 is already in use on 127.0.0.1, so Cellucid could not start its server there. Another program is listening on that port, often an earlier Cellucid that was never stopped. Serve on a different port with --port 9000, or let the operating system pick a free one with --port 0 and use the Viewer URL it prints.
```

The file was read and checked; only the socket failed. Nothing was wrong with
your data — another process holds the port.

### A launch choice names something the store does not have

```text
$ cellucid serve pbmc_demo.h5ad \
      --dataset-name "PBMC demo" --dataset-id pbmc-demo --no-browser --knn connectivities
cellucid v0.9.5

[1/4] Detecting format...
      Path: /path/to/pbmc_demo.h5ad
      Format: AnnData (.h5ad)
      ✓ Format detected

[2/4] Reading metadata...
      Cells: 4,096
      Genes: 512
      obsm: X_umap
      obs fields: 3
      obsp: none
      layers: none
      ✓ Metadata read

[3/4] Checking launch choices...
Error: knn_key names obsp['connectivities'], which this store does not have. obsp holds: []. Compute the neighbour graph with sc.pp.neighbors(adata), or leave knn_key unset.
```

Every launch choice is checked the same way before the socket opens, and the
message lists what the store does hold:

```text
Error: embedding_keys['2d'] names obsm['X_tsne'], which this store does not have. obsm holds: ['X_umap'].
Error: expression_source 'layer:counts' names a layer this store does not have. layers holds: [].
Error: vector_field_default names field 'velocity', which this store does not declare. Its vector fields (obsm '<field>_<basis>' as wide as 'X_<basis>'): [].
```

Fix the name, or drop the flag and let the viewer's default rule choose.

### Nothing in `obsm` can be plotted

The server checks names, not values, so a store with no usable embedding still
starts. The viewer then reports it when it opens the store, listing every
`obsm` entry and why it was not used:

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell (4,096) and at least one numeric column. obsm holds: X_df (table), spatial (4,096 × 2, NaN at cell 17). Add a 2-D embedding to obsm (for example with sc.tl.umap in Scanpy), then open the data again.
```

Here a table and an array with a `NaN` coordinate. Repair the array, compute an
embedding, or name a usable one with `--embedding-2d KEY`, and serve again.

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
    example `http://127.0.0.1:64088/?anndata=true`. Copy that URL.

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
  - The viewer draws `obsp['connectivities']` by default. A graph stored under
    another key (for example `scvi_connectivities`) is drawn when you name it
    with `--knn KEY` (CLI) or `knn_key="KEY"` (Python), or pick it under
    **Session › Data mapping › Neighbour graph**.

- **Windows firewall prompt**:
  - If you allow public network access accidentally, other machines may reach your server.

- **Large in-memory AnnData input**:
  - An object is written once to a temporary `.h5ad` before the server binds,
    so it needs that much free disk space in the temporary directory; serving
    the file it came from skips the write.

- **Mixed content**:
  - Pointing an HTTPS Cellucid page at an `http://` server is refused by
    Cellucid itself, before the browser gets a chance to block the fetch:
    `An HTTPS Cellucid page requires an explicit HTTPS remote server URL`.
  - Always open the exact local Viewer URL printed by the server
    (`/?anndata=true`), which serves the verified UI generation.
  - Prefer an SSH tunnel when the kernel/server is remote.

- **Vector fields only exist in one dimension**:
  - If you have `velocity_umap` (for `X_umap`) but no `velocity_umap_3d` (for `X_umap_3d`), the overlay dropdown will be empty in 3D.
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
2) You pasted the **Viewer URL** into the **Link or server** bar. Its
   `?anndata=true` query makes the bar open it as a link to data, not as a
   server; type `127.0.0.1:8765`, the port alone, or the **Local URL**
   (`http://127.0.0.1:8765`).
3) The server is bound to `127.0.0.1` but you are trying to access it from another machine.
4) The Cellucid page is on HTTPS and the server address is `http://`, which
   Cellucid refuses outright.

**How to confirm**
- Open the server URL directly in your browser:

  ```text
  http://127.0.0.1:8765/_cellucid/health
  ```

  (If the server is running locally, you should get a small JSON response.)
- The dataset the server names is at
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
- Or write the file with `prepare(var_gene_id_column="...")`.

---

### Symptom: “Vector field overlay toggle is disabled / no fields appear”

**Likely causes (ordered)**
1) The dataset contains no vectors (common).
2) Vector fields exist, but are not named `<field>_<basis>` after an `X_<basis>` embedding the store holds (scVelo's convention: `velocity_umap` for `X_umap`), or do not have as many columns as it. A `velocity_umap_2d` beside a bare `X_umap` (no `X_umap_2d`) is not bound: rename it `velocity_umap`.
3) Dimension mismatch: vectors exist for 2D but you’re viewing 3D (or vice versa).

**How to confirm**
- Read step 2 of the startup output: its `obsm:` line lists every `obsm` key.
  Look for `<field>_<basis>` entries such as `velocity_umap` beside `X_umap`.

**Fix**
- Rename or regenerate vector fields using the exact keys documented in
  {doc}`../i_vector_field_velocity/01_what_vector_fields_are_user_facing`, or
  bind any `obsm` array that fits an embedding under
  **Session › Data mapping › Vector fields** ({doc}`13_data_mapping`).
- Switch the viewer to the dimension that has vectors.

For overlay UI behavior and deeper debugging, see:
- {doc}`../i_vector_field_velocity/index`
- {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

---

### Symptom: “The connectivity/edge overlay has nothing to draw”

**Likely causes (ordered)**
1) The object has no `obsp['connectivities']`, and no other graph was named.
2) The graph is stored under another key and was not named.
3) The file was written by `prepare()` without a `connectivities=` argument.

**How to confirm**
- Read step 2 of the startup output: its `obsp:` line lists every `obsp` key.

**Fix**
- Cause 1: compute the graph first, for example with `sc.pp.neighbors(adata)`.
- Cause 2: restart with `--knn KEY`, or `knn_key="KEY"` from Python, or pick it
  under **Session › Data mapping › Neighbour graph**. A matrix that is not
  symmetric, such as `distances`, is refused when the viewer reads it.
- Cause 3: write it again with `prepare(connectivities=...)`.

---

### Symptom: “It’s extremely slow”

**Likely causes**
- A high-latency link between the browser and the server, such as an SSH
  tunnel across continents: the viewer reads a store with many range requests,
  so each round trip adds up.
- The first gene of a row-major (CSR) expression matrix: the browser builds a
  gene index over the matrix once, and every later gene reads from it.

**Fix**
- Run the server close to the data (on the same machine as the file) and
  tunnel to it.
- For a matrix you open often, store it column-major (CSC), so one gene is one
  contiguous read. An in-memory object is written that way automatically.

## Next Steps

- Dataset identity (sessions/sharing): {doc}`06_dataset_identity_why_it_matters`
- Full troubleshooting matrix: {doc}`08_troubleshooting_data_loading`
- Want notebook embedding + programmatic control? → {doc}`05_jupyter_tutorial`
- Want browser-only loading without any server? → {doc}`03_browser_file_picker_tutorial`
- Want to publish a static public catalog instead? →
  {doc}`11_custom_dataset_repository`

Will it be usable at your scale? These two answer that, and are worth reading
before you hand the URL to a lab:

- What is read when, what it costs, and measured timings →
  {doc}`../../python_package/d_viewing_apis/14_performance_scaling_and_lazy_loading`
- What to turn on and off in the viewer at tens of millions of cells →
  {doc}`../n_benchmarking_performance/03_large_dataset_best_practices`
