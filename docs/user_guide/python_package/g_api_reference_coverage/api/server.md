# Server (browser tab + local HTTP server)

```{eval-rst}
.. currentmodule:: cellucid
```

This page documents **server-mode** usage: open Cellucid in a browser tab while a local server serves the data.

You have two equivalent entry points:
- **CLI:** `cellucid serve …` (recommended for most users; detects the format)
- **Python:** {func}`~cellucid.serve_anndata` / {class}`~cellucid.AnnDataServer`

Both take any file or folder Cellucid reads — an `.h5ad` file, an AnnData
Zarr store (folder or zipped), an `.h5mu` or MuData Zarr store, a `.loom`
file, a 10x `.h5` file or Cell Ranger `outs/` folder, a MatrixMarket folder
(also GEO's `GSM…_matrix.mtx.gz` naming), an `.rds`/`.rda`/`.RData` (gzipped
or not) — or, in Python, an in-memory `AnnData` or `MuData`.

---

## Fast path (beginner-friendly)

### A) Serve from the terminal

```bash
# .h5ad, Zarr store folder, or zipped Zarr store:
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
cellucid serve /path/to/data.zarr --dataset-name "My dataset" --dataset-id my-dataset
cellucid serve /path/to/data.zarr.zip --dataset-name "My dataset" --dataset-id my-dataset
# MuData, Loom, 10x, MatrixMarket, R objects:
cellucid serve /path/to/cite.h5mu --dataset-name "CITE-seq" --dataset-id cite --modality prot
cellucid serve /path/to/GSM4567_seurat.rds.gz --dataset-name "GSM4567" --dataset-id gsm4567
cellucid serve /path/to/outs --dataset-name "PBMC 1k" --dataset-id pbmc-1k
```

`--dataset-name` and `--dataset-id` default to the store's own
`uns['cellucid']['name']` and `uns['cellucid']['id']` (every file
{func}`~cellucid.prepare` writes declares both); name the dataset when the store
declares none.

Then open the printed URL (or let it auto-open).

### B) Serve from Python

```python
from cellucid import serve_anndata

serve_anndata("my_study.h5ad")  # a prepared file names itself; blocks until Ctrl+C
```

---

## Practical path (deployment patterns)

### 1) A prepared file or any AnnData

**A Cellucid-ready AnnData** (recommended for speed + reproducibility):
- Create with {func}`~cellucid.prepare` (see {doc}`prepare`)
- Serve with {func}`~cellucid.serve_anndata` or `cellucid serve my_study.h5ad`

**Any other AnnData**:
- Serve with {func}`~cellucid.serve_anndata` or
  `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`
- The viewer reads the data itself — in place, with HTTP range requests, where
  the format allows it, once from the start for an R object or a MatrixMarket
  folder. A file in priority order (`prepare()`'s) is served byte for byte; an
  in-memory AnnData, and an AnnData file in any other order, are served as a
  Zarr store in priority order whose chunks are made as they are asked for
  (`cellucid._served_store`), so a subsample transfers only its cells; a
  MuData object is first written once to a temporary Cellucid-ready `.h5mu`.
- Python reads the identity (`uns['cellucid']`) only of an AnnData or MuData,
  so an R object, a Loom, 10x or MatrixMarket dataset always needs
  `--dataset-name` and `--dataset-id`.

### 2) Local machine vs remote server (HPC / cloud VM)

**Local machine**:
- Keep default host `127.0.0.1` (only accessible on your machine).

**Remote server + SSH tunnel** (recommended for HPC):
1. On the remote machine (where the data lives):
   ```bash
   cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset --host 127.0.0.1 --port 8765 --no-browser
   ```
2. On your laptop:
   ```bash
   ssh -L 8765:127.0.0.1:8765 user@remote
   ```
3. In your laptop browser:
   - Open the printed Viewer URL:
     `http://127.0.0.1:8765/?anndata=true`

**Direct LAN access** (use carefully):
- Bind to `0.0.0.0` to make the server accessible to other machines on the network:
  ```bash
  cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset --host 0.0.0.0
  ```
- This is also what a compute node needs when the tunnel has to terminate on a
  different machine, such as a cluster login node — see
  {doc}`../../d_viewing_apis/17_hpc_slurm_and_compute_node_serving`.
- A non-loopback bind has **no authentication**: every user who can reach the
  port reads the dataset. Do this only if you understand your network/security
  posture.

### 3) Launch choices

Every launch option is optional: `embedding_keys`, `knn_key`,
`expression_source`, `gene_id_column`, `obs_keys`, `latent_key`,
`vector_field_default`, `normalize_embeddings`, `centroid_min_points`, and
`modality`. For an AnnData or MuData, each one you set is checked against the
store's metadata before the server binds (a named `obsm`, `obsp`, layer,
column or modality must exist); for a format only the browser reads (an R
object, Loom, 10x, MatrixMarket), only its form is checked here and the viewer
checks it against the data when it opens it. It is handed to the viewer as a
launch default. One you leave unset is
taken from the store's own `uns['cellucid']['defaults']`, else the viewer
applies its default rule; see
{doc}`../../../web_app/b_data_loading/13_data_mapping`.

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --knn connectivities
```

```python
from cellucid import serve_anndata

server = serve_anndata(
    "data.h5ad",
    dataset_name="My dataset",
    dataset_id="my-dataset",
    knn_key="connectivities",
)
```

No choice costs startup time: the viewer reads the neighbour graph the first
time a view draws it, and a vector field when the overlay shows it. Naming an
element the store does not have is an error before the bind, for example
`knn_key names obsp['connectivities'], which this store does not have.`

A choice you set overrides the store's own `uns['cellucid']['defaults']` for
this server only; the file is not changed.

### 4) Which URL the banner prints

A server bound to a wildcard address — `--host 0.0.0.0`, `--host ::`, or an
empty host — accepts connections on every interface, but the wildcard is not an
address any browser can open. So the banner prints the loopback origin as
`Local URL` and `Viewer URL` (that is what a browser on the serving machine
uses, and what an SSH tunnel terminates on), and names the addresses other
machines use separately:

A wildcard bind prints the loopback origin plus the machine's own name; see {doc}`../../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`.

- The `From another machine` block appears only when this machine has a usable
  name of its own. A hostname that resolves back to loopback names nothing
  outside the machine, so it is reported as no origin at all rather than as a
  URL that cannot travel.
- An IPv6 literal bind is bracketed as a URL requires: `http://[::1]:8765`.
- Open the exact Viewer URL the banner printed. `0.0.0.0` is never printed as a
  URL because nothing can connect to it.

---

## Deep path (how server mode works)

### What the server provides

- The data's bytes under `/_cellucid/anndata/<dataset_id><ending>` for one
  file (`.h5ad`, `.zarr.zip`, `.h5mu`, `.loom`, `.h5`, `.rds.gz`, …; the
  ending its format names it by), `<dataset_id>.zarr/<key>` for a Zarr store
  folder, or `<dataset_id>/<path>` for a folder of files, with `Range`,
  `206 Partial Content`, `Content-Range`, `HEAD`, and `ETag` answered exactly
  and the bytes sent by the kernel (`sendfile`). A Zarr store folder is served
  key by key (404 for an absent key), with consolidated metadata built in
  memory when the store has none. The route is named by the dataset id, so the
  local file name is never published.
- CORS headers so the viewer UI can request data
- `Host` validation ahead of routing: on a loopback bind (the default), only
  one well-formed authority naming `localhost` or a loopback IP literal on the
  bound port is served (421 otherwise), which is what defeats DNS rebinding.
  Behind a reverse proxy, declare its host name with `--allowed-host HOST` or
  `allowed_hosts=[...]`; declaring names also applies the check to a
  non-loopback bind, which otherwise answers any `Host`; see
  {doc}`../../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`
- Health/info endpoints:
  - `/_cellucid/health`
  - `/_cellucid/info`
  - `/_cellucid/datasets`
  - `/_cellucid/protocol` (where the notebook viewer channel is, and what it carries)
- The notebook viewer channel `/_cellucid/ws` (a WebSocket; opened only by a
  registered notebook viewer with its token, see
  {doc}`../../e_jupyter_hooks/05_architecture_viewer_channel`). The server
  answers no `POST` (405).
- `/_cellucid/datasets` is
  `{"datasets": [{"id", "name", "anndata": {"file": {"url", "kind", "layout", "size", "etag"}, "defaults": {...}}}]}`.
  `url` is the data's route relative to the server's base URL, `kind` is the
  input format registry's id (`h5ad`, `zarr`, `h5mu`, `mudata-zarr`, `loom`,
  `tenx-h5`, `tenx-outs`, `mtx`, `rds`, `rda`), `layout` is `file`, `zip`, or
  `folder`, `size` and `etag` are given for a file only (`etag` is the strong
  validator its bytes are answered with, read from the file at each request:
  the version the browser files what it keeps of the file under), a folder of
  files (`tenx-outs`, `mtx`) adds `files: [{"path", "size", "etag"}, …]`
  because a browser cannot list a folder, and `defaults` holds only the
  launch choices you set.
  `/_cellucid/health` is
  `{"status": "ok", "type": "anndata", "version", "format", "n_cells", "n_genes"}`
  (the counts `null` for a format Python does not read: an R object, Loom,
  10x or MatrixMarket), and `/_cellucid/info` adds the bound `host` and
  `port`, and `runtime`: the program serving the viewer, which the viewer
  names in its Dataset box (`"python"`; `"jupyter"` for a viewer in a notebook
  kernel; the R package's viewer says `"r"`). The server is
  threaded and keeps connections alive (HTTP/1.1), because a browser reads one
  file with many range requests. The full route table is in
  {doc}`../../../web_app/b_data_loading/04_server_tutorial`.
- `/_cellucid/info` always carries the package `version`.

### Viewer UI hosting and source requirements

By default, Cellucid establishes and serves the exact verified **viewer UI**
generation from the same server origin. This avoids mixed-content and
cross-origin issues for the browser tab. (A notebook widget, `show_anndata`,
uses no server: the browser loads the app from `web_source_url`.)

If the runtime cannot fetch and verify the exact source generation, startup
raises before binding the server. Use `--web-cache-dir PATH` or the
`web_cache_dir=` argument to select its publication directory.

### Startup progress (direct AnnData: four steps)

Startup prints four numbered steps: `[1/4] Detecting format`,
`[2/4] Reading metadata`, `[3/4] Checking launch choices`, `[4/4] Starting
server`. Steps 1 to 3 read metadata only, never an array. The transcript is in
{doc}`../../../web_app/b_data_loading/04_server_tutorial`.

- Step 2 prints the cell and gene counts and the `obsm`, `obs`, `obsp`, and
  `layers` entries an AnnData or MuData holds; for another format it says the
  browser reads the contents when it opens the data.
- Step 3 prints every launch choice it accepted, and
  `Consolidated metadata: built in memory (…)` for a Zarr store folder that
  has none.
- A store holding an array the browser cannot decode (an HDF5 filter other than
  gzip, shuffle, Fletcher-32, or LZF, such as Blosc, Zstandard, szip, or
  scaleoffset; a Zarr codec such as bz2, lzma, or zfpy) prints
  `Browser cannot decode:` with the arrays, and is served through the store in
  priority order, which decodes each chunk on the kernel as it is asked for
  (h5py or zarr with their codec plugins). MuData, Loom and 10x files are
  read on the kernel the same way.
- A MuData file not written by `prepare()`, a Loom file and a 10x Cell Ranger
  `outs/` folder print `Serving: a Zarr store in priority order, …` (or `a
  MuData Zarr store …`): the store holds what the browser presents from the
  file itself, cells in priority order (`cellucid._served_formats` mirrors the
  browser's layout views). A Loom file with no embedding, and a 10x `.h5`
  alone (its projections are in its `outs/` folder), are refused in the
  browser's words.

---

## API reference

### Functions

```{eval-rst}
.. autofunction:: serve_anndata
```

The viewer's web build is cached locally (see
{doc}`../../d_viewing_apis/09_server_mode_advanced`):

```{eval-rst}
.. autofunction:: get_web_cache_dir
.. autofunction:: clear_web_cache
```

### Classes


```{eval-rst}
.. autoclass:: AnnDataServer
   :members:
   :show-inheritance:
```

---

## Edge cases (do not skip)

### “It works locally but not on the remote server”
- If the server is on a remote machine, your browser cannot directly reach its localhost.
- Use SSH tunneling (preferred) or bind to an external interface + open firewall rules (higher risk).

### “I bound to 0.0.0.0 and now anyone can access it”
- Binding `--host 0.0.0.0` exposes the server to your network.
- There is no authentication on a non-loopback bind: every user who can reach
  the port reads the dataset. `--host`'s own help text states this.
- The banner still prints the loopback origin as `Local URL` / `Viewer URL`,
  plus a `Bound to every network interface. From another machine:` block when
  this machine has a routable name. `0.0.0.0` itself is never rendered inside a
  URL.
- Cellucid servers are designed for local/private use (not hardened for public internet exposure).

### “The graph is there, but the viewer will not draw edges”
- The viewer draws `obsp['connectivities']` by default. A graph
  stored under another key is drawn when you name it with `--knn KEY` (or
  `knn_key="KEY"`), or pick it under **Session › Data mapping**.
- A file written by {func}`~cellucid.prepare` holds a graph only when
  `connectivities=` was passed.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Port already in use”
Fix:
- Use `--port 0` if supported by your workflow (otherwise choose a free port).
- Or pick another port manually: `--port 9000`.

---

### Symptom: “Data directory not found”
Fix:
- Confirm the path exists on the machine where you ran the server.
- If you’re on a remote server, remember your local machine path is different.

---

### Symptom: “I can open `/_cellucid/health` but the dataset won’t load”
Likely causes:
- The store holds no plottable embedding in `obsm`, or a shape mismatch.
- A Zarr store folder's consolidated metadata is stale.

How to confirm:
- Visit `/_cellucid/datasets` and confirm the store route and launch choices.
- Read the viewer's own error; `viewer.debug_connection()` probes the store.

Fix:
- Name the embedding with `embedding_keys`, or write the file again with
  {func}`~cellucid.prepare`.

---

### Symptom: “CORS blocked” (browser console error)
Likely causes:
- A corporate environment / browser extension blocks cross-origin requests.

Fix:
- Prefer opening the viewer URL served by the same origin (default).
- Disable the conflicting extension for the site, or verify the request in a
  clean profile with extensions disabled.

---

### Symptom: web-generation startup failure
Likely causes:
- The configured source inventory or one of its declared objects could not be
  fetched or verified.

Fix:
- Ensure source access at startup and pass a writable `web_cache_dir`.

---

### Symptom: “Browser can’t connect / connection refused”
Fix:
- Confirm server is running and you are opening the correct URL.
- If remote: confirm your SSH tunnel is active and mapped to the correct port.

---

### Symptom: “It works for small datasets but gets slow for large ones”
Likely causes:
- A high-latency link between browser and server: the viewer reads a store
  with many range requests.
- The first gene of a row-major (CSR) expression matrix: the browser builds a
  gene index over the matrix once.

Fix:
- Run the server on the machine that holds the data and tunnel to it.
- Store a matrix you open often column-major (CSC), so one gene is one
  contiguous read.

---

### Symptom: `cellucid serve` stops on a module it could not import
What’s happening:
- An import failure is not one condition, and the message names which of three
  it is:
  - **The distribution is not installed.** The message names the package and the
    exact `pip install <name>` that supplies it.
  - **It is installed, but too old to hold the name Cellucid asked of it.**
    Installing the same version again repeats the failure, so the message names
    the file it imported and asks for `pip install --upgrade <name>`.
  - **It is a standard-library module this interpreter cannot provide.** No
    package installs it. The message names the running Python version, that
    interpreter’s own path, and the version range Cellucid supports.

Fix:
- Follow the exact remedy in the message; for the third case, run Cellucid on a
  Python inside the supported range instead of installing anything.

---

## See also

- {doc}`cli` for CLI options and examples
- {doc}`prepare` for `prepare(...)` and the Cellucid-ready AnnData
- {doc}`jupyter` for notebook embedding
