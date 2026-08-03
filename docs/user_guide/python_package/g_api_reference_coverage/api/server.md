# Server (browser tab + local HTTP server)

```{eval-rst}
.. currentmodule:: cellucid
```

This page documents **server-mode** usage: open Cellucid in a browser tab while a local server serves the data.

You have three equivalent entry points:
- **CLI:** `cellucid serve …` (recommended for most users; auto-detects format)
- **Python (exported data):** {func}`~cellucid.serve` / {class}`~cellucid.CellucidServer`
- **Python (AnnData):** {func}`~cellucid.serve_anndata` / {class}`~cellucid.AnnDataServer`

---

## Fast path (beginner-friendly)

### A) Serve anything (auto-detected) from the terminal

```bash
# Exported folder, .h5ad, or .zarr are all supported:
cellucid serve /path/to/data
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
cellucid serve /path/to/data.zarr --dataset-name "My dataset" --dataset-id my-dataset
```

Then open the printed URL (or let it auto-open).

### B) Serve an exported dataset from Python

```python
from cellucid import serve

serve("./my_export")  # blocks until Ctrl+C
```

---

## Practical path (deployment patterns)

### 1) Exported dataset vs AnnData mode

**Exported dataset** (recommended for speed + reproducibility):
- Create with {func}`~cellucid.prepare` (see {doc}`export`)
- Serve with {func}`~cellucid.serve` or `cellucid serve ./my_export`

**AnnData mode** (convenient, often slower):
- Serve with {func}`~cellucid.serve_anndata` or
  `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`
- Useful during exploratory analysis when you don’t want to write an export yet

### 2) Local machine vs remote server (HPC / cloud VM)

**Local machine**:
- Keep default host `127.0.0.1` (only accessible on your machine).

**Remote server + SSH tunnel** (recommended for HPC):
1. On the remote machine (where the data lives):
   ```bash
   cellucid serve /path/to/data --host 127.0.0.1 --port 8765 --no-browser
   ```
2. On your laptop:
   ```bash
   ssh -L 8765:127.0.0.1:8765 user@remote
   ```
3. In your laptop browser:
   - Open the printed prepared-export Viewer URL:
     `http://127.0.0.1:8765/?source=remote`

**Direct LAN access** (use carefully):
- Bind to `0.0.0.0` to make the server accessible to other machines on the network:
  ```bash
  cellucid serve /path/to/data --host 0.0.0.0
  ```
- This is also what a compute node needs when the tunnel has to terminate on a
  different machine, such as a cluster login node — see
  {doc}`../../d_viewing_apis/17_hpc_slurm_and_compute_node_serving`.
- A non-loopback bind has **no authentication**: every user who can reach the
  port reads the dataset. Do this only if you understand your network/security
  posture.

### 3) Connectivity is asked for, not assumed (AnnData mode)

Serving an AnnData object directly leaves the `obsp['connectivities']` neighbor
graph **off by default**. Ask for it explicitly:

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --connectivity
```

```python
from cellucid import serve_anndata

server = serve_anndata(
    "data.h5ad",
    dataset_name="My dataset",
    dataset_id="my-dataset",
    serve_connectivity=True,
)
```

**Off (the default).** `adata.obsp['connectivities']` is never read. The dataset
identity reports `has_connectivity: false` with `n_edges: null`, and
`connectivity_manifest.json` and the three `connectivity/edges.*` routes —
`edges.src.bin`, `edges.dst.bin`, and `edges.weights.f64.bin` (Float64 weights)
— answer 404.

**On.** The whole graph is read and validated *before the server binds*, because
the manifest the viewer fetches first declares the edge count and the neighbor
maximum. Asking for connectivity on an object whose `obsp` holds no
`connectivities` matrix is an error at startup, not silence.

**Why off is the default.** On a large object, building the edge list is the
single longest part of startup — a 50-neighbor graph over millions of cells is
hundreds of millions of stored neighbors — and most sessions never draw the
graph.

`--connectivity` and `serve_connectivity=` apply to direct AnnData input only.
Passing `--connectivity` with a prepared export directory is rejected, like
every other AnnData-only flag: an export already holds the connectivity
artifacts or it does not, and {func}`~cellucid.serve` publishes exactly what is
there.

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
- A prepared-export server prints the same shape with its own viewer query,
  `/?source=remote`.
- An IPv6 literal bind is bracketed as a URL requires: `http://[::1]:8765`.
- Open the exact Viewer URL the banner printed. `0.0.0.0` is never printed as a
  URL because nothing can connect to it.

---

## Deep path (how server mode works)

### What the server provides

- Static or virtual endpoints for dataset files (exported data or AnnData-backed “virtual files”)
- CORS headers so the viewer UI can request data
- `Host` validation ahead of routing: only one well-formed authority naming
  `localhost` or a loopback IP literal on the bound port is served, which is what
  defeats DNS rebinding. Behind a reverse proxy, declare its host name with
  `--allowed-host HOST` or `allowed_hosts=[...]`; see
  {doc}`../../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`
- Health/info endpoints:
  - `/_cellucid/health`
  - `/_cellucid/info`
  - `/_cellucid/datasets`
  - `/_cellucid/protocol` (which wire capabilities this `cellucid` accepts)

### Viewer UI hosting and source requirements

By default, Cellucid establishes and serves the exact verified **viewer UI**
generation from the same server origin. This avoids:
- mixed-content issues (HTTPS notebook + HTTP localhost)
- cross-origin iframe restrictions

If the runtime cannot fetch and verify the exact source generation, startup
raises before binding the server. Use `--web-cache-dir PATH` or the
`web_cache_dir=` argument to select its publication directory.

### Startup progress (direct AnnData: five steps)

Opening a large object spends minutes before the socket is bound, so each phase
reports as it happens:

Startup prints five numbered steps; the transcript is in {doc}`../../../web_app/b_data_loading/04_server_tutorial`.

- Step 2 reports the adapter build itself: the obs columns being classified, the
  `obsm` key resolution and the exact key resolved for each dimension, the
  vector-field scan, and — only when you asked for it — the connectivity read.
- Step 4 is the manifest and centroid phase. It has its own step because it is
  real work on a large object, not a formality after the load succeeds.
- `Connectivity` in step 3 has three states, not two:
  - `yes (N edges)` — the graph was asked for, read, and validated;
  - `not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)` — a graph exists and was not asked for;
  - `no` — the object has no graph. Checking `obsp` membership is a key lookup,
    not a read of the matrix, so naming the middle state costs nothing that
    opting out was meant to save.

A **prepared export** runs three steps instead, unchanged:
`[1/3] Validating dataset`, `[2/3] Loading dataset info`, `[3/3] Starting
server`.

---

## API reference

### Functions

```{eval-rst}
.. autofunction:: serve
```

```{eval-rst}
.. autofunction:: serve_anndata
```

### Classes

```{eval-rst}
.. autoclass:: CellucidServer
   :members:
   :show-inheritance:
```

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
- In AnnData mode the neighbor graph is opt-in. Without `--connectivity` (or
  `serve_connectivity=True`) the server never reads `obsp['connectivities']`,
  and the connectivity manifest and edge routes answer 404.
- Step 3 of startup says so explicitly: `Connectivity: not served
  (obsp['connectivities'] is present; pass serve_connectivity=True, or
  --connectivity, to draw it)`.
- A prepared export is different: it serves the connectivity artifacts it
  contains, so re-export with `connectivities=` if they are missing.

### “The dataset directory contains multiple datasets”
- {class}`~cellucid.CellucidServer` supports serving a directory that contains multiple exported datasets as subfolders.
- `server.viewer_url` opens that served catalog without embedding an arbitrary
  dataset id. One unique entry is selected automatically; multiple entries
  require an exact dataset-id selection.
- Confirm the exact ids and paths by visiting `/_cellucid/datasets`.
- Once any immediate child is recognized as a dataset candidate, **every**
  immediate subdirectory must be one complete current export. A stray or
  partial subdirectory rejects the root instead of being silently omitted.

### “Exported folder missing some files”
- `CellucidServer` validates the complete declared prepared-artifact inventory
  during construction.
- Missing required metadata, point files, or any artifact declared by a
  current manifest raises before a socket is bound.

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
- You served a directory that is not a valid exported dataset (missing manifests/points files).
- You served a parent folder with multiple datasets but opened the wrong path.

How to confirm:
- Visit `/_cellucid/datasets` and confirm the dataset path you intended exists.
- Try fetching `dataset_identity.json` and `obs_manifest.json` in the browser.

Fix:
- Serve the correct dataset directory (the folder produced by `prepare(...)`).
- Or re-export with {func}`~cellucid.prepare`.

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

### Symptom: “It works for small datasets but crashes/gets slow for large ones”
Likely causes:
- AnnData mode is generating virtual files on demand and may hit memory/CPU limits.

Fix:
- Export with {func}`~cellucid.prepare` (quantization + compression) and serve the exported directory.
- Use read-only-backed `.h5ad` input when staying in AnnData mode; Zarr input
  is loaded eagerly.

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

### Symptom: direct AnnData startup says “Importing dependencies…” and feels stuck
What’s happening:
- The first import may be slow in fresh environments because large scientific dependencies are loaded.

Fix:
- Wait once (subsequent runs are faster in the same process).
- If it truly hangs, run with `--verbose` and check for import errors.

---

## See also

- {doc}`cli` for CLI options and examples
- {doc}`export` for `prepare(...)` and export layout
- {doc}`jupyter` for notebook embedding
