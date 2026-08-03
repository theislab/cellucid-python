# Troubleshooting (viewing)

This page is symptom-driven troubleshooting for:
- `cellucid serve …` (CLI server mode),
- `serve(...)` / `serve_anndata(...)` (Python server mode),
- `show(...)` / `show_anndata(...)` (notebook embedding).

If you’re new, skim the “Quick triage” section first.

---

## Quick triage (do these in order)

### 1) What workflow are you using?

- Terminal/browser → you’re in **server mode** (`cellucid serve …`)
- Notebook iframe → you’re in **notebook mode** (`show` / `show_anndata`)

This matters because notebook mode adds extra failure modes (HTTPS mixed content, notebook proxies).

### 2) Is the server alive?

Open (replace `<port>`):

```text
http://127.0.0.1:<port>/_cellucid/health
```

If this doesn’t load, the viewer can’t load either.

### 3) Is the exact viewer generation being served?

Open the exact **Viewer URL** printed by the server:

```text
Prepared export: http://127.0.0.1:<port>/?source=remote
Direct AnnData:  http://127.0.0.1:<port>/?anndata=true
```

A wildcard bind prints the loopback origin plus the machine's own name; see
{doc}`12_remote_servers_ssh_tunneling_and_cloud`.

If this fails:
- exported folders may be incomplete, or
- AnnData server may have crashed while generating metadata.

### 4) Did one field fail, or did everything fail?

A viewer that loads, draws points, and refuses exactly one gene or one
continuous `obs` column is not a broken server. The direct-AnnData server checks
each continuous column for finiteness before it sends a byte, and answers
**HTTP 422** with a JSON diagnosis when the column cannot be drawn. Fetch the
route by hand and read the body:

```bash
curl -s "http://127.0.0.1:<port>/var/0.values.f32" | python -m json.tool
```

Use `GET`, not `curl -I`: a `HEAD` returns the `422` status with no body, and
the body is the whole diagnosis. See
[a gene or column that cannot be drawn](#symptom-a-gene-or-column-is-found-but-cannot-be-drawn-http-422).

### 5) If you are in a notebook: run a structured debug report

```python
viewer.debug_connection()
```

---

## Symptom: “Port already in use” / server won’t start

### Likely causes (ordered)
- another Cellucid server is still running
- another program is using the port (often `8765`)

### How to confirm
- startup raises an address-in-use `OSError`

### Fix
- choose a port explicitly: `cellucid serve ... --port 9000`
- in Python, pass `port=0` to request one operating-system-assigned port and
  then use the actual `viewer.server_url`
- in notebooks, call `viewer.stop()` on old viewers (see {doc}`11_viewer_lifecycle_cleanup_ports_and_multiple_viewers`)

### Prevention
- stop servers when done
- use fixed ports only when you need them (SSH tunnels)

---

## Symptom: viewer startup cannot establish the web generation

This failure occurs before the viewer server is published.

### Likely causes (ordered)
- your environment cannot reach the configured `web_source_url`
- the inventory or a declared object has an invalid MIME type, length, hash,
  file set, or build identity
- the selected `web_cache_dir` is not writable

### How to confirm
- read the exception raised by `cellucid serve`, `show(...)`, or
  `show_anndata(...)`
- exercise the same establishment contract directly with
  `ensure_web_ui_cached(force=True)`

### Fix
- allow access to the configured source and correct the exact response failure
- pass a writable generation location:

```bash
cellucid serve /path/to/data --web-cache-dir /path/to/cache
```

### Prevention
- verify source access and destination permissions before a workshop or job
- do not plan on a previous local generation being substituted when the source
  is unavailable

---

## Symptom: Viewer loads but dataset is blank / stuck loading

### Likely causes (ordered)
- export directory is missing required files (manifests or points binaries)
- you passed the wrong directory (parent folder vs dataset folder)
- AnnData mode: no readable UMAP embedding in `obsm`, or a shape mismatch

### How to confirm
- open `/_cellucid/datasets` and see what it lists
- open `/dataset_identity.json` and check:
  - `stats.n_cells`
  - `embeddings.available_dimensions`
- in AnnData mode: check `adata.obsm.keys()` and shapes

### Fix
- for exports: validate your folder layout (see {doc}`07_exported_directory_mode_show_and_serve`)
- for AnnData: confirm `obsm` carries an embedding Cellucid can read — a bare
  `X_umap` of 1, 2, or 3 columns, or an explicit `X_umap_1d` / `X_umap_2d` /
  `X_umap_3d` (see {doc}`08_anndata_mode_show_anndata_and_serve_anndata`)

---

## Symptom: my UMAP is in `X_umap` and it used to be refused

### What happens now
- an object that declares none of `X_umap_1d`, `X_umap_2d`, or `X_umap_3d` but
  carries the bare `X_umap` that `sc.tl.umap()` writes is served as written
- that one array is read at the dimension its own column count states: one
  column is 1D, two columns are 2D, three columns are 3D
- nothing in your object is renamed, assigned, or written back — resolution is a
  read, and a served `.h5ad` is opened read-only-backed

### How to confirm
- the startup transcript names the array it resolved:

```text
      Embeddings: resolving obsm keys
      Embeddings: 2D from obsm['X_umap']
```

- or open `/dataset_identity.json` and read `embeddings.available_dimensions`

### Good to know
- one array yields one dimension. To offer 2D *and* 3D in the same session,
  declare `X_umap_2d` and `X_umap_3d`
- a dimensional key anywhere in `obsm` decides the set on its own, and the bare
  key then stays out: an object holding both `X_umap_3d` and `X_umap` serves 3D
  only
- the same rule applies to `show_anndata(...)`, `serve_anndata(...)`,
  `AnnDataViewer(...)`, `AnnDataAdapter(...)` including `AnnDataAdapter.from_file`,
  `add_transition_drift_to_obsm(...)`, `cellucid serve` on a path it detects as
  direct AnnData input, and the browser’s own `.h5ad` and `.zarr` readers
- `prepare()` shares none of it: coordinates arrive there as arrays in its
  `X_umap_1d=`, `X_umap_2d=`, and `X_umap_3d=` arguments, and it never looks at
  `obsm` keys

Full rules: {doc}`08_anndata_mode_show_anndata_and_serve_anndata`.

---

## Symptom: “No UMAP embedding Cellucid can read was found in adata.obsm”

### Likely causes (ordered)
- `X_umap` is present but is some other width — a latent space named after a
  plot, and no rename makes it drawable
- UMAP was never computed, so `obsm` carries no embedding at all
- a dimensional key exists but does not carry the column count its name
  declares, or its rows do not match `adata.n_obs`

### How to confirm

```python
print(list(adata.obsm.keys()))
for k in ["X_umap", "X_umap_1d", "X_umap_2d", "X_umap_3d"]:
    if k in adata.obsm:
        print(k, getattr(adata.obsm[k], "shape", None))
print("n_cells:", adata.n_obs)
```

The message itself lists the `obsm` keys the object has, and names the shape of
an `X_umap` whose width no viewer can draw:

```text
No UMAP embedding Cellucid can read was found in adata.obsm. Cellucid reads the exact keys 'X_umap_1d', 'X_umap_2d', and 'X_umap_3d', and reads a bare 'X_umap' as the dimension its own column count states. Available obsm keys: ['X_pca', 'X_umap'].
  adata.obsm['X_umap'] has shape (3696, 10), and Cellucid draws 1, 2, or 3 dimensions, so that array names a dimension no viewer renders. Assign the columns you mean to draw under the key for their own count.
```

### Fix
- assign the columns you actually mean to draw under the key for their own
  count — when some other `obsm` array has `n_obs` rows and 1, 2, or 3 columns,
  the message prints one copy-pasteable statement per candidate:

```python
adata.obsm["X_umap_2d"] = umap_coordinates_2d
```

- or compute UMAP and store it under `X_umap_2d` / `X_umap_3d` directly
- a key whose array contradicts its own name is an error, not a
  reinterpretation: `X_umap_3d` holding `(n_obs, 2)` fails construction rather
  than being served as 2D
- serving a `.h5ad` or `.zarr` path? Re-save the file after the assignment

---

## Symptom: “obs field '…' has unsupported dtype”

### Likely causes (ordered)
- `adata.obs` carries a `datetime64`, timedelta, or period column; Cellucid
  serves categorical, boolean, numeric, string, and object columns only

### How to confirm

```python
print(adata.obs.dtypes)
```

### Fix
- convert the column: `adata.obs["collection_date"] = adata.obs["collection_date"].astype(str)`
- or leave it out by naming the columns to serve:

```python
show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
    obs_keys=["cell_type", "leiden", "total_counts"],
)
```

- from the CLI, repeat `--obs-key` once per column

---

## Symptom: Gene search returns nothing / “Gene 'X' not found”

### Likely causes (ordered)
- your gene IDs are not in `var.index` (default)
- wrong `gene_id_column`
- the chosen gene names are duplicated, or carry invisible characters; AnnData
  construction rejects them

### How to confirm

```python
print(adata.var.index[:5])
print(adata.var.columns)
```

### Fix
- pass `gene_id_column="..."` to `show_anndata(...)` / `serve_anndata(...)`
- ensure the chosen gene ID field is unique

This symptom is about a gene the server never found. A gene it *did* find and
then refused is the next section.

---

## Symptom: a gene or column is found but cannot be drawn (HTTP 422)

The gene is in `var.index`, the search offers it, you select it — and the viewer
answers with a notification titled `Gene cannot be shown`. Everything else keeps
working.

### What happened

Cellucid publishes continuous values as finite float32, because a colour scale
has no position for an infinity: one `+Inf` makes a field's range infinite and
collapses every other cell onto one colour. The direct-AnnData server checks each
continuous column before it sends a byte and answers **HTTP 422** with a JSON
diagnosis instead of serving it.

On this path neither an infinity nor a `NaN` is published, and both produce this
same 422. `NaN` is drawable only when the browser reads your `.h5ad` or `.zarr`
itself, where Cellucid renders it as the neutral grey meaning "not measured
here"; the Python server publishes finite float32 and nothing else. If the body
below counts `nan` and nothing else, check that column first — an all-`NaN`
column usually means it was never computed rather than that it broke.

### How to confirm

Fetch the route with `GET` and read the body; `curl -I` sends a `HEAD`, which
returns the status with no body:

```bash
curl -s "http://127.0.0.1:<port>/var/0.values.f32" | python -m json.tool
```

The `application/json` body names the field, counts the offenders by kind
(`nan`, `positive_infinity`, `negative_infinity`, `float32_overflow`,
`float32_underflow`), gives `total` and `offending`, and lists the first
offending cell indices in `examples`:

```json
{
  "error": "non_finite_continuous_payload",
  "kind": "gene",
  "name": "A2ML1",
  "counts": { "total": 18142044, "offending": 12, "positive_infinity": 12 },
  "examples": [4127, 88301, 250994, 1004112, 7730015],
  "message": "Gene 'A2ML1' cannot be published: of 18,142,044 cells, 12 +Inf. …"
}
```

Read `counts` before acting: twelve cells out of eighteen million and every cell
in the column arrive as the same status and want opposite responses. The two
float32 counts are values that are finite in your object and stop being finite —
or stop being nonzero — once narrowed to float32. `kind` is `"gene"` for a
`/var/…` route and `"field"` for `/obs/…`.

### Fix

Repair the values in the object, then reload:

```python
import numpy as np

adata.X.data[~np.isfinite(adata.X.data)] = 0          # sparse expression
np.nan_to_num(adata.X, copy=False)                    # dense expression
adata.obs["score"] = adata.obs["score"].replace([np.inf, -np.inf], np.nan).fillna(0)
```

Serving a `.h5ad` or `.zarr` path rather than an in-memory object? Re-save after
the repair.

:::{note}
A prepared export cannot reach this state: `prepare()` refuses the same values
with the same counted diagnosis. Exporting once is also how to find every
affected gene at once — the export scans them all, while the server refuses one
per request as you select it.
:::

---

## Symptom (notebook): iframe is blank, unreachable, or mixed-content blocked

### Likely causes (ordered)
- notebook served over HTTPS
- kernel/server is remote (JupyterHub/cloud)
- the browser-facing proxy/tunnel base was not supplied as `client_server_url=`

### How to confirm
- compare `viewer.viewer_url` with the URL that the browser can actually reach
- inspect the browser console for mixed-content or connection errors

### Fix
- expose the selected port through an HTTPS proxy or tunnel
- pass that exact browser-reachable base as `client_server_url=` when creating
  the viewer

Deep dive: {doc}`10_notebook_widget_mode_advanced`.

---

## Symptom: Vector field / velocity overlay is missing

### Likely causes (ordered)
- vectors not present (exports: missing `vector_fields` metadata; AnnData: missing `obsm` keys)
- naming mismatch (e.g. `velocity_umap2d` vs `velocity_umap_2d`)
- dimension mismatch (you have 2D vectors but are viewing 3D)

### How to confirm
- exports: open `/dataset_identity.json` and search for `vector_fields`
- AnnData: inspect `adata.obsm.keys()` for `*_umap_2d` / `*_umap_3d`

### Fix
- rename keys to the supported convention
- ensure you have vectors for the dimension you’re viewing

---

## Symptom: the neighbor-graph overlay is missing

### Likely causes (ordered)
- connectivity is **opt-in and off by default** on the direct-AnnData path, so a
  `.h5ad` straight out of `sc.pp.neighbors()` serves without a graph
- exports: the export was written without `connectivities=`, so it holds no
  connectivity artifacts to serve

### How to confirm
- read the `Connectivity` line of the step 3 dataset report, which has three
  states:

```text
      Connectivity: yes (33,476 edges)
      Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)
      Connectivity: no
```

- or open `/dataset_identity.json` and read `stats.has_connectivity` and
  `stats.n_edges`; with connectivity off they are `false` and `null`
- `connectivity_manifest.json`, `connectivity/edges.src.bin`,
  `connectivity/edges.dst.bin`, and `connectivity/edges.weights.f64.bin` answer
  404 while it is off

### Fix
- in Python, pass `serve_connectivity=True`:

```python
viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    serve_connectivity=True,
)
```

- from the terminal, add `--connectivity`:

```bash
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1 \
  --connectivity
```

- the keyword is accepted by `show_anndata(...)`, `serve_anndata(...)`,
  `AnnDataViewer(...)`, `AnnDataServer(...)`, `AnnDataAdapter(...)`, and
  `AnnDataAdapter.from_file`; `--connectivity` is AnnData-only and is rejected
  for a prepared export, like every other AnnData-only flag
- asking for a graph the object does not carry is now an error, not silence:
  `Connectivity was asked for, and adata.obsp has no 'connectivities' matrix to
  serve.` Compute one with `sc.pp.neighbors(adata)` first

### Expect it to cost time
Reading the graph is the single longest part of starting a server on a large
object — a 50-neighbor graph over millions of cells is hundreds of millions of
stored neighbors — which is why it is off unless asked for. When asked for, the
whole graph is read and validated **before** the server binds, because the first
document the viewer fetches is the manifest, and the manifest declares the exact
edge count and neighbor maximum.

`prepare()` is unchanged and was already opt-in through its own
`connectivities=` argument, and a prepared export serves whatever artifacts it
holds.

---

## Symptom: Performance is extremely slow / browser crashes

### Likely causes (ordered)
- using browser `.h5ad` loading instead of Python-backed mode
- loading a large in-memory AnnData or Zarr store
- huge categorical fields or heavy overlays

### Fix
- use `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`
- export once (`prepare`) and use export mode
- reduce or remove pathological fields (very high cardinality)

Performance guide: {doc}`14_performance_scaling_and_lazy_loading`.

---

## Reporting a bug (what to include)

Please capture:

- the exact command you ran (`cellucid serve ...`), including flags
- whether the input was export vs `.h5ad` vs `.zarr`
- the full startup transcript: the numbered steps, the `Embeddings: … from
  obsm[…]` line, the `Connectivity` line, and the banner
- the last `Error:` line verbatim — it is written after everything else, and it
  names the one condition to correct
- the output of:
  - `/_cellucid/health`
  - `/_cellucid/info`
- if one gene or one continuous field failed: the full JSON body of the `422`
  response for its `/var/<index>.values.f32` or `/obs/<index>.values.f32` route,
  fetched with `GET` — it already names the field, counts the NaN, `+Inf`,
  `-Inf`, float32-overflow, and float32-underflow values, and lists the first
  offending cell indices
- in notebooks: `viewer.debug_connection()` output
- if safe: `dataset_identity.json` (or at least its `stats` block)

This will save you (and maintainers) a lot of back-and-forth.

## Server-side symptoms

A banner URL that will not open, a startup that prints nothing for minutes, and
an import Cellucid says is missing are covered in
{doc}`../i_troubleshooting_index/04_server_mode_issues`. Cluster-specific
failures are in {doc}`../i_troubleshooting_index/07_hpc_and_remote_access_issues`.
