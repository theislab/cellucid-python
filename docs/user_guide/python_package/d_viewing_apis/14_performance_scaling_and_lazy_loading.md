# Performance, scaling, and lazy loading

This page explains what makes Cellucid fast (or slow) and how to choose a workflow that scales.

If you are new, you can treat this as “rules of thumb + common failure modes”.

## At a glance

**Audience**
- Computational users working with large datasets
- Anyone hitting slow loads, memory spikes, or browser crashes

## The three biggest performance drivers

1) **Number of cells (`n_cells`)**
   - affects point rendering, selection, and most UI interactions
2) **Number of genes (`n_genes`)**
   - affects gene search lists and the cost of gene expression queries
3) **Data format and loading mode**
   - exported vs `.h5ad` vs `.zarr` vs in-memory

## Exports vs AnnData: performance tradeoffs

### Exported directory (best performance)

Best for:
- repeated viewing,
- sharing,
- very large datasets.

Why it’s fast:
- data is already in the viewer’s on-disk format,
- manifests are compact,
- and gene/field data is served in a viewer-friendly layout.

### AnnData direct mode (best convenience)

Best for:
- interactive notebook exploration,
- quick sanity checks.

Tradeoffs:
- the server has to adapt your AnnData into Cellucid’s format on the fly,
- some operations (especially gene column reads) can be slower than exported mode.

## Lazy loading (what it means in practice)

### The problem

Gene expression is huge (`n_cells × n_genes`). Loading it all into browser memory is often impossible.

### The solution

Cellucid uses **lazy gene loading** where possible:
- the viewer fetches only the genes you request,
- and the server responds with one gene’s values at a time.

### Which modes are lazy?

- Exported directory: ✅ lazy by design
- `.h5ad` served by Python: ✅ lazy (backed mode by default)
- `.zarr` served by Python: the Zarr store is loaded eagerly, while browser gene
  requests remain on demand
- Browser file picker for `.h5ad`: ❌ usually not truly lazy (browser holds the file)

See the full matrix: {doc}`02_the_14_loading_options_breakdown`.

## Exact `.h5ad` loading mode

A `.h5ad` path is always opened with `anndata.read_h5ad(path, backed="r")` —
read-only, HDF5-backed, no option and no flag. The file is never read into
memory as a whole, so startup time and resident memory are set by the metadata
(`obs`, `var`, `obsm`), not by the size of the expression matrix.

```bash
cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
```

Startup prints five numbered steps; the transcript is in
{doc}`../../web_app/b_data_loading/04_server_tutorial`.

What backed mode costs is per-gene speed. For an in-memory CSR matrix the
adapter builds a CSC copy on the first gene request, which doubles the matrix's
memory but makes every later column read cheap. For a backed `.h5ad` it
deliberately does **not**: building that cache would pull the whole matrix into
memory and defeat the reason for backed mode. Backed columns are sliced straight
out of the file instead — slower per gene, flat in memory, which is the right
trade at multi-million-cell scale.

```{note}
`.zarr` is the other direction: `anndata.read_zarr()` loads the store eagerly at
startup, and only the browser's gene requests stay on demand.
```

## The neighbor graph is the expensive optional part

Connectivity is off unless you ask for it, in Python with
`serve_connectivity=True` and on the command line with `--connectivity`:

```bash
cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --connectivity
```

`--connectivity` is an AnnData-only flag. Supplying it for a prepared export is
an error, exactly like `--dataset-name` and `--obs-key`: an export either holds
the connectivity artifacts already or it does not.

**With connectivity off (the default)**

- `adata.obsp['connectivities']` is never read.
- `stats.has_connectivity` in `dataset_identity.json` is `false` and
  `stats.n_edges` is `null`.
- `connectivity_manifest.json` and the three payload routes —
  `connectivity/edges.src.bin`, `connectivity/edges.dst.bin`, and
  `connectivity/edges.weights.f64.bin` — answer 404.
- The viewer has no graph to draw, and startup skips the single most expensive
  thing it could do.

**With connectivity on**

- The whole graph is read and validated *before* the server binds its port,
  because `connectivity_manifest.json` is one of the first documents the viewer
  fetches and it declares both the edge count and the neighbor maximum. There
  is no way to serve a manifest that has not been computed.
- Asking for connectivity when `adata.obsp` holds no `connectivities` matrix is
  an error naming `sc.pp.neighbors(adata)`, not a silent absence.

The dataset report in step 3 therefore has three states, not two:

| Reported | Meaning |
| --- | --- |
| `Connectivity: yes (N edges)` | The graph was asked for, validated, and is served |
| `Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)` | A graph exists on the object and this server was not asked for it |
| `Connectivity: no` | The object carries no graph at all |

Distinguishing the middle state costs nothing: membership in `obsp` is a key
lookup, not a read of the matrix.

`prepare(...)` is unaffected. It was always opt-in through its own
`connectivities=` argument, and a prepared export served with `serve(...)`
publishes whichever connectivity artifacts it holds.

### What a neighbor graph actually costs

`sc.pp.neighbors(adata, n_neighbors=50)` writes a symmetric
`obsp['connectivities']`. Take five million cells:

- **Stored neighbors.** Fifty neighbors per cell is 250,000,000 directed
  neighbor relations. Symmetry stores the reverse of each relation that is not
  already mutual, so the matrix ends up holding between 250,000,000 and
  500,000,000 stored neighbors. `adata.obsp['connectivities'].nnz` tells you
  which, in constant time.
- **Deduplicated edges.** What gets served is the strict upper triangle, which
  is exactly half the stored count: 125,000,000 to 250,000,000 undirected
  edges.
- **The symmetry check.** Before deduplication the validator proves the matrix
  is exactly symmetric in topology *and* in weight, which means materializing
  the transpose and comparing it against the original — a second and a third
  copy of an object that already costs roughly 12 bytes per stored neighbor as
  CSR with Float64 data.
- **The coordinate pass.** Duplicate coordinates are detected by sorting one
  64-bit key per stored neighbor: 2 GB at 250,000,000. The row, column, and
  weight arrays it validates are another 24 bytes per stored neighbor, or 6 GB
  at the same size.
- **The ordering pass.** Edges are ordered by source and then destination. The
  permutation is one 64-bit index per edge — 1 GB at 125,000,000 edges — and
  the reordered arrays exist alongside the originals while the reorder runs,
  for roughly 3 GB more.
- **What survives.** The payload the server holds and the browser eventually
  downloads is 4 bytes of source index, 4 bytes of destination index, and 8
  bytes of Float64 weight per edge: 2 GB at 125,000,000 edges.

Above **50,000,000 stored neighbors** — exactly the shape of a 50-neighbor
graph over one million cells — the `cellucid.connectivity` logger names the
cost before the validation starts, so the wait is announced rather than
discovered:

```text
Reading a neighbor graph of 250,000,000 stored neighbors over 5,000,000 cells. Validating it is symmetric, deduplicating it into edges, and ordering them takes minutes and several times the graph's own memory at this size. Serve without connectivity to skip this work entirely.
```

The degree count that produces `max_neighbors` for the manifest uses
`np.bincount` rather than an unbuffered `np.add.at` scatter-add; on the
hundreds of millions of elements a graph this size brings to that line, the
difference was measured at roughly twenty-fold on that step. It is now among
the cheapest parts of the pass — the symmetry check and the ordering pass are
what dominate.

```{note}
None of this arithmetic applies when connectivity is off, which is the default.
If you are not going to draw the graph, none of it is paid.
```

## What the terminal shows while it works

The five steps exist so that a long phase is recognisable as work rather than a
hang. A default startup over a large backed `.h5ad`, with connectivity off:

```text
$ cellucid serve atlas.h5ad \
      --dataset-name "Cell atlas" --dataset-id cell-atlas-v1 --no-browser
cellucid v0.9.1

Importing dependencies (anndata, numpy, scipy)... done

[1/5] Detecting format...
      Path: atlas.h5ad
      Format: h5ad
      ✓ Format detected

[2/5] Loading AnnData...
      Mode: read-only backed h5ad
      Obs columns: classifying 18
      Embeddings: resolving obsm keys
      Embeddings: 2D from obsm['X_umap']
      ✓ File opened

[3/5] Analyzing dataset...
      Cells: 5,000,000
      Genes: 24,516
      Embeddings: 2D
      Obs fields: 15 categorical, 3 continuous
      Vector fields: no
      Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)
      ✓ Analysis complete

[4/5] Building manifests...
      Centroids: one per categorical field and dimension
      ✓ Manifests built

[5/5] Starting server...
      Viewer UI generation: establishing exact configured source
Cellucid web UI cache: 100%|██████████████████| 486/486 [00:13<00:00, 36.66file/s]
      ✓ Viewer UI generation established
      ✓ Server ready

════════════════════════════════════════════════════════════
  CELLUCID SERVER RUNNING
════════════════════════════════════════════════════════════

  Local URL:    http://127.0.0.1:8765
  Viewer URL:   http://127.0.0.1:8765/?anndata=true

  Press Ctrl+C to stop

════════════════════════════════════════════════════════════
```

The same object with `--connectivity` adds two lines to step 2, and the pause
between them is where the minutes go:

```text
[2/5] Loading AnnData...
      Mode: read-only backed h5ad
      Obs columns: classifying 18
      Embeddings: resolving obsm keys
      Embeddings: 2D from obsm['X_umap']
      Connectivity: reading obsp['connectivities']
      Connectivity: 125,481,922 edges, 214 neighbors at most
      ✓ File opened

[3/5] Analyzing dataset...
      Cells: 5,000,000
      Genes: 24,516
      Embeddings: 2D
      Obs fields: 15 categorical, 3 continuous
      Vector fields: no
      Connectivity: yes (125,481,922 edges)
      ✓ Analysis complete
```

Reading the transcript top to bottom tells you which phase you are waiting in:

- stuck after `Connectivity: reading obsp['connectivities']` — the graph
  validation described above, which is why it is off by default;
- stuck after `Obs columns: classifying N` — one validation pass per continuous
  obs column, usually storage-bound on a network filesystem;
- stuck inside `[4/5] Building manifests` — centroids, one per category per
  available dimension;
- stuck inside `[5/5] Starting server` — the viewer generation is being fetched
  and verified.

A prepared export served with `cellucid serve` still runs three numbered steps,
not five. It has no adapter to build.

## Practical recommendations by dataset size (rule of thumb)

```{note}
Exact cutoffs depend on your machine, browser, and GPU. These are starting points.
```

- **Small (≤50k–100k cells):**
  - Direct AnnData viewing is fine and very convenient.
- **Medium (100k–500k cells):**
  - prefer `show_anndata("data.h5ad", dataset_name="My dataset", dataset_id="my-dataset")`.
  - exports become attractive if you view often.
- **Large (500k–millions):**
  - prefer exports (`prepare` → `serve/show`) or a read-only-backed `.h5ad`.
  - avoid in-memory AnnData.
  - leave connectivity off unless the session will actually draw the graph.

For repeated sessions over a very large object, run `prepare(...)` once and
serve the export every time after. The adapter build in step 2 and the centroid
work in step 4 happen on every single direct-AnnData startup; an export pays
both once, at export time, and the prepared server then runs three steps that
read finished artifacts. The break-even is low — two or three sessions over a
multi-million-cell object is usually enough. Connectivity sharpens this
further: an export validates and writes the edge list once, while the direct
path revalidates the whole graph on every startup that asks for it.

## Browser/GPU considerations (web app side)

Cellucid’s rendering load is ultimately in the browser (WebGL).

Common symptoms of GPU pressure:
- “WebGL context lost”
- extreme slowdowns when switching to high-detail modes
- browser tab crashes

For web-app-specific performance tuning (render modes, LOD, etc.), see:
{doc}`../../web_app/n_benchmarking_performance/index`.

## Troubleshooting performance

If performance is unexpectedly bad:

- confirm you’re not using browser `.h5ad` loading (#4) for a huge dataset
- confirm direct H5AD input is being served through Python
- read the numbered steps and identify which phase the time is in
- confirm connectivity is off unless you need it — it is the longest phase of a
  large startup, and it is off unless you passed `serve_connectivity=True` or
  `--connectivity`
- consider exporting once (`prepare`) and using export mode

More symptom-based help: {doc}`15_troubleshooting_viewing`.
