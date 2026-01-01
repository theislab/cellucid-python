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
- `.zarr` served by Python: ✅ lazy-ish (chunked access)
- Browser file picker for `.h5ad`: ❌ usually not truly lazy (browser holds the file)

See the full matrix: {doc}`02_the_14_loading_options_breakdown`.

## `.h5ad` backed mode vs in-memory mode

When you run:

```bash
cellucid serve data.h5ad
```

Cellucid uses backed mode (lazy) by default.

If you run:

```bash
cellucid serve data.h5ad --no-backed
```

Cellucid loads the whole AnnData into memory.

Use `--no-backed` only when:
- the dataset is small enough that RAM is not a concern, and
- you want faster random access patterns at the cost of memory.

## Practical recommendations by dataset size (rule of thumb)

```{note}
Exact cutoffs depend on your machine, browser, and GPU. These are starting points.
```

- **Small (≤50k–100k cells):**
  - `show_anndata(adata)` is fine and very convenient.
- **Medium (100k–500k cells):**
  - prefer `show_anndata("data.h5ad")` (backed) or `"data.zarr"`.
  - exports become attractive if you view often.
- **Large (500k–millions):**
  - prefer exports (`prepare` → `serve/show`) or `.h5ad` backed / `.zarr` served by Python.
  - avoid in-memory AnnData.

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
- confirm you didn’t force in-memory mode (`--no-backed`)
- consider exporting once (`prepare`) and using export mode

More symptom-based help: {doc}`15_troubleshooting_viewing`.
