# Performance profiling and scaling

This page explains where Cellucid Python performance bottlenecks typically come from and how to profile them.

It covers:
- export-time performance (`prepare`)
- server performance (exported vs AnnData)
- notebook performance (embedding + hooks)

---

## Performance “north star”: what users experience

Users mostly care about:

- time to first view (TTFV): “how long until I can rotate and color-by?”
- interaction latency: “does selecting/filtering feel instant?”
- memory footprint: “does this crash my laptop/server?”

Developer work should aim to improve one without catastrophically harming the others.

---

## Export-time performance (`prepare`)

### Main cost drivers

1) **Gene export is O(n_cells × n_genes)** (and writes one file per gene).
2) **Connectivity export can be expensive** if the KNN graph is large.
3) **Compression trades CPU for disk**.
4) **Quantization trades precision for disk** (usually worth it for visualization).

### Practical tuning knobs

- `var_quantization=8` is usually the best size/quality tradeoff for colormaps.
- `compression=6` is a reasonable gzip default.
- `gene_identifiers=[...]` dramatically reduces export size/time if you only need a gene panel.

### Memory footguns

- Dense gene expression matrices can explode memory.
- Sparse matrices are preferred; `prepare` will convert to CSC for per-gene access.

---

## AnnData server performance

AnnData mode is designed for convenience and lazy loading, but it has real bottlenecks:

- gene expression access is column-oriented (better with CSC)
- caches exist (LRU for recently accessed genes) but can still grow
- the HTTP server is single-process and largely sequential (standard library `HTTPServer`)

Where to look:
- `cellucid-python/src/cellucid/anndata_adapter.py`:
  - CSR→CSC conversion cache
  - gene expression LRU cache (default size 100)
  - centroid/outlier caches

Practical guidance:
- If users complain “gene coloring is slow”, check whether the backend is repeatedly converting formats or thrashing the LRU cache.

---

## Exported server performance

Exported mode is “as fast as your static files”:

- file size and compression dominate
- browser-side decompression of `.gz` files can become noticeable for huge arrays

If exported mode is slow:
- check export settings (quantization + compression),
- check file sizes,
- check client hardware limits (browser memory/GPU).

---

## Notebook embedding performance

Notebook mode adds overhead:

- iframe rendering and notebook layout
- proxy/tunnel URLs (remote notebooks)
- large hook payloads (huge selections)

Guidance:
- keep hook payloads lightweight when possible,
- avoid sending large “command” payloads frequently (e.g., repeated highlight updates every frame).

---

## Profiling playbook (developer steps)

### Step 1 — Measure first

Before optimizing, capture:
- dataset size (`n_cells`, `n_genes`)
- export settings (quantization/compression)
- mode (exported vs AnnData)
- rough timings (“prepare took 20 minutes”, “gene color-by takes 3 seconds”)

### Step 2 — Profile the Python side

Useful built-in tools:

```bash
python -m cProfile -o profile.out -m your_script_that_calls_prepare
python -c "import pstats; p=pstats.Stats('profile.out'); p.sort_stats('cumtime').print_stats(30)"
```

For memory:
- prefer tooling like `tracemalloc` in targeted code paths (avoid broad guesses).

### Step 3 — Profile HTTP behavior

In browser DevTools Network:
- check which endpoints dominate (var files? obs files? points?),
- check whether requests are repeated (caching issues).

### Step 4 — Decide optimization strategy

Common strategies:
- reduce data size (quantize/subset)
- reduce recomputation (cache)
- reduce per-request overhead (batch requests, avoid repeated conversions)

---

## Edge cases

### Millions of cells

Expect:
- large point arrays
- heavy browser memory pressure
- slow decompression if exported `.gz` is huge

Mitigation:
- encourage pre-export with quantization
- consider downsampling workflows in documentation

### Extremely high category counts

Even if export succeeds, categories with tens of thousands of levels can become unusable in the UI.

Mitigation:
- document recommended practice (collapse categories, use higher-level labels).

---

## Troubleshooting

### Symptom: “Export folder is enormous”

Likely causes:
- no quantization,
- no compression,
- exporting all genes.

Fix:
- enable `var_quantization` + `compression`,
- export fewer genes.

### Symptom: “AnnData server uses too much memory”

Likely causes:
- conversion to CSC doubles sparse matrix memory,
- caches storing large arrays.

Fix:
- prefer pre-export for large datasets,
- reduce gene interactions, or restart kernel to clear caches after exploration.
