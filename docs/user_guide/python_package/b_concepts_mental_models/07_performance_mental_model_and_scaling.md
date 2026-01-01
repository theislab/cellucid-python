# Performance mental model and scaling

**Audience:** everyone (computational users get the most value)  
**Time:** 15–30 minutes  
**What you’ll learn:**
- What costs happen at export time vs view time
- Why `prepare(...)` exists (and when it’s worth it)
- Where bottlenecks usually come from (CPU, disk, network, GPU)
- Practical scaling recommendations and edge cases

---

## Mental model (one sentence)

Cellucid performance is best when you do expensive work **once** (export/quantize/compress) and keep view-time work mostly to **streaming + GPU rendering**.

---

## Two performance regimes: export-first vs direct AnnData

### Export-first (`prepare(...)` → `show(...)` / `serve(...)`)

What gets “paid once”:
- converting arrays to viewer-optimized binary layouts,
- quantizing values (optional),
- gzip compression (optional),
- writing manifests/metadata.

What gets “paid every time you view”:
- reading and streaming those binaries to the browser,
- GPU rendering in the viewer.

This is the recommended path for:
- large datasets,
- sharing/hosting,
- reproducible viewing.

### Direct AnnData (`show_anndata(...)`)

What gets “paid at view time”:
- dynamic conversion from AnnData to the viewer’s resource format,
- server-side work per gene query / field access,
- potential extra memory structures for fast gene access.

This is the recommended path for:
- quick exploration,
- iterating on preprocessing before committing to an export.

---

## The bottleneck map (what can be slow)

### 1) Disk I/O

- Export mode: reading `.bin(.gz)` files from disk.
- AnnData mode: reading from `.h5ad` (HDF5) or `.zarr` chunks.

Symptoms:
- slow first load,
- stutter when switching fields/gene expression.

### 2) CPU (serialization + decompression + preprocessing)

- gzip decompression in the browser and/or server
- quantization/dequantization (mostly at export time; minimal at view time)
- centroid/outlier computations (export time)

Symptoms:
- high CPU usage during load,
- UI responsiveness drops when scrubbing controls quickly.

### 3) Network transfer (even on localhost)

The viewer loads data via HTTP requests.

For large datasets, file size + compression level matter:
- higher compression → smaller transfer but more CPU,
- lower compression → larger transfer but less CPU.

### 4) GPU / WebGL

Rendering millions of points is GPU-bound.

Symptoms:
- low FPS,
- “WebGL context lost”,
- smoke/volumetric modes are especially heavy.

Web app performance guidance lives here:
- {doc}`../../web_app/n_benchmarking_performance/index`

---

## Practical recommendations by dataset size (rule of thumb)

```{note}
These are guidelines, not hard limits. Hardware, browser, and field complexity matter.
```

### Small (≤ ~50k cells)

- `show_anndata(...)` is usually fine.
- Export is optional unless you want sharing/reproducibility.

### Medium (~50k–500k cells)

- Export-first is strongly recommended for repeated viewing.
- Use quantization + compression if you will share/host:
  - `var_quantization=8`
  - `obs_continuous_quantization=8`
  - `compression=6`

### Large (≥ ~500k cells)

- Prefer export-first for reliable UX.
- Avoid browser file picker for huge folders unless you’ve tested the browser/OS combination.
- Consider server mode + SSH tunneling for remote workflows.

---

## Performance knobs you can control (Python-side)

### Export knobs (`prepare(...)`)

Common “good defaults” for large-ish datasets:

```python
prepare(
    ...,
    compression=6,
    var_quantization=8,
    obs_continuous_quantization=8,
)
```

Tradeoffs:
- More compression → smaller files, slower export, potentially slower load on weak CPUs.
- Less compression → larger files, faster export, potentially faster load on strong networks/SSDs.

### AnnData mode knobs (`show_anndata(...)`)

AnnData mode optimizes for interactive convenience, but there are real costs:
- the adapter may build extra structures for fast gene access (e.g., CSC caches),
- gene expression results are cached in an LRU (helps repeated gene queries),
- backed mode (`.h5ad`, `.zarr`) can reduce RAM usage compared to fully in-memory data.

If you hit memory ceilings, export-first is usually the right move.

---

## Edge cases (performance-specific)

### “It’s fast on my laptop but slow on the cluster”

Remote notebooks add complexity:
- proxy/tunnel overhead,
- weaker CPU per user,
- shared disks.

Mitigations:
- use a fixed port and SSH tunnel,
- export to a local SSD on the compute node if possible,
- prefetch the web UI cache once (see {doc}`06_privacy_security_and_offline_vs_online`).

### “My export folder is huge”

Most often:
- you did not quantize gene expression (`var_quantization=None`),
- you exported many genes (or didn’t subset),
- you used no compression.

Fix:
- export fewer genes if appropriate,
- use `var_quantization=8` and `compression=6`,
- consider excluding rarely used fields.

### “Gene expression is slow in AnnData mode”

Likely causes:
- huge sparse matrix with expensive column access patterns,
- cold cache (first query for a gene is always slower than repeated queries),
- running on network-mounted storage.

Fix:
- export-first for repeated viewing,
- ensure the `.h5ad` is on fast storage,
- keep gene queries focused (avoid rapid “scrolling” through thousands of genes).

---

## Troubleshooting

### Symptom: “The first load takes forever”

Likely causes:
- web UI assets are being fetched for the first time (online + cache miss),
- very large points file (millions of cells),
- slow disk.

How to confirm:
- Run `viewer.debug_connection()` and check `web_ui.cache` and `viewer_index_probe`.
- Watch network requests in the browser devtools network tab.

Fix:
- prefetch the web UI cache once,
- export-first with compression/quantization,
- move data to faster storage.

### Symptom: “Interactions lag when I scrub sliders”

Likely causes:
- expensive recomputation in the viewer (filters, outlier thresholds),
- CPU-bound rendering options.

Fix:
- use coarser steps (avoid continuous scrubbing),
- reduce simultaneous overlays,
- switch to a lighter render mode.

---

## Next steps

- Debugging checklist: {doc}`08_debugging_mental_model_where_to_look`
- Export-time performance tuning: {doc}`../c_data_preparation_api/10_performance_tuning_guide_prepare_export`
