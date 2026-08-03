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

The important word is **every**. Nothing in that list is written down between
sessions, which is the whole argument for exporting once.

Two of those costs are worth carrying as a mental model, because they are the
ones that scale:

- **Building the adapter** — classifying obs columns, resolving the embedding,
  computing one centroid per category per dimension. Proportional to the object,
  redone on every direct-AnnData startup, and absent from the prepared-export
  path entirely.
- **Reading the neighbor graph** — proportional to stored neighbors, not to
  cells, so it grows fastest of all. It is therefore asked for rather than
  assumed: off unless you pass `serve_connectivity=True`, or `--connectivity`.

Startup reports five numbered steps so a long pause has a label above it rather
than looking like a hang; the second is the adapter build and the fourth is the
centroid and manifest pass. A prepared export runs three steps, because the
artifacts already exist.

For the step-by-step cost breakdown, the graph arithmetic, and the transcript,
see {doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading`.

---

## The bottleneck map (what can be slow)

### 1) Disk I/O

- Export mode: reading `.bin(.gz)` files from disk.
- AnnData mode: on-demand reads from read-only-backed `.h5ad`, or reads from
  the in-memory AnnData produced by eager `anndata.read_zarr`.

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
- Leave connectivity off unless the session will draw the graph.

Direct AnnData is the right tool while preprocessing is still moving. Once the
object stops changing, export it: an export pays the adapter build and the
centroids once instead of once per session.

---

## Performance knobs you can control (Python-side)

### Export knobs (`prepare(...)`)

Common “good defaults” for large-ish datasets:

```python
prepare(
    ...,
    dataset_name="My study",
    dataset_id="my-study-v1",
    obs_categorical_dtype="uint16",
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
- gene expression results are cached in an LRU bounded by a 256 MiB budget (helps
  repeated gene queries without growing without limit),
- `.h5ad` paths are opened read-only-backed; `.zarr` paths are loaded eagerly.

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
- ensure the configured viewer source and selected generation directory have
  adequate throughput (see
  {doc}`06_privacy_security_and_offline_vs_online`).

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
- the complete configured web generation is being fetched and verified for
  this startup,
- very large points file (millions of cells),
- slow disk.

How to confirm:
- Run `viewer.debug_connection()` and check `web_ui.cache` and `viewer_index_probe`.
- Watch network requests in the browser devtools network tab.

Fix:
- diagnose source download or generation-directory throughput,
- export-first with compression/quantization,
- move data to faster storage.

### Symptom: “Startup sat there for minutes with nothing happening”

Likely causes:
- connectivity was asked for on a very large object, and the graph is being
  validated,
- centroids are being computed for many categorical fields across several
  dimensions,
- the viewer generation is being fetched on a cold cache.

How to confirm:
- read the last line printed. Each of the five steps names the phase it is in,
  and step 2 names each part of the adapter build as it starts.

Fix:
- drop `serve_connectivity=True` / `--connectivity` if the session will not draw
  the graph,
- narrow the served obs columns with `obs_keys=` / `--obs-key`,
- export once and serve the export for every session after.

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
