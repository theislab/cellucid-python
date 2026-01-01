# Troubleshooting (prepare/export)

**Audience:** everyone  
**Time:** varies (treat this as a reference index)  
**Goal:** diagnose export failures quickly and fix them without guessing

This page is the deep troubleshooting guide for {func}`~cellucid.prepare`.

Use it when:
- `prepare()` raises an exception,
- export “succeeds” but the viewer looks wrong,
- the web app refuses to load the export folder,
- or exports are too slow/too big.

---

## How to use this page (recommended workflow)

1) Find the symptom that matches what you see.
2) Run the “confirm” checks (they are designed to be quick).
3) Apply the “safe fix” first (usually: re-export to a fresh folder, or `force=True`).
4) Only then try advanced fixes (quantization changes, vector scaling changes, etc.).

```{important}
The most common root cause is **row-order mismatch** across inputs.
If things look “plausible but wrong”, assume alignment is wrong until proven otherwise.
```

---

## Symptom: `prepare()` errors immediately (“missing required …”)

### Likely causes

- You did not pass `latent_space` (required).
- You did not pass `obs` (required).
- You did not pass any of `X_umap_1d/2d/3d` (at least one required).

### How to confirm

Read the exact exception message; the exporter raises clear errors such as:
- `ValueError: latent_space is required for outlier quantile calculation.`
- `ValueError: obs DataFrame is required.`
- `ValueError: At least one dimensional embedding must be provided.`

### Fix

Provide the missing arguments. A minimal working call is:

```python
prepare(
    latent_space=latent,   # (n_cells, n_latent_dims)
    obs=obs_df,            # n_cells rows
    X_umap_2d=X_umap_2d,   # (n_cells, 2)  (or X_umap_3d)
    out_dir="./export",
    force=True,
)
```

### Prevention

Use the preflight checklist:
- {doc}`02_input_requirements_global`

---

## Symptom: shape mismatch / `n_cells` mismatch errors

### Likely causes (ordered)

1) You filtered `adata` but reused an old embedding array (or vice versa).
2) You concatenated datasets and did not reorder all arrays consistently.
3) `gene_expression` is not aligned to the same `n_cells` as the embeddings.
4) `connectivities` was computed on a different subset.
5) Vector fields came from a different cell ordering than the embeddings.

### How to confirm

Print shapes for every cell-aligned input:

```python
print("emb:", X_umap_2d.shape if X_umap_2d is not None else None, X_umap_3d.shape if X_umap_3d is not None else None)
print("latent:", latent.shape)
print("obs:", obs.shape)
print("X:", None if X is None else X.shape)
print("conn:", None if conn is None else conn.shape)
```

Then assert they all match the same `n_cells` (rows).

### Fix (safe)

- Rebuild all inputs from the same source object (AnnData) after any filtering:

```python
adata2 = adata[mask].copy()
prepare(
    latent_space=adata2.obsm["X_pca"],
    obs=adata2.obs,
    X_umap_2d=adata2.obsm["X_umap_2d"],
    gene_expression=adata2.X,
    var=adata2.var,
    connectivities=adata2.obsp.get("connectivities"),
    out_dir="./export",
    force=True,
)
```

### Prevention

- Treat AnnData as the canonical alignment object.
- If you must use raw arrays, enforce a canonical cell id ordering and reorder everything explicitly before export.

---

## Symptom: `prepare()` succeeds but “nothing changed” (stale outputs)

### Likely causes

- You reused an `out_dir` and `force=False`, so existing files were skipped.
- You changed inputs/parameters but the manifest files already exist and were not rewritten.

### How to confirm

Look for console messages like:

```text
⚠ Skipping obs manifest: ... already exists (use force=True to overwrite)
⚠ Skipping points_3d.bin.gz: ... already exists (use force=True to overwrite)
```

### Fix (safe)

Pick one:

1) Use `force=True` while iterating, or
2) Export to a fresh `out_dir` each run.

### Prevention

Adopt an iteration convention:
- `exports/<dataset_id>/v001/`, `v002/`, etc. for major changes
- or `exports/tmp/` for scratch and `exports/final/` for publish

---

## Symptom: export folder is huge

### Likely causes

1) You exported too many genes (one file per gene).
2) `var_quantization` is disabled (float32).
3) `compression` is disabled.
4) You exported large string-heavy obs fields (large category lists).
5) You exported a dense connectivity graph (many edges).

### How to confirm

- Check if `var_manifest.json` exists; if yes, you exported genes.
- Count files in `var/` (rough proxy for number of genes exported).
- Inspect `dataset_identity.json["export_settings"]`.

### Fix

Safe, high-impact fixes:
- Export a curated gene list (`gene_identifiers`).
- Enable `var_quantization=8`.
- Enable `compression=6`.
- Remove large/high-cardinality obs fields from `obs_keys`.

### Prevention

- Start with a minimal export and grow it gradually.
- For large datasets needing full gene access, use server mode instead of static gene exports.

See:
- {doc}`06_gene_expression_matrix`
- {doc}`10_performance_tuning_guide_prepare_export`

---

## Symptom: export is extremely slow

### Likely causes

1) Many genes exported (file count dominates).
2) High gzip compression level on a slow CPU.
3) Exporting connectivities on a huge graph (edge extraction loops).
4) Writing to a slow filesystem (network mount, HDD).

### How to confirm

- If the progress bar says `Exporting genes`, gene export is the bottleneck.
- If it hangs at connectivities, the graph extraction is the bottleneck.

### Fix (in order)

1) Export fewer genes (`gene_identifiers`).
2) Disable compression during iteration.
3) Skip connectivities unless required.
4) Export to a local SSD.

### Prevention

- Use the iteration workflow in the performance guide.
- Precompute a stable gene panel for exports.

See: {doc}`10_performance_tuning_guide_prepare_export`

---

## Symptom: fields are missing in the UI after export

### Likely causes

1) You passed `obs_keys` and did not include the field.
2) `obs_manifest.json` was skipped due to `force=False` (stale manifest).
3) Two field keys collided after safe filename sanitization and overwrote files.

### How to confirm

- Open `<out_dir>/obs_manifest.json` and search for the key.
- Check whether a corresponding file exists under `<out_dir>/obs/`.
- Check the exporter logs for “Skipping obs manifest”.

### Fix

- Re-export to a fresh folder or set `force=True`.
- Make field names unique after sanitization.

### Prevention

- Export a curated `obs_keys` list and treat it as part of your pipeline spec.

See: {doc}`04_obs_cell_metadata`

---

## Symptom: gene expression missing / gene search disabled

### Likely causes

- You did not export `gene_expression` + `var`.
- `var_manifest.json` exists but is stale (skipped export).
- You exported a restricted `gene_identifiers` list and the requested gene isn’t included.

### How to confirm

- Does `<out_dir>/var_manifest.json` exist?
- Open it and inspect `fields` for the gene ID you expect.

### Fix

- Re-export genes with `force=True`.
- Ensure `var_gene_id_column` is what you think it is.
- Export the genes you need (or use server mode for full access).

See:
- {doc}`05_var_gene_metadata`
- {doc}`06_gene_expression_matrix`

---

## Symptom: gene expression looks wrong (flat / inverted / “wrong gene”)

### Likely causes

1) `var` rows do not match `gene_expression` columns (alignment bug).
2) You exported a different expression representation than you think (counts vs log1p).
3) Quantization + extreme outliers compress contrast.

### How to confirm

- Pick a gene and compare:
  - values in Python (`adata[:, gene].X`) vs
  - values loaded by the viewer (visually, or by inspecting the exported file if you know the format).
- Check that `var_gene_id_column` matches the IDs you’re using.

### Fix

- Fix alignment: rebuild `var` and matrix together (AnnData is recommended).
- Decide the expression representation explicitly and export that representation.
- If contrast is poor, try 16-bit or float32 export for that gene panel.

---

## Symptom: vector fields not detected / overlay missing

### Likely causes

1) You didn’t export `vector_fields` at all.
2) You exported vectors for 3D but didn’t export `X_umap_3d` (vectors skipped).
3) Vector field ids contain unsafe characters (export would error).
4) Vector arrays have wrong shape (not `(n_cells, dim)`).

### How to confirm

- Open `<out_dir>/dataset_identity.json` and search for `"vector_fields"`.
- Confirm `<out_dir>/vectors/` exists and contains `*.bin(.gz)` files.
- Check exporter logs for “Skipping vector field … embedding points_<dim>d not provided”.

### Fix

- Export vectors with explicit keys like `velocity_umap_2d`, `velocity_umap_3d`.
- Export matching embeddings.
- Re-export with `force=True`.

See: {doc}`08_vector_fields_velocity_displacement`

---

## Symptom: connectivities missing / graph-based features disabled

### Likely causes

- `connectivities` was not provided to export.
- Shape mismatch prevented export.
- `connectivity_manifest.json` is missing or stale.

### How to confirm

- Check for `<out_dir>/connectivity_manifest.json` and `<out_dir>/connectivity/edges.src.bin(.gz)`.

### Fix

- Export connectivities and re-export with `force=True`.

See: {doc}`07_connectivities_knn_graph`

---

## Symptom: export succeeds but the web app refuses to load the folder

### Likely causes

1) You selected the wrong folder in the file picker (one level too high/low).
2) `dataset_identity.json` is missing or malformed.
3) Export is partial (interrupted run; missing points file).

### How to confirm

Check that the folder root contains:
- `dataset_identity.json`
- at least one points file (`points_2d.bin(.gz)` or `points_3d.bin(.gz)`)

Web app loading docs:
- {doc}`../../web_app/b_data_loading/03_browser_file_picker_tutorial`
- {doc}`../../web_app/b_data_loading/08_troubleshooting_data_loading`

<!-- SCREENSHOT PLACEHOLDER
ID: invalid-export-folder-error
Where it appears: User Guide → Python Package → Data Preparation API → Troubleshooting
Capture:
  - Try to load an invalid folder (missing dataset_identity.json) in the web app
  - Capture the exact on-screen error message/toast/notification
Crop:
  - Include: the error message + the loader panel context
Redact:
  - Remove: any private paths/dataset names
Alt text:
  - Error message indicating the selected folder is not a valid Cellucid export.
Caption:
  - If the export folder is missing `dataset_identity.json`, the viewer cannot load it and shows an invalid export error.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for an invalid export folder error.
:width: 100%

If the export folder is missing `dataset_identity.json`, the viewer cannot load it and shows an invalid export error.
```

### Fix

- Re-export to a clean folder with `force=True`.
- Ensure you load the dataset folder root (the folder that contains `dataset_identity.json`).

### Prevention

- Add an automated “export validation” step to your pipeline (check for identity + points file).

---

## Reference: collecting debug info for maintainers

When reporting a bug, include:

### 1) Minimal reproducible inputs

- dataset scale (`n_cells`, `n_genes`)
- which inputs you provided (embeddings, obs, genes, connectivities, vectors)
- the exact `prepare()` call (parameters, especially quantization/compression)

### 2) The exact error message + stack trace

Copy/paste the full exception text.

### 3) Export folder tree (redacted)

Example:

```text
my_export/
  dataset_identity.json
  obs_manifest.json
  points_3d.bin.gz
  obs/
    leiden.codes.u8.gz
    leiden.outliers.u8.gz
    pct_mito.values.u8.gz
```

### 4) Key JSON files

Attach (or paste) the contents of:
- `dataset_identity.json`
- `obs_manifest.json`
- `var_manifest.json` (if present)
- `connectivity_manifest.json` (if present)

### 5) Environment

- OS + filesystem type (local SSD vs network mount)
- Python version
- `cellucid` version

---

## Next steps

- Re-read the global invariants: {doc}`02_input_requirements_global`
- Performance tuning: {doc}`10_performance_tuning_guide_prepare_export`
