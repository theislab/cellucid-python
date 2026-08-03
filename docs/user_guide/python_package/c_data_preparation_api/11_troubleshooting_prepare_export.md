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
- `TypeError: obs must be a pandas DataFrame, got dict.` (also raised for
  `var`, naming whatever type you actually passed)

### Fix

Provide the missing arguments. A minimal working call is:

```python
prepare(
    latent_space=latent,   # (n_cells, n_latent_dims)
    obs=obs_df,            # n_cells rows
    X_umap_2d=X_umap_2d,   # (n_cells, 2)  (or X_umap_3d)
    dataset_name="My study",
    dataset_id="my-study-v1",
    obs_categorical_dtype="uint16",
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
    dataset_name="My study",
    dataset_id="my-study-v1",
    obs_categorical_dtype="uint16",
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

## Symptom: `prepare()` refuses a gene or a continuous `obs` column as non-finite

### Likely causes

1) The expression matrix carries `NaN` or `±Inf` — usually from a division, a
   log of zero, or a merge that left cells unmeasured.
2) A continuous `obs` column carries `NaN` where a score was never computed.
3) A value is finite in your object but not once narrowed to `float32`: either
   too large (it becomes an infinity) or nonzero but smaller than the smallest
   `float32` (it becomes exactly zero).

Cellucid publishes continuous values as finite `float32` only, because a colour
scale has no position for an infinity: one `+Inf` makes the whole field's range
infinite and collapses every other cell onto a single colour.

### How to confirm

Read the exception. `NonFinitePayloadError` (a `ValueError`, so existing
`except ValueError` handlers still catch it) carries a **counted** diagnosis
rather than a bare "must be finite":

```text
Gene 'A2ML1' cannot be published: of 18,142,044 cells, 12 +Inf. First affected
cells: 4,127, 88,301, 250,994, 1,004,112, 7,730,015, ... Cellucid publishes
finite float32 only, because a colour scale has no position for an infinity.
Repair the matrix before serving or exporting it -- ...
```

The counts are the point. Twelve cells out of eighteen million and every cell in
the column are the same sentence in a message that only says "not finite", and
they want opposite responses: the first is a cell to look at, the second is the
wrong matrix. The message counts `NaN`, `+Inf`, `-Inf`, values beyond the
`float32` range, and values below the smallest `float32` separately, then names
the first five offending cell indices. A continuous `obs` column reports the
same shape under `Continuous obs field '<key>'`.

Every other numeric input `prepare()` takes — `X_umap_*d`, `latent_space`, a
vector field, `transition_matrix` — is refused with the same counts. Those have
axes rather than one value per cell, so their positions read as `[row, column]`,
zero-based and row-major, the way you would index the array you passed. A sparse
input is counted over its stored non-zero values, which have no position the
caller would recognise, so those messages give the counts and no positions.
`cellucid-r` prints the same list one-based and column-major, the way its caller
indexes; see
{doc}`../../r_package/g_api_reference_coverage/02_error_messages_and_exceptions_document_patterns`.

### Fix

Repair the values, then export again:

```python
import numpy as np

adata.X.data[~np.isfinite(adata.X.data)] = 0          # sparse expression
np.nan_to_num(adata.X, copy=False)                    # dense expression
adata.obs["score"] = adata.obs["score"].replace([np.inf, -np.inf], np.nan).fillna(0)
```

Deciding what the repair should be is yours: `0` is a value with a position on
the colour scale, and it is not the same claim as "not measured".

### Prevention

Export is the cheapest way to find *every* affected gene, because it scans them
all in one pass. The direct-AnnData server checks the same values but only for
the column being requested, so it refuses one gene per selection with an HTTP
422 — see {doc}`../d_viewing_apis/15_troubleshooting_viewing`. If you plan to
serve a large object repeatedly, run `prepare()` once first and treat the export
as the finiteness audit.

---

## Symptom: `prepare()` refuses a categorical field's outlier quantiles

### Likely causes

- `centroid_min_points` (default `10`) is higher than the size of every category
  in that field, so no category qualifies for a centroid and **every** generated
  quantile is missing. A set with no quantile at all has no bound to publish and
  nothing for the viewer to decode, so it cannot be encoded into the quantized
  payload.

Small, many-category fields (per-sample IDs, plate wells, fine-grained
sub-clusters) are the usual origin.

### How to confirm

The check runs per field but reports once, so the call names every affected
field at the same time, before anything is written:

```text
ValueError: 1 generated categorical outlier quantile set(s) cannot be encoded:
a set with no quantile at all has no value to publish.
  1 set(s) with obs_continuous_quantization=8:
    'sample_id': no category holds at least centroid_min_points=50 cells, so every generated quantile is missing
  Fix: ...
```

Check the named field with `adata.obs["sample_id"].value_counts()`.

### Fix

The message ends with the three options, and any one of them works:

1) lower `centroid_min_points` so a category qualifies,
2) drop the field from `obs_keys=`, or
3) pass `obs_continuous_quantization=None` to export the generated quantiles as
   full-precision `float32` (this is already the default; the failure only
   arises when you asked for 8- or 16-bit quantization).

### Prevention

Match `centroid_min_points` to the smallest category you actually want an
outlier filter for.

---

## Symptom: “Field 'x' has N categories, but uint8 supports at most 255.”

### Likely causes

- `obs_categorical_dtype="uint8"` with a field that has more than 255
  categories. `uint8` spends code 255 on the missing marker, so it carries 255
  categories (codes `0…254`); `uint16` spends 65535 the same way and carries
  65,535.

### How to confirm

```python
print(obs["batch"].astype("category").cat.categories.size)
```

### Fix

- Pass `obs_categorical_dtype="uint16"`, or
- drop the field from `obs_keys=` if it is a per-cell identifier rather than a
  grouping.

### Prevention

A field with thousands of categories is legal but not usable in a legend. Keep
cell IDs out of the export and map selection indices back in Python instead —
see {doc}`04_obs_cell_metadata`.

---

## Symptom: a non-empty target is rejected

### Likely causes

- You reused a non-empty `out_dir` with `force=False`.

### How to confirm

The call raises:

```text
FileExistsError: ... Use force=True to replace it.
```

### Fix (safe)

Pick one:

1) Use `force=True` for an intentional complete atomic replacement, or
2) Export to a fresh `out_dir` each run.

### Prevention

Adopt an iteration convention:
- `exports/<dataset_id>/v001/`, `v002/`, etc. for major changes
- or `exports/tmp/` for scratch and `exports/final/` for publish

---

## Symptom: `out_dir` is rejected as not a dedicated dataset output directory

### Likely causes

- You passed the directory you are working in (`.`, `./`, or `os.getcwd()`),
  your home directory (`~`), the filesystem root, or the directory holding
  every home (`..`, `/Users`, `/home`) — directly or through a symbolic link.

### How to confirm

The call raises, before anything is created or removed:

```text
ValueError: out_dir must name a dedicated dataset output directory, not /home/you.
```

### Fix (safe)

Name a child directory of your own, such as
`prepare(..., out_dir="./exports/my_dataset")`.

### Prevention

Publishing replaces the whole `out_dir`, and `force=True` removes everything the
previous generation held, so `out_dir` should never be a directory you keep
anything else in. Keep one directory per dataset under a single `exports/` root.

---

## Symptom: a category label or an identity string is rejected as “not shown as written”

### Likely causes

- A source annotation carries padding: `"Liver "` instead of `"Liver"`. Reading
  a fixed-width table, a hand-edited spreadsheet column, or a CSV exported from
  Excel are the usual origins.
- A byte-order mark (`U+FEFF`) leads the first label of a UTF-8 CSV column.
- A label is empty (`""`), often from `.fillna("")`.
- Two labels differ only by a run of whitespace: `"T cell"` and `"T  cell"`.

### How to confirm

```python
labels = adata.obs["organ"].astype("string").dropna().unique()
print([label for label in labels if label != label.strip()])
print([label for label in labels if any(ord(c) < 32 or 0x7f <= ord(c) <= 0x9f for c in label)])
```

### Fix

```python
adata.obs["organ"] = adata.obs["organ"].astype("string").str.strip()
adata.obs["organ"] = adata.obs["organ"].astype("category")
```

Check the result before re-exporting: if the column previously held both
`"Liver"` and `"Liver "`, stripping merges them into one category and moves
cells between them. That may be exactly what you want — but it is a change to
your annotation, which is why `prepare()` will not make it for you.

### Prevention

- Strip and re-factor categorical columns as part of the analysis pipeline, not
  at export time.
- Read CSVs with `encoding="utf-8-sig"` so a byte-order mark never reaches a
  label.

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
3) Exporting connectivities on a huge graph. The extraction is fully vectorised,
   but the symmetry proof materializes the transpose and compares it against the
   original, and the deduplicated edges are then sorted by source and
   destination — several times the graph's own memory, and minutes at hundreds
   of millions of stored neighbors.
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
2) The candidate export failed validation and was never published.
3) A field identifier was duplicated, or carried characters with no glyph
   (control, zero-width, or leading/trailing whitespace), so publication was
   rejected.

### How to confirm

- Open `<out_dir>/obs_manifest.json` and search for the key. Each entry is
  `[index, key, ...]`, whose first member names the file the field owns.
- Check whether `<out_dir>/obs/<index>.*` exists.
- Read the exact exporter exception.

### Fix

- Rename a duplicate, and strip invisible characters:
  `obs.columns = obs.columns.str.strip()`.
- Re-export to a fresh folder or set `force=True` for an intentional complete
  replacement.

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
2) You declared vectors for 3D but didn’t provide `X_umap_3d` (validation fails).
3) A declaration key does not match `<field>_umap_<1|2|3>d` exactly — the
   `_umap` tail and the dimension suffix are both required.
4) A vector field id is empty, duplicated, or carries characters the viewer
   cannot draw (control, zero-width, or leading/trailing whitespace). Punctuation
   and spaces are fine: an id is a name, not a filename.
5) Vector arrays have wrong shape (not `(n_cells, dim)`).

### How to confirm

- Open `<out_dir>/dataset_identity.json` and search for `"vector_fields"`.
- Confirm `<out_dir>/vectors/` exists and contains `*.bin(.gz)` files.
- Read the validation error for the vector key and its required embedding dimension.

### Fix

- Export vectors with required dimension-suffixed keys such as
  `velocity_umap_2d` and `velocity_umap_3d`.
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

- Check for `<out_dir>/connectivity_manifest.json` and every manifest-declared
  source, destination, and Float64 weight payload.

### Fix

- Export connectivities and re-export with `force=True`.

See: {doc}`07_connectivities_knn_graph`

---

## Symptom: “manifest does not describe the payloads that were written”

### What it means

Before publishing, `prepare()` re-expands every payload path the emitted
manifest points a reader at and compares that set against the directory it
actually wrote. A mismatch either way aborts the candidate, so the previous
generation stays untouched:

```text
RuntimeError: Observation manifest does not describe the payloads that were
written. Declared but absent: ['obs/7.values.u8.gz']. Written but undeclared: [].
```

A payload the manifest does not declare is invisible to the viewer; a payload it
declares but that was never written fails the dataset in the browser, long after
the export "succeeded". Both are caught here instead.

### Likely causes

- The export directory was written to by something else while `prepare()` ran.
- A filesystem that dropped or renamed a file (network mounts under pressure).
- A defect in Cellucid — this check exists to make it loud rather than silent.

### Fix

Re-export to a fresh `out_dir` on a local disk with nothing else writing to it.
If it reproduces, the message lists both sides of the mismatch exactly; include
it verbatim in a bug report along with the `prepare()` call.

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
    0.codes.u8.gz
    0.outliers.u8.gz
    1.values.u8.gz
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
