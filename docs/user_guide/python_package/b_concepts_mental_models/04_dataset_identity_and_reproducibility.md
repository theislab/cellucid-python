# Dataset identity and reproducibility

**Audience:** everyone (especially computational users and collaborators)  
**Time:** 15–30 minutes  
**What you’ll learn:**
- What “dataset identity” means in Cellucid (and why it matters)
- How `dataset_identity.json` is generated and used
- How session bundles decide whether state is safe to restore/apply
- A practical versioning strategy for “paper-ready” exports

---

## Mental model (one sentence)

If you want sessions, share links, and collaborative annotation to work reliably, treat your dataset export as an **immutable, versioned artifact** with a **stable dataset ID**.

---

## What is “dataset identity”?

In Cellucid, “dataset identity” is the answer to:

> “Is the dataset loaded right now the *same* dataset that this session / link / annotation was created for?”

Cellucid needs a dataset identity for:
- **Session restore** (web app): restore against the exact matching dataset or
  reject the bundle.
- **Session → AnnData application** (Python): mutate only the exact matching
  AnnData or raise before mutation.
- **Multi-dataset servers**: how the viewer lists and labels datasets behind one server URL.
- **Collaboration**: if you send someone a session bundle, you need a stable reference to the dataset it belongs to.

---

## The concrete identity file: `dataset_identity.json`

When you export with `prepare(...)`, Cellucid writes a `dataset_identity.json` into the export folder.

What it contains (high level):
- `id`: dataset ID (string)
- `name`: human-readable dataset name
- `description`: optional description
- `stats`: lightweight counts (cells/genes/fields) and connectivity summary
- `embeddings`: available dimensions + default dimension
- `export_settings`: compression/quantization knobs used for the export
- optional `source` metadata (name/url/citation)

```{note}
`dataset_identity.json` is a **metadata contract**, not a cryptographic hash of all content.
If you need tamper-proof guarantees, add your own checksums/version control around the export folder.
```

---

## How the dataset ID is chosen (important)

### Export mode (`prepare(...)`)

`prepare(...)` requires both `dataset_name` (human-friendly) and `dataset_id`
(stable, machine-friendly). It does not derive either value or rewrite an
unsafe ID.

Example:

```python
prepare(
    ...,
    out_dir="./exports/pbmc_v1",
    dataset_name="PBMC 10k (Seurat v3 tutorial)",
    dataset_id="pbmc_10k_v1",
    created_at="2026-07-26T12:34:56Z",
    dataset_description="PBMC demo export used for tutorial figures",
    source_name="Seurat PBMC tutorial",
    source_url="https://satijalab.org/seurat/articles/pbmc3k_tutorial.html",
    obs_categorical_dtype="uint16",
)
```

Omit `created_at` for an ordinary export: Cellucid writes the current UTC time.
Use an explicit, provenance-fixed value only when a canonical builder must
reproduce the complete export directory byte-for-byte. Its exact accepted form
is `YYYY-MM-DDTHH:MM:SSZ`.

### AnnData direct mode (`show_anndata(...)`)

Direct AnnData mode also requires an explicit stable name and ID:

```python
from cellucid import show_anndata
viewer = show_anndata("data.h5ad", dataset_name="My study (raw)", dataset_id="my_study_raw_v1")
```

---

## How sessions decide “same dataset or not”

A `.cellucid-session` bundle contains a lightweight **dataset fingerprint** in its manifest (not the full dataset).

In Python, `apply_cellucid_session_to_anndata(...)` uses this fingerprint as a safety guard:
- compares `cellCount` to `adata.n_obs`,
- compares `varCount` to `adata.n_vars`,
- compares `datasetId` to the required `expected_dataset_id`.

Any mismatch raises `ValueError` before AnnData is mutated. Application is
strict and atomic.

```{important}
These guards prevent the most dangerous failure mode: silently applying highlights/labels to the wrong cells.
They do not guarantee perfect identity (e.g., two datasets can have the same shape).
Use stable IDs + versioning to make identity meaningful.
```

---

## A versioning strategy that works (practical and boring — on purpose)

### Rule 1: treat export folders as immutable

Once you share an export folder (or publish it), do not modify it in place.
If you need to change anything, create a new folder and bump the version.

Good:
- `exports/pbmc_10k_v1/`
- `exports/pbmc_10k_v2/`

Risky:
- overwriting `exports/pbmc_10k/` while reusing the same dataset ID.

### Rule 2: version the dataset ID (not just the folder)

Make `dataset_id` reflect the version:
- `he_cortex_v2026_01_01`
- `hlca_lung_atlas_v3`
- `projectX_qcfiltered_v1`

### Rule 3: record provenance in `dataset_description` and/or `source_*`

If you want paper-grade reproducibility, include:
- which pipeline produced the AnnData,
- which commit/tag produced the export,
- which parameters (neighbors/umap) were used,
- any filtering decisions (QC thresholds).

Even a short description helps collaborators interpret sessions correctly.

### Rule 4: keep cell order stable across “same dataset”

Many Cellucid interactions use **row index** as identity.

If you reorder cells, previously saved sessions/highlights can map to different biological cells even if the shape is unchanged.

If you must reorder/subset:
- treat it as a new dataset version (new `dataset_id`),
- or do not reuse sessions across those changes.

---

## How this relates to collaboration tools

### Sessions (shareable state)

If you send a collaborator:
- a `.cellucid-session` file **without** the matching dataset export folder (or matching hosted dataset),
they may not be able to restore the state meaningfully.

### Community annotation (`cellucid-annotation`)

Community annotation workflows expect:
- stable dataset identity, so annotations can be applied consistently across contributors.

If you plan to use `cellucid-annotation`, decide your dataset ID/versioning early and do not change it casually.

---

## Edge cases (common identity mistakes)

### “I exported twice to the same folder with `force=True`”

This can produce a folder whose metadata says one thing while some files are from another export run.
If you must overwrite, do it only while iterating privately — and treat the result as not shareable until you export cleanly into a fresh folder.

### “Two exports have the same `dataset_id`”

This is a recipe for confusion:
- sessions may “match” when they shouldn’t,
- collaborators may restore the wrong state onto the wrong dataset.

Fix: choose unique IDs (include a version suffix).

### “My AnnDataViewer identity includes local file paths”

Direct AnnData identity may include source info (like the path).
Be mindful when sharing debug reports or screenshots (privacy/security covered in {doc}`06_privacy_security_and_offline_vs_online`).

---

## Troubleshooting

### Symptom: “Session restore says dataset mismatch”

Likely causes:
- different number of cells/genes,
- different dataset ID,
- you are loading a different export folder than you think.

How to confirm:
- inspect `dataset_identity.json` in the export folder you loaded,
- compare counts (`n_cells`, `n_genes`) to the session’s fingerprint (web app warning or Python logs).

Fix:
- load the correct dataset export folder,
- or create a fresh session from the dataset generation you intend to use.

The current web and Python readers do not accept a layout-only partial
restore. Identity or format mismatch rejects the complete operation.

### Symptom: “Applying a session raises a dataset fingerprint mismatch”

Likely causes:
- `expected_dataset_id` does not match the bundle's `datasetId`,
- cell or variable counts differ,
- you applied to an AnnData with different row order or subset.

Fix:
- pass the exact ID used to create the matching viewer as
  `expected_dataset_id=...`,
- apply only to the same dataset generation and row order, or regenerate the
  session from the dataset version you intend to use.

---

## Next steps

- Capture state into Python and apply safely: {doc}`05_sessions_to_anndata_bridge`
- Understand what persists where: {doc}`03_state_persistence_and_scope`
