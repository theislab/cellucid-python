# Dataset identity (why it matters)

**Audience:** everyone (wet lab, computational, developers)  
**Time:** 5–10 minutes  
**What you’ll learn:** what Cellucid treats as “the same dataset”, what must stay stable, and how to choose/verify a dataset id  
**Prerequisites:** none (helpful: you know which loading path you use)

---

## The one-sentence definition

In Cellucid, **dataset identity** is the app’s “primary key” for a dataset: a stable identifier (plus a few stability assumptions) that the UI uses to decide whether saved state, caches, and annotations belong to *this* dataset or a *different* one.

If dataset identity changes, Cellucid will behave as if you opened a new dataset—even if the points “look similar”.

---

## Fast path (non-technical): what you should do

1) Pick a **stable dataset id** once (e.g. `pbmc_10k`, `mouse_brain_v1`).
2) Keep it stable across re-exports *unless you intentionally want a new dataset*.
3) If you share datasets (GitHub exports, sessions, community annotation), treat the dataset id like a **permalink**: changing it breaks links between artifacts and the dataset.

---

## Where the dataset id comes from (by loading method)

### 1) Pre-exported folders (`prepare(...)`)

For exported datasets, the dataset id is:

- `dataset_identity.json["id"]`

This file is written by `cellucid.prepare(...)` into the export directory.

**Key property:** because it’s *file-backed metadata*, it stays stable as long as you keep exporting with the same `dataset_id`.

### 2) Public GitHub-hosted exports (no server)

Same as above: each dataset folder contains `dataset_identity.json`, and the id comes from `dataset_identity.json["id"]`.

Additionally:
- `datasets.json` is a *manifest* that lists datasets and their paths.
- The UI uses this to populate the dataset dropdown.

See {doc}`02_local_demo_tutorial` for the exact required layout.

### 3) Server mode / Jupyter (`.h5ad` / `.zarr` / in-memory AnnData)

When you serve or display AnnData directly, Cellucid still generates a `dataset_identity.json` payload—dynamically—via the AnnData adapter.

Default behavior:
- `dataset_name` defaults to the file name (e.g. `mydata.h5ad` → `mydata`) or “AnnData Dataset”.
- `dataset_id` defaults to a URL-safe slug derived from `dataset_name`.

Best practice (recommended when you care about persistence/sharing):
- Pass an explicit `dataset_id=...` (and optionally `dataset_name=...`) when using:
  - Python server mode (`serve_anndata(..., dataset_id="...")`)
  - Jupyter (`show_anndata(..., dataset_id="...")`)

Note: `cellucid serve ...` currently derives dataset identity from the path/filename; if you need a stable dataset id for sharing/sessions, prefer exports (`prepare`) or the Python API.

Why this matters:
- If you rename a `.h5ad` file, the default dataset id changes.
- Your saved sessions / caches / annotation links may no longer match.

---

## What must stay stable (and what can change)

### Must stay stable if you want things to “carry over”

**1) Dataset id**
- Used throughout the web app as the primary dataset key (URL selection, session fingerprints, caches).

**2) Cell ordering / row alignment**
- Cellucid assumes that arrays like embeddings, `obs` fields, and gene expression all refer to the same cell order.
- If you reorder cells and reuse the same dataset id, you risk “silent mismatch” bugs:
  - selections/highlights point to the wrong biological cells
  - cached results apply to the wrong rows

**3) Field keys**
- `obs` column names and gene identifier conventions should remain stable if you want saved “color-by”, filters, or analysis inputs to restore.

### Can change safely (usually)

- `dataset_name` and `description` (UI labels)
- export timestamp (`created_at`)
- minor export settings (compression level, quantization), as long as the dataset semantics are unchanged

Rule of thumb:
- If changes would alter the meaning of an existing selection/filter (different cells, different gene set, different IDs), treat it as a new dataset (new dataset id).

---

## Choosing a good dataset id (computational + lab-friendly guidance)

Good dataset ids are:
- **short**: `pbmc_10k`, not a full sentence
- **stable**: not “run_2025_12_30_19_04”
- **URL-safe**: lowercase letters, numbers, `_` or `-`
- **versioned only when needed**: add `v2` only when it truly becomes a new dataset

Common conventions:
- `project_sample` (e.g. `crc_patient12`)
- `study_dataset_v1` (e.g. `miller_atlas_v1`)
- `assay_tissue` (e.g. `scatac_mouse_brain`)

Avoid:
- spaces and special characters
- putting private identifiers into the dataset id (it often ends up in URLs)

---

## How to verify the dataset id (quick checks)

### Exported dataset folder

Open your export directory and inspect:

```text
<export_dir>/dataset_identity.json
```

The dataset id is:

```json
{
  "id": "pbmc_demo"
}
```

### AnnData server / Jupyter

If you have a Python handle, print the viewer URL and/or inspect the server’s `dataset_identity.json` endpoint.

Example (server/Jupyter):

```python
from cellucid import show_anndata

viewer = show_anndata("data.h5ad", dataset_id="pbmc_10k", dataset_name="PBMC 10k")
print(viewer.viewer_url)
```

---

## Edge cases (read if you share sessions or annotations)

- **Two different datasets with the same dataset id**:
  - The app may treat them as “the same dataset” for some caches/sessions.
  - At minimum, you risk confusing collaborators.
  - Use unique dataset ids for different datasets, even if their names are similar.

- **Dataset folder name vs dataset id** (GitHub exports):
  - The folder name can be arbitrary (e.g. `exports/pbmc_demo/`).
  - The dataset id is inside `dataset_identity.json["id"]`.
  - Keep them aligned to reduce confusion, but know they are not the same thing.

- **Changing `obs` column names**:
  - Anything that refers to field keys (color-by, filters, analysis inputs) can fail to restore.
  - Prefer additive changes (new fields) over renaming/removing fields.

---

## Troubleshooting (identity-related)

### Symptom: “My session won’t restore / says dataset mismatch”

**Likely causes**
- Dataset id changed (common when the `.h5ad` filename changed in AnnData mode).
- You re-exported a different dataset but reused the same dataset id.

**How to confirm**
- Compare the current dataset id to the saved one:
  - exported datasets: check `dataset_identity.json["id"]`
  - AnnData mode: check the file stem / passed `dataset_id`

**Fix**
- If you intended to keep identity: use an explicit `dataset_id=...` (server/Jupyter) or `prepare(..., dataset_id=...)` (exports).
- If you intended a new dataset: accept the mismatch and treat it as a new dataset (new session/annotation scope).

### Symptom: “Community annotation shows ‘no votes’ / ‘wrong dataset’”

**Likely causes**
- The annotation repo is keyed to a dataset id that no longer matches the dataset you opened.

**Fix**
- Ensure the dataset id you open in Cellucid matches the dataset id expected by the annotation repo.
- See {doc}`../j_community_annotation/index`.

---

## Next steps

- Picking a loading path: {doc}`01_loading_options_overview`
- Sharing datasets via GitHub exports: {doc}`02_local_demo_tutorial`
- Understanding export folder contents: {doc}`07_folder_file_format_expectations_high_level_link_to_spec`
- When loading fails: {doc}`08_troubleshooting_data_loading`
