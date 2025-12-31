# Dataset Identity and Reproducibility

**Audience:** computational users, maintainers, anyone sharing datasets  
**Time:** 10–15 minutes  
**Goal:** understand what makes an export “the same dataset” and how to make exports reproducible.

Every export folder includes:

- `dataset_identity.json`

This file is the viewer’s “source of truth” for:
- how many cells/genes exist,
- which embeddings are available (1D/2D/3D),
- what obs fields exist (categorical vs continuous),
- whether connectivity/vector fields are present,
- and what export settings were used (compression/quantization).

## What counts as “dataset identity” in Cellucid?

In practice, identity is a combination of:

1) **The dataset ID** (`id` in `dataset_identity.json`)
2) **The dataset structure** (n_cells, embeddings, field keys, etc.)
3) **The stable meaning of row order** (cell ordering across files)

If you want share links, session bundles, or community annotation to behave predictably, you need to treat identity as a first-class concept.

## How `cellucid-r` chooses `dataset_id`

If you do not provide `dataset_id=...`, `cellucid_prepare()` derives it from:
- `dataset_name` (if provided), else
- the folder name (`basename(out_dir)`),

and then sanitizes it to a filesystem-safe identifier.

### Recommendation

- **Always set `dataset_id` explicitly** for datasets you plan to share or publish.
- Use something stable and versioned, like:
  - `pbmc3k_v1`
  - `mouse_brain_atlas_2024-08`
  - `projectX_patient12_v2`

## Reproducibility checklist (practical)

### 1) Record your export script

Store the exact R code used to produce the export folder alongside the export (or in a repo).

### 2) Freeze the inputs

If your source object (Seurat/SCE) is mutable, make sure you know:
- which assay/slot you exported (counts vs logcounts vs scaled)
- which reduction you used (UMAP vs tSNE; 2D vs 3D)
- which metadata columns you exported

### 3) Make “row order” explicit

When exporting from Seurat/SCE:
- choose a cell ID vector (e.g. `cells <- colnames(seu)`)
- explicitly reorder every input to match that order

This prevents “it works on my machine” bugs when objects are subset/reordered.

### 4) Version your exports

If you re-export the same biological dataset after changes (filters, normalization, different embedding):
- keep the old folder (don’t overwrite silently),
- bump `dataset_id` or encode a version in `out_dir`,
- and write a short changelog for collaborators.

## Where dataset identity shows up in the UI

Cellucid uses identity information to:
- show dataset names and stats,
- decide what fields/embeddings exist,
- and (in broader workflows) reason about session compatibility.

For the web app’s perspective on sessions and compatibility, see:
- {doc}`../../web_app/l_sessions_sharing/07_versioning_compatibility_and_dataset_identity`

## Troubleshooting

### Symptom: “I exported again but the viewer still shows old fields”

**Likely cause**
- You reused the same `out_dir` and export skipped writing because `force=FALSE` (default).

**Fix**
- Re-run with `force=TRUE`, or delete the export folder and re-export.

### Symptom: “My collaborator’s session doesn’t apply to my export”

**Likely causes**
- Different `dataset_id`
- Different cell order / different `n_cells`
- Different field keys (renamed columns)

**Fix**
- Align on a stable dataset identity plan; don’t overwrite exports without versioning.
