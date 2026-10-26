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
2) Write it into the file once: `prepare(..., dataset_id=...)` stores it as
   `uns['cellucid']['id']`, and every route then opens the file under it.
3) Keep it stable when you write the file again *unless you intentionally want a new dataset*.
4) If you share datasets (catalogs, sessions, community annotation), treat the dataset id like a **permalink**: changing it breaks links between artifacts and the dataset.

---

## Where the dataset id comes from (by loading method)

### 1) The file itself: `uns['cellucid']['id']`

A Cellucid-ready AnnData (written by `cellucid.prepare(...)` in Python or
`cellucid_prepare()` in R) states its id in `uns['cellucid']['id']`. Any
AnnData file can carry it: `adata.uns["cellucid"] = {"id": "pbmc_10k"}`. A
Seurat object carries it in `misc$cellucid` and a SingleCellExperiment in
`metadata$cellucid`, both read as `uns['cellucid']`.

**Key property:** because it is stored in the file, the same file has the same
id from every route: opened from your disk, from a link, from a catalog, or
served by Python. A session saved on one route restores on the others.

### 2) A file that states no id

| Route | Dataset id |
|---|---|
| **Open…** or a drop (`.h5ad`) | `h5ad_<file name without .h5ad>` |
| **Open…** or a drop (`.zarr.zip`, `.zarr` folder) | `zarr_<store name without .zarr / .zip>` |
| **Open…** or a drop, any other format | `<format id>_<name without its ending>`: `loom_atlas`, `h5mu_citeseq`, `rds_pbmc`, `tenx-h5_sample`, `mtx_<folder name>` |
| A link in the **Link or server** bar, one file or store | the link's own stable id (from the address or the record) |
| A dataset catalog (`datasets.json`) | the catalog entry's `id` |

A file without an id therefore has a different id on each route; write one into
it to make sessions and annotations portable.

### 3) A dataset catalog (sample datasets, linked folders)

`datasets.json` lists each dataset under an `id`. When the file states its own
id, the two must be equal: a file listed under another id is refused when it is
chosen, naming both. See {doc}`07_file_format_expectations`.

### 4) Server mode / Jupyter (`.h5ad` / Zarr store / in-memory AnnData)

When you serve or display AnnData directly, the server names the dataset in
`/_cellucid/datasets` (its `id` and `name`) and serves the store under a route
named by the dataset id, so the local file name is never published. The
browser reads the store under that id and name.

Exact behavior:
- `dataset_name` is non-empty, unpadded, and free of control characters;
  human-readable Unicode is preserved.
- `dataset_id` uses the same portable identifier as `uns['cellucid']['id']`:
  1–180 ASCII characters, beginning with a letter or digit, followed only by
  letters, digits, `.`, `_`, or `-`, with no trailing `.` or reserved Windows
  device name.
- Both default to what the store declares, `uns['cellucid']['name']` and
  `uns['cellucid']['id']`. A store that declares neither needs both
  `dataset_name` and `dataset_id` (`--dataset-name` and `--dataset-id` for
  `cellucid serve <path>`); without them startup stops with
  `dataset_name and dataset_id are required: this AnnData declares no Cellucid identity of its own (uns['cellucid']['id'] and ['name']).`
- The server's id and name are the dataset's: a server that names an id the
  file does not state serves the file under the server's id.
- `dataset_name`, `dataset_description`, and the `source_*` strings are shown
  to the reader exactly as stored, so none of them may carry a control
  character, a zero-width character (`U+200B`, `U+2060`, `U+FEFF`), or leading
  or trailing whitespace of any kind. Category labels obey the same rule. Both
  writers reject the value and name the offending character rather than
  trimming it, because trimming rewrites text nobody asked to change.

Use the same exact `dataset_id` whenever you reopen the same dataset
generation. Changing it intentionally creates a different dataset identity,
which prevents a session for the old identity from being applied implicitly.

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

**3) Field names**
- `obs` column names and the names you publish your genes under should remain stable if you want saved “color-by”, filters, or analysis inputs to restore. A session records field names, so a field that is renamed between writes of a file is a field a saved session can no longer find.

### Can change safely (usually)

- `name` and `description` (UI labels), citation and source
- how the file is stored (`.h5ad` or Zarr, chunking, compression), as long as the dataset's values are unchanged

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
- `study_dataset_v1` (e.g. `smith_atlas_v1`)
- `assay_tissue` (e.g. `scatac_mouse_brain`)

Avoid:
- spaces and special characters
- putting private identifiers into the dataset id (it often ends up in URLs)

---

## How to verify the dataset id (quick checks)

### In the file

```python
import anndata as ad

adata = ad.read_zarr("pbmc_demo.zarr")   # or ad.read_h5ad("pbmc_demo.h5ad")
print(adata.uns.get("cellucid", {}).get("id"))
```

In the web app, the **Dataset** box under **Session** shows the dataset that
is open; a saved session records its id in its fingerprint.

### AnnData server / Jupyter

If you have a Python handle, print the viewer URL and/or inspect the server’s `/_cellucid/datasets` endpoint, which names the dataset's `id` and `name`.

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

- **File name vs dataset id**:
  - The file name can be anything (e.g. `pbmc_demo_2026.zarr`).
  - The dataset id is `uns['cellucid']['id']` (or the catalog's id for it).
  - Keep them aligned to reduce confusion, but know they are not the same thing.

- **Changing `obs` column names**:
  - Anything that refers to field keys (color-by, filters, analysis inputs) can fail to restore.
  - Prefer additive changes (new fields) over renaming/removing fields.

---

## Troubleshooting (identity-related)

### Symptom: “My session won’t restore / says dataset mismatch”

**Likely causes**
- The `dataset_id` given at launch (or declared in the store's
  `uns['cellucid']['id']`) changed.
- You wrote a different dataset but reused the same dataset id.
- The file states no id, and you opened it on another route than the one the
  session was saved on (each route then gives it its own id).

**How to confirm**
- Compare the current dataset id to the saved one:
  - the file: check `uns['cellucid']['id']`
  - a server: check the `dataset_id` passed to `show_anndata()`,
    `serve_anndata()`, or `cellucid serve`, else the store's
    `uns['cellucid']['id']`; `server_url + "/_cellucid/datasets"` shows the one
    in use

**Fix**
- If you intended to keep identity: write the id into the file (`prepare(..., dataset_id=...)`, or `adata.uns["cellucid"]["id"]`), or pass the same `dataset_id=...` to the server.
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
- Sharing datasets via a catalog: {doc}`02_local_demo_tutorial`
- Publishing a complete custom catalog: {doc}`11_custom_dataset_repository`
- What a file and a catalog contain: {doc}`07_file_format_expectations`
- When loading fails: {doc}`08_troubleshooting_data_loading`
- How the writers set the identity: {doc}`../../python_package/b_concepts_mental_models/04_dataset_identity_and_reproducibility` (Python) and {doc}`../../r_package/b_concepts_mental_models/03_dataset_identity_and_reproducibility` (R)
