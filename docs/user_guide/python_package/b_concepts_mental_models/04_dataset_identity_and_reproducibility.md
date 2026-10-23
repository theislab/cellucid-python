# Dataset identity and reproducibility

**Audience:** everyone (especially computational users and collaborators)  
**Time:** 15–30 minutes  
**What you’ll learn:**
- What “dataset identity” means in Cellucid (and why it matters)
- Where a dataset's id and name come from (`uns['cellucid']`, or the server's arguments)
- How session bundles decide whether state is safe to restore/apply
- A practical versioning strategy for “paper-ready” datasets

---

## Mental model (one sentence)

If you want sessions, share links, and collaborative annotation to work reliably, treat your dataset file as an **immutable, versioned artifact** with a **stable dataset ID**.

---

## What is “dataset identity”?

In Cellucid, “dataset identity” is the answer to:

> “Is the dataset loaded right now the *same* dataset that this session / link / annotation was created for?”

Cellucid needs a dataset identity for:
- **Session restore** (web app): restore against the exact matching dataset or
  reject the bundle.
- **Session → AnnData application** (Python): mutate only the exact matching
  AnnData or raise before mutation.
- **Dataset catalogs**: how the viewer lists and labels datasets published together.
- **Collaboration**: if you send someone a session bundle, you need a stable reference to the dataset it belongs to.

---

## Where the identity lives: `uns['cellucid']`

A Cellucid-ready AnnData carries its identity in the file itself:

- `uns['cellucid']['id']`: the dataset ID (stable, machine-friendly)
- `uns['cellucid']['name']`: the human-readable name (also `uns['title']`)
- optional `description`, `citation` and `source`, and the launch `defaults`

```{note}
The id is a **metadata contract**, not a cryptographic hash of the content.
If you need tamper-proof guarantees, add your own checksums or version control
around the file.
```

---

## How the dataset ID is chosen (important)

### Writing it into the file (`prepare(...)`)

`prepare(...)` requires both `dataset_name` (human-friendly) and `dataset_id`
(stable, machine-friendly). It does not derive either value or rewrite an
unsafe ID.

Example:

```python
prepare(
    "pbmc_10k_v1.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    dataset_id="pbmc_10k_v1",
    dataset_name="PBMC 10k (Seurat v3 tutorial)",
    dataset_description="PBMC demo dataset used for tutorial figures",
    source_name="Seurat PBMC tutorial",
    source_url="https://satijalab.org/seurat/articles/pbmc3k_tutorial.html",
)
```

### Naming it when you open it (`show_anndata(...)`, `cellucid serve`)

A file that declares `uns['cellucid']['id']` and `['name']` is opened under
them. Any other AnnData needs a stable name and ID, passed explicitly:

```python
from cellucid import show_anndata
viewer = show_anndata("data.h5ad", dataset_name="My study (raw)", dataset_id="my_study_raw_v1")
```

An explicit `dataset_name` / `dataset_id` overrides what the file declares.
Data that declares neither (and is given neither) is refused before the server
starts:

```text
dataset_name and dataset_id are required: this AnnData declares no Cellucid identity of its own (uns['cellucid']['id'] and ['name']). Name the dataset, for example dataset_name="My study", dataset_id="my-study-v1" (--dataset-name/--dataset-id on the command line).
```

Sessions and community annotations are keyed by this dataset ID.

---

## How sessions decide “same dataset or not”

A `.cellucid-session` bundle contains a lightweight **dataset fingerprint** in
its manifest (not the full dataset). In the current session format (format 2) it has exactly six required fields: three scalars — `datasetId`,
`cellCount`, `varCount` — the dataset's display `name` (or null; used in
messages only, never compared), `cellOrder`, recording *which ordering* of the
cells the session was saved against, and `dataMapping`, recording *which arrays*
it was showing. Beside the fingerprint, the manifest's `requires` lists every
obs column and gene the session uses. It does not record how the dataset was opened, so a session
applies to the same dataset whether it was opened from a sample, a link, a
server, a notebook or a local folder.

In Python, `apply_cellucid_session_to_anndata(...)` uses this fingerprint as a safety guard:
- compares `datasetId` to the required `expected_dataset_id`,
- compares `cellCount` to `adata.n_obs`,
- checks that every key `dataMapping` names exists on `adata` (its `obsm`
  embeddings and their columns, `obsp` graph, layer or `raw`, latent space, and
  gene-name column),
- requires every obs column and gene the manifest's `requires` lists, naming
  all that are missing,
- compares `varCount` to the gene count of the axis that mapping reads:
  `adata.raw.n_vars` for `expression: "raw"`, `0` for a mapping without
  expression, `adata.n_vars` otherwise.

Any mismatch raises `ValueError` before AnnData is mutated. Application is
strict and atomic, and a bundle without `cellOrder` or `dataMapping` is refused.

### `cellOrder` is shape-checked here, its digest is verified in the viewer

`cellOrder` is exactly `embedding`, `dimension` (1, 2, or 3), and a
16-character lowercase hexadecimal `digest` of the coordinates of that
embedding. `embedding` is `{"key", "columns"}`: an `obsm` key with the 0-based
columns that were read, one per dimension. It closes a hole the three scalars
cannot: sessions store selections as row numbers, and a dataset rewritten at
the same id from re-sorted input matches all three scalars while every row
number now denotes a different cell.

`dataMapping` records which stored arrays the session was showing: `kind`
(`"anndata"`), `embeddings` (`enabled`, the embeddings that are turned on, and
`defaults`, the `1d`/`2d`/`3d` default of each dimension), `knn`, `geneNames`
(`"index"` or `"column:<name>"`), `expression` (`"X"`, `"raw"`,
`"layer:<name>"`, or null), `latent`, and `vectorFields`.

Both readers require both records. Python validates `cellOrder`'s shape and
never its digest: the digest is taken over normalized coordinates the viewer
drew, and `apply_cellucid_session_to_anndata(...)` targets an `AnnData` whose
row order the caller controls, so there is nothing here to check it against.
Enforcing cell-order identity is the viewer's job.

The consequence runs one way: **a bundle Python accepts may still be refused by
the web app**, because only the viewer compares the digest. Do not read that as a disagreement between the two readers.

```{important}
These guards prevent the most dangerous failure mode: silently applying highlights/labels to the wrong cells.
They do not guarantee perfect identity (e.g., two datasets can have the same shape), and on the Python side the row-order assumption is yours to keep — see Rule 4 below.
Use stable IDs + versioning to make identity meaningful.
```

---

## A versioning strategy that works (practical and boring — on purpose)

### Rule 1: treat published files as immutable

Once you share or publish a dataset file, do not modify it in place.
If you need to change anything, write a new file and bump the version.

Good:
- `pbmc_10k_v1.h5ad`
- `pbmc_10k_v2.h5ad`

Risky:
- overwriting `pbmc_10k.h5ad` while reusing the same dataset ID.

### Rule 2: version the dataset ID (not just the folder)

Make `dataset_id` reflect the version:
- `he_cortex_v2026_01_01`
- `hlca_lung_atlas_v3`
- `projectX_qcfiltered_v1`

### Rule 3: record provenance in `dataset_description` and/or `source_*`

If you want paper-grade reproducibility, include:
- which pipeline produced the AnnData,
- which commit/tag produced the file,
- which parameters (neighbors/umap) were used,
- any filtering decisions (QC thresholds).

Even a short description helps collaborators interpret sessions correctly.

### Rule 4: keep cell order stable across “same dataset”

Many Cellucid interactions use **row index** as identity.

If you reorder cells, previously saved sessions/highlights can map to different biological cells even if the shape is unchanged.

The web app catches that for a rewritten dataset, because `cellOrder` digests
the coordinates it drew. **Python does not**, and cannot: applying a
session to an `AnnData` you reordered in memory is checked only against
`n_obs`, `n_vars`, and `datasetId`, all of which a permutation leaves intact.
Keeping row order stable remains your responsibility on this side.

If you must reorder/subset:
- treat it as a new dataset version (new `dataset_id`),
- or do not reuse sessions across those changes.

---

## How this relates to collaboration tools

### Sessions (shareable state)

If you send a collaborator:
- a `.cellucid-session` file **without** the matching dataset file (or matching hosted dataset),
they may not be able to restore the state meaningfully.

### Community annotation (`cellucid-annotation`)

Community annotation workflows expect:
- stable dataset identity, so annotations can be applied consistently across contributors.

If you plan to use `cellucid-annotation`, decide your dataset ID/versioning early and do not change it casually.

---

## Edge cases (common identity mistakes)

### “I rewrote the same file with `force=True`”

`prepare(..., force=True)` replaces the whole file at once, so it is never a
mix of two runs. But anyone holding a session or annotation for the old content
now restores it onto the new one: change `dataset_id` whenever the cells
change.

### “Two datasets have the same `dataset_id`”

This is a recipe for confusion:
- sessions may “match” when they shouldn’t,
- collaborators may restore the wrong state onto the wrong dataset.

Fix: choose unique IDs (include a version suffix).

### “Does my AnnDataViewer publish local file paths?”

No. The server publishes the dataset ID and name, and serves the store under a
route named by the dataset ID (`_cellucid/anndata/<dataset_id>.h5ad`); the local
file name is never part of it. The terminal output still prints the resolved
path, so be mindful when sharing logs or screenshots (privacy/security covered
in {doc}`06_privacy_security_and_offline_vs_online`).

---

## Troubleshooting

### Symptom: “Session restore says dataset mismatch”

Likely causes:
- different number of cells/genes,
- different dataset ID,
- you are loading a different file than you think,
- the dataset was rewritten at the same ID from re-sorted input, so the
  counts match but the cell order does not,
- the session names an embedding, graph, layer, or gene-name column the dataset
  does not have,
- the session file predates the current fingerprint (no `cellOrder` with an
  `embedding`, or no `dataMapping`), so it is refused outright.

How to confirm:
- inspect `uns['cellucid']` and the shape of the file you loaded,
- compare counts (`n_obs`, `n_vars`) to the session’s fingerprint (web app warning or Python logs).

Fix:
- load the correct dataset file,
- or create a fresh session from the dataset generation you intend to use.

The current web and Python readers do not accept a layout-only partial
restore. Identity or format mismatch rejects the complete operation.

### Symptom: “Applying a session raises a dataset fingerprint mismatch”

Likely causes:
- `expected_dataset_id` does not match the bundle's `datasetId`,
- cell or variable counts differ,
- the bundle's `dataMapping` names an `obsm`, `obsp`, layer, `raw`, or `var`
  column the target `AnnData` does not have,
- `varCount` disagrees with the gene axis the mapping reads,
- the bundle's `cellOrder` or `dataMapping` record is missing or malformed.

A different row order is *not* on that list. Python compares counts and the
dataset ID, and a permutation changes neither, so it will apply the session
without complaint. Only apply to the same dataset generation and row order.

Fix:
- pass the exact ID used to create the matching viewer as
  `expected_dataset_id=...`,
- apply only to the same dataset generation and row order, or regenerate the
  session from the dataset version you intend to use.

---

## Next steps

- Capture state into Python and apply safely: {doc}`05_sessions_to_anndata_bridge`
- Understand what persists where: {doc}`03_state_persistence_and_scope`
