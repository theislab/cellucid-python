# Sessions → AnnData (No-Download Bridge)

This page explains the “session bundle → AnnData mutation” workflow in `cellucid-python`.

Goal:
- Interact in the viewer.
- Pull the current session state into Python **as an object** (no manual browser download).
- Apply the meaningful parts onto an `anndata.AnnData` safely and efficiently.

---

## The two state layers (mental model)

Cellucid exposes two complementary layers of notebook state:

1) **`viewer.state` (live snapshot)**  
   A small “latest event” snapshot updated continuously from UI events (selection/hover/click/ready).

2) **Session bundle (`.cellucid-session`) (durable artifact)**  
   A file-format snapshot of meaningful UI state (highlights, user-defined fields, overlay registries, layout state).
   In Jupyter you can request it as a `CellucidSessionBundle` object.

If you want to mutate an `AnnData` reproducibly, use the **session bundle**.

---

## Workflow: Jupyter → bundle → AnnData

```python
from cellucid import show_anndata

viewer = show_anndata("data.h5ad", height=650)
viewer.wait_for_ready(timeout=60)

# Pull the current viewer session bundle into Python (no download UI step)
bundle = viewer.get_session_bundle(timeout=60)

# Apply to AnnData (non-destructive by default)
adata2 = bundle.apply_to_anndata(adata, inplace=False)
```

Convenience one-liner (requests bundle + applies + cleans up the temp bundle file):

```python
adata2 = viewer.apply_session_to_anndata(adata, inplace=False)
```

---

## What gets applied to AnnData (current behavior)

`apply_cellucid_session_to_anndata(...)` currently materializes:

- **Highlights → `adata.obs`**  
  One boolean `.obs` column per highlight group (`cellucid_highlight__<groupId>` by default).

- **User-defined categorical fields → `adata.obs` / `adata.var`**  
  Categorical columns decoded from `user-defined/codes/<fieldId>` (categories come from `core/field-overlays`).

It also stores metadata under:
- `adata.uns["cellucid"]["session"]` (manifest, dataset fingerprint, decoded JSON chunks)

---

## Safety: dataset mismatch + index-based identity

### Dataset mismatch policy

Sessions are treated as untrusted input. Before applying dataset-dependent chunks, Python checks the session fingerprint against the AnnData shape.

You control behavior with `dataset_mismatch`:
- `"error"`: raise on mismatch
- `"warn_skip"`: warn and skip dataset-dependent chunks
- `"skip"`: silently skip dataset-dependent chunks

### Critical constraint: cell identity is index-based (today)

Session highlights and codes are aligned to **cell indices** (row position).

Only apply a session to an `AnnData` whose row order matches the dataset that produced the session.

Recommended extra guard:
- pass `expected_dataset_id=...` when you have a stable dataset identity string

Roadmap direction:
- move toward stable per-cell identifiers (e.g. mapping by `obs_names`) so sessions can be applied safely after reordering/subsetting.

---

## Advanced: apply from a saved file

```python
from cellucid import CellucidSessionBundle

bundle = CellucidSessionBundle("my.cellucid-session")
adata2 = bundle.apply_to_anndata(adata, inplace=False)
```

