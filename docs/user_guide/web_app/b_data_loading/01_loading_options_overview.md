# Data Loading Overview (All 14 Options)

This page is the **map** of how Cellucid can load data.

If you are new, start here. If you already know what you want, jump to the tutorial that matches your workflow:

- `02_local_demo_tutorial`: publish/share a dataset without running a server (GitHub-hosted exports)
- `03_browser_file_picker_tutorial`: load data from your computer using the browser file picker
- `04_server_tutorial`: run a local/remote Python server for large datasets
- `05_jupyter_tutorial`: embed Cellucid inside a notebook and interact programmatically

## At A Glance

**Audience**
- Wet lab / non-technical: choose a safe workflow and avoid common traps.
- Computational: understand formats, lazy loading, and performance tradeoffs.
- Power users: understand what happens under the hood so you can debug failures.

**Time**
- Fast path decision: ~5 minutes
- Full read (with edge cases): ~20–30 minutes

**What you’ll learn**
- What “exported”, “h5ad”, and “zarr” mean in Cellucid
- What “lazy loading” is and why it matters for large datasets
- How to choose between **File Picker**, **Server Mode**, **GitHub-hosted exports**, and **Jupyter**
- The top failure modes (and how to quickly diagnose them)

## The Three Data Formats (What You Actually Load)

Cellucid understands your data in three high-level ways:

1) **Pre-exported folder (recommended)**
   - Produced by `cellucid.prepare(...)` in Python.
   - Output is a directory containing binary files + `dataset_identity.json` + compact manifests.
   - **Fastest** in the browser and the most stable.
   - Can optionally include **vector fields** (e.g. RNA velocity) for the overlay.

2) **AnnData `.h5ad` file**
   - A single HDF5 file.
   - Can be loaded in several ways:
     - directly in the browser (works, but has hard performance limits)
     - via a Python server (recommended for large datasets)
     - inside Jupyter (recommended for analysis workflows)

3) **AnnData `.zarr` directory**
   - A directory-based format (chunked arrays).
   - Typically better than `.h5ad` for lazy loading.
   - Can be loaded:
     - directly in the browser (often feasible; still depends on dataset size)
     - via a Python server (best reliability)
     - inside Jupyter

**Rule of thumb:** if you care about smooth interaction and reproducibility, **export once** using `prepare()`.

For deeper context:
- why stable IDs matter: {doc}`06_dataset_identity_why_it_matters`
- what files/keys are expected: {doc}`07_folder_file_format_expectations_high_level_link_to_spec`
- when loading fails: {doc}`08_troubleshooting_data_loading`

## What “Lazy Loading” Means (In Plain Language)

When you load a dataset, Cellucid needs several kinds of data:

- **Embeddings** (UMAP coordinates): small-ish (n_cells × 2/3) → loaded early
- **Obs fields** (metadata for coloring/filtering): usually manageable → loaded early
- **Gene expression** (n_cells × n_genes): huge → ideally loaded **on demand**

**Lazy loading** means:
- Cellucid loads just enough to render and interact.
- When you search a gene (or request some data), it fetches only what it needs.

Why you care:
- For large datasets, loading everything up front can crash the browser or take minutes.
- Lazy loading keeps memory bounded and makes the UI responsive sooner.

Important nuance:
- Some formats/modes are truly lazy.
- Some are “lazy-ish” (metadata up front; some chunks later).
- Some are **not lazy** in practice (browser `.h5ad` requires loading the whole file).

## Vector fields (velocity/drift overlay): what you need for it to appear

Vector fields are optional per-cell vectors (e.g. RNA velocity, drift/displacement) that Cellucid can render as an animated overlay on top of your embedding.

Important rules:
- Vector fields are **dimension-specific**: you need 2D vectors for a 2D embedding, 3D vectors for 3D, etc.
- Vector fields must match the **same row order** as your points/cells.
- The overlay only shows fields available for the **current dimension**.

Where vector fields come from depends on your loading format:

- **Exports (`prepare`)**: `dataset_identity.json` contains a `vector_fields` block, and binary vectors live under `vectors/`.
- **AnnData (`.h5ad` / `.zarr` / in-memory)**: vector fields are discovered from `adata.obsm` keys like:
  - `velocity_umap_2d`, `velocity_umap_3d`
  - `T_fwd_umap_2d`

If you expect vector fields but don’t see the overlay toggle or dropdown:
- first check your data format/naming: {doc}`07_folder_file_format_expectations_high_level_link_to_spec`
- then check the overlay section: {doc}`../i_vector_field_velocity/index`

Internal reference (naming conventions):
- `cellucid/markdown/VECTOR_FIELD_OVERLAY_CONVENTIONS.md`

## Fast Path: Choose A Workflow (Decision Tree)

Pick the first row that matches you.

| Your situation | Recommended workflow | Why |
|---|---|---|
| “I just want to look at my data quickly, no Python.” | Browser File Picker | Zero setup; good for quick preview |
| “My dataset is big (hundreds of thousands to millions of cells).” | Server Mode or pre-export + Server Mode | True lazy loading; avoids browser memory limits |
| “I’m already working in a notebook.” | Jupyter (`show_anndata`, `show`) | Tight analysis loop; programmatic control |
| “I want to share a dataset publicly without running a server.” | GitHub-hosted exports | Shareable URL; no server |
| “I need the fastest possible web experience.” | Pre-exported folder | Best performance and stability |

If you’re unsure, start with **Server Mode** for `.h5ad`/`.zarr`, or **File Picker** for exported folders.

## The 14 Loading Options (Complete Matrix)

This is the canonical list used throughout the documentation.

**Legend**
- **Exported** = folder created by `cellucid.prepare()`
- **Lazy genes**: whether gene expression is fetched on demand (best) vs effectively “load all” (worst)
- **Vector fields**: supported in all loading options; the overlay appears only if your dataset includes vectors for the current dimension.

| # | Where you run things | How you point Cellucid to the data | Data format | Lazy genes | Best for |
|---:|---|---|---|---|---|
| 1 | Cellucid web app (demo mode) | Choose a built-in demo dataset | Exported | ✅ | Learning the UI with known-good data |
| 2 | Cellucid web app (public GitHub) | Connect to a public repo (or `?github=...`) | Exported | ✅ | Sharing a dataset publicly, no server |
| 3 | Cellucid web app | Browser **Folder** picker | Exported | ✅ | Fast local viewing of prepared exports |
| 4 | Cellucid web app | Browser **.h5ad** picker | `.h5ad` | ❌* | Quick preview of small `.h5ad` |
| 5 | Cellucid web app | Browser **.zarr** folder picker | `.zarr` | ✅† | Quick preview of `.zarr` without Python |
| 6 | Terminal CLI | `cellucid serve <export_dir>` | Exported | ✅ | Reliable viewing of large exports |
| 7 | Terminal CLI | `cellucid serve <data.h5ad>` | `.h5ad` | ✅ | Large `.h5ad` with true lazy loading |
| 8 | Terminal CLI | `cellucid serve <data.zarr>` | `.zarr` | ✅ | Large `.zarr` with true lazy loading |
| 9 | Python | `cellucid.serve(<export_dir>)` | Exported | ✅ | Scripting server startup |
| 10 | Python | `cellucid.serve_anndata(<data.h5ad>)` | `.h5ad` | ✅ | Scripting server startup |
| 11 | Python | `cellucid.serve_anndata(<data.zarr>)` | `.zarr` | ✅ | Scripting server startup |
| 12 | Jupyter | `cellucid.show(<export_dir>)` | Exported | ✅ | Notebook-based exploration of exports |
| 13 | Jupyter | `cellucid.show_anndata(<data.h5ad>)` | `.h5ad` | ✅ | Notebook-based exploration of `.h5ad` |
| 14 | Jupyter | `cellucid.show_anndata(<data.zarr or AnnData>)` | `.zarr` / in-memory | ✅ | Notebook-based exploration of `.zarr` or in-memory |

\* Browser `.h5ad` loading is **not truly lazy**: the whole file is loaded into browser memory before use.

† Browser `.zarr` loading can be **effectively lazy** for gene expression (chunked storage), but still needs to read metadata up front.

## Minimal Commands (Copy/Paste)

### Install

```bash
pip install cellucid
```

### 1) Export once (recommended)

```python
from cellucid import prepare

# prepare(...) writes an export folder that loads fast in the browser.
# You will typically pass:
# - obs (cell metadata)
# - var (gene metadata)
# - gene_expression (adata.X)
# - embeddings (e.g., X_umap_2d / X_umap_3d)
# - vector_fields (optional: velocity/drift overlays)
# See the prepare() documentation for full details.
```

### 2) Serve any data (auto-detected)

```bash
cellucid serve /path/to/data.h5ad
cellucid serve /path/to/data.zarr
cellucid serve /path/to/export_dir
```

### 3) Notebook quick start

```python
from cellucid import show_anndata

viewer = show_anndata("data.h5ad")
```

```python
# Optional: sanity-check your AnnData quickly
#
# This cell is intentionally lightweight; it's here so you can quickly
# confirm that your dataset has the minimum requirements before you pick a workflow.

from __future__ import annotations

from pathlib import Path

import numpy as np

try:
    import anndata as ad
except Exception as e:  # pragma: no cover
    raise RuntimeError(
        "anndata is required for this cell. Install with: pip install anndata"
    ) from e


def summarize_anndata(path: str | Path) -> dict:
    path = Path(path)
    adata = ad.read_h5ad(path, backed="r")  # doesn't load full X into memory
    # Vector fields (if present) are stored in obsm and are NOT prefixed with X_.
    # Common naming pattern: <field>_umap_<dim>d (e.g., velocity_umap_2d).
    import re

    vector_field_re = re.compile(r"^(?!X_).+_umap(?:_[123]d)?$")
    vector_field_obsm_keys = sorted(
        [k for k in adata.obsm.keys() if vector_field_re.match(k)]
    )
    summary = {
        "path": str(path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obsm_keys": list(adata.obsm.keys()),
        "vector_field_obsm_keys": vector_field_obsm_keys,
        "has_X": adata.X is not None,
        "X_type": type(adata.X).__name__,
    }
    return summary


# Example usage (uncomment):
# print(summarize_anndata("/path/to/data.h5ad"))
```

## Screenshot Placeholders (Add Later)

The next notebooks are extremely step-by-step. This overview page still benefits from a few “orientation” images.
For a capture checklist, see {doc}`09_screenshots`.

<!-- SCREENSHOT PLACEHOLDER
ID: data-loading-overview-loader-panel
Suggested filename: data_loading/02_overview-loader-panel.png
Where it appears: Data Loading → 01_loading_options_overview → Screenshot Placeholders
Capture:
  - UI location: left sidebar → Dataset Connections / Load Data area
  - State prerequisites: Cellucid open; no dataset loaded
  - Action to reach state: open Cellucid, then ensure you are in the default/empty state
Crop:
  - Include: left sidebar + top-left of the canvas (so readers can orient)
  - Exclude: personal bookmarks/extensions
Redact:
  - Remove: private dataset names and local paths
Annotations:
  - Callouts: highlight the 3 main paths: Local, Remote server, GitHub repo
Alt text:
  - Left sidebar with the data loading controls highlighted.
Caption:
  - Explain what this panel is for.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Cellucid data loading controls.
:width: 100%

The left sidebar is where you choose how to load data (local files, a server, or a public GitHub repo).
```

<!-- SCREENSHOT PLACEHOLDER
ID: data-loading-overview-success-state
Suggested filename: data_loading/03_overview-success-state.png
Where it appears: Data Loading → 01_loading_options_overview → Screenshot Placeholders
Capture:
  - UI location: main canvas + dataset name/summary area
  - State prerequisites: a dataset successfully loaded
  - Action to reach state: load any small demo dataset
Crop:
  - Include: dataset name + point count + one visible embedding
  - Exclude: unrelated browser chrome
Redact:
  - Remove: private dataset identifiers
Alt text:
  - Cellucid canvas showing a loaded embedding with a visible legend.
Caption:
  - Describe what “loaded successfully” looks like.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a successful dataset load.
:width: 100%

A successful load shows a populated embedding and a field/legend you can color by.
```

## Edge Cases (Read This If Something Feels “Weird”)

### Data edge cases
- **No UMAP / embedding keys**: your AnnData is missing required `obsm` entries (`X_umap`, `X_umap_2d`, or `X_umap_3d`).
- **NaN/Inf in embeddings**: points may disappear, render at the origin, or crash WebGL shaders.
- **Duplicate cell IDs**: selection/highlighting can behave unpredictably.
- **Huge categorical fields** (e.g., 100k unique categories): legends and color assignment become unusable.

### Scale edge cases
- **Browser `.h5ad`**: if the file is large, the browser may freeze or crash due to memory pressure.
- **Millions of cells + dense gene matrix**: even with lazy loading, gene queries can be slow or memory-heavy.

### Environment edge cases
- **Safari / older browsers**: directory picking and File System Access APIs may be limited.
- **Corporate environments**: content blockers can break `raw.githubusercontent.com` (GitHub-hosted exports).

## Troubleshooting (Quick Index)

This is a compact troubleshooting index. Each follow-up notebook contains a **much longer** troubleshooting section.

### Symptom: “Nothing loads / it spins forever”
**Likely causes**
- Large `.h5ad` opened directly in the browser
- CORS/network blocked (server mode or GitHub mode)
- Browser GPU/WebGL context lost due to memory pressure

**How to confirm**
- Open Developer Tools → Console → look for network errors / memory warnings
- Try loading a built-in demo dataset: if demos load, your environment is probably fine

**Fix**
- Prefer `cellucid serve data.h5ad` (server mode) for large `.h5ad`
- Prefer exports created by `prepare()` for best performance

---

### Symptom: “I don’t see any embedding / it says no UMAP”
**Likely causes**
- AnnData missing `obsm['X_umap']` (or 2d/3d variants)

**How to confirm**
- In Python: `print(adata.obsm.keys())`

**Fix**
- Compute UMAP and store it in `obsm` before viewing

---

### Symptom: “GitHub connect says datasets.json not found”
**Likely causes**
- The repo path points at the wrong folder
- You uploaded only a single dataset folder without a top-level `datasets.json`

**Fix**
- Create an exports root directory containing `datasets.json` and one or more dataset subfolders.
- See {doc}`02_local_demo_tutorial` for the exact layout and a helper function to generate `datasets.json`.

## Next Steps

Choose your path:

- Want to share a dataset publicly (no server)? → {doc}`02_local_demo_tutorial`
- Want to load from your own machine right now? → {doc}`03_browser_file_picker_tutorial`
- Working with a big `.h5ad`/`.zarr`? → {doc}`04_server_tutorial`
- Already in notebooks? → {doc}`05_jupyter_tutorial`
