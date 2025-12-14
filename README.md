<p>
  <img src="https://raw.githubusercontent.com/theislab/cellucid-python/main/cellucid-logo.svg" alt="Cellucid logo" width="360">
</p>

[![PyPI version](https://img.shields.io/pypi/v/cellucid.svg)](https://pypi.org/project/cellucid/)
[![Python versions](https://img.shields.io/pypi/pyversions/cellucid.svg)](https://pypi.org/project/cellucid/)
[![Documentation Status](https://readthedocs.org/projects/cellucid/badge/?version=latest)](https://cellucid.readthedocs.io/en/latest/)

# Cellucid

Python package for visualizing single-cell data with the [Cellucid WebGL viewer](https://cellucid.com). Cellucid renders millions of cells in real-time 3D, enabling interactive exploration of UMAP/tSNE embeddings with gene expression overlays, filtering, and KNN connectivity visualization.

## Features

- **Interactive Single-Cell Data Visualization** of UMAP/tSNE embeddings
- **Gene expression overlays** with fast queries
- **Cell metadata filtering** and categorical coloring
- **KNN connectivity** edge visualization
- **Multi-dimensional** support (1D timelines, 2D, 3D)
- **Jupyter integration** with bidirectional communication
- **Scales to millions** of cells with adaptive LOD

## Installation

```bash
pip install cellucid
```

All features (Jupyter, AnnData, server) are included in the standard install.

## Quick Start

### Option 1: Direct AnnData Visualization (No Export)

The fastest way to get started - visualize your AnnData directly:

```python
from cellucid import show_anndata

# In-memory AnnData
show_anndata(adata)

# From h5ad file (lazy loading)
show_anndata("/path/to/data.h5ad")

# From zarr store (lazy loading)
show_anndata("/path/to/data.zarr")
```

### Option 2: Pre-export for Best Performance

For production use or sharing, export to optimized binary format:

```python
from cellucid import prepare, show

# Export once (quantized + compressed)
prepare(adata, output_dir="./my_export", compress=True)

# View anytime (fastest loading)
show("./my_export")
```

## AnnData Requirements

Your AnnData needs UMAP coordinates in `obsm`:

```python
# Required: one of these
adata.obsm['X_umap_3d']  # shape (n_cells, 3) - recommended
adata.obsm['X_umap_2d']  # shape (n_cells, 2)
adata.obsm['X_umap']     # shape (n_cells, 2 or 3)

# Optional
adata.obs                # Cell metadata (categorical/continuous)
adata.X                  # Gene expression (dense or sparse)
adata.obsp['connectivities']  # KNN graph edges
```

Computing 3D UMAP with scanpy:

```python
import scanpy as sc
sc.pp.neighbors(adata)
sc.tl.umap(adata, n_components=3)
adata.obsm['X_umap_3d'] = adata.obsm['X_umap']
```

## All 14 Loading Options

Cellucid supports 6 deployment modes, each with support for pre-exported binary data, h5ad files, and zarr stores:

| # | Method | Exported | h5ad | zarr | Python | Lazy Load | Performance |
|---|--------|----------|------|------|--------|-----------|-------------|
| 1 | Local Demo (GitHub) | ✅ | - | - | No* | Yes | Best |
| 2 | Remote Demo (GitHub) | ✅ | - | - | No* | Yes | Best |
| 3 | Browser File Picker | ✅ | - | - | No | Yes | Best |
| 4 | Browser File Picker | - | ✅ | - | No | **No** | Slower |
| 5 | Browser File Picker | - | - | ✅ | No | **No** | Slower |
| 6 | Server CLI | ✅ | - | - | Yes | Yes | Best |
| 7 | Server CLI | - | ✅ | ✅ | Yes | Yes | Good |
| 8 | Python serve() | ✅ | - | - | Yes | Yes | Best |
| 9 | Python serve_anndata() | - | ✅ | ✅ | Yes | Yes | Good |
| 10 | Jupyter show() | ✅ | - | - | Yes | Yes | Best |
| 11 | Jupyter show_anndata() | - | ✅ | ✅ | Yes | Yes | Good |

\* Python required for initial export, not for viewing

**Summary by method:**
| Method | Exported | h5ad | zarr | Total |
|--------|----------|------|------|-------|
| Local/Remote Demo | ✅ | - | - | 2 |
| Browser File Picker | ✅ | ✅ | ✅ | 3 |
| Server CLI | ✅ | ✅ | ✅ | 3 |
| Python serve | ✅ | ✅ | ✅ | 3 |
| Jupyter | ✅ | ✅ | ✅ | 3 |
| **Total** | | | | **14** |

### Key Notes:
- **Browser h5ad/zarr**: Entire file loaded into memory - no lazy loading due to JavaScript limitations
- **Python h5ad/zarr modes**: True lazy loading via AnnData backed mode (h5ad) or zarr's native chunked access
- **Pre-exported data**: Always fastest - use for production and sharing
- **zarr stores**: Can be a directory (.zarr) or a file - the Python server auto-detects the format

## Deployment Modes

| Mode | Best For | Command/Function |
|------|----------|------------------|
| **Jupyter** | Interactive analysis | `show_anndata(adata)` |
| **Local server** | Development | `cellucid serve ./export` |
| **Remote + SSH** | Team access | `cellucid serve-anndata data.h5ad` |
| **Browser** | Quick preview | [cellucid.com](https://cellucid.com) file picker |
| **Static hosting** | Public sharing | GitHub Pages |

### CLI Commands

```bash
# Serve h5ad file directly (lazy loading, no export needed)
cellucid serve-anndata /path/to/data.h5ad

# Serve zarr store directly (lazy loading, no export needed)
cellucid serve-anndata /path/to/data.zarr

# Serve pre-exported data (fastest)
cellucid serve /path/to/export

# With options
cellucid serve-anndata data.h5ad --port 9000 --no-browser

# Show version
cellucid --version
```

### Remote Server Access

```bash
# On remote server
cellucid serve-anndata /data/cells.h5ad --no-browser

# On local machine (SSH tunnel)
ssh -L 8765:localhost:8765 user@server
# Then open: https://cellucid.com?remote=http://localhost:8765
```

## Jupyter Features

```python
from cellucid import show_anndata

viewer = show_anndata(adata, height=600)

# Programmatic control
viewer.highlight_cells([1, 2, 3], color="#ff0000")
viewer.set_color_by("cell_type")
viewer.set_visibility([0, 1, 2], visible=False)

# Cleanup
viewer.stop()
```

## Performance Guide

| Method | Load Time | Gene Queries | Recommended For |
|--------|-----------|--------------|-----------------|
| Pre-exported (`show`) | Fast | Fast | Production, sharing |
| Direct AnnData (`show_anndata`) | Medium | Medium | Exploration |
| Browser file picker | Slower | Slower | Quick preview |

For datasets > 500k cells, pre-export with `prepare()` is recommended.

## Documentation

- **[Full documentation](https://cellucid.readthedocs.io)** - API reference and tutorials
- **[Web viewer](https://github.com/theislab/cellucid)** - JavaScript repository

## License

BSD-3-Clause
