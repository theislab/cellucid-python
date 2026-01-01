![Cellucid](https://raw.githubusercontent.com/theislab/cellucid-python/main/cellucid-logo.svg)

[![PyPI version](https://img.shields.io/pypi/v/cellucid.svg)](https://pypi.org/project/cellucid/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/cellucid.svg)](https://anaconda.org/conda-forge/cellucid)
[![Python versions](https://img.shields.io/pypi/pyversions/cellucid.svg)](https://pypi.org/project/cellucid/)
[![Documentation Status](https://readthedocs.org/projects/cellucid/badge/?version=latest)](https://cellucid.readthedocs.io/en/latest/)
[![License](https://img.shields.io/pypi/l/cellucid.svg)](https://pypi.org/project/cellucid/)

# Cellucid

**See every cell. Query any gene. Fly through millions. Interactive, GPU-accelerated single-cell visualization in the browser.**

Cellucid is a browser-first viewer for exploring large single-cell datasets in real time: fly through 2D/3D embeddings (UMAP/t-SNE/PCA), color by genes or metadata, filter and compare populations, and share reproducible views with collaborators.

## Highlights

- **GPU-accelerated WebGL rendering** with adaptive detail for millions of cells
- **AnnData-first workflow** for the Scanpy ecosystem (`.h5ad` and Zarr supported)
- **Shareable exports** you can host locally, on GitHub, or behind a server
- **Genes + metadata overlays** optimized for interactive querying
- **Connectivity + dynamics**: KNN edges and animated vector fields (RNA velocity / drift)
- **Collaboration**: community annotation voting with optional GitHub sync
- **Publication export**: SVG (vector) and high-DPI PNG figures
- **Works everywhere**: web app, local/remote server, and Jupyter notebooks

## Install

```bash
pip install cellucid
```

## Quickstart

Try the web app (no setup):

1. Open https://cellucid.com
2. Load a pre-exported folder, `.h5ad`, or `.zarr`

Or visualize an AnnData from Python/Jupyter:

```python
from cellucid import show_anndata

show_anndata(adata)  # or: show_anndata("dataset.h5ad")
```

## Links

- Web app: https://cellucid.com
- Documentation: https://cellucid.readthedocs.io
- Community annotation voting: https://cellucid.readthedocs.io/en/latest/user_guide/web_app/j_community_annotation/index.html
- Source: https://github.com/theislab/cellucid-python
- Web viewer: https://github.com/theislab/cellucid
- Annotation template: https://github.com/theislab/cellucid-annotation
- Issues: https://github.com/theislab/cellucid-python/issues

## Community

- Contributing: [CONTRIBUTING.md](CONTRIBUTING.md)
- Code of Conduct: [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)
- Security: [SECURITY.md](SECURITY.md)
- Support: [SUPPORT.md](SUPPORT.md)
- Citation: [CITATION.cff](CITATION.cff)

## License

BSD-3-Clause
