![Cellucid](https://raw.githubusercontent.com/theislab/cellucid-python/main/cellucid-logo.svg)

[![PyPI version](https://img.shields.io/pypi/v/cellucid.svg)](https://pypi.org/project/cellucid/)
[![Python versions](https://img.shields.io/pypi/pyversions/cellucid.svg)](https://pypi.org/project/cellucid/)
[![Documentation Status](https://readthedocs.org/projects/cellucid/badge/?version=latest)](https://cellucid.readthedocs.io/en/latest/)
[![License](https://img.shields.io/pypi/l/cellucid.svg)](https://pypi.org/project/cellucid/)

# Cellucid

**See every cell. Query any gene. Explore in 2D and 3D. Interactive, GPU-accelerated single-cell visualization in the browser.**

Cellucid is a browser-first viewer for exploring large single-cell datasets in real time: fly through 2D/3D embeddings (UMAP/t-SNE/PCA), color by genes or metadata, filter and compare populations, and share reproducible views with collaborators.

## Highlights

- **GPU-accelerated WebGL rendering** with level-of-detail culling; a built-in
  synthetic benchmark measures what your own GPU sustains
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

1. Open [cellucid.com](https://www.cellucid.com)
2. Load a prepared folder, one `.h5ad`, or one `.zarr.zip` archive

Or visualize an AnnData from Python/Jupyter:

```python
from cellucid import show_anndata

show_anndata(
    adata,
    dataset_name="My dataset",
    dataset_id="my-dataset",
)
```

## Documentation and ecosystem

- [Python package guide](https://cellucid.readthedocs.io/en/latest/user_guide/python_package/index.html)
- [Complete Cellucid documentation](https://cellucid.readthedocs.io/en/latest/)
- [Live web application](https://www.cellucid.com) and
  [web viewer source](https://github.com/theislab/cellucid)
- [R package](https://github.com/theislab/cellucid-r)
- [Official public demo datasets](https://github.com/theislab/cellucid-datasets)
- [Three custom dataset repository examples](https://github.com/theislab/cellucid-demo-custom-datasets)
- [Community annotation repository](https://github.com/theislab/cellucid-annotation)
  and [annotation guide](https://cellucid.readthedocs.io/en/latest/user_guide/web_app/j_community_annotation/index.html)
- [Python package issues](https://github.com/theislab/cellucid-python/issues)

## Community

- [Contributing](https://github.com/theislab/cellucid-python/blob/main/CONTRIBUTING.md)
- [Code of Conduct](https://github.com/theislab/cellucid-python/blob/main/CODE_OF_CONDUCT.md)
- [Security](https://github.com/theislab/cellucid-python/blob/main/SECURITY.md)
- [Support](https://github.com/theislab/cellucid-python/blob/main/SUPPORT.md)
- [Citation](https://github.com/theislab/cellucid-python/blob/main/CITATION.cff)

## License

BSD-3-Clause
