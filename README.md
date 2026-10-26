![Cellucid](https://raw.githubusercontent.com/theislab/cellucid-python/main/cellucid-logo.svg)

[![PyPI version](https://img.shields.io/pypi/v/cellucid.svg)](https://pypi.org/project/cellucid/)
[![Python versions](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue.svg)](https://docs.cellucid.com/en/latest/user_guide/python_package/installation.html)
[![Documentation Status](https://readthedocs.org/projects/cellucid/badge/?version=latest)](https://docs.cellucid.com/en/latest/)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# Cellucid

**See every cell. Query any gene. Explore in 2D and 3D. Interactive, GPU-accelerated single-cell visualization in the browser.**

Cellucid is a browser-first viewer for exploring large single-cell datasets in real time: fly through 2D/3D embeddings (UMAP/t-SNE/PCA), color by genes or metadata, filter and compare populations, and share reproducible views with collaborators.

## Highlights

- **GPU-accelerated WebGL rendering** with level-of-detail culling; a built-in
  synthetic benchmark measures what your own GPU sustains
- **AnnData-first workflow** for the Scanpy ecosystem (`.h5ad` and Zarr), plus MuData (`.h5mu`, MuData Zarr), Loom, 10x Cell Ranger `.h5` and `outs/`, MatrixMarket folders and R objects (`.rds`, `.rda`), read by the browser itself
- **Cellucid-ready AnnData** (`prepare()`): one standard `.h5ad`, `.zarr` or `.zarr.zip`, laid out for fast browser reads, that you can host locally, on GitHub, on any web host, or behind a server
- **Genes + metadata overlays** optimized for interactive querying
- **Connectivity + dynamics**: KNN edges and animated vector fields (RNA velocity / drift)
- **Collaboration**: community annotation voting with optional GitHub sync
- **Publication export**: SVG (vector) and high-DPI PNG figures
- **Works everywhere**: web app, local/remote server, and Jupyter notebooks

## Install

Cellucid requires Python 3.11, 3.12, 3.13, or 3.14. Use the dedicated
[Installation guide](https://docs.cellucid.com/en/latest/user_guide/python_package/installation.html)
to choose PyPI or GitHub based on current registry availability. To install the
active source directly from the official GitHub repository:

```bash
pip install "cellucid @ git+https://github.com/theislab/cellucid-python"
```

From PyPI:

```bash
pip install cellucid
```

The [PyPI package page](https://pypi.org/project/cellucid/) is authoritative for
registry availability. Until it lists the active version, `pip install cellucid`
resolves only to the versions PyPI already holds, and those are earlier
generations whose API differs from the one this README and the documentation
describe.

## Quickstart

Try the web app (no setup):

1. Open [cellucid.com](https://www.cellucid.com)
2. Open a dataset file or folder: an `.h5ad`, a Zarr store (folder or
   `.zarr.zip`), an `.h5mu`, a `.loom`, a 10x `.h5` or `outs/` folder, a
   MatrixMarket folder, or an `.rds`/`.rda`

Or visualize an AnnData from Python/Jupyter:

```python
from cellucid import show_anndata

show_anndata(
    adata,
    dataset_name="My dataset",
    dataset_id="my-dataset",
)
```

An `adata` straight out of Scanpy works as written, and so does a `MuData`.
The viewer reads the data itself, as each view needs it: `cellucid serve`
hands it any file or folder in place (`.h5ad`, a Zarr store, `.h5mu`, `.loom`,
a 10x `.h5` or `outs/` folder, a MatrixMarket folder, an `.rds`/`.rda`), and
an in-memory object is written once to a temporary `.h5ad` (or `.h5mu`).
Python and the notebook viewer talk over one WebSocket to the local server,
so hooks and commands work wherever the notebook's browser can reach that
server: Jupyter and JupyterLab on the same machine, an SSH-forwarded port, or
a proxy that passes WebSockets (VS Code's notebook webview and Google Colab
have not been tested).

To share or host a dataset, write a Cellucid-ready AnnData once. It is still a
standard AnnData, with the expression matrix gene-major and the dataset's name,
citation and launch defaults in `uns['cellucid']`, so the browser reads one
gene with one round trip and the file names itself (on a 562k-cell dataset
over a simulated 80 ms link: opened in 0.25 s, the first gene in 0.09 s; see
the specification's measured choices):

```python
from cellucid import prepare

prepare(
    "my_dataset.h5ad",
    format="h5ad",
    data=adata,
    dataset_id="my-dataset",
    dataset_name="My dataset",
)
```

`data=` takes an AnnData or a MuData whole (a MuData is written as `.h5mu`
with `format="h5mu"`, or as MuData Zarr with `"zarr"`/`"zarr.zip"`; `pip
install "cellucid[mudata]"`); any
component passed explicitly (`embeddings=`, `obs=`, `gene_expression=`, …)
replaces that one component.

`format` is required: `"zarr"` (a folder of small files) for GitHub, which
refuses files over 100 MB and serves many small files at once over HTTP/2;
`"h5ad"` (or `"zarr.zip"`, the Zarr store in one file) for everything else —
single-file hosts such as Zenodo, Figshare, Hugging Face, S3 or a lab server,
local files, `cellucid serve` and Jupyter. The documentation's "Which format?"
guide has the measurements.

Then `cellucid serve my_dataset.h5ad`, `show_anndata("my_dataset.h5ad")`, or
open the file in the web app.

## Documentation and ecosystem

- [Python package guide](https://docs.cellucid.com/en/latest/user_guide/python_package/index.html)
- [Complete Cellucid documentation](https://docs.cellucid.com/en/latest/)
- [Live web application](https://www.cellucid.com) and
  [web viewer source](https://github.com/theislab/cellucid)
- [R package](https://github.com/theislab/cellucid-r)
- [Official public demo datasets](https://github.com/theislab/cellucid-datasets)
- [Custom dataset repository examples](https://github.com/theislab/cellucid-demo-custom-datasets)
- [Community annotation repository](https://github.com/theislab/cellucid-annotation)
  and [annotation guide](https://docs.cellucid.com/en/latest/user_guide/web_app/j_community_annotation/index.html)
- [Python package issues](https://github.com/theislab/cellucid-python/issues)

## Citation

<!-- BEGIN cellucid-citation (generated by cellucid-private/scripts/sync_citation.py from CITATION.cff; do not edit) source-sha256=2486463f2b7c8cc656ae41e2df0261c1eb79780c4d635be147215f409b89aa2a body-sha256=5186fc40b9eddfb904fab52bee5362bf27a8d6c3713f9b67b1defa1c68c54d3f -->
> **Placeholder citation.** The Cellucid preprint is in preparation for bioRxiv and has no DOI yet; the reference below is a placeholder and will be replaced when the preprint is published.

If you use Cellucid, please cite:

Inecik, K. PLACEHOLDER: Cellucid preprint, title to be announced (forthcoming on bioRxiv). *bioRxiv* (in preparation).

BibTeX and other formats: [Citing Cellucid](https://docs.cellucid.com/en/latest/citing.html).
<!-- END cellucid-citation -->

## Community

- [Contributing](https://github.com/theislab/cellucid-python/blob/main/CONTRIBUTING.md)
- [Code of Conduct](https://github.com/theislab/cellucid-python/blob/main/CODE_OF_CONDUCT.md)
- [Security](https://github.com/theislab/cellucid-python/blob/main/SECURITY.md)
- [Support](https://github.com/theislab/cellucid-python/blob/main/SUPPORT.md)
- [Citation](https://github.com/theislab/cellucid-python/blob/main/CITATION.cff)

## License

BSD-3-Clause
