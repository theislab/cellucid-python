# Cellucid and other single-cell viewers

**Audience:** anyone choosing an interactive viewer for single-cell data  
**Time:** 3 minutes

Several good tools let you explore single-cell data interactively. They were
built for different jobs, and the right one depends on yours. This page says
what each project describes itself as, from its own documentation, and what
Cellucid was built for, so you can choose. It does not rank them, and it states
nothing about another tool that its own documentation does not.

## The tools

**CZ CELLxGENE.** CZ CELLxGENE Discover is a platform of curated, interoperable
single-cell data from the Chan Zuckerberg Initiative; each dataset page opens
in its Explorer ([Discover paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11701654/)).
For your own files, CELLxGENE Annotate is "an interactive data explorer for
single-cell datasets", installed with `pip install cellxgene` and started with
`cellxgene launch` on an AnnData file
([README](https://github.com/chanzuckerberg/cellxgene)).

**UCSC Cell Browser.** "A fast, lightweight viewer for single-cell data"; its
public site hosts datasets that labs submit
([documentation](https://cellbrowser.readthedocs.io/en/master/)). To build your
own, `cbBuild` converts an expression matrix and metadata into JSON and binary
files that any web server can host, with importers for Scanpy and Seurat
objects ([setup guide](https://cellbrowser.readthedocs.io/en/master/basic_usage.html)).

**Vitessce.** A framework for "integrative visualization of multimodal and
spatially resolved single-cell data", with JavaScript, Python and R packages
([vitessce.io](https://vitessce.io/)). It reads AnnData-Zarr, MuData-Zarr,
SpatialData, OME-TIFF and OME-Zarr, among other formats
([data file types](https://vitessce.io/docs/data-file-types/)).

## What Cellucid was built for

- **Very large datasets on ordinary hardware.** Points are drawn on the GPU
  with WebGL2, with optional level of detail and frustum culling; see
  {doc}`../n_benchmarking_performance/01_performance_mental_model`.
- **Any number of 1D, 2D and 3D embeddings**, one per view and any mix side
  by side, with orbit, planar and free-fly navigation per view; see
  {doc}`../c_core_interactions/05_dimension_switching_1d_2d_3d`.
- **No server needed to share.** An AnnData file, or a folder of them with a
  `datasets.json` catalog, is static: host it on GitHub and anyone opens it
  from a link
  ({doc}`../b_data_loading/11_custom_dataset_repository`). A server is there
  when you want one: `cellucid serve` in Python, or a notebook.
- **Cellucid-ready AnnData from Python and R** with the same layout:
  `cellucid.prepare` and `cellucid_prepare()`.
- **Community annotation**: several people vote on cell-type labels for one
  dataset through a GitHub repository; see {doc}`../j_community_annotation/index`.
- **In-browser analysis and figures**: compare selections, find marker genes,
  and export vector figures.

## Choosing

- **Browsing curated public atlases**: CELLxGENE Discover gives you its corpus
  ready to explore.
- **Spatial and multimodal data with images**: Vitessce is designed for it.
- **Publishing a lab's datasets as static files**: both the UCSC Cell Browser
  and Cellucid produce them. Cellucid opens AnnData files and Zarr folders from a GitHub repository in
  the hosted app, so you need no web server or site of your own.
- **Millions of cells in 2D or 3D, on a laptop or a phone, in your own browser
  with nothing uploaded**: that is what Cellucid is for.

Most of these tools read AnnData, so trying more than one on the same dataset
is cheap. To start with Cellucid, see
{doc}`06_how_to_visualize_single_cell_data_interactively`.
