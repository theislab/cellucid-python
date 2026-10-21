# How to visualize single-cell data interactively

**Audience:** anyone with a single-cell dataset — an `.h5ad` from Scanpy, a
Seurat or SingleCellExperiment object, or just a file someone sent you  
**Time:** 5 minutes to a first view  
**What you'll learn:** the three ways to open your own data in Cellucid, and
which one to pick

---

## Pick the route that matches what you have

| You have | Route | Install |
| --- | --- | --- |
| An `.h5ad` file, of any size | [Open it in the browser](#route-1-open-a-file-in-the-browser) | Nothing |
| An AnnData in Python (Scanpy), or an `.h5ad` you open often | [Serve it from Python](#route-2-serve-it-from-python) | `pip install cellucid` |
| A Seurat or SingleCellExperiment object in R | [Export it from R](#route-3-export-it-from-r) | the `cellucid` R package |

Every route ends in the same viewer, with millions of cells drawn on the GPU in
1D, 2D or 3D. Your data stays on your machine: nothing is uploaded.

## Route 1: open a file in the browser

1. Open [www.cellucid.com](https://www.cellucid.com) and choose **Choose a
   dataset** on the welcome card.
2. In the **Session** panel, under **Open data**, press **H5AD** and pick your
   file.

The file needs at least one embedding in `obsm`. The `X_umap` that
`sc.tl.umap()` writes is read as it is; see
{doc}`../b_data_loading/03_browser_file_picker_tutorial` for the exact rule and
how a large file is read.

## Route 2: serve it from Python

From a terminal, for an `.h5ad` file of any size:

```bash
pip install cellucid
cellucid serve my_dataset.h5ad --dataset-name "My dataset" --dataset-id my-dataset
```

Open the viewer URL it prints. The file is opened in backed mode and read on
demand, and Python answers each gene the browser asks for.

Inside a Jupyter notebook, show an AnnData you already have in memory:

```python
from cellucid import show_anndata

viewer = show_anndata(adata, dataset_name="My dataset", dataset_id="my-dataset")
```

The viewer opens in the output cell and can be driven from Python
(`viewer.set_color_by("leiden")`). See
{doc}`../../python_package/a_landing_pages/04_quick_start_3_levels`.

## Route 3: export it from R

`cellucid_prepare()` writes your object's embeddings, metadata and (optionally)
gene expression into an export folder. In the viewer, press **Prepared** under
**Open data** and choose that folder. The recipes extract everything from a
Seurat or SingleCellExperiment object:

- {doc}`../../r_package/e_integrations_recipes/01_seurat_recipe`
- {doc}`../../r_package/e_integrations_recipes/02_singlecellexperiment_recipe`

## What to do once it is open

- **Colour by anything**: a cluster column, a QC metric, or a gene, from
  **Coloring & filtering**.
- **Filter** to a subset of categories or a value range, and the plot follows.
- **Compare** two colourings, two subsets or the 2D and 3D layout side by side
  with **Keep view** under **Compare views**.
- **Select cells** with lasso, neighbourhood or proximity tools and compare the
  selection in **Analysis** (differential expression, marker genes, gene-set
  scores).
- **Export a figure** as vector SVG or high-resolution PNG.

{doc}`03_quick_tour_60_seconds` walks through each of these.

## Share it

A dataset exported to a public GitHub repository opens from a link, a README
badge or an embedded frame on a web page; see
{doc}`../l_sessions_sharing/13_share_links_badges_and_embeds`. To send your exact
view — camera, colours, filters and selections — save a session file.

See also {doc}`07_cellucid_and_other_viewers` for how Cellucid relates to other
interactive single-cell viewers.
