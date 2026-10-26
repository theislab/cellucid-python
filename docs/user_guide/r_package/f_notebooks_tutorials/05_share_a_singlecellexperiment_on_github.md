# Share a SingleCellExperiment on GitHub

**Audience:** anyone who wants others to open their dataset from one link
**Time:** 15 minutes
**Goal:** write a SingleCellExperiment as a Cellucid-ready AnnData Zarr
folder, check it, publish it in a public GitHub repository, and share the
link

You need the `cellucid` and `SingleCellExperiment` packages, `git`, and a
GitHub account. A Zarr folder is many small files, so GitHub's 100 MB limit
per file never applies and no Git LFS is needed
({doc}`../../web_app/b_data_loading/15_host_on_github`).

## Step 1: the object

Use your own SingleCellExperiment. Here one is built from SeuratObject's
80-cell example, with its counts, log-normalised values, cell metadata, PCA
and t-SNE:

```r
library(SingleCellExperiment)
library(SeuratObject)
data("pbmc_small")

sce <- SingleCellExperiment(
  assays = list(
    counts = LayerData(pbmc_small, layer = "counts"),
    logcounts = LayerData(pbmc_small, layer = "data")
  ),
  colData = pbmc_small[[]],
  reducedDims = list(
    PCA = Embeddings(pbmc_small, "pca"),
    TSNE = Embeddings(pbmc_small, "tsne")
  )
)
sce
```

```text
class: SingleCellExperiment 
dim: 230 80 
metadata(0):
assays(2): counts logcounts
rownames(230): MS4A1 CD79B ... SPON2 S100B
rowData names(0):
colnames(80): ATGCCAGAACGACT CATGGCCTGTGCAT ... GGAACACTTCAGAC
  CTTGATTGATCTTC
colData names(7): orig.ident nCount_RNA ... groups RNA_snn_res.1
reducedDimNames(2): PCA TSNE
mainExpName: NULL
altExpNames(0):
```

## Step 2: name it and write it

The dataset's id is what links, saved sessions and community annotations use:
choose it once and keep it. The object can carry it, with its name and
description, in `metadata(sce)$cellucid`, so every later write and every
viewer uses the same:

```r
library(cellucid)
metadata(sce)$cellucid <- list(
  id = "pbmc-small-sce",
  name = "PBMC small (SingleCellExperiment)",
  description = "SeuratObject's 80-cell example, as a SingleCellExperiment."
)
cellucid_prepare("pbmc_small.zarr", format = "zarr", data = sce)
```

`data = sce` takes everything the object holds: `logcounts` as the expression
matrix (gene-major, so the browser reads one gene at a time), `colData` as the
cell metadata, each reduced dimension as an embedding (`X_PCA`, `X_TSNE`),
`colPairs` as the neighbour graph. Any argument you pass replaces one of them
({doc}`../c_data_preparation_api/01_cellucid_prepare_overview`).

## Step 3: check it before publishing

```r
files <- list.files("pbmc_small.zarr", recursive = TRUE, all.files = TRUE)
length(files)
max(file.size(file.path("pbmc_small.zarr", files)))
```

```text
[1] 105
[1] 8937
```

105 files, the largest 9 KB. Look at it in the viewer exactly as others will:

```r
viewer <- cellucid_view("pbmc_small.zarr")
```

The 80 cells open on the t-SNE, named *PBMC small (SingleCellExperiment)*.
`cellucid_stop(viewer)` when you are done. In Python,
`anndata.read_zarr("pbmc_small.zarr")` reads the same folder.

## Step 4: publish it

Create a **public** repository on GitHub (here `OWNER/REPOSITORY`), then, in a
terminal, in the folder that holds `pbmc_small.zarr`:

```sh
git init
git add pbmc_small.zarr
git commit -m "Add the PBMC small dataset"
git branch -M main
git remote add origin https://github.com/OWNER/REPOSITORY.git
git push -u origin main
```

Use Git from the command line: GitHub's upload page takes at most 100 files at
a time.

## Step 5: share the link

Open [Cellucid](https://www.cellucid.com), and in **Session › Link or server**
enter `OWNER/REPOSITORY/pbmc_small.zarr` and click **Open**. The address bar
now carries the link to share:

```text
https://www.cellucid.com/?data=gh:OWNER/REPOSITORY/pbmc_small.zarr
```

Anyone who opens it reads only what each view needs, straight from GitHub.
Several datasets in one repository can be listed in a catalog
(`datasets.json`) and picked from the **Sample dataset** list
({doc}`../../web_app/b_data_loading/11_custom_dataset_repository`); the
repository [theislab/cellucid-demo-custom-datasets](https://github.com/theislab/cellucid-demo-custom-datasets)
is a complete example.

## Next

- {doc}`../d_viewing_loading/02_host_stores_for_sharing`: Hugging Face,
  Zenodo, or any web host
- {doc}`06_citeseq_seurat_to_mudata`: an object with several assays
