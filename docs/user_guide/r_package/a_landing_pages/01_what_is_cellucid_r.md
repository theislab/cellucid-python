# What is cellucid-r?

**Audience:** wet lab scientists, computational biologists, engineers  
**Time:** 5 minutes  
**Goal:** understand what the R package does, and what workflow it supports

`cellucid-r` (package name: `cellucid`) does two things from R:

- it **opens** a Seurat, SingleCellExperiment or SummarizedExperiment object
  (or a file) in the **Cellucid viewer**, in the RStudio or Positron Viewer
  pane or a browser, served from memory with no file written and no Python;
- it **writes** the object, or matrices and data frames, as a
  **Cellucid-ready AnnData** store that the **Cellucid web app** opens
  directly (several assays: a **Cellucid-ready MuData**).

If you remember one sentence:

> `cellucid_view(seu)` shows your object in Cellucid; `cellucid_prepare()`
> writes it as a standard *AnnData Zarr store* that anyone opens in a browser —
> or in Python with `anndata`.

```{note}
Cellucid itself is a **web app** (the viewer UI). `cellucid-python` and
`cellucid-r` are peer packages that prepare data and serve it to the viewer
from Python and R. `cellucid-annotation` is for community annotation
workflows.
```

## What `cellucid-r` does

### Fast path (wet lab / non-technical)
- Look at an object: `viewer <- cellucid_view(seu)`.
- Share it: one R call, `cellucid::cellucid_prepare("pbmc.zarr", "zarr", data = seu, ...)`,
  writes one **`.zarr` folder** (or, with `format = "zarr.zip"`, one **`.zarr.zip` file**).
- Anyone opens that folder or file in the Cellucid web app and explores it:
  embeddings, metadata colouring, gene expression, and more.

### Practical path (computational users)
`cellucid_prepare()` writes a standard AnnData store (AnnData Zarr, format 2,
or 3 with `zarr_format = 3`). With `data =` every component comes from the
object; any argument you pass replaces that component:

| What you pass | Where it lands |
|---|---|
| `data` (Seurat, SingleCellExperiment, SummarizedExperiment) | every element below; several assays as modalities (`mod/<name>`) |
| `embeddings` (named list of matrices) | `obsm[<key>]`, float32, every column |
| `obs` (data frame) | `obs`, every column, with anndata's column encodings |
| `var` (data frame) | `var`, every column |
| `gene_expression` (cells × genes) | `X`, gene-major float32 (float64 for integer counts above 2^24; `csc_matrix` from a sparse matrix, stored zeros dropped) |
| `connectivities` (cells × cells) | `obsp['connectivities']`, float32 `csr_matrix` |
| `vector_fields` | `obsm['<field>_<basis>']` for each embedding `X_<basis>` (`velocity_umap`) |
| dataset name, id, citation, source, launch defaults | `uns['cellucid']` and `uns['title']` |

Its dependencies are light:
- Hard dependencies: `jsonlite`, and for the viewer `httpuv`, `later` and `curl`
- Optional: `Matrix` (sparse expression matrices and neighbour graphs;
  recommended), `SeuratObject`, `SingleCellExperiment` /
  `SummarizedExperiment` (each only to read its own objects)
- A few small C routines (checksums, the outlier quantiles, one-pass matrix
  checks, the decoders of a saved session), so a source install needs the
  platform's R build toolchain.

### Deep path (expert / maintainer)
- Expression is stored **gene-major** and in chunks sized for range reads, so
  the viewer fetches one gene with one small read.
- Every chunk is byte-shuffled where that helps and **deflate**-compressed,
  which browsers decompress natively.
- Consolidated metadata (`.zmetadata`) lets a store on a web host open with one request.
- Embeddings are stored at their **own coordinates**; the viewer scales them to
  its view.
- With `latent_key`, each categorical `obs` column's **outlier quantiles** in
  that latent space are computed in C and stored in
  `uns['cellucid']['outlierQuantiles']`, the same bytes Python's
  `cellucid.prepare()` stores, so the viewer's outlier filter works as soon as
  the file opens.
- `cellucid_view()` serves exactly the store `cellucid_prepare(format = "zarr")`
  would write, held in memory and encoded chunk by chunk as the viewer asks.

## What `cellucid-r` does not do

- It writes **no HDF5**: no `.h5ad` and no `.h5mu` (Python's
  `cellucid.prepare()` writes both), so it depends on no HDF5 library.
  `cellucid_view()` still opens an `.h5ad` or `.h5mu` by path; the browser
  reads the file itself.
- It has no Shiny or htmlwidgets integration; the viewer is a page served by a
  small local web server (`httpuv`), shown in the IDE's Viewer pane or a
  browser.
- It writes one expression matrix per modality: `counts` beside `data`, other
  layers and `scale.data` are not written.

## The core mental model

1) **Row order defines cell identity.**  
   Row `i` of every embedding, of `obs`, of `gene_expression`, of the graph,
   and of every vector field is the same cell. Where several inputs carry row
   names, they must agree.
2) **The store is the dataset.**  
   The browser reads a Cellucid-ready AnnData, from a file, a web host, or the
   store `cellucid_view()` holds in memory for an R object.
3) **The web app is the UI.**  
   R prepares, serves and talks to the viewer; the browser is where you
   interact.

## Common misconceptions

- “Can I write from Seurat/SCE directly?”  
  Yes: `cellucid_prepare(out_path, format, data = seu)`. See {doc}`../e_integrations_recipes/index`.
- “Do I need gene expression to use Cellucid?”  
  No. Embeddings alone are enough; expression enables gene colouring and search.
- “Is this a Cellucid-only format?”  
  No. It is plain AnnData: `anndata.read_zarr("my_dataset.zarr")` reads it in
  Python (`mudata.read_zarr()` for several assays).

## Next steps

- Copy/paste minimal example: {doc}`04_quick_start_3_levels`
- Seurat or SingleCellExperiment recipes: {doc}`../e_integrations_recipes/index`
- The viewer from R: {doc}`../d_viewing_loading/04_view_in_rstudio_positron`
- Detailed API and store layout: {doc}`../c_data_preparation_api/index`
