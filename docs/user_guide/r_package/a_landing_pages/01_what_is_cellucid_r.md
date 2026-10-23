# What is cellucid-r?

**Audience:** wet lab scientists, computational biologists, engineers  
**Time:** 5 minutes  
**Goal:** understand what the R package does, and what workflow it supports

`cellucid-r` (package name: `cellucid`) is a **minimal R package** with one job:
write your single-cell data as a **Cellucid-ready AnnData** store that the
**Cellucid web app** opens directly.

If you remember one sentence:

> `cellucid_prepare()` turns *R matrices and data frames* into a standard
> *AnnData Zarr store* that you open in your browser (Cellucid) — or in Python
> with `anndata`.

```{note}
Cellucid itself is a **web app** (the viewer UI). `cellucid-python` and
`cellucid-r` are **helper packages** that prepare data and (in Python) provide
notebook/server integrations. `cellucid-annotation` is for community annotation
workflows.
```

## What `cellucid-r` does

### Fast path (wet lab / non-technical)
- You (or a collaborator) run one R function: `cellucid::cellucid_prepare(...)`.
- It writes one **`.zarr` folder** (or, with `format = "zarr.zip"`, one **`.zarr.zip` file**).
- You open the Cellucid web app and open that folder or file.
- You explore interactively: embeddings, metadata colouring, gene expression, and more.

### Practical path (computational users)
`cellucid_prepare()` writes a standard AnnData store (AnnData Zarr, format 2,
or 3 with `zarr_format = 3`):

| What you pass | Where it lands |
|---|---|
| `embeddings` (named list of matrices) | `obsm[<key>]`, float32, every column |
| `obs` (data frame) | `obs`, every column, with anndata's column encodings |
| `var` (data frame) | `var`, every column |
| `gene_expression` (cells × genes) | `X`, gene-major float32 (float64 for integer counts above 2^24; `csc_matrix` from a sparse matrix, stored zeros dropped) |
| `connectivities` (cells × cells) | `obsp['connectivities']`, float32 `csr_matrix` |
| `vector_fields` | `obsm['<field>_<basis>']` for each embedding `X_<basis>` (`velocity_umap`) |
| dataset name, id, citation, source, launch defaults | `uns['cellucid']` and `uns['title']` |

It is intentionally dependency-light:
- Hard dependency: `jsonlite`
- Optional (recommended): `Matrix`, for sparse expression matrices and neighbour graphs
- One small C routine (the CRC-32 checksum of the `.zarr.zip` writer), so a
  source install needs the platform's R build toolchain.

### Deep path (expert / maintainer)
- Expression is stored **gene-major** and in chunks sized for range reads, so
  the viewer fetches one gene with one small read.
- Every chunk is **zlib**-compressed, which browsers decompress natively.
- `.zmetadata` (consolidated metadata) lets a store on a web host open with one request.
- Embeddings are stored at their **own coordinates**; the viewer scales them to
  its view.
- With `latent_key`, each categorical `obs` column's **outlier quantiles** in
  that latent space are computed in C and stored in
  `uns['cellucid']['outlierQuantiles']`, the same bytes Python's
  `cellucid.prepare()` stores, so the viewer's outlier filter works as soon as
  the file opens.

## What `cellucid-r` does NOT do (yet)

- It does **not** embed Cellucid inside RStudio / Shiny (no widget yet).
- It does **not** run a server; the viewer reads the store's files directly.
- It does **not** read sessions back into your data. A `.cellucid-session`
  is checked against, and applied to, an AnnData — its dataset id in
  `uns['cellucid']`, its `obs`, `var`, `obsm` and `obsp` names — and this
  package writes AnnData stores but never reads one, so a Seurat or
  SingleCellExperiment object has nothing to check it against. Open the store
  in Python and use `cellucid.apply_cellucid_session_to_anndata()`
  ({doc}`../../python_package/e_jupyter_hooks/11_session_bundles_apply_session_to_anndata`).
- It does **not** depend on Seurat or SingleCellExperiment; you **extract
  matrices/data frames** yourself (see the recipes).

If you want a “viewer object” you can control from a notebook, see {doc}`../../python_package/index`.

## The core mental model

1) **Row order defines cell identity.**  
   Row `i` of every embedding, of `obs`, of `gene_expression`, of the graph,
   and of every vector field is the same cell. Where several inputs carry row
   names, they must agree.
2) **The store is the dataset.**  
   The browser does not load an R object — it reads the store's files.
3) **The web app is the UI.**  
   R is for preparation. The browser is for interaction.

## Common misconceptions

- “Can I write from Seurat/SCE directly?”  
  Not as a single call (by design), but extraction takes a few lines. Start with {doc}`../e_integrations_recipes/index`.
- “Do I need gene expression to use Cellucid?”  
  No. Embeddings alone are enough; expression enables gene colouring and search.
- “Is this a Cellucid-only format?”  
  No. It is plain AnnData: `anndata.read_zarr("my_dataset.zarr")` reads it in Python.

## Next steps

- Copy/paste minimal example: {doc}`04_quick_start_3_levels`
- Seurat or SingleCellExperiment recipes: {doc}`../e_integrations_recipes/index`
- Detailed API and store layout: {doc}`../c_data_preparation_api/index`
