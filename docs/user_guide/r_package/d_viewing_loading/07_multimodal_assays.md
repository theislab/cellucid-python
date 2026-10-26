# Multimodal Assays (RNA + ADT, RNA + ATAC)

**Audience:** CITE-seq and multiome users
**What you'll do:** view and write an object with several feature spaces.

A Seurat object with several assays, or a SingleCellExperiment with alternative
experiments, is viewed and written as a **Cellucid-ready MuData**: one
modality per assay, the cells' metadata, embeddings and neighbour graph shared.
The viewer adds a **Modality** choice to Data mapping; gene search, colour by
gene, markers and dot plots then read the chosen modality.

```r
viewer <- cellucid_view(cite)                 # RNA and ADT assays
cellucid_prepare("cite.zarr", "zarr", data = cite,
                 dataset_id = "cite", dataset_name = "CITE-seq")
```

## Which modalities, in which order

- By default every assay is written, the **active assay first** (Seurat), or
  the **main experiment first** (SingleCellExperiment, named by
  `mainExpName()`, else `main`).
- `modalities = c("RNA", "ADT")` chooses and orders them; one name writes a
  plain AnnData of that assay.
- The viewer opens on `modality = "ADT"` when you name one, else on the first
  modality named like gene expression (`rna`, `gex`, `gene expression`,
  `gene_expression`, `sct`, any case), else the first.

## What is written

| Element | Holds |
| --- | --- |
| `obs`, `obsm`, `obsp/connectivities`, `uns['cellucid']` | the cells: shared by every modality |
| `obsm/<modality>` | which cells each modality measured (mudata's masks) |
| `mod/<name>` | each modality: gene-major `X` and its `var` |
| `obsmap/<name>` | each cell's row in the modality, `0` where it was not measured |
| `mod` attribute `mod-order` | the modality order |

It is MuData's own Zarr layout: `mudata.read_zarr("cite.zarr")` reads it in
Python, and `cellucid.prepare(data = mdata)` writes the same bytes from the same
data. R writes MuData as Zarr only (no `.h5mu`; this package depends on no HDF5
library). The contract is
{ref}`the Cellucid-ready MuData section <r-cellucid-ready-mudata>` of the
specification.

The gene arguments (`gene_expression`, `var`, `gene_identifiers`)
describe one feature space, so with several they are refused: change the assay
in the object first, or write one modality.
