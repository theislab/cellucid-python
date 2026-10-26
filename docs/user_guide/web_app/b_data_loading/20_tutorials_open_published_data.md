# Tutorials: open published data as it is

**Audience:** anyone who found a dataset in a paper's data availability
section and wants to look at it  
**Time:** 5 minutes each  
**What you'll do:** open three kinds of published single-cell data with no
conversion and nothing installed: a GEO sample's count matrix by its link, a
Seurat object on Zenodo by its link, and a 10x Cell Ranger output folder from
your computer  
**Prerequisites:** a current desktop Chrome, Edge, Firefox or Safari

Each tutorial uses a real public dataset. The numbers it states (cells, what
was read, how long it took) were measured on 2026-10-24 in Chromium against
the current web app; your connection decides the times.

---

## 1. A GEO sample's MatrixMarket files, by link

GEO keeps many count matrices as three files per sample: the matrix
(`…matrix.mtx.gz`), its features and its barcodes. Some submitters add a
table of each cell's UMAP coordinates beside them, and then Cellucid can
draw the sample. Sample GSM9843410 of series
[GSE336826](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336826)
(zebrafish embryos at 12 hours post fertilization, 10x Multiome) has five
files:

| File | What it holds |
|---|---|
| `GSM9843410_12hpf_scRNA_matrix.mtx.gz` (161 MB) | the counts |
| `GSM9843410_12hpf_scRNA_features.tsv.gz` | the features: 32,058 genes and 239,160 chromatin peaks |
| `GSM9843410_12hpf_scRNA_barcodes.tsv.gz` | the 4,278 cells |
| `GSM9843410_12hpf_umap_projection.csv.gz` | each cell's `UMAP-1`, `UMAP-2` |
| `GSM9843410_12hpf_clusters_ID.csv.gz` | each cell's `Cluster` |

1. Open [Cellucid](https://www.cellucid.com) and expand **Session**.
2. Paste this into the **Link or server** bar under **Open data** and click
   **Open**:

   ```text
   geo:GSE336826/GSM9843410_12hpf_scRNA_matrix.mtx.gz
   ```

   or open [this link](https://www.cellucid.com/?data=geo:GSE336826/GSM9843410_12hpf_scRNA_matrix.mtx.gz),
   which carries the same address.
3. Cellucid lists the series' files, pairs the matrix with the features and
   barcodes that share its name, and finds the two cell tables whose names
   start with the same sample prefix (`GSM9843410_12hpf_`). The notification
   shows the matrix downloading: a compressed text file has no random access,
   so it is **read once** from the start, into this browser's storage
   ({doc}`17_open_r_objects`). Nothing is saved to your Downloads folder.
4. After about 170 MB and a minute and a half, the sample shows its
   **4,278 cells** on the submitter's UMAP (`X_umap`, from `UMAP-1` and
   `UMAP-2`), with the field **Cluster** (each cell's cluster, as categories:
   a cell-table column named for a clustering is categorical even when its
   labels are numbers).

Search a gene (`sox17`) to colour by it. The features of a Multiome sample
like this one are of two types, so they open as two modalities, as a 10x
`.h5` of several feature types does: **Gene Expression** (the 32,058 genes,
opened first) and **Peaks** (the 239,160 peaks), chosen under
**Session › Data mapping › Genes › Modality** ({doc}`19_multimodal_mudata`).

A GEO sample without such a table is refused with the lines that compute an
embedding in Python and in R, because Cellucid draws cells at embedding
coordinates ({ref}`no-cell-embedding`). A series that holds several samples
names each one's link when you paste the series alone (`geo:GSE336826`).

---

## 2. A Seurat object on Zenodo, by link

Many studies share their analysis as a Seurat object saved with `saveRDS()`.
Zenodo record [8071138](https://zenodo.org/records/8071138) (Tauber et al.,
mast cells across organs, CC0) holds one of 371 mouse mast cells,
`Mice_Mast_cells.Rds`, 20 MB.

1. Paste this into the **Link or server** bar and click **Open**:

   ```text
   zenodo:8071138/Mice_Mast_cells.Rds
   ```

   or open [this link](https://www.cellucid.com/?data=zenodo:8071138/Mice_Mast_cells.Rds).
2. The `.rds` is read once as it downloads (19 MB, about 17 seconds): Cellucid
   reads R's own file format itself, so no R is needed.
3. The object opens on its UMAP (`X_umap`, from the reduction `umap`) with
   `X_pca` beside it, and every column of its `meta.data` as a field
   (`seurat_clusters`, `Group`, `percent.mt`, …), plus `active_ident`.
4. The object has three assays, so **Session › Data mapping › Genes** shows a
   **Modality** select: `integrated` (the object's active assay, 2,000
   genes), `RNA` (28,936 genes) and `decontX`. Choose **RNA** and click
   **Apply** to search every gene ({doc}`19_multimodal_mudata`).

Cellucid passes over what it does not use (each assay's `scale.data`, the
PCA's feature loadings, the object's command log) as the file streams. The
same works for a SingleCellExperiment `.rds` and for an `.rda` or
`.RData` file ({doc}`17_open_r_objects`).

:::{note}
A GEO supplementary file named `.rds` (not `.rds.gz`) saved with
`saveRDS()`'s default gzip compression is labelled gzip-encoded by NCBI's
file host even in parts of it, and no browser can decode a part. Cellucid
notices and downloads the whole file once instead (the notification says
why), then reads it as usual.
:::

---

## 3. A Cell Ranger `outs/` folder, from your computer

A Cell Ranger count matrix holds counts and no embedding; Cell Ranger writes
its UMAP or t-SNE, PCA and clusterings beside it, in `analysis/`. Opening the
whole `outs/` folder shows them. 10x Genomics' own site does not let other
websites read its files, so download the files first. This uses the public
1k PBMC dataset (Cell Ranger 3.0):

```sh
curl -LO https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
curl -LO https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_analysis.tar.gz
mkdir outs
cp pbmc_1k_v3_filtered_feature_bc_matrix.h5 outs/filtered_feature_bc_matrix.h5
tar -xzf pbmc_1k_v3_analysis.tar.gz -C outs
```

(A folder Cell Ranger wrote itself is already laid out this way.)

1. In the Session panel, under **Open data**, click **Open…**, then
   **Folder…**, and choose the `outs` folder itself. Or drag the folder onto
   the viewer.
2. The matrix is read in place, like an `.h5ad`: opening reads only its
   structure and barcodes.
3. **1,222 cells** appear on Cell Ranger's t-SNE (`X_tsne`; Cell Ranger 3.1
   and later also write a UMAP, which is then the default), with `X_pca`
   beside it and each clustering as a field: `graphclust` and
   `kmeans_2_clusters` to `kmeans_10_clusters`.

The `.h5` on its own is refused, naming the folder to open instead. A folder
with `filtered_feature_bc_matrix/` (the `matrix.mtx.gz`, `features.tsv.gz`
and `barcodes.tsv.gz` form) in place of the `.h5` opens the same way, read
once rather than in place ({doc}`17_open_r_objects`).

---

## Related pages

- {doc}`14_open_from_link`: every provider and link form
- {doc}`17_open_r_objects`: R objects and MatrixMarket folders, what is read
- {doc}`18_open_loom_and_10x`: Loom files and Cell Ranger output
- {doc}`../../python_package/f_notebooks_tutorials/25_serve_loom_10x_and_matrixmarket`:
  the same files served from Python
