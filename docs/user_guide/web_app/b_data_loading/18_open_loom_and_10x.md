# Open Loom and 10x Cell Ranger files

**Audience:** anyone with a `.loom` file or a Cell Ranger output  
**Time:** 5 minutes  
**What you’ll learn:** how to open a Loom file, a 10x Cell Ranger `.h5` and a
Cell Ranger `outs/` folder directly, what Cellucid reads from each, and what to
do with a count matrix that has no embedding  
**Prerequisites:** none

Cellucid opens these files as they are: there is no conversion step. Each is
read in place, like an `.h5ad`: opening reads the file's structure and the
values the first view needs, and a gene or an embedding is read when you show
it ({doc}`07_file_format_expectations`).

---

## Loom (`.loom`)

**Open… › File…**, drop the file on the viewer, or paste a link to it
({doc}`14_open_from_link`). Loom 2.0.1 and 3.0.0 files open, whichever tool
wrote them: loompy (`loompy.create`, `adata.write_loom`), velocyto, SCENIC and
SCope, LoomExperiment, loomR and SeuratDisk.

| In the Loom file | In Cellucid |
|---|---|
| `/matrix` (genes × cells) | the expression matrix `X` |
| `/layers/<name>` (`spliced`, `unspliced`, …) | the layers, listed under Data mapping › Genes › Expression |
| a cell attribute with one value per cell (`ClusterName`, `nUMI`) | an obs field |
| a cell attribute with one row per cell (`X_umap`, n × 2) | an embedding |
| coordinate columns: `_X` / `_Y`, `_tSNE1` / `_tSNE2`, `umap_1` / `umap_2` | one embedding each: `X_embedding`, `X_tSNE`, `X_umap` |
| SCope's compound `Embedding`, `Embeddings_X` / `Embeddings_Y` | `X_embedding`, `X_embedding_<id>` |
| `CellID` | the cell names |
| `Gene` (repeats allowed) | the gene names, made unique as for any file |
| `/col_graphs/<name>` | a neighbour graph |
| `/attrs` (Loom 3) or the file's attributes (Loom 2) | `uns`; `title` names the dataset |

Loom 2.0.1 stores text as ASCII with character references (`&#233;` for
`é`); Cellucid shows the text it stands for. A file written with the matrix
as cells × genes (loomR's order) is recognised by the length of its
attributes and read without transposing.

```{figure} ../../../_static/screenshots/data_loading/loom-dataset-box.png
:alt: The Dataset box after a Loom file was opened: the title stored in the file, its cells and genes, "Loom (.loom) · this computer", and Details open with the file, its size, how it is read, and what it holds: four obs fields and a UMAP.
:width: 488px

A Loom file opened from this computer, with **Details** open. It is read in
place, as an `.h5ad` is: no *Read once* line and nothing kept. Its `title`
attribute names the dataset, and its `X_umap` cell attribute is the embedding.
```

**A Loom file with no embedding is refused**: Cellucid draws cells at
embedding coordinates, and a velocyto output, for example, has only counts.
The message gives the lines that compute one; with scanpy:

```python
import scanpy as sc

adata = sc.read_loom("cells.loom")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata.write_h5ad("cells.h5ad")
```

---

## 10x Cell Ranger: `filtered_feature_bc_matrix.h5` and `outs/`

A Cell Ranger count matrix holds counts and nothing to draw them with. Cell
Ranger computes a UMAP, a t-SNE, a PCA and clusterings beside it, in the
`analysis/` folder of `outs/`, and **opening the `outs/` folder** shows them:

1. In the Session panel, under **Open data**, click **Open…**, then
   **Folder…** (or drag the folder onto the viewer).
2. Choose the `outs/` folder itself (the one holding
   `filtered_feature_bc_matrix.h5` and `analysis/`).
3. The dataset opens on Cell Ranger's UMAP (`X_umap`, written by Cell Ranger
   3.1 and later), else its t-SNE (`X_tsne`); the PCA (`X_pca`) is listed
   under Data mapping, and each clustering (`graphclust`,
   `kmeans_<k>_clusters`) is a categorical field.

```{figure} ../../../_static/screenshots/data_loading/tenx-outs-window.png
:alt: The whole Cellucid window after a Cell Ranger outs/ folder was opened with Open… › Folder…: the Dataset box reads "10x Cell Ranger outs folder · this computer", and the canvas shows Cell Ranger's UMAP coloured by its clustering, gene_expression_clusters.
:width: 1440px

A Cell Ranger `outs/` folder opened with **Open… › Folder…**: the Dataset box
reads *Cell Ranger outs/ folder · this computer*, the view is Cell Ranger's
UMAP, and it is coloured by Cell Ranger's own clustering,
`gene_expression_clusters`, chosen under **Categorical obs**.
```

| In `outs/` | In Cellucid |
|---|---|
| `filtered_feature_bc_matrix.h5` (Cell Ranger 3 and later) or `filtered_gene_bc_matrices_h5.h5` (1 and 2) | the expression matrix, read in place |
| `analysis/umap/<…>_components/projection.csv` | `X_umap` |
| `analysis/tsne/<…>_components/projection.csv` | `X_tsne` |
| `analysis/pca/<…>_components/projection.csv` | `X_pca` |
| `analysis/clustering/<name>/clusters.csv` | the field `<name>` |
| barcodes | the cell names |
| feature names, ids, types and genome | the gene names; `var` columns `gene_ids`, `feature_types`, `genome` (as `scanpy.read_10x_h5` names them) |

Projections and clusterings are matched to the matrix by barcode. A cell a
projection does not list cannot be drawn in it (the embedding is refused with
the cell named); a cell a clustering does not list has no category.

A matrix of several **feature types** (Gene Expression with Antibody Capture
from CITE-seq, CRISPR Guide Capture, Multiplexing Capture) opens on Gene
Expression and shows a **Modality** select under Data mapping › Genes, one
modality per feature type ({doc}`19_multimodal_mudata`). A Cell Ranger 2 file
of several genomes has one modality per genome.

**The `.h5` alone** (**Open… › File…**, a drop or a link) is refused with this
message, because it has no embedding:

```text
filtered_feature_bc_matrix.h5 is a 10x Cell Ranger count matrix that has no
cell embedding, and Cellucid draws cells at embedding coordinates. Cell Ranger
writes UMAP and t-SNE projections beside it: open the whole outs/ folder
(Open… › Folder…), which holds this file and analysis/. Or compute one and open
the result: …
```

Opened from a link (Google Drive, OneDrive, a web address), the advice is what
a link allows: a link opens one file, never a folder, and no provider shares a
folder Cellucid could read, so the message says to download the whole `outs/`
folder and open it from this computer (**Open… › Folder…**), or to compute an
embedding.

To compute an embedding of your own instead, in Python:

```python
import scanpy as sc

adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata.write_h5ad("pbmc.h5ad")
```

or in R, with Seurat, then `cellucid::cellucid_prepare("pbmc.zarr", data = seu)`
after `RunUMAP()`.

`molecule_info.h5` is a table of molecules, not a count matrix: it is refused,
naming the matrix to open. The 10x download site does not let other websites
read its files, so a link to a file there is refused; download it, then open
the file or folder.

---

## Where to go next

- {doc}`19_multimodal_mudata`: the Modality select.
- {doc}`13_data_mapping`: which embedding, graph and genes a view shows.
- {doc}`08_troubleshooting_data_loading`: every refusal and what to do.
