# Serve Loom, 10x and MatrixMarket files

**Audience:** anyone with a Cell Ranger output, a MatrixMarket folder or a
Loom file, on this computer or on a cluster
**Time:** 10 minutes
**Goal:** view each one with `cellucid serve` (or `show_anndata`), with no
conversion to AnnData first

`cellucid serve` and `show_anndata()` hand a file or folder to the viewer. A
Loom file and a Cell Ranger `outs/` folder are served as a store with their
cells in priority order, made from the file on the kernel (h5py) as the
viewer reads it, holding exactly what the viewer shows from the file itself:
a subsample the viewer opens transfers only its cells. A MatrixMarket matrix
is served as it is and read once from the start. Serving needs no scanpy and
no loompy. Each needs a cell embedding, because Cellucid draws cells at
embedding coordinates.

The data is 10x Genomics' public 1k PBMC dataset (Cell Ranger 3.0).

## A Cell Ranger `outs/` folder

Cell Ranger writes its t-SNE (and from 3.1 a UMAP), PCA and clusterings in
`analysis/`, beside the count matrix. A folder Cell Ranger wrote is already
laid out this way; this one is put together from the dataset's downloads:

```sh
curl -LO https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
curl -LO https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_analysis.tar.gz
mkdir outs
cp pbmc_1k_v3_filtered_feature_bc_matrix.h5 outs/filtered_feature_bc_matrix.h5
tar -xzf pbmc_1k_v3_analysis.tar.gz -C outs

cellucid serve outs --dataset-name "PBMC 1k (Cell Ranger outs)" --dataset-id pbmc-1k-outs
```

```text
[1/4] Detecting format...
      Format: Cell Ranger outs/ folder, a folder of 27 files
      ✓ Format detected
...
```

The browser opens on the 1,222 cells' t-SNE (`X_tsne`), with `X_pca` and the
clusterings `graphclust` and `kmeans_2_clusters` to `kmeans_10_clusters` as
fields. The terminal says `Serving: a Zarr store in priority order, made as
the browser reads it`; the matrix is read from the file one gene at a time,
for the cells open, as you show it. A matrix of several feature types (Gene
Expression with Antibody Capture) is served as a MuData store, one modality
per type, as the viewer shows the folder itself.

A dataset that is not a Cellucid-ready AnnData does not name itself, so it
is named after its file (its id the file's name without the extension,
printed when the server starts); `--dataset-name` and `--dataset-id` name it
otherwise. The id is what sessions are kept under.

## A MatrixMarket folder

The same counts as `matrix.mtx.gz`, `features.tsv.gz` and `barcodes.tsv.gz`,
with the same `analysis/` beside them for the embedding:

```sh
curl -LO https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz
mkdir pbmc_mtx
tar -xzf pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz -C pbmc_mtx
tar -xzf pbmc_1k_v3_analysis.tar.gz -C pbmc_mtx

cellucid serve pbmc_mtx --dataset-name "PBMC 1k (MatrixMarket)" --dataset-id pbmc-1k-mtx
```

The format line says `MatrixMarket folder, a folder of 29 files`. A
compressed text matrix has no random access, so the browser reads it once
when it opens, into its own storage, and every view reads from there. Any
`.csv` or `.tsv` table of the cells in the folder (a `metadata.tsv` with
`UMAP_1`/`UMAP_2` columns) gives the embedding the same way `analysis/` does
({doc}`../../web_app/b_data_loading/17_open_r_objects`).

## A Loom file

A Loom file opens when it carries an embedding, as a two-column cell
attribute (`X_umap`) or as coordinate columns (`_X`/`_Y`). This one is
written with scanpy and loompy (`pip install scanpy loompy`) from the same
`.h5`, to have one:

```python
import scanpy as sc

adata = sc.read_10x_h5("pbmc_1k_v3_filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
adata.write_loom("pbmc_1k.loom", write_obsm_varm=True)
```

```sh
cellucid serve pbmc_1k.loom --dataset-name "PBMC 1k (Loom)" --dataset-id pbmc-1k-loom
```

The file opens on its UMAP, with `X_pca` and the field `leiden`; it is served
in priority order like the `outs/` folder, the matrix read from the file
(transposed: Loom stores genes by cells) for the cells the viewer opens.

## From a notebook

`show_anndata()` takes the same paths and embeds the viewer:

```python
import cellucid

viewer = cellucid.show_anndata("outs", dataset_name="PBMC 1k (Cell Ranger outs)", dataset_id="pbmc-1k-outs")
```

For any of these, the hooks, selections and session bundles of
{doc}`23_programmatic_highlighting_and_selection_callbacks` work as for an
AnnData.

## When something is refused

- `filtered_feature_bc_matrix.h5` served alone has no embedding, and is
  refused in step 1 saying so: serve the `outs/` folder that holds it and
  `analysis/`.
- A folder with neither a Cell Ranger matrix, a MatrixMarket matrix nor Zarr
  metadata: `Unable to detect format for: <path>`.
- A Loom file with no embedding is refused in step 1 with the lines that
  compute one; a MatrixMarket dataset with none, by the viewer ({doc}`../../web_app/b_data_loading/08_troubleshooting_data_loading`).

## Next

- {doc}`../d_viewing_apis/02_loading_routes`: every format by path
- {doc}`../../web_app/b_data_loading/20_tutorials_open_published_data`: the
  same files opened in the browser with no Python
