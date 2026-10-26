# CITE-seq as a MuData: prepare and show

**Audience:** anyone with RNA and surface-protein (ADT) counts of the same cells
**Time:** 10 minutes
**Goal:** build a MuData from a 10x CITE-seq matrix, write it as a
Cellucid-ready `.h5mu`, and view it in a notebook with its two modalities

You need `pip install "cellucid[mudata]" scanpy` (scanpy for the embedding;
`sc.tl.leiden(flavor="igraph")` also needs `igraph`). The data is 10x Genomics'
public 1k PBMC CITE-seq run (Cell Ranger 3.0): RNA and 17 antibodies.

```sh
curl -LO https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5
```

## Step 1: one matrix, two feature types

```python
import scanpy as sc
import mudata as md
import cellucid

adata = sc.read_10x_h5("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5", gex_only=False)
adata.var_names_make_unique()
print(adata)
print(adata.var["feature_types"].value_counts())
```

```text
AnnData object with n_obs × n_vars = 713 × 33555
    var: 'gene_ids', 'feature_types', 'genome', 'pattern', 'read', 'sequence'
feature_types
Gene Expression     33538
Antibody Capture       17
```

## Step 2: a MuData of two modalities

Split the matrix by feature type, compute the embedding on RNA, and keep the
shared parts (the cells' clusters, the UMAP and the neighbour graph) at the
MuData's root, where every modality shares them:

```python
rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
prot = adata[:, adata.var["feature_types"] == "Antibody Capture"].copy()

sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.pca(rna)
sc.pp.neighbors(rna)
sc.tl.umap(rna)
sc.tl.leiden(rna, flavor="igraph", n_iterations=2)
sc.pp.log1p(prot)

mdata = md.MuData({"rna": rna, "prot": prot})
mdata.obs["leiden"] = rna.obs["leiden"]
mdata.obsm["X_umap"] = rna.obsm["X_umap"]
mdata.obsp["connectivities"] = rna.obsp["connectivities"]
print(mdata)
```

```text
MuData object with n_obs × n_vars = 713 × 33555
  obs:	'leiden'
  obsm:	'X_umap'
  obsp:	'connectivities'
  2 modalities
    rna:	713 × 33538
      ...
    prot:	713 × 17
      ...
```

## Step 3: write a Cellucid-ready `.h5mu`

```python
cellucid.prepare(
    "pbmc_1k_citeseq.h5mu",
    format="h5mu",
    data=mdata,
    dataset_id="pbmc-1k-citeseq",
    dataset_name="PBMC 1k CITE-seq (10x)",
    modality="rna",
    embedding_defaults={"2d": "X_umap"},
    obs_keys=["leiden"],
)
```

`data=mdata` takes every component the MuData holds; each modality is written
as a Cellucid-ready AnnData (gene-major, chunked for the browser), and
`modality="rna"` is the one the viewer opens on. `format="zarr"` writes the
same as a MuData Zarr folder, the form to publish on GitHub
({doc}`../c_data_preparation_api/12_prepare_from_anndata_or_mudata`). The file
is an ordinary MuData:

```python
back = md.read_h5mu("pbmc_1k_citeseq.h5mu")
print(back.uns["cellucid"]["defaults"])
```

```text
{'embedding2d': 'X_umap', 'modality': 'rna', 'obsKeys': array(['leiden'], dtype=object)}
```

## Step 4: show it

```python
viewer = cellucid.show_anndata("pbmc_1k_citeseq.h5mu")
```

The 713 cells open on the UMAP, with the field `leiden`. In the viewer's
**Session › Data mapping › Genes**, **Modality** lists `rna` (33,538
features) and `prot` (17). Choose **prot**, click **Apply**, and search a
protein (`CD3_TotalSeqB`, `CD19_TotalSeqB`) to colour the cells by it; marker
genes and the dot plot then compare proteins
({doc}`../../web_app/b_data_loading/19_multimodal_mudata`).

The in-memory MuData shows the same way, written once to a temporary `.h5mu`;
it names itself only once prepared, so give it a name and an id. `modality=`
opens it on the proteins:

```python
viewer = cellucid.show_anndata(
    mdata,
    dataset_name="PBMC 1k CITE-seq (in memory)",
    dataset_id="pbmc-1k-citeseq-live",
    modality="prot",
)
```

`cellucid serve pbmc_1k_citeseq.h5mu` serves the file from a terminal, and
the web app opens it from **Open… › File…** with no Python at all.

## Next

- {doc}`25_serve_loom_10x_and_matrixmarket`: Loom, 10x and MatrixMarket files from Python
- {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`
  (§5: the Cellucid-ready MuData)
- R: a Seurat object with RNA and ADT assays is written as the same MuData
  ({doc}`../../r_package/f_notebooks_tutorials/06_citeseq_seurat_to_mudata`)
