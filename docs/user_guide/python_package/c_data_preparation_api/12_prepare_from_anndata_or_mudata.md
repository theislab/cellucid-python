# Prepare an AnnData or a MuData whole (`data=`)

**Audience:** anyone holding an AnnData or a MuData in Python
**Time:** 5 minutes
**Goal:** write a Cellucid-ready file from the object itself, and know what is taken from it

{func}`~cellucid.prepare` takes the object with `data=`. Every component the
object holds becomes the default of the matching argument, and an argument you
pass replaces that one component. The file is byte for byte what the same
components passed one by one write. The R package's `cellucid_prepare(data =
...)` follows the same rule for Seurat and SingleCellExperiment objects.

## From an AnnData

```python
import anndata as ad
import cellucid

adata = ad.read_h5ad("pbmc.h5ad")
cellucid.prepare("pbmc_cellucid.h5ad", format="h5ad", data=adata,
                 dataset_id="pbmc-3k", dataset_name="PBMC 3k")
```

| The AnnData holds | `prepare()` takes it as |
|---|---|
| `obs` | `obs` |
| `uns['<column>_colors']` of each `obs` column written (Scanpy's colors) | `obs_colors`: the viewer opens the column in them |
| the matrix the viewer opens the AnnData with: `raw.X` and `raw.var` when `raw` holds more genes than `X` or `X` is scaled (Scanpy's `adata.raw = adata`, then the highly variable genes kept and `sc.pp.scale`), else `X` and `var`; without `X`, `raw`, else a layer | `gene_expression` and `var`, a message saying which and why |
| every numeric `obsm` array | `embeddings`, each under its own key |
| `obsm['<field>_<basis>']` as wide as `obsm['X_<basis>']` | `vector_fields` (`velocity_umap` is `velocity` over `X_umap`) |
| `obsp['connectivities']` | `connectivities` |
| `uns['cellucid']` | `dataset_id`, `dataset_name`, `dataset_description`, the citation and source, and the launch defaults (`embedding_defaults`, `latent_key`, `obs_keys`, `var_gene_id_column`, `vector_field_default`, `normalize_embeddings`, `centroid_min_points`) |

Colors that cannot color the column written (not a list of text, fewer
colors than its categories, a column not written as a categorical) are left
out, and a message says why. The other matrices (layers, `raw` or `X`), other
`obsp` matrices, `obsm` tables and the rest of `uns` are not written: a Cellucid-ready AnnData holds `X` alone
({doc}`09_cellucid_ready_anndata_specification`). A file `prepare()` wrote
before carries its own `uns['cellucid']`, so preparing it again needs no
identity:

```python
cellucid.prepare("pbmc_v2.zarr", format="zarr", data=ad.read_h5ad("pbmc_cellucid.h5ad"))
```

To choose a component, pass it: `embeddings={"X_umap": adata.obsm["X_umap"]}`
writes that embedding alone, `gene_expression=adata.layers["log1p"]` writes
that layer as `X` (with `adata.var` as its genes), and
`gene_expression=adata.X, var=adata.var` writes `X` where the object's genes
are read from `raw`.

A Scanpy object is read as Scanpy and the viewer read it: after
`adata.raw = adata`, `sc.pp.highly_variable_genes(..., subset=True)` and
`sc.pp.scale`, its genes are the log-normalized `raw` of every gene, not the
scaled `X` of the highly variable ones, so `prepare(data=adata)` writes
`raw.X` (with `raw.var`) and the file opens on the values the object itself
opens on. To leave out something the
object holds, remove it from a copy first (`adata.obsp.pop("connectivities")`).

## From a MuData

A MuData holds several feature spaces over the same cells: RNA and surface
proteins (CITE-seq), RNA and chromatin. `prepare()` writes it as a
**Cellucid-ready MuData** ({ref}`cellucid-ready-mudata`): `format="h5mu"` for
one `.h5mu` file, `"zarr"` for a MuData Zarr folder, `"zarr.zip"` for that
folder in one file. Reading and building a MuData needs the `mudata` package
(`pip install "cellucid[mudata]"`); writing needs nothing more.

```python
import mudata as md
import cellucid

mdata = md.read_h5mu("citeseq.h5mu")          # modalities "rna" and "prot"
cellucid.prepare("citeseq_cellucid.h5mu", format="h5mu", data=mdata,
                 dataset_id="citeseq-pbmc", dataset_name="CITE-seq PBMC",
                 modality="rna")
```

| The MuData holds | Written as |
|---|---|
| `mdata.obs`, `mdata.obsm`, `mdata.obsp['connectivities']`, `mdata.uns['cellucid']` | the root: the cells, the shared embeddings and vector fields, the shared graph, the identity and launch defaults (as for an AnnData) |
| `mdata.uns['<column>_colors']`, a modality's `uns['<column>_colors']` | the colors of the root's columns at the root, of a modality's own columns in its `mod/<name>/uns` (`obs_colors` names them `<name>:<column>`) |
| each modality's matrix the viewer opens it with (as for an AnnData: `raw.X` with `raw.var` when `raw` holds more genes than `X` or `X` is scaled, else `X`; without `X`, `raw` or a layer), `obs`, numeric `obsm` arrays, `obsp['connectivities']` | `mod/<name>`, a Cellucid-ready AnnData: gene-major `X`, every array laid out for the browser; a modality written from `raw` makes the root `var` and its `varmap` follow `raw`'s genes |
| `mdata.obsmap`, `mdata.varmap`, the per-modality masks in `mdata.obsm` | as mudata keeps them: which cells and features each modality has |
| `mdata.mod` order | `mod-order` |

The viewer opens a Cellucid-ready MuData with one more select under
**Data mapping**, **Modality** ({doc}`../../web_app/b_data_loading/19_multimodal_mudata`).
`modality=` chooses the one it opens on; without it the viewer opens on the
first gene-expression modality (`rna`, `gex`, `gene expression`,
`gene_expression`, `sct`, any case), else the first. `var_gene_id_column`
names a column of that modality's `var`. Each modality's own `obs` columns and
`obsm` arrays are shown as `<modality>:<name>` (`rna:leiden`,
`rna:X_umap`), missing where the modality did not measure the cell; columns
mudata copied to the root under those names are not written twice. A
modality's embedding can be a launch default (`embedding_defaults={"2d":
"rna:X_umap"}`); `latent_key` names a root embedding or one of a modality that
measured every cell.

The genes of a MuData are each modality's own, so `gene_expression`, `var` and
`gene_identifiers` are refused with a MuData; change the modality first
(`mdata.mod["rna"] = mdata.mod["rna"][:, genes].copy(); mdata.update()`), or
write that one modality with `modalities=`.

`modalities=` chooses the modalities written, in that order (the R package's
`cellucid_prepare(modalities = )` takes the same names). Several are a
Cellucid-ready MuData of those modalities alone. One is a plain Cellucid-ready
AnnData: the root's cells, `obs`, embeddings, graph and `uns['cellucid']`, with
that modality's `X`, `var` and `obsm` arrays placed in the MuData's cell order
-- it must have measured every cell -- and `gene_expression`, `var` and
`gene_identifiers` then apply to it:

```python
prepare("rna.h5ad", format="h5ad", data=mdata, modalities=["rna"],
        gene_identifiers=marker_genes)
prepare("cite.zarr", format="zarr", data=mdata, modalities=["rna", "prot"])
```

The modalities are written one after the other, each matrix converted only
while it is written, so memory beyond the MuData itself holds one modality's
converted matrix.

### What a MuData must be

| `prepare()` refuses | Do instead |
|---|---|
| modalities that share features, not cells (`mdata.axis` 1 or -1) | prepare each modality as an AnnData: `data=mdata.mod["rna"]` |
| a repeated cell name, in the MuData or in one modality | `mdata.obs_names_make_unique(); mdata.update()` (or the modality's `obs_names_make_unique()`) |
| a modality with no cells or no features | `del mdata.mod["<name>"]; mdata.update()` |
| a modality without `X` | give it one: `mdata.mod["<name>"].X = mdata.mod["<name>"].layers["counts"]` |
| no embedding anywhere (root or modality `obsm`) | compute one (`sc.pp.neighbors(mdata.mod["rna"]); sc.tl.umap(mdata.mod["rna"])`) or pass `embeddings=` |
| `format="h5ad"` | `"h5mu"`, `"zarr"` or `"zarr.zip"`; an `.h5ad` holds one AnnData |

## Read it back

```python
mdata = md.read_h5mu("citeseq_cellucid.h5mu")      # md.read_zarr for a folder
adata = ad.read_h5ad("pbmc_cellucid.h5ad")
```

A session saved in the viewer applies to the object it was saved on, a MuData
included: {func}`~cellucid.apply_cellucid_session_to_anndata` checks it against
the cells and the modality its data mapping names, plans the change, and
writes the labels, highlight groups and user-defined fields to `mdata.obs`
and the record of the apply to `mdata.uns`.

## Next

- {doc}`09_cellucid_ready_anndata_specification` (§5: the Cellucid-ready MuData contract)
- {doc}`../g_api_reference_coverage/api/prepare`
