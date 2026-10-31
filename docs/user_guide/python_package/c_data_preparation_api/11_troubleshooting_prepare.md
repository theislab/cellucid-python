# Troubleshooting `prepare()`

**Audience:** everyone
**Time:** look up your message
**Goal:** symptom → cause → fix for every refusal `prepare()` makes

Every input is checked before anything is written, so a refusal leaves nothing
on disk. The message names the argument, what was found, and usually the fix.

## The output path

| Message | Cause and fix |
|---|---|
| `prepare() needs format=, the kind of output to write. format must be "h5ad", "zarr", "zarr.zip" or "h5mu": ...` | `format` is required. `"zarr"`: a folder of small files, for GitHub (no file over 100 MB, no LFS). `"h5ad"`: one file, for single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files, `cellucid serve` and Jupyter. `"zarr.zip"`: the Zarr store in one file, wherever an `.h5ad` goes. `"h5mu"`: one file for a MuData (`data=mdata`). See {ref}`which-format`. |
| `format must be "h5ad", "zarr", "zarr.zip" or "h5mu": ... Got '...'.` | Pass one of the four exactly, in lower case. |
| `format="zarr" writes a .zarr folder, and out_path 'data.h5ad' ends in .h5ad, which format="h5ad" writes: write 'data.zarr' instead.` | `out_path` ends in the suffix of `format`. Write the path named, or pass the `format` the suffix belongs to. |
| `... already exists; pass force=True to replace it.` | Pass `force=True` to replace an earlier output of the same kind. |
| `... exists and is not a Zarr store, so force=True does not replace it` | A `.zarr` path names a folder holding other things. Choose another path, or remove the folder yourself. |
| `... is a folder, and format="h5ad" writes a file` (or `"zarr.zip"`) | A folder sits where the file would go; remove it yourself. |
| `zarr_format applies to format="zarr" or "zarr.zip", not format="h5ad".` (or `"h5mu"`) | Drop `zarr_format` for an `.h5ad` or `.h5mu`. |

## Identity and text

| Message | Cause and fix |
|---|---|
| `dataset_id ... must be 1-180 ASCII bytes, start with an ASCII letter or digit ...` | Use letters, digits, `.`, `_` and `-` only, e.g. `"pbmc-3k-v1"`. |
| `dataset_name is displayed verbatim, so it must not carry characters that have no glyph ...` | Remove leading/trailing spaces or invisible characters from the text named. |
| `citation_text is required whenever citation_url is supplied.` (or `source_name` / `source_url`) | A link needs its text. |
| `... must be an absolute http(s) address` | Use a full `https://...` address. |

## Shapes and values

| Message | Cause and fix |
|---|---|
| `All embeddings must have one row per cell: ...` / `obs has N rows, but the embeddings have M.` | Inputs come from different cell sets or orders; take every array from one AnnData. |
| `gene_expression has N rows, but the embeddings have M cells. gene_expression is cells x genes; transpose ...` | Pass cells × genes (`matrix.T` for a genes × cells matrix). |
| `... must hold finite numbers within the float32 range. ...: of N values, 3 NaN. First affected positions: ...` | Remove or replace the NaN/infinite values at the positions named, or cast values beyond float32's range. |
| `obs column 'x' holds N infinite values (first at row r)` | Replace infinities with NaN, which is stored as missing. |
| `obs column 'x' mixes text with int values` | `obs["x"] = obs["x"].astype(str)`. |
| `obs column 'x' is a datetime64[ns] column, which has no AnnData encoding Cellucid writes` | Convert it, e.g. `.astype(str)`. |

## Names and keys

| Message | Cause and fix |
|---|---|
| `Embedding 'X/umap' contains a character that cannot appear in a file name on every platform ...` | Rename the key or column; see the path rules in {doc}`02_input_requirements_global`. |
| `... differ only in letter case` | Two keys would be one folder on macOS and Windows; rename one. |
| `obs column '_index' is reserved` | Rename the column. |
| `obs row names must be non-empty text` / `var row names ...` | `obs.index = obs.index.astype(str)`. |
| `The genes have no names: ...` | Give `var` an index of gene ids, or pass `var_gene_id_column`. |
| `var describes the genes of gene_expression, and no gene_expression was supplied.` | Pass `gene_expression`, or drop `var`, `gene_identifiers` and `var_gene_id_column`. |
| `gene_identifiers contains identifiers not found in var: [...]` | Identifiers match the viewer's gene names (`var_gene_id_column`, else the index) after duplicates were made unique. |

## The graph and vector fields

| Message | Cause and fix |
|---|---|
| `Connectivity matrix must be exactly symmetric: a directed neighbour graph ...` | A directed graph (Seurat's `RNA_nn`, Scanpy's `distances`): pass a symmetric one (`obsp['connectivities']`, Seurat's `RNA_snn`), or symmetrize: `(g + g.T) / 2`. A diagonal (each cell its own neighbour, as in `RNA_snn`) is left out, never refused. |
| `Connectivity weights must all be non-negative.` | Pass connectivities, not signed or distance-like values. |
| `Vector field 'velocity' on embedding 'pca': the viewer binds a vector field to the embedding 'X_<basis>' it is named for ('velocity_umap' to 'X_umap'), so the embedding needs an 'X_<basis>' key. ...` | Name the embedding `X_<basis>` (`X_pca`) ({doc}`08_vector_fields_velocity_displacement`). |
| `Vector field 'velocity' on embedding 'X_umap' has 3 columns, and 'X_umap' has 2: a vector field has one component per embedding column.` | Pass an array exactly as wide as its embedding. |
| `Vector field 'T_fwd' on embedding 'X_umap' would be stored as obsm['T_fwd_umap'], which the viewer reads as the field 'T' on 'X_fwd_umap'. Rename the field.` | Another embedding's longer basis claims the key (here `X_fwd_umap`); rename the field ({doc}`09_cellucid_ready_anndata_specification`). |
| `vector_field_default must name one of vector_fields: [...]` | Use the field name you passed (`"velocity"`), not `"velocity_umap"`. |

## Launch defaults

| Message | Cause and fix |
|---|---|
| `embedding_defaults['2d'] names obsm['X_tsne'], which this store does not have.` | Pass that array in `embeddings`, or name one that is. |
| `latent_key must name one of the embeddings: [...]` | Pass the latent space in `embeddings` and name its key. |
| `obs_keys names obs['x'], which this store does not have.` | Name existing `obs` columns, each once. |
| `var_gene_id_column names var['x'], which var does not have.` | Name an existing `var` column. |

## `data=` and MuData

| Message | Cause and fix |
|---|---|
| `data must be an anndata.AnnData or a mudata.MuData, got ...` | Pass the object, or its parts as `embeddings=`, `obs=`, ... |
| `prepare() needs dataset_id: data's uns['cellucid'] has no 'id', so pass dataset_id=.` | The object does not name itself; pass `dataset_id` and `dataset_name`. |
| `data holds no embedding: its obsm has no numeric array.` | Compute one (`sc.pp.neighbors(adata); sc.tl.umap(adata)`) or pass `embeddings=`. |
| `format="h5mu" writes a MuData: ...` / `A MuData is written as format="h5mu", "zarr" or "zarr.zip"; ...` | `.h5mu` is for a MuData, `.h5ad` for an AnnData; for one modality alone pass `data=mdata.mod["rna"]`. |
| `This MuData shares features, not cells (axis=1); ...` | Prepare each modality as an AnnData: `data=mdata.mod["<name>"]`. |
| `... repeat (...). MuData places each modality's cells among the global ones by name, ...` | `mdata.obs_names_make_unique(); mdata.update()` (or the modality's). |
| `Modality '<name>' has no cells` / `has no features` | `del mdata.mod["<name>"]; mdata.update()`. |
| `Modality '<name>' has no X. ...` | Give it one: `mdata.mod["<name>"].X = mdata.mod["<name>"].layers["counts"]`. |
| `This MuData holds no embedding: ...` | Compute one in a modality or at the root, or pass `embeddings=`. |
| `gene_expression describes the genes of one matrix, and a MuData's genes are each modality's own.` (also `var`, `gene_identifiers`) | Change the modality first: `mdata.mod["rna"] = mdata.mod["rna"][:, genes].copy(); mdata.update()`. |
| `modality names '<name>', which this data does not have. Its modalities: [...]` | Name one of `mdata.mod`. |
| `modality names '<name>', and this data has one feature space: ...` | With an AnnData, leave `modality` unset. |
| `modalities names [...], which data does not hold. Its modalities: [...]` | Name modalities of `mdata.mod`. |
| `latent_key names '<mod>:<key>', and modality '<mod>' did not measure every cell, ...` | Name a root embedding, or one of a modality that measured every cell. |

All of it, with examples: {doc}`12_prepare_from_anndata_or_mudata`.

## The file opens but looks wrong

- **The wrong embedding opens first**: set `embedding_defaults`
  ({doc}`03_embeddings_and_coordinates`).
- **Gene search shows ids instead of symbols**: set `var_gene_id_column`
  ({doc}`05_var_gene_metadata`).
- **No velocity overlay**: the view must show the embedding the field was
  written for (`velocity` on `X_umap` is `obsm['velocity_umap']`, drawn over
  `X_umap`). An `obsm['velocity_umap_2d']` written by hand beside a bare
  `X_umap` (no `X_umap_2d`) no longer binds: name it `velocity_umap`, scVelo's
  own name, or pass it through `vector_fields={"velocity": {"X_umap": v}}`.
- **Slow over the internet**: host an `.h5ad` or `.zarr.zip` rather than a
  `.zarr` folder ({doc}`10_performance_and_storage`).

For loading problems after the file is written, see
{doc}`../d_viewing_apis/15_troubleshooting_viewing`.

The same page for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/11_troubleshooting_prepare`.
