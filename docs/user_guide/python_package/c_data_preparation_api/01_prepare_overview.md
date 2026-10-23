# `prepare()` overview

**Audience:** everyone who will share, host, or repeatedly open a dataset
**Time:** 10 minutes to a first file
**Goal:** write a Cellucid-ready AnnData and open it

{func}`~cellucid.prepare` takes arrays and data frames — embeddings, cell and
gene metadata, an expression matrix, a neighbour graph, vector fields — checks
every one of them, and writes **one standard AnnData** that the Cellucid web
app reads fast:

| You pass | It is stored as |
|---|---|
| `embeddings={"X_umap": ..., "X_pca": ...}` | `obsm[key]`, float32, every column |
| `obs=` (a DataFrame) | `obs`; repeated text becomes categorical |
| `var=` (a DataFrame) | `var` |
| `gene_expression=` (sparse or dense) | `X`, gene-major (CSC when sparse), float32 |
| `connectivities=` | `obsp['connectivities']`, symmetric CSR float32 |
| `vector_fields={"velocity": {"X_umap": ...}}` | `obsm['velocity_umap']` (`<field>_<basis>` for `X_<basis>`, scVelo's name) |
| `dataset_id`, `dataset_name`, citation, source, launch defaults | `uns['cellucid']` (and `uns['title']`) |

The required `format` names the container, and the output path ends in its
suffix:

| `format` | Writes | Choose it for |
|---|---|---|
| `"zarr"` | a `.zarr` store folder of small files (Zarr format 2 unless `zarr_format=3`) | GitHub: no file over 100 MB, no LFS, many small files at once over HTTP/2 |
| `"h5ad"` | one HDF5 `.h5ad` file | everything else: single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files, `cellucid serve` and Jupyter |
| `"zarr.zip"` | the store folder in one `.zarr.zip` file, entries stored uncompressed | wherever an `.h5ad` goes, as Zarr |

There is no default: a call without `format` is refused with these three
choices. Why each goes where, measured: {ref}`which-format`. Every layout
detail is in {doc}`09_cellucid_ready_anndata_specification`.

Values and coordinates are stored as given (float32), and the browser
normalizes embeddings and computes centroids itself. The one thing computed
here is the outlier filter's quantiles, when `latent_key` names a latent
space: every categorical `obs` column's, stored in
`uns['cellucid']['outlierQuantiles']` ({doc}`03_embeddings_and_coordinates`).

---

## Fast path

### From an AnnData object

```python
import anndata as ad
import cellucid

adata = ad.read_h5ad("pbmc.h5ad")

cellucid.prepare(
    "pbmc_cellucid.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    dataset_id="pbmc-3k",
    dataset_name="PBMC 3k",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
)
```

`prepare()` returns the path it wrote. It writes beside the target under a
temporary hidden name and moves the finished file into place, so a failed or
interrupted call never leaves a half-written file. An existing output is
replaced only with `force=True`.

### Open it

- In Python: `cellucid.serve_anndata("pbmc_cellucid.h5ad")` or, in a notebook,
  `cellucid.show_anndata("pbmc_cellucid.h5ad")`. The file names itself
  (`uns['cellucid']`), so no `dataset_name`/`dataset_id` is passed.
- On the command line: `cellucid serve pbmc_cellucid.h5ad`.
- In the browser: the web app's **Files** button, or **Open from link** for a
  hosted file ({doc}`../../web_app/b_data_loading/index`).

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Cellucid Session panel: the Dataset summary, the Open data group with Sample dataset, Files, Server and Open from link, the Session file group with its Save and Load buttons, and the Data mapping group.
:width: 524px

The Session panel gathers every loading path under **Open data** and keeps the **Save** and **Load** buttons under **Session file**, directly below it.
```

### What success looks like

- The dataset opens under the name you gave it.
- Every `obs` column is listed (or the `obs_keys` you chose, in that order).
- Typing a gene name colours the cells by its expression.
- The embedding you named in `embedding_defaults` is the one shown first.

---

## Required and optional inputs

Required: `out_path`, `embeddings` (at least one), `dataset_id`,
`dataset_name`. Everything else is optional:

| Argument | Without it |
|---|---|
| `obs` | Cells are named `"0"` … `"n - 1"`; no metadata fields. |
| `var`, `gene_expression` | No `X`; gene search is unavailable. `var` without `gene_expression` is refused. |
| `connectivities` | No neighbour graph. |
| `vector_fields` | No velocity overlay. |
| `dataset_description`, `citation_text`/`citation_url`, `source_name`/`source_url` | Not shown. |
| `embedding_defaults`, `latent_key`, `obs_keys`, `var_gene_id_column`, `vector_field_default`, `normalize_embeddings`, `centroid_min_points` | The viewer's own choice; nothing is written to `uns['cellucid']['defaults']`. Without `latent_key` no outlier quantiles are stored, and the outlier filter waits for a latent space chosen in the viewer. |

`dataset_id` is the dataset's stable key — sessions and community annotations
are keyed by it — so keep it when you rewrite the same dataset, and change it
when the cells change. See
{doc}`../b_concepts_mental_models/04_dataset_identity_and_reproducibility`.

## `prepare()` or open the AnnData directly?

| You want to | Use |
|---|---|
| look at an AnnData now, on this machine or in a notebook | {func}`~cellucid.show_anndata` / {func}`~cellucid.serve_anndata` / `cellucid serve`, no preparation |
| share or host a dataset, or open a large one repeatedly | `prepare()` once, then open the result anywhere |
| ship a dataset that opens with a name, citation and chosen defaults | `prepare()` |

A file `prepare()` did not write still opens; only a row-major (CSR)
expression matrix makes the browser build a gene index first, which on a
562k-cell dataset took 4.1 s over a slow link where the prepared file took
0.09 s ({doc}`10_performance_and_storage`).

## Next

- {doc}`02_input_requirements_global`
- {doc}`11_troubleshooting_prepare`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/01_cellucid_prepare_overview`
