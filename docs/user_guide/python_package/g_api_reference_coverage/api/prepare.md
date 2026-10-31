# Data preparation (`prepare`)

```{eval-rst}
.. currentmodule:: cellucid
```

This page documents {func}`~cellucid.prepare`, which writes a **Cellucid-ready
AnnData**, and {func}`~cellucid.restore_original_order`, which reads one back
in its cells' original order. A Cellucid-ready AnnData is one standard AnnData (`.h5ad`, `.zarr` folder, or `.zarr.zip`)
laid out for the web app's range reads and naming itself in `uns['cellucid']`.
Given a MuData it writes a **Cellucid-ready MuData** (`.h5mu`, or MuData Zarr
as a folder or `.zarr.zip`): each modality a Cellucid-ready AnnData.
The result can be:
- opened in the web app (**Open…**, or the **Link or server** bar once hosted),
- served with `cellucid serve <file>`,
- embedded in notebooks with {func}`~cellucid.show_anndata`,
- read by any AnnData tool (`anndata.read_h5ad`, `anndata.read_zarr`).

Long-form walkthroughs: {doc}`../../c_data_preparation_api/index`. The written
contract: {doc}`../../c_data_preparation_api/09_cellucid_ready_anndata_specification`.

---

## Fast path (copy/paste)

```python
from cellucid import prepare

# Everything the AnnData holds (obs, X with var, obsm, obsp['connectivities']):
path = prepare("my_study.h5ad", format="h5ad", data=adata, dataset_id="my-study-v1", dataset_name="My study")

# The same, component by component; an argument given with data= replaces that component.
path = prepare(
    "my_study.h5ad",
    format="h5ad",
    force=True,
    embeddings={"X_umap": adata.obsm["X_umap"], "X_pca": adata.obsm["X_pca"]},
    dataset_id="my-study-v1",
    dataset_name="My study",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp["connectivities"],
    latent_key="X_pca",
)
```

`prepare()` returns the {class}`pathlib.Path` it wrote.

---

## Signature

```text
prepare(
    out_path,
    *,
    format,
    data=None,
    embeddings=None,
    dataset_id=None,
    dataset_name=None,
    obs=None,
    var=None,
    gene_expression=None,
    connectivities=None,
    vector_fields=None,
    modalities=None,
    gene_identifiers=None,
    dataset_description=None,
    citation_text=None,
    citation_url=None,
    source_name=None,
    source_url=None,
    embedding_defaults=None,
    latent_key=None,
    obs_keys=None,
    var_gene_id_column=None,
    vector_field_default=None,
    normalize_embeddings=None,
    centroid_min_points=None,
    modality=None,
    priority_seed=0,
    zarr_format=None,
    force=False,
) -> Path
```

Every argument after `out_path` is keyword-only. `format` is required and
has no default. `embeddings`, `dataset_id` and `dataset_name` are required
unless `data` holds them (`obsm`, `uns['cellucid']['id']`, `['name']`). The R package's `cellucid_prepare()` takes the same arguments
and writes the same Zarr layout (`format = "zarr"` or `"zarr.zip"`; R writes no
`.h5ad` or `.h5mu`, so a MuData from R is always MuData Zarr).

### Which `format`?

| `format` | Writes | Choose it for |
|---|---|---|
| `"zarr"` | a `.zarr` folder of small files | GitHub: no file over 100 MB, no LFS, many small files at once over HTTP/2 |
| `"h5ad"` | one `.h5ad` file | everything else: single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files, `cellucid serve` and Jupyter |
| `"zarr.zip"` | the Zarr store in one `.zarr.zip` file | wherever an `.h5ad` goes, as Zarr |
| `"h5mu"` | one `.h5mu` file (a MuData only) | a MuData wherever an `.h5ad` goes; `"zarr"` and `"zarr.zip"` write a MuData as MuData Zarr |

The measurements behind this are in {ref}`which-format`.

`out_path` must end in the suffix of `format`; a mismatch is refused with the
path to write instead.

| Argument | Stored as |
|---|---|
| `out_path` | the output; ends in `.h5ad`, `.zarr`, `.zarr.zip` or `.h5mu`, matching `format` |
| `format` | the container: `"h5ad"`, `"zarr"`, `"zarr.zip"` or `"h5mu"` (required) |
| `data` | an AnnData or a MuData whose components are the defaults of the arguments below ({doc}`../../c_data_preparation_api/12_prepare_from_anndata_or_mudata`) |
| `embeddings` | `obsm[key]`, float32, every column |
| `dataset_id`, `dataset_name` | `uns['cellucid']['id']`, `['name']`; `uns['title']` |
| `obs`, `var` | `obs`, `var` (repeated text as categorical) |
| `gene_expression` | `X`: CSC when sparse, float32 (float64 for counts above 2^24) |
| `connectivities` | `obsp['connectivities']`, symmetric CSR float32 |
| `vector_fields` | `obsm['<field>_<basis>']` for each embedding `X_<basis>` (`velocity_umap`) |
| `modalities` | with a MuData, the modalities written, in order: several as a MuData of those, one as an AnnData |
| `gene_identifiers` | keeps those genes, in that order |
| `dataset_description`, `citation_text`/`citation_url`, `source_name`/`source_url` | `uns['cellucid']['description']`, `['citation']`, `['source']` |
| `embedding_defaults` … `centroid_min_points`, `modality` | `uns['cellucid']['defaults']` (only the ones given); `modality` (a MuData's) is the modality the viewer opens on |
| `latent_key` | also `uns['cellucid']['outlierQuantiles']`: every categorical `obs` column's outlier quantiles in that latent space |
| `priority_seed` | the cells' order: by ascending priority, the subsampling `uniform` rank under this seed (default 0), recorded in `uns['cellucid']['priority']`, each row's original position in `obs['cellucid_original_position']`, each gene's values below each prefix level in `varm['cellucid_prefix_counts']`; `None` stores the cells as given ({ref}`priority-order`) |
| `zarr_format` | 2 (default) or 3, for `format="zarr"` or `"zarr.zip"`; refused with `"h5ad"` and `"h5mu"` |
| `force` | replace an existing output of the same kind |

---

## Interface reference

```{figure} ../../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: The whole Cellucid window with the Pancreas sample coloured by cell_type: the sidebar shows Coloring and filtering with the categorical field selected and its eight-row legend, and the 3D embedding fills the canvas in the matching colours.
:width: 1440px

A loaded dataset in Cellucid: the sidebar controls the active view while the categorical legend maps directly to the colored points.
```

---

## API reference

```{eval-rst}
.. autofunction:: prepare

.. autofunction:: restore_original_order
```

---

## Edge cases

- **Validation first.** Every input is checked before anything is written; a
  refusal leaves nothing on disk.
- **Atomic placement.** The output is written beside `out_path` under a
  temporary hidden name and moved into place only once complete.
- **Existing outputs.** Replaced only with `force=True`; a folder that is not a
  Zarr store is never replaced.
- **Names.** Repeated cell or gene names are made unique (`MATR3`, `MATR3-1`)
  and reported; keys and columns must be portable path components
  ({doc}`../../c_data_preparation_api/02_input_requirements_global`).
- **Values.** Expression and the graph must be finite and fit in float32;
  embeddings and vectors fit float32, a row of NaN being a cell the array does
  not place (written as given, not drawn in that embedding); `obs`/`var`
  numbers may hold NaN but not infinities.
- **Cell order.** The cells are stored in priority order (`priority_seed`), so
  `anndata.read_h5ad` / `anndata.read_zarr` give them in that order;
  {func}`~cellucid.restore_original_order` puts them back in the order given.
  `uns["cellucid"]["priority"]` records the seed and the cell count (`nObs`)
  it was written for, and is verified before it is used: an object subset,
  combined or reordered since is taken as cells Cellucid did not write, with a
  warning (the specification's "Verified before use"). An object in a
  verified stored order is put back in its original order before it is
  prepared again, whatever `priority_seed` is.
  Every result Cellucid hands back (sessions, selections, exports) names
  cells by their original positions and obs names.
- **Vector fields** bind to an embedding named `X_<basis>`, with exactly as
  many columns as it, under a key that reads back as that field and embedding
  ({doc}`../../c_data_preparation_api/08_vector_fields_velocity_displacement`).

## Troubleshooting

Every refusal, its cause and its fix:
{doc}`../../c_data_preparation_api/11_troubleshooting_prepare`.

## See also

- {doc}`server` (`serve_anndata`, `cellucid serve`)
- {doc}`jupyter` (`show_anndata`)
- {doc}`../../c_data_preparation_api/10_performance_and_storage`
