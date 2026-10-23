# Data preparation (`prepare`)

```{eval-rst}
.. currentmodule:: cellucid
```

This page documents {func}`~cellucid.prepare`, which writes a **Cellucid-ready
AnnData**: one standard AnnData (`.h5ad`, `.zarr` folder, or `.zarr.zip`)
laid out for the web app's range reads and naming itself in `uns['cellucid']`.
The result can be:
- opened in the web app (**Files**, or **Open from link** once hosted),
- served with `cellucid serve <file>`,
- embedded in notebooks with {func}`~cellucid.show_anndata`,
- read by any AnnData tool (`anndata.read_h5ad`, `anndata.read_zarr`).

Long-form walkthroughs: {doc}`../../c_data_preparation_api/index`. The written
contract: {doc}`../../c_data_preparation_api/09_cellucid_ready_anndata_specification`.

---

## Fast path (copy/paste)

```python
from cellucid import prepare

path = prepare(
    "my_study.h5ad",
    format="h5ad",
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
    embeddings,
    dataset_id,
    dataset_name,
    obs=None,
    var=None,
    gene_expression=None,
    connectivities=None,
    vector_fields=None,
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
    zarr_format=None,
    force=False,
) -> Path
```

Every argument after `out_path` is keyword-only. `format` is required and
has no default. The R package's `cellucid_prepare()` takes the same arguments
and writes the same Zarr layout (`format = "zarr"` or `"zarr.zip"`; R writes no
`.h5ad`).

### Which `format`?

| `format` | Writes | Choose it for |
|---|---|---|
| `"zarr"` | a `.zarr` folder of small files | GitHub: no file over 100 MB, no LFS, many small files at once over HTTP/2 |
| `"h5ad"` | one `.h5ad` file | everything else: single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files, `cellucid serve` and Jupyter |
| `"zarr.zip"` | the Zarr store in one `.zarr.zip` file | wherever an `.h5ad` goes, as Zarr |

The measurements behind this are in {ref}`which-format`.

`out_path` must end in the suffix of `format`; a mismatch is refused with the
path to write instead.

| Argument | Stored as |
|---|---|
| `out_path` | the output; ends in `.h5ad`, `.zarr` or `.zarr.zip`, matching `format` |
| `format` | the container: `"h5ad"`, `"zarr"` or `"zarr.zip"` (required) |
| `embeddings` | `obsm[key]`, float32, every column |
| `dataset_id`, `dataset_name` | `uns['cellucid']['id']`, `['name']`; `uns['title']` |
| `obs`, `var` | `obs`, `var` (repeated text as categorical) |
| `gene_expression` | `X`: CSC when sparse, float32 (float64 for counts above 2^24) |
| `connectivities` | `obsp['connectivities']`, symmetric CSR float32 |
| `vector_fields` | `obsm['<field>_<basis>']` for each embedding `X_<basis>` (`velocity_umap`) |
| `gene_identifiers` | keeps those genes, in that order |
| `dataset_description`, `citation_text`/`citation_url`, `source_name`/`source_url` | `uns['cellucid']['description']`, `['citation']`, `['source']` |
| `embedding_defaults` … `centroid_min_points` | `uns['cellucid']['defaults']` (only the ones given) |
| `latent_key` | also `uns['cellucid']['outlierQuantiles']`: every categorical `obs` column's outlier quantiles in that latent space |
| `zarr_format` | 2 (default) or 3, for `format="zarr"` or `"zarr.zip"`; refused with `"h5ad"` |
| `force` | replace an existing output of the same kind |

---

## Interface reference

```{figure} ../../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: Cellucid web app with the sidebar open and a single-cell embedding colored by cell type.
:width: 1440px

A loaded dataset in Cellucid: the sidebar controls the active view while the categorical legend maps directly to the colored points.
```

---

## API reference

```{eval-rst}
.. autofunction:: prepare
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
- **Values.** Embeddings, expression, graph and vectors must be finite and fit
  in float32; `obs`/`var` numbers may hold NaN but not infinities.
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
