# Open a Store in the Cellucid Web App

**Audience:** everyone
**Time:** 5 minutes
**Goal:** see the cells you wrote with `cellucid_prepare()` in the browser.

`cellucid_prepare()` writes a standard AnnData Zarr store. The web app reads it
directly, in your browser: nothing is uploaded, and only the parts a view needs
(an embedding, the obs columns, one gene at a time) are read.

## From your computer: a `.zarr.zip`

A browser folder picker may leave out hidden files, and every Zarr store keeps
its metadata in hidden files (`.zgroup`, `.zattrs`, `.zarray`, `.zmetadata`).
To open a store from your own disk, write it as one file:

```r
library(cellucid)

set.seed(1)
umap <- matrix(rnorm(200), ncol = 2)

cellucid_prepare(
  out_path = file.path(tempdir(), "my_dataset.zarr.zip"),
  format = "zarr.zip",
  embeddings = list(X_umap = umap),
  dataset_id = "my-dataset",
  dataset_name = "My dataset",
  force = TRUE
)
```

Then:

1. Open <https://www.cellucid.com>.
2. In the **Session** panel, under **Files**, click **Zarr ZIP**.
3. Select `my_dataset.zarr.zip`.

The dataset opens under the name you passed as `dataset_name`, with the launch
defaults you set (for example `embedding_defaults` or `obs_keys`) applied.

## From a web address: a `.zarr` folder or a `.zarr.zip`

Once a store is on a web host, paste its address into **Open from link** in the
**Session** panel (or put it in the page's `?data=` parameter):

- a `.zarr` folder: `https://example.org/data/my_dataset.zarr/` — the store's
  consolidated metadata (`.zmetadata`, which `cellucid_prepare()` always writes)
  lets the app open it with one request, and every chunk is read as a view needs
  it;
- a `.zarr.zip` file: `https://example.org/data/my_dataset.zarr.zip` — read in
  place with HTTP range requests.

GitHub, Hugging Face, Zenodo and Dropbox links are accepted too; see
{doc}`02_host_stores_for_sharing` and
{doc}`../../web_app/b_data_loading/14_open_from_link`.

## What success looks like

- the number of cells matches your input;
- the embedding you named in `embedding_defaults` (or the app's own choice)
  is drawn;
- your `obs` columns are listed as fields (in `obs_keys` order, if you set it);
- searching a gene name returns its expression, if you passed `gene_expression`.

If something is missing or the app refuses the store, see
{doc}`03_validate_stores_and_debug_loading` and
{doc}`../i_troubleshooting_index/04_web_app_loading_issues`.
