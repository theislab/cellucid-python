# Web App Loading Issues

**Audience:** everyone  
**Time:** 10 minutes  
**What you'll get:** why the browser will not open a store that exists

Use this page when your `.zarr` folder or `.zarr.zip` file exists but Cellucid
will not open it, or opens it wrongly. The web app's own guides start at
{doc}`../../web_app/index`.

## Symptom: “the folder I picked does not hold a Zarr store”

**Likely cause**
- The browser's folder picker can leave out hidden files, and a Zarr store's
  metadata files (`.zgroup`, `.zattrs`, `.zarray`, `.zmetadata`) start with a
  dot.
- You selected the folder above the store, or a folder inside it.

**Fix**
- Select the `.zarr` folder itself.
- Or write the store as one file, `cellucid_prepare(..., format = "zarr.zip")` with an
  `out_path` ending in `.zarr.zip`, and open that file.

## Symptom: “a link to my store is refused or never loads”

**Likely causes**
- The host does not allow cross-origin requests (CORS) or HTTP range requests.
- The link points at a folder listing page, not at the store or file itself.

**Fix**
- Follow {doc}`../d_viewing_loading/02_host_stores_for_sharing` and
  {doc}`../../web_app/b_data_loading/14_open_from_link`, which lists the hosts
  and link forms that work.
- A `.zarr` folder on a web host needs its `.zmetadata`, which
  `cellucid_prepare()` always writes; upload the hidden files too.

## Symptom: “the store opens but shows the wrong embedding, genes or fields”

**Likely causes**
- No launch defaults were stored, so the viewer guessed.
- Cells were in different orders in the inputs (only caught when the inputs
  carry row names).

**Fix**
- Pass `embedding_defaults`, `latent_key`, `var_gene_id_column`, `obs_keys` as
  needed; or change the choice in the viewer under **Session › Data mapping**
  ({doc}`../../web_app/b_data_loading/13_data_mapping`).
- Align every input to one cell order and write the store again; see
  {doc}`../c_data_preparation_api/02_input_requirements_global`.

## Symptom: “the vector field overlay is not offered”

The viewer binds `obsm['<field>_<basis>']` to the embedding `X_<basis>` it is
named for (scVelo's convention: `velocity_umap` to `X_umap`) when both have the
same number of columns, and draws it only while a view shows that embedding.
Show that embedding, or bind the field under **Session › Data mapping › Vector
fields**. An `obsm[['velocity_umap_2d']]` beside a bare `X_umap` (no
`X_umap_2d`) no longer binds: name it `velocity_umap`, or write it with
`vector_fields = list(velocity = list(X_umap = v))`. See
{doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`.

## Symptom: “the canvas is blank / WebGL context lost”

**Likely causes**
- GPU/WebGL problems, or a dataset too large for the GPU's memory.

**Fix**
- Confirm hardware acceleration and WebGL2 are enabled, and update the browser
  and GPU driver.
- Try a smaller subset of cells to confirm the store itself is fine.

## Next steps

- Validate the store: {doc}`../d_viewing_loading/03_validate_stores_and_debug_loading`
