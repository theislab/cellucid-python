# Current format and dataset identity

**Audience:** anyone reopening or sharing sessions

**Time:** 20–35 minutes
**What you’ll learn:**
- what Cellucid means by “the same dataset”;
- why the current reader rejects rather than partially restores; and
- how to make a session reproducible across machines.

## The safety rule

A Cellucid session is cell-indexed scientific state. Applying it to the wrong
dataset could make a filter, highlight, categorical code, or saved analysis
result look plausible while referring to different cells. Cellucid therefore uses one
strict rule:

> The complete current bundle and the complete current dataset fingerprint
> must validate, or no restore succeeds.

There is no layout-only partial restore, unknown-chunk skip, compatibility
reader, or migration path. If any chunk, field, byte count, or identity check
fails, the transaction rolls back and reports the exact failure.

## Format version

Every session file starts with the signature `CELLUCID_SESSION/2`: the file
kind and the session format version. This release reads format 2 only. Any
file that does not start with that signature is refused before anything is
read:

> This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2).

## Exact fingerprint

The session stores six fields:

- `datasetId`: the dataset's identity id (`uns['cellucid']['id']` in the file,
  else the id its route gives it, such as `cellxgene:<dataset id>` for a link,
  or `<format>_<name>` for a file or folder, `h5ad_atlas` for `atlas.h5ad`);
- `name`: the dataset's display name, or `null`. It is not compared; messages
  use it to name the dataset (`"<name>" (id "<id>")`);
- `cellCount`: the number of observations;
- `varCount`: the number of variables;
- `cellOrder`: a record of *which ordering* of those cells the session was saved
  against, with exactly three members; and
- `dataMapping`: which stored arrays the dataset was shown with.

`cellOrder.digest` is 16 lowercase hexadecimal characters computed over every
cell coordinate of one embedding in dataset row order. `cellOrder.embedding`
names that embedding — the embedding the live view showed when the session was
saved — as an `obsm` key with the 0-based columns plotted
(`{"key": "X_pca", "columns": [0, 1]}`). `cellOrder.dimension` is its dimension,
exactly 1, 2, or 3. The embedding is part of the record because every embedding
is a separate payload normalized on its own, so a digest only means anything
within one of them.

`dataMapping` names the embeddings that are turned on and each dimension's
default, the vector-field bindings in use, the neighbour graph, the modality (a
MuData dataset only), the gene names, the expression matrix and its `log1p`
transform, the highly variable genes, and the latent space. They start at the file's defaults
(`uns['cellucid']['defaults']`, or the best guess), however the file was
opened, and every group can be changed under **Session › Data mapping**.

The three scalars and `cellOrder` must equal the currently loaded dataset. The
data mapping need not equal what is on screen: a session that turns on other
embeddings, names other defaults, or names other choices (graph, modality, gene
names, matrix, …) restores them, embeddings before its views, provided this
dataset can show or use them; otherwise it is refused and names what is
missing. A group the open dataset cannot change must be exactly the one it is
shown with. For a dataset that names its own id, the route is not part of the identity: the same files
opened from a GitHub repository, a web address, a server or a local folder are
the same dataset, and `cellOrder` proves they hold their cells in the same
order. Moving a dataset to a different published identifier does require a
newly saved session.

The first three fields say *which* dataset. They cannot say which ordering,
because all three survive writing the same dataset again from re-sorted input.
Sessions store selections as row numbers, so without `cellOrder` a republication
under an unchanged identity would silently re-point every lasso, KNN, proximity,
and annotation selection at different cells. That is the failure `cellOrder`
exists to catch.

Restoring is *verify, then apply*: the counts and identity, then the data
mapping, then the digest of the very embedding the session names are all
checked before anything changes on screen. The named embedding is loaded for the
check whichever one is on screen, so a session saved in 2D opens while 3D is
shown.

## Fields the session uses

Beside the fingerprint, the session lists every observation field and gene its
views show or filter on, or a highlight group is bound to (`requires`), by the
dataset's own key. Under the session's own data mapping, and before any part of
the session is applied, Cellucid checks that the open dataset has all of them.
A dataset that lacks any is refused with one message naming every missing key,
for example: *This session uses the observation field "cell_type" and the genes
"CD3E", "MS4A1", which the dataset that is open now does not have.* A field
added to the dataset since the session was saved does not matter.

## Outcome matrix

| Situation | Result |
|---|---|
| Current bundle, exact fingerprint, complete valid chunks | eager and lazy schedules finish, the transaction commits, then Cellucid reports **Session fully restored** |
| Different id, cell count, or variable count | complete rejection; the message names the saved and current cell and gene counts and asks you to open the dataset the session was saved on |
| The session turns on other embeddings, or names other defaults, that this dataset can show | the embeddings are turned on, then the views are restored, then the session's data mapping is applied |
| The session turns on an embedding that this dataset does not have (or cannot plot) | complete rejection; the message names it, for example "This session shows obsm['X_tsne'] columns [0, 1], which this dataset does not have." |
| The session names another graph, modality, gene names, expression matrix, latent space or vector bindings, ones this dataset can use | the choice is restored, then the session |
| The session names a graph, modality, gene-name column, expression matrix, latent space or vector binding this dataset does not have (or cannot use) | complete rejection; the message names it |
| The embedding the session recorded its cell order on does not exist in this dataset | complete rejection; the message names the embedding |
| Same id and same cell and gene counts, but the dataset stores its cells in a different order | complete rejection; the message says the coordinates differ from the ones the session was saved against, names both causes a digest cannot tell apart — a different order, or the same cells written again from a re-computed embedding — and asks you to load the version the session was saved on, or re-create the selections here |
| The dataset lacks an observation field or gene the session uses | complete rejection; the message names every missing key |
| A highlight saved as a whole category, or part of one, whose category no longer holds the same cells | complete rejection: "The category this selection was saved relative to no longer holds the same cells." |
| A file that does not start with `CELLUCID_SESSION/2` | complete rejection before anything is read, with the message above |
| Missing, extra, duplicate, aliased, reordered, or noncanonical chunk | complete rejection before success |
| Corrupt/truncated framing, dishonest byte count, or invalid gzip/JSON/binary payload | complete rejection, bounded decode, rollback |
| Newer load starts or the user presses **Cancel** | older restore aborts and its progress is dismissed without a red product error |
| Bundle made by a different schema/build | accepted only if it is byte-for-byte valid under the one current contract; otherwise rejected |
| Bundle saved in a different browser (for example Safari, opened in Chrome or Firefox) | restores like any other. A saved camera carries its orbit and free-fly forms; restore compares exactly the one the navigation mode drives and ignores the derived one, whose last bits differ between browser engines |

Cellucid reads one session format, the one in the signature, and does not
guess how to translate another. For reproducible work, record the app build
shown in the footer and keep the dataset publication revision with the session.

## Identity is a guard, not a content hash

The fingerprint is still deliberately small. `cellOrder` digests the coordinates
of one embedding; it does not cover expression values, observation or variable
fields, categories, cell identifiers, or the embeddings that are not on screen.
You can therefore still defeat the guard by republishing changed *values* under
the same dataset id and sizes while leaving the coordinates alone.

What you can no longer do silently is re-order. A dataset written again at the
same id from re-sorted input now fails the check instead of restoring selections
that mark the wrong cells.

Treat a dataset id as a content identity:

1. keep cell order, variable order, and the published file immutable;
2. publish changed filtering, ordering, embeddings, fields, or categories as a
   new dataset id or immutable repository revision;
3. save a new session after the changed dataset is loaded; and
4. tell collaborators the exact loading route and dataset revision.

## The Python reader

The `cellucid` Python package reads the same six-field fingerprint and requires
every field: a file without `cellOrder.embedding` or `dataMapping` is refused.
It checks the data mapping against the `AnnData` it applies the session to — the
`obsm`, `obsp`, layer and `var` keys it names must exist — and compares
`varCount` on the mapped gene axis. The exact shapes are in
{doc}`../../python_package/e_jupyter_hooks/11_session_bundles_apply_session_to_anndata`.

Python validates `cellOrder`'s shape and never its digest: the digest is taken
over coordinates the viewer drew, and `apply_cellucid_session_to_anndata()`
targets an `AnnData` whose row order the caller controls, so there is nothing on
that side to check it against. Enforcing cell-order identity is the web app's
job.

The practical consequence is one-directional: a bundle Python accepts may still
be refused by the app. That is not an inconsistency to report. See
{doc}`../../python_package/b_concepts_mental_models/04_dataset_identity_and_reproducibility`.

## Reliable collaboration recipes

### Shared GitHub publication

Give collaborators the exact `owner/repository`, branch or immutable revision,
and dataset id. Everyone should open it with the **Link or server** bar, then
choose **Load**.

### Shared dataset file

Distribute one unchanged file or folder (an `.h5ad`, `.zarr.zip` or `.zarr`
folder, or any other format Cellucid reads) and the session together, with a stable id in `uns['cellucid']['id']`. Everyone opens
it with **Open…** or by dropping it; how it was opened does
not matter. Re-zipping without changing the store is fine; writing it again
from changed data, or changing its id, is a new publication.

### Shared remote or Jupyter server

Keep the dataset identity the server gives stable for that exact file.
If you publish a new generation at the same service, give it a new dataset id
and save a new session.

## Ordinary sessions versus catalog default sessions

A dataset opened from a catalog — the official samples included — may apply a
small SHA-256-pinned default session automatically. That seven-chunk file is
named by the catalog's `session` entry and validated on its own path; it is not
a file for an ordinary **Load**.

To share your own state for Suo, Garcia, He, Kanemaru, or Pancreas:

1. choose the official sample;
2. wait for the verified default view;
3. make any changes you want; and
4. use **Save**.

The downloaded result has the full current user-session inventory and can be
restored manually against the same official dataset identity. See
{doc}`04_official_sample_states`.

## If restore is rejected

Do not rename the file, edit its manifest, remove a chunk, or retry against a
“similar” dataset. Instead:

1. confirm the expected dataset is fully loaded;
2. confirm the exact dataset id (and, for a dataset without its own id, the
   same loading route);
3. compare the app build in the footer with the sender's recorded build;
4. obtain the original unchanged session again if transfer corruption is
   possible; and
5. ask the sender to load the intended dataset and create a fresh bundle with
   **Save** in the current app.

See {doc}`10_troubleshooting_sessions` for symptom-led checks and
{doc}`12_reference` for the exact container contract.
