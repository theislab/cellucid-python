# Current format and dataset identity

**Audience:** anyone reopening or sharing sessions

**Time:** 20–35 minutes
**What you’ll learn:**
- what Cellucid means by “the same dataset”;
- why the current reader rejects rather than partially restores; and
- how to make a session reproducible across machines.

## The safety rule

A Cellucid session is cell-indexed scientific state. Applying it to the wrong
dataset could make a filter, highlight, categorical code, or cached analysis
look plausible while referring to different cells. Cellucid therefore uses one
strict rule:

> The complete current bundle and the complete current dataset fingerprint
> must validate, or no restore succeeds.

There is no layout-only partial restore, unknown-chunk skip, compatibility
reader, or migration path. If any chunk, field, byte count, or identity check
fails, the transaction rolls back and reports the exact failure.

## Exact fingerprint

The session stores six fields:

- `sourceType`: the loading route, such as `local-demo`, `github-repo`,
  `local-user`, `remote`, or `jupyter`;
- `datasetId`: the route's exact dataset identifier;
- `cellCount`: the number of observations;
- `varCount`: the number of variables;
- `cellOrder`: a record of *which ordering* of those cells the session was saved
  against, with exactly three members; and
- `dataMapping`: which stored arrays the dataset was shown with.

`cellOrder.digest` is 16 lowercase hexadecimal characters computed over every
cell coordinate of one embedding in dataset row order. `cellOrder.embedding`
names that embedding — the embedding the live view showed when the session was
saved — as an exported embedding id (`{"key": "X_umap", "columns": null}`) or,
for an AnnData source, an `obsm` key with the 0-based columns plotted
(`{"key": "X_pca", "columns": [0, 1]}`). `cellOrder.dimension` is its dimension,
exactly 1, 2, or 3. The embedding is part of the record because every embedding
is a separate payload normalized on its own, so a digest only means anything
within one of them.

`dataMapping` names the embeddings that are turned on and each dimension's
default, the vector-field bindings in use, and for an AnnData source also the
neighbour graph, the gene names, the expression matrix, and the latent space.
For a prepared export these are fixed when the export is written, except which
of its embeddings are turned on, which is each dimension's default, which
exported gene-name set names the genes and which of its vector bindings are in
use; for an `.h5ad` or Zarr file opened in the browser, or an object served by
the Python package, they start at the defaults (or the best guess), and every
group can be changed under **Session › Data mapping**.

The four scalars and `cellOrder` must equal the currently loaded dataset. The
data mapping must too, with one exception: its embeddings. A session that turns
on other embeddings, or names other defaults, than the ones on screen restores
them,
before its views, provided this dataset can show them; otherwise it is refused
and names the embedding. This means the same files
loaded through a GitHub repository and through a local folder are intentionally
different identities. Likewise, moving a remote dataset to a different
published identifier requires a newly saved session.

The first four fields say *which* dataset. They cannot say which ordering,
because all four survive a re-export of the same dataset from re-sorted input.
Sessions store selections as row numbers, so without `cellOrder` a republication
under an unchanged identity would silently re-point every lasso, KNN, proximity,
and annotation selection at different cells. That is the failure `cellOrder`
exists to catch.

Restoring is *verify, then apply*: the counts and identity, then the data
mapping, then the digest of the very embedding the session names are all
checked before anything changes on screen. The named embedding is loaded for the
check whichever one is on screen, so a session saved in 2D opens while 3D is
shown.

## Outcome matrix

| Situation | Result |
|---|---|
| Current bundle, exact fingerprint, complete valid chunks | eager and lazy schedules finish, the transaction commits, then Cellucid reports **Session fully restored** |
| Different source type, id, cell count, or variable count | complete rejection; the message names the saved and current cell and gene counts and asks you to open the dataset the session was saved on |
| The session turns on other embeddings, or names other defaults, that this dataset can show | the embeddings are turned on, then the views are restored, then the session's data mapping is applied |
| The session turns on an embedding that this dataset does not have (or cannot plot) | complete rejection; the message names it, for example "This session shows obsm['X_tsne'] columns [0, 1], which this dataset does not have." |
| The session names another graph, gene names, expression matrix, latent space or vector bindings, ones this dataset can use | the choice is restored, then the session |
| The session names a graph, gene-name column, expression matrix, latent space or vector binding this dataset does not have (or cannot use) | complete rejection; the message names it |
| The session was saved with gene names a prepared dataset with one name set does not show | complete rejection; the message names the names on each side |
| The embedding the session recorded its cell order on does not exist in this dataset | complete rejection; the message names the embedding |
| Same id and same cell and gene counts, but the dataset stores its cells in a different order | complete rejection; the message says the coordinates differ from the ones the session was saved against, names both causes a digest cannot tell apart — a different order, or the same cells re-exported from a re-computed embedding — and asks you to load the version the session was saved on, or re-create the selections here |
| A session file saved before Cellucid recorded cell order — a fingerprint with only the four scalars | complete rejection with its own message: such a file can never be confirmed, so re-create the selections and save a new session file |
| A session file saved by an earlier release — a fingerprint with `cellOrder` but no `dataMapping` | complete rejection with its own message: that release read an export format this one does not open, so re-create the view on the current dataset and save a new session file |
| Missing, extra, duplicate, aliased, reordered, or noncanonical chunk | complete rejection before success |
| Corrupt/truncated framing, dishonest byte count, or invalid gzip/JSON/binary payload | complete rejection, bounded decode, rollback |
| Newer load starts or the user presses **Cancel** | older restore aborts and its progress is dismissed without a red product error |
| Bundle made by a different schema/build | accepted only if it is byte-for-byte valid under the one current contract; otherwise rejected |
| Bundle saved in a different browser (for example Safari, opened in Chrome or Firefox) | restores like any other. A saved camera carries its orbit and free-fly forms; restore compares exactly the one the navigation mode drives and ignores the derived one, whose last bits differ between browser engines |

Cellucid does not store a format-version field and does not guess how to
translate another schema. For reproducible work, record the app build shown in
the footer and keep the dataset publication revision with the session.

## Identity is a guard, not a content hash

The fingerprint is still deliberately small. `cellOrder` digests the coordinates
of one embedding; it does not cover expression values, observation or variable
fields, categories, cell identifiers, or the embeddings that are not on screen.
You can therefore still defeat the guard by republishing changed *values* under
the same source type, dataset id, and sizes while leaving the coordinates alone.

What you can no longer do silently is re-order. A dataset re-exported at the
same id from re-sorted input now fails the check instead of restoring selections
that mark the wrong cells.

Treat a dataset id as a content identity:

1. keep cell order, variable order, and prepared files immutable;
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
and dataset id. Everyone should use the **GitHub** loading route, then choose
**Load State**.

### Shared prepared folder

Distribute one unchanged folder and the session together. Everyone should use
the **Prepared** browser picker. Re-zipping without changing files is fine;
re-preparing or renaming the dataset identity is a new publication.

### Shared remote or Jupyter server

Keep the endpoint's dataset identity stable for that exact prepared generation.
If you publish a new generation at the same service, give it a new dataset id
and save a new session.

## Ordinary sessions versus official starting states

Official catalog samples may apply a small SHA-256-pinned state automatically.
That internal five-chunk file is validated through a separate advertised
capability; it is not an ordinary **Load State** artifact.

To share your own state for Suo, Garcia, He, Kanemaru, or Pancreas:

1. choose the official sample;
2. wait for the verified starting view;
3. make any changes you want; and
4. use **Save State**.

The downloaded result has the full current user-session inventory and can be
restored manually against the same official dataset identity. See
{doc}`04_official_sample_states`.

## If restore is rejected

Do not rename the file, edit its manifest, remove a chunk, or retry against a
“similar” dataset. Instead:

1. confirm the expected dataset is fully loaded;
2. confirm the same loading route and exact dataset id;
3. compare the app build in the footer with the sender's recorded build;
4. obtain the original unchanged session again if transfer corruption is
   possible; and
5. ask the sender to load the intended dataset and create a fresh **Save
   State** bundle with the current app.

See {doc}`10_troubleshooting_sessions` for symptom-led checks and
{doc}`12_reference` for the exact container contract.
