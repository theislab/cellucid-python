# Current format and dataset identity

**Audience:** anyone reopening or sharing sessions

**Time:** 20–35 minutes
**What you'll learn:**
- which session files this release reads;
- what Cellucid means by "the same dataset";
- why a load onto another dataset is refused whole, while a file from another version of Cellucid is restored part by part; and
- how to make a session reproducible across machines.

## The safety rule

A Cellucid session names cells by their position in the dataset. Applying it
to the wrong dataset could make a filter, highlight, field you made, or saved
analysis result look plausible while referring to different cells. Cellucid
therefore uses one strict rule:

> The whole session file and the whole dataset fingerprint must check out, or
> nothing is loaded.

There is no partial load, no skipping of a part Cellucid does not know, no
reader for an older format, and no migration. If any check fails, what was
shown before stays, and the exact reason is reported.

## Format version

Every session file starts with the signature `CELLUCID_SESSION/4`: the file
kind and the session format version. This release reads format 4 only, and
format 4 holds the **history**: the trackers, their steps, and the step shown
(or, as the person sharing it chose, the current state only, or the steps
only).

A session of another format version is refused before anything is read,
naming its version. One saved by an earlier Cellucid (format 3, which held one
saved state rather than the history of the steps, or formats 1 and 2 before
it) is refused with:

> This session was saved by an earlier version of Cellucid (session format 3), which this version no longer reads: a session file now holds the history of every step (session format 4). It was not opened. Open it in the version of Cellucid that saved it, or re-create the session in this version.

One saved by a newer Cellucid is refused with:

> This session was saved by a newer version of Cellucid (session format 5), which this version cannot read. Open it in that version of Cellucid.

(with its own format number), and any other file with:

> This file is not a Cellucid session file of format 4 (it does not start with CELLUCID_SESSION/4).

The Python package (`cellucid.read_session()`) and the R package
(`cellucid_read_session()`) read the same format and refuse another one with
the same sentences. A file written by a newer Cellucid than the package, with
records the package does not know, is refused as newer, never as damaged:

> This session was written by a newer Cellucid (X) than this package reads (Y); update the cellucid package to read it.

(X the version that wrote it, Y the package's).

A file with the right signature that breaks the format is refused with a
sentence that says what is wrong: a missing or extra key, a chunk listed
twice, or a step referring to a value the file lacks is said as a damaged file,
its **Details** starting *This session file is damaged:*; a file cut short or with bytes after its last chunk says that it
is an invalid session file.

## Exact fingerprint

The session stores six fields:

- `datasetId`: the dataset's identity id (`uns['cellucid']['id']` in the file,
  else the id its route gives it, such as `cellxgene:<dataset id>` for a link,
  or, for a file or folder on your computer, `<format>_<name>_<structure>`:
  its format, its name and eight hexadecimal digits of its cell and gene counts
  and gene names, `h5ad_atlas_1bd42f4e` for an `atlas.h5ad`);
- `name`: the dataset's display name, or `null`. It is not compared; messages
  use it to name the dataset;
- `cellCount`: the number of observations;
- `varCount`: the number of variables;
- `cellOrder`: a record of *which ordering* of those cells the session was saved
  against, with exactly three members; and
- `dataMapping`: which stored arrays the dataset was shown with.

`cellOrder.digest` is 16 lowercase hexadecimal characters computed over the
stored coordinates of one embedding, each cell named by its original position
(its row in the dataset's own order), for the cells the session was saved on.
How those cells were opened does not change it: whole, grown step by step, a
share raised or lowered with **Apply**, or a file that stores its cells in
priority order. `cellOrder.embedding` names that embedding (the embedding the
live view showed when the session was saved) as an `obsm` key with the
0-based columns plotted (`{"key": "X_pca", "columns": [0, 1]}`).
`cellOrder.dimension` is its dimension, exactly 1, 2, or 3. The embedding is
part of the record because a digest only means anything within one
embedding's columns. The digest is specified in
{doc}`../p_developer_docs/10_sessions_persistence_and_serialization`.

`dataMapping` names the embeddings that are turned on and each dimension's
default, the vector-field bindings in use, the neighbor graph, the modality (a
MuData dataset only), the gene names, the expression matrix and its `log1p`
transform, the highly variable genes, and the latent space. They start at the
file's defaults (`uns['cellucid']['defaults']`, or the best guess), however the
file was opened, and every group can be changed under **Session › Data
mapping**; each change is a step of the history.

The three scalars and `cellOrder` must equal the currently loaded dataset. The
data mapping need not equal what is on screen: a session that turns on other
embeddings, names other defaults, or names other choices (graph, modality, gene
names, matrix, …) puts them back, right after the cells and before every view,
provided this dataset can show or use them; otherwise it is refused and names
what is missing. A group the open dataset cannot change must be exactly the one
it is shown with. For a dataset that names its own id, the route is not part of
the identity: the same files opened from a GitHub repository, a web address, a
server or a local folder are the same dataset, and `cellOrder` proves they hold
their cells in the same order. Moving a dataset to a different published
identifier does require a newly saved session.

The first three fields say *which* dataset. They cannot say which ordering,
because all three survive writing the same dataset again from re-sorted input.
Sessions name cells by their original position, so without `cellOrder` a
republication under an unchanged identity would silently re-point every lasso,
KNN, proximity, and annotation selection at different cells. That is the
failure `cellOrder` exists to catch.

Restoring is *verify, then apply*: the cells the session was saved on are
opened, then the counts and identity, the data mapping, and the digest of the
very embedding the session names are all checked before anything changes on
screen. The named embedding is checked whichever one is on screen (its stored
values; it is not loaded for display), so a session saved in 2D opens while 3D
is shown.

## Fields the session uses

Beside the fingerprint, the session lists every observation field and gene the
step it shows uses (`requires`): what its views are colored by, filter on or
have changed the colors of, and the fields its highlight groups are bound to,
by the dataset's own key.
Under the session's own data mapping, and before any part of the session is
put back, Cellucid checks that the open dataset has all of them. A dataset that
lacks any is refused with one message naming every missing key, for example:

> This session uses the observation field "cell_type" and the genes "CD3E", "MS4A1", which the dataset that is open now does not have. The session was not opened. Open the dataset the session was saved on, then load the session again.

A field added to the dataset since the session was saved does not matter.
Genes are compared here when the session names genes the way the open dataset
does (the same gene-name column and modality); otherwise the part that holds a
missing gene names it when it is put back, and the load is refused then. A
field that only another step uses is checked when you show that step.

## Outcome matrix

| Situation | Result |
|---|---|
| A format-4 file, exact fingerprint, every field present, every chunk valid | the step shown is put back, the other steps are read, then Cellucid reports **Session fully restored.** |
| Different cell count or gene count | complete refusal; the message names the saved and current cell and gene counts and asks you to open the dataset the session was saved on |
| Same cell and gene counts, different id | complete refusal; the message says the session was saved on a different dataset with the same number of cells and genes, naming both datasets by the names they show |
| No dataset is open | complete refusal: "No dataset is open. Open the dataset the session was saved on, then load the session again." (with the dataset's name when the session gives one) |
| The session turns on other embeddings, or names other defaults, that this dataset can show | the embeddings are turned on, then the views are put back |
| The session turns on an embedding that this dataset does not have (or cannot plot) | complete refusal; the message names it, for example "This session shows obsm['X_tsne'] columns [0, 1], which this dataset does not have." |
| The session names another graph, modality, gene names, expression matrix, latent space or vector bindings, ones this dataset can use | the choice is put back with the rest of the session |
| The session names a graph, modality, gene-name column, expression matrix, latent space or vector binding this dataset does not have (or cannot use) | complete refusal; the message names it |
| The embedding the session recorded its cell order on does not exist in this dataset | complete refusal; the message names the embedding |
| Same id and same cell and gene counts, but the dataset stores its cells in a different order | complete refusal; the message says the coordinates differ from the ones the session was saved against, names both causes a digest cannot tell apart (a different order, or the same cells written again from a re-computed embedding), and asks you to load the version the session was saved on, or re-create the selections here |
| The dataset lacks an observation field or gene the session uses | complete refusal; the message names every missing key |
| The file holds a part of the state this version of Cellucid does not have | restored part by part: that part is left out of every step and named ("A part of the state this version of Cellucid does not have ("…") is left out: it was saved by a newer version of Cellucid (…)."); everything else is restored |
| One of the file's histories (trackers) has records this version cannot read | restored part by part: that history is left out and named, the others are restored; when it is the one shown, what it showed is kept as a new history of its name |
| A session saved by an earlier Cellucid (session format 3 and earlier) or a newer one | complete refusal before anything is read, naming its format, with the messages above |
| A file that does not start with `CELLUCID_SESSION/4` | complete refusal before anything is read, with the message above |
| Missing, extra, duplicate or reordered chunk; a value or array that does not match its hash; a step tree that does not hold together | complete refusal: "This session file is damaged: …" |
| Corrupt or truncated framing, dishonest byte count, or invalid gzip/JSON | complete refusal, bounded decode, nothing changed |
| A newer load starts or you press **Cancel** | the older load stops, what was shown comes back, and its progress is dismissed without a red error |
| A session saved in a different browser (for example Safari, opened in Chrome or Firefox) | loads like any other. A saved camera carries its orbit and free-fly forms; the load compares the one the navigation mode drives and ignores the derived one, whose last bits differ between browser engines |

Cellucid reads one session format, the one in the signature, and does not
guess how to translate another. For reproducible work, record the app build
shown in the footer and keep the dataset publication revision with the session.

## Identity is a guard, not a content hash

The fingerprint is still deliberately small. `cellOrder` digests the stored
coordinates of one embedding; it does not cover expression values, observation
or variable fields, categories, cell identifiers, or the embeddings that are not
on screen. You can therefore still defeat the guard by republishing changed
*values* under the same dataset id and sizes while leaving the coordinates
alone.

What you cannot do silently is re-order. A dataset written again at the same
id from re-sorted input fails the check instead of loading selections that mark
the wrong cells.

Treat a dataset id as a content identity:

1. keep cell order, variable order, and the published file immutable;
2. publish changed filtering, ordering, embeddings, fields, or categories as a
   new dataset id or immutable repository revision;
3. save a new session after the changed dataset is loaded; and
4. tell collaborators the exact loading route and dataset revision.

## The Python and R readers

The `cellucid` Python package reads the same format and the same six-field
fingerprint, and requires every field: a file without `cellOrder.embedding` or
`dataMapping` is refused. `cellucid.read_session()` reads the trackers,
branches, steps, notes and turns and gives the state at any step;
`apply_cellucid_session_to_anndata(adata, session, tracker=..., branch=..., step=...)`
plans, then writes, one step's labels, highlight groups, fields you made and
gene sets to an `AnnData` (by default the step the file shows, of the tracker
in use). It checks only what that step needs, before anything is written: the
cell count, the dataset id when `uns['cellucid']['id']` names one, the data
mapping (the `obsm`, `obsp`, layer and `var` keys it names must exist) and
`varCount` on the mapped gene axis, the cell order (when the `AnnData` has
the embedding `cellOrder` names, the cell-order digest of its stored values
over the session's cells must be the session's), and the fields the step
uses. The details are in
{doc}`../../python_package/g_api_reference_coverage/api/sessions` and
{doc}`../../python_package/e_jupyter_hooks/11_session_bundles_apply_session_to_anndata`.
R's `cellucid_read_session()` and `cellucid_apply_session()` check the same
({doc}`../../r_package/d_viewing_loading/06_hooks_selection_and_sessions`).

The viewer checks every part of the fingerprint when it loads a session;
an apply checks the parts the chosen step needs, and the cell order only when
the object has the embedding the session names. So a file an apply accepts
may still be refused by the app. That is not an inconsistency to report. See
{doc}`../../python_package/b_concepts_mental_models/04_dataset_identity_and_reproducibility`.

## Reliable collaboration recipes

### Shared GitHub publication

Give collaborators the exact `owner/repository`, branch or immutable revision,
and dataset id. A session saved on it records the link and opens it when it
is loaded; anyone can also open it first with the **Link or server** bar.

### Shared dataset file

Distribute one unchanged file or folder (an `.h5ad`, `.zarr.zip` or `.zarr`
folder, or any other format Cellucid reads) and the session together, with a
stable id in `uns['cellucid']['id']` (`prepare()` and `cellucid_prepare()`
write one). Everyone opens it with **Open…** or by dropping it; how it was
opened does not matter. Re-zipping without changing the store is fine; writing
it again from changed data, or changing its id, is a new publication. A file
without its own id is named by its format, its name, its cell and gene counts
and its gene names: the same data written again (an observation column added)
keeps its id, and two files of one name holding other data are two datasets.

### Shared remote or Jupyter server

Keep the dataset identity the server gives stable for that exact file.
If you publish a new generation at the same service, give it a new dataset id
and save a new session.

## Your sessions versus a catalog's starting views

A dataset opened from a catalog (the official samples included) may put back
a small SHA-256-pinned starting view by itself. That file is named by the
catalog's `session` entry: a format-4 file of one tracker of one root and no
step, without the cells to open or the layout. It goes through the same
fingerprint and field checks as your own sessions; on a subsample of its
dataset the cell-order digest is not compared, because the catalog pins the
file by its SHA-256.

To share your own work on Suo, Garcia, He, Kanemaru, or Pancreas:

1. choose the official sample;
2. wait for its starting view;
3. make any changes you want; and
4. use **Share…**.

The file holds your whole history from the starting view on, and loads onto
the same official dataset anywhere. See {doc}`04_official_sample_states`.

## If a load is refused

Do not rename the file, edit it, or retry against a "similar" dataset.
Instead:

1. confirm the expected dataset is fully loaded;
2. confirm the exact dataset id (and, for a dataset without its own id, the
   same file name and loading route);
3. compare the app build in the footer with the sender's recorded build;
4. obtain the original unchanged session again if the transfer may have
   damaged it; and
5. ask the sender to load the intended dataset and save a fresh file with
   **Share…** in the current app.

See {doc}`10_troubleshooting_sessions` for symptom-led checks and
{doc}`12_reference` for the exact file format.
