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

The session stores five fields:

- `sourceType`: the loading route, such as `local-demo`, `github-repo`,
  `local-user`, `remote`, or `jupyter`;
- `datasetId`: the route's exact dataset identifier;
- `cellCount`: the number of observations;
- `varCount`: the number of variables; and
- `cellOrder`: a record of *which ordering* of those cells the session was saved
  against, with exactly two members.

`cellOrder.digest` is 16 lowercase hexadecimal characters computed over every
cell coordinate of the displayed embedding in dataset row order.
`cellOrder.dimension` is the embedding dimension that was on screen when the
session was saved, exactly 1, 2, or 3; it is part of the record because the 1D,
2D, and 3D embeddings are separate files normalized independently, so a digest
only means anything within one of them.

All five must equal the currently loaded dataset — six compared values, since
`cellOrder` contributes two. This means the same files loaded through a GitHub
repository and through a local folder are intentionally different identities.
Likewise, moving a remote dataset to a different published identifier requires a
newly saved session.

The first four fields say *which* dataset. They cannot say which ordering,
because all four survive a re-export of the same dataset from re-sorted input.
Sessions store selections as row numbers, so without `cellOrder` a republication
under an unchanged identity would silently re-point every lasso, KNN, proximity,
and annotation selection at different cells. That is the failure `cellOrder`
exists to catch.

## Outcome matrix

| Situation | Result |
|---|---|
| Current bundle, exact fingerprint, complete valid chunks | eager and lazy schedules finish, the transaction commits, then Cellucid reports **Session fully restored** |
| Different source type, id, cell count, or variable count | complete rejection; the message names the saved and current cell and gene counts and asks you to open the dataset the session was saved on |
| A different embedding dimension is on screen than the one the session was saved in | complete rejection; the message names both dimensions and asks you to switch back, then load again. This is the recoverable case: nothing is wrong with the data |
| Same id and same cell and gene counts, but the dataset stores its cells in a different order | complete rejection; the message says every saved selection would mark the wrong cells and asks you to load the version the session was saved on, or re-create the selections here |
| A session file saved before Cellucid recorded cell order — a fingerprint with only the four scalars | complete rejection with its own message: such a file can never be confirmed, so re-create the selections and save a new session file |
| Missing, extra, duplicate, aliased, reordered, or noncanonical chunk | complete rejection before success |
| Corrupt/truncated framing, dishonest byte count, or invalid gzip/JSON/binary payload | complete rejection, bounded decode, rollback |
| Newer load starts or the user presses **Cancel** | older restore aborts and its progress is dismissed without a red product error |
| Bundle made by a different schema/build | accepted only if it is byte-for-byte valid under the one current contract; otherwise rejected |

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

## The Python reader is deliberately more permissive

The `cellucid` Python package reads the same `.cellucid-session` files, and it
accepts a fingerprint **with or without** `cellOrder`, validating only the
record's shape when it is present. It cannot do more: the digest is taken over
exported coordinates, and
`apply_cellucid_session_to_anndata()` targets an `AnnData` whose row order the
caller controls, so there is nothing on that side to check the digest against.
Enforcing cell-order identity is the web app's job.

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

To share your own state for Suo, Garcia, He, or Pancreas:

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
