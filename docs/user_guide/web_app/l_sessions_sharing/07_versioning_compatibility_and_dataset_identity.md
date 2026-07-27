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

The session stores four fields:

- `sourceType`: the loading route, such as `local-demo`, `github-repo`,
  `local-user`, `remote`, or `jupyter`;
- `datasetId`: the route's exact dataset identifier;
- `cellCount`: the number of observations; and
- `varCount`: the number of variables.

All four must equal the currently loaded dataset. This means the same files
loaded through a GitHub repository and through a local folder are intentionally
different identities. Likewise, moving a remote dataset to a different
published identifier requires a newly saved session.

## Outcome matrix

| Situation | Result |
|---|---|
| Current bundle, exact fingerprint, complete valid chunks | eager and lazy schedules finish, the transaction commits, then Cellucid reports **Session fully restored** |
| Different source type, id, cell count, or variable count | complete rejection; current state remains or is restored by rollback |
| Missing, extra, duplicate, aliased, reordered, or noncanonical chunk | complete rejection before success |
| Corrupt/truncated framing, dishonest byte count, or invalid gzip/JSON/binary payload | complete rejection, bounded decode, rollback |
| Newer load starts or the user presses **Cancel** | older restore aborts and its progress is dismissed without a red product error |
| Bundle made by a different schema/build | accepted only if it is byte-for-byte valid under the one current contract; otherwise rejected |

Cellucid does not store a format-version field and does not guess how to
translate another schema. For reproducible work, record the app build shown in
the footer and keep the dataset publication revision with the session.

## Identity is a guard, not a content hash

The fingerprint is deliberately small. It does not hash every coordinate,
field, category, or cell identifier. You can defeat the guard by publishing
changed content under the same source type, dataset id, and sizes.

Treat a dataset id as a content identity:

1. keep cell order, variable order, and prepared files immutable;
2. publish changed filtering, ordering, embeddings, fields, or categories as a
   new dataset id or immutable repository revision;
3. save a new session after the changed dataset is loaded; and
4. tell collaborators the exact loading route and dataset revision.

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
