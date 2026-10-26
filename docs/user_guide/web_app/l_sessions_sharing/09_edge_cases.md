# Edge cases: sessions

**Audience:** power users and anyone building repeatable workflows

**Time:** 15–30 minutes
**What you’ll learn:**
- which unusual outcomes are expected under the strict current reader;
- how cancellation, large bundles, and different screens behave; and
- when to create a fresh current session.

## Different identity, even with similar files

A session saved on a file opened from GitHub restores on the same file opened
with **Open…** (or a drop) from your computer, from a web address,
or from a server: the route is not part of a dataset's identity, as long as the
file states its own id in `uns['cellucid']['id']`. A file without one gets an
id from how it was opened (its file name, its link), and then the route does
matter. If a publication changes
its dataset id, `datasetId` differs, and that rejects the complete restore.

This is intentional. Cellucid will not extract the floating layout or skip only
cell-indexed chunks. Load the dataset identity used by the sender, or create a
fresh session after loading it.

See {doc}`07_versioning_compatibility_and_dataset_identity`.

## Same fingerprint, changed content

The fingerprint checks dataset id, cell count, variable count, the
data mapping, and a `cellOrder` record — a digest of the coordinates of the
embedding the live view showed, plus that embedding and its dimension. It is
not a hash of expression values, obs/var fields, category inventories, cell
identifiers, or the embeddings that were not displayed.

So the residual gap is *values*, not order: republishing changed field values,
renamed or re-ordered categories, or a different gene table under the same
identity and the same sizes, while leaving the displayed coordinates untouched,
still passes the guard. Re-ordering cells or recomputing the embedding no longer
does — both change the digest and are refused with their own message.

Never reuse an identity for changed content. Publish the changed file under a
new immutable id and save a new session.

## Eager and lazy scheduling

Large highlight memberships, inactive categorical code columns, and analysis
artifacts are scheduled lazily. This means progress can remain active longer,
not that Cellucid has partly succeeded. The terminal **Session fully restored**
notification appears only after every lazy chunk has applied and the
transaction has committed.

If a late chunk is invalid, earlier camera, fields, filters, highlights,
windows, playback, and cache ownership roll back. If you press **Cancel** or
start another restore, the older progress is dismissed and neither success nor
a red product failure is published.

## Empty state replaces existing state

An ordinary current bundle always carries:

- `cinematic/camera`, even when there are no keyframes; and
- `analysis/windows` and `analysis/docked`, even when no analysis is open.

Restoring those empty payloads clears a destination Camera Path and closes
floating analysis windows the session did not have. It
does not preserve stale state from whatever you explored before loading the
session. Saved playback state is restored exactly; autoplay begins only at a
successful commit when the saved settings require it.

## Analysis results in sessions

Differential expression, marker-gene, signature, correlation, dot plot and
composition results are saved with their analysis (in the sidebar panel or a floating window) and come
back without being run again. A session is still not the canonical
scientific-results archive: export result tables for publication or downstream
analysis. A saved result that names a gene or category the open dataset lacks
is refused.

## Different screens and browsers

Floating panels are clamped so they remain reachable on the current viewport.
A layout saved on a large monitor can therefore move slightly on a laptop or
mobile-sized window. Reposition and save a new device-specific session if exact
panel placement matters.

Restore does not promise canvas focus, hover state, or pointer-lock state.
Click the canvas once before using keyboard shortcuts.

## Large bundles

Sessions grow with:

- many highlight groups containing millions of indices;
- many categorical user-defined fields;
- many saved views; and
- many saved analysis results.

The loader enforces stored and decoded limits, checks gzip output length before
native decompression, and yields during heavy preflight/restore work so Cancel
and newer loads remain responsive. You can keep ordinary workflows lighter by
saving meaningful milestones and deleting obsolete highlights. A highlight
that is a whole category costs a few bytes whatever its size: it is stored
relative to the category.

## Missing field or changed schema

A matching fingerprint does not make a missing field acceptable. If a saved
active field, category inventory, code column, view descriptor, or analysis
result no longer matches the loaded dataset, the exact contributor rejects
the bundle and the transaction rolls back. Cellucid does not quietly keep the
current field or apply the remaining chunks.

Load the unchanged dataset file, or configure the changed dataset and
use **Save** to author a new current bundle.

## Catalog default sessions are not manual sessions

A catalog's default session has an intentionally smaller, integrity-pinned
profile that is applied only when its dataset is opened from that catalog. Downloading it and choosing **Load** is rejected because the manual
reader requires the explicit camera-path chunk, which a default does not carry.

Choose the dataset from the catalog first, then use **Save** if you need a shareable manual
session. See {doc}`04_official_sample_states`.

## Next steps

- {doc}`10_troubleshooting_sessions` (symptom → diagnosis → fix)
- {doc}`12_reference` (exact framing and chunk profiles)
