# Edge cases: sessions

**Audience:** power users and anyone building repeatable workflows

**Time:** 15–30 minutes
**What you’ll learn:**
- which unusual outcomes are expected under the strict current reader;
- how cancellation, large bundles, and different screens behave; and
- when to create a fresh current session.

## Different identity, even with similar files

If a session was saved through GitHub and the same prepared files are later
opened with the local **Prepared** picker, `sourceType` differs. If a remote
publication changes its dataset id, `datasetId` differs. Either case rejects
the complete restore.

This is intentional. Cellucid will not extract the floating layout or skip only
cell-indexed chunks. Load the exact route and dataset identity used by the
sender, or create a fresh session after loading through the intended route.

See {doc}`07_versioning_compatibility_and_dataset_identity`.

## Same lightweight fingerprint, changed content

The fingerprint checks source type, dataset id, cell count, and variable count;
it is not a hash of every cell, gene, field, or coordinate. Reordering cells or
recomputing an embedding under the same identity and sizes can therefore pass
the identity guard while changing scientific meaning.

Never reuse an identity for changed content. Publish a new prepared generation
under a new immutable id and save a new session.

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
- `analysis/cache-inventory`, even when there are no artifacts.

Restoring those empty payloads clears a destination Camera Path or cache. It
does not preserve stale state from whatever you explored before loading the
session. Saved playback state is restored exactly; autoplay begins only at a
successful commit when the saved settings require it.

## Analysis results versus saved cache artifacts

Analysis window descriptors and the exact cache inventory can restore, but a
session is not the canonical scientific-results archive. Export result tables
for publication or downstream analysis. After restore, recompute any result
that the session did not explicitly cache.

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
- cached analysis artifacts.

The loader enforces stored and decoded limits, checks gzip output length before
native decompression, and yields during heavy preflight/restore work so Cancel
and newer loads remain responsive. You can keep ordinary workflows lighter by
saving meaningful milestones, deleting obsolete highlights, and retaining only
the analysis caches you expect collaborators to use.

## Missing field or changed schema

A matching fingerprint does not make a missing field acceptable. If a saved
active field, category inventory, code column, view descriptor, or analysis
artifact no longer matches the loaded dataset, the exact contributor rejects
the bundle and the transaction rolls back. Cellucid does not quietly keep the
current field or apply the remaining chunks.

Load the unchanged prepared generation, or configure the changed dataset and
use **Save State** to author a new current bundle.

## Official sample defaults are not manual sessions

The catalog's `default.cellucid-session` has an intentionally smaller,
integrity-pinned profile that is applied only by the official starting-state
path. Downloading it and choosing **Load State** is rejected because the manual
reader requires explicit camera-path and cache-inventory chunks.

Choose the sample first, then use **Save State** if you need a shareable manual
session. See {doc}`04_official_sample_states`.

## Next steps

- {doc}`10_troubleshooting_sessions` (symptom → diagnosis → fix)
- {doc}`12_reference` (exact framing and chunk profiles)
