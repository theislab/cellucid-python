# Reference: current session bundles

**Audience:** power users and developers  
**Time:** 20–40 minutes
**What you’ll learn:**
- the exact `.cellucid-session` framing and manifest;
- which chunks a current user-authored bundle contains; and
- how ordinary sessions differ from integrity-pinned official sample states.

## File and framing

**Save State** downloads one binary file with the
`.cellucid-session` extension. Its framing is:

1. ASCII magic bytes `CELLUCID_SESSION\n`;
2. `manifestByteLength` as an unsigned 32-bit little-endian integer;
3. exactly that many UTF-8 JSON manifest bytes; and
4. for each manifest entry, a 32-bit little-endian stored length followed by
   exactly that chunk's stored bytes.

There may be no duplicate chunks, unlisted payload, truncated payload, or
trailing bytes. Cellucid treats every bundle as untrusted input and validates
the complete framing, inventory, byte lengths, gzip stream, and dataset
identity before reporting success.

## Exact manifest

The root has exactly `createdAt`, `datasetFingerprint`, and `chunks`:

```json
{
  "createdAt": "2026-07-27T12:00:00.000Z",
  "datasetFingerprint": {
    "sourceType": "local-demo",
    "datasetId": "suo",
    "cellCount": 1000,
    "varCount": 2000
  },
  "chunks": []
}
```

- `createdAt` is a canonical ISO timestamp.
- `datasetFingerprint` has exactly `sourceType`, `datasetId`, `cellCount`, and
  `varCount`.
- `chunks` is the ordered inventory.

Every chunk entry has exactly these nine keys:

```json
{
  "id": "core/state",
  "contributorId": "core-state",
  "priority": "eager",
  "kind": "json",
  "codec": "gzip",
  "label": "Core state",
  "datasetDependent": true,
  "storedBytes": 1234,
  "uncompressedBytes": 5678
}
```

`priority` is either `eager` or `lazy`, `kind` is `json` or `binary`, and
`codec` is `none` or `gzip`. The current chunk object has exactly the nine keys
shown above; an extra key such as `dependsOn` is invalid.

Eager and lazy are deterministic scheduling classes inside one awaited restore:
all eager chunks precede all lazy chunks, and contributor groups retain their
registered order. **Session fully restored** is emitted only after both classes
have validated, applied, and committed.

## Exact dataset identity

All four fingerprint fields must match the currently loaded dataset. A
different source route, dataset id, cell count, or variable count rejects the
whole restore. Cellucid does not apply a layout-only subset and does not skip
dataset-dependent chunks.

This strict equality prevents a session's cell-indexed filters, categorical
codes, highlights, or cached analysis values from being attached to a
different dataset. It is still your responsibility not to reuse the same
identity for reordered or scientifically changed content.

## Current user-authored inventory

Every ordinary **Save State** bundle contains all current singleton chunks,
including explicit empty camera-path and analysis-cache inventories:

| Chunk id | Priority | Dataset dependent | Contains |
|---|---:|---:|---|
| `core/field-overlays` | eager | yes | field rename/delete registries and user-defined field metadata |
| `core/state` | eager | yes | camera, navigation, dimension, filters, active fields, and multiview |
| `ui/dockable-layout` | eager | no | non-analysis floating-panel state and geometry |
| `analysis/windows` | eager | yes | analysis window state and geometry |
| `highlights/meta` | eager | yes | highlight pages and group metadata |
| `analysis/cache-inventory` | eager | yes | exact ordered analysis-artifact IDs, including an empty list |
| `cinematic/camera` | eager | yes | exact keyframes, playback settings, and explicit empty path |

Dynamic binary chunks are closed by those singleton inventories:

| Chunk family | Priority | Required relationship |
|---|---:|---|
| `user-defined/codes/<fieldId>` | eager or lazy | exactly one for every categorical user-defined field; eager exactly when a target live/snapshot active field needs it |
| `highlights/cells/<groupId>` | lazy | exactly one membership payload for every advertised highlight group |
| `analysis/artifacts/bulk-gene/<cacheKey>` | lazy | exact ordered match to `analysis/cache-inventory` |

The built-in singleton profiles and dynamic metadata are exact. Missing,
duplicate, reordered, aliased, unknown, or dishonestly labelled chunks reject
the operation. Empty replacement state is represented by a real singleton
chunk; omission never means “keep whatever is already open.”

## Atomic restore and cancellation

Restore is one transaction:

1. validate framing, manifest, profiles, dataset identity, and declared bounds;
2. capture the small reversible owner state;
3. apply eager and lazy chunks in order;
4. prepare and commit all contributors; and
5. refresh the final UI before publishing success.

Any validation, decode, contributor, commit, or final-refresh failure rolls
back registered state, including categorical code references, camera-path
playback state, panels, analysis windows, and analysis-cache ownership.
Starting a newer restore or pressing **Cancel** aborts the older operation and
dismisses its progress without a false failure or success.

## Official sample state is a separate profile

An official sample may advertise an integrity-pinned
`default.cellucid-session`. Cellucid fetches and applies it automatically only
through the catalog capability described in
{doc}`04_official_sample_states`. Its exact small profile contains these five
gzip JSON chunks, in order:

1. `core/field-overlays`
2. `core/state`
3. `ui/dockable-layout`
4. `analysis/windows`
5. `highlights/meta`

It deliberately contains no cinematic or analysis-cache chunk. Do not download
that internal default and choose **Load State**: the ordinary loader requires
the complete user-authored profile and rejects it. To create a portable manual
session for an official sample, choose the sample, wait for its verified
starting view, then use **Save State**.

## Implementation pointers

- Orchestration and exact profiles:
  `cellucid/assets/js/app/session/session-serializer.js`
- Bundle framing:
  `cellucid/assets/js/app/session/bundle/format.js`
- Feature contributors:
  `cellucid/assets/js/app/session/contributors/`
- Maintainer inventory:
  `cellucid/assets/js/app/state-serializer/README.md`
