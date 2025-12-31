# Reference (session bundles)

**Audience:** power users and developers  
**Time:** 15–40 minutes  
**What you’ll learn:**
- The `.cellucid-session` container format (high level)
- The manifest schema (what metadata is stored)
- Chunk inventory (what chunk IDs you may see)
- The `state-snapshots.json` manifest schema used for auto-restore

---

## File extension

Cellucid session bundles use:

- `.cellucid-session`

They are downloaded from the web app via **Save State** and restored via **Load State**.

---

## Container format (high level)

A `.cellucid-session` file is a single binary container:

1) ASCII MAGIC bytes: `CELLUCID_SESSION\\n`
2) `manifestByteLength` (u32 little-endian)
3) `manifest` JSON bytes (UTF-8)
4) repeated payload chunks:
   - `chunkByteLength` (u32 little-endian)
   - `chunkBytes` (raw stored payload; may be gzip-compressed; decoded based on manifest metadata)

Dev-phase constraints:
- no backward compatibility guarantees
- no migration/version fields
- sessions treated as untrusted input (strict bounds + size guards)

---

## Manifest schema (practical)

The manifest is a JSON object that includes:

- `createdAt`: ISO timestamp
- `dataSource`: a snapshot of the active data source selection (debug/provenance; not used to auto-load the dataset)
- `datasetFingerprint`: used to detect dataset mismatches on restore
- `chunks`: array of chunk metadata entries

Example shape (simplified):

```json
{
  "createdAt": "2025-12-30T02:00:14.000Z",
  "dataSource": { "sourceType": "local-demo", "datasetId": "suo", "userPath": null },
  "datasetFingerprint": { "sourceType": "local-demo", "datasetId": "suo", "cellCount": 12345, "varCount": 20000 },
  "chunks": [
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
  ]
}
```

Important fields per chunk:
- `priority`: `eager` vs `lazy` (progressive restore)
- `kind`: `json` vs `binary`
- `codec`: `none` vs `gzip`
- `datasetDependent`: whether the chunk is skipped on dataset mismatch
- `dependsOn` (optional): dependency ordering hints for contributors

---

## Dataset fingerprint (what “same dataset” means)

The dataset fingerprint is intentionally small and fast to compute.

At minimum it includes:
- `sourceType`
- `datasetId`

It may also include lightweight size guards:
- `cellCount`
- `varCount`

If the fingerprint does not match, dataset-dependent chunks are skipped and only dataset-agnostic layout restores.

Deep dive: {doc}`07_versioning_compatibility_and_dataset_identity`

---

## Chunk inventory (what you may see)

These are the chunk IDs used by the current session system:

| Chunk id | Priority | Dataset dependent | Contains |
|---|---:|---:|---|
| `core/field-overlays` | eager | yes | rename/delete registries + user-defined field definitions (metadata only) |
| `core/state` | eager | yes | camera + UI controls + dimension + filters + active fields + multiview descriptors |
| `ui/dockable-layout` | eager | no | floating non-analysis panels geometry + open/closed |
| `analysis/windows` | eager | yes | open analysis windows descriptors (settings + geometry) |
| `highlights/meta` | eager | yes | highlight pages + group shells (no cell indices) |
| `user-defined/codes/<fieldId>` | eager/lazy | yes | user-defined categorical codes (binary) |
| `highlights/cells/<groupId>` | lazy | yes | highlight group membership indices (binary) |
| `analysis/artifacts/*` | lazy | yes | analysis caches/artifacts (binary) |

Notes:
- the presence and exact contents of chunks are dev-phase and can change across builds
- missing contributors are skipped (best-effort robustness)

---

## `state-snapshots.json` (auto-restore manifest)

Auto-restore reads `state-snapshots.json` from the dataset exports directory and finds entries whose filename ends with `.cellucid-session`.

Supported shapes:

- `{ "states": [ "file.cellucid-session", ... ] }` (recommended)
- `[ "file.cellucid-session", ... ]` (also accepted)

The **last** matching entry is auto-restored on app startup.

Deep dive: {doc}`04_auto_restore_latest_from_dataset_exports`

---

## Implementation pointers (for developers)

Session bundle orchestrator:
- `cellucid/assets/js/app/session/session-serializer.js`

Container format framing:
- `cellucid/assets/js/app/session/bundle/format.js`

What is saved/restored (developer-facing list):
- `cellucid/assets/js/app/state-serializer/README.md`

Design plan:
- `cellucid/markdown/session-serializer-plan.md`
