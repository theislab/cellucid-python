# Sessions: persistence and serialization

This page documents Cellucid’s **session bundle** system: what gets saved, how it is encoded, and how restore is staged to get “first pixels” quickly while heavier artifacts load in the background.

It is written for:
- contributors changing state fields/UI controls,
- anyone debugging “why didn’t my session restore X?”,
- anyone hosting `.cellucid-session` bundles alongside dataset exports.

## At a glance

**Audience**
- Wet lab / non-technical: read “What sessions are (and aren’t)” + “Troubleshooting”.
- Computational users: read “Coverage: what is saved” + “Dataset mismatch behavior”.
- Developers: read everything; this is a common source of subtle bugs.

**Time**
- 30–60 minutes

**Prerequisites**
- {doc}`05_app_architecture_overview`
- {doc}`06_state_datastate_and_events` (for understanding what state is)

---

## What sessions are (and aren’t)

### Sessions are…

- A **portable snapshot** of UI + app state for a dataset.
- Downloadable as a single `.cellucid-session` file.
- Restorable later to recover:
  - camera/view layout,
  - active fields and filters,
  - highlights and highlight pages (and optionally heavy artifacts),
  - some UI control state.

### Sessions are NOT…

- The dataset itself (no data files are embedded).
- A guarantee of long-term backward compatibility (dev-phase constraint).
- A copy of network/auth state (community annotation tokens are not included).

---

## Where the session system lives (code map)

Session orchestrator:
- `cellucid/assets/js/app/session/session-serializer.js`

Bundle container framing:
- `cellucid/assets/js/app/session/bundle/format.js`

Feature contributors (what chunks are emitted/restored):
- `cellucid/assets/js/app/session/contributors/`

Helpers for capturing/restoring state:
- `cellucid/assets/js/app/state-serializer/`
- `cellucid/assets/js/app/state-serializer/README.md` (detailed coverage list)

UI wiring:
- `cellucid/assets/js/app/ui/modules/session-controls.js`

---

## Bundle format (how `.cellucid-session` is encoded)

Cellucid session bundles are a single binary container with:

1) **MAGIC** bytes (ASCII): `CELLUCID_SESSION\n`
2) `manifestByteLength` (u32 little-endian)
3) manifest JSON bytes (UTF‑8)
4) repeated chunks:
   - `chunkByteLength` (u32 little-endian)
   - `chunkBytes...`

Code:
- `cellucid/assets/js/app/session/bundle/format.js`

### Chunk metadata (what the manifest contains)

The manifest contains a `chunks` array.
Each chunk describes:

- `id`: unique identifier
- `contributorId`: which contributor handles it
- `priority`: `eager` or `lazy`
- `kind`: `json` or `binary`
- `codec`: `none` or `gzip`
- `label`: display label for UI/progress reporting
- `datasetDependent`: whether it should be skipped on dataset mismatch
- `storedBytes`: size after codec (on disk)
- `uncompressedBytes`: expected decoded size (guard)
- `dependsOn` (optional): chunk ordering constraints

Design constraints (dev phase):
- No version fields, no migrations.
- Session files are treated as **untrusted input** (strict bounds checks; size guards).

---

## Restore staging (eager vs lazy)

Restoring a session is intentionally staged:

- **Eager chunks** restore “first pixels + UI-ready” state quickly:
  - camera, layout, active fields, filters, core UI controls.
- **Lazy chunks** restore heavier state in the background:
  - highlight memberships (large per-cell index lists),
  - user-defined field codes (potentially large),
  - analysis caches/artifacts.

Why this matters:
- A session restore should feel fast even for large datasets.
- Heavy artifacts can be cancelled if the user switches datasets or closes the tab.

Implementation:
- `cellucid/assets/js/app/session/session-serializer.js` processes eager chunks first, then schedules lazy processing with yielding between chunks.

---

## Coverage (what is saved/restored)

The authoritative list is in:
- `cellucid/assets/js/app/state-serializer/README.md`

High-level summary:

### Core visualization + UI (“first pixels”)

Saved/restored (eager):
- camera state
  - locked cameras: one global camera
  - unlocked cameras: per-view cameras
- live vs snapshot view layout
- active view id
- per-view dimension levels
- active field selection (obs/var) per view
- active filters (modified-only)
- generic sidebar controls (by DOM id), with explicit exclusions
- floating panels layout (non-analysis)

### Registries and overlays

Saved/restored (eager):
- field rename registry
- field delete/purge registry
- user-defined fields metadata (but large codes may be lazy)

### Highlights and analysis artifacts (heavier)

Saved/restored (often lazy):
- highlight pages/groups and memberships
- user-defined field codes (if large)
- analysis caches/artifacts

### Exclusions (intentional)

Not saved/restored:
- the dataset itself
- dataset selection/connection UI inputs
- figure export UI state
- benchmark UI state
- community annotation state (auth, votes, drafts, moderation UI)
- DOM/WebGL runtime objects (only declarative state is stored)

---

## Auto-load on startup from dataset exports

Cellucid supports a legacy-but-useful workflow:
- datasets can ship “latest session snapshot” pointers in their exports folder

Mechanism:
- On startup, `main.js` calls `sessionSerializer.restoreLatestFromDatasetExports()`.
- The serializer reads `state-snapshots.json` in the current dataset base URL and loads the last `.cellucid-session` entry.

Expected files:
- `<dataset_base_url>/state-snapshots.json`
- `<dataset_base_url>/cellucid-session-....cellucid-session`

Supported `state-snapshots.json` shapes (dev-phase):
- `{ "states": ["file.cellucid-session", ...] }` (recommended)
- `["file.cellucid-session", ...]` (also accepted)

Security model:
- Still treated as untrusted input (bounds checks apply).

---

## Dataset mismatch behavior

Sessions contain a dataset fingerprint.
If the current dataset does not match the session’s dataset:

- dataset-dependent chunks are skipped (highlights/caches/core state)
- dataset-agnostic layout can still restore (e.g. floating panels)

This is a safety feature:
- restoring highlights/caches against the wrong dataset can silently corrupt user interpretation.

---

## Adding a new persisted feature (developer workflow)

If you add a new feature and want it to persist:

1) Decide what kind of state it is:
   - “first pixels” (eager, small, required for UI coherence)
   - “heavy artifact” (lazy, large, optional)

2) Add or extend a contributor in:
   - `cellucid/assets/js/app/session/contributors/`

3) Ensure the serialized payload is:
   - bounded in size (guarded)
   - JSON-safe (or encode as binary)
   - deterministic (avoid including transient runtime-only values)

4) Update exclusions if needed:
   - UI controls: `data-state-serializer-skip="true"` (DOM-level exclusion)
   - Analysis windows: handled by analysis module exclusions

5) Update docs and add a troubleshooting entry (symptoms users will see).

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Load session does nothing (no UI change)”

Likely causes (ordered):
1) User cancelled file picker.
2) Session file is corrupt or truncated.
3) Restore was cancelled because a dataset reload started.

How to confirm:
- DevTools → Console: look for session serializer warnings/errors.
- DevTools → Network: if auto-load, check whether `state-snapshots.json` and the `.cellucid-session` file were fetched.

Fix:
- Try a known-good small session file.
- If hosting sessions, ensure the server serves the binary file correctly (no HTML fallback).

### Symptom: “Session loads, but highlights are missing”

Likely causes:
- Highlights are stored in lazy chunks and are still loading.
- Dataset mismatch caused dataset-dependent highlight chunks to be skipped.

How to confirm:
- Watch notifications (lazy chunk progress).
- Check console logs for “dataset mismatch” skips.

Fix:
- Ensure you loaded the same dataset id/fingerprint.
- Wait for lazy chunks (or inspect which chunks are present in the manifest).

### Symptom: “Invalid chunk length … (session file truncated?)”

Likely cause:
- The `.cellucid-session` response is incomplete/corrupt.

Fix:
- Re-download the session.
- If serving from a host/proxy, ensure it is not truncating large binary responses.

---

Next: {doc}`11_analysis_architecture` (analysis uses sessions for reopening windows and optionally restoring artifacts).
