# Edge cases (sessions)

**Audience:** power users and anyone building repeatable workflows  
**Time:** 15–35 minutes  
**What you’ll learn:**
- Common edge cases that make sessions restore “imperfectly”
- What behavior is expected vs a bug
- How to design workflows that avoid the worst cases

**Prerequisites:**
- Basic familiarity with sessions ({doc}`01_session_mental_model`)

---

## Dataset identity and mismatch edge cases

### Restoring into a different dataset (mismatch)

**Expected behavior**
- Cellucid warns about a dataset mismatch.
- Dataset-dependent state is skipped (filters/highlights/active fields).
- Only dataset-agnostic layout may restore.

**Why this exists**
- Safety: applying highlights/filters to the wrong dataset is worse than “restore didn’t work”.

**What to do**
- Load the correct dataset first, then load the session again.
- Or create a new session for the dataset you actually want.

See {doc}`07_versioning_compatibility_and_dataset_identity`.

### “Same dataset content, different load method”

Example:
- you saved a session while loading from GitHub,
- your collaborator opens the same export folder locally.

**Expected behavior**
- session may mismatch because the source type and dataset id differ.

Fix:
- align the load method, or
- create a new session under the collaborator’s load method.

---

## Progressive restore edge cases (eager vs lazy)

### “My highlight groups are there, but nothing is highlighted yet”

This can happen when:
- highlight meta restores eagerly, and
- large highlight memberships are still restoring lazily.

Expected behavior:
- memberships appear later (watch the session progress notifications).

If it never resolves, see {doc}`10_troubleshooting_sessions`.

### “Session restore finished, but analysis results are missing”

Expected behavior (dev-phase):
- analysis windows restore (settings + geometry),
- analysis results are not treated as authoritative artifacts and may need recomputation,
- some caches may restore lazily to accelerate recomputation.

---

## UI layout edge cases (different screens, different browsers)

### Floating panels appear “clamped” or moved

Expected behavior:
- when restoring on a different screen size, floating panel geometry is clamped/cascaded so panels remain visible.

What to do:
- manually reposition panels and resave the session if you want “perfect layout” on that machine.

### Keyboard shortcuts / focus feels different after restore

Expected behavior:
- restore does not guarantee canvas focus or pointer-lock state.

Fix:
- click the canvas once, then retry shortcuts.

---

## Large session edge cases (performance + file size)

Sessions can become large when you have:
- huge highlight groups (hundreds of thousands to millions of indices),
- many groups/pages,
- many snapshots,
- many user-defined categorical fields (per-cell codes),
- analysis caches/artifacts.

Symptoms:
- Save State takes longer.
- Load State shows progress for longer.
- Browser may become memory-constrained.

Mitigations:
- save milestone sessions (not every minute),
- keep highlight groups coarse and meaningful,
- reduce snapshots before saving if you don’t need them,
- consider stripping analysis caches (future/advanced; not always available).

---

## Field/feature availability edge cases

### A saved active field no longer exists

If you load a session into an export where:
- the active field was deleted/purged,
- a gene is missing,
- obs/var schema changed,

then active field restore may fall back to a default field.

Expected behavior:
- best-effort restore; missing fields cannot be magically restored.

Fix:
- load the export version that matches the session,
- or update the session (set a new active field and resave).

---

## Auto-restore edge cases (`state-snapshots.json`)

### `state-snapshots.json` exists but is ignored

Common causes:
- it contains no `.cellucid-session` filenames,
- filenames don’t end with `.cellucid-session`,
- it is served as HTML instead of JSON (SPA fallback).

See {doc}`04_auto_restore_latest_from_dataset_exports`.

---

## Next steps

- {doc}`10_troubleshooting_sessions` (the full symptom → diagnosis → fix map)
- {doc}`12_reference` (format and implementation notes)
