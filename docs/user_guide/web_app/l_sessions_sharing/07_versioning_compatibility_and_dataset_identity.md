# Versioning, compatibility, and dataset identity

**Audience:** computational users, power users, and anyone sharing sessions across machines  
**Time:** 25–50 minutes  
**What you’ll learn:**
- Why “same dataset” is not a vague concept in Cellucid (it’s a concrete fingerprint check)
- How sessions behave when the dataset or app version changes
- How to build collaboration workflows that avoid silent mis-restores
- What to do when a restore warns “dataset mismatch”

**Prerequisites:**
- You understand the session mental model ({doc}`01_session_mental_model`)

---

## The big idea

Cellucid sessions are intentionally conservative:

- If Cellucid is confident you loaded the **same dataset**, it will restore dataset-dependent state.
- If it is not confident, it will **skip dataset-dependent chunks** and restore only safe layout.

This is designed to prevent the most dangerous failure mode:

> applying highlights/filters from one dataset onto a different dataset.

---

## “Same dataset” is determined by a dataset fingerprint

Each session bundle stores a small **dataset fingerprint**.

At a high level, the fingerprint includes:

- `sourceType` (how the dataset was loaded)
  - examples: local demo, local folder, GitHub, remote server
- `datasetId` (the dataset identifier within that source)
- lightweight size guards like `cellCount` and `varCount` (when available)

On restore, Cellucid compares:
- fingerprint in the session file vs
- fingerprint of the currently loaded dataset.

If they differ, you get a warning like:

> “Session dataset mismatch (…) Restoring only dataset-agnostic layout.”

And then Cellucid skips dataset-dependent state.

---

## Compatibility matrix (what to expect)

### Case 1: same app build + same dataset fingerprint

Expected outcome:
- full restore (eager + lazy), subject to normal “not saved” exclusions.

### Case 2: same dataset fingerprint, but the dataset content actually changed

This is the dangerous case.

Because the fingerprint is intentionally lightweight, it may not detect “subtle” changes like:
- cell order changes,
- different filtering during export,
- recomputed embeddings with the same shape,
- different gene list ordering with the same length.

If you reuse the same dataset id for changed content, a session can restore “successfully” but be semantically wrong.

**Best practice:**
- Treat dataset ids as content ids.
- If you change anything that would reorder or remove/add cells, treat it as a new dataset id/version.

### Case 3: different dataset fingerprint (mismatch)

Expected outcome:
- dataset-dependent chunks are skipped,
- only dataset-agnostic layout restores (e.g., some floating panel geometry),
- your highlights/filters/active fields won’t apply.

Fix:
- load the correct dataset (same source type + dataset id), then restore again,
- or create a new session for the dataset you actually want to use.

### Case 4: different app build/version

Sessions are dev-phase artifacts:
- no backward compatibility guarantees,
- no migrations,
- chunks are owned by features and may change across builds.

Expected outcome:
- best-effort restore, but any of the following can happen:
  - some chunks are skipped because the contributor doesn’t exist anymore,
  - restore errors for certain chunks (isolated; session may partially restore),
  - a session fails to load with a manifest/format error.

If you need long-term stability:
- treat exported figures and exported data as the archival artifacts,
- and treat sessions as “working state” artifacts.

---

## Why sessions can mismatch across “equivalent” dataset loads

It’s common to think:

> “It’s the same export folder, just loaded a different way.”

But the fingerprint includes the **source type**.

Practical examples:

- If you save a session while loading from **GitHub**, then later load the export folder from your **local disk**, the session can mismatch.
- If you save a session while loading from a **remote server URL**, then later load the same files through a different host/path, the dataset id may differ.

**Collaboration rule:**
When you share sessions, agree on one dataset access method:
- “Everyone loads from GitHub path X”, or
- “Everyone uses the shared export folder Y locally”, or
- “Everyone uses the same hosted export URL”.

If you can’t align access methods, create a session bundle using the same method the recipient will use.

---

## How to design dataset identity so sessions are safe

If you are the person exporting datasets (often the computational user), you control most of this story.

Recommended approach:

1) Assign a stable dataset id for a stable dataset content snapshot.
2) When content meaningfully changes (cell order/selection, embedding recompute, gene list changes), create a new dataset id/version.
3) Keep old exports and sessions if reproducibility matters.

This aligns with:
- Community Annotation (which also keys off dataset identity),
- long-term provenance for figures.

---

## “Fields missing” and “non-restorable state” pitfalls

Even if the dataset fingerprint matches, a restore can still look incomplete if:

- a field referenced by the session is missing in the export,
- a field was renamed/deleted/purged since the session was saved,
- gene expression is unavailable (var fields missing).

Some UI actions can make state intentionally non-restorable:
- permanently purging a field,
- changing field schema or category encoding upstream.

When this happens you may see:
- a default active field instead of the expected one,
- filters silently not applying (because the field doesn’t exist),
- highlight groups present but meaningless (if they were based on missing derived fields).

Fix strategy:
- restore the exact export folder version used to create the session, or
- accept that the session is not portable across those changes and re-create it.

---

## Next steps

- {doc}`05_share_workflows_links_bundles_exports` (collaboration-safe sharing patterns)
- {doc}`10_troubleshooting_sessions` (dataset mismatch and partial restore fixes)
