# What gets saved and restored

**Audience:** computational users, power users, and anyone who has been surprised by “restore didn’t restore X”  
**Time:** 20–40 minutes  
**What you’ll learn:**
- The explicit “saved vs not saved” contract for `.cellucid-session` bundles
- Which parts restore immediately (eager) vs later (lazy)
- Which parts are dataset-dependent (and skipped on dataset mismatch)
- How to verify a restore and diagnose partial restores quickly

**Prerequisites:**
- A dataset loaded in the web app
- Optional: a saved `.cellucid-session` bundle you can test with

---

## Quick summary (read this if you only have 60 seconds)

Saved/restored (typical):
- camera + navigation mode
- view layout (live + snapshots)
- active fields (obs/var) per view
- filter state per view (category visibility, continuous ranges, outliers)
- highlight pages/groups, plus memberships (large groups may restore later)
- some UI control values (generic controls by DOM id)
- floating panel layout (geometry/open/closed)
- analysis windows (settings + geometry), plus some caches (background)

Not saved/restored (typical):
- the dataset itself (points, obs, gene expression)
- dataset selection / connection controls (file picker, remote URL, GitHub path)
- Community Annotation state (votes, repo auth)
- Figure Export UI state (intentionally excluded)
- Benchmarking state
- ephemeral UI like toasts, in-progress selections, hover state

If you see a “dataset mismatch” warning during restore, only dataset-agnostic layout restores; see {doc}`07_versioning_compatibility_and_dataset_identity`.

---

## The contract: “session bundles restore state, not data”

The mental model:

- A session bundle restores **how the app is configured**.
- The dataset loader restores **what data the app is configured on**.

This is intentional because:
- datasets can be huge (sessions should be shareable and fast),
- dataset access methods differ (local folder vs server vs GitHub),
- sessions are treated as untrusted input and should not be a “data injection” path.

---

## Saved + restored (organized by feature area)

This section is the best “source-of-truth” checklist for what you should expect.

:::{note}
Some features restore only if they exist in the current build.
For example, if the Analysis module is not initialized, “analysis windows” cannot restore.
:::

### 1) Core visualization + UI (“first pixels + UI-ready”) — eager

Restored early so you see the correct view quickly:

- **Camera state**
  - camera position/orbit target/navigation mode
  - locked cameras (shared camera) vs unlocked cameras (per view)
- **Dimension levels**
  - which dimension (1D/2D/3D) is active per view
- **Views / snapshots**
  - layout mode (single vs grid)
  - active view selection
  - snapshot descriptors and replay plan (so snapshots can be rebuilt)
- **Active field selection**
  - whether you’re coloring by an obs field or a var field (gene expression)
  - active field keys per view
- **Filter state**
  - per-view filters (see below)
- **Generic sidebar UI controls**
  - many checkboxes/selects/sliders are captured “by id”
  - accordion open/closed state for many sections

Practical expectation:
after eager restore, you should already recognize “where you were”.

### 2) Filtering state (per view) — eager

Filters are restored as “modified-only” deltas to keep session files small.

Common things that restore:
- categorical category visibility toggles (legend checkboxes)
- continuous filter ranges (Min/Max sliders)
- continuous color ranges and log-scale toggles (when present)
- outlier filter enable/threshold (when present)
- “filter enabled/disabled” state for fields

Because filtering is view-scoped, restore replays the correct filter state into:
- the live view, and
- each snapshot view.

### 3) Field overlays and field registry state — eager

Sessions persist “field overlays” that affect what fields exist and how they appear:

- **Rename registry**
  - field display renames
  - category label renames
- **Delete registry**
  - soft-deleted fields
  - permanently purged fields (non-restorable)
- **User-defined fields (metadata)**
  - definitions/provenance for derived fields (but not necessarily the heavy codes yet)

This is one reason sessions can restore an experience even if you created derived fields in the UI.

### 4) Floating panels layout (non-analysis) — eager, dataset-agnostic

If you float/collapse panels, sessions can restore:
- which accordion sections were floated,
- their geometry (position/size),
- open/closed state.

This chunk is considered **dataset-agnostic** and can restore even when the dataset mismatches.

### 5) Highlights — meta eager, memberships lazy

Highlights restore in two layers:

1) **Highlight meta (eager)**
   - highlight pages
   - highlight groups “shells” (names, enabled/disabled, etc.)

2) **Highlight memberships (lazy, per group)**
   - actual arrays of cell indices can be large
   - they restore in the background

Practical implication:
- after eager restore you may see your pages/groups quickly,
- but large group memberships (and any derived counts) may finish restoring later.

### 6) User-defined field codes (derived categorical columns) — eager/lazy split

User-defined categorical fields can require storing a code per cell (potentially millions of values).

To keep sessions responsive:
- codes needed for “initial view correctness” can restore eagerly (active coloring fields + snapshot actives),
- other codes restore lazily.

### 7) Analysis windows + caches — windows eager, caches lazy

Sessions can persist:
- **open floating analysis windows** (which analysis modes were open, their settings, and geometry) — eager
- **some analysis caches/artifacts** to speed up later analysis — lazy

Important nuance:
- analysis **results** are generally not treated as “authoritative” artifacts inside sessions.
- you should expect to recompute results after restore (but caches can make it faster).

---

## Not saved + not restored (intentional exclusions)

This is the “don’t expect it to come back” list.

### Data and data loading

- the dataset itself (points/obs/var/expression matrices)
- dataset selection controls / connection UI
  - local file picker state
  - remote server URL field
  - GitHub repo path field
  - dataset dropdown selection (as a UI action)

### Modules intentionally excluded from sessions

- **Figure Export module UI state**
  - by design: sessions should not “silently change export settings”
- **Benchmarking module state**
- **Community Annotation state**
  - votes, suggestions, moderation
  - GitHub auth/session state

### Ephemeral UI state

- an in-progress selection (candidate set) you haven’t confirmed
- hover state
- toast/notification history
- transient loading spinners

### Low-level rendering buffers

Sessions do **not** store GPU buffers like:
- per-cell color buffers
- transparency buffers
- other derived rendering caches

Those are recomputed from higher-level state after restore.

---

## Dataset mismatch behavior (critical)

On restore, Cellucid compares a **dataset fingerprint** saved in the session with the currently loaded dataset.

If they differ, Cellucid:
- warns about the mismatch, and
- restores only dataset-agnostic chunks (mostly layout).

This is a safety feature (it prevents applying filters/highlights to the wrong dataset).

See {doc}`07_versioning_compatibility_and_dataset_identity` for:
- what counts as “the same dataset”,
- common ways you accidentally create a mismatch,
- collaboration-safe workflows.

---

## How to verify a restore (practical checklist)

After you load a session, verify in this order:

1) **Did you load the correct dataset?**
   - If you saw a dataset mismatch warning, stop and fix that first.
2) **Eager state looks right**
   - camera/framing
   - active field
   - filters (Active filters panel)
   - snapshot layout
3) **Wait for lazy restore (if large)**
   - look for session progress notifications
   - confirm highlight memberships appear
4) **Spot-check “derived” consequences**
   - highlight counts
   - analysis window presence and settings

If any part is wrong, jump to {doc}`10_troubleshooting_sessions`.

---

## Next steps

- {doc}`03_save_restore_ux` (how to use Save State / Load State safely)
- {doc}`04_auto_restore_latest_from_dataset_exports` (how sessions can auto-load from an export folder)
- {doc}`12_reference` (chunk inventory, `state-snapshots.json` schema, and implementation notes)
