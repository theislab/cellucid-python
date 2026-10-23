# What gets saved and restored

**Audience:** computational users, power users, and anyone who has been surprised by “restore didn’t restore X”  
**Time:** 20–40 minutes  
**What you’ll learn:**
- The explicit “saved vs not saved” contract for `.cellucid-session` bundles
- How eager and lazy chunks are ordered inside one atomic restore
- Which parts are dataset-dependent and why identity mismatch rejects everything
- How to verify a complete restore quickly

**Prerequisites:**
- A dataset loaded in the web app
- Optional: a saved `.cellucid-session` bundle you can test with

---

## Quick summary (read this if you only have 60 seconds)

Saved/restored (typical):
- every view's camera + navigation mode, and whether cameras are linked
- every view's embedding
- view layout (live + kept views)
- active fields (obs/var) per view
- filter state per view (category visibility, continuous ranges, outliers)
- centroid points and labels, per view
- highlight pages/groups and their exact memberships
- sidebar control values (generic controls by DOM id), including the
  background, renderer settings and the app theme (light/dark)
- floating panel layout (geometry/open/closed)
- the analysis panel in the sidebar (which analysis is open, each analysis's
  settings and result) and floating analysis windows (settings, result,
  geometry)
- Figure Export settings
- Camera Path settings and keyframes, including an explicitly empty path

Not saved/restored (typical):
- the dataset itself (points, obs, gene expression)
- dataset selection / connection controls (file picker, remote URL, GitHub path)
- Community Annotation state (votes, repo auth)
- Antialiasing (a property of the machine drawing the view)
- Benchmarking state
- ephemeral UI like toasts, in-progress selections, hover state

If the dataset fingerprint differs, the complete operation is rejected and
rolled back; see {doc}`07_versioning_compatibility_and_dataset_identity`.

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
The current bundle contract is closed. A missing required feature owner,
unknown contributor, or noncanonical chunk rejects the restore; Cellucid does
not silently omit that feature.
:::

### 1) Core visualization + UI (“first pixels + UI-ready”) — eager

Restored early so you see the correct view quickly:

- **Camera state**
  - every view's camera position/orbit target and navigation mode
  - whether **Link cameras** is on (views of one embedding share a camera)
  - Camera Path keyframes, timing, interpolation, navigation interaction
    controls, Loop playback, and Autoplay for the matching dataset
  - restored Camera Paths remain stopped unless their saved Autoplay setting
    is enabled; enabled autoplay begins only after the complete restore commits
- **Embeddings**
  - the embedding each view shows, as an embedding reference (an export's
    embedding id, or an `obsm` key with its columns); a view's dimension is
    that embedding's
  - the data mapping: the embeddings that are turned on and the default of
    each dimension, as set under **Session › Data mapping** (turned on before
    the views are restored; see {doc}`../b_data_loading/13_data_mapping`)
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
  - many checkboxes/selects/sliders are captured “by id”, including the
    background, point and lighting settings, and the app theme (light/dark),
    which a restore applies exactly as a choice made by hand
  - accordion open/closed state for many sections
- **Centroids**
  - whether each view shows centroid points and centroid labels

Eager describes ordering, not an intermediate success state. Treat the view as
restored only after the terminal **Session fully restored** notification.

:::{note}
“Captured by id” means the control's DOM identifier is the key — not the label
you read beside it. Renaming a control's visible text does not change what a
saved session restores; changing its id does, and would make older bundles stop
carrying it.

Restoring a control writes its value and then fires the event its owner listens
for, and the owner's own handler is what publishes the value onward. This is why
a restored control and the thing it drives must never disagree: if the owner
rejects the value, the mismatch fails the whole restore instead of leaving the
panel showing something the viewer never accepted.
:::

:::{note}
**`Antialiasing` is the one exception.** It is not
recorded in a session, and loading one **never changes it**, for two reasons.

It describes the machine you are sitting at rather than the analysis you are
doing, and your choice is stored in this browser rather than in the file.

And it has three states while a checkbox has two. The stored preference is
`auto` — the default, and what no stored preference means — or `on`, or `off`.
Under `auto`, Cellucid decides from the size of the dataset you have open and
the dimension of the view you are working in: on below 5,000,000 cells in 3D and
below 200,000 in a 1D or 2D view, off at or above. It re-decides on every
dataset publication, because datasets are switched in place, and when the view
switches between 2D and 3D. Clicking the box ends
automatic selection permanently, in both directions. A captured tick could only
ever replay as `on` or `off`, which would turn “let the dataset decide” into a
choice you never made.

So a session you open — including the starting view every built-in sample
publishes — keeps your choice, or the dataset's choice when you have not made
one, rather than the setting whoever saved the file was using. The status line
under the control tells you which you are on: it reads `Chosen automatically:
<n> cells in <d>D.` while the dataset is deciding, and nothing at all once you have
clicked the box. (It reads `This browser is not providing antialiasing for this
view.` when the device cannot multisample.)

The control itself is live — a change applies to the next frame, and no reload
is involved anywhere — and your click persists across reloads because the
preference is stored in this browser. Nothing else in a session behaves this
way.
:::

The four renderer settings are the easiest instance of this to check yourself,
because you can see both the control and its effect:

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-06-renderer-controls-changed.png
:alt: The renderer settings block with Shader quality set to Light, Level-of-Detail checked, Force LOD level set to 4, and Frustum culling checked.
:width: 480px

Four sidebar controls set to non-default values before **Save** (under **Session file**).
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-08-renderer-controls-restored.png
:alt: The renderer settings block after a restore, showing the same four values as the saved state.
:width: 480px

The same four after **Load** into a freshly started app that began at the
defaults. Walked through step by step in {doc}`03_save_restore_ux`.
```

### 2) Filtering state (per view) — eager

Filters are restored as “modified-only” deltas to keep session files small.

Common things that restore:
- categorical category visibility toggles (legend checkboxes)
- continuous filter ranges (Min/Max sliders)
- continuous color ranges and log-scale toggles (when present)
- outlier filter enable/threshold on categorical fields, restored by getting
  their quantiles (stored in the file, or computed with progress and Cancel);
  a notification names each threshold that could not be restored
- “filter enabled/disabled” state for fields

Because filtering is view-scoped, restore replays the correct filter state into:
- the live view, and
- each snapshot view.

Each view is saved from its own state, whichever view is selected when you
save: saving with a kept view selected still records the live view's own
coloring field and filters as the live view's.

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
  - a categorical field's category centroids are not stored: they are derived
    from its codes and the coordinates of each embedding, and recomputed when
    the codes are restored, so its labels come back exactly where they were

This is one reason sessions can restore an experience even if you created derived fields in the UI.

### 4) Floating panels layout (non-analysis) — eager

If you float/collapse panels, sessions can restore:
- which accordion sections were floated,
- their geometry (position/size),
- open/closed state.

Its manifest profile records `datasetDependent: false`, but the enclosing
bundle still requires an exact dataset fingerprint. Cellucid never extracts
this chunk from a mismatched session.

### 5) Highlights — metadata eager, memberships lazy

Highlights restore in two layers:

1) **Highlight meta (eager)**
   - highlight pages
   - highlight groups “shells” (names, enabled/disabled, etc.)

2) **Highlight memberships (lazy, per group)**
   - actual arrays of cell indices can be large
   - exactly one membership chunk is required for each advertised group

The scheduler yields between lazy chunks so the browser remains responsive,
but public completion waits for all membership arrays. A missing or malformed
group rolls the transaction back.

### 6) User-defined field codes (derived categorical columns) — eager/lazy split

User-defined categorical fields can require storing a code per cell (potentially millions of values).

To keep sessions responsive:
- codes needed for “initial view correctness” can restore eagerly (active coloring fields + snapshot actives),
- other codes restore lazily.

Every categorical user-defined field still requires exactly one codes chunk.
Priority is part of the current contract: a code column is eager exactly when a
target live or snapshot active field needs it, and lazy otherwise.

### 7) Analysis — eager

Sessions persist:
- **the analysis panel in the sidebar**: which analysis is open, and every
  analysis's settings and result
- **open floating analysis windows**: which analysis each shows, its settings,
  its result, and its geometry

Results that are expensive to compute — differential expression, marker
genes, gene signatures, correlation — are stored with the session (genes and
categories by name, numeric columns as their raw values), so they come back
without being run again. A restored result naming a gene or category the open
dataset lacks is refused. Both reopen only after the highlight pages they refer
to have been restored, so an analysis comparing saved pages opens on those
pages even in a fresh tab.

The per-cell gene values analyses read are not stored: they are re-read from
the dataset when needed.

### 8) Figure Export — eager

Every Figure Export setting is saved: plot size and preset, the frame,
background, font and text sizes, title and axis labels, axes, legend and its
position, orientation, depth sorting, reference grid, selection emphasis,
colorblind preview, download format, DPI, SVG strategy, and **Export all
views**. After **Load**, the panel holds the settings the figure was exported with.

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

### Modules and settings intentionally excluded from sessions

- **Antialiasing** — a property of the device (see the note above)
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

## Dataset identity behavior (critical)

On restore, Cellucid compares a **dataset fingerprint** saved in the session with the currently loaded dataset.

If they differ, Cellucid rejects the complete restore before success. Any
registered mutation is rolled back; no layout, filters, fields, highlights, or
analysis subset is committed. A dataset that lacks a field or gene the session
uses is refused the same way, with every missing key named.

This is a safety feature (it prevents applying filters/highlights to the wrong dataset).

See {doc}`07_versioning_compatibility_and_dataset_identity` for:
- what counts as “the same dataset”,
- common ways you accidentally create a mismatch,
- collaboration-safe workflows.

---

## How to verify a restore (practical checklist)

After you load a session, verify in this order:

1) **Did you load the correct dataset?**
   - A mismatch is a terminal rejection. Open the dataset the message names,
     then use the unchanged bundle again.
2) **Did the terminal success appear?**
   - Require **Session fully restored**.
   - A progress card is not success; eager and lazy scheduling are still
     inside the same operation.
3) **State looks right**
   - camera/framing
   - active field
   - filters (Filters panel)
   - snapshot layout
4) **Spot-check exact replacement**
   - an empty saved Camera Path is empty after restore
   - a saved analysis result is shown without running it again
5) **Spot-check derived consequences**
   - highlight counts
   - analysis panel and window presence, settings and results

If any part is wrong, jump to {doc}`10_troubleshooting_sessions`.

---

## Next steps

- {doc}`03_save_restore_ux` (how to use **Save** / **Load** safely)
- {doc}`12_reference` (chunk inventory and implementation notes)
