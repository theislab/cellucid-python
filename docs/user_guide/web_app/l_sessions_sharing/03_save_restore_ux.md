# Save/restore UX (manual save + restore)

**Audience:** everyone (wet lab + computational)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- Where the session controls live in the UI
- The exact steps to Save State and Load State
- How to confirm the one complete restore outcome
- Best practices for naming, organizing, and sharing session files

**Prerequisites:**
- A dataset loaded in the web app

---

## Where the session controls live

In the left sidebar:

1) Open **Session** (the accordion that also contains dataset loading controls).
2) Find the block labeled **Session state**.
3) Use:
   - **Save State** (downloads a `.cellucid-session` file)
   - **Load State** (opens a file picker and restores a `.cellucid-session` file)

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The Session accordion showing a dataset summary table, the sample dataset picker, local, remote-server and GitHub loading controls, and a Session state block containing Save State and Load State buttons.
:width: 524px

The **Session** accordion. Everything that brings a dataset in sits above, and
**Session state** sits at the bottom — the order matches the rule that the
dataset must be loaded before a session is restored onto it.
```

---

## Fast path (wet lab / non-technical): “save my work”

### Save State

1) Click **Save State**.
2) Wait for a “Session saved successfully” notification.
3) Find the downloaded file in your browser’s downloads.
   - It will end in `.cellucid-session`.
   - The default name includes a timestamp (so files don’t overwrite each other).

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-02-save-state-action.png
:alt: The Session state block with a mouse pointer over the Save State button, beside the Load State button.
:width: 488px

**Save State** and **Load State** sit side by side. Save downloads a file
immediately; nothing is uploaded anywhere.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-03-save-confirmation.png
:alt: A notification reading "Session saved successfully".
:width: 610px

The confirmation to look for. If you do not see it, the save did not happen —
check the browser's download permissions for the site.
```

The default filename looks like
`cellucid-session-2026-08-01T15-18-20.cellucid-session`: a fixed prefix and the
save time, so repeated saves never overwrite each other.

### Load State

1) Load the same dataset first.
2) Click **Load State**.
3) Choose the `.cellucid-session` file.
4) Wait for the download progress to finish with **Session fully restored** and
   for the panel status to say **Session loaded successfully**.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-04-restore-confirmation.png
:alt: A stack of notifications showing a download entry that ends in "Session fully restored.", a processing entry for LOD node mappings, a "Session fully restored." entry, and a final "Session loaded successfully" entry.
:width: 592px

A successful restore reports twice, and both matter. **Session fully restored.**
is the serializer confirming that every chunk applied and every owner committed.
**Session loaded successfully** is the panel confirming the operation returned.
Anything less than both is not a completed restore.
```

If the view already matched the file, the visual change may be subtle; use the
checks below.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-07-restored-view.png
:alt: The full application window immediately after a restore, with the dataset rendered in the viewer and the sidebar showing the restored control values.
:width: 1440px

The application immediately after a restore. The dataset was already open; the
session put the camera, the fields, the filters, and the sidebar controls back
where they were when it was saved.
```

:::{note}
This file-picker workflow is for a state you or a collaborator saved. Official
built-in samples are different: their catalog-advertised, integrity-verified
static starting state applies automatically after the scientific dataset is
published. See {doc}`04_official_sample_states`.

The catalog's internal `default.cellucid-session` is not an ordinary manual
bundle. If you download it and choose **Load State**, the strict manual loader
rejects its intentionally smaller five-chunk profile. Choose the sample and
then use **Save State** to create a portable manual bundle.
:::

---

## Practical path (computational): what happens under the hood

### Save State is a file download

When you click **Save State**:
- Cellucid captures state from multiple subsystems (camera, filters, highlights, etc.).
- It writes a single-file bundle and triggers a browser download.

There is no server upload by default.

### Load State is ordered and atomic

When you click **Load State** and pick a file:

1) Cellucid validates the whole container, current chunk inventory, exact
   dataset fingerprint, and declared byte bounds.
2) **Eager** chunks run before dependent **lazy** chunks.
3) The loader waits for large highlight memberships, analysis artifacts, and
   inactive user-defined codes too.
4) All feature owners prepare and commit, then the final UI refreshes.
5) Only then does the operation report success.

The scheduler yields between heavy chunks to keep Cancel and newer loads
responsive. Eager is not a partial-success boundary: a failure in any later
chunk rolls back earlier state.

---

## How to confirm a restore worked (quick checklist)

After loading:

1) Look for a session notification:
   - require **Session fully restored**, followed by the panel's **Session
     loaded successfully** status.
2) If the operation was rejected, read the exact error.
   - A dataset mismatch rejects everything; see
     {doc}`07_versioning_compatibility_and_dataset_identity`.
3) Confirm the high-signal state after success:
   - active field (Coloring & Filtering section)
   - Active filters panel (should match what you remember)
   - snapshot layout (grid vs single)
   - Camera Path and playback state
   - highlight memberships and analysis windows

If you press **Cancel** or start a newer restore, the older progress is
dismissed. That intentional cancellation is neither success nor a product
failure.

If anything is off, jump to {doc}`10_troubleshooting_sessions`.

---

## Worked example: renderer settings survive the round trip

The sidebar controls are restored generically, by writing each control's saved
value and then telling whatever owns that control that it changed. The four
renderer settings in **Visualization** are a good example to check against,
because they are visible, they are easy to set to something unmistakable, and
they drive the GPU rather than the UI.

Start from the defaults:

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-05-renderer-controls-default.png
:alt: The renderer settings block at its defaults, with Shader quality set to Full (lighting plus fog), Level-of-Detail unchecked, and Frustum culling unchecked.
:width: 480px

Defaults: **Shader quality** `Full (lighting + fog)`, **Level-of-Detail (LOD)**
off, **Frustum culling** off. The **Force LOD level** slider is hidden while LOD
is off, because it has nothing to act on.
```

Change all four, then **Save State**:

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-06-renderer-controls-changed.png
:alt: The renderer settings block with Shader quality set to Light (circular, no lighting), Level-of-Detail checked, Force LOD level set to 4, and Frustum culling checked.
:width: 480px

The saved state: shader quality `Light (circular, no lighting)`, LOD on, forced
to level `4`, frustum culling on. Enabling LOD is what reveals the **Force LOD
level** slider.
```

Reload the app with the same dataset — which puts all four back to the defaults
above — and **Load State** the file:

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-08-renderer-controls-restored.png
:alt: The renderer settings block after restoring the session, showing Shader quality Light, Level-of-Detail checked, Force LOD level 4, and Frustum culling checked, matching the saved state exactly.
:width: 480px

After the restore. All four match what was saved, and the viewer is drawing with
them — the panel and the renderer agree.
```

:::{important}
That last sentence is the point of the check, and it is worth doing once
yourself. A control counts as restored only when the thing it drives has
followed it: a panel showing `Light` in the dropdown while the viewer keeps
drawing at `Full` is a failed restore, not a cosmetic mismatch. Each of these
four controls belongs to the module that owns it, and that module is built
early enough in startup to hear the restore and act on it.

If you ever see a restored control whose effect has not followed it, that is a
defect worth reporting, not something to work around — the restore is supposed
to be all-or-nothing, and a control that refuses its restored value is required
to fail the whole restore rather than diverge quietly.
:::

---

## Best practices (these save you hours later)

### 1) Treat the dataset + session as a pair

Because sessions do not contain the dataset:
- store the `.cellucid-session` file next to the dataset export folder,
- or store a clear pointer to where the dataset is hosted (URL/GitHub path).

### 2) Use a naming convention that encodes intent

Recommended naming format:

`<dataset>__<goal>__<who>__<YYYY-MM-DD>__v<N>.cellucid-session`

Examples:
- `pbmc_demo__qc-filtered__kelly__2025-12-30__v1.cellucid-session`
- `suo__figure1_candidate__team__2025-12-30__v3.cellucid-session`

Why:
- session files are binary; the filename is your metadata unless you add a separate README.

### 3) Save “milestones”, not every minute

Sessions can become large when they include:
- many highlight groups with huge memberships,
- many user-defined categorical fields,
- analysis caches.

Save milestones like:
- “after QC filters + first labeling”
- “after final highlight pages”
- “before figure export settings”

### 4) Save before you load someone else’s session

Loading is destructive to the current UI state.

Workflow:
1) Save your current session.
2) Load the other session.
3) If needed, switch back by loading your saved session.

---

## Interface reference

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: Cellucid Session panel showing sample, local-file, remote-server, GitHub, and session-state controls.
:width: 524px

The Session panel presents each loading path separately and keeps Save State and Load State beside the dataset controls.
```

---

## Next steps

- {doc}`04_official_sample_states` (understand automatic, integrity-verified static starting states for official samples)
- {doc}`05_share_workflows_links_bundles_exports` (how to send sessions to humans safely)
