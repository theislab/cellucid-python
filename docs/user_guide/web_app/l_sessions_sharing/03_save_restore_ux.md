# Save/restore UX (manual save + restore)

**Audience:** everyone (wet lab + computational)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- Where the session controls live in the UI
- The exact steps to save and load a session file
- How to confirm the one complete restore outcome
- Best practices for naming, organizing, and sharing session files

**Prerequisites:**
- A dataset loaded in the web app

---

## Where the session controls live

In the left sidebar:

1) Open **Session** (the accordion that also contains dataset loading controls).
2) Find the group labeled **Session file**.
3) Use:
   - **Save** (downloads a `.cellucid-session` file)
   - **Load** (opens a file picker and restores a `.cellucid-session` file)

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The Session accordion showing the Dataset box for the loaded sample, the Open data group with the Sample dataset dropdown, the Open… button and the Link or server bar, the Session file group with Save and Load below them, and the collapsed Data mapping group.
:width: 540px

The **Session** accordion. Everything that brings a dataset in sits under
**Open data**, and **Session file** sits below it — a session is restored onto
a dataset, which **Load** opens first when the session records where it is.
```

---

## Fast path (wet lab / non-technical): “save my work”

### Save

1) Click **Save**.
2) Wait for a “Session saved successfully” notification.
3) Find the downloaded file in your browser’s downloads.
   - It will end in `.cellucid-session`.
   - The default name includes a timestamp (so files don’t overwrite each other).

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-02-save-state-action.png
:alt: The Session file group with its (i), Save and Load side by side, and the pointer on Save.
:width: 504px

**Save** and **Load** sit side by side. **Save** downloads a file
immediately; nothing is uploaded anywhere.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-03-save-confirmation.png
:alt: A notification reading "Session saved successfully" after Save produced a session bundle.
:width: 520px

The confirmation to look for. If you do not see it, the save did not happen —
check the browser's download permissions for the site.
```

The default filename looks like
`cellucid-session-2026-08-01T15-18-20.cellucid-session`: a fixed prefix and the
save time, so repeated saves never overwrite each other.

### Load

1) Click **Load**.
2) Choose the `.cellucid-session` file.
3) Wait for the restore's progress notification to finish as **Loaded** with
   **Session fully restored.**

A session records where its dataset was opened from: a sample dataset, a data
link (the **Link or server** bar or `?data=`), or a Cellucid server. When another
dataset — or none — is open, **Load** opens that one first ("Opening the
session's dataset from …") and then applies the session. A session saved on a
file or folder from your own computer records only its name, since a page
cannot open a file on its own, and one saved in a notebook records no source:
open that dataset yourself before loading the session (a refusal names the
file it needs).

A session saved on a link records the link. A link that carries its own access
key — a Dropbox share link, a signed bucket address — is then in the session
file too, so share such a session only with people who may open the data.

If a dataset is still loading when you choose the file, a notification says
the session is waiting for it, and the session is applied as soon as that
dataset has loaded. Choosing another dataset in the meantime means the session
is applied to that one instead, and refused with an explanation if it was saved
from a different dataset. See
{ref}`changing-dataset-while-loading`.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-04-restore-confirmation.png
:alt: The one notification a completed restore publishes: the restore's progress entry, now titled "Loaded" and reading "Session fully restored.".
:width: 520px

A successful restore reports once. **Session fully restored.** is the
restore's own progress entry completing: every chunk applied and every owner
committed. A restore that fails ends the same entry as **Loading failed**
instead, naming the error.
```

If the view already matched the file, the visual change may be subtle; use the
checks below.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-07-restored-view.png
:alt: The whole application window immediately after a session restore, with the dataset drawn in the viewer and the sidebar showing the restored control values.
:width: 1440px

The application immediately after a restore. The dataset was already open; the
session put the camera, the fields, the filters, and the sidebar controls back
where they were when it was saved.
```

:::{note}
This file-picker workflow is for a state you or a collaborator saved. A
dataset opened from a catalog is different: the default session its catalog
names (the built-in samples have one) is integrity-verified and applies
automatically after the scientific dataset is published. See
{doc}`04_official_sample_states`.

A catalog's default session is not an ordinary manual bundle. If you download it and choose **Load**, the strict manual loader
rejects its intentionally smaller seven-chunk profile. Choose the sample and
then use **Save** to create a portable manual bundle.
:::

---

## Practical path (computational): what happens under the hood

### Save is a file download

When you click **Save**:
- Cellucid captures state from multiple subsystems (camera, filters, highlights, etc.).
- It writes a single-file bundle and triggers a browser download.

There is no server upload by default.

### Load is ordered and atomic

When you click **Load** and pick a file:

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
   - require **Session fully restored.** on the restore's progress entry.
2) If the operation was rejected, read the exact error.
   - A dataset mismatch rejects everything; see
     {doc}`07_versioning_compatibility_and_dataset_identity`.
3) Confirm the high-signal state after success:
   - active field (Coloring & filtering section)
   - Filters panel (should match what you remember)
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
:alt: The Cells and Performance controls at their defaults: Shader quality Full (lighting + fog), Point size, the Depth cues sliders, and under Performance Level of detail (LOD) and Frustum culling unchecked and Antialiasing checked.
:width: 496px

Defaults: **Shader quality** `Full (lighting + fog)`, **Level of detail (LOD)**
off, **Frustum culling** off. The **Force LOD level** slider is hidden while LOD
is off, because it has nothing to act on.
```

Change all four, then click **Save**:

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-06-renderer-controls-changed.png
:alt: The same controls changed: Shader quality Light, Level of detail (LOD) checked with Force LOD level at 4, and Frustum culling checked.
:width: 496px

The saved state: shader quality `Light (circular, no lighting)`, LOD on, forced
to level `4`, frustum culling on. Enabling LOD is what reveals the **Force LOD
level** slider.
```

Reload the app with the same dataset — which puts all four back to the defaults
above — and **Load** the file:

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-08-renderer-controls-restored.png
:alt: The same controls after a saved session was loaded into a fresh page: Shader quality Light, Level of detail (LOD) checked with Force LOD level 4, and Frustum culling checked, as they were saved.
:width: 496px

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

Sessions do not contain the dataset. A session saved on a sample, a data link
or a server points to it and opens it; one saved in a notebook does not, and
one saved on a local file names the file. In Chrome and Edge, **Load** in the
browser that opened that file reopens it from the handle kept of it, with one
click on **Reopen <name>** ({ref}`reopen-local-file`); anywhere else the file
has to be opened first. So:
- store the `.cellucid-session` file next to the dataset file,
- or host the dataset (see {doc}`../b_data_loading/14_open_from_link`) and save
  the session on the hosted copy, so the session opens it by itself.

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
- many analysis results.

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
:alt: The Session panel with the Pancreatic endocrinogenesis sample loaded: the Dataset box with its name, cells and genes, format and host, description and a collapsed Details; Open data with the Sample dataset dropdown, the Open… button and the Link or server bar; Save and Load under Session file; and the collapsed Data mapping group.
:width: 532px

The Session panel gathers every loading path under **Open data** and keeps the **Save** and **Load** buttons under **Session file**, directly below it.
```

---

## Next steps

- {doc}`04_official_sample_states` (understand automatic, integrity-verified default sessions named by a catalog)
- {doc}`05_share_workflows_links_sessions_data_files` (how to send sessions to humans safely)
