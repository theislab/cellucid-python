# Save/restore UX (manual save + restore)

**Audience:** everyone (wet lab + computational)  
**Time:** 10–20 minutes  
**What you'll learn:**
- Where the session controls live in the UI
- The exact steps to share (save) and load a session file
- How to pick up where you left off without a file
- How to confirm a load finished
- Best practices for naming, organizing, and sharing session files

**Prerequisites:**
- A dataset loaded in the web app

---

## Where the session controls live

In the left sidebar:

1) Open **State** (the section below **Session**).
2) Find the group labeled **Session**.
3) Use:
   - **Share…** (chooses what the file holds, then downloads a `.cellucid-session` file of the dataset shown)
   - **Load…** (every history this browser keeps, and **Open a session file…**; see
     {doc}`../m_state_history/03_saved_states_and_autosave`)

Below **Session**, the same section shows the **History**: the steps as a
graph, the branch shown and the tracker in use, what a session file holds;
see {doc}`../m_state_history/index`.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The State section with the Pancreatic endocrinogenesis sample open: Session with Share… and Load…, and the History with the recording switch, the branch shown and the graph of the steps taken so far, the step shown marked.
:width: 540px

The **State** section holds the session's **Share…** and **Load…** above the
history: **Share…** writes those steps to a file, **Load…** puts a file's
steps back beside yours.
```

The **Session** section above it is where data comes in: **Open data**,
**Subsampling** (**Cells to open**) and **Data mapping**. A session is loaded
onto a dataset, which **Load…** opens first when the session records where it
is.

---

## Fast path (wet lab / non-technical): “save my work”

### Save

1) Under **State › Session**, click **Share…**.
2) Leave **Share** on **History** (the trackers and branches ticked, with
   every value); **What's inside** lists exactly what the file will hold.
   The other choices are in {doc}`../m_state_history/06_sharing_and_replay`.
3) Click **Save file**, and wait for the notification *Saved <name>.cellucid-session.*
4) Find the downloaded file in your browser's downloads.
   - It ends in `.cellucid-session`.
   - It is named after the dataset (`Pancreatic endocrinogenesis.cellucid-session`);
     your browser adds a number when a file of that name is already there.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-02-save-state-action.png
:alt: The Session group of the State section with its (i), Share… and Load… side by side, and the pointer on Share….
:width: 504px

**Share…** and **Load…** sit side by side. **Share…** chooses what the file
holds and downloads it; nothing is uploaded anywhere.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-03-save-confirmation.png
:alt: A notification reading "Saved Pancreatic endocrinogenesis (scVelo).cellucid-session." after Share… › Save file wrote the session file.
:width: 760px

The confirmation to look for. If you do not see it, the file was not written:
check the browser's download permissions for the site.
```

**Share…** stays greyed out until a dataset is open. Changes made while
recording was paused become one step first. To keep a version in this
browser instead of a file, type a name in **Load…** and **Save**
({doc}`../m_state_history/03_saved_states_and_autosave`).

### Load

1) Under **State › Session**, click **Load…**.
2) Under **From a file**, click **Open a session file…** and choose the
   `.cellucid-session` file.
3) Wait for the progress notification to finish as **Loaded** with
   **Session fully restored.**

The file's trackers are added beside yours, and the file's active tracker is
shown at the step it showed when it was saved; **Undo** walks back through its
steps (see {ref}`loading-adds-beside` below). A camera path saved with
**Autoplay** on starts playing.

A session records where its dataset was opened from: a sample dataset, a data
link (the **Link or server** bar or `?data=`), or a Cellucid server. When another
dataset, or none, is open, **Load…** opens that one first ("Opening the
session's dataset from …") and then loads the session. It also opens the
cells the session was saved on, reading the dataset again when it is open with
other cells. A session saved on a file or folder from your own computer records
only its name, since a page cannot open a file on its own: in Chrome and Edge
the browser that opened it may have kept it, and **Load…** then reopens it, or
offers **Reopen <name>** ({ref}`reopen-local-file`). One saved in a notebook
records no source. Otherwise, open that dataset yourself before loading the
session (a refusal names the dataset it needs).

A session saved on a link records the link in its public form: a key,
signature or token the link carried (a Dropbox share key, a signed bucket
address) is never written into the session. Whoever has access to the data
opens it from the session; anyone else opens the link you sent them first,
then loads the session ({ref}`link-access-keys`).

If a dataset is still loading when you choose the file, a notification says
the session is waiting for it, and the session is loaded as soon as that
dataset has loaded. Choosing another dataset in the meantime means the session
is loaded onto that one instead, and refused with an explanation if it was
saved on a different dataset. See {ref}`changing-dataset-while-loading`.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-04-restore-confirmation.png
:alt: The two notifications a completed restore publishes: "Loaded "Tracker 1 (shared)" from Pancreatic endocrinogenesis (scVelo).cellucid-session beside yours." with Undo, and the restore's progress entry, now titled "Loaded" and reading "Session fully restored.".
:width: 760px

A successful load says what it added (**Loaded "…" from <file> beside
yours.**, with **Undo**, which takes the loaded trackers out again), and its
own progress entry completes as **Session fully restored.**: every step and
every value read and checked, and the step shown put back. A load that fails ends the same entry as **Loading
failed** instead, naming the reason.
```

If the view already matched the file, the visual change may be subtle; use the
checks below.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-07-restored-view.png
:alt: The whole application window immediately after a session file was loaded, with the dataset drawn in the viewer and the State section showing the loaded history, the step shown marked.
:width: 1440px

The application immediately after a load. The dataset was already open; the
session put its history back, and with it the camera, the fields, the filters
and the sidebar controls of the step that was shown when it was saved.
```

:::{note}
This file-picker workflow is for a history you or a collaborator saved. A
dataset opened from a catalog is different: the starting view its catalog
names (the built-in samples have one) is integrity-verified and applies
automatically after the dataset is published, before its history starts. See
{doc}`04_official_sample_states`.

To share your work on a sample, choose the sample, work from its starting
view, and use **Share…**.
:::

---

## Continue where you left off (no file needed)

You do not need a file to pick up your own work on the same computer:

- Every step is kept in this browser as you make it (**Autosave in this
  browser**, on by default, in **History settings**).
- When you open a dataset again whose history this browser kept, a
  notification asks **Continue where you left off?** with **Continue** and
  **Start fresh**. Nothing is loaded until you choose **Continue**; **Start
  fresh** stops the question for that history.
- From **Session › Open data › My datasets › All…**, a dataset's **Continue**
  opens it with its latest saved state.
- **State › Session › Load…** lists every history kept, for every dataset:
  open one, replay its steps here, rename it, or delete it.

These stay in this browser only. A session file is how a history goes to
another browser, computer or person. Details:
{doc}`../m_state_history/03_saved_states_and_autosave`.

---

## Practical path (computational): what happens under the hood

### Share is a file download

When you click **Save file** in **Share…**:
- Cellucid waits for actions under way to finish, and turns changes made
  while recording was paused into one step, so what is saved is what you see.
- It takes what you chose (the trackers and branches ticked with every value
  their steps name, the current state only, or the steps only; notes,
  conversations and Recently deleted when ticked; results, times and how the
  dataset is named as chosen) together with the dataset's fingerprint, where
  the dataset opens from (unless you chose its name only), and the cells that
  are open.
- It writes a single file and starts a browser download.

There is no server upload.

### Load is ordered and all or nothing

When you choose **Open a session file…** and pick a file:

1) Cellucid opens the session's dataset if needed, with the cells it was saved
   on.
2) It checks the file's signature and manifest, the dataset fingerprint, and
   the fields and genes the shown step uses, before anything changes.
3) It reads the steps and the values of the step shown, adds the file's
   trackers beside yours and shows that step.
4) It reads the values of every other step while that step is shown; a new
   step waits until they are read.
5) Only then does it report **Session fully restored.**

A failure at any point, or **Cancel**, puts back what was shown before and
keeps the trackers you had.

---

## How to confirm a restore worked (quick checklist)

After loading:

1) Look for the session notification:
   - **Loaded** with **Session fully restored.** on the load's progress entry.
2) If the load was refused, read the exact message.
   - A dataset mismatch refuses everything; see
     {doc}`07_versioning_compatibility_and_dataset_identity`.
3) Confirm the high-signal state after success:
   - **State › History**: the steps, with the step shown marked, and the line
     above the graph saying where the tracker comes from
   - active field (Coloring & filtering section)
   - filters (legend and Active filters)
   - kept views and their layout (grid vs single)
   - camera path and playback state
   - highlight groups and analysis windows

If you press **Cancel** or start a newer load, the older progress is
dismissed. That cancellation is neither success nor a failure.

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

Change all four, then **Share…** and **Save file**:

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-06-renderer-controls-changed.png
:alt: The same controls changed: Shader quality Light, Level of detail (LOD) checked with Force LOD level at 4, and Frustum culling checked.
:width: 496px

The saved state: shader quality `Light (circular, no lighting)`, LOD on, forced
to level `4`, frustum culling on. Enabling LOD is what reveals the **Force LOD
level** slider.
```

Open the same dataset in a fresh page (do not choose **Continue** if it offers
to continue), which shows all four at the defaults above, and load the file
(**Load… › Open a session file…**):

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-08-renderer-controls-restored.png
:alt: The same controls after a saved session was loaded into a fresh page: Shader quality Light, Level of detail (LOD) checked with Force LOD level 4, and Frustum culling checked, as they were saved.
:width: 496px

After the load. All four match what was saved, and the viewer is drawing with
them: the panel and the renderer agree.
```

:::{important}
That last sentence is the point of the check, and it is worth doing once
yourself. A control counts as restored only when the thing it drives has
followed it: a panel showing `Light` in the dropdown while the viewer keeps
drawing at `Full` is a failed restore, not a cosmetic mismatch.

If you ever see a restored control whose effect has not followed it, that is a
defect worth reporting, not something to work around: a control that refuses
its restored value is required to fail the whole load rather than diverge
quietly.
:::

**Undo** shows the same thing without a file: each of the four changes is a
step, and undoing them puts the defaults back, in the panel and in the
renderer.

---

## Best practices (these save you hours later)

### 1) Treat the dataset + session as a pair

Sessions do not contain the dataset. A session saved on a sample, a data link
or a server points to it and opens it; one saved in a notebook does not, and
one saved on a local file names the file. In Chrome and Edge, **Load…** in the
browser that opened that file reopens it from the copy kept of it
({ref}`reopen-local-file`); anywhere else the file has to be opened first. So:
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
- session files are binary; the filename is your metadata unless you add a
  separate README. Inside Cellucid, rename trackers and branches to say what
  each line of work is ({doc}`../m_state_history/02_branches_trackers_and_comparing`).

### 3) Share a file when the work leaves this browser

This browser already keeps every step. Share a file when you want the work
somewhere else:
- to send it to someone,
- to keep it with a paper, a figure or a lab notebook entry,
- before clearing this browser's data or moving to another computer.

One file holds every step, so you do not need one file per step. Files grow
with:
- many highlight groups with huge numbers of cells,
- many fields you made,
- many analysis results.

(loading-adds-beside)=
### 4) Loading adds the file's trackers beside yours

**Open a session file…** adds the file's trackers beside the ones on screen and shows the
file's active one; nothing of yours is replaced. A loaded tracker whose name
you already use is named apart (`"Tracker 1 (shared)"`). A notice says what was
loaded, with **Undo**: it moves the loaded trackers to **Recently deleted**,
from where they come back like any removed tracker. Each load is kept in
**Recently deleted** (**⋯ › Recently deleted…**) as "Loaded … from <file>", so it can be undone
later too. A file of the steps only opens as a replay instead: a tracker of its
own, its steps applied again on this dataset. A file of another dataset can be
replayed on the dataset shown with **Load… › Replay a session file here…**
({doc}`../m_state_history/06_sharing_and_replay`).

Opening one of the saved states in **Load…** replaces the history shown; the
history you had stays in its own saved state, listed there.

---

## Next steps

- {doc}`04_official_sample_states` (starting views a catalog names, verified and put back by themselves)
- {doc}`05_share_workflows_links_sessions_data_files` (how to send sessions to humans safely)
- {doc}`../m_state_history/03_saved_states_and_autosave` (the history kept in this browser)
