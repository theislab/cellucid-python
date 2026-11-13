# What gets saved and restored

**Audience:** computational users, power users, and anyone who has been surprised by "restore didn't restore X"  
**Time:** 20–40 minutes  
**What you'll learn:**
- The exact "saved vs not saved" list for `.cellucid-session` files
- The sixteen parts of the state each step records
- Which parts depend on the dataset and why an identity mismatch refuses everything
- How to check a load quickly

**Prerequisites:**
- A dataset loaded in the web app
- Optional: a saved `.cellucid-session` file you can test with

---

## Quick summary (read this if you only have 60 seconds)

A session file holds the **history** of the dataset shown, as far as **Share…**
includes it ({doc}`../m_state_history/06_sharing_and_replay`):
- the trackers and branches ticked (by default every tracker with its main
  branch and the branch it shows), with all of their steps, including steps
  you undid, and their groups;
- the step each tracker shows, and which tracker is in use;
- for each step: what was done (the action and its settings, by name), who did
  it (you, an agent by its name and model, or a notebook), when (exact, dates
  only or no times, as chosen), and the value of every part it changed, before
  and after;
- notes and bookmarks (unless left out), conversations with agents and
  Recently deleted (only when ticked), analysis results (or only their
  settings, to compute on open: once the file opens, each analysis shown at
  its step is run again as a step by *Session*, and you are told why when its
  inputs differ from the shared result's; one not shown is left with its
  settings and tried again at the next opening).

**Current state** writes one step holding what is shown, with no notes,
conversations, agents or times; **Steps only** writes the steps without cells
or values, a plan to replay on any dataset.

At every step, these sixteen parts are kept:
- **cells to open**: every cell or a subsample, and how many cells genes are
  read for ({doc}`../b_data_loading/21_subsampled_loading`);
- **data mapping**: which stored embeddings and arrays are plotted and read
  ({doc}`../b_data_loading/13_data_mapping`);
- **expression values**: whether the genes are shown and analyzed normalized
  to 10,000 counts per cell and/or `log1p` (Data mapping › Expression values);
- **fields and labels**: renamed fields and categories, deleted fields, and
  fields you made or a notebook pushed, with their values;
- **gene sets**: made in the app or pushed by a notebook;
- **views**: the kept views, their layout, the focused view, linked cameras,
  the embedding each view shows, centroid points and labels;
- **active field**: what each view is colored by, a field or a gene;
- **filters**, per view: hidden categories, continuous ranges, outlier and
  missing-value filters, and **filter by cell count**;
- **colors**, per view: category colors, color maps, color ranges, log
  scale, Gray for zero;
- **hidden cells**, per view: cells hidden one by one;
- **highlights**: highlight pages and their groups of cells;
- **cameras**, per view, with the navigation mode;
- **camera path**: keyframes and playback settings, an empty path included;
- **settings**: every stable sidebar control (background, point and lighting
  settings, overlays, render mode, …) and every Figure Export setting;
- **analysis**: the analysis panel and floating analysis windows, with their
  settings and results;
- **layout**: docked and floating panels and where they float, open sidebar
  sections, the sidebar, and the theme.

Also in the file: where the dataset opens from (unless you chose **Name
only**), the cells it was saved on, the dataset fingerprint, every observation
field and gene the shown step uses, and the selection tool chosen (put back
after a load; choosing a tool is never a step).

Not saved:
- the dataset itself (cells, fields, gene expression);
- dataset selection and connection controls (the Open… menu, the **Link or
  server** bar, the **Sample dataset** list); the session records where its
  dataset came from instead (see {doc}`03_save_restore_ux`);
- per-cell gene values analyses read (read again from the dataset);
- Community Annotation (votes, comments, moderation, GitHub sync and sign-in,
  drafts);
- **Antialiasing** (a property of the machine drawing the view);
- benchmarking state;
- notifications, and anything in progress: a selection not yet confirmed, a
  pointer lock, hover;
- which floating panel is in front.

If the dataset fingerprint differs, the whole load is refused and nothing
changes; see {doc}`07_versioning_compatibility_and_dataset_identity`.

---

## The contract: "a session file keeps the history of the state, not the data"

The mental model:

- A session file restores **what you did and how the app is set up**.
- The dataset loader restores **what data the app is set up on**.

This is intentional because:
- datasets can be huge (sessions should be shareable and fast),
- dataset access methods differ (local folder vs server vs GitHub),
- sessions are treated as untrusted input and should not be a "data injection" path.

---

## Saved and restored, part by part

This section is the checklist of what you should expect. Each part is one
value per step; a load puts back the value of each part at the step shown,
and **Undo** puts back the values of the step before.

:::{note}
A session holds every part's value as it is: there is no separate format per
feature, and nothing is left out because a feature "was not ready". A file
holding a part this version of Cellucid does not have is refused, naming it.
:::

### Cells to open

- every cell, or the subsample's specification (method, share, seed and so
  on) and its exact cells;
- how many of the open cells genes were read for.

Changing **Cells to open** (Subsampling › **Apply**) is a step like any other:
**Undo** reads the dataset again with the cells before. Loading a session opens the
dataset with the cells the session was saved on first, reading it again if it
is open with other cells
({doc}`../b_data_loading/21_subsampled_loading`).

### Data mapping

The embeddings that are turned on and the default of each dimension, and every
other choice under **Session › Data mapping** (modality, neighbor graph, gene
names, expression matrix and its transform, highly variable genes, latent
space, vector fields). It is put back right after the cells and before the
views, so every view the session names can be shown; see
{doc}`../b_data_loading/13_data_mapping`.

### Expression values

Whether the genes are normalized to 10,000 counts per cell and log-transformed
(**Session › Data mapping › Expression values**), as asked. Put back right after
the data mapping and before the views, so a view colored by a gene, its color
range and its filter are restored in the values they were set in; restoring a
normalized state reads each cell's total first (an obs column, else every gene
once). A session saved before this part existed opens with the values as
stored.

### Views and cameras

- the kept views and their layout (single or grid), the focused view, whether
  the live view is hidden;
- whether **Link cameras** is on (views of one embedding share a camera);
- the embedding each view shows, as an embedding reference (an `obsm` key with
  its columns); a view's dimension is that embedding's;
- whether each view shows centroid points and centroid labels;
- every view's camera position and navigation mode.

Each view is saved from its own state, whichever view is selected when you
save: saving with a kept view selected still records the live view's own
coloring field and filters as the live view's.

### Active field, filters, colors and hidden cells (per view)

- what each view is colored by: an observation field or a gene;
- filters: hidden categories, continuous ranges, the outlier filter of a
  categorical field and its threshold, missing values hidden, a field's filter
  switched off, and **filter by cell count** (a category keeping only some of
  its cells shown);
- colors: category colors, color maps, color ranges, log scale, Gray for
  zero;
- cells hidden one by one.

Only what differs from how the field opened is recorded, so a session stays
small however many fields you looked at. A color range is measured against
the data, not against the filter.

### Fields and labels

- renamed fields and renamed categories;
- deleted fields;
- fields you made (or a notebook pushed): their definition, and their values
  per cell by original position.

A categorical field's category centroids are not stored: they are derived from
its values and the coordinates of each embedding, and computed again when the
field comes back, so its labels sit exactly where they were.

### Gene sets and highlights

- gene sets, with their names, descriptions and genes;
- highlight pages (name, color) and their groups (label, on or off, and the
  exact cells, by original position).

### Camera path

Keyframes with their holds, timing, interpolation, **Loop playback** and **Autoplay**, and an
explicitly empty path. A camera path saved with **Autoplay** on starts playing
once the load is done; otherwise it stays stopped.

### Settings

Every stable control of the sidebar and its floating panels is captured by its
DOM id, including the background, point and lighting settings, the render mode
and overlays. A restore applies each exactly as a choice made by hand.

- The connectivity **Max edges** slider records the number of edges you asked
  for (250,000 until you move it), not the count it shows: its range is the
  graph's visible edges, read only once edges are shown, so a restore keeps
  your number and the slider shows as much of it as the graph allows.
- Every **Figure Export** setting is saved: plot size and preset, the frame,
  background, font and text sizes, title and axis labels, axes, legend and its
  position, orientation, depth sorting, reference grid, selection emphasis,
  colorblind preview, download format, DPI, SVG strategy, and **Export all
  views**.
- Saved settings are read one at a time. A value this version of Cellucid
  cannot use — a control removed since the file was saved, a choice the
  control no longer offers, a Figure Export value no version has — keeps
  the setting as it is shown, and one notice says which (`2 saved settings
  could not be used and stay as shown`, each with its reason under
  **Details**); every other setting is restored. A setting added since the
  file was saved keeps its value. The Figure Export settings are one record
  with a version: a record of any version but today's is refused whole (one
  notice line, the reason under **Details**), the panel keeps its settings, and
  everything else in the session is restored.

:::{note}
"Captured by id" means the control's DOM identifier is the key, not the label
you read beside it. Renaming a control's visible text does not change what a
session restores.

Restoring a control writes its value and then fires the event its owner
listens for, and the owner's own handler is what passes the value on. This is
why a restored control and the thing it drives must never disagree: if the
owner refuses a value it was handed, the whole load fails instead of leaving
the panel showing something the viewer never accepted. (A value the control
does not offer is never handed to it: it is one of the settings the notice
names, above.)
:::

:::{note}
**`Antialiasing` is the one exception.** It is not recorded in a session, and
loading one **never changes it**, for two reasons.

It describes the machine you are sitting at rather than the analysis you are
doing, and your choice is stored in this browser rather than in the file.

And it has three states while a checkbox has two. The stored preference is
`auto` (the default, and what no stored preference means), or `on`, or `off`.
Under `auto`, Cellucid decides from the size of the dataset you have open and
the dimension of the view you are working in: on below 5,000,000 cells in 3D and
below 200,000 in a 1D or 2D view, off at or above. It decides again on every
dataset publication and when the view switches between 2D and 3D. Clicking
the box ends automatic selection for good, in both directions.

So a session you open, including the starting view a built-in sample
publishes, keeps your choice, or the dataset's choice when you have not made
one. The first line of the control's (i) tells you which you are on: it reads
`Chosen automatically: <n> cells in <d>D.` while the dataset is deciding, and
is gone once you have clicked the box. (It reads `This browser is not
providing antialiasing for this view.` when the device cannot multisample.)
:::

The four renderer settings are the easiest instance of a restored setting to
check yourself, because you can see both the control and its effect:

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-06-renderer-controls-changed.png
:alt: The same controls changed: Shader quality Light, Level of detail (LOD) checked with Force LOD level at 4, and Frustum culling checked.
:width: 496px

Four sidebar controls set to non-default values before the session was shared.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-08-renderer-controls-restored.png
:alt: The same controls after a saved session was loaded into a fresh page: Shader quality Light, Level of detail (LOD) checked with Force LOD level 4, and Frustum culling checked, as they were saved.
:width: 496px

The same four after the session was loaded into a freshly started app that began at the
defaults. Walked through step by step in {doc}`03_save_restore_ux`.
```

### Analysis

- **the analysis panel in the sidebar**: which analysis is open, and every
  analysis's settings and result;
- **floating analysis windows**: which analysis each shows, its settings, its
  result, and its place and size.

Results that are expensive to compute (differential expression, marker genes,
gene signatures, correlation, dot plot, composition) are stored with the
session, genes and categories by name, so they come back without being run
again. A result is saved with what it was computed on, so it comes back
current, or marked **Stale · Run** when that has changed on the page it is
loaded into. A result naming a gene or category the open dataset lacks is
refused. A result's plot is drawn when its panel is shown; the load never
waits for that.

The per-cell gene values analyses read are not stored: they are read again
from the dataset when needed.

### Layout

Docked and floating panels and where they float, which sidebar sections are
open, the sidebar (open or closed, and a width you dragged it to), and the
theme. Which floating panel is in front is not saved.

---

## Not saved, not restored (intentional exclusions)

This is the "don't expect it to come back" list.

### Data and data loading

- the dataset itself (cells, observation and gene tables, expression matrices)
- dataset selection and connection controls
  - the Open… menu (a file or folder from this computer)
  - the **Link or server** bar
  - the **Sample dataset** and **My datasets** lists

### Kept elsewhere, or not at all

- **Antialiasing**: a property of the device (see the note above)
- **Benchmarking** state
- **Community Annotation** state
  - votes, comments, moderation, drafts
  - GitHub sync and sign-in

### In-progress and on-screen only

- a selection you haven't confirmed yet
- a pointer lock, hover state
- notifications
- loading spinners

### Rendering buffers

Sessions do **not** store GPU buffers, such as:
- per-cell color buffers
- transparency buffers
- other derived rendering caches

Those are computed again from the parts above after a load.

---

## Dataset identity behavior (critical)

On load, Cellucid compares a **dataset fingerprint** saved in the session with
the dataset that is open.

If they differ, Cellucid refuses the whole load before anything changes: no
layout, filters, fields, highlights or analyses are put back on their own. A
dataset that lacks a field or gene the shown step uses is refused the same
way, with every missing key named.

This is a safety feature: it keeps filters and highlights from landing on the
wrong cells.

See {doc}`07_versioning_compatibility_and_dataset_identity` for:
- what counts as "the same dataset",
- common ways you accidentally create a mismatch,
- collaboration-safe workflows.

---

## How to verify a restore (practical checklist)

After you load a session, check in this order:

1) **Did you load the correct dataset?**
   - A mismatch refuses the load. Open the dataset the message names, then load
     the unchanged file again.
2) **Did the load finish?**
   - The load must end in a notice reading **Session fully restored.**
   - **Session restored but for N parts** names each part that could not be
     put back (left as it was shown), one sentence each, with the failure's
     words under **Details**; a part saved by another version of Cellucid
     that this one cannot read says so.
   - A progress card still running is not success.
3) **Is the history there?**
   - **State › History** lists the steps, with the step shown marked;
     **Tracker** lists every tracker the file holds.
4) **Does the state look right?**
   - camera and framing
   - active field
   - filters (the legend and Active filters)
   - kept views and their layout
5) **Spot-check exact replacement**
   - an empty saved camera path is empty after the load
   - a saved analysis result is shown without running it again
6) **Spot-check derived consequences**
   - highlight counts
   - analysis panel and windows, their settings and results

If any part is wrong, jump to {doc}`10_troubleshooting_sessions`.

---

## Next steps

- {doc}`03_save_restore_ux` (how to use **Share…** / **Load…** safely)
- {doc}`12_reference` (the file format)
- {doc}`../m_state_history/01_the_history` (steps, undo and redo)
