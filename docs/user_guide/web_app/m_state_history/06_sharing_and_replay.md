# Sharing, loading and replaying

**Audience:** everyone
**Time:** 10 minutes
**What you'll learn:**
- What **Share…** writes, and how to choose exactly what a file holds
- How to save the dataset itself, with what you did written in, as a new AnnData file
- What happens when you load a history someone shared
- How to replay a history on another dataset ("run my workflow on my next sample")
- How to apply one branch's steps again, and how to verify that a branch reproduces

---

## Share a history

**State › Session › Share…** opens **Share this history**. **Save file** writes
a `.cellucid-session` file of the dataset shown with what you chose; nothing
is uploaded, and the dataset itself is never in it.

**Share** decides what kind of file it is:

| Choice | What the file holds | Opens as |
|---|---|---|
| **History** (the default) | The trackers and branches ticked, with every value their steps name | The same steps, on the same dataset |
| **Current state** | Only what is shown now, as one step ("Shared state"): no notes, conversations, agents or times | The state, without how you got there |
| **Steps only** | The steps without cells or values: a replay plan | A replay on any dataset ({ref}`history-replay-another-dataset`) |

**Include** (for **History** and **Steps only**):

- the trackers and branches: every tracker with its main branch and the
  branch it shows is ticked; tick or untick others (archived branches are
  not ticked);
- **Notes and bookmarks** (ticked);
- **Conversations** with agents, with how many there are (not ticked);
- **Recently deleted** (not ticked).

Then **Results** (**Include** the analysis results, or **Compute on open** from
their settings; always computed again for **Steps only**), **Times**
(**Exact**, **Dates**, or **None**; **None** by default for **Current state**)
and **Dataset** (**Name and where** it opens, the file name, link or server, or
**Name only**).

**Results left out** (**Compute on open**) are computed again once the file is
opened, for each analysis shown at the step it opens on: each run is an
ordinary step by *Session*, which Undo takes back. A notification says what
was computed again (*Computed again on open (left out of the file): …*), and
says why when its inputs differ from the shared result's (the cells open, the
filters, the genes). An analysis that is not shown there cannot draw, so it
is not run: it says so, keeps its settings without a result, and is tried
again the next time the file is opened. Run it yourself to see its result.

The line above the buttons counts what the file will hold (`14 steps · 2
branches · 3 results`). **What's inside** lists each kind of content with its
count: branches, cell sets and values, analysis results, notes and bookmarks,
conversations, agents and models, the names of people, Recently deleted,
times, and what is said of the dataset. Read it before you send a file.

Changes made while recording is paused become one step first
({ref}`history-recording-paused`). A conversation, a note or your name
({ref}`history-settings`) is in the file only when you include it.

---

(history-save-anndata)=
## Save AnnData with this state

**⋯ › Save AnnData with this state…**, in the **History** header, writes a
**new** file of the dataset shown: a Cellucid-ready AnnData `.zarr.zip`
({doc}`../../python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification`)
with the state shown written in. Python, R and any AnnData tool read it, and
Cellucid opens it like any other file. The dataset's own file is never
written over.

What is shown is what is saved. Changes made while recording is paused become
one step first, as for **Share…**. The state is then applied to the dataset
by the same rules as `apply_cellucid_session_to_anndata()` in Python and
`cellucid_apply_session()` in R
({doc}`../../python_package/e_jupyter_hooks/11_session_bundles_apply_session_to_anndata`):

| In the state | In the new file |
|---|---|
| A renamed field or category | The obs column or category renamed in place |
| A field you made (labels, a copy, a pushed column) | One obs column, named by the field's key |
| A highlight group, on any page | One true/false obs column, `cellucid_highlight__<page>__<group>` |
| A gene set | Kept in `uns['cellucid']['gene_sets']`, and one true/false var column, `cellucid_gene_set__<name>`, true on the genes the set lists |
| Filters, colours, views, the camera, analyses | Not written: they are how the data was shown, not the data |
| The save itself | One record added to `uns['cellucid']['history_applied']`, with the package `web`, the step applied and the SHA-256 of the "Current state only" session file of that state |

Everything else is copied as the dataset stores it: every cell, whichever
cells are open, every layer, `raw`, every embedding and graph, every value
with its type. A dataset of several modalities (a MuData, or an R object with
several assays) is saved as a Cellucid-ready MuData: the labels and
highlights go to its global obs, and the gene set columns to the var of the
modality the data mapping reads. It works for every format Cellucid opens.

A few things to know:

- **Cells not opened.** When only part of the cells are open (a subsample),
  a cell that was never opened has no known value: it is missing in a
  field's column and `<NA>` in a highlight's, never "not selected".
- **Gene set columns.** They go on the var whose genes the data mapping
  reads (`X` or a layer). When it reads the raw matrix's genes, the sets are
  kept in `uns['cellucid']['gene_sets']` only. A gene the dataset does not
  have is left out of the column. A var column of that name that Cellucid
  did not write is the dataset's own and is never replaced.
- **The cells' order.** The cells stay in the order the file stores them. A
  file Cellucid prepared keeps its priority order, so put the cells back in
  their original order after reading it (below).
- **A copy of a copy.** What the file holds already is not written again. A
  save that adds nothing says *nothing new applied* and adds no record.

**Where the file goes.** When the browser offers a save dialog (Chrome and
Edge do), you choose the name; it starts as "<dataset> with this
state.zarr.zip", and the file is written as it is made. Otherwise, and
always when an agent or a notebook asks for the copy, the file is
downloaded once it is whole. Closing the save dialog without a name saves
nothing.

**While it runs**, a notification shows how much is read, with **Cancel**.
It ends with one line, such as *Saved PBMC with this state.zarr.zip (412 MB):
5 changes applied, 2 items left out (see Details).* **Details** lists
everything left out, one line each: filters, colours and the other parts of
the state that are not data, a gene set's genes missing from the dataset, and
anything of the dataset the browser cannot copy.

**Cancel** stops it at once. The file you named is abandoned, never finished,
so no partial file is left under that name. A copy meant for a download is
removed from the browser's storage, and nothing is downloaded.

**It is refused, before anything is written**, with one sentence, when:

- no dataset is open;
- a gene set or a field has a name a file cannot hold (a `/` in it, for
  instance): rename it, then save;
- the dataset's matrix is still being read (a MatrixMarket folder streaming
  in): save once every cell's values are in;
- you choose the file the dataset is open from: a copy is always a new file;
- a change would write, rename or replace `obs['cellucid_original_position']`
  (each row's original position in a file stored in priority order), for
  instance a field renamed to that name: the copy keeps that column as the
  dataset holds it;
- another copy is being saved: wait for it, or cancel it.

**Read it in Python.** `anndata.read_zarr` reads a `.zarr.zip` through
zarr's `ZipStore` (a MuData: `mudata.read_zarr` the same way):

```python
import anndata, cellucid, zarr

adata = anndata.read_zarr(zarr.storage.ZipStore("PBMC with this state.zarr.zip", mode="r"))
adata = cellucid.restore_original_order(adata)   # a file stored in priority order: cells back in their order
adata.write_h5ad("PBMC with this state.h5ad")   # the same data as one .h5ad
cellucid.history_applied(adata)[-1]["package"]   # 'web'
```

**Show it from R** with `cellucid_view("PBMC with this state.zarr.zip")`.

Agents, Python and R save a copy with the action `history.anndata.save`
(`viewer.history_anndata_save(file=...)` in Python,
`cellucid_history_anndata_save()` in R). An agent asks you first: it writes a
file ({doc}`../g_agents/03_permissions_and_what_agents_see`).

---

## Load a shared history

**State › Session › Load… › Open a session file…** reads the file and adds its
trackers **beside yours**: nothing you have is replaced, and a tracker whose
name you already use gets "(shared)". Its step is shown as soon as its values
are read. A notification says *Loaded "T cell analysis" from analysis.cellucid-session
beside yours.* with **Undo**, which removes the loaded trackers to Recently
deleted.

A line above the graph says where the tracker in use comes from, until you
dismiss it (`✕`):

- `Shared · "T cell analysis" · 3 branches`, with **▶ Tour** and **Methods**;
- `Shared as a state, without its history · "T cell analysis"` for a file of
  the current state only.

A file saved on a sample, a link or a server opens its dataset first. A file
of another dataset is refused with the reason
({doc}`../l_sessions_sharing/10_troubleshooting_sessions`); its steps can
still be replayed on the dataset shown (below). A file of steps only opens as
a replay. The saving and loading details, refusals and limits are in
{doc}`../l_sessions_sharing/index`.

---

(history-replay-another-dataset)=
## Replay a history on another dataset

A replay applies the steps of a saved history again, **by name**, on the
dataset shown, in a new tracker: colour by the same field, hide the same
categories, keep the same genes. It is how you run the same workflow on your
next sample. Start one from:

- **Load… › Replay here** on a saved state of this browser (of this dataset
  or another);
- **Load… › Replay a session file here…** for a session file;
- **⋯ › Replay a saved history…** in the **History** header, which opens
  **Load…**.

The branch replayed is the history's main branch. A notification shows the
progress, with **Cancel**; a stopped replay keeps the steps applied before it.
Then one sentence says how it went (`Replayed 16 steps: 14 applied; 2 steps
did not match.`), and the report opens by itself when a step did not match.

**What a step does on the other dataset:**

| Step | On the other dataset |
|---|---|
| An action whose names exist there | Applied again (`✓`), or changed nothing because the state already was so (`=`) |
| An action naming a field, category or gene the dataset lacks, or cells by position | Kept in the new tracker, flagged with the reason (`✗`, *Did not match here: …*); it changes nothing |
| A step that holds only values (changes made while paused, merged steps, a fresh start, a shared state) | Flagged: its values belong to the dataset it was recorded on |
| A change made outside the history | Skipped (`–`) |
| A step recorded by a newer Cellucid | Flagged (*Recorded with a newer Cellucid*) |
| A step flagged where it was recorded | Kept, flagged |

Nothing is put back from the source: every step is applied again from what
the dataset shows. The new tracker's line above the graph says
`Replayed from Sample A: every step matched`, or how many did not match, with
**Report…**.

**The replay report** lists each step replayed with its mark (`✓` applied,
`=` changed nothing here, `✗` not applied with why, `–` skipped, `■` stopped);
clicking a step that was not applied shows it in the graph. The report is kept
with the tracker, so a saved and loaded replay reports the same.

**Use other names.** A step that did not apply because a field, a category or
a gene it names is not here offers **Use instead ▾**: the fields, the field's
categories, or the genes that match what you type, as this dataset has them.
A name missing in several steps is chosen once. Then:

- **Replay again with these names** shows where the replay started and
  replays the same steps from the same source (the saved state or the session
  file) in a new tracker, which keeps the names chosen (with any it was given
  before) and reports again. The first replay's tracker stays beside it:
  remove it when you no longer need it;
- **Keep as is** closes the report and keeps the replay as it is.

Where the replay came from (the saved state's id, or the session file's
name) is kept with the tracker. **Replay again** reads the steps from that
source: a saved state of this browser works any time; a session file only
while this page still holds it, so after a reload replay the file again from
**Load… › Replay a session file here…**. Only fields, categories and genes
are offered; a view named by a step is not mapped yet. The names chosen are
used for this replay only: the next sample starts without them.

From Python and R, `viewer.history_replay(saved=...)` and
`cellucid_history_replay()` replay a saved state too, and take a `mapping` of
other names to use for fields, categories and genes.

---

(history-apply-branch-here)=
## Apply a branch's steps here

**Apply its steps here**, in a branch's menu, applies that branch's own steps
(since it parted) again, by name, on the step shown, as a new branch: a
cherry-pick of a line of work onto another. It is a replay within the same
dataset, with the same marks for steps that do not apply.

---

(history-verify-branch)=
## Verify a branch

**Verify**, in a branch's menu, replays the branch from its start, on this
dataset, into a scratch tracker that is not kept, and compares every part of
the state at its end with the branch's last step (times and generated ids
aside). It says either `“Main” reproduced exactly (14 steps).` or which parts
differ and from which step (`Not reproduced: Highlights differs from step 9
(“…”).`). Use it before you publish a workflow.

The replay runs on what you see, and at its end you are back on the step you
were on. A step that could not run at all, because its field could not be read
(the host is failing, or the store changed on its host), gives no verdict: it
says `“Main” was not verified, because step 3 (“Coloured by sex”) could not
run: …`, and nothing else changes. Verify again once the data can be read.

---

## Next

- {doc}`../l_sessions_sharing/index`: session files in depth: what they hold,
  identity checks, troubleshooting.
- {doc}`04_tour_methods_and_report`: tour, pseudocode, methods paragraph,
  report.
