# Session mental model

**Audience:** everyone (wet lab + computational + power users)  
**Time:** 15–30 minutes  
**What you'll learn:**
- What a Cellucid session file is (and what it is not)
- The layers that keep your work: the page, this browser, the session file, browser preferences
- Why a load is all or nothing
- How to choose between sharing a link, sharing a session, and sharing the dataset file

**Prerequisites:**
- A dataset loaded in the web app (demo datasets are fine)

---

## Mental model (one sentence)

A Cellucid session file is the **history of one dataset** (every tracker, every step of each, its branches, and the step shown) saved as a single `.cellucid-session` file, which **Load…** puts back beside your own trackers, at the step it showed, with every step still there to undo to.

---

## What a session is (plain language)

If you're not technical, think:

- **Dataset** = the "slide under the microscope" (the actual data).
- **Session** = the lab notebook of everything you did at the microscope, in
  order: where you looked, what you colored by, what you filtered and
  highlighted, and which page of the notebook is open.

Every change you make is one **step** of the history (see
{doc}`../m_state_history/01_the_history`). When you click **Share…** under
**State › Session**, then **Save file**, Cellucid downloads the steps you chose
as a file named after the dataset, like
`Pancreatic endocrinogenesis.cellucid-session`.

When you choose that file in **Load… › Open a session file…**, Cellucid adds
its history beside yours: the view is the one that was shown when the file
was saved, and **Undo** walks back through the steps that led to it.

:::{important}
A session file does **not** contain the dataset itself.

It can be loaded only onto the dataset it was saved on. A session saved on a
sample, a data link or a Cellucid server opens that dataset by itself.
:::

---

## The layers that keep your work (this prevents most confusion)

Cellucid keeps your work in several places. They are easy to mix up:

1) **The page (volatile)**
   - What is on screen now, with its history.
   - Lost from the page when you refresh, crash or close the tab, but see 2.

2) **The history kept in this browser**
   - Every step is written to this browser as you make it (**Autosave in this
     browser**, on by default, in **History settings**).
   - Reopening the dataset offers **Continue where you left off?**, and
     **State › Session › Load…** lists every history kept.
   - It stays in this browser and this computer. See
     {doc}`../m_state_history/03_saved_states_and_autosave`.

3) **Session file (`.cellucid-session`)**
   - The file you download with **Share…**: the history, the current state
     only, or the steps only, as you choose.
   - Intended for:
     - keeping your work outside this browser,
     - reproducibility ("these exact steps, views, filters and highlights"),
     - sharing with collaborators.

4) **Browser preferences (small conveniences)**
   - Antialiasing, and the theme and background new pages open with.
   - Not a reproducibility mechanism, and not shared.

Practical implication:
- To **share** or **reproduce**, use a **session file** (and the dataset).
- To pick up your own work on the same computer, the history kept in this
  browser is enough.

---

## Fast path (wet lab / non-technical): save your work and send it to someone

1) Load your dataset.
2) Explore until you have what you want (color by a field, filters,
   highlights, kept views).
3) Open the **State** section of the sidebar.
4) Under **Session**, click **Share…**, then **Save file**.
5) Find the downloaded file (usually in your Downloads folder).
6) Send the file to your collaborator. A session saved on a sample, a data link
   or a server records where the dataset is and opens it; for a dataset opened
   from your own computer, send the dataset file with it.

To reopen it later, click **Load… › Open a session file…** and pick your
`.cellucid-session` file. Its dataset opens first if it is not open already; a
session saved on a local file needs that file opened again first (in Chrome
and Edge, **Load…** offers
**Reopen <name>** when this browser kept the file).

---

## Practical path (computational users): what "state" actually means

What you can change is split into sixteen parts. A step records which parts it
changed, and their value before and after. A session file holds those values
for every step.

- **Parts of the whole dataset view:**
  - the cells to open and the genes read for them,
  - the data mapping (which stored arrays are plotted and read),
  - the expression values (genes normalized and/or log1p, or as stored),
  - fields and labels (renames, deleted fields, fields you made),
  - gene sets,
  - the views (kept views, their layout, the focused view, linked cameras,
    the embedding each view shows),
  - highlights (pages and groups of cells),
  - the camera path,
  - settings (every stable sidebar control, Figure Export),
  - analysis windows and the analysis panel, with their results,
  - layout (docked and floating panels, open sections, sidebar, theme).
- **Parts each view has its own of** (the live view and each kept view):
  - the field or gene it is colored by,
  - its filters,
  - its colors,
  - cells hidden one by one,
  - its camera and navigation mode.

The full list, with what each part holds, is
{doc}`02_what_gets_saved_and_restored`.

Two rules make a session apply to the right things:

- **Things are named, not numbered.** A field is saved by the key the dataset
  gives it, a category by its label, a gene by its name, a view by its place
  (the live view first, then the kept views in order).
- **Cells are named by original position**: the cell's row in the dataset as
  it was written. That is the same however the cells were opened (all of them,
  a share, grown step by step), which is why the dataset's cell order is
  checked before anything is loaded (below).

### One load, one outcome

A load is all or nothing:

1. The session's dataset is opened if it is not open, with the cells the
   session was saved on.
2. The file and the dataset are checked (below) before anything changes.
3. As soon as the steps and the values of the step shown are read, that step is
   shown.
4. The other steps' values are read while it is shown; a new step waits for
   them.
5. The load ends in a notice reading **Session fully restored.**

If anything fails, or you press **Cancel**, Cellucid puts back what was shown
before and keeps the trackers you had. Loading another session while one is
loading cancels the first.

### Dataset identity (the first restore check)

Each session carries a **dataset fingerprint**: the dataset id (and its
display name, used only in messages), cell count, gene count, a
**cell-order record**, and the **data mapping**, which stored embeddings and
arrays the dataset was shown with. Beside it the session lists every
observation field and gene the shown step uses, so a dataset missing one is
refused before anything changes, naming every missing key.

The first three say *which* dataset. The cell-order record says *which ordering
of its cells*, which the first three cannot: a session names cells by their
position in the dataset's own order, and a dataset written again under the same
name with the same number of cells and genes but its rows in a different order
matches all three while every position now points at a different cell. The
record is a short digest of the stored coordinates of one embedding, each cell
named by its position, plus that embedding and its dimension (a digest only
means anything within one embedding). Cellucid checks the digest on that
embedding whichever one is on screen when you load; how the cells were opened
does not change it.

If any part of the fingerprint does not match the open dataset:
- Cellucid refuses the whole load;
- no part of the session is put back; and
- what was shown stays as it was.

You will see one of these messages, each naming what actually happened:

| What you see | What it means | What to do |
|---|---|---|
| The session was saved on a different dataset (with both cell and gene counts named) | it is not this dataset | open the dataset the session was saved on |
| The session was saved on a different dataset with the same number of cells and genes (naming both datasets by the names they show) | the sizes agree exactly and the identity does not: two datasets, or one dataset reached two ways under two ids | open the dataset the session was saved on |
| The session uses an embedding (or neighbor graph, modality, gene names, expression matrix, expression transform, highly variable genes, latent space, or vector field) this dataset does not have or cannot use, naming it | it was saved on another version of the dataset (a saved choice this dataset does have is put back, not refused) | open the dataset version the session was saved on |
| The dataset has the same name and counts but its cell coordinates differ | either the cells are listed in a different order, or the same cells were exported again from a re-computed embedding; a one-way digest cannot say which | load the version the session was saved on, or re-create the selections here |
| The session uses observation fields or genes the dataset does not have, naming every one | it was saved on another version of the dataset | open the dataset version the session was saved on |
| The session was saved by an earlier version of Cellucid (session format 3, or older), which this version no longer reads | a session file now holds the history of every step (session format 4); an earlier one held one saved state | open it in the version of Cellucid that saved it, or re-create the session in this version |
| The file is not a Cellucid session file of format 4 | it does not start with the format-4 signature `CELLUCID_SESSION/4` | check that you picked a session file saved by this release |

The format version is part of every session file's signature, and this
release reads format 4 only. The Python and R packages read the same format
and refuse an earlier one with the same sentence; a file written by a newer
Cellucid than the package, with records it does not know, is refused as newer
(update the package), never as damaged.

See {doc}`07_versioning_compatibility_and_dataset_identity` for the exact
matching rules and safe collaboration recipes.

---

## Deep path (power users / developers): what's in the file (high level)

A `.cellucid-session` is a single file that contains:

- a JSON **manifest**: the dataset fingerprint, where the dataset opens from,
  the cells it was saved on, the fields the shown step uses, which tracker is
  active, and each tracker's root, step shown, branches and names; and
- length-prefixed, gzip-compressed **chunks**: the steps of each tracker, the
  values of the parts at those steps, and the typed arrays those values hold.

A value or an array is stored once, however many steps refer to it: a
highlight group's cells are written once even when fifty later steps keep it.
There is no separate encoder per feature: a part's value is stored as it is.

Sessions are treated as **untrusted input**:
- the manifest, every chunk and every step tree are checked exactly;
- every chunk is decompressed to exactly its declared size, within bounds;
- every array is checked against its digest and every value against its
  hash; and
- a failure anywhere puts back what was shown.

If you want the exact reference, see {doc}`12_reference`.

---

## When to use sessions vs other artifacts

Use this decision rule:

- **A session file**: when you want to keep the interactive work, with its steps
  (views, filters, highlights, analyses).
- **Exported figure**: when you want a publication artifact that does not
  depend on the app to interpret.
- **Methods paragraph or printable report**: when you want the steps in words
  for a lab notebook or a paper ({doc}`../m_state_history/04_tour_methods_and_report`).
- **Community annotation**: when you want many people to label and vote
  together without conflicts (GitHub-backed).

Sessions are great for:
- "Please look at this exact selection and camera angle."
- "Here are the highlight groups I think are interesting, and how I made them."
- "Here's the state right before I ran DE; can you verify?"

Sessions are not ideal for:
- multi-user, long-running labeling workflows (use community annotation),
- distributing the dataset itself (share the dataset file or hosting link).

---

## What changes when you load a session

Loading a session is intentionally "strong":

- It **adds the file's trackers beside yours** (a name you already use gets
  "(shared)") and shows the step the file showed; the notice after it has
  **Undo**, which removes them again. Opening a saved state from **Load…**
  instead replaces the history shown.
- Every part is put back exactly, empty ones included: an empty camera path in
  the file clears the one on screen, and analysis windows the file does not
  have are closed.
- The cells the session was saved on are opened first; a dataset open with
  other cells is read again.
- A camera path saved with **Autoplay** on starts playing once the load is
  done.
- Loading another session while one is loading cancels the first.

Nothing you had is lost: your trackers stay beside the loaded ones, and while
**Autosave in this browser** is on the history is kept in this browser as
before. To keep it outside this browser, **Share…** it first.

---

## Interface reference

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The State section with the Pancreatic endocrinogenesis sample open: Session with Share… and Load…, and the History with the recording switch, the branch shown and the graph of the steps taken so far, the step shown marked.
:width: 540px

The **State** section holds the session's **Share…** and **Load…** above the
history: the steps a session file keeps, with the step shown marked.
```

---

## Next steps

- {doc}`02_what_gets_saved_and_restored` (exact list of what is and is not saved)
- {doc}`03_save_restore_ux` (**Share…** / **Load…** step by step)
- {doc}`05_share_workflows_links_sessions_data_files` (how to share with collaborators)
- {doc}`10_troubleshooting_sessions` (if anything goes wrong)
- {doc}`../m_state_history/index` (the history: steps, undo, branches, trackers)
