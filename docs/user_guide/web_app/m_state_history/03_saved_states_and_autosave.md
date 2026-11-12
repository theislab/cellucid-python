# Saved states and autosave

**Audience:** everyone
**Time:** 10 minutes
**What you'll learn:**
- How every step is kept in this browser as you work
- How to continue where you left off
- How to keep a version under a name, and find it again in **Load…**
- How to turn autosave off, and what this browser then keeps

---

## Every step is kept in this browser

While **Autosave in this browser** is on (it is on by default; the switch is in
**History settings**, {ref}`history-settings`), each change of the history is
written to this browser's storage a moment after you make it: a new step, an
undo or a jump, a renamed branch, a new tracker, a deleted step, a note, a
conversation with an agent, a session file you loaded. A crash, a reload or a
closed tab loses nothing that was written. What it can lose is the change of
its last moment: one still on its way into storage, a fraction of a second
after you made it (tens of milliseconds on a laptop). What is written is
always whole: a write a closing tab cuts short is not kept at all, so a saved
state never holds half a step or half a conversation. Changes made while
recording is paused are kept too, and come back as one step.

- Only what changed is written, and a state that several steps share is
  written once.
- A history is kept once it holds something of yours (a step, a note, a
  conversation, changes made while paused, or a history you loaded): opening a
  dataset and looking at it writes nothing.
- Each opening of a dataset keeps its own history. Continuing an earlier one
  (below) keeps writing into that one.
- One tab writes a history at a time. A second tab that continues the same one
  opens a copy of it.
- Nothing is uploaded. The history is kept in this site's storage in this
  browser profile only (the IndexedDB database `cellucid-history`; see
  {doc}`../o_accessibility_privacy_security/02_privacy_model`).

---

## Continue where you left off

When you open a dataset whose earlier history this browser kept, a
notification asks:

> Continue where you left off? Kept in this browser from 3 Nov 2026, 14:05 on this dataset: 12 steps and 1 conversation.

- **Continue** puts that history back, at the step it showed, with every step
  to undo to, and keeps writing into it. Changes made while recording was
  paused come back as one step, and the notification says so.
- **Start fresh** keeps the history you have now; you are not asked about that
  earlier history again. It stays in **Load…**.

Nothing is put back unless you choose **Continue**. The offer names the most
recent kept history of this dataset that holds something of yours. When that
history is open in another tab, the offer is **Open a copy here** instead.

**A part that cannot be put back.** A saved state (and a session file) is put
back part by part: the views, the fields, the highlights, the filters, the
Expression values and each other part of the state, the history's steps,
and its conversations. A part that cannot be put back is left as it was
shown, and the rest is restored. The notification says *Session restored but
for 1 part* and names each such part in a sentence of its own, with the
failure's own words under **Details**. A part saved by an earlier (or a
newer) version of Cellucid that this version cannot read is said as that:
*Expression values could not be restored: it was saved by an earlier
version of Cellucid (0.9.7), and this version cannot read this part.*
The history then holds one more step, *Not restored from the saved state:
… (kept as shown)*, so what you see and the history agree; undoing it tries
the saved part again. A tracker whose steps this version cannot read is left
out and named; when it is the one that was shown, what it showed is
restored as the start of a new tracker of that name, with its
conversations. A conversation this version cannot show is named the same
way and stays with the history. The saved state you continued is then not
written into: it keeps every part as it was, and your work goes on as a new
saved state beside it. An agent opening a saved state (`history_saved_open`)
is told the same sentences (`not_restored`).

```{figure} ../../../_static/screenshots/state_history/continue-offer.png
:alt: The notification Cellucid shows when a dataset whose history this browser kept is opened again: "Continue where you left off? Kept in this browser from … on this dataset: 3 steps.", with Continue and Start fresh.
:width: 760px

The offer to continue, shown when the Pancreas sample was opened again after three steps.
```

You can also continue from the list of your datasets: **Session › Open data ›
My datasets › All…** shows each dataset's saved states, and its **Continue**
opens the dataset with its latest one
({doc}`../b_data_loading/01_loading_options_overview`).

---

## The Load dialog

**State › Session › Load…** lists, under **Saved in this browser**, every
history this browser keeps, for every dataset, newest first. Each row shows:

- its name: the name you gave it, or, for a history kept as you worked, what
  its latest step did ("Colored by gene CD33"), so autosaves read apart; with
  a badge (**Auto** or **Named**) and its size;
- its dataset, how many steps, and when it last changed (hover for when it was
  first saved).

| To… | Do this |
|---|---|
| Open one | **Open**. If it is of another dataset, that dataset opens first, as for a session file (a file from your computer is opened through the handle this browser kept, or you are asked to choose it). The history shown now is marked **Shown**. |
| Run its steps on the dataset shown | **Replay here**: its steps are applied again, by name, in a new tracker ({doc}`06_sharing_and_replay`) |
| Rename one | The edit button |
| Delete one | The bin button, then **Delete** in its row (**Keep** changes your mind). The line above says how many steps of which dataset would go. It cannot be undone. |
| Keep the history shown under a name | Type a name in *Name for the history shown*, choose what to keep (**History**, **Current state only** or **Steps only**), then **Save** (or `Enter`) |

Agents, Python and R reach the same list: `history.saved.list` lists the
saved states with their ids, and `history.saved.delete` deletes one by its id
(`viewer.history_saved_delete(saved=...)` in Python,
`cellucid_history_saved_delete()` in R). The history shown stays open, and a
deleted saved state cannot be brought back, so an agent asks you first
(**Deleting what was saved**,
{doc}`../g_agents/03_permissions_and_what_agents_see`).

An agent asked to continue your saved work opens it as **Open** does
(`history.saved.open`, by its id): its dataset first when another is shown,
then its history, which goes on where it left off. A whole-file download of
that dataset is asked about first (**Large downloads**). When it cannot open
it, it says why: a file of your computer the browser lets only you read
again (the page offers **Reopen** *file*, after which it opens), a saved
state this version cannot open, or what failed. Applying its steps to the
dataset shown instead is **Replay here** (`history.replay`). A notebook's
viewer shows the dataset its notebook serves, so it refuses
`history_saved_open()` from Python (`cellucid_history_saved_open()` in R), as
it refuses opening another dataset.

Below the list, **From a file**: **Open a session file…** adds a session file's
trackers beside yours, and **Replay a session file here…** applies its steps
on the dataset shown ({doc}`06_sharing_and_replay`).

A **named** save is the history as it was when you saved it; it does not
change afterwards. Opening it later continues from it as a new history, so the
named save stays as it was, and the history you had before stays in its own
autosave, listed beside it. Opening an **Autosave** keeps writing into it.

A saved state of the **current state only** or of **the steps only** keeps what
the choice holds: the current state only is the state shown as one step,
without how you got there; the steps only keeps no cells or values at all and
opens as a replay. Each row's size is what that saved state holds, each value,
part and array counted once, and it goes down when steps are deleted.

Saved states **imported** from another browser (**My datasets › All… ›
Import…**, a `.cellucid-library` file) are added as named saves, each called
"<name> (imported)", beside the ones kept here; nothing kept here is
replaced, and importing the same file twice keeps two copies. Open one as
any named save.

```{figure} ../../../_static/screenshots/state_history/saved-states-dialog.png
:alt: The Load dialog: under Saved in this browser, a state saved under the name "Without Ductal" and the autosave of the Pancreas history, named by its latest step "Colored by clusters_coarse", each with its kind (Named, Auto), size, dataset, steps and when; Open and Replay here on the named state, Shown on the autosave, and rename and delete on both; the name field, the History / Current state only / Steps only choice and Save; and under From a file, Open a session file… and Replay a session file here….
:width: 960px

**Load…** with the autosave of the history shown and a state saved under a name.
```

Where the browser gives the site no storage for this (some private
windows), nothing is kept: the list says *Nothing is kept in this browser.*,
and a session file can still be opened.

---

## Turn autosave off

The switch **Autosave in this browser**, in **History settings** (the sliders
button in the **History** header), is on by default.

- **Off:** nothing more is written. Saving under a name in **Load…** still
  works, and what was kept before stays until you delete it.
- **Work not kept is never dropped without a question.** While the history
  shown holds steps made since you last saved it (a named save in **Load…**,
  or **Save file** of the whole history) and this browser is not keeping them
  (Autosave off, a browser that keeps nothing for the site, or no room left),
  opening another dataset (the sample list, My datasets, a file, a link, a
  server, a saved state or session file of another dataset, **None**, the
  benchmark) first asks *Your work is not kept*: **Save…** opens **Save
  file** (open the dataset again once it is saved), **Open anyway** drops the
  steps, **Cancel** keeps the dataset shown. Steps that only arranged the
  page (a section opened, a panel moved or docked) are not counted: losing
  them is no loss to ask about. An agent opening a dataset is
  asked about through its permission card instead (**Replacing work that is
  not kept**, {doc}`../g_agents/03_permissions_and_what_agents_see`). Reloading
  or closing the tab gets the browser's own leave-page question, in the same
  case only. Reading the same dataset on other cells (**Cells to open**) keeps
  the history and asks nothing.
- **On again:** the whole history shown is written at once, then every step as
  before.

The switch is a setting of this browser, not of your work: it is not saved in
session files and does not travel with them.

---

## When the browser runs out of room

If the browser refuses to keep more, Cellucid stops keeping the history and
says so once:

> This browser has no room left to keep the history, so the steps from now on are not kept. Remove saved states you no longer need (Load… in State › Session, or Stored in Session › Stored); keeping resumes then.

Delete saved states you no longer need in **Load…**: keeping starts again, in
every tab, and writes everything it missed. **Session › Open data › Stored ›
Manage…** lists the saved states with their sizes beside everything else this
browser keeps for Cellucid, which helps you see what takes the room; its
**Keep it when space is low** asks the browser not to clear the site's data
when the disk runs low.

---

## Saved states or a session file?

| | Saved state | Session file |
|---|---|---|
| Made | Automatically at every step, or with **Save** under a name in **Load…** | With **State › Session › Share…** |
| Where | This browser only | A `.cellucid-session` file you keep and send |
| Opens on | This browser | Any computer and browser, with the same dataset |
| Holds | The whole history (or the state or steps only, when saved so) | What you chose in **Share…**: the trackers and branches, notes, conversations, Recently deleted, results and times |
| Removed | When you delete it, or clear the site's data | When you delete the file |

Use saved states to never lose work, and a session file to share it or keep it
with your project. See {doc}`06_sharing_and_replay` and
{doc}`../l_sessions_sharing/index`. To take saved states to another browser of
your own, export them with your datasets from **My datasets › All… ›
Export…** ({doc}`../b_data_loading/01_loading_options_overview`).

---

## Next

- {doc}`04_tour_methods_and_report`: tour, pseudocode, methods paragraph, report.
- {doc}`06_sharing_and_replay`: sharing, loading and replaying a history.
