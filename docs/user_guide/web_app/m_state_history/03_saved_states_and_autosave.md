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

- its name: **Autosave** for a history kept as you worked, or the name you gave
  it, with a badge (**Auto** or **Named**) and its size;
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
:alt: The Load dialog: under Saved in this browser, a state saved under the name "Without Ductal" and the autosave of the Pancreas history, each with its kind, size, dataset, steps and when, Open (or Shown), Replay here, rename and delete; the name field, the History / Current state only / Steps only choice and Save; and under From a file, Open a session file… and Replay a session file here….
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
