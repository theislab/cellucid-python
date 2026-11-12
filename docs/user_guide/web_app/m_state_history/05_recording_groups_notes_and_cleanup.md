# Recording, groups, notes and cleanup

**Audience:** everyone
**Time:** 10 minutes
**What you'll learn:**
- How to pause and resume recording, and what happens to the changes made meanwhile
- What History settings decide
- Groups and merging steps into one, and when to use which
- Notes, bookmarks and search
- Recently deleted: undo for the history itself
- Clearing the history while keeping the state, and showing the state at a time

---

(history-recording-paused)=
## Pause and resume recording

The **● Recording** switch at the start of the status line pauses the
recording. While it reads **⏸ Paused**:

- every change still applies to the view, but makes no step;
- the graph shows a **Paused** row after the step shown, and the status line
  says since when (`since 10:42`);
- the changes made meanwhile are kept in this browser a moment after each is
  made, so a crash or a reload keeps them (they come back with **Continue**,
  {doc}`03_saved_states_and_autosave`);
- an agent asks before every change, whatever your permission settings allow,
  because such a change cannot be undone on its own.

**Resume** (the same switch) turns everything that differs from the step
shown into **one step**, marked `◐`, that stands for "Changes made while
recording was paused", with everyone who acted meanwhile. Saving or sharing a
session file, loading one and opening a saved state do the same first, so a
file never holds unfinished changes.

**Going back while paused.** Undo, redo, showing another step, switching
branch or tracker and the tour would lose those changes, so Cellucid asks
first:

> **Keep your changes first?** You changed Filters, Colors while recording was paused. Keep and go makes them one step, then moves. Discard and go drops them.

**Keep and go** (the default) records them as one step, then moves; **Discard
and go** drops them; closing the question stays where you are. Tick **Don't
ask again while paused** to keep the answer until you resume, or turn **Ask
before going back while paused** off in History settings to always keep them.

With **Record history for new datasets** off (History settings), a dataset you
open starts paused.

---

(history-settings)=
## History settings

The sliders button in the **History** header opens **Record**, the history's
settings in this browser. Each row has an **i** button saying what it does,
and the **History** header's own **i** tip links here (**History settings…**).

| Setting | On (the default) | Off |
|---|---|---|
| **Camera moves** | Every camera move (and the navigation mode) is a step | Camera moves make no step and undo never moves the camera; the steps made meanwhile still keep the camera, and their rows dim |
| **Panels, windows and theme** | Docking, floating and resizing panels and analysis windows, opening sidebar sections and the theme are steps | They make no step and undo never moves them |
| **Record history for new datasets** | A dataset you open starts recording | It starts paused |
| **Ask before going back while paused** | The question above | Changes made while paused are kept as one step without asking |
| **Autosave in this browser** | Every step is kept in this browser as you make it ({doc}`03_saved_states_and_autosave`) | Nothing more is written; named saves still are |
| **Your name in shared histories** | Empty: nobody is named | Your name, recorded with your steps and notes and shown in a shared file |

Everything else is always recorded: it is the analysis itself. These settings
are kept in this browser (`cellucid:history-settings`), never in a session
file or a step. The picture switch of the printed report is in the report's
own dialog ({doc}`04_tour_methods_and_report`).

---

## Groups and merging steps into one

Two ways to tidy a run of consecutive steps of one branch. Select the first
row, `Shift`+click the last, then right-click (or press `Ctrl+G` for a group):

| | **Group…** (`▣`) | **Merge into one step…** (`⧉`) |
|---|---|---|
| What it does | Names the steps and shows them as one node you can expand | Replaces them with one step, from the state before the first to the state after the last |
| The steps | All kept, each still undoable and showable | The steps in between go to **Recently deleted**, where you can restore them |
| Notes and bookmarks | Stay on their steps; the group has notes of its own | Move to the merged step |
| On another dataset | Each step replays | The merged step brings its values only, flagged |
| In the pseudocode and methods | A section with the group's name | One line that holds only values |

Groups are the safe choice: nothing changes or is lost. A group's row says
where the step shown is among its steps (`at 4 of 7`). A group holds whole
turns of an agent (a selection that would split one is widened to the turn,
and the question says so) and never another group. Its details show its steps,
when, and branches that start inside it, with **Show last step**,
**Rename…**, **Ungroup** and **Merge into one step…**. **Add to group "QC"**,
in the menu of a step after a group, grows the nearest group before it to
that step.

**Merge into one step…** asks first and offers **Group instead**. It is
refused, with the reason, when a branch starts inside the run (except at its
last step), and when the run would cut a group (merge inside it, or all of
it). A turn's and a group's details offer it for their steps.

Every group edit and every merge comes with a notification whose **Undo**
reverses it, and each is listed in **Recently deleted**.

---

## Notes and bookmarks

- **Notes.** A step, a branch, a group, a turn and a tracker each take notes
  (**Add note** in the details; an agent's turn calls its first note its
  **Summary**). A note is plain text with a little Markdown (`*emphasis*`,
  `` `code` ``, `[links](https://…)`); `Ctrl+Enter` saves it. Each note keeps
  who wrote it and when, and can be edited or removed. The tracker's first
  note is the abstract of the methods paragraph and the report; a step's
  first note is its caption's second line in the tour.
- **Bookmarks.** **Bookmark…** in a step's menu names the step (`★ name`), to
  find it by: a bookmark is unique on its tracker, a jump can name it, and the
  printed report pictures bookmarked steps by default. An empty name removes
  it.

Notes and bookmarks are shared in a session file only when you tick **Notes
and bookmarks** in **Share…** (it is ticked by default).

---

## Search and filter

The search button in the **History** header (or `/`) opens **Steps, notes,
bookmarks**: the steps whose words, notes or bookmark match stay, the others
fold into runs (`12 steps hidden`), the step shown always stays, and the
status line counts the matches (`8 of 214 steps`). A branch that holds a match
is expanded.

**Filter** narrows further, by who acted (**By you**, **By an agent**, **From
Python**, **From R**, **Replayed**), by kind (**Analysis changes**, **Views and
camera**, **Panels and layout**), by **Part changed…** (each part of the state
some step changed, with its count of steps) and **Day…** (each day a step was
taken, newest first), and **With a note**, **Bookmarked** or **Stale or
flagged**. `Escape` closes the search.

---

(history-recently-deleted)=
## Recently deleted

Ctrl+Z undoes changes to the data and the view. **⋯ › Recently deleted…** is
undo for the history itself: every edit of the history that can still be
undone, newest first, across the trackers of the dataset (deleted steps and
branches, removed trackers, merges, groups, renames, notes, loaded files),
each with who did it and when.

- **Restore** puts it back exactly. When later edits changed what it made, it
  comes back as new steps or branches named "… (restored)", applied again by
  name.
- `✕` empties one entry, and **Empty…** all of them, after a question: they can
  no longer be restored, and the memory they hold is freed.

Each history edit also shows a notification with **Undo**, the same as
**Restore** on its entry; a new edit's notification replaces the previous
one. Recently deleted is kept in this browser with the history, and shared in
a session file only when you tick **Recently deleted** in **Share…**.

---

## Clear the history, keep the state

**⋯ › Clear history, keep the state…** starts a new branch, **Fresh start**,
holding what is shown as one step from the dataset as opened (`⧉`). Every
branch stays (archive them when you like), and the clear is undone from
Recently deleted. Use it before sharing when the way you got there does not
matter, or share the current state only ({doc}`06_sharing_and_replay`).

---

## Show the state at a time

**⋯ › Show the state at…** asks for a time (`10:30` today, or a date and time
such as `2026-11-04 10:30`) and shows the last step on the branch shown at or
before it.

---

## Next

- {doc}`04_tour_methods_and_report`: tour, pseudocode, methods paragraph, report.
- {doc}`06_sharing_and_replay`: sharing, loading and replaying a history.
