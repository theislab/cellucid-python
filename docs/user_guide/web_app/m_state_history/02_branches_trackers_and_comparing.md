# Branches, trackers and comparing

**Audience:** everyone
**Time:** 10 minutes
**What you'll learn:**
- How a branch starts, and how to go back to one
- The branch picker, the branch map, and what you can do with a branch
- How to look at two branches side by side, and what differs between any two steps
- When to start another tracker

---

## Branches

You never lose work by undoing. If you undo a few steps and then do something
new, the steps you undid stay on a branch of their own, and the new steps go
on a new branch.

Example: you colour by `cell_type`, hide the Ductal cells and keep 500 Beta
cells. You undo the last two steps, then hide the Alpha cells instead. The
history now has two branches:

```text
Opening ─ Coloured by cell_type ─┬─ Hid Ductal of cell_type ─ Kept 500 cells of Beta in cell_type   (Main)
                                 └─ Hid Alpha of cell_type                                           (Branch 2)
```

In the graph, the line on the left is the path to the step shown. Where
another branch leaves it, a row marked `⑂` names that branch, its number of
steps and how long ago it last changed; `▸` (or `→`) expands it to show its
steps in a lane of its own, and double-clicking it shows its last step.

- The first branch is **Main**. Later ones are **Branch 2**, **Branch 3**, …;
  a branch an agent starts is named after it and what it did
  (**Claude: Coloured by Ins1**).
- **Branch here** (in a step's details or menu, `Ctrl+Shift+B`) starts the
  next change on a new branch from that step, even at the newest step. The
  status line says so (*Next change starts "Branch 3" from step 4.*), and
  **Cancel** beside it takes it back.

A branch's name goes into the pseudocode, the methods paragraph and the
report, so a short, telling name ("Without Ductal") helps a reader later.

### The branch picker

**Main ▾** on the status line opens the branch picker under it:

- the tracker in use, as a button (**Tracker 1 ▾**: switch, start, rename or
  remove trackers; {ref}`history-trackers`), and **Map ⤢**, the branch map (below);
- every branch with its colour, `★` for the main branch and `🔒` for a locked
  one, where it starts (`from step 3`), its steps and how long ago it last
  changed; a click shows its last step;
- `⋯` beside each branch: its menu (below);
- **Find a branch** when there are many, and the archived branches on request;
- **⑂ Branch from the step shown**: the next change starts a new branch.

```{figure} ../../../_static/screenshots/state_history/history-branches.png
:alt: The State section after two undos and a new step: the branch button reads Branch 1, and the branch picker under it shows the tracker button, Map, Main (starred as the main branch) and Branch 1 (from step 1, the branch shown), each with its ⋯ menu, and Branch from the step shown; below, the graph with Main folded into a stub where Branch 1 leaves it.
:width: 540px

After two undos and a new step: the branch picker with **Main** and **Branch 2**, the branch shown marked.
```

### What you can do with a branch

The branch's menu (`⋯` in the picker or the map, or a right-click on its `⑂`
row):

| Item | What it does |
|---|---|
| **Show** | Shows its last step |
| **Rename…** (`F2`) | A new name, unique on this tracker |
| **Compare with the shown branch** | Its last step beside the shown one, in kept views (below) |
| **Make main** | The main branch is listed first, shared by default, and is the one a replay takes; it cannot be deleted or archived |
| **Lock** / **Unlock** | A locked branch cannot be renamed or deleted, and a change made at its last step starts a new branch instead of growing it (*"Main" is locked; your change started "Branch 4".*) |
| **Archive** / **Unarchive** | Hides it from the graph and the picker without deleting it; the map lists archived branches when asked |
| **Copy pseudocode** | Its steps as pseudocode ({doc}`04_tour_methods_and_report`) |
| **Apply its steps here** | Applies its own steps again, by name, on the step shown, as a new branch ({ref}`history-apply-branch-here`) |
| **Verify** | Replays it and checks that every part of the state comes out exactly as it was kept ({ref}`history-verify-branch`) |
| **Delete…** | Deletes it (asks first, below) |

**Deleting a branch** asks first and says what goes: its own steps since it
parted from the others go to **Recently deleted**, where you can restore them;
the results computed on it go with them; a branch that starts from it keeps
the steps they share (or, with the box ticked, goes too); kept views made from
it keep their picture but can no longer be replayed; and you see where it
started if you were on it. The main branch and a locked branch cannot be
deleted. The notification after it has **Undo**.

### The branch map

**Map ⤢** in the picker, or **⋯ › Branch map…** in the **History** header,
opens every branch at once: on the left the graph with every branch expanded
(the same rows, keys and row rule as the sidebar: click selects, double-click
shows, `›` details), on the right a table of the branches with where each
starts, its steps, its last change and its menu. **Find a branch** narrows the
table, **Archived** lists the archived ones, and the boxes choose several for
**Archive chosen** or **Delete chosen…** (one question, one notification whose
**Undo** restores them all). On a phone the table alone is shown.

---

## Compare two branches side by side

**Compare with the shown branch** (or **Compare** in a branch's details) puts
the two branches next to each other in the multiview: it keeps one view
showing the branch shown at its last step, one view showing the other branch
at its last step, and lays them out side by side. A notification says which
views were kept (`Kept 2 views: Main and Branch 2, side by side.`).

- Each kept view shows what the selected view shows at that branch's last
  step: its colouring, filters and camera.
- The step shown does not move; the history stays on the branch it showed.
- The kept views and the side-by-side layout are steps of the history like
  any kept view, so **Undo** takes them away again.

Kept views are explained in
{doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`.

---

## Compare two steps

To see what differs between two points of your work, open **⋯ › Compare two
steps…** in the **History** header, **Compare…** in a step's details or menu,
or **What changed ›** beside a step's **Changed** (that step against the one
before it).

- **From** and **To** list every step of the tracker in use under its branch,
  each once, by its number (`3 · Hid Alpha of cell_type`).
- **⋯ › Compare two steps…** starts with the step shown against the step
  before it; a step's **Compare…** starts with that step against the step
  shown.
- The line under the choices says whether the two steps are on one line of
  steps or where their branches part, and how many parts of the state differ
  (`“Main” and “Branch 2” part after step 1. 2 parts of the state differ.`),
  or that both show the same state.
- Below it, each part that differs (*Filters*, *Colors*, *Fields and labels*,
  *Highlights*, *Views*, *Camera*, *Settings*, *Layout*, …) is listed with
  short lines saying what changed.

Comparing changes nothing and shows nothing: it reads the two steps as they
were kept.

```{figure} ../../../_static/screenshots/state_history/compare-two-steps.png
:alt: The Compare two steps dialog: From and To each name a step, the line under them says the steps are on one line of steps and how many parts of the state differ, and each part that differs is listed with what changed.
:width: 960px

**Compare two steps**: the step before the one shown against the step shown, with each part of the state that differs.
```

---

(history-trackers)=
## Trackers

A tracker is one whole, independent history of the dataset, with its own steps
and branches. One tracker is enough for most work. Start another when you do
something separate on the same dataset (for example, annotating clusters in
one tracker and preparing figures in another), so the two never mix.

The tracker button at the top of the branch picker (**Tracker 1 ▾**) lists
the trackers of the dataset:

| To… | Do this |
|---|---|
| Switch tracker | Choose it in the list: Cellucid shows that tracker's step, as it was |
| Start another | **Start another tracker…** (also in **⋯** of the **History** header): **From what is shown**, or **From the dataset as opened**. It is named **Tracker 2**, **Tracker 3**, … |
| Rename it | **Rename tracker…** |
| Remove it | **Remove tracker…**, then **Remove tracker** |

A removed tracker, with its steps and branches, goes to **Recently deleted**,
where it can be restored; the tracker used before it is shown. The only
tracker of a dataset cannot be removed (start another first).

A tracker belongs to its dataset. Every tracker is kept in this browser with
the rest of the history and can be shared in a session file (you choose which
trackers and branches, {doc}`06_sharing_and_replay`); the one in use is the
one shown when the history is opened again. A tracker loaded from a file, and
one a replay made, says where it comes from in a line above the graph.

---

## Next

- {doc}`03_saved_states_and_autosave`: what this browser keeps, Load…, and
  Continue.
- {doc}`05_recording_groups_notes_and_cleanup`: the recording switch, groups,
  notes, search and Recently deleted.
- {doc}`06_sharing_and_replay`: sharing, loading and replaying a history.
