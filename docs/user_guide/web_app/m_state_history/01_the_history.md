# The history

**Audience:** everyone
**Time:** 10 minutes
**What you'll learn:**
- What a step is, and what is not one
- How to read the history graph: its rows, marks and numbers
- How to show any earlier step, and what its details tell you
- Undo and redo, with the buttons and the keyboard
- How to delete, revert or leave out one step and keep the rest
- How steps taken by an agent or a notebook look

---

## What a step is

Every change you make to the dataset shown becomes one step of its history:

- colouring by a field or a gene, and every colour setting (a category's
  colour, a colour map, log scale, the colour range);
- a filter: hiding categories, a range, an outlier filter, missing values,
  keeping only some cells of a category, hiding the categories with fewer
  cells than a count ({doc}`../e_filtering/index`);
- labels: renaming or merging categories, renaming or deleting a field, a new
  field made from highlights or written by a notebook;
- highlights: a confirmed selection, a highlight page made, renamed or
  recoloured;
- views: a kept view, the multiview layout, the selected view, linked cameras,
  the embedding a view shows;
- the camera: one step per move, so a drag, a scroll and a pan in a row on the
  same view become one step;
- settings in the sidebar (point size, background, navigation, overlays): one
  step per setting, so dragging a slider is one step;
- the figure export settings, and the camera path (its keyframes, duration,
  loop and playback);
- the layout: docking or floating a panel, an analysis window's size and
  place, opening or closing a sidebar section, the theme;
- an analysis: opening its window, its options, its view and plot, running
  it, closing it;
- **Cells to open** in **Session**: reading the dataset again on another share
  of the cells ({doc}`../b_data_loading/21_subsampled_loading`). Undoing it
  reads the dataset again on the cells open before.

The first row, the **opening**, is the dataset as it opened, with its starting
view if its catalog names one.

Not a step:

- something that changes nothing (choosing what is already chosen);
- a selection while you draw it: the selection tools keep their own undo until
  you confirm, and the confirmed selection is the step. Which selection tool
  is chosen is never a step either; a session file keeps it and puts it back
  when it is loaded;
- opening another dataset: that dataset gets a history of its own;
- saving or sharing a session file, copying pseudocode, or anything else that
  only reads;
- a change made while recording is paused: it applies, and becomes one step
  when you resume ({doc}`05_recording_groups_notes_and_cleanup`);
- a camera move, or a panel, window or theme change, while its switch in
  **History settings** is off ({ref}`history-settings`).

:::{note}
Labels follow the history. Going back to a step shows the category names,
merged categories, colours and fields as they were at that step, and each
branch has its own. An analysis result comes back as it was when you go back
to it; it is not computed again.
:::

---

## The history graph

**State › History** draws the history of the dataset shown as a graph, one row
per step, from the opening at the top to the newest step at the bottom. The
line down the left is the path to the step shown; the steps of other branches
hang off it ({doc}`02_branches_trackers_and_comparing`).

Each step's row has a mark on the line, one short sentence saying what the
step did ("Coloured by cell_type", "Hid Ductal and Alpha of cell_type",
"Kept 500 cells of Beta in cell_type"), a badge naming who took it when it was
not you (an agent's name and model such as **Claude · Sonnet 4.6**,
**Python**, **R**), and how long ago (`now`, `4 min`, `2 h`, `3 d`, then the
date). Hover a row for the exact time; nothing is cut short, a long label
wraps.

| Mark | Row |
|---|---|
| `●` | a step |
| `◉` | the step shown |
| `○` | a step you undid (faded; redo brings it back) |
| `⑂` | where another branch leaves the line (expand it to see its steps) |
| `◆` | an agent's or a notebook's turn: the steps it took for one request, under it |
| `▣` | a group of steps you named ({doc}`05_recording_groups_notes_and_cleanup`) |
| `⋯` | a run of camera, view or layout steps folded into one row |
| `⏸` | recording is paused here |
| `◐` | the changes made while recording was paused, kept as one step |
| `⧉` | steps merged into one, a fresh start after the history was cleared, or a state shared without its history |

**Folded runs.** Consecutive steps that only present the analysis (camera
moves, panels, open sections, colours that do not change what is coloured by
or filtered) fold into one row such as `Camera ×6 · Layout ×2`. Folding is for
reading only: undo still goes step by step, and expanding the row shows each
step. A step with a note, a bookmark, a branch starting from it or a stale
result never folds.

**Numbers.** Steps are numbered in the order they were made, and named by
their number everywhere: in the details, the menus, Compare, the branch
picker and the search ("at step 12"). The opening is not counted.

**The status line** under the header says where you are:

| It reads | When |
|---|---|
| `● Recording  Main ▾  14 steps` | the newest step of branch **Main** is shown |
| `14 steps · 5 undone` | five steps after the one shown are undone |
| `8 of 214 steps` | a search shows 8 of the steps |
| `⏸ Paused  Main ▾  since 10:42` | recording has been paused since 10:42 |
| `⚠ 2 stale` | two steps' results no longer match their inputs (click it to see them) |
| `⚠ 1.2 GB` | the history holds more memory than this device has room for now (click it to see where it goes) |

With no step yet it says *Every change you make appears here. Undo with
Ctrl+Z.* (`⌘Z` on a Mac). **Recording** is the recording switch, and **Main ▾**
opens the branch picker ({doc}`02_branches_trackers_and_comparing`).

**Memory.** The history keeps every value of every step in memory, and
never drops one by itself. Once what Cellucid draws (the open cells, the
history included) is more than this device has room for, the status line
shows how much the history holds (`⚠ 1.2 GB`) and says once which parts take
the most (*Largest: Highlights at step 31: 8 MB; …*). The chip opens
**History memory**: every part that holds memory, largest first (a step's
values, what Recently deleted keeps, changes made while paused), each step's
with **Show**, and **Recently deleted…** to empty what it keeps
({ref}`history-recently-deleted`). It is a warning, not a limit: you keep
working.

The graph draws only the rows in view, so a history of ten thousand steps
opens and scrolls like a short one. When the step shown is scrolled out of
sight, a pill such as **Shown: step 12 ↑** brings it back.

---

## Show any step

The rows follow one rule, in the sidebar, in the branch map and on a phone:

| To… | Do this |
|---|---|
| Select a row | Click it (`Shift`+click selects every step between) |
| Show the step | Double-click it (double-tap on a phone), or select it and press `Enter` |
| Open its details | `›` at the end of the row, or `Space` |
| Expand a branch, turn, group or run | `▸`, or `→` |
| Open the row's menu | Right-click it (the menu opens at the pointer), `Shift+F10`, or `⋯` in its details |

Showing a step puts back the dataset as it was right after that step. Nothing
is lost: every later step stays in the graph, faded, to show again, and redo
walks back to them one by one. If you act from there, a new branch starts
({doc}`02_branches_trackers_and_comparing`).

**A step's details** say:

| Detail | What it says |
|---|---|
| **By** | Who took it: *You* (or your name, {ref}`history-settings`), an agent's name with its model, *Python*, *R*; **turn ↗** opens the agent's conversation |
| **When** | The time it was taken, how many changes it merged (a drag is one step), how long it took when that was long, and its branch |
| ⚠ **Stale since…** | For a step whose result no longer matches its inputs: which later step changed them, and, for an analysis or a highlight of a category or a range, **Run again** |
| ⚠ **No longer applies:** | Why a step stopped applying (a step it needed was deleted) |
| **Stands for** | For a step that holds only values: what it stands for (changes made while recording was paused, steps merged into one, a state shared without its history) |
| **Code**, or the parameters | The step's line of pseudocode (`⧉` copies it), or what it was given, by name: the field, the categories, the gene, the number of cells, the seed |
| **Changed** | The parts of the state it changed (*Filters*, *Colors*, *Camera*, *Fields and labels*, …); **What changed ›** compares the state before and after it |
| **Bookmark** | Its bookmark, `★ name` |
| **View** | **Show this step's view**: each view's camera put where the step kept it, nothing else changed (a camera step while camera moves are recorded; applied even while they are not, {ref}`history-settings`) |
| **Notes** | Its notes, with **+ Note** to add one |

Under them: **Show**, **Branch here**, **Compare…** and `⋯` (the row's menu:
**Bookmark…**, **Revert**, **Branch without this step**, **Copy code**,
**Group…**, **Merge into one step…**, **Delete…**).

```{figure} ../../../_static/screenshots/state_history/history-step-details.png
:alt: The History block of the State section after three steps, with the second step shown by a double-click: the line under the recording switch and the branch Main reads 3 steps · 1 undone, the step shown is marked in the graph, and its details, opened with ›, say who took it, when and on which branch, its code (keep_cells with the field, category, count and seed), what it changed with What changed ›, + Note, and Shown, Branch here, Compare… and ⋯.
:width: 496px

The second of three steps shown, with its details open: who took it, when, its parameters, what it changed, and its buttons.
```

**Keyboard.** The graph is one stop of `Tab`; the **i** button of **History**
lists the same keys.

| Keys | Do |
|---|---|
| `↑` `↓` | move between rows (`Shift`: select several) |
| `→` `←` | open or close a branch, turn, group or run |
| `Enter`, or double-click | show the step |
| `Space`, or `›` | open its details |
| `Alt+↑` `Alt+↓` | show the step before or after |
| `F2` | rename (a branch, a group, a bookmark) |
| `Delete` | delete (asks first) |
| `Shift+F10` | the row's menu |
| `[` `]` | the previous or next fork |
| `Ctrl+G` (`⌘G`) | group the steps selected |
| `Ctrl+Shift+B` (`⇧⌘B`) | branch from the step selected |
| `/` | search ({doc}`05_recording_groups_notes_and_cleanup`) |
| `Ctrl+Z`, `Ctrl+Shift+Z` | undo, redo |

---

## Undo and redo

| | Button | Keyboard |
|---|---|---|
| Undo the step shown | **Undo** (curved arrow left) in the **History** header | `Ctrl+Z`, or `⌘Z` on a Mac |
| Redo the next step | **Redo** (curved arrow right) | `Ctrl+Shift+Z` (`⇧⌘Z` on a Mac), or `Ctrl+Y` on Windows and Linux |

Hover a button to see which step it undoes or redoes ("Undo “Hid Ductal of
cell_type” (Ctrl+Z)"). A screen reader hears "Undid “…”" or "Redid “…”" once
the history has moved.

The keys are left to whatever owns them first:

- while you type in a text field, the field keeps its own undo (a slider, a
  checkbox or a menu has none, so the keys undo the history there too);
- an open dialog or menu keeps its keys;
- while a selection is in progress (the selection bar shows its steps), the
  keys undo the selection's own last step, as the bar's **Undo** does.

Undo and redo are off while a tour plays. While recording is paused with
changes, going back first asks whether to keep them
({ref}`history-recording-paused`).

**Acting after an undo starts a branch.** The steps you undid are not thrown
away: they stay on a branch of their own, which you can show again from the
graph or the branch picker. See {doc}`02_branches_trackers_and_comparing`.

---

## Stale results

A step's result is kept with the inputs it read: an analysis's cells,
filters and genes, the cells a highlight of a category or a range took from
the view, the values a copied field or a kept view was made from. When a
later step changes one of them, the result is **stale**: the status line
shows `⚠ 1 stale`, and the step's details say which step changed what
(`⚠ Stale since step 7 (“Hid Alpha of cell_type”): Filters changed.`).
Renaming a category or a field never makes a result stale. Click `⚠ 1 stale`
to see every stale step.

**Run again**, in the details, brings the result up to date as a new step.
It is offered only where doing the same again does that:

- an analysis is run again on what it reads now;
- a highlight made from a category or a range of values (from the legend, an
  agent or a notebook) gives its group, already on the page, the cells the
  view shows now. The group is updated in place, not added a second time, and
  the new step reads *Updated the cells of …*. If the view shows the same
  cells, nothing changes; if it shows none of them, it is refused and the
  group keeps its cells.

Other stale steps (a copied field, a kept view, a keyframe) offer no **Run
again**: doing them again would add a second one of what they made. Make it
again yourself if you want it up to date.

---

## Delete one step

To take out one step but keep everything after it:

1. Open the step's menu (right-click its row) and choose **Delete…**, or
   select the row and press `Delete`.
2. The question says what will happen: the step goes to **Recently deleted**,
   where it can be restored; how many later steps changed or read what it
   changed and are applied again without it (and roughly how long that takes
   when it is long); which branches start from it (they will start from the
   step before); and which step you will see if it was the step shown.
3. Confirm **Delete step**.

Only the later steps that changed or read what the deleted step changed are
applied again, in order, each on what is now before it. If one cannot be (for
example, it hid a category the deleted step had created), it is not dropped:
it stays in the graph, marked **No longer applies:** with the reason, and its
details open. A notification says what happened, with **Undo**
(`Deleted “Hid Ductal of cell_type”; 1 step no longer applies`).

The opening cannot be deleted. **Undo** in the notification, or **Restore** in
**⋯ › Recently deleted…**, puts the step back exactly where it was
({ref}`history-recently-deleted`).

Two other ways to take a step's effect out:

- **Revert** adds a new step that puts back what the step changed; the step
  itself stays. It is greyed out, saying why, when a later step on the path
  shown changed the same part of the state, when the step changed nothing, or
  when it is not on the path to the step shown.
- **Branch without this step** starts a new branch from the step before it and
  applies the later steps again without it; the branch you were on stays as
  it was.

---

## Steps of agents and notebooks

An AI agent you allow to work in Cellucid (in the chat bar, or connected from
Claude, ChatGPT or a local tool such as Claude Code), and a Python or R
notebook that drives the viewer, add their steps to the same history,
labelled with who took them:

- An agent's steps from one request are one **turn** (`◆`): a row with the
  agent's name and model and how many steps it took, its steps under it. Its
  details say what was asked and when, hold its summary note, and offer
  **Undo turn**, **Merge into one step…** and, for a turn of the chat bar,
  **Show in the chat**. A message in the chat shows its steps here in turn.
- **Undo turn** goes back to the step before the turn when nothing came after
  it. When you or another agent acted after the turn began, it asks how:
  **Only Claude's 3 steps** (they go to Recently deleted, and the later steps
  that needed them are applied again without them; one that cannot be stays,
  marked as no longer applying), or **All 5 steps since** (back to before the
  turn, kept for redo).
- While recording is paused, an agent asks before every change, whatever your
  permission settings allow, because a change made while paused cannot be
  undone on its own, and its activity line starts with `⏸`.
- When History settings leave camera moves (or panels, windows and theme) out
  of the history, every agent's tool list says so ("Camera moves are not
  recorded in this history"), so it knows such a change makes no step and
  Undo does not take it back.
- A notebook's steps carry **Python** or **R**. A notebook reads the history
  and moves it the way the buttons do: `viewer.history_get()`,
  `viewer.history_undo()`, `viewer.history_redo()` and
  `viewer.history_jump(step=...)` in Python
  ({doc}`../../python_package/e_jupyter_hooks/06_python_to_frontend_commands`),
  `cellucid_history_get()`, `cellucid_history_undo()`,
  `cellucid_history_redo()` and `cellucid_history_jump()` in R
  ({doc}`../../r_package/d_viewing_loading/06_hooks_selection_and_sessions`).
  These move the history and add no step. Every other operation of this
  chapter (branches, groups, notes, Recently deleted, sharing, replay) is a
  `history_*` method in Python and a `cellucid_history_*()` function in R
  too.

A history you load from a session file keeps who took each step.

---

## Edge cases

- **A step marked "Not recorded".** A row badged *Not recorded*, reading
  *Changed outside the history*, means a part of the state changed without a
  step. Cellucid keeps the change as its own step so undo still works. It
  should not happen: if you see it, please report it with what you did just
  before.
- **Long histories.** The graph draws only the rows in view, and each distinct
  state is kept once, however many steps show it. The **History memory**
  dialog lists what it holds but cannot yet free a step's results to compute
  later, or merge steps from there: use **Merge into one step…** or empty
  Recently deleted ({doc}`05_recording_groups_notes_and_cleanup`).
- **Another dataset.** Opening another dataset ends this history in the tab.
  It is kept in this browser while **Autosave in this browser** is on, and
  reopening the dataset offers to continue it
  ({doc}`03_saved_states_and_autosave`).
- **The same dataset on other cells.** Changing **Cells to open** keeps the
  history: it is a step of it.

---

## Next

- {doc}`02_branches_trackers_and_comparing`: branches, the branch map,
  comparing, trackers.
- {doc}`03_saved_states_and_autosave`: what this browser keeps, Load…, and
  Continue.
- {doc}`05_recording_groups_notes_and_cleanup`: the recording switch, History
  settings, groups, notes, search and Recently deleted.
