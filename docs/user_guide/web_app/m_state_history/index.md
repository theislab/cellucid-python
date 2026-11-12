# State and History

The **State** section of the sidebar keeps the history of everything you do
with the dataset shown. Every change you make is one **step**: what you color
by, a filter, a renamed category, a highlight, a kept view, a camera move, a
panel you dock. You can go back to any step, try something else from there,
keep several independent lines of work, and come back to all of it later.

These pages are written for everyone. You do not need to know anything about
files or code to use the history; the last page shows how to turn it into
pseudocode, a methods paragraph or a printed report for a paper or a lab book.

:::{important}
The history belongs to the dataset shown. Open another dataset and a new
history starts; the one you leave is kept in this browser (see
{doc}`03_saved_states_and_autosave`), a session file carries it to another
computer (see {doc}`06_sharing_and_replay`), and its steps can be replayed on
another dataset.
:::

---

## Mental model (one sentence)

Cellucid writes down every change you make as a step in a tree of steps, so
you can show any step again, branch off from it, keep it in this browser,
share it as a file, replay it on another dataset, and say exactly what you
did.

---

## Fast path (choose what you want to do)

| You want to… | Do this | Read |
|---|---|---|
| Undo the last thing you did | Press `Ctrl+Z` (`⌘Z` on a Mac), or **Undo** in the **History** header | {doc}`01_the_history` |
| Go back to an earlier point | Double-click that step in the graph | {doc}`01_the_history` |
| See what a step did | `›` at the end of its row: its details list who took it, its parameters and what it changed | {doc}`01_the_history` |
| Remove one step but keep what came after | Right-click the step, **Delete…** | {doc}`01_the_history` |
| Try an alternative without losing your work | Undo, then act: the steps you undid stay on their own branch | {doc}`02_branches_trackers_and_comparing` |
| Find and manage branches | **Main ▾** on the status line (the branch picker), or **⋯ › Branch map…** | {doc}`02_branches_trackers_and_comparing` |
| Look at two alternatives side by side | A branch's `⋯`, **Compare with the shown branch** | {doc}`02_branches_trackers_and_comparing` |
| See what changed between two points | **⋯ › Compare two steps…** | {doc}`02_branches_trackers_and_comparing` |
| Start separate work on the same dataset | **⋯ › Start another tracker…** | {doc}`02_branches_trackers_and_comparing` |
| Pick up tomorrow where you stopped | Open the dataset again and choose **Continue** | {doc}`03_saved_states_and_autosave` |
| Keep a version under a name | **State › Session › Load…**, type a name, **Save** | {doc}`03_saved_states_and_autosave` |
| Try things without recording them | Click **● Recording** to pause; resume to keep the changes as one step | {doc}`05_recording_groups_notes_and_cleanup` |
| Tidy a long run of steps | Select it (`Shift`+click), then **Group…** or **Merge into one step…** | {doc}`05_recording_groups_notes_and_cleanup` |
| Get back a deleted step or branch | **Undo** in its notification, or **⋯ › Recently deleted…** | {doc}`05_recording_groups_notes_and_cleanup` |
| Send your work to someone | **State › Session › Share…**, **Save file**, and send the file | {doc}`06_sharing_and_replay` |
| Save the dataset with your labels, highlights and gene sets written in | **⋯ › Save AnnData with this state…**: a new `.zarr.zip` that Python and R read | {doc}`06_sharing_and_replay` |
| Run the same workflow on another sample | **Load… › Replay here** (or **Replay a session file here…**) with the other sample open | {doc}`06_sharing_and_replay` |
| Show your analysis in a talk | **⋯ › Play as a tour** | {doc}`04_tour_methods_and_report` |
| Write the methods section | **⋯ › Copy methods** | {doc}`04_tour_methods_and_report` |
| Keep a printed record | **⋯ › Print report…** | {doc}`04_tour_methods_and_report` |

---

## The State section at a glance

The **State** section sits in the sidebar right below **Session**. From top to
bottom:

| Part | What it holds |
|---|---|
| **Session** | **Share…** writes a `.cellucid-session` file of the history, the current state only or the steps only; **Load…** lists every history this browser keeps (open, replay here, rename, delete, save the history shown under a name) and opens or replays a session file. |
| **History** header | Search (`/`), **Undo**, **Redo**, **History settings** (what is recorded, autosave, your name) and **⋯**: **Play as a tour**, **Compare two steps…**, **Branch map…**, **Recently deleted…**, **Copy pseudocode**, **Copy methods**, **Print report…**, **Save AnnData with this state…**, **Replay a saved history…**, **Start another tracker…**, **Show the state at…**, **Clear history, keep the state…**. |
| Status line | **● Recording** (the recording switch), the branch shown (**Main ▾**, the branch picker with trackers and the branch map), the step count (`14 steps · 5 undone`), and `⚠ stale` when a step's result (an analysis, a highlight of a category, a copied field, …) no longer matches its inputs. |
| The graph | One row per step on the path to the step shown, with the branches that leave it, turns of agents, groups and folded runs. Click selects, double-click shows, `›` opens the details. |

Each part has an `i` button that says the same in a few lines.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The State section with the Pancreatic endocrinogenesis sample open: Session with Share… and Load…, and the History with the recording switch, the branch shown and the graph of the steps taken so far, the step shown marked.
:width: 540px

The **State** section after three steps on the Pancreas sample: **Session** with **Share…** and **Load…**, and the **History** graph, the step shown marked.
```

---

## Where your history is kept

| Where | What it holds | How long |
|---|---|---|
| This tab | The history you are working on: every step of every tracker of the dataset shown | Until you open another dataset or close the tab |
| This browser | Every step, written as you make it (**Autosave in this browser**, on by default), and the states you saved under a name | Until you delete them (**State › Session › Load…**, or **Session › Open data › Stored › Manage…**), or clear the site's data |
| A session file | What you chose in **Share…**: the history of the trackers and branches ticked, the current state only, or the steps only | As long as you keep the file; it opens on any computer and browser with the same dataset, and its steps replay on others |

Nothing is uploaded. The dataset itself is never part of the history or the
session file. What this browser keeps, and how to remove it, is listed in
{doc}`../o_accessibility_privacy_security/02_privacy_model`.

---

## Settings for the heavier parts

The **History settings** box (the sliders button in the **History** header)
holds the history's switches, remembered by this browser; everything is on
by default ({ref}`history-settings`):

| Switch | Off means |
|---|---|
| **Camera moves**, **Panels, windows and theme** | Those changes make no step and the history never moves them |
| **Record history for new datasets** | A dataset you open starts with recording paused |
| **Ask before going back while paused** | Changes made while paused are kept as one step without asking |
| **Autosave in this browser** | Nothing more is written to this browser; saving under a name still works ({doc}`03_saved_states_and_autosave`) |

**Pictures of the steps chosen**, in **⋯ › Print report…**, decides whether
the report draws pictures ({doc}`04_tour_methods_and_report`).

Going back to a step never computes anything again: every state a step showed
is kept, so there is nothing else to switch.

---

## Who can add steps

Steps are taken by you, and also by:

- an AI agent working in Cellucid for you, in the chat bar or connected from
  Claude, ChatGPT or a local tool (each of its steps carries its name and
  model, such as **Claude · Sonnet 4.6**, and **Undo turn** undoes its whole
  turn);
- a Python or R notebook that drives the viewer (its steps carry **Python** or
  **R**).

They all add to the same history, so you always see who did what. A notebook
reads and moves the same history: in Python `viewer.history_get()`,
`viewer.history_undo()`, `viewer.history_redo()` and
`viewer.history_jump(step=...)` ({doc}`../../python_package/e_jupyter_hooks/06_python_to_frontend_commands`);
in R `cellucid_history_get()`, `cellucid_history_undo()`,
`cellucid_history_redo()` and `cellucid_history_jump()`
({doc}`../../r_package/d_viewing_loading/06_hooks_selection_and_sessions`).

---

## Pages in this chapter

1) {doc}`01_the_history`: what a step is, the history graph, showing a step, undo and redo, deleting a step, steps of agents and notebooks.
2) {doc}`02_branches_trackers_and_comparing`: branches, the branch picker and map, comparing branches side by side, comparing two steps, trackers.
3) {doc}`03_saved_states_and_autosave`: every step kept in this browser, Continue where you left off, the Load dialog, turning autosave off.
4) {doc}`04_tour_methods_and_report`: the guided tour, pseudocode, the methods paragraph and the printable report.
5) {doc}`05_recording_groups_notes_and_cleanup`: pausing recording, History settings, groups and merging, notes and bookmarks, search, Recently deleted.
6) {doc}`06_sharing_and_replay`: Share…, saving an AnnData with this state, loading a shared history, replaying on another dataset, applying a branch, Verify.

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
