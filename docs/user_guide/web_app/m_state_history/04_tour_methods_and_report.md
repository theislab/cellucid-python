# Tour, pseudocode, methods and report

**Audience:** everyone; especially anyone presenting, teaching or writing up
**Time:** 10 minutes
**What you'll learn:**
- How to play your history back as a guided tour
- How to copy it as pseudocode or as a methods paragraph
- How to print a report with every step, its parameters and pictures

All four live in the **⋯** menu of the **State › History** header. Each uses
the steps of the branch shown, from the opening of the dataset to the step
shown. To use another branch, show it first, or use **Copy pseudocode** in that
branch's own menu ({doc}`02_branches_trackers_and_comparing`).

```{figure} ../../../_static/screenshots/state_history/history-menu.png
:alt: The History block with its ⋯ menu open: Play as a tour, Compare two steps…, Branch map…, Recently deleted…, Copy pseudocode, Copy methods, Print report…, Save AnnData with this state…, Replay a saved history…, Start another tracker…, Show the state at… and Clear history, keep the state….
:width: 496px

The **⋯** menu of **History**.
```

---

## Play as a tour

**⋯ › Play as a tour** plays the steps back one at a time, from the start, for
a talk or a class:

- a bar over the top of the view shows each step's sentence as a caption,
  the step's first note as a second line, the name of a group as a chapter
  title on the group's first step, who took the step when it was not you, and
  where you are (`3 / 12`);
- a camera step is flown to smoothly, from the camera before it; every other
  step is shown at once and held for a moment. While camera moves are not
  recorded ({ref}`history-settings`), the tour still flies to the camera each
  step kept;
- **Pause** holds the step on screen (press it again to play on),
  **Previous** and **Next** step by hand, and **Stop** (or `Escape`) ends the
  tour.

After the last step the bar says *End of the tour*, and Cellucid returns to the
step you were on. If you (or an agent) act during the tour, by a change, an
undo or a jump, the tour ends there and stays where you acted, so the new step
is yours to keep; a note or a bookmark added meanwhile lets it play on. While recording is paused with changes, the tour
first asks whether to keep them ({ref}`history-recording-paused`). A tracker
loaded from a file or made by a replay also offers **▶ Tour** in the line
above the graph.

With **reduce motion** turned on in your system settings, nothing flies: each
step is shown at once and held. Undo and redo are off while the tour plays.

An agent can play and stop the tour too (`history.tour.start`,
`history.tour.stop`; {doc}`../g_agents/index`).

For a scripted camera flight of your own, use the camera path
({doc}`../c_core_interactions/08_camera_path`).

---

## Copy pseudocode

**⋯ › Copy pseudocode** copies the steps as short lines, one per step, to
paste into a message, a notebook or a lab book. For example:

```text
# Cellucid 0.9.7 · Pancreatic endocrinogenesis
# Tracker "Tracker 1", branch "Main": 6 steps
color_by("cell_type")
# --- QC ---
hide_categories("cell_type", ["Ductal", "Ngn3 low EP"])
keep_cells("cell_type", "Beta", 500, seed=0)
# Bookmark: Beta cells only
# --- end of QC ---
reset_camera("live")
# Changes made with recording paused: Filters, Colors (values only)
color_by(gene="Ins1")  # Claude
```

- The first lines name the Cellucid version, the dataset, the tracker, the
  branch and the number of steps.
- Each step is the line it was recorded with. Steps name fields, categories
  and genes by name, never by position.
- A step taken by an agent or a notebook says so at the end (`# Claude`,
  `# Python`).
- The steps of a group sit between `# --- QC ---` and `# --- end of QC ---`.
  A bookmarked step is followed by `# Bookmark: <name>`. When steps were
  merged into one, the merged step keeps their bookmarks, each followed by
  `(on a step merged into this one)`.
- A step that no longer applies, a change made outside the history, and a
  step that holds only values are kept as `#` comments, so nothing is left out
  silently. A step that holds only values has no line of its own: the changes
  made while recording was paused, several steps merged into one, a fresh
  start after the history was cleared, or a state shared without its history.

The lines describe what was done; they are not a script to run. To repeat the
work on another dataset, replay the history there ({doc}`06_sharing_and_replay`),
keep the session file (it reopens every step), or drive the viewer from
Python or R ({doc}`../../python_package/e_jupyter_hooks/06_python_to_frontend_commands`).

---

## Copy methods paragraph

**⋯ › Copy methods** copies a paragraph written the way a methods section
is: past tense, in order, with the software, its version and the dataset. For
example:

> Single-cell data were explored interactively in Cellucid (version 0.9.7) [citation]. The dataset “Pancreatic endocrinogenesis” (… cells and … genes; its citation) was opened in full. The following 4 steps were applied, in this order: (1) colored by cell_type; (2) hid Ductal and Ngn3 low EP of cell_type; (3) kept 500 cells of Beta in cell_type; (4) colored by gene Ins1. Step 4 was taken by an AI assistant (Claude, claude-sonnet-4-6) at the user's request. 1 further step only changed the camera, the panel layout, colors or display settings and is not listed. The complete history, with every parameter, can be reopened from the accompanying Cellucid session file.

What it holds:

- the tracker's first note, its abstract, as the opening paragraph, when the
  tracker has one;
- how to cite Cellucid, and the dataset's name, size and citation when the
  dataset gives one;
- the cells opened: in full, or how many of how many, chosen how and with
  which seed ({doc}`../b_data_loading/21_subsampled_loading`);
- every step that changed the data, its labels, what is shown of it or what
  was computed (fields and labels, filters, hidden cells, highlights, the data
  mapping, gene sets, analyses), numbered in order;
- which steps an agent (with its model, when known) or a notebook took;
- a step that holds only values, said as what it is, with the parts of the
  state it changed: "changes made with recording paused (filters, colors)",
  "3 steps merged into one (highlights)", "a fresh start from the state
  reached before the history was cleared (…)" or "a shared state, received
  without its history (…)". A sentence says that only the state it reached
  was kept, not the operations;
- how many steps only changed the presentation (camera, panels, colors,
  display settings): counted, not listed;
- how many steps no longer applied after a step was deleted: counted, left
  out.

When a step it lists is in a group, the text is split into sections. The first
paragraph says how many steps were applied and holds the counts; then each
group's steps follow under the group's name on a line of its own, and steps in
no group under **Before *group*** or **After *group***:

```text
Single-cell data were explored interactively in Cellucid (version 0.9.7). The dataset … was opened in full. The 4 steps applied are listed below, in order. 1 further step only changed the camera, … The complete history, with every parameter, can be reopened from the accompanying Cellucid session file.

Before QC
The following 1 step was applied: (1) colored by cell_type.

QC
The following 2 steps were applied, in this order: (2) hid Ductal and Ngn3 low EP of cell_type; (3) kept 500 cells of Beta in cell_type.

After QC
The following 1 step was applied: (4) colored by gene Ins1. Step 4 was taken by an AI assistant (Claude, claude-sonnet-4-6) at the user's request.
```

Notes on single steps are not part of it; they are in the report. Read it
before you paste it: it says what was done, not why.

---

## Print report

**⋯ › Print report…** opens a dialog:

- **Pictures of the steps chosen** (on by default): a picture of the selected
  view after each step you tick. The bookmarked steps, the last step of each
  group and the step shown are ticked; tick others to show how the view
  changed. This browser remembers the switch for the next report.
- **Open report** makes the report in a new window, to read, print or save as
  PDF. Allow pop-ups for the site if the browser blocks it.

```{figure} ../../../_static/screenshots/state_history/print-report-dialog.png
:alt: The Print report dialog: the switch Pictures of the steps chosen (on), the steps of the history with the newest ticked, and Open report.
:width: 960px

**Print report** with **Pictures of the steps chosen** on and the step shown ticked.
```

The report holds:

| Section | What is in it |
|---|---|
| Header | *History report*, the dataset, the Cellucid version and when the report was made |
| **Abstract** | The tracker's first note, when it has one |
| **Notes** | The tracker's other notes and the branch's notes, each with who wrote it and when |
| **Dataset** | Its name, identity, cells and genes, the cells open (when a subsample), its citation, the tracker, the branch, the number of steps, the Cellucid version |
| **Methods** | The methods text above, with its group headings |
| **Steps** | Every step, under its group's name (with the group's notes) when the branch has groups: its sentence, its bookmarks (those it kept from steps merged into it marked *merged*), who took it (everyone who acted, for changes made while recording was paused), when, its parameters, the parts of the state it changed, a flag when it no longer applies, was made outside the history or holds only values, its notes (those kept from steps merged into it say so), and its picture when one was chosen |
| **Pseudocode** | The pseudocode above |

Pictures are never stored in the history: each one is drawn when the report is
made, and only for the steps picked.

Each picture is made by showing that step, drawing the selected view the way
the figure export draws it ({doc}`../k_figure_export/index`), then returning
to the step you were on. Pictures take a moment per step; turn the switch off
for a report of the steps alone. Where the page cannot draw pictures, the
report opens without the dialog and says so.

**An agent asks for the report** with `history.report`. A browser opens a
new window only from a click, so an agent gets the same document as an
`.html` file downloaded to your computer (asked first: it writes a file), or
as its text, without pictures, to read itself
({doc}`../g_agents/03_permissions_and_what_agents_see`).

---

## What to keep with a paper

For a figure or an analysis you report:

1. A session file (**State › Session › Share…**, **History**): the branches
   you report, which reopen every step on the same dataset
   ({doc}`06_sharing_and_replay`).
2. The methods paragraph, edited into your methods section.
3. The printed report (as PDF) in your supplementary material or lab book.

---

## Next

- {doc}`01_the_history`: what a step is, undo and redo.
- {doc}`../l_sessions_sharing/05_share_workflows_links_sessions_data_files`: sharing with collaborators.
